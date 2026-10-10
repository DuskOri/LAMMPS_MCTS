"""LAMMPS_MCTS 本地前端接口服务。

这个服务只面向本机使用，负责三件事：
1. 为 frontend/ 提供静态页面；
2. 修改 config.yaml 中常用运行参数；
3. 在后台按阶段执行主流程，并把状态返回给前端。

启动方式：
    python web_server.py

浏览器访问：
    http://localhost:8010/frontend/
"""

import csv
import json
import os
import shutil
import threading
import time
import secrets
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse, unquote

import main
from md_engine import terminate_active_lammps
from md_engine.process_control import PipelineCancelled, reset_cancellation, cancel_processes
from md_engine.thermal_conductivity_writer import _normalise_params
from mcts.fragment_registry import FragmentRegistry
from run_storage import atomic_file, write_json, new_run_id


PROJECT_ROOT = Path(__file__).resolve().parent
HOST = "127.0.0.1"
PORT = int(os.environ.get("LAMMPS_MCTS_PORT", "8010"))
SESSION_TOKEN = secrets.token_urlsafe(32)

THERMAL_PROFILES = {
    "quick": {
        "system": {"molecule_count": 4, "box_size": 45.0},
        "thermal_conductivity": {
            "density_equilibration": True,
            "melt_temperature": 550.0,
            "melt_heating_steps": 2000,
            "melt_hold_steps": 2000,
            "precompression_density": 0.70,
            "precompression_steps": 5000,
            "density_block_steps": 2000,
            "density_sample_every": 100,
            "density_max_blocks": 10,
            "density_plateau_tolerance": 0.05,
            "cooling_steps": 3000,
            "nvt_steps": 2000,
            "npt_steps": 8000,
            "production_steps": 100000,
            "nemd_hot_temperature": 450.0,
            "nemd_cold_temperature": 150.0,
            "nemd_thermostat_damp": 1.0,
            "nemd_steady_steps": 20000,
            "nemd_production_block_steps": 10000,
            "nemd_min_production_blocks": 5,
            "nemd_heat_flux_tolerance": 0.15,
            "nemd_profile_nevery": 10,
            "nemd_profile_repeat": 100,
            "nemd_min_gradient_r2": 0.70,
            "thermo_every": 500,
            "dump_every": 5000,
        },
        "lammps": {
            "timeout_seconds": 1800,
            "compression_timeout_seconds": 1800,
            "nemd_timeout_seconds": 1800,
        },
    },
    "standard": {
        "system": {"molecule_count": 10, "box_size": 60.0},
        "thermal_conductivity": {
            "density_equilibration": True,
            "melt_temperature": 550.0,
            "melt_heating_steps": 10000,
            "melt_hold_steps": 20000,
            "precompression_density": 0.70,
            "precompression_steps": 50000,
            "density_block_steps": 10000,
            "density_sample_every": 100,
            "density_max_blocks": 20,
            "density_plateau_tolerance": 0.02,
            "cooling_steps": 50000,
            "nvt_steps": 10000,
            "npt_steps": 50000,
            "production_steps": 500000,
            "nemd_hot_temperature": 450.0,
            "nemd_cold_temperature": 150.0,
            "nemd_thermostat_damp": 1.0,
            "nemd_steady_steps": 50000,
            "nemd_production_block_steps": 25000,
            "nemd_min_production_blocks": 5,
            "nemd_heat_flux_tolerance": 0.08,
            "nemd_profile_nevery": 10,
            "nemd_profile_repeat": 100,
            "nemd_min_gradient_r2": 0.85,
            "thermo_every": 1000,
            "dump_every": 10000,
        },
        "lammps": {
            "timeout_seconds": 7200,
            "compression_timeout_seconds": 3600,
            "nemd_timeout_seconds": 7200,
        },
    },
}

RUN_LOCK = threading.Lock()
RUN_CANCEL_EVENT = threading.Event()
RUN_SERIES_FILES = {"density": "", "heat_flux": ""}


class RunCancelled(PipelineCancelled):
    """Raised inside the background pipeline after a user cancellation."""

RUN_STATE = {
    "running": False,
    "stage": "idle",
    "message": "等待运行",
    "started_at": "",
    "finished_at": "",
    "error": "",
    "log": [],
    "console": [],
    "iteration": 0,
    "total_iterations": 0,
    "evaluation_count": 0,
    "top_k": 5,
    "live_candidates": [],
    "telemetry": {"evaluation": 0, "density": [], "heat_flux": []},
}

STAGE_LABELS = {
    "idle": "等待运行",
    "search": "MCTS 搜索",
    "mcts_lammps": "MCTS 热导反馈",
    "build": "RDKit 建链",
    "packmol": "Packmol 初始体系",
    "write_lammps": "写出 LAMMPS 输入",
    "compress": "高温熔融与体系压缩",
    "run_lammps": "NEMD 热导计算",
    "export_build": "导出 Top K 结构",
    "export_packmol": "导出 Top K 体系",
    "export_lammps": "导出 Top K 输入",
    "save": "保存结果",
    "done": "完成",
    "failed": "失败",
    "stopping": "正在终止",
    "stopped": "已终止",
}


class LammpsMctsHandler(SimpleHTTPRequestHandler):
    """处理静态页面和 API 请求。"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(PROJECT_ROOT / "frontend"), **kwargs)

    def end_headers(self):
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def _allowed_request(self, mutation=False):
        host = self.headers.get("Host", "")
        allowed = {f"127.0.0.1:{self.server.server_port}", f"localhost:{self.server.server_port}"}
        origin = self.headers.get("Origin")
        if host not in allowed or (origin is not None and origin != f"http://{host}"):
            self.send_error(HTTPStatus.FORBIDDEN, "Local same-origin requests required")
            return False
        if mutation and not secrets.compare_digest(self.headers.get("X-CSRF-Token", ""), SESSION_TOKEN):
            self.send_error(HTTPStatus.FORBIDDEN, "Session token required")
            return False
        return True

    def do_OPTIONS(self):
        if not self._allowed_request():
            return
        self.send_response(HTTPStatus.NO_CONTENT)
        self.end_headers()

    def do_GET(self):
        if not self._allowed_request():
            return
        path = urlparse(self.path).path
        if path == "/api/session":
            self._send_json({"token": SESSION_TOKEN})
            return
        if path == "/":
            self.send_response(HTTPStatus.FOUND)
            self.send_header("Location", "/frontend/")
            self.end_headers()
            return
        if path == "/api/config":
            self._send_json(load_public_config())
            return
        if path == "/api/graphs":
            self._send_json(load_graph_info())
            return
        if path == "/api/status":
            refresh_telemetry()
            self._send_json({**RUN_STATE, "can_stop": True, "counts": collect_counts()})
            return
        if path == "/api/results":
            self._send_json(load_results())
            return
        decoded = unquote(path)
        if not decoded.startswith("/frontend/"):
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        relative = decoded[len("/frontend/"):]
        asset = (PROJECT_ROOT / "frontend" / relative).resolve()
        root = (PROJECT_ROOT / "frontend").resolve()
        if not asset.is_relative_to(root) or (relative and asset.is_dir()):
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        self.path = "/" + relative
        super().do_GET()

    def do_HEAD(self):
        self.send_error(HTTPStatus.METHOD_NOT_ALLOWED)

    def do_POST(self):
        if not self._allowed_request(mutation=True):
            return
        path = urlparse(self.path).path
        try:
            payload = self._read_json_body()
            if path == "/api/config":
                with RUN_LOCK:
                    config = update_config(payload)
                self._send_json({"ok": True, "config": public_config_view(config)})
                return
            if path == "/api/custom-graph":
                with RUN_LOCK:
                    saved = save_custom_graph(payload)
                self._send_json({"ok": True, "path": str(saved)})
                return
            if path == "/api/run":
                started = start_background_run(payload)
                self._send_json({"ok": started, "status": RUN_STATE})
                return
            if path == "/api/stop":
                result = stop_background_run()
                self._send_json({"ok": result["requested"], **result, "status": RUN_STATE})
                return
            if path == "/api/clear-data":
                with RUN_LOCK:
                    result = clear_output_data(payload)
                self._send_json({"ok": True, **result, "counts": collect_counts()})
                return
            self.send_error(HTTPStatus.NOT_FOUND, "API endpoint not found")
        except Exception as exc:
            self._send_json({"ok": False, "error": str(exc)})

    def _read_json_body(self):
        length = int(self.headers.get("Content-Length", "0") or 0)
        if length > 1024 * 1024:
            raise ValueError("request body is too large")
        if length <= 0:
            return {}
        raw = self.rfile.read(length).decode("utf-8")
        return json.loads(raw) if raw.strip() else {}

    def _send_json(self, data):
        body = json.dumps(data, ensure_ascii=False, indent=2).encode("utf-8")
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def load_public_config():
    """读取配置并返回前端需要的精简视图。"""
    config = main.load_config(PROJECT_ROOT / "config.yaml")
    return public_config_view(config)


def public_config_view(config):
    """把完整配置压缩成前端表单字段。"""
    mcts_config = config.get("mcts", {})
    polymer_config = config.get("polymer", {})
    chemistry_config = config.get("chemistry", {})
    tc_config = config.get("thermal_conductivity", {})
    lammps_config = config.get("lammps", {})
    try:
        effective = _normalise_params(tc_config)
    except ValueError:
        effective = dict(tc_config)
    return {
        "effective_parameters": {
            "力场": config.get("system", {}).get("force_field", "uff_screening"),
            "链数": config.get("system", {}).get("molecule_count", 4),
            "热浴 / K": effective.get("nemd_hot_temperature", ""),
            "冷浴 / K": effective.get("nemd_cold_temperature", ""),
            "热浴阻尼 / fs": effective.get("nemd_thermostat_damp", ""),
            "时间步 / fs": effective.get("timestep", ""),
            "稳态步数": effective.get("nemd_steady_steps", ""),
            "采样块步数": effective.get("nemd_production_block_steps", ""),
            "最低采样块数": effective.get("nemd_min_production_blocks", ""),
            "实际生产步数上限": effective.get("production_steps", ""),
        },
        "iterations": int(mcts_config.get("iterations", 300)),
        "max_steps": int(mcts_config.get("max_steps", 5)),
        "top_k": int(mcts_config.get("top_k", 5)),
        "exploration_weight": float(mcts_config.get("exploration_weight", 1.41421356237)),
        "target_length_angstrom": float(
            polymer_config.get("target_chain_length_angstrom", 120.0)
        ),
        "active_families": ",".join(chemistry_config["active_families"]) if isinstance(chemistry_config.get("active_families"), list) else str(chemistry_config.get("active_families", "polyimide")),
        "run_lammps": bool(mcts_config.get("feedback_run_lammps", False)),
        "feedback_mode": str(mcts_config.get("feedback_mode", "heuristic")),
        "feedback_run_lammps": bool(mcts_config.get("feedback_run_lammps", False)),
        "feedback_every_iteration": bool(mcts_config.get("feedback_every_iteration", True)),
        "feedback_max_evaluations": int(mcts_config.get("feedback_max_evaluations", 0) or 0),
        "thermal_profile": str(tc_config.get("profile", "quick")),
        "density_equilibration": bool(tc_config.get("density_equilibration", True)),
        "precompression_density": float(tc_config.get("precompression_density", 0.70)),
        "density_plateau_percent": float(
            tc_config.get("density_plateau_tolerance", 0.05)
        ) * 100.0,
        "density_max_blocks": int(tc_config.get("density_max_blocks", 10)),
        "compression_timeout_seconds": int(
            lammps_config.get(
                "compression_timeout_seconds",
                lammps_config.get("timeout_seconds", 1800),
            ) or 0
        ),
        "nemd_timeout_seconds": int(
            lammps_config.get(
                "nemd_timeout_seconds",
                lammps_config.get("timeout_seconds", 1800),
            ) or 0
        ),
        "production_steps": int(tc_config.get("production_steps", 100000)),
        "heat_flux_tolerance_percent": float(
            tc_config.get("nemd_heat_flux_tolerance", 0.15)
        ) * 100.0,
    }


def update_config(payload):
    """更新常用配置并写回 config.yaml。"""
    if RUN_STATE["running"]:
        raise RuntimeError("Configuration cannot be changed while a task is running")
    config = main.load_config(PROJECT_ROOT / "config.yaml")
    mcts_config = config.setdefault("mcts", {})
    polymer_config = config.setdefault("polymer", {})
    chemistry_config = config.setdefault("chemistry", {})
    tc_config = config.setdefault("thermal_conductivity", {})
    lammps_config = config.setdefault("lammps", {})

    if "iterations" in payload:
        mcts_config["iterations"] = int(payload["iterations"])
    if "exploration_weight" in payload:
        mcts_config["exploration_weight"] = float(payload["exploration_weight"])
    if "max_steps" in payload:
        mcts_config["max_steps"] = int(payload["max_steps"])
    if "top_k" in payload:
        mcts_config["top_k"] = int(payload["top_k"])
    if "target_length_angstrom" in payload:
        polymer_config["target_chain_length_angstrom"] = float(
            payload["target_length_angstrom"]
        )
    if "active_families" in payload:
        chemistry_config["active_families"] = str(payload["active_families"])
    if "feedback_mode" in payload:
        mcts_config["feedback_mode"] = str(payload["feedback_mode"])
    if "feedback_run_lammps" in payload:
        mcts_config["feedback_run_lammps"] = bool(payload["feedback_run_lammps"])
    if "feedback_every_iteration" in payload:
        mcts_config["feedback_every_iteration"] = bool(payload["feedback_every_iteration"])
    if "feedback_max_evaluations" in payload:
        mcts_config["feedback_max_evaluations"] = int(payload["feedback_max_evaluations"])
    if "thermal_profile" in payload:
        requested_profile = str(payload["thermal_profile"])
        if requested_profile != tc_config.get("profile", "quick") or payload.get("apply_profile") is True:
            apply_thermal_profile(config, requested_profile)
    if "density_equilibration" in payload:
        tc_config["density_equilibration"] = bool(payload["density_equilibration"])
    if "precompression_density" in payload:
        tc_config["precompression_density"] = max(
            0.05, float(payload["precompression_density"])
        )
    if "density_plateau_percent" in payload:
        tc_config["density_plateau_tolerance"] = max(
            0.0001, float(payload["density_plateau_percent"]) / 100.0
        )
    if "density_max_blocks" in payload:
        tc_config["density_max_blocks"] = max(3, int(payload["density_max_blocks"]))
    if "compression_timeout_seconds" in payload:
        lammps_config["compression_timeout_seconds"] = max(
            0, int(payload["compression_timeout_seconds"])
        )
    if "nemd_timeout_seconds" in payload:
        lammps_config["nemd_timeout_seconds"] = max(
            0, int(payload["nemd_timeout_seconds"])
        )
    if "production_steps" in payload:
        tc_config["production_steps"] = int(payload["production_steps"])
    if "heat_flux_tolerance_percent" in payload:
        tolerance_percent = min(
            100.0, max(0.01, float(payload["heat_flux_tolerance_percent"]))
        )
        tc_config["nemd_heat_flux_tolerance"] = tolerance_percent / 100.0

    if mcts_config.get("feedback_mode") == "thermal":
        mcts_config["feedback_every_iteration"] = True
        mcts_config["feedback_max_evaluations"] = max(
            int(mcts_config.get("feedback_max_evaluations", 0) or 0),
            int(mcts_config.get("iterations", 0) or 0),
        )
        if "run_lammps" in payload:
            mcts_config["feedback_run_lammps"] = bool(payload["run_lammps"])

    for key in ("iterations", "max_steps", "top_k"):
        if int(mcts_config[key]) <= 0:
            raise ValueError(f"{key} must be positive")
    if polymer_config["target_chain_length_angstrom"] <= 0:
        raise ValueError("target chain length must be positive")
    _normalise_params(tc_config)
    write_config(config, PROJECT_ROOT / "config.yaml")
    return config


def apply_thermal_profile(config, profile_name):
    """应用快速筛选或标准复核参数。"""
    if profile_name not in THERMAL_PROFILES:
        raise ValueError(f"Unknown thermal profile: {profile_name}")
    name = profile_name
    profile = THERMAL_PROFILES[name]
    for section, values in profile.items():
        config.setdefault(section, {}).update(values)
    config.setdefault("thermal_conductivity", {})["profile"] = name


def write_config(config, filename):
    """JSON is valid YAML and preserves nested graphs without optional parsers."""
    write_json(filename, config)


def load_graph_info():
    """返回前端的有向图选项和当前自定义片段。"""
    custom_path = PROJECT_ROOT / "custom_fragments.json"
    custom_data = {}
    if custom_path.exists():
        custom_data = json.loads(custom_path.read_text(encoding="utf-8"))
    return {
        "presets": [
            {"key": "polyimide", "label": "聚酰亚胺"},
            {"key": "polyurethane", "label": "聚氨酯"},
            {"key": "phenolic", "label": "酚醛树脂"},
            {"key": "custom", "label": "自定义片段"},
        ],
        "custom": custom_data,
    }


def save_custom_graph(payload):
    """保存前端创建的自定义有向图片段。"""
    if RUN_STATE["running"]:
        raise RuntimeError("Graph cannot be changed while a task is running")
    fragments = payload.get("fragments", [])
    transitions = payload.get("transitions", {})
    if not isinstance(fragments, list) or not isinstance(transitions, dict):
        raise ValueError("custom graph payload must contain fragments and transitions")
    data = {"fragments": fragments, "transitions": transitions}
    FragmentRegistry().configure(active_families="custom", custom_fragments=fragments, custom_transitions=transitions)
    path = PROJECT_ROOT / "custom_fragments.json"
    write_json(path, data)
    return path


def start_background_run(payload=None):
    """启动后台运行线程。"""
    with RUN_LOCK:
        if RUN_STATE["running"]:
            return False
        if payload:
            update_config(payload)
        reset_run_state()
        thread = threading.Thread(target=run_pipeline, daemon=True)
        thread.start()
        return True


def stop_background_run():
    """Request cancellation and terminate the active LAMMPS child process."""
    with RUN_LOCK:
        if not RUN_STATE["running"]:
            return {"requested": False, "process_terminated": False}
        RUN_CANCEL_EVENT.set()
        RUN_STATE["stage"] = "stopping"
        RUN_STATE["message"] = "正在终止当前任务"
        RUN_STATE["log"].append(
            {
                "time": time.strftime("%H:%M:%S"),
                "stage": STAGE_LABELS["stopping"],
                "message": RUN_STATE["message"],
            }
        )

    process_terminated = terminate_active_lammps()
    return {"requested": True, "process_terminated": process_terminated}


def reset_run_state():
    """重置运行状态。"""
    RUN_CANCEL_EVENT.clear()
    reset_cancellation()
    RUN_SERIES_FILES.update({"density": "", "heat_flux": ""})
    RUN_STATE.update(
        {
            "running": True,
            "stage": "search",
            "message": "准备开始 MCTS 搜索",
            "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "finished_at": "",
            "error": "",
            "log": [],
            "console": [],
            "iteration": 0,
            "total_iterations": 0,
            "evaluation_count": 0,
            "top_k": 5,
            "live_candidates": [],
            "telemetry": {"evaluation": 0, "density": [], "heat_flux": []},
        }
    )


def set_stage(stage, message, details=None):
    """更新当前阶段，并接收 MCTS 的轮次和实时候选事件。"""
    details = dict(details or {})
    series_files = details.get("series_files") or {}
    if series_files:
        for name in ("density", "heat_flux"):
            if name in series_files:
                RUN_SERIES_FILES[name] = str(series_files[name] or "")
        RUN_STATE["telemetry"]["evaluation"] = int(
            details.get("telemetry_evaluation", RUN_STATE["evaluation_count"] + 1)
        )
        refresh_telemetry()
    if (
        RUN_CANCEL_EVENT.is_set()
        and stage not in ("stopping", "stopped")
        and not details.get("console_line")
    ):
        raise RunCancelled("run cancelled by user")
    console_line = str(details.get("console_line", "")).rstrip()
    if console_line:
        RUN_STATE["console"].append(
            {
                "time": time.strftime("%H:%M:%S"),
                "text": console_line,
            }
        )
        RUN_STATE["console"] = RUN_STATE["console"][-200:]

    if not console_line:
        label = STAGE_LABELS.get(stage, stage)
        print(f"[PIPELINE][{label}] {message}", flush=True)

    if details.get("console_only"):
        return

    if "iteration" in details:
        RUN_STATE["iteration"] = int(details["iteration"])
    if "total_iterations" in details:
        RUN_STATE["total_iterations"] = int(details["total_iterations"])

    candidate = details.get("candidate")
    if candidate:
        upsert_live_candidate(candidate)

    if stage == "candidate":
        RUN_STATE["evaluation_count"] = max(
            RUN_STATE["evaluation_count"], len(RUN_STATE["live_candidates"])
        )
        return
    if stage == "evaluation":
        RUN_STATE["evaluation_count"] = int(
            details.get("evaluation", RUN_STATE["evaluation_count"])
        )
        RUN_STATE["stage"] = "mcts_lammps"
        RUN_STATE["message"] = message
        return

    RUN_STATE["stage"] = stage
    RUN_STATE["message"] = message
    RUN_STATE["log"].append(
        {
            "time": time.strftime("%H:%M:%S"),
            "stage": STAGE_LABELS.get(stage, stage),
            "message": message,
        }
    )
    RUN_STATE["log"] = RUN_STATE["log"][-50:]


def refresh_telemetry(max_points=240):
    """读取当前候选的密度和热流文件，供前端轮询绘图。"""
    telemetry = RUN_STATE.setdefault(
        "telemetry", {"evaluation": 0, "density": [], "heat_flux": []}
    )
    for name in ("density", "heat_flux"):
        filename = RUN_SERIES_FILES.get(name)
        telemetry[name] = read_lammps_series(filename, max_points=max_points)


def read_lammps_series(filename, max_points=240):
    """读取 fix ave/time 的两列输出，返回时间步和值。"""
    if not filename:
        return []
    path = Path(filename)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    if not path.exists():
        return []

    points = []
    try:
        for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            try:
                points.append({"step": float(fields[0]), "value": float(fields[-1])})
            except ValueError:
                continue
    except OSError:
        return []

    limit = max(int(max_points), 2)
    if len(points) <= limit:
        return points
    last_index = len(points) - 1
    return [
        points[round(index * last_index / (limit - 1))]
        for index in range(limit)
    ]


def upsert_live_candidate(candidate):
    """按序列更新实时候选，并始终保持 reward 从高到低排列。"""
    sequence = list(candidate.get("sequence") or [])
    if not sequence:
        return

    sequence_text = " -> ".join(str(item) for item in sequence)
    rows = RUN_STATE["live_candidates"]
    row = next((item for item in rows if item["sequence_text"] == sequence_text), None)
    if row is None:
        row = {"sequence": sequence, "sequence_text": sequence_text}
        rows.append(row)

    for key, value in candidate.items():
        if key != "sequence":
            row[key] = value
    row["iteration"] = int(candidate.get("iteration", RUN_STATE["iteration"]))

    rows.sort(key=lambda item: float(item.get("score", float("-inf"))), reverse=True)
    for rank, item in enumerate(rows, start=1):
        item["rank"] = rank


def run_pipeline():
    """Save completed evaluations before optional structure exports."""
    config = None
    status = "failed"
    try:
        config = main.prepare_run_config(main.load_config(PROJECT_ROOT / "config.yaml"))
        RUN_STATE["total_iterations"] = int(config["mcts"]["iterations"])
        RUN_STATE["top_k"] = int(config["mcts"]["top_k"])
        RUN_STATE["run_id"] = config["runtime"]["run_id"]
        set_stage("search", "正在执行 MCTS 搜索")
        candidates, feedback_records = main.run_search(config, progress_callback=set_stage)
        main.save_outputs(candidates, feedback_records, [], [], [], config)
        builds, systems, inputs = [], [], []
        if config["output"].get("export_final_structures", True):
            set_stage("export_build", "正在导出 Top K 聚合物结构")
            builds = main.build_candidates(candidates, config)
            set_stage("export_packmol", "正在导出 Top K 体系")
            systems = main.prepare_systems(builds, config)
            set_stage("export_lammps", "正在导出 Top K 输入")
            inputs = main.write_thermal_inputs(systems, config)
        set_stage("save", "正在保存结果")
        main.save_outputs(candidates, feedback_records, builds, systems, inputs, config)
        status = "done"
        set_stage("done", "运行完成")
    except PipelineCancelled:
        status = "cancelled"
        RUN_STATE["error"] = ""
        set_stage("stopped", "任务已由用户终止，已完成的记录已保存")
    except Exception as exc:
        RUN_STATE["error"] = str(exc)
        try:
            set_stage("failed", f"运行失败：{exc}")
        except PipelineCancelled:
            status = "cancelled"
            set_stage("stopped", "任务已终止")
    finally:
        try:
            if config:
                main.finish_run(config, status, RUN_STATE["error"])
        finally:
            with RUN_LOCK:
                RUN_STATE["running"] = False
                RUN_STATE["finished_at"] = time.strftime("%Y-%m-%d %H:%M:%S")


def collect_counts():
    results = load_results()
    return {key: len(results[key]) for key in ("candidates", "build_results", "thermal_inputs", "low_k_database")}


def result_paths():
    config = main.load_config(PROJECT_ROOT / "config.yaml")
    root = Path(config["output"]["directory"])
    if not root.is_absolute():
        root = PROJECT_ROOT / root
    directory = root
    pointer = root / "latest_run.json"
    if pointer.exists():
        manifest = json.loads(pointer.read_text(encoding="utf-8"))
        candidate = Path(manifest["directory"]).resolve()
        if candidate.is_relative_to((root / "runs").resolve()):
            directory = candidate
    return {key: directory / config["output"][key + "_csv"]
            for key in ("candidates", "build_results", "system_results", "thermal_inputs", "low_k_database", "feedback_records")}


def load_results():
    """Return the latest run using configured directories and CSV names."""
    paths = result_paths()
    results = {key: read_csv(path) for key, path in paths.items()}
    results["thermal_results"] = results["feedback_records"]
    for row in results["thermal_results"] + results["low_k_database"]:
        row["measurement_status"] = "已通过校验" if row.get("run_id") and row.get("metadata_path") and row.get("success", "").lower() == "true" else "未通过新校验"
    return results


def clear_output_data(payload):
    """按前端请求清理输出数据。"""
    if RUN_STATE.get("running"):
        raise RuntimeError("pipeline is running, data cannot be cleared")

    mode = str(payload.get("mode", "")).strip()
    paths = result_paths()
    output_dir = paths["candidates"].parent
    output_dir.mkdir(parents=True, exist_ok=True)
    if mode == "selected_candidates":
        keys = set(str(item) for item in payload.get("candidate_keys", []))
        return clear_selected_candidates(paths["candidates"], keys)
    if mode == "all_candidates":
        return clear_csv_rows(paths["candidates"])
    if mode == "thermal_results":
        cleared = []
        for name in ["low_k_database.csv", "feedback_records.csv"]:
            result = clear_csv_rows(paths[Path(name).stem])
            cleared.extend(result["cleared_files"])
        return {"cleared_files": cleared, "message": "thermal feedback tables cleared"}
    if mode == "all_tables":
        cleared = []
        for name in [
            "candidates.csv",
            "build_results.csv",
            "system_results.csv",
            "thermal_inputs.csv",
            "low_k_database.csv",
            "feedback_records.csv",
        ]:
            result = clear_csv_rows(paths[Path(name).stem])
            cleared.extend(result["cleared_files"])
        return {"cleared_files": cleared, "message": "all output tables cleared"}
    raise ValueError(f"unknown clear mode: {mode}")


def clear_selected_candidates(filename, candidate_keys):
    """删除 candidates.csv 中选中的候选行。"""
    path = Path(filename)
    rows = read_csv(path)
    if not rows or not candidate_keys:
        return {"cleared_files": [], "removed": 0, "message": "no selected candidates"}

    backup_file(path)
    remaining = []
    removed = 0
    for row in rows:
        key = candidate_row_key(row)
        if key in candidate_keys:
            removed += 1
        else:
            remaining.append(row)
    write_csv_rows(path, remaining, rows[0].keys())
    return {
        "cleared_files": [str(path)],
        "removed": removed,
        "message": f"{removed} selected candidates removed",
    }


def clear_csv_rows(filename):
    """清空 CSV 数据行并保留表头。"""
    path = Path(filename)
    rows = read_csv(path)
    fieldnames = []
    if path.exists() and path.stat().st_size > 0:
        with path.open("r", encoding="utf-8", newline="") as file_obj:
            reader = csv.reader(file_obj)
            fieldnames = next(reader, [])
    if path.exists():
        backup_file(path)
    write_csv_rows(path, [], fieldnames)
    return {"cleared_files": [str(path)], "removed": len(rows)}


def backup_file(path):
    """清理前备份输出表。"""
    if not path.exists():
        return None
    archive_dir = path.parent / "archive"
    archive_dir.mkdir(parents=True, exist_ok=True)
    stamp = new_run_id()
    target = archive_dir / f"{path.stem}_{stamp}{path.suffix}"
    shutil.copy2(path, target)
    return target


def write_csv_rows(filename, rows, fieldnames):
    """写回 CSV。"""
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fieldnames or [])
    if not names and rows:
        names = list(rows[0].keys())
    with atomic_file(path, newline="") as file_obj:
        if not names:
            return
        writer = csv.DictWriter(file_obj, fieldnames=names)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in names})


def candidate_row_key(row):
    """生成前后端共用的候选行标识。"""
    rank = str(row.get("rank", "")).strip()
    sequence = str(row.get("sequence_text", "")).strip()
    return f"{rank}|{sequence}"


def read_csv(filename):
    """读取 CSV 为字典列表。"""
    path = Path(filename)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as file_obj:
        return list(csv.DictReader(file_obj))


def main_server():
    """启动本地服务。"""
    server = ThreadingHTTPServer((HOST, PORT), LammpsMctsHandler)
    print(f"LAMMPS_MCTS web server: http://localhost:{PORT}/frontend/")
    try:
        server.serve_forever()
    finally:
        RUN_CANCEL_EVENT.set()
        cancel_processes()
        server.server_close()


if __name__ == "__main__":
    main_server()
