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
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

import main


PROJECT_ROOT = Path(__file__).resolve().parent
HOST = "127.0.0.1"
PORT = 8010

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
            "production_steps": 30000,
            "sample_nevery": 10,
            "correlation_samples": 300,
            "thermo_every": 500,
            "dump_every": 5000,
        },
        "lammps": {"timeout_seconds": 1800},
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
            "production_steps": 100000,
            "sample_nevery": 10,
            "correlation_samples": 1000,
            "thermo_every": 1000,
            "dump_every": 10000,
        },
        "lammps": {"timeout_seconds": 7200},
    },
}

RUN_LOCK = threading.Lock()
RUN_STATE = {
    "running": False,
    "stage": "idle",
    "message": "等待运行",
    "started_at": "",
    "finished_at": "",
    "error": "",
    "log": [],
    "iteration": 0,
    "total_iterations": 0,
    "evaluation_count": 0,
    "top_k": 5,
    "live_candidates": [],
}

STAGE_LABELS = {
    "idle": "等待运行",
    "search": "MCTS 搜索",
    "mcts_lammps": "MCTS 热导反馈",
    "build": "RDKit 建链",
    "packmol": "Packmol 初始体系",
    "write_lammps": "写出 LAMMPS 输入",
    "compress": "高温熔融与体系压缩",
    "run_lammps": "Green-Kubo 计算",
    "export_build": "导出 Top K 结构",
    "export_packmol": "导出 Top K 体系",
    "export_lammps": "导出 Top K 输入",
    "save": "保存结果",
    "done": "完成",
    "failed": "失败",
}


class LammpsMctsHandler(SimpleHTTPRequestHandler):
    """处理静态页面和 API 请求。"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(PROJECT_ROOT), **kwargs)

    def end_headers(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(HTTPStatus.NO_CONTENT)
        self.end_headers()

    def do_GET(self):
        path = urlparse(self.path).path
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
            self._send_json({**RUN_STATE, "counts": collect_counts()})
            return
        if path == "/api/results":
            self._send_json(load_results())
            return
        super().do_GET()

    def do_POST(self):
        path = urlparse(self.path).path
        try:
            payload = self._read_json_body()
            if path == "/api/config":
                config = update_config(payload)
                self._send_json({"ok": True, "config": public_config_view(config)})
                return
            if path == "/api/custom-graph":
                saved = save_custom_graph(payload)
                self._send_json({"ok": True, "path": str(saved)})
                return
            if path == "/api/run":
                if payload:
                    update_config(payload)
                started = start_background_run()
                self._send_json({"ok": started, "status": RUN_STATE})
                return
            if path == "/api/clear-data":
                result = clear_output_data(payload)
                self._send_json({"ok": True, **result, "counts": collect_counts()})
                return
            self.send_error(HTTPStatus.NOT_FOUND, "API endpoint not found")
        except Exception as exc:
            self._send_json({"ok": False, "error": str(exc)})

    def _read_json_body(self):
        length = int(self.headers.get("Content-Length", "0") or 0)
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
    return {
        "iterations": int(mcts_config.get("iterations", 300)),
        "max_steps": int(mcts_config.get("max_steps", 5)),
        "top_k": int(mcts_config.get("top_k", 5)),
        "exploration_weight": float(mcts_config.get("exploration_weight", 1.41421356237)),
        "target_length_angstrom": float(
            polymer_config.get("target_chain_length_angstrom", 120.0)
        ),
        "active_families": str(chemistry_config.get("active_families", "polyimide")),
        "run_lammps": bool(mcts_config.get("feedback_run_lammps", False)),
        "feedback_mode": str(mcts_config.get("feedback_mode", "heuristic")),
        "feedback_run_lammps": bool(mcts_config.get("feedback_run_lammps", False)),
        "feedback_every_iteration": bool(mcts_config.get("feedback_every_iteration", True)),
        "feedback_max_evaluations": int(mcts_config.get("feedback_max_evaluations", 0) or 0),
        "gk_profile": str(tc_config.get("profile", "quick")),
        "density_equilibration": bool(tc_config.get("density_equilibration", True)),
        "precompression_density": float(tc_config.get("precompression_density", 0.70)),
        "density_plateau_percent": float(
            tc_config.get("density_plateau_tolerance", 0.05)
        ) * 100.0,
        "density_max_blocks": int(tc_config.get("density_max_blocks", 10)),
        "timeout_seconds": int(lammps_config.get("timeout_seconds", 1800) or 0),
    }


def update_config(payload):
    """更新常用配置并写回 config.yaml。"""
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
    if "gk_profile" in payload:
        apply_thermal_profile(config, str(payload["gk_profile"]))
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
    if "timeout_seconds" in payload:
        lammps_config["timeout_seconds"] = max(0, int(payload["timeout_seconds"]))

    if mcts_config.get("feedback_mode") == "thermal":
        mcts_config["feedback_every_iteration"] = True
        mcts_config["feedback_max_evaluations"] = max(
            int(mcts_config.get("feedback_max_evaluations", 0) or 0),
            int(mcts_config.get("iterations", 0) or 0),
        )
        if "run_lammps" in payload:
            mcts_config["feedback_run_lammps"] = bool(payload["run_lammps"])

    write_config(config, PROJECT_ROOT / "config.yaml")
    return config


def apply_thermal_profile(config, profile_name):
    """应用快速筛选或标准复核参数。"""
    name = profile_name if profile_name in THERMAL_PROFILES else "quick"
    profile = THERMAL_PROFILES[name]
    for section, values in profile.items():
        config.setdefault(section, {}).update(values)
    config.setdefault("thermal_conductivity", {})["profile"] = name


def write_config(config, filename):
    """写出项目当前使用的二级 YAML 配置。"""
    lines = []
    for section, values in config.items():
        if not isinstance(values, dict):
            continue
        lines.append(f"{section}:")
        for key, value in values.items():
            if isinstance(value, (list, tuple)):
                value = ",".join(str(item) for item in value)
            if isinstance(value, dict):
                continue
            if isinstance(value, bool):
                value = "true" if value else "false"
            elif isinstance(value, str) and not value:
                value = '""'
            lines.append(f"  {key}: {value}")
        lines.append("")
    Path(filename).write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


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
    fragments = payload.get("fragments", [])
    transitions = payload.get("transitions", {})
    if not isinstance(fragments, list) or not isinstance(transitions, dict):
        raise ValueError("custom graph payload must contain fragments and transitions")
    data = {"fragments": fragments, "transitions": transitions}
    path = PROJECT_ROOT / "custom_fragments.json"
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


def start_background_run():
    """启动后台运行线程。"""
    with RUN_LOCK:
        if RUN_STATE["running"]:
            return False
        reset_run_state()
        thread = threading.Thread(target=run_pipeline, daemon=True)
        thread.start()
        return True


def reset_run_state():
    """重置运行状态。"""
    RUN_STATE.update(
        {
            "running": True,
            "stage": "search",
            "message": "准备开始 MCTS 搜索",
            "started_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "finished_at": "",
            "error": "",
            "log": [],
            "iteration": 0,
            "total_iterations": 0,
            "evaluation_count": 0,
            "top_k": 5,
            "live_candidates": [],
        }
    )


def set_stage(stage, message, details=None):
    """更新当前阶段，并接收 MCTS 的轮次和实时候选事件。"""
    details = dict(details or {})
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
    """按阶段执行主流程。"""
    old_cwd = Path.cwd()
    os.chdir(PROJECT_ROOT)
    try:
        config = main.load_config(PROJECT_ROOT / "config.yaml")
        RUN_STATE["total_iterations"] = int(config["mcts"]["iterations"])
        RUN_STATE["top_k"] = int(config["mcts"]["top_k"])
        set_stage("search", "正在执行 MCTS 搜索")
        candidates, feedback_records = main.run_search(config, progress_callback=set_stage)

        set_stage("export_build", "正在为最终 Top K 候选导出聚合物结构")
        build_results = main.build_candidates(candidates, config)

        set_stage("export_packmol", "正在为最终 Top K 候选导出 Packmol 体系")
        system_results = main.prepare_systems(build_results, config)

        set_stage("export_lammps", "正在为最终 Top K 候选导出 LAMMPS 输入")
        thermal_input_results = main.write_thermal_inputs(system_results, config)

        set_stage("save", "正在保存 CSV 和说明文件")
        main.save_outputs(
            candidates,
            feedback_records,
            build_results,
            system_results,
            thermal_input_results,
            config,
        )
        set_stage("done", "运行完成")
    except Exception as exc:
        RUN_STATE["error"] = str(exc)
        set_stage("failed", f"运行失败：{exc}")
    finally:
        RUN_STATE["running"] = False
        RUN_STATE["finished_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
        os.chdir(old_cwd)


def collect_counts():
    """统计前端状态卡片需要的数量。"""
    return {
        "candidates": len(read_csv(PROJECT_ROOT / "outputs" / "candidates.csv")),
        "build_results": len(read_csv(PROJECT_ROOT / "outputs" / "build_results.csv")),
        "thermal_inputs": len(read_csv(PROJECT_ROOT / "outputs" / "thermal_inputs.csv")),
        "low_k_database": len(read_csv(PROJECT_ROOT / "outputs" / "low_k_database.csv")),
    }


def load_results():
    """读取数据库表格。"""
    output_dir = PROJECT_ROOT / "outputs"
    return {
        "candidates": read_csv(output_dir / "candidates.csv"),
        "build_results": read_csv(output_dir / "build_results.csv"),
        "thermal_inputs": read_csv(output_dir / "thermal_inputs.csv"),
        "thermal_results": [],
        "low_k_database": read_csv(output_dir / "low_k_database.csv"),
    }


def clear_output_data(payload):
    """按前端请求清理输出数据。"""
    if RUN_STATE.get("running"):
        raise RuntimeError("pipeline is running, data cannot be cleared")

    mode = str(payload.get("mode", "")).strip()
    output_dir = PROJECT_ROOT / "outputs"
    output_dir.mkdir(exist_ok=True)
    if mode == "selected_candidates":
        keys = set(str(item) for item in payload.get("candidate_keys", []))
        return clear_selected_candidates(output_dir / "candidates.csv", keys)
    if mode == "all_candidates":
        return clear_csv_rows(output_dir / "candidates.csv")
    if mode == "thermal_results":
        cleared = []
        for name in ["low_k_database.csv", "feedback_records.csv"]:
            result = clear_csv_rows(output_dir / name)
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
            result = clear_csv_rows(output_dir / name)
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
        "cleared_files": [str(path.relative_to(PROJECT_ROOT))],
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
    return {"cleared_files": [str(path.relative_to(PROJECT_ROOT))], "removed": len(rows)}


def backup_file(path):
    """清理前备份输出表。"""
    if not path.exists():
        return None
    archive_dir = PROJECT_ROOT / "outputs" / "archive"
    archive_dir.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
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
    with path.open("w", encoding="utf-8", newline="") as file_obj:
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
    server.serve_forever()


if __name__ == "__main__":
    main_server()
