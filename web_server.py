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

RUN_LOCK = threading.Lock()
RUN_STATE = {
    "running": False,
    "stage": "idle",
    "message": "等待运行",
    "started_at": "",
    "finished_at": "",
    "error": "",
    "log": [],
}

STAGE_LABELS = {
    "idle": "等待运行",
    "search": "MCTS 搜索",
    "build": "RDKit 建链",
    "packmol": "Packmol 初始体系",
    "write_lammps": "写出 LAMMPS 输入",
    "run_lammps": "运行 LAMMPS",
    "analyze": "热导率后处理",
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
        self.send_error(HTTPStatus.NOT_FOUND, "API endpoint not found")

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
    lammps_config = config.get("lammps", {})
    return {
        "iterations": int(mcts_config.get("iterations", 300)),
        "max_steps": int(mcts_config.get("max_steps", 6)),
        "top_k": int(mcts_config.get("top_k", 5)),
        "dp": int(polymer_config.get("degree_of_polymerization", 10)),
        "active_families": str(chemistry_config.get("active_families", "polyimide,polyurethane,phenolic")),
        "run_lammps": bool(lammps_config.get("run_enabled", False)),
        "feedback_mode": str(mcts_config.get("feedback_mode", "heuristic")),
    }


def update_config(payload):
    """更新常用配置并写回 config.yaml。"""
    config = main.load_config(PROJECT_ROOT / "config.yaml")
    mcts_config = config.setdefault("mcts", {})
    polymer_config = config.setdefault("polymer", {})
    chemistry_config = config.setdefault("chemistry", {})
    lammps_config = config.setdefault("lammps", {})

    if "iterations" in payload:
        mcts_config["iterations"] = int(payload["iterations"])
    if "max_steps" in payload:
        mcts_config["max_steps"] = int(payload["max_steps"])
    if "top_k" in payload:
        mcts_config["top_k"] = int(payload["top_k"])
    if "dp" in payload:
        polymer_config["degree_of_polymerization"] = int(payload["dp"])
    if "active_families" in payload:
        chemistry_config["active_families"] = str(payload["active_families"])
    if "run_lammps" in payload:
        lammps_config["run_enabled"] = bool(payload["run_lammps"])
    if "feedback_mode" in payload:
        mcts_config["feedback_mode"] = str(payload["feedback_mode"])

    write_config(config, PROJECT_ROOT / "config.yaml")
    return config


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
            {"key": "polyimide,polyurethane,phenolic", "label": "全部内置"},
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
        }
    )


def set_stage(stage, message):
    """更新当前阶段。"""
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


def run_pipeline():
    """按阶段执行主流程。"""
    old_cwd = Path.cwd()
    os.chdir(PROJECT_ROOT)
    try:
        config = main.load_config(PROJECT_ROOT / "config.yaml")
        set_stage("search", "正在执行 MCTS 搜索")
        candidates, feedback_records = main.run_search(config)

        set_stage("build", "正在使用 RDKit 构建聚合物")
        build_results = main.build_candidates(candidates, config)

        set_stage("packmol", "正在准备 Packmol 初始体系")
        system_results = main.prepare_systems(build_results, config)

        set_stage("write_lammps", "正在写出 LAMMPS 热导率输入")
        thermal_input_results = main.write_thermal_inputs(system_results, config)

        set_stage("run_lammps", "正在按配置运行 LAMMPS")
        lammps_run_results = main.run_lammps_jobs(thermal_input_results, config)

        set_stage("analyze", "正在解析热导率结果")
        thermal_analysis_results = main.analyze_thermal_results(thermal_input_results, config)

        set_stage("save", "正在保存 CSV 和说明文件")
        main.save_outputs(
            candidates,
            feedback_records,
            build_results,
            system_results,
            thermal_input_results,
            lammps_run_results,
            thermal_analysis_results,
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
        "thermal_results": len(read_csv(PROJECT_ROOT / "outputs" / "thermal_results.csv")),
    }


def load_results():
    """读取数据库表格。"""
    output_dir = PROJECT_ROOT / "outputs"
    return {
        "candidates": read_csv(output_dir / "candidates.csv"),
        "build_results": read_csv(output_dir / "build_results.csv"),
        "thermal_inputs": read_csv(output_dir / "thermal_inputs.csv"),
        "thermal_results": read_csv(output_dir / "thermal_results.csv"),
        "low_k_database": read_csv(output_dir / "low_k_database.csv"),
    }


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

