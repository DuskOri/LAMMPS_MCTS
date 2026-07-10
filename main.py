"""项目主入口：执行 MCTS 搜索并生成候选聚合物结构。"""

import json
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None

from generator import build_polymer_from_sequence, prepare_initial_system
from md_engine import run_lammps_input, write_fast_tc_input
from mcts import MCTSEngine, ThermalFeedbackEvaluator, configure_fragment_registry
from mcts.ai_fragment_api import generate_ai_fragment_file
from post_process import (
    analyze_thermal_outputs,
    write_build_results_csv,
    write_candidates_csv,
    write_lammps_run_results_csv,
    write_mcts_feedback_records_csv,
    write_system_results_csv,
    write_thermal_analysis_results_csv,
    write_thermal_input_results_csv,
)
from post_process.explain_model import summarize_candidates


PROJECT_ROOT = Path(__file__).resolve().parent


DEFAULT_CONFIG = {
    "mcts": {
        "iterations": 300,
        "max_steps": 6,
        "start_fragment": "Start",
        "top_k": 5,
        "exploration_weight": 1.41421356237,
        "rollout_limit": 32,
        "random_seed": 7,
        "feedback_mode": "heuristic",
        "feedback_max_evaluations": 10,
        "feedback_output_dir": "outputs/mcts_feedback",
        "feedback_run_lammps": False,
        "feedback_reward_offset": 0.05,
        "feedback_fallback_to_heuristic": True,
        "feedback_failure_reward": -1.0,
    },
    "chemistry": {
        "active_families": ["polyimide", "polyurethane", "phenolic"],
        "custom_fragment_file": "custom_fragments.json",
        "ai_enabled": False,
        "ai_api_url": "https://api.openai.com/v1/chat/completions",
        "ai_api_key_env": "OPENAI_API_KEY",
        "ai_model": "gpt-4.1-mini",
        "ai_output_file": "ai_fragments.json",
        "ai_refresh": False,
        "ai_temperature": 0.2,
        "ai_max_fragments": 8,
        "ai_timeout_seconds": 60,
        "ai_validate_rdkit": True,
        "ai_prompt": "Generate low thermal conductivity polymer fragments for MCTS screening.",
        "custom_fragments": [],
        "custom_transitions": {},
        "allow_custom_self_transitions": True,
    },
    "polymer": {
        "degree_of_polymerization": 10,
        "build_pdb": True,
    },
    "system": {
        "prepare_initial_system": True,
        "molecule_count": 10,
        "box_size": 60.0,
        "tolerance": 2.0,
        "packmol_executable": "packmol",
    },
    "thermal_conductivity": {
        "write_fast_input": True,
        "temperature": 300.0,
        "hot_temperature": 450.0,
        "cold_temperature": 150.0,
        "timestep": 0.5,
        "tdamp": 50.0,
        "cutoff": 12.0,
        "equil_steps": 5000,
        "gradient_steps": 20000,
        "dump_every": 5000,
        "analyze_outputs": True,
        "trim_fraction": 0.20,
    },
    "lammps": {
        "run_enabled": False,
        "executable": "lmp",
        "timeout_seconds": 120,
    },
    "output": {
        "directory": "outputs",
        "candidates_csv": "candidates.csv",
        "build_results_csv": "build_results.csv",
        "system_results_csv": "system_results.csv",
        "thermal_inputs_csv": "thermal_inputs.csv",
        "lammps_runs_csv": "lammps_runs.csv",
        "thermal_results_csv": "thermal_results.csv",
        "low_k_database_csv": "low_k_database.csv",
        "feedback_records_csv": "feedback_records.csv",
    },
}


def load_config(filename="config.yaml"):
    """读取配置文件；没有配置时使用默认参数。"""
    config_path = Path(filename)
    if not config_path.is_absolute():
        config_path = PROJECT_ROOT / config_path

    if not config_path.exists() or config_path.stat().st_size == 0:
        return DEFAULT_CONFIG

    if yaml is not None:
        with config_path.open("r", encoding="utf-8") as file_obj:
            user_config = yaml.safe_load(file_obj) or {}
    else:
        user_config = _load_simple_yaml(config_path)

    return _merge_config(DEFAULT_CONFIG, user_config)


def run_search(config):
    """运行 MCTS 搜索，返回按分数排序的候选序列。"""
    configure_chemistry(config)
    mcts_config = config["mcts"]
    reward_fn = None
    feedback_records = []

    if mcts_config.get("feedback_mode", "heuristic") == "thermal":
        evaluator = ThermalFeedbackEvaluator(
            config=config,
            output_dir=mcts_config["feedback_output_dir"],
            max_evaluations=mcts_config["feedback_max_evaluations"],
            run_lammps=mcts_config["feedback_run_lammps"],
            reward_offset=mcts_config["feedback_reward_offset"],
            fallback_to_heuristic=mcts_config["feedback_fallback_to_heuristic"],
            failure_reward=mcts_config["feedback_failure_reward"],
        )
        reward_fn = evaluator
    else:
        evaluator = None

    engine = MCTSEngine(
        max_steps=mcts_config["max_steps"],
        start_fragment=mcts_config["start_fragment"],
        exploration_weight=mcts_config["exploration_weight"],
        rollout_limit=mcts_config["rollout_limit"],
        reward_fn=reward_fn,
        random_seed=mcts_config.get("random_seed"),
    )
    candidates = engine.ranked_candidates(
        iterations=mcts_config["iterations"],
        top_k=mcts_config["top_k"],
    )
    if evaluator is not None:
        feedback_records = evaluator.records
    return candidates, feedback_records


def build_candidates(candidates, config):
    """对候选序列调用 RDKit 建模，LAMMPS 模拟暂不在这里执行。"""
    polymer_config = config["polymer"]
    output_dir = Path(config["output"]["directory"]) / "pdb"

    if not polymer_config.get("build_pdb", True):
        return []

    results = []
    for index, candidate in enumerate(candidates, start=1):
        result = build_polymer_from_sequence(
            candidate["sequence"],
            dp=polymer_config["degree_of_polymerization"],
            output_dir=output_dir,
            name=f"candidate_{index:03d}",
        )
        results.append(result)
    return results


def prepare_systems(build_results, config):
    """把 PDB 转成 LAMMPS data，并尝试用 Packmol 生成多链初始体系。"""
    system_config = config["system"]
    output_dir = Path(config["output"]["directory"]) / "systems"

    if not system_config.get("prepare_initial_system", True):
        return []

    results = []
    for result in build_results:
        if not result.success or not result.pdb_path:
            continue

        system_result = prepare_initial_system(
            template_pdb=result.pdb_path,
            output_dir=output_dir,
            molecule_count=system_config["molecule_count"],
            box_size=system_config["box_size"],
            tolerance=system_config["tolerance"],
            packmol_executable=system_config["packmol_executable"],
        )
        results.append(system_result)

    return results


def write_thermal_inputs(system_results, config):
    """为已经生成的体系 data 文件输出快速热导率 LAMMPS 输入脚本。"""
    tc_config = config["thermal_conductivity"]
    output_dir = Path(config["output"]["directory"]) / "lammps"

    if not tc_config.get("write_fast_input", True):
        return []

    params = {
        key: value
        for key, value in tc_config.items()
        if key != "write_fast_input"
    }

    results = []
    for index, system_result in enumerate(system_results, start=1):
        if not system_result.success or not system_result.system_data_path:
            continue

        stem = Path(system_result.system_data_path).stem.replace("_system", "")
        input_file = output_dir / f"{stem}_fast_tc.in"
        output_prefix = output_dir / stem
        result = write_fast_tc_input(
            data_file=system_result.system_data_path,
            input_file=input_file,
            output_prefix=output_prefix,
            params=params,
        )
        results.append(result)

    return results


def run_lammps_jobs(thermal_input_results, config):
    """按配置执行 LAMMPS 快速热导率脚本。"""
    lammps_config = config["lammps"]
    if not lammps_config.get("run_enabled", False):
        return []

    results = []
    for result in thermal_input_results:
        if not result.success:
            continue
        run_result = run_lammps_input(
            input_path=result.input_path,
            lammps_executable=lammps_config["executable"],
            timeout=lammps_config.get("timeout_seconds"),
        )
        results.append(run_result)
    return results


def analyze_thermal_results(thermal_input_results, config):
    """解析热流和温度剖面，估算快速热导率。"""
    tc_config = config["thermal_conductivity"]
    if not tc_config.get("analyze_outputs", True):
        return []

    results = []
    for result in thermal_input_results:
        if not result.success:
            continue
        analysis = analyze_thermal_outputs(
            input_path=result.input_path,
            flux_path=result.flux_output,
            temp_path=result.temp_output,
            trim_fraction=tc_config.get("trim_fraction", 0.20),
        )
        results.append(analysis)
    return results


def save_outputs(
    candidates,
    feedback_records,
    build_results,
    system_results,
    thermal_input_results,
    lammps_run_results,
    thermal_analysis_results,
    config,
):
    """保存候选序列、结构生成结果和解释信息。"""
    output_dir = Path(config["output"]["directory"])
    output_dir.mkdir(parents=True, exist_ok=True)

    write_candidates_csv(candidates, output_dir / config["output"]["candidates_csv"])
    write_mcts_feedback_records_csv(
        feedback_records,
        output_dir / config["output"]["feedback_records_csv"],
    )
    write_build_results_csv(
        build_results,
        output_dir / config["output"]["build_results_csv"],
    )
    write_system_results_csv(
        system_results,
        output_dir / config["output"]["system_results_csv"],
    )
    write_thermal_input_results_csv(
        thermal_input_results,
        output_dir / config["output"]["thermal_inputs_csv"],
    )
    write_lammps_run_results_csv(
        lammps_run_results,
        output_dir / config["output"]["lammps_runs_csv"],
    )
    write_thermal_analysis_results_csv(
        thermal_analysis_results,
        output_dir / config["output"]["thermal_results_csv"],
    )
    write_thermal_analysis_results_csv(
        thermal_analysis_results,
        output_dir / config["output"]["low_k_database_csv"],
    )
    _write_explain_text(summarize_candidates(candidates), output_dir / "candidate_notes.txt")


def _write_explain_text(explanations, filename):
    """把候选序列的简单解释写成文本，方便人工快速查看。"""
    lines = []
    for index, item in enumerate(explanations, start=1):
        lines.append(f"[{index}] {item['sequence']}")
        lines.append(f"score: {item['score']}")
        for note in item["notes"]:
            lines.append(f"- {note}")
        lines.append("")

    Path(filename).write_text("\n".join(lines), encoding="utf-8")


def configure_chemistry(config):
    """根据配置刷新片段库。

    自定义化合物片段可以直接写入配置，也可以放在 custom_fragment_file 指向的
    JSON 文件中。JSON 入口只依赖 Python 标准库，适合没有 PyYAML 的运行环境。
    """
    chemistry_config = dict(config.get("chemistry", {}))
    ai_enabled = _as_config_bool(chemistry_config.get("ai_enabled", False))
    if ai_enabled:
        ai_result = generate_ai_fragment_file(chemistry_config, project_root=PROJECT_ROOT)
        if ai_result.success:
            chemistry_config["ai_generated_file"] = ai_result.output_path
        else:
            print(f"AI fragment generation skipped: {ai_result.message}")

    ai_file = chemistry_config.get("ai_generated_file") or chemistry_config.get("ai_output_file")
    if ai_enabled and ai_file:
        loaded = _load_custom_fragment_file(ai_file)
        if loaded:
            file_fragments = loaded.get("custom_fragments", loaded.get("fragments", []))
            file_transitions = loaded.get("custom_transitions", loaded.get("transitions", {}))
            chemistry_config["custom_fragments"] = list(
                chemistry_config.get("custom_fragments") or []
            ) + list(file_fragments or [])
            chemistry_config["custom_transitions"] = {
                **(file_transitions or {}),
                **(chemistry_config.get("custom_transitions") or {}),
            }

    custom_file = chemistry_config.get("custom_fragment_file")
    if custom_file:
        loaded = _load_custom_fragment_file(custom_file)
        if loaded:
            file_fragments = loaded.get("custom_fragments", loaded.get("fragments", []))
            file_transitions = loaded.get("custom_transitions", loaded.get("transitions", {}))
            chemistry_config["custom_fragments"] = list(
                chemistry_config.get("custom_fragments") or []
            ) + list(file_fragments or [])
            chemistry_config["custom_transitions"] = {
                **(file_transitions or {}),
                **(chemistry_config.get("custom_transitions") or {}),
            }
    configure_fragment_registry(chemistry_config)


def _as_config_bool(value):
    """把配置中的布尔值统一转换为 Python bool。"""
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


def _load_custom_fragment_file(filename):
    """读取自定义片段 JSON 文件；没有该文件时使用内置片段库。"""
    path = Path(filename)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as file_obj:
        return json.load(file_obj)


def _merge_config(default_config, user_config):
    """递归合并配置，用户只需要写想改的字段。"""
    merged = {}
    for key, value in default_config.items():
        if isinstance(value, dict):
            merged[key] = _merge_config(value, user_config.get(key, {}))
        else:
            merged[key] = user_config.get(key, value)
    return merged


def _load_simple_yaml(filename):
    """读取本项目当前使用的二级 YAML 配置。"""
    data = {}
    current_section = None

    for raw_line in Path(filename).read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line:
            continue

        if not raw_line.startswith(" ") and line.endswith(":"):
            current_section = line[:-1].strip()
            data[current_section] = {}
            continue

        if current_section is None or ":" not in line:
            continue

        key, value = line.strip().split(":", 1)
        data[current_section][key.strip()] = _parse_config_value(value.strip())

    return data


def _parse_config_value(value):
    """把配置文件里的字符串转换为常用 Python 类型。"""
    if value.lower() in ("true", "false"):
        return value.lower() == "true"

    try:
        if "." in value:
            return float(value)
        return int(value)
    except ValueError:
        return value.strip("\"'")


def main():
    """执行完整的非 LAMMPS 流程。"""
    config = load_config()
    candidates, feedback_records = run_search(config)
    build_results = build_candidates(candidates, config)
    system_results = prepare_systems(build_results, config)
    thermal_input_results = write_thermal_inputs(system_results, config)
    lammps_run_results = run_lammps_jobs(thermal_input_results, config)
    thermal_analysis_results = analyze_thermal_results(thermal_input_results, config)
    save_outputs(
        candidates,
        feedback_records,
        build_results,
        system_results,
        thermal_input_results,
        lammps_run_results,
        thermal_analysis_results,
        config,
    )

    print(f"候选序列数量: {len(candidates)}")
    print(f"MCTS反馈记录数量: {len(feedback_records)}")
    print(f"结构生成数量: {len(build_results)}")
    print(f"初始体系数量: {len(system_results)}")
    print(f"热导率输入脚本数量: {len(thermal_input_results)}")
    print(f"LAMMPS运行数量: {len(lammps_run_results)}")
    print(f"热导率结果数量: {len(thermal_analysis_results)}")
    print(f"结果目录: {config['output']['directory']}")


if __name__ == "__main__":
    main()
