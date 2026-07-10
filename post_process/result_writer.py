"""候选序列和构建结果的 CSV 写出工具。"""

import csv
from pathlib import Path

from .feature_extractor import extract_sequence_features


def write_candidates_csv(candidates, filename):
    """保存 MCTS 候选序列及其统计特征。"""
    rows = []
    for index, candidate in enumerate(candidates, start=1):
        sequence = candidate["sequence"]
        features = extract_sequence_features(sequence)
        row = {
            "rank": index,
            "sequence": " ".join(sequence),
            "score": candidate.get("score", ""),
            "visits": candidate.get("visits", ""),
            "average_reward": candidate.get("average_reward", ""),
        }
        row.update(features)
        rows.append(row)

    _write_rows(rows, filename)


def write_build_results_csv(results, filename):
    """保存聚合物生成结果，便于后续接热导率计算。"""
    rows = []
    for index, result in enumerate(results, start=1):
        features = extract_sequence_features(result.sequence)
        row = {
            "index": index,
            "sequence": " ".join(result.sequence),
            "dp": result.dp,
            "success": result.success,
            "smiles": result.smiles,
            "pdb_path": result.pdb_path,
            "message": result.message,
        }
        row.update(features)
        rows.append(row)

    _write_rows(rows, filename)


def write_system_results_csv(results, filename):
    """保存 LAMMPS data 和 Packmol 初始体系的准备结果。"""
    rows = []
    for index, result in enumerate(results, start=1):
        rows.append(
            {
                "index": index,
                "template_pdb": result.template_pdb,
                "single_data_path": result.single_data_path,
                "packmol_input_path": result.packmol_input_path,
                "packed_pdb_path": result.packed_pdb_path,
                "system_data_path": result.system_data_path,
                "molecule_count": result.molecule_count,
                "success": result.success,
                "message": result.message,
            }
        )

    _write_rows(rows, filename)


def write_thermal_input_results_csv(results, filename):
    """保存快速热导率 LAMMPS 输入脚本的生成结果。"""
    rows = []
    for index, result in enumerate(results, start=1):
        rows.append(
            {
                "index": index,
                "data_path": result.data_path,
                "input_path": result.input_path,
                "flux_output": result.flux_output,
                "temp_output": result.temp_output,
                "dump_output": result.dump_output,
                "success": result.success,
                "message": result.message,
            }
        )

    _write_rows(rows, filename)


def write_lammps_run_results_csv(results, filename):
    """保存 LAMMPS 执行结果。"""
    rows = []
    for index, result in enumerate(results, start=1):
        rows.append(
            {
                "index": index,
                "input_path": result.input_path,
                "log_path": result.log_path,
                "success": result.success,
                "returncode": result.returncode,
                "message": result.message,
            }
        )
    _write_rows(rows, filename)


def write_thermal_analysis_results_csv(results, filename):
    """保存热导率后处理结果。"""
    rows = []
    sorted_results = sorted(
        results,
        key=lambda item: (not item.success, item.conductivity_w_mk),
    )
    for index, result in enumerate(sorted_results, start=1):
        rows.append(
            {
                "low_k_rank": index,
                "input_path": result.input_path,
                "flux_path": result.flux_path,
                "temp_path": result.temp_path,
                "heat_flux": result.heat_flux,
                "temperature_gradient": result.temperature_gradient,
                "conductivity_w_mk": result.conductivity_w_mk,
                "success": result.success,
                "message": result.message,
            }
        )
    _write_rows(rows, filename)


def write_mcts_feedback_records_csv(records, filename):
    """保存 MCTS 奖励回传记录。"""
    rows = []
    sorted_records = sorted(
        records,
        key=lambda item: (not item.success, item.conductivity_w_mk),
    )
    for index, record in enumerate(sorted_records, start=1):
        rows.append(
            {
                "low_k_rank": index,
                "sequence": " ".join(record.sequence),
                "reward": record.reward,
                "conductivity_w_mk": record.conductivity_w_mk,
                "success": record.success,
                "message": record.message,
                "input_path": record.input_path,
            }
        )
    _write_rows(rows, filename)


def _write_rows(rows, filename):
    """根据字典列表写 CSV，空结果也会生成一个文件。"""
    output_file = Path(filename)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        output_file.write_text("", encoding="utf-8")
        return

    with output_file.open("w", newline="", encoding="utf-8") as file_obj:
        writer = csv.DictWriter(file_obj, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
