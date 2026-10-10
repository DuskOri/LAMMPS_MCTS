"""把热导率计算结果接回 MCTS 奖励函数。"""

from dataclasses import dataclass, field
from pathlib import Path
import shutil
from run_storage import new_run_id, file_hash, write_json
from md_engine.process_control import PipelineCancelled, raise_if_cancelled

from generator import build_polymer_from_sequence, prepare_initial_system
from md_engine import (
    run_lammps_input,
    write_adaptive_nemd_input,
    write_compression_input,
    write_direct_nemd_input,
)
from mcts.mcts_engine import default_reward
from post_process import analyze_direct_nemd_outputs
from post_process.result_writer import write_mcts_feedback_records_csv
from md_engine.thermal_conductivity_writer import _normalise_params


@dataclass
class ThermalFeedbackRecord:
    """记录一次序列热导率评估过程。"""

    sequence: list
    reward: float
    conductivity_w_mk: float
    success: bool
    message: str
    input_path: str = ""
    run_id: str = ""
    evaluation_id: str = ""
    metadata_path: str = ""
    relative_error: float = 0.0
    energy_balance_error: float = 0.0
    flux_agreement_error: float = 0.0


@dataclass
class ThermalFeedbackEvaluator:
    """把聚合物序列转换为热导率奖励值。"""

    config: dict
    output_dir: str = "data"
    max_evaluations: object = 10
    run_lammps: bool = False
    reward_scale: float = 0.20
    fallback_to_heuristic: bool = False
    failure_reward: float = -1.0
    progress_callback: object = None
    cache: dict = field(default_factory=dict)
    records: list = field(default_factory=list)
    run_id: str = ""

    def __post_init__(self):
        self.run_id = self.config.get("runtime", {}).get("run_id") or new_run_id()
        if not self.config.get("runtime", {}).get("run_id"):
            self.output_dir = str(Path(self.output_dir) / "runs" / self.run_id)
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def __call__(self, sequence):
        """作为 reward_fn 被 MCTS 调用。"""
        key = tuple(sequence)
        raise_if_cancelled()
        if key in self.cache:
            return self.cache[key]

        self._stages = {"compression": "not_run", "nemd": "not_run", "measurement": "not_run"}
        try:
            record = self._failed(sequence, "thermal evaluation limit reached") if self._limit_reached() else self._evaluate(sequence)
        except PipelineCancelled:
            write_json(Path(self.output_dir) / f"eval_{len(self.records) + 1:04d}" / "metadata.json", {"run_id": self.run_id, "sequence": list(sequence), "status": "cancelled", "config": self.config})
            raise
        except (OSError, ValueError, RuntimeError) as error:
            record = self._failed(sequence, f"candidate evaluation failed: {error}")
        record.run_id = self.run_id
        record.evaluation_id = f"eval_{len(self.records) + 1:04d}"
        base_dir = Path(self.output_dir) / record.evaluation_id
        record.metadata_path = str(base_dir / "metadata.json")
        artifacts = []
        if base_dir.exists():
            for path in sorted(base_dir.rglob("*")):
                if path.is_file() and path.name != "metadata.json":
                    item = {"path": str(path.relative_to(base_dir)), "bytes": path.stat().st_size}
                    if path.suffix in (".in", ".data", ".mol", ".dat"):
                        item["sha256"] = file_hash(path)
                    artifacts.append(item)
        write_json(record.metadata_path, {"schema_version": 2, "record": record, "config": self.config, "stages": self._stages, "artifacts": artifacts})
        self.cache[key] = record.reward
        self.records.append(record)
        write_json(Path(self.output_dir) / "records.json", self.records)
        write_mcts_feedback_records_csv(self.records, Path(self.config["output"]["directory"]) / self.config["output"]["feedback_records_csv"])
        self._notify(
            "evaluation",
            f"第 {len(self.cache)} 个候选已完成热导评价，reward 已回传",
            {
                "evaluation": len(self.cache),
                "candidate": {
                    "sequence": list(record.sequence),
                    "score": record.reward,
                    "reward": record.reward,
                    "conductivity_w_mk": record.conductivity_w_mk,
                    "success": record.success,
                    "message": record.message,
                    "input_path": record.input_path,
                },
            },
        )
        return record.reward

    def _limit_reached(self):
        """判断热导反馈次数是否已经用完。"""
        if self.max_evaluations is None:
            return False
        return len(self.cache) >= int(self.max_evaluations)

    def _evaluate(self, sequence):
        """执行一次从序列到热导率的快速评估。"""
        index = len(self.records) + 1
        base_dir = Path(self.output_dir) / f"eval_{index:04d}"
        base_dir.mkdir(parents=True, exist_ok=False)
        write_json(base_dir / "metadata.json", {"run_id": self.run_id, "sequence": list(sequence), "status": "running", "config": self.config})
        pdb_dir = base_dir / "polymer"
        system_dir = base_dir / "packmol"

        polymer_config = self.config["polymer"]
        system_config = self.config["system"]
        tc_config = self.config["thermal_conductivity"]
        lammps_config = self.config["lammps"]
        self._notify(
            "mcts_lammps",
            f"MCTS thermal feedback {index}: building polymer structure",
        )

        build_result = build_polymer_from_sequence(
            sequence,
            target_length_angstrom=polymer_config["target_chain_length_angstrom"],
            min_dp=polymer_config["min_degree_of_polymerization"],
            max_dp=polymer_config["max_degree_of_polymerization"],
            output_dir=pdb_dir,
            name="candidate",
        )
        if not build_result.success:
            return self._failed(sequence, build_result.message)

        self._notify(
            "mcts_lammps",
            f"MCTS thermal feedback {index}: preparing Packmol system",
        )
        system_result = prepare_initial_system(
            template_pdb=build_result.pdb_path,
            output_dir=system_dir,
            molecule_count=system_config["molecule_count"],
            box_size=system_config["box_size"],
            tolerance=system_config["tolerance"],
            packmol_executable=system_config["packmol_executable"],
            template_mol=build_result.mol_path,
            force_field=system_config.get("force_field", "uff_screening"),
            timeout=system_config.get("packmol_timeout_seconds", 600),
        )
        if not system_result.success:
            return self._failed(sequence, system_result.message)

        initial_data = base_dir / "initial.data"
        compressed_data = base_dir / "compressed.data"
        shutil.copy2(system_result.system_data_path, initial_data)
        thermal_params = {
            key: value
            for key, value in tc_config.items()
            if key not in ("write_fast_input", "analyze_outputs", "profile")
        }

        self._notify(
            "mcts_lammps",
            f"MCTS thermal feedback {index}: writing LAMMPS input",
        )

        if self.run_lammps:
            compression_result = write_compression_input(
                data_file=initial_data,
                input_file=base_dir / "compression.in",
                compressed_data_file=compressed_data,
                output_prefix=base_dir / "compression",
                params=thermal_params,
            )
            if not compression_result.success:
                return self._failed(sequence, compression_result.message)

            self._notify(
                "compress",
                f"MCTS thermal feedback {index}: melting and compressing the system",
                {
                    "telemetry_evaluation": index,
                    "series_files": {
                        "density": str(base_dir / "compression_density.dat"),
                        "heat_flux": "",
                    },
                },
            )

            def report_compression_line(line):
                text = str(line).strip()
                self._notify("compress", text, {"console_line": text})

            compression_run = run_lammps_input(
                input_path=compression_result.input_path,
                lammps_executable=lammps_config["executable"],
                log_path=base_dir / "compression.log",
                timeout=lammps_config.get(
                    "compression_timeout_seconds",
                    lammps_config.get("timeout_seconds"),
                ) or None,
                line_callback=report_compression_line,
            )
            self._stages["compression"] = {"returncode": compression_run.returncode, "success": compression_run.success}
            self._stages["density"] = "passed" if compression_run.success and tc_config.get("density_equilibration", True) else "not_validated"
            if not compression_run.success:
                return self._failed(
                    sequence, compression_run.message, compression_result.input_path
                )
            if not compressed_data.exists():
                return self._failed(
                    sequence,
                    "compression finished without writing compressed.data",
                    compression_result.input_path,
                )

            input_result = write_adaptive_nemd_input(
                data_file=compressed_data,
                input_file=base_dir / "nemd.in",
                output_prefix=base_dir / "nemd",
                params=thermal_params,
            )
            if not input_result.success:
                return self._failed(sequence, input_result.message)

            self._notify(
                "run_lammps",
                f"MCTS thermal feedback {index}: running adaptive NEMD",
                {
                    "telemetry_evaluation": index,
                    "series_files": {
                        "density": str(base_dir / "compression_density.dat"),
                        "heat_flux": str(base_dir / "nemd_heat_flux.dat"),
                    },
                },
            )

            def report_nemd_line(line):
                text = str(line).strip()
                message = text
                if text.startswith("Heat flux convergence check"):
                    message = f"NEMD production sampling: {text}"
                elif text.startswith("Heat flux converged"):
                    message = "NEMD heat flux converged; preparing result"
                self._notify("run_lammps", message, {"console_line": text})

            nemd_run = run_lammps_input(
                input_path=input_result.input_path,
                lammps_executable=lammps_config["executable"],
                log_path=base_dir / "nemd.log",
                timeout=lammps_config.get(
                    "nemd_timeout_seconds",
                    lammps_config.get("timeout_seconds"),
                ) or None,
                line_callback=report_nemd_line,
            )
            self._stages["nemd"] = {"returncode": nemd_run.returncode, "success": nemd_run.success}
            self._stages["heat_flux_convergence"] = "passed" if nemd_run.success else "not_validated"
            if not nemd_run.success:
                return self._failed(sequence, nemd_run.message, input_result.input_path)
        else:
            input_result = write_direct_nemd_input(
                data_file=initial_data,
                input_file=base_dir / "nemd.in",
                output_prefix=base_dir / "nemd",
                params=thermal_params,
            )
            if not input_result.success:
                return self._failed(sequence, input_result.message)
            return self._failed(sequence, "input prepared; LAMMPS was not executed", input_result.input_path)

        self._notify(
            "mcts_lammps",
            f"MCTS thermal feedback {index}: analyzing thermal output",
        )
        analysis = analyze_direct_nemd_outputs(
            input_path=input_result.input_path,
            temperature_profile_path=input_result.correlation_output,
            heat_flux_path=input_result.conductivity_output,
            tail_fraction=tc_config.get("trim_fraction", 0.50),
            min_gradient_r2=tc_config.get("nemd_min_gradient_r2", 0.70),
            min_temperature_span=tc_config.get("nemd_min_temperature_span", 5.0),
            bath_energy_path=base_dir / "nemd_bath_energy.dat",
            timestep_fs=tc_config.get("timestep", 0.5),
            block_steps=_normalise_params(thermal_params)["nemd_production_block_steps"],
            max_relative_error=tc_config.get("nemd_max_relative_error", 0.20),
            energy_balance_tolerance=tc_config.get("nemd_energy_balance_tolerance", 0.20),
            flux_agreement_tolerance=tc_config.get("nemd_flux_agreement_tolerance", 0.30),
        )
        if not analysis.success:
            self._stages["measurement"] = "failed"
            record = self._failed(sequence, analysis.message, input_result.input_path)
            record.relative_error = analysis.relative_error
            record.energy_balance_error = analysis.energy_balance_error
            record.flux_agreement_error = analysis.flux_agreement_error
            return record

        analysis_message = (
            f"ok; gradient={analysis.gradient_k_per_a:.6f} K/A; "
            f"R2={analysis.gradient_r2:.4f}"
        )
        self._stages["measurement"] = "passed"
        return ThermalFeedbackRecord(
            sequence=list(sequence),
            reward=self._conductivity_to_reward(analysis.conductivity_w_mk),
            conductivity_w_mk=analysis.conductivity_w_mk,
            success=True,
            message=analysis_message,
            input_path=input_result.input_path,
            relative_error=analysis.relative_error,
            energy_balance_error=analysis.energy_balance_error,
            flux_agreement_error=analysis.flux_agreement_error,
        )

    def _conductivity_to_reward(self, conductivity_w_mk):
        """将低热导目标映射到 0 到 1，使奖励项和 UCB 探索项尺度接近。"""
        conductivity = max(conductivity_w_mk, 0.0)
        scale = max(float(self.reward_scale), 1.0e-12)
        return 1.0 / (1.0 + conductivity / scale)

    def _failed(self, sequence, message, input_path=""):
        """失败时返回固定惩罚或启发式分数。"""
        reward = default_reward(sequence) if self.fallback_to_heuristic else self.failure_reward
        return ThermalFeedbackRecord(
            sequence=list(sequence),
            reward=reward,
            conductivity_w_mk=0.0,
            success=False,
            message=message,
            input_path=input_path,
        )

    def _notify(self, stage, message, details=None):
        """把搜索阶段的反馈进度交给外层界面或命令行。"""
        if self.progress_callback is None:
            return
        try:
            self.progress_callback(stage, message, details or {})
        except TypeError:
            self.progress_callback(stage, message)


def _is_lammps_thermo_row(text):
    """Return whether a streamed line starts with numeric thermo data."""
    fields = str(text).split()
    if len(fields) < 2:
        return False
    try:
        int(fields[0])
        float(fields[1])
    except ValueError:
        return False
    return True
