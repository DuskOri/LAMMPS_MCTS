"""把热导率计算结果接回 MCTS 奖励函数。"""

from dataclasses import dataclass, field
from pathlib import Path

from generator import build_polymer_from_sequence, prepare_initial_system
from md_engine import run_lammps_input, write_direct_nemd_input
from mcts.mcts_engine import default_reward
from post_process import analyze_direct_nemd_outputs


@dataclass
class ThermalFeedbackRecord:
    """记录一次序列热导率评估过程。"""

    sequence: list
    reward: float
    conductivity_w_mk: float
    success: bool
    message: str
    input_path: str = ""


@dataclass
class ThermalFeedbackEvaluator:
    """把聚合物序列转换为热导率奖励值。"""

    config: dict
    output_dir: str = "outputs/mcts_feedback"
    max_evaluations: object = 10
    run_lammps: bool = False
    reward_scale: float = 0.20
    fallback_to_heuristic: bool = False
    failure_reward: float = -1.0
    progress_callback: object = None
    cache: dict = field(default_factory=dict)
    records: list = field(default_factory=list)

    def __call__(self, sequence):
        """作为 reward_fn 被 MCTS 调用。"""
        key = tuple(sequence)
        if key in self.cache:
            return self.cache[key]

        if self._limit_reached():
            reward = default_reward(sequence) if self.fallback_to_heuristic else self.failure_reward
            self.cache[key] = reward
            self.records.append(
                ThermalFeedbackRecord(
                    sequence=list(sequence),
                    reward=reward,
                    conductivity_w_mk=0.0,
                    success=False,
                    message="thermal evaluation limit reached",
                )
            )
            return reward

        record = self._evaluate(sequence)
        self.cache[key] = record.reward
        self.records.append(record)
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
        index = len(self.cache) + 1
        base_dir = Path(self.output_dir) / f"eval_{index:04d}"
        pdb_dir = base_dir / "pdb"
        system_dir = base_dir / "systems"
        lammps_dir = base_dir / "lammps"

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
        )
        if not system_result.success:
            return self._failed(sequence, system_result.message)

        self._notify(
            "mcts_lammps",
            f"MCTS thermal feedback {index}: writing LAMMPS input",
        )
        input_stem = "direct_nemd"
        input_result = write_direct_nemd_input(
            data_file=system_result.system_data_path,
            input_file=lammps_dir / f"{input_stem}.in",
            output_prefix=lammps_dir / input_stem,
            params={
                key: value
                for key, value in tc_config.items()
                if key not in ("write_fast_input", "analyze_outputs", "profile")
            },
        )
        if not input_result.success:
            return self._failed(sequence, input_result.message)

        if self.run_lammps:
            density_equilibration = bool(tc_config.get("density_equilibration", True))
            initial_stage = "compress" if density_equilibration else "run_lammps"
            initial_message = (
                f"MCTS thermal feedback {index}: melting and compressing the system"
                if density_equilibration
                else f"MCTS thermal feedback {index}: running thermal simulation"
            )
            self._notify(initial_stage, initial_message)

            nemd_started = False
            nemd_run_index = 0
            nemd_run_target = int(tc_config.get("nemd_steady_steps", 0) or 0)

            def report_lammps_line(line):
                nonlocal nemd_started, nemd_run_index, nemd_run_target
                text = str(line).strip()
                if text.startswith("Density plateau check"):
                    self._notify("compress", text, {"console_line": text})
                elif text.startswith("Final equilibrated density"):
                    nemd_started = True
                    self._notify(
                        "run_lammps",
                        f"{text}; starting thermal-conductivity sampling",
                        {"console_line": text},
                    )
                elif nemd_started and text.startswith("run "):
                    fields = text.split()
                    if len(fields) >= 2 and fields[1].isdigit():
                        nemd_run_index += 1
                        nemd_run_target = int(fields[1])
                    phase = "establishing steady state" if nemd_run_index <= 1 else "sampling"
                    self._notify(
                        "run_lammps",
                        f"NEMD {phase}: 0/{nemd_run_target} steps",
                        {"console_line": text},
                    )
                elif nemd_started and _is_lammps_thermo_row(text):
                    current_step = int(text.split()[0])
                    phase = "establishing steady state" if nemd_run_index <= 1 else "sampling"
                    self._notify(
                        "run_lammps",
                        f"NEMD {phase}: {current_step}/{nemd_run_target} steps",
                        {"console_line": text},
                    )
                else:
                    self._notify(
                        "run_lammps" if nemd_started else "compress",
                        text,
                        {"console_line": text, "console_only": True},
                    )

            run_result = run_lammps_input(
                input_path=input_result.input_path,
                lammps_executable=lammps_config["executable"],
                timeout=lammps_config.get("timeout_seconds") or None,
                line_callback=report_lammps_line,
            )
            if not run_result.success:
                return self._failed(sequence, run_result.message, input_result.input_path)

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
        )
        if not analysis.success:
            return self._failed(sequence, analysis.message, input_result.input_path)

        analysis_message = (
            f"ok; gradient={analysis.gradient_k_per_a:.6f} K/A; "
            f"R2={analysis.gradient_r2:.4f}"
        )
        return ThermalFeedbackRecord(
            sequence=list(sequence),
            reward=self._conductivity_to_reward(analysis.conductivity_w_mk),
            conductivity_w_mk=analysis.conductivity_w_mk,
            success=True,
            message=analysis_message,
            input_path=input_result.input_path,
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
