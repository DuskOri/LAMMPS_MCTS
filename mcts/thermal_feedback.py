"""把热导率计算结果接回 MCTS 奖励函数。"""

from dataclasses import dataclass, field
from pathlib import Path

from generator import build_polymer_from_sequence, prepare_initial_system
from md_engine import run_lammps_input, write_fast_tc_input
from mcts.mcts_engine import default_reward
from post_process import analyze_thermal_outputs


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
    max_evaluations: int = 10
    run_lammps: bool = False
    reward_offset: float = 0.05
    fallback_to_heuristic: bool = True
    failure_reward: float = -1.0
    cache: dict = field(default_factory=dict)
    records: list = field(default_factory=list)

    def __call__(self, sequence):
        """作为 reward_fn 被 MCTS 调用。"""
        key = tuple(sequence)
        if key in self.cache:
            return self.cache[key]

        if len(self.cache) >= self.max_evaluations:
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
        return record.reward

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

        build_result = build_polymer_from_sequence(
            sequence,
            dp=polymer_config["degree_of_polymerization"],
            output_dir=pdb_dir,
            name="candidate",
        )
        if not build_result.success:
            return self._failed(sequence, build_result.message)

        system_result = prepare_initial_system(
            template_pdb=build_result.pdb_path,
            output_dir=system_dir,
            molecule_count=system_config["molecule_count"],
            box_size=system_config["box_size"],
            tolerance=system_config["tolerance"],
            packmol_executable=system_config["packmol_executable"],
        )
        if not system_result.success:
            return self._failed(sequence, system_result.message)

        input_result = write_fast_tc_input(
            data_file=system_result.system_data_path,
            input_file=lammps_dir / "fast_tc.in",
            output_prefix=lammps_dir / "fast_tc",
            params={
                key: value
                for key, value in tc_config.items()
                if key not in ("write_fast_input", "analyze_outputs")
            },
        )

        if self.run_lammps:
            run_result = run_lammps_input(
                input_path=input_result.input_path,
                lammps_executable=lammps_config["executable"],
                timeout=lammps_config.get("timeout_seconds"),
            )
            if not run_result.success:
                return self._failed(sequence, run_result.message, input_result.input_path)

        analysis = analyze_thermal_outputs(
            input_path=input_result.input_path,
            flux_path=input_result.flux_output,
            temp_path=input_result.temp_output,
            trim_fraction=tc_config.get("trim_fraction", 0.20),
        )
        if not analysis.success:
            return self._failed(sequence, analysis.message, input_result.input_path)

        return ThermalFeedbackRecord(
            sequence=list(sequence),
            reward=self._conductivity_to_reward(analysis.conductivity_w_mk),
            conductivity_w_mk=analysis.conductivity_w_mk,
            success=True,
            message="ok",
            input_path=input_result.input_path,
        )

    def _conductivity_to_reward(self, conductivity_w_mk):
        """热导率越低，奖励越高。"""
        conductivity = max(conductivity_w_mk, 0.0)
        return 1.0 / (conductivity + self.reward_offset)

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
