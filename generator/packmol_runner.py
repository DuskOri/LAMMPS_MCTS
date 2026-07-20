"""调用 Packmol 生成多链初始体系。"""

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess

from md_engine import (
    read_pdb_topology,
    write_lammps_data,
    write_lammps_data_from_template,
    write_uff_lammps_data_from_template,
)


@dataclass
class InitialSystemResult:
    """记录单个候选结构的初始体系准备结果。"""

    template_pdb: str
    single_data_path: str
    packmol_input_path: str
    packed_pdb_path: str
    system_data_path: str
    molecule_count: int
    success: bool
    message: str = ""
    force_field: str = "topology_only"
    force_field_ready: bool = False


def prepare_initial_system(
    template_pdb,
    output_dir,
    molecule_count=10,
    box_size=60.0,
    tolerance=2.0,
    packmol_executable="packmol",
    template_mol=None,
    force_field="uff_screening",
):
    """由单链 PDB 生成单链 data，并尝试通过 Packmol 生成多链初始体系。"""
    template_path = Path(template_pdb)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    stem = template_path.stem
    single_data = output_path / f"{stem}_single.data"
    packmol_input = output_path / f"{stem}_packmol.inp"
    packed_pdb = output_path / f"{stem}_system.pdb"
    system_data = output_path / f"{stem}_system.data"

    topology = read_pdb_topology(template_path)
    write_lammps_data(topology, single_data, comment=f"single chain from {template_path.name}")
    _write_packmol_input(
        template_path,
        packmol_input,
        packed_pdb,
        molecule_count=molecule_count,
        box_size=box_size,
        tolerance=tolerance,
    )

    executable = _resolve_packmol(packmol_executable)
    if executable is None:
        return InitialSystemResult(
            template_pdb=str(template_path),
            single_data_path=str(single_data),
            packmol_input_path=str(packmol_input),
            packed_pdb_path="",
            system_data_path="",
            molecule_count=molecule_count,
            success=False,
            message="Packmol executable not found.",
        )

    run_result = _run_packmol(executable, packmol_input)
    if run_result.returncode != 0:
        return InitialSystemResult(
            template_pdb=str(template_path),
            single_data_path=str(single_data),
            packmol_input_path=str(packmol_input),
            packed_pdb_path="",
            system_data_path="",
            molecule_count=molecule_count,
            success=False,
            message=run_result.stderr.strip() or run_result.stdout.strip(),
        )

    try:
        if force_field == "uff_screening" and template_mol:
            write_uff_lammps_data_from_template(
                template_mol=template_mol,
                coordinate_pdb=packed_pdb,
                filename=system_data,
                molecule_count=molecule_count,
                box_size=box_size,
            )
            force_field_ready = True
        else:
            write_lammps_data_from_template(
                template_pdb=template_path,
                coordinate_pdb=packed_pdb,
                filename=system_data,
                molecule_count=molecule_count,
                box_size=box_size,
            )
            force_field_ready = False
    except (OSError, ValueError) as error:
        return InitialSystemResult(
            template_pdb=str(template_path),
            single_data_path=str(single_data),
            packmol_input_path=str(packmol_input),
            packed_pdb_path=str(packed_pdb),
            system_data_path="",
            molecule_count=molecule_count,
            success=False,
            message=f"force-field parameterization failed: {error}",
            force_field=str(force_field),
            force_field_ready=False,
        )

    return InitialSystemResult(
        template_pdb=str(template_path),
        single_data_path=str(single_data),
        packmol_input_path=str(packmol_input),
        packed_pdb_path=str(packed_pdb),
        system_data_path=str(system_data),
        molecule_count=molecule_count,
        success=True,
        message="ok",
        force_field=str(force_field),
        force_field_ready=force_field_ready,
    )


def _write_packmol_input(template_pdb, filename, output_pdb, molecule_count, box_size, tolerance):
    """写 Packmol 输入文件，约束所有链在立方盒子内随机排布。"""
    input_text = "\n".join(
        [
            f"tolerance {tolerance}",
            "filetype pdb",
            f"output {_packmol_path(output_pdb)}",
            "",
            f"structure {_packmol_path(template_pdb)}",
            f"  number {molecule_count}",
            f"  inside box 0.0 0.0 0.0 {box_size} {box_size} {box_size}",
            "end structure",
            "",
        ]
    )
    Path(filename).write_text(input_text, encoding="utf-8")


def _resolve_packmol(packmol_executable):
    """定位 Packmol 可执行文件。"""
    executable = Path(packmol_executable)
    if executable.exists():
        return str(executable)
    return shutil.which(packmol_executable)


def _run_packmol(executable, input_file):
    """运行 Packmol，优先使用新版支持的 -i 参数。"""
    return subprocess.run(
        [executable, "-i", str(input_file)],
        capture_output=True,
        check=False,
        text=True,
    )


def _packmol_path(path):
    """Packmol 输入中统一使用正斜杠路径。"""
    return Path(path).resolve().as_posix()
