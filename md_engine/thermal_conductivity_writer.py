"""生成快速 Green-Kubo 热导率计算输入文件。"""

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass
class ThermalConductivityInputResult:
    """记录一个候选体系的 Green-Kubo 输入和输出路径。"""

    data_path: str
    input_path: str
    correlation_output: str
    conductivity_output: str
    dump_output: str
    method: str
    success: bool
    message: str = ""

def write_rapid_gk_input(data_file, input_file, output_prefix, params=None):
    """生成忽略分子间热流交叉相关项的快速 Green-Kubo 输入。"""
    params = _merge_params(_default_params(), params or {})
    data_path = Path(data_file)
    input_path = Path(input_file)
    output_prefix = Path(output_prefix)
    correlation_file = Path(f"{output_prefix}_intra_hfacf.dat")
    conductivity_file = Path(f"{output_prefix}_kappa.dat")
    dump_file = Path(f"{output_prefix}_dump.lammpstrj")

    missing = validate_force_field_data(
        data_path,
        force_field_include=params.get("force_field_include", ""),
    )
    if missing:
        return ThermalConductivityInputResult(
            data_path=str(data_path),
            input_path="",
            correlation_output="",
            conductivity_output="",
            dump_output="",
            method="rapid_green_kubo",
            success=False,
            message=(
                "LAMMPS data is topology-only; missing force-field sections: "
                + ", ".join(missing)
            ),
        )

    molecule_count = max(int(params.get("molecule_count", 1)), 1)
    sample_nevery = max(int(params["sample_nevery"]), 1)
    correlation_samples = max(int(params["correlation_samples"]), 2)
    correlation_every = sample_nevery * correlation_samples
    params["sample_nevery"] = sample_nevery
    params["correlation_samples"] = correlation_samples
    params["correlation_every"] = correlation_every
    params["production_steps"] = max(
        int(params["production_steps"]), correlation_every
    )

    input_path.parent.mkdir(parents=True, exist_ok=True)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    script = _render_rapid_gk_script(
        data_file=data_path.as_posix(),
        correlation_file=correlation_file.as_posix(),
        conductivity_file=conductivity_file.as_posix(),
        dump_file=dump_file.as_posix(),
        molecule_count=molecule_count,
        params=params,
    )
    input_path.write_text(script, encoding="utf-8")

    return ThermalConductivityInputResult(
        data_path=str(data_path),
        input_path=str(input_path),
        correlation_output=str(correlation_file),
        conductivity_output=str(conductivity_file),
        dump_output=str(dump_file),
        method="rapid_green_kubo",
        success=True,
        message="ok",
    )


def write_fast_tc_input(data_file, input_file, output_prefix, params=None):
    """保留旧函数名，实际生成快速 Green-Kubo 输入。"""
    return write_rapid_gk_input(data_file, input_file, output_prefix, params=params)


def validate_force_field_data(data_file, force_field_include=""):
    """确认 data 或额外 include 中存在运行 MD 所需的力场系数。"""
    data_path = Path(data_file)
    if not data_path.exists():
        return ["data file"]

    if force_field_include:
        include_path = Path(force_field_include)
        if include_path.exists():
            return []
        return ["force-field include file"]

    text = data_path.read_text(encoding="utf-8", errors="ignore")
    required = []
    if _type_count(text, "atom") > 0 and not _has_any_section(
        text, ("Pair Coeffs", "PairIJ Coeffs")
    ):
        required.append("Pair Coeffs")
    for type_name, section in (
        ("bond", "Bond Coeffs"),
        ("angle", "Angle Coeffs"),
        ("dihedral", "Dihedral Coeffs"),
        ("improper", "Improper Coeffs"),
    ):
        if _type_count(text, type_name) > 0 and not _has_any_section(text, (section,)):
            required.append(section)
    return required


def _type_count(text, type_name):
    pattern = rf"^\s*(\d+)\s+{re.escape(type_name)}\s+types\s*$"
    match = re.search(pattern, text, flags=re.MULTILINE | re.IGNORECASE)
    return int(match.group(1)) if match else 0


def _has_any_section(text, names):
    return any(
        re.search(rf"^\s*{re.escape(name)}\s*$", text, flags=re.MULTILINE)
        for name in names
    )


def _default_params():
    """快速筛选默认值；正式复核时应延长平衡和采样时间。"""
    return {
        "temperature": 300.0,
        "pressure": 1.0,
        "timestep": 0.5,
        "tdamp": 50.0,
        "pdamp": 500.0,
        "cutoff": 12.0,
        "nvt_steps": 10000,
        "npt_steps": 50000,
        "production_steps": 100000,
        "sample_nevery": 10,
        "correlation_samples": 1000,
        "thermo_every": 1000,
        "dump_every": 10000,
        "random_seed": 87287,
        "pair_style": "lj/cut/coul/long 12.0",
        "bond_style": "harmonic",
        "angle_style": "harmonic",
        "dihedral_style": "harmonic",
        "improper_style": "cvff",
        "kspace_style": "pppm 1.0e-4",
        "special_bonds": "amber",
        "force_field_include": "",
        "molecule_count": 1,
    }


def _merge_params(default_params, user_params):
    merged = dict(default_params)
    merged.update(user_params)
    return merged


def _style_line(keyword, value):
    value = str(value or "").strip()
    return f"{keyword:<16}{value}" if value and value.lower() != "none" else ""


def _render_rapid_gk_script(
    data_file,
    correlation_file,
    conductivity_file,
    dump_file,
    molecule_count,
    params,
):
    """渲染按分子求和的热流自相关计算脚本。"""
    style_lines = [
        _style_line("pair_style", params["pair_style"]),
        _style_line("bond_style", params["bond_style"]),
        _style_line("angle_style", params["angle_style"]),
        _style_line("dihedral_style", params["dihedral_style"]),
        _style_line("improper_style", params["improper_style"]),
    ]
    style_block = "\n".join(line for line in style_lines if line)
    include_value = params.get("force_field_include", "")
    include_file = "" if include_value is None else str(include_value).strip()
    include_line = f"include         {Path(include_file).as_posix()}" if include_file else ""
    kspace_line = _style_line("kspace_style", params.get("kspace_style", ""))
    special_line = _style_line("special_bonds", params.get("special_bonds", ""))

    molecule_blocks = []
    heat_flux_values = []
    for molecule_id in range(1, molecule_count + 1):
        molecule_blocks.extend(
            [
                f"group           mol_{molecule_id} molecule {molecule_id}",
                f"compute         ke_{molecule_id} mol_{molecule_id} ke/atom",
                f"compute         pe_{molecule_id} mol_{molecule_id} pe/atom",
                f"compute         stress_{molecule_id} mol_{molecule_id} stress/atom NULL virial",
                (
                    f"compute         flux_{molecule_id} mol_{molecule_id} heat/flux "
                    f"ke_{molecule_id} pe_{molecule_id} stress_{molecule_id}"
                ),
            ]
        )
        heat_flux_values.extend(
            f"c_flux_{molecule_id}[{axis}]" for axis in (1, 2, 3)
        )

    molecule_block = "\n".join(molecule_blocks)
    correlation_values = " ".join(heat_flux_values)
    first_column = 3
    last_column = first_column + len(heat_flux_values) - 1
    trap_sum = "+".join(f"trap(f_HFACF[{column}])" for column in range(first_column, last_column + 1))

    return f"""# 快速 Green-Kubo 热导率计算
# 分子间作用力保留，只在热流相关函数中省略不同分子之间的交叉项。
units           real
dimension       3
atom_style      full
boundary        p p p

variable        T equal {params["temperature"]}
variable        P equal {params["pressure"]}
variable        dt equal {params["timestep"]}
variable        tdamp equal {params["tdamp"]}
variable        pdamp equal {params["pdamp"]}
variable        s equal {params["sample_nevery"]}
variable        p equal {params["correlation_samples"]}
variable        d equal {params["correlation_every"]}

{style_block}
read_data       {data_file}
{include_line}
{special_line}
{kspace_line}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        ${{dt}}
thermo          {params["thermo_every"]}
thermo_style    custom step temp press density vol pe ke etotal
thermo_modify   flush yes

minimize        1.0e-4 1.0e-6 500 5000
velocity        all create ${{T}} {params["random_seed"]} mom yes rot yes dist gaussian

fix             eq_nvt all nvt temp ${{T}} ${{T}} ${{tdamp}}
run             {params["nvt_steps"]}
unfix           eq_nvt

fix             eq_npt all npt temp ${{T}} ${{T}} ${{tdamp}} iso ${{P}} ${{P}} ${{pdamp}}
run             {params["npt_steps"]}
unfix           eq_npt
reset_timestep  0

{molecule_block}

# type auto 只保留每个分子三个热流分量各自的自相关，不生成分子间交叉相关。
fix             HFACF all ave/correlate ${{s}} ${{p}} ${{d}} {correlation_values} type auto file {correlation_file} ave running

# real 单位转换：kcal/(mol*Angstrom*fs) -> W/m，并除以 kB*T^2*V。
variable        kB equal 1.380649e-23
variable        kcal2J equal 4184.0/6.02214076e23
variable        A2m equal 1.0e-10
variable        fs2s equal 1.0e-15
variable        convert equal v_kcal2J*v_kcal2J/v_fs2s/v_A2m
variable        scale equal v_convert/v_kB/${{T}}/${{T}}/vol*${{s}}*${{dt}}
variable        kappa equal ({trap_sum})*v_scale/3.0

fix             kappa_out all ave/time ${{d}} 1 ${{d}} v_kappa file {conductivity_file}
fix             production all nve
dump            trajectory all custom {params["dump_every"]} {dump_file} id mol type q x y z vx vy vz
thermo_style    custom step temp press density etotal v_kappa
thermo_modify   flush yes
run             {params["production_steps"]}

undump          trajectory
unfix           production
unfix           kappa_out
unfix           HFACF
"""
