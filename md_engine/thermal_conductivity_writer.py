"""生成直接 NEMD 热导率计算输入文件。"""

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass
class ThermalConductivityInputResult:
    """记录一个候选体系的热导率输入和输出路径。"""

    data_path: str
    input_path: str
    correlation_output: str
    conductivity_output: str
    dump_output: str
    method: str
    success: bool
    message: str = ""


def write_direct_nemd_input(data_file, input_file, output_prefix, params=None):
    """生成冷热 Langevin 恒温浴驱动的直接 NEMD 输入。"""
    params = _merge_params(_default_params(), params or {})
    data_path = Path(data_file)
    input_path = Path(input_file)
    output_prefix = Path(output_prefix)
    profile_file = Path(f"{output_prefix}_temperature_profile.dat")
    transfer_file = Path(f"{output_prefix}_heat_flux.dat")
    dump_file = Path(f"{output_prefix}_nemd_dump.lammpstrj")
    density_file = Path(f"{output_prefix}_density.dat")

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
            method="direct_nemd",
            success=False,
            message=(
                "LAMMPS data is topology-only; missing force-field sections: "
                + ", ".join(missing)
            ),
        )

    profile_nevery = max(int(params.get("nemd_profile_nevery", 100)), 1)
    profile_repeat = max(int(params.get("nemd_profile_repeat", 20)), 1)
    profile_every = profile_nevery * profile_repeat
    production_steps = max(int(params["production_steps"]), profile_every)
    params["nemd_profile_nevery"] = profile_nevery
    params["nemd_profile_repeat"] = profile_repeat
    params["nemd_profile_every"] = profile_every
    params["production_steps"] = (
        (production_steps + profile_every - 1) // profile_every
    ) * profile_every

    density_sample_every = max(int(params["density_sample_every"]), 1)
    density_block_steps = max(int(params["density_block_steps"]), density_sample_every)
    density_window_samples = max(density_block_steps // density_sample_every, 1)
    params["density_sample_every"] = density_sample_every
    params["density_window_samples"] = density_window_samples
    params["density_block_steps"] = density_sample_every * density_window_samples
    params["density_max_blocks"] = max(int(params["density_max_blocks"]), 3)
    params["density_plateau_tolerance"] = max(
        float(params["density_plateau_tolerance"]), 1.0e-6
    )

    input_path.parent.mkdir(parents=True, exist_ok=True)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    script = _render_direct_nemd_script(
        data_file=data_path.as_posix(),
        profile_file=profile_file.as_posix(),
        transfer_file=transfer_file.as_posix(),
        dump_file=dump_file.as_posix(),
        density_file=density_file.as_posix(),
        params=params,
    )
    input_path.write_text(script, encoding="utf-8")

    return ThermalConductivityInputResult(
        data_path=str(data_path),
        input_path=str(input_path),
        correlation_output=str(profile_file),
        conductivity_output=str(transfer_file),
        dump_output=str(dump_file),
        method="direct_nemd",
        success=True,
        message="ok",
    )


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
        "nvt_steps": 10000,
        "npt_steps": 50000,
        "production_steps": 100000,
        "nemd_hot_temperature": 350.0,
        "nemd_cold_temperature": 250.0,
        "nemd_thermostat_damp": 100.0,
        "nemd_steady_steps": 100000,
        "nemd_profile_nevery": 100,
        "nemd_profile_repeat": 20,
        "nemd_min_gradient_r2": 0.70,
        "nemd_min_temperature_span": 5.0,
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
    }


def _merge_params(default_params, user_params):
    merged = dict(default_params)
    merged.update(user_params)
    return merged


def _style_line(keyword, value):
    value = str(value or "").strip()
    return f"{keyword:<16}{value}" if value and value.lower() != "none" else ""


def _as_bool(value):
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "on")
    return bool(value)


def _render_density_equilibration(density_file, params):
    """生成熔融、预压缩和三窗口密度平台判断流程。"""
    if not _as_bool(params.get("density_equilibration", True)):
        return f"""fix             eq_nvt all nvt temp ${{T}} ${{T}} ${{tdamp}}
run             {params["nvt_steps"]}
unfix           eq_nvt

fix             eq_npt all npt temp ${{T}} ${{T}} ${{tdamp}} iso ${{P}} ${{P}} ${{pdamp}}
run             {params["npt_steps"]}
unfix           eq_npt"""

    remaining_blocks = max(int(params["density_max_blocks"]) - 2, 1)
    return f"""# 高温熔融，释放初始构象中的局部应力。
variable        Tmelt equal {params["melt_temperature"]}
fix             melt_heat all nvt temp ${{T}} ${{Tmelt}} ${{tdamp}}
run             {params["melt_heating_steps"]}
unfix           melt_heat

fix             melt_hold all nvt temp ${{Tmelt}} ${{Tmelt}} ${{tdamp}}
run             {params["melt_hold_steps"]}
unfix           melt_hold

# 预压缩只消除 Packmol 大空隙，不规定最终平衡密度。
variable        rho_system equal (mass(all)/0.602214076)/vol
variable        rho_pre equal {params["precompression_density"]}
variable        pre_scale equal ternary(v_rho_system<v_rho_pre,(v_rho_system/v_rho_pre)^(1.0/3.0),1.0)
print           "Density before precompression: ${{rho_system}} g/cm^3"
print           "Precompression box scale: ${{pre_scale}}"
fix             pre_nvt all nvt temp ${{Tmelt}} ${{Tmelt}} ${{tdamp}}
fix             pre_box all deform 100 x scale ${{pre_scale}} y scale ${{pre_scale}} z scale ${{pre_scale}} remap x units box
run             {params["precompression_steps"]}
unfix           pre_box
unfix           pre_nvt

# 高温、目标压力下分块运行 NPT。连续三个窗口的平均密度进入平台后退出。
variable        rho_now equal density
variable        rho_tol equal {params["density_plateau_tolerance"]}
fix             density_npt all npt temp ${{Tmelt}} ${{Tmelt}} ${{tdamp}} iso ${{P}} ${{P}} ${{pdamp}}
fix             density_window all ave/time {params["density_sample_every"]} {params["density_window_samples"]} {params["density_block_steps"]} v_rho_now ave one file {density_file}
thermo          {params["density_block_steps"]}

run             {params["density_block_steps"]}
variable        rho_a equal $(f_density_window)
run             {params["density_block_steps"]}
variable        rho_b equal $(f_density_window)
variable        density_loop loop {remaining_blocks}

label           density_plateau_loop
run             {params["density_block_steps"]}
variable        rho_c equal $(f_density_window)
variable        rho_spread equal sqrt(((v_rho_a-v_rho_b)^2+(v_rho_b-v_rho_c)^2+(v_rho_a-v_rho_c)^2)/3.0)/((v_rho_a+v_rho_b+v_rho_c)/3.0)
print           "Density plateau check: rho=${{rho_c}} g/cm^3, relative spread=${{rho_spread}}"
if              "${{rho_spread}} <= ${{rho_tol}}" then "jump SELF density_plateau_done"
variable        rho_a delete
variable        rho_a equal ${{rho_b}}
variable        rho_b delete
variable        rho_b equal ${{rho_c}}
variable        rho_c delete
variable        rho_spread delete
next            density_loop
jump            SELF density_plateau_loop
jump            SELF density_plateau_limit

label           density_plateau_done
print           "Density plateau reached within tolerance."
jump            SELF density_plateau_cleanup

label           density_plateau_limit
print           "WARNING: Density plateau was not reached before the configured window limit."

label           density_plateau_cleanup
unfix           density_window
unfix           density_npt
print           "Density plateau stage finished at: ${{rho_system}} g/cm^3"

# 在目标压力下降温，再完成常温 NPT 与 NVT 平衡。
fix             cool_npt all npt temp ${{Tmelt}} ${{T}} ${{tdamp}} iso ${{P}} ${{P}} ${{pdamp}}
run             {params["cooling_steps"]}
unfix           cool_npt

fix             final_npt all npt temp ${{T}} ${{T}} ${{tdamp}} iso ${{P}} ${{P}} ${{pdamp}}
run             {params["npt_steps"]}
unfix           final_npt

fix             final_nvt all nvt temp ${{T}} ${{T}} ${{tdamp}}
run             {params["nvt_steps"]}
unfix           final_nvt
print           "Final equilibrated density: ${{rho_system}} g/cm^3"
"""


def _render_direct_nemd_script(
    data_file,
    profile_file,
    transfer_file,
    dump_file,
    density_file,
    params,
):
    """按照参考算例渲染冷热区直接 NEMD 脚本。"""
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
    equilibration_block = _render_density_equilibration(density_file, params)

    return f"""# 冷热 Langevin 恒温浴直接 NEMD 热导率计算
units           real
dimension       3
atom_style      full
boundary        p p p

variable        T equal {params["temperature"]}
variable        P equal {params["pressure"]}
variable        dt equal {params["timestep"]}
variable        tdamp equal {params["tdamp"]}
variable        pdamp equal {params["pdamp"]}

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

{equilibration_block}
reset_timestep  0

# 沿 x 方向依次划分固定端、热源、传热主体、冷源和固定端。
variable        fix_thick equal 0.05*lx
variable        source_thick equal 0.10*lx
variable        flux_thick equal 0.10*lx
region          bottomfix block $(xlo) $(xlo+v_fix_thick) EDGE EDGE EDGE EDGE units box
region          hot block $(xlo+v_fix_thick) $(xlo+v_fix_thick+v_source_thick) EDGE EDGE EDGE EDGE units box
region          cold block $(xhi-v_fix_thick-v_source_thick) $(xhi-v_fix_thick) EDGE EDGE EDGE EDGE units box
region          topfix block $(xhi-v_fix_thick) $(xhi) EDGE EDGE EDGE EDGE units box
region          flux block $(0.5*(xlo+xhi)-0.5*v_flux_thick) $(0.5*(xlo+xhi)+0.5*v_flux_thick) EDGE EDGE EDGE EDGE units box
group           bottomfix region bottomfix
group           topfix region topfix
group           hot region hot
group           cold region cold
group           flux region flux

compute         Thot all temp/region hot
compute         Tcold all temp/region cold
fix             fhot all langevin {params["nemd_hot_temperature"]} {params["nemd_hot_temperature"]} {params["nemd_thermostat_damp"]} {params["random_seed"]}
fix             fcold all langevin {params["nemd_cold_temperature"]} {params["nemd_cold_temperature"]} {params["nemd_thermostat_damp"]} {int(params["random_seed"]) + 104729}
fix_modify      fhot temp Thot
fix_modify      fcold temp Tcold
fix             hold_bottom bottomfix setforce 0.0 0.0 0.0
fix             hold_top topfix setforce 0.0 0.0 0.0
group           mobile subtract all bottomfix topfix
fix             nemd_nve mobile nve
thermo_style    custom step temp press density c_Thot c_Tcold
run             {params["nemd_steady_steps"]}

reset_timestep  0
compute         myKE flux ke/atom
compute         myPE flux pe/atom
compute         myStress flux stress/atom NULL virial
compute         heat_flux flux heat/flux myKE myPE myStress
variable        flux_volume equal ly*lz*v_flux_thick
variable        Jx equal c_heat_flux[1]/v_flux_volume
fix             average_flux flux ave/time {params["nemd_profile_nevery"]} {params["nemd_profile_repeat"]} {params["nemd_profile_every"]} v_Jx ave one file {transfer_file}

compute         layers all chunk/atom bin/1d x lower 0.5 units box
fix             temperature_profile all ave/chunk {params["nemd_profile_nevery"]} {params["nemd_profile_repeat"]} {params["nemd_profile_every"]} layers temp norm sample file {profile_file}
dump            trajectory all custom {params["dump_every"]} {dump_file} id mol type q x y z vx vy vz
thermo          {params["nemd_profile_every"]}
thermo_style    custom step temp press density c_Thot c_Tcold v_Jx
thermo_modify   flush yes
run             {params["production_steps"]}

undump          trajectory
unfix           nemd_nve
unfix           fhot
unfix           fcold
unfix           hold_bottom
unfix           hold_top
unfix           average_flux
unfix           temperature_profile
"""
