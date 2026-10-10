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


@dataclass
class LammpsStageInputResult:
    """Paths produced while preparing one standalone LAMMPS stage."""

    data_path: str
    input_path: str
    output_data_path: str
    success: bool
    message: str = ""


def write_compression_input(
    data_file,
    input_file,
    compressed_data_file,
    output_prefix,
    params=None,
):
    """Write the melt, density plateau and cooling stage as a separate job."""
    params = _normalise_params(params or {})
    data_path = Path(data_file)
    input_path = Path(input_file)
    compressed_path = Path(compressed_data_file)
    density_file = Path(f"{output_prefix}_density.dat")
    missing = validate_force_field_data(
        data_path,
        force_field_include=params.get("force_field_include", ""),
    )
    if missing:
        return LammpsStageInputResult(
            str(data_path), "", str(compressed_path), False,
            "LAMMPS data is topology-only; missing force-field sections: "
            + ", ".join(missing),
        )

    input_path.parent.mkdir(parents=True, exist_ok=True)
    compressed_path.parent.mkdir(parents=True, exist_ok=True)
    density_file.parent.mkdir(parents=True, exist_ok=True)
    script = _render_compression_script(
        data_path.as_posix(),
        compressed_path.as_posix(),
        density_file.as_posix(),
        params,
    )
    input_path.write_text(script, encoding="utf-8")
    return LammpsStageInputResult(
        str(data_path), str(input_path), str(compressed_path), True, "ok"
    )


def write_adaptive_nemd_input(data_file, input_file, output_prefix, params=None):
    """Write a NEMD job whose production loop stops after heat-flux convergence."""
    params = _normalise_params(params or {})
    data_path = Path(data_file)
    input_path = Path(input_file)
    output_prefix = Path(output_prefix)
    profile_file = Path(f"{output_prefix}_temperature_profile.dat")
    transfer_file = Path(f"{output_prefix}_heat_flux.dat")
    dump_file = Path(f"{output_prefix}_nemd_dump.lammpstrj")
    final_data_file = Path(f"{output_prefix}_final.data")
    missing = validate_force_field_data(
        data_path,
        force_field_include=params.get("force_field_include", ""),
    )
    if missing:
        return ThermalConductivityInputResult(
            str(data_path), "", "", "", "", "direct_nemd", False,
            "LAMMPS data is topology-only; missing force-field sections: "
            + ", ".join(missing),
        )

    input_path.parent.mkdir(parents=True, exist_ok=True)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    script = _render_adaptive_nemd_script(
        data_path.as_posix(),
        profile_file.as_posix(),
        transfer_file.as_posix(),
        dump_file.as_posix(),
        final_data_file.as_posix(),
        params,
    )
    input_path.write_text(script, encoding="utf-8")
    return ThermalConductivityInputResult(
        str(data_path), str(input_path), str(profile_file), str(transfer_file),
        str(dump_file), "direct_nemd", True, "ok",
    )


def write_direct_nemd_input(data_file, input_file, output_prefix, params=None):
    """生成冷热 Langevin 恒温浴驱动的直接 NEMD 输入。"""
    params = _normalise_params(params or {})
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
        re.search(
            rf"^\s*{re.escape(name)}(?:\s+#.*)?\s*$",
            text,
            flags=re.MULTILINE,
        )
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
        "nemd_steady_steps": 20000,
        "nemd_production_block_steps": 10000,
        "nemd_min_production_blocks": 5,
        "nemd_heat_flux_tolerance": 0.10,
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


def _normalise_params(user_params):
    """Merge defaults and align averaging, density and production windows."""
    params = _merge_params(_default_params(), user_params)
    profile_nevery = max(int(params.get("nemd_profile_nevery", 100)), 1)
    profile_repeat = max(int(params.get("nemd_profile_repeat", 20)), 1)
    profile_every = profile_nevery * profile_repeat
    block_steps = max(
        int(params.get("nemd_production_block_steps", profile_every)),
        profile_every,
    )
    block_steps = ((block_steps + profile_every - 1) // profile_every) * profile_every
    min_blocks = max(int(params.get("nemd_min_production_blocks", 5)), 3)
    production_steps = int(params.get("production_steps", block_steps))
    max_blocks = production_steps // block_steps
    if max_blocks < min_blocks:
        raise ValueError(f"production_steps must be at least {min_blocks * block_steps}; configured {production_steps}")
    if not (float(params["nemd_hot_temperature"]) > float(params["nemd_cold_temperature"]) > 0):
        raise ValueError("NEMD hot temperature must exceed cold temperature, both positive")
    for key in ("timestep", "tdamp", "pdamp", "nemd_thermostat_damp"):
        if float(params[key]) <= 0:
            raise ValueError(f"{key} must be positive")
    density_sample_every = max(int(params["density_sample_every"]), 1)
    density_block_steps = max(int(params["density_block_steps"]), density_sample_every)

    params.update(
        {
            "nemd_profile_nevery": profile_nevery,
            "nemd_profile_repeat": profile_repeat,
            "nemd_profile_every": profile_every,
            "nemd_production_block_steps": block_steps,
            "nemd_min_production_blocks": min_blocks,
            "nemd_max_production_blocks": max_blocks,
            "production_steps": max_blocks * block_steps,
            "nemd_heat_flux_tolerance": max(
                float(params.get("nemd_heat_flux_tolerance", 0.10)), 1.0e-6
            ),
            "density_sample_every": density_sample_every,
            "density_window_samples": max(
                density_block_steps // density_sample_every, 1
            ),
            "density_max_blocks": max(int(params["density_max_blocks"]), 3),
            "density_plateau_tolerance": max(
                float(params["density_plateau_tolerance"]), 1.0e-6
            ),
        }
    )
    params["density_block_steps"] = (
        params["density_sample_every"] * params["density_window_samples"]
    )
    return params


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
fix             density_window all ave/time {params["density_sample_every"]} {params["density_window_samples"]} {params["density_block_steps"]} v_rho_now ave one file "{density_file}"
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
quit            2

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


def _render_lammps_header(data_file, params):
    style_lines = [
        _style_line("pair_style", params["pair_style"]),
        _style_line("bond_style", params["bond_style"]),
        _style_line("angle_style", params["angle_style"]),
        _style_line("dihedral_style", params["dihedral_style"]),
        _style_line("improper_style", params["improper_style"]),
    ]
    style_block = "\n".join(line for line in style_lines if line)
    include_file = str(params.get("force_field_include", "") or "").strip()
    include_line = f'include         "{Path(include_file).as_posix()}"' if include_file else ""
    kspace_line = _style_line("kspace_style", params.get("kspace_style", ""))
    special_line = _style_line("special_bonds", params.get("special_bonds", ""))
    return f"""units           real
dimension       3
atom_style      full
boundary        p p p

variable        T equal {params["temperature"]}
variable        P equal {params["pressure"]}
variable        dt equal {params["timestep"]}
variable        tdamp equal {params["tdamp"]}
variable        pdamp equal {params["pdamp"]}

{style_block}
read_data       "{data_file}"
{include_line}
{special_line}
{kspace_line}

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
timestep        ${{dt}}
thermo          {params["thermo_every"]}
thermo_style    custom step temp press density vol pe ke etotal
thermo_modify   flush yes"""


def _render_compression_script(data_file, compressed_data_file, density_file, params):
    header = _render_lammps_header(data_file, params)
    equilibration = _render_density_equilibration(density_file, params)
    return f"""# Melt, density equilibration and room-temperature relaxation
{header}

minimize        1.0e-4 1.0e-6 500 5000
velocity        all create ${{T}} {params["random_seed"]} mom yes rot yes dist gaussian

{equilibration}
write_data      "{compressed_data_file}"
print           "Compressed structure written: {compressed_data_file}"
"""


def _render_flux_window_update(block_steps):
    return f"""run             {block_steps}
variable        flux_a delete
variable        flux_a equal ${{flux_b}}
variable        flux_b delete
variable        flux_b equal ${{flux_c}}
variable        flux_c delete
variable        flux_c equal $(f_convergence_flux)"""


def _render_flux_check(params):
    return f"""variable        flux_mean equal (v_flux_a+v_flux_b+v_flux_c)/3.0
variable        flux_floor equal ternary(abs(v_flux_mean)>1.0e-30,abs(v_flux_mean),1.0e-30)
variable        flux_spread equal sqrt(((v_flux_a-v_flux_b)^2+(v_flux_b-v_flux_c)^2+(v_flux_a-v_flux_c)^2)/3.0)/v_flux_floor
print           "Heat flux convergence check: relative spread=${{flux_spread}}, tolerance={params["nemd_heat_flux_tolerance"]}"
if              "${{flux_spread}} <= {params["nemd_heat_flux_tolerance"]} && ${{flux_a}} > 0 && ${{flux_b}} > 0 && ${{flux_c}} > 0" then "jump SELF heat_flux_converged"
variable        flux_mean delete
variable        flux_floor delete
variable        flux_spread delete"""


def _render_adaptive_production(params):
    block_steps = params["nemd_production_block_steps"]
    min_blocks = params["nemd_min_production_blocks"]
    max_blocks = params["nemd_max_production_blocks"]
    warmup_updates = "\n".join(
        _render_flux_window_update(block_steps) for _ in range(max(min_blocks - 3, 0))
    )
    remaining = max(max_blocks - min_blocks, 0)
    continuation = ""
    if remaining:
        continuation = f"""
variable        production_loop loop {remaining}
label           heat_flux_production_loop
{_render_flux_window_update(block_steps)}
{_render_flux_check(params)}
next            production_loop
jump            SELF heat_flux_production_loop
jump            SELF heat_flux_limit
"""
    return f"""variable        heat_flux_ok internal 0
print           "NEMD production starts: block={block_steps} steps, minimum={min_blocks} blocks, maximum={max_blocks} blocks"
run             {block_steps}
variable        flux_a equal $(f_convergence_flux)
run             {block_steps}
variable        flux_b equal $(f_convergence_flux)
run             {block_steps}
variable        flux_c equal $(f_convergence_flux)
{warmup_updates}
{_render_flux_check(params)}
{continuation}
jump            SELF heat_flux_limit

label           heat_flux_converged
variable        heat_flux_ok delete
variable        heat_flux_ok internal 1
print           "Heat flux converged; production sampling stopped early."
jump            SELF heat_flux_cleanup

label           heat_flux_limit
print           "WARNING: Heat flux did not converge before the production limit."

label           heat_flux_cleanup"""


def _render_adaptive_nemd_script(
    data_file,
    profile_file,
    transfer_file,
    dump_file,
    final_data_file,
    params,
    equilibration="",
):
    header = _render_lammps_header(data_file, params)
    production = _render_adaptive_production(params)
    return f"""# Direct NEMD with adaptive heat-flux convergence
{header}

velocity        all create ${{T}} {params["random_seed"]} mom yes rot yes dist gaussian
{equilibration}
reset_timestep  0

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
velocity        bottomfix set 0.0 0.0 0.0
velocity        topfix set 0.0 0.0 0.0

compute         Thot all temp/region hot
compute         Tcold all temp/region cold
fix             fhot all langevin {params["nemd_hot_temperature"]} {params["nemd_hot_temperature"]} {params["nemd_thermostat_damp"]} {params["random_seed"]} tally yes
fix             fcold all langevin {params["nemd_cold_temperature"]} {params["nemd_cold_temperature"]} {params["nemd_thermostat_damp"]} {int(params["random_seed"]) + 104729} tally yes
fix_modify      fhot temp Thot
fix_modify      fcold temp Tcold
fix             hold_bottom bottomfix setforce 0.0 0.0 0.0
fix             hold_top topfix setforce 0.0 0.0 0.0
group           mobile subtract all bottomfix topfix
fix             nemd_nve mobile nve
thermo_style    custom step temp press density c_Thot c_Tcold
print           "NEMD steady state starts: {params["nemd_steady_steps"]} steps"
run             {params["nemd_steady_steps"]}

reset_timestep  0
compute         myKE all ke/atom
compute         myPE all pe/atom
compute         stressPair all stress/atom NULL pair kspace
compute         stressBond all centroid/stress/atom NULL bond angle dihedral improper
# real units: stress (atm A^3) / 68568.415 = energy (kcal/mol).
# rmask selects the current occupants, including atoms crossing slab boundaries.
variable        local_Jx atom "rmask(flux)*((c_myKE+c_myPE)*vx-((c_stressPair[1]+c_stressBond[1])*vx+(c_stressPair[4]+c_stressBond[4])*vy+(c_stressPair[5]+c_stressBond[5])*vz)/68568.415)"
compute         regional_flux all reduce sum v_local_Jx
variable        flux_volume equal ly*lz*v_flux_thick
variable        Jx equal c_regional_flux/v_flux_volume
variable        cross_section equal ly*lz
fix             average_flux all ave/time {params["nemd_profile_nevery"]} {params["nemd_profile_repeat"]} {params["nemd_profile_every"]} v_Jx ave one file "{transfer_file}"
fix             convergence_flux all ave/time {params["nemd_profile_nevery"]} {params["nemd_production_block_steps"] // params["nemd_profile_nevery"]} {params["nemd_production_block_steps"]} v_Jx ave one
fix             bath_energy all ave/time {params["nemd_profile_nevery"]} 1 {params["nemd_profile_every"]} f_fhot f_fcold v_cross_section file "{transfer_file.replace('_heat_flux.dat', '_bath_energy.dat')}"

compute         layers all chunk/atom bin/1d x lower 0.5 units box
fix             temperature_profile all ave/chunk {params["nemd_profile_nevery"]} {params["nemd_profile_repeat"]} {params["nemd_profile_every"]} layers temp norm sample file "{profile_file}"
dump            trajectory all custom {params["dump_every"]} "{dump_file}" id mol type q x y z vx vy vz
thermo          {params["nemd_profile_every"]}
thermo_style    custom step temp press density c_Thot c_Tcold v_Jx
thermo_modify   flush yes

{production}

write_data      "{final_data_file}"
undump          trajectory
unfix           convergence_flux
unfix           average_flux
unfix           bath_energy
unfix           temperature_profile
unfix           nemd_nve
unfix           fhot
unfix           fcold
unfix           hold_bottom
unfix           hold_top
if              "${{heat_flux_ok}} == 0" then "quit 2"
"""


def _render_direct_nemd_script(
    data_file, profile_file, transfer_file, dump_file, density_file, params,
):
    """Use the same validated measurement stage for combined and split jobs."""
    equilibration = "minimize 1.0e-4 1.0e-6 500 5000\n"
    equilibration += _render_density_equilibration(density_file, params)
    return _render_adaptive_nemd_script(
        data_file, profile_file, transfer_file, dump_file,
        transfer_file.replace("_heat_flux.dat", "_final.data"), params,
        equilibration=equilibration,
    )
