"""快速热导率计算输入脚本生成模块。"""

from dataclasses import dataclass
from pathlib import Path


@dataclass
class ThermalConductivityInputResult:
    """记录一个 LAMMPS 热导率输入脚本的生成结果。"""

    data_path: str
    input_path: str
    flux_output: str
    temp_output: str
    dump_output: str
    success: bool
    message: str = ""


def write_fast_tc_input(data_file, input_file, output_prefix, params=None):
    """生成一个忽略分子间作用力的快速 NEMD 热导率计算脚本。"""
    params = _merge_params(_default_params(), params or {})
    input_path = Path(input_file)
    input_path.parent.mkdir(parents=True, exist_ok=True)

    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    flux_file = f"{output_prefix}_flux.profile"
    temp_file = f"{output_prefix}_temp.profile"
    dump_file = f"{output_prefix}_dump.lammpstrj"

    script = _render_fast_tc_script(
        data_file=Path(data_file).as_posix(),
        flux_file=Path(flux_file).as_posix(),
        temp_file=Path(temp_file).as_posix(),
        dump_file=Path(dump_file).as_posix(),
        params=params,
    )
    input_path.write_text(script, encoding="utf-8")

    return ThermalConductivityInputResult(
        data_path=str(data_file),
        input_path=str(input_path),
        flux_output=str(flux_file),
        temp_output=str(temp_file),
        dump_output=str(dump_file),
        success=True,
        message="ok",
    )


def _default_params():
    """默认参数偏向快速跑通流程，不追求收敛精度。"""
    return {
        "temperature": 300.0,
        "hot_temperature": 450.0,
        "cold_temperature": 150.0,
        "timestep": 0.5,
        "tdamp": 50.0,
        "cutoff": 12.0,
        "fix_fraction": 0.05,
        "source_fraction": 0.10,
        "flux_fraction": 0.10,
        "equil_steps": 5000,
        "gradient_steps": 20000,
        "sample_nevery": 10,
        "sample_nrepeat": 100,
        "sample_nfreq": 1000,
        "temp_bin_width": 1.0,
        "thermo_every": 1000,
        "dump_every": 5000,
        "random_seed": 87287,
        "langevin_hot_seed": 699483,
        "langevin_cold_seed": 399481,
    }


def _merge_params(default_params, user_params):
    """合并配置参数，保留默认值作为兜底。"""
    merged = dict(default_params)
    merged.update(user_params)
    return merged


def _render_fast_tc_script(data_file, flux_file, temp_file, dump_file, params):
    """按示例脚本的思路渲染快速热导率计算输入文件。"""
    return f"""# 快速热导率计算输入脚本
# 用于跑通 NEMD 热流和温度梯度统计流程，并忽略分子间非键作用力。
# 当前模式适合筛选和流程验证，正式计算需要补充完整力场参数。

units           real
dimension       3
atom_style      full
boundary        p p p

variable        T_equil   equal {params["temperature"]}
variable        T_hot     equal {params["hot_temperature"]}
variable        T_cold    equal {params["cold_temperature"]}
variable        dt        equal {params["timestep"]}
variable        tdamp     equal {params["tdamp"]}
variable        rc        equal {params["cutoff"]}

# 当前 data 文件只有拓扑和坐标，没有完整非键/键参数。
# pair_style zero 用来忽略所有 pair 力；neigh_modify 行保留“忽略分子间作用力”的显式设置。
pair_style      zero ${{rc}}
bond_style      harmonic

read_data       {data_file}

# 将所有 bond type 的力常数置零，使脚本可以在缺少真实力场时快速跑通。
pair_coeff      * *
bond_coeff      * 0.0 1.0
special_bonds   lj/coul 0.0 0.0 0.5
neigh_modify    exclude molecule/inter all

neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

timestep        ${{dt}}
thermo          {params["thermo_every"]}
thermo_style    custom step temp press vol etotal
thermo_modify   flush yes

velocity        all create ${{T_equil}} {params["random_seed"]} dist gaussian

# 先做短时间 NVE，让速度场稳定；快速模式不做长时间压缩和退火。
fix             eq all nve
run             {params["equil_steps"]}
unfix           eq

variable        xlo_now equal xlo
variable        xhi_now equal xhi
variable        ylo_now equal ylo
variable        yhi_now equal yhi
variable        zlo_now equal zlo
variable        zhi_now equal zhi
variable        Lx equal v_xhi_now-v_xlo_now
variable        Ly equal v_yhi_now-v_ylo_now
variable        Lz equal v_zhi_now-v_zlo_now

variable        fix_thick  equal {params["fix_fraction"]}*v_Lx
variable        src_thick  equal {params["source_fraction"]}*v_Lx
variable        flux_thick equal {params["flux_fraction"]}*v_Lx

variable        bottom_fix_bottom equal v_xlo_now
variable        bottom_fix_top    equal v_xlo_now+v_fix_thick
variable        hot_bottom        equal v_bottom_fix_top
variable        hot_top           equal v_hot_bottom+v_src_thick
variable        cold_top          equal v_xhi_now-v_fix_thick
variable        cold_bottom       equal v_cold_top-v_src_thick
variable        top_fix_bottom    equal v_cold_top
variable        top_fix_top       equal v_xhi_now

region          bottomfix block v_bottom_fix_bottom v_bottom_fix_top INF INF INF INF
region          topfix    block v_top_fix_bottom v_top_fix_top INF INF INF INF
region          hot       block v_hot_bottom v_hot_top INF INF INF INF
region          cold      block v_cold_bottom v_cold_top INF INF INF INF

group           bottomfix region bottomfix
group           topfix    region topfix
group           hot       region hot
group           cold      region cold
group           mobile    subtract all bottomfix topfix

variable        flux_center equal (v_xlo_now+v_xhi_now)/2.0
variable        flux_bottom equal v_flux_center-v_flux_thick/2.0
variable        flux_top    equal v_flux_center+v_flux_thick/2.0
region          flux block v_flux_bottom v_flux_top INF INF INF INF
group           flux region flux

compute         Thot all temp/region hot
compute         Tcold all temp/region cold

fix             freeze1 bottomfix setforce 0.0 0.0 0.0
fix             freeze2 topfix setforce 0.0 0.0 0.0
velocity        bottomfix set 0.0 0.0 0.0
velocity        topfix set 0.0 0.0 0.0

fix             int mobile nve
fix             fhot all langevin ${{T_hot}} ${{T_hot}} ${{tdamp}} {params["langevin_hot_seed"]} tally yes
fix             fcold all langevin ${{T_cold}} ${{T_cold}} ${{tdamp}} {params["langevin_cold_seed"]} tally yes
fix_modify      fhot temp Thot
fix_modify      fcold temp Tcold

# 用热源和冷源能量交换的平均值估算输入热流。
variable        area equal v_Ly*v_Lz
variable        elapsed_time equal step*v_dt
variable        Qin equal abs(f_fhot)
variable        Qout equal abs(f_fcold)
variable        Jsource equal 0.5*(v_Qin+v_Qout)/(v_area*v_elapsed_time)

compute         myKE flux ke/atom
compute         myPE flux pe/atom
compute         myStress flux stress/atom NULL virial
compute         flux_vector flux heat/flux myKE myPE myStress
variable        Volume_flux equal v_Ly*v_Lz*v_flux_thick
variable        Jx equal c_flux_vector[1]/v_Volume_flux

fix             ave_flux all ave/time {params["sample_nevery"]} {params["sample_nrepeat"]} {params["sample_nfreq"]} v_Jsource v_Jx file {flux_file}
compute         layers all chunk/atom bin/1d x lower {params["temp_bin_width"]} units box
fix             ave_temp all ave/chunk {params["sample_nevery"]} {params["sample_nrepeat"]} {params["sample_nfreq"]} layers temp norm sample file {temp_file}

dump            traj all custom {params["dump_every"]} {dump_file} id mol type x y z vx vy vz
thermo_style    custom step temp c_Thot c_Tcold v_Jsource v_Jx f_fhot f_fcold
run             {params["gradient_steps"]}

undump          traj
unfix           ave_temp
unfix           ave_flux
unfix           fhot
unfix           fcold
unfix           int
unfix           freeze1
unfix           freeze2
"""
