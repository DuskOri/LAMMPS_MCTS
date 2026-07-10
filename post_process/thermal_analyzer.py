"""热流和温度剖面的后处理计算。"""

from dataclasses import dataclass
from pathlib import Path


REAL_TO_W_MK = 4184.0 / 6.02214076e23 / 1.0e-25


@dataclass
class ThermalConductivityResult:
    """记录一条候选体系的热导率估算结果。"""

    input_path: str
    flux_path: str
    temp_path: str
    heat_flux: float
    temperature_gradient: float
    conductivity_w_mk: float
    success: bool
    message: str = ""


def analyze_thermal_outputs(input_path, flux_path, temp_path, trim_fraction=0.20):
    """根据热流文件和温度剖面估算热导率。"""
    try:
        heat_flux = read_average_heat_flux(flux_path)
        gradient = read_temperature_gradient(temp_path, trim_fraction=trim_fraction)
    except OSError as error:
        return ThermalConductivityResult(
            input_path=str(input_path),
            flux_path=str(flux_path),
            temp_path=str(temp_path),
            heat_flux=0.0,
            temperature_gradient=0.0,
            conductivity_w_mk=0.0,
            success=False,
            message=str(error),
        )
    except ValueError as error:
        return ThermalConductivityResult(
            input_path=str(input_path),
            flux_path=str(flux_path),
            temp_path=str(temp_path),
            heat_flux=0.0,
            temperature_gradient=0.0,
            conductivity_w_mk=0.0,
            success=False,
            message=str(error),
        )

    if gradient == 0.0:
        conductivity = 0.0
        success = False
        message = "temperature gradient is zero"
    else:
        conductivity = abs(heat_flux / gradient) * REAL_TO_W_MK
        success = True
        message = "ok"

    return ThermalConductivityResult(
        input_path=str(input_path),
        flux_path=str(flux_path),
        temp_path=str(temp_path),
        heat_flux=heat_flux,
        temperature_gradient=gradient,
        conductivity_w_mk=conductivity,
        success=success,
        message=message,
    )


def read_average_heat_flux(filename):
    """读取 flux.profile，返回后半段平均热流。"""
    rows = _read_numeric_rows(filename, min_columns=3)
    if not rows:
        raise ValueError("flux profile has no numeric rows")

    start = len(rows) // 2
    values = [row[1] for row in rows[start:]]
    return sum(values) / len(values)


def read_temperature_gradient(filename, trim_fraction=0.20):
    """读取 temp.profile 最后一帧，拟合中间区域温度梯度。"""
    rows = _read_last_temp_block(filename)
    if len(rows) < 3:
        raise ValueError("temperature profile has too few rows")

    rows.sort(key=lambda row: row[0])
    trim_count = int(len(rows) * trim_fraction)
    if trim_count > 0 and len(rows) - 2 * trim_count >= 3:
        rows = rows[trim_count:-trim_count]

    xs = [row[0] for row in rows]
    temps = [row[1] for row in rows]
    return _linear_slope(xs, temps)


def _read_numeric_rows(filename, min_columns):
    """读取普通数字表，跳过注释和空行。"""
    rows = []
    for raw_line in Path(filename).read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) < min_columns:
            continue

        try:
            rows.append([float(part) for part in parts])
        except ValueError:
            continue
    return rows


def _read_last_temp_block(filename):
    """读取 ave/chunk 温度文件中的最后一个时间块。"""
    blocks = []
    current = []

    for raw_line in Path(filename).read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        if len(parts) == 3:
            if current:
                blocks.append(current)
                current = []
            continue

        if len(parts) < 4:
            continue

        try:
            coord = float(parts[1])
            temp = float(parts[3])
        except ValueError:
            continue
        current.append((coord, temp))

    if current:
        blocks.append(current)

    if not blocks:
        raise ValueError("temperature profile has no data block")
    return blocks[-1]


def _linear_slope(xs, ys):
    """最小二乘直线斜率。"""
    n = len(xs)
    x_mean = sum(xs) / n
    y_mean = sum(ys) / n
    denominator = sum((x - x_mean) ** 2 for x in xs)
    if denominator == 0.0:
        return 0.0
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    return numerator / denominator
