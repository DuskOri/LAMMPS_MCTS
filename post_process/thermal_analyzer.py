"""直接 NEMD 热导率输出读取。"""

from dataclasses import dataclass
from pathlib import Path
import math


@dataclass
class ThermalConductivityResult:
    """记录一次热导率分析结果。"""

    input_path: str
    correlation_path: str
    conductivity_path: str
    conductivity_w_mk: float
    success: bool
    message: str = ""
    method: str = "direct_nemd"
    gradient_k_per_a: float = 0.0
    gradient_r2: float = 0.0
    heat_flux_w_m2: float = 0.0


def analyze_direct_nemd_outputs(
    input_path,
    temperature_profile_path,
    heat_flux_path,
    tail_fraction=0.50,
    min_gradient_r2=0.80,
    min_temperature_span=5.0,
):
    """用中央区热流和温度梯度计算直接 NEMD 热导率。"""
    try:
        coordinates, temperatures = read_temperature_profile(
            temperature_profile_path,
            tail_fraction=tail_fraction,
        )
        gradient, gradient_r2, temperature_span = fit_temperature_gradient(
            coordinates,
            temperatures,
        )
        native_flux = read_average_heat_flux(
            heat_flux_path,
            tail_fraction=tail_fraction,
        )
    except (OSError, ValueError) as error:
        return ThermalConductivityResult(
            input_path=str(input_path),
            correlation_path=str(temperature_profile_path),
            conductivity_path=str(heat_flux_path),
            conductivity_w_mk=0.0,
            success=False,
            message=str(error),
            method="direct_nemd",
        )

    if temperature_span < float(min_temperature_span):
        message = (
            f"NEMD temperature span {temperature_span:.3f} K is below "
            f"{float(min_temperature_span):.3f} K"
        )
        return _failed_nemd_result(
            input_path, temperature_profile_path, heat_flux_path, message,
            gradient, gradient_r2,
        )
    if gradient_r2 < float(min_gradient_r2):
        message = (
            f"NEMD temperature-gradient R2 {gradient_r2:.4f} is below "
            f"{float(min_gradient_r2):.4f}"
        )
        return _failed_nemd_result(
            input_path, temperature_profile_path, heat_flux_path, message,
            gradient, gradient_r2,
        )
    if gradient <= 0.0:
        return _failed_nemd_result(
            input_path, temperature_profile_path, heat_flux_path,
            "NEMD temperature gradient is zero", gradient, gradient_r2,
        )

    # real 单位下 Jx 为 kcal/(mol A^2 fs)，温度梯度为 K/A。
    kcal_per_mol_to_joule = 4184.0 / 6.02214076e23
    heat_flux_w_m2 = abs(native_flux) * kcal_per_mol_to_joule / 1.0e-35
    conductivity = heat_flux_w_m2 / (gradient * 1.0e10)
    if not math.isfinite(conductivity) or conductivity <= 0.0:
        return _failed_nemd_result(
            input_path, temperature_profile_path, heat_flux_path,
            "NEMD conductivity is not a positive finite value",
            gradient, gradient_r2, heat_flux_w_m2,
        )

    return ThermalConductivityResult(
        input_path=str(input_path),
        correlation_path=str(temperature_profile_path),
        conductivity_path=str(heat_flux_path),
        conductivity_w_mk=conductivity,
        success=True,
        message="ok",
        method="direct_nemd",
        gradient_k_per_a=gradient,
        gradient_r2=gradient_r2,
        heat_flux_w_m2=heat_flux_w_m2,
    )


def read_temperature_profile(filename, tail_fraction=0.50):
    """读取 ave/chunk 完整数据块，并平均末段温度剖面。"""
    blocks = _read_chunk_blocks(filename)
    if not blocks:
        raise ValueError("NEMD temperature profile has no complete data blocks")

    fraction = min(max(float(tail_fraction), 0.0), 1.0)
    tail_count = max(1, int(len(blocks) * fraction)) if fraction else 1
    selected = blocks[-tail_count:]
    reference = selected[-1]
    averaged = []
    for row_index, reference_row in enumerate(reference):
        values = []
        for block in selected:
            if len(block) != len(reference):
                continue
            row = block[row_index]
            if (
                int(row[0]) == int(reference_row[0])
                and row[2] > 0.0
                and math.isfinite(row[-1])
            ):
                values.append(row[-1])
        if values:
            averaged.append((reference_row[1], sum(values) / len(values)))
    if len(averaged) < 4:
        raise ValueError("NEMD temperature profile has too few populated layers")
    return [item[0] for item in averaged], [item[1] for item in averaged]


def fit_temperature_gradient(coordinates, temperatures):
    """在冷热源之间的中央 60% 区域拟合温度梯度。"""
    if len(coordinates) != len(temperatures) or len(coordinates) < 4:
        raise ValueError("NEMD temperature profile is incomplete")
    x_min = min(coordinates)
    x_max = max(coordinates)
    span = x_max - x_min
    lower = x_min + 0.20 * span
    upper = x_min + 0.80 * span
    points = [
        (x, temperature)
        for x, temperature in zip(coordinates, temperatures)
        if lower <= x <= upper and math.isfinite(temperature)
    ]
    if len(points) < 4:
        raise ValueError("NEMD central fitting region has too few layers")
    slope, r_squared = _linear_regression(points)
    temperature_span = max(value for _, value in points) - min(
        value for _, value in points
    )
    return abs(slope), r_squared, temperature_span


def read_average_heat_flux(filename, tail_fraction=0.50):
    """读取 ave/time 输出，返回末段 Jx 的有符号平均值。"""
    rows = _read_numeric_rows(filename, min_columns=2)
    values = [row[-1] for row in rows if math.isfinite(row[-1])]
    if not values:
        raise ValueError("NEMD heat-flux output has no finite numeric rows")
    fraction = min(max(float(tail_fraction), 0.0), 1.0)
    tail_count = max(1, int(len(values) * fraction)) if fraction else 1
    return sum(values[-tail_count:]) / tail_count


def _read_chunk_blocks(filename):
    lines = Path(filename).read_text(encoding="utf-8", errors="ignore").splitlines()
    blocks = []
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        index += 1
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            chunk_count = int(float(parts[1]))
        except ValueError:
            continue
        block = []
        while index < len(lines) and len(block) < chunk_count:
            row_line = lines[index].strip()
            index += 1
            if not row_line or row_line.startswith("#"):
                continue
            try:
                row = [float(value) for value in row_line.split()]
            except ValueError:
                break
            if len(row) >= 4:
                block.append(row)
        if len(block) == chunk_count:
            blocks.append(block)
    return blocks


def _linear_regression(points):
    count = len(points)
    mean_x = sum(x for x, _ in points) / count
    mean_y = sum(y for _, y in points) / count
    denominator = sum((x - mean_x) ** 2 for x, _ in points)
    if denominator <= 0.0:
        raise ValueError("NEMD temperature coordinates do not span a distance")
    slope = sum((x - mean_x) * (y - mean_y) for x, y in points) / denominator
    intercept = mean_y - slope * mean_x
    residual = sum((y - (slope * x + intercept)) ** 2 for x, y in points)
    total = sum((y - mean_y) ** 2 for _, y in points)
    r_squared = 1.0 if total <= 0.0 else max(0.0, 1.0 - residual / total)
    return slope, r_squared


def _failed_nemd_result(
    input_path,
    profile_path,
    flux_path,
    message,
    gradient=0.0,
    gradient_r2=0.0,
    heat_flux_w_m2=0.0,
):
    return ThermalConductivityResult(
        input_path=str(input_path),
        correlation_path=str(profile_path),
        conductivity_path=str(flux_path),
        conductivity_w_mk=0.0,
        success=False,
        message=message,
        method="direct_nemd",
        gradient_k_per_a=gradient,
        gradient_r2=gradient_r2,
        heat_flux_w_m2=heat_flux_w_m2,
    )

def _read_numeric_rows(filename, min_columns):
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
