"""快速 Green-Kubo 热导率输出读取。"""

from dataclasses import dataclass
from pathlib import Path
import math


@dataclass
class ThermalConductivityResult:
    """记录 Green-Kubo 热导率分析结果。"""

    input_path: str
    correlation_path: str
    conductivity_path: str
    conductivity_w_mk: float
    success: bool
    message: str = ""
    method: str = "rapid_green_kubo"

def analyze_rapid_gk_outputs(
    input_path,
    correlation_path,
    conductivity_path,
    tail_fraction=0.20,
):
    """读取热导率运行值，并以末段平均作为快速筛选结果。"""
    try:
        conductivity = read_plateau_conductivity(
            conductivity_path,
            tail_fraction=tail_fraction,
        )
    except (OSError, ValueError) as error:
        return ThermalConductivityResult(
            input_path=str(input_path),
            correlation_path=str(correlation_path),
            conductivity_path=str(conductivity_path),
            conductivity_w_mk=0.0,
            success=False,
            message=str(error),
        )

    if not math.isfinite(conductivity) or conductivity <= 0.0:
        return ThermalConductivityResult(
            input_path=str(input_path),
            correlation_path=str(correlation_path),
            conductivity_path=str(conductivity_path),
            conductivity_w_mk=0.0,
            success=False,
            message="Green-Kubo conductivity is not a positive finite value",
        )

    return ThermalConductivityResult(
        input_path=str(input_path),
        correlation_path=str(correlation_path),
        conductivity_path=str(conductivity_path),
        conductivity_w_mk=conductivity,
        success=True,
        message="ok",
    )


def read_plateau_conductivity(filename, tail_fraction=0.20):
    """读取 fix ave/time 输出，返回末段有限正值的平均值。"""
    rows = _read_numeric_rows(filename, min_columns=2)
    values = [row[-1] for row in rows if math.isfinite(row[-1])]
    if not values:
        raise ValueError("Green-Kubo conductivity output has no finite numeric rows")

    fraction = min(max(float(tail_fraction), 0.0), 1.0)
    tail_count = max(1, int(len(values) * fraction)) if fraction else 1
    tail = values[-tail_count:]
    return sum(tail) / len(tail)


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
