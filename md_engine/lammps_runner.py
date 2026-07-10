"""LAMMPS 脚本运行接口。"""

from dataclasses import dataclass
from pathlib import Path
import subprocess


@dataclass
class LammpsRunResult:
    """记录一次 LAMMPS 运行的基本结果。"""

    input_path: str
    log_path: str
    success: bool
    returncode: int
    message: str = ""


def run_lammps_input(input_path, lammps_executable="lmp", log_path=None, timeout=None, work_dir=None):
    """执行一个 LAMMPS 输入脚本，并把标准输出保存成日志。"""
    input_file = Path(input_path)
    if log_path is None:
        log_file = input_file.with_suffix(".log")
    else:
        log_file = Path(log_path)

    log_file.parent.mkdir(parents=True, exist_ok=True)
    command = [str(lammps_executable), "-in", str(input_file)]

    try:
        result = subprocess.run(
            command,
            cwd=str(work_dir or Path.cwd()),
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
    except FileNotFoundError:
        return LammpsRunResult(
            input_path=str(input_file),
            log_path=str(log_file),
            success=False,
            returncode=127,
            message="LAMMPS executable not found.",
        )
    except subprocess.TimeoutExpired as error:
        log_file.write_text(error.stdout or "", encoding="utf-8", errors="ignore")
        return LammpsRunResult(
            input_path=str(input_file),
            log_path=str(log_file),
            success=False,
            returncode=124,
            message="LAMMPS run timeout.",
        )

    log_text = (result.stdout or "") + ("\n" + result.stderr if result.stderr else "")
    log_file.write_text(log_text, encoding="utf-8", errors="ignore")
    return LammpsRunResult(
        input_path=str(input_file),
        log_path=str(log_file),
        success=result.returncode == 0,
        returncode=result.returncode,
        message="ok" if result.returncode == 0 else _last_error_line(log_text),
    )


def _last_error_line(log_text):
    """提取日志末尾的错误信息。"""
    lines = [line.strip() for line in log_text.splitlines() if line.strip()]
    for line in reversed(lines):
        if "ERROR" in line.upper() or "WARNING" in line.upper():
            return line
    return lines[-1] if lines else "LAMMPS run failed."
