"""LAMMPS 脚本运行接口。"""

from dataclasses import dataclass
from pathlib import Path
import subprocess
import threading


@dataclass
class LammpsRunResult:
    """记录一次 LAMMPS 运行的基本结果。"""

    input_path: str
    log_path: str
    success: bool
    returncode: int
    message: str = ""


def run_lammps_input(
    input_path,
    lammps_executable="lmp",
    log_path=None,
    timeout=None,
    work_dir=None,
    stream_output=True,
):
    """执行 LAMMPS，持续写日志并把关键进度实时输出到命令行。"""
    input_file = Path(input_path)
    if log_path is None:
        log_file = input_file.with_suffix(".log")
    else:
        log_file = Path(log_path)

    log_file.parent.mkdir(parents=True, exist_ok=True)
    command = [str(lammps_executable), "-in", str(input_file)]
    if stream_output:
        print(f"[LAMMPS] input: {input_file.resolve()}", flush=True)
        print(f"[LAMMPS] log: {log_file.resolve()}", flush=True)

    try:
        process = subprocess.Popen(
            command,
            cwd=str(work_dir or Path.cwd()),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
    except FileNotFoundError:
        return LammpsRunResult(
            input_path=str(input_file),
            log_path=str(log_file),
            success=False,
            returncode=127,
            message="LAMMPS executable not found.",
        )

    timed_out = threading.Event()

    def stop_after_timeout():
        if process.poll() is None:
            timed_out.set()
            process.kill()

    timer = None
    if timeout is not None and float(timeout) > 0:
        timer = threading.Timer(float(timeout), stop_after_timeout)
        timer.daemon = True
        timer.start()

    lines = []
    try:
        with log_file.open("w", encoding="utf-8", errors="replace") as log_obj:
            for line in iter(process.stdout.readline, ""):
                lines.append(line)
                log_obj.write(line)
                log_obj.flush()
                if stream_output and _should_stream_line(line):
                    print(f"[LAMMPS] {line.rstrip()}", flush=True)
        process.stdout.close()
        returncode = process.wait()
    finally:
        if timer is not None:
            timer.cancel()

    log_text = "".join(lines)
    if timed_out.is_set():
        return LammpsRunResult(
            input_path=str(input_file),
            log_path=str(log_file),
            success=False,
            returncode=124,
            message="LAMMPS run timeout.",
        )

    return LammpsRunResult(
        input_path=str(input_file),
        log_path=str(log_file),
        success=returncode == 0,
        returncode=returncode,
        message="ok" if returncode == 0 else _last_error_line(log_text),
    )


def _should_stream_line(line):
    """筛选适合实时展示的 LAMMPS 进度行。"""
    text = line.strip()
    if not text:
        return False

    upper = text.upper()
    if upper.startswith(("ERROR", "WARNING")):
        return True
    if text.startswith(("LAMMPS (", "Reading data file", "Loop time of", "Total wall time")):
        return True
    if text.startswith(("minimize", "run ", "fix             eq_", "fix             production")):
        return True
    if text.startswith("Step"):
        return True

    columns = text.split()
    if len(columns) < 2:
        return False
    try:
        int(columns[0])
        float(columns[1])
    except ValueError:
        return False
    return True


def _last_error_line(log_text):
    """提取日志末尾的错误信息。"""
    lines = [line.strip() for line in log_text.splitlines() if line.strip()]
    for line in reversed(lines):
        if "ERROR" in line.upper() or "WARNING" in line.upper():
            return line
    return lines[-1] if lines else "LAMMPS run failed."
