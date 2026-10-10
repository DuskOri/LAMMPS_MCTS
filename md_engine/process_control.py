"""Cancellation and child-process ownership shared by simulation stages."""
import threading
import subprocess

_CANCEL = threading.Event()
_LOCK = threading.Lock()
_PROCESSES = set()


class PipelineCancelled(RuntimeError):
    pass


def reset_cancellation():
    _CANCEL.clear()


def raise_if_cancelled():
    if _CANCEL.is_set():
        raise PipelineCancelled("run cancelled by user")


def register_process(process):
    with _LOCK:
        _PROCESSES.add(process)
    if _CANCEL.is_set():
        close_process(process)
        if process.stdout is not None:
            process.stdout.close()
        unregister_process(process)
        raise_if_cancelled()


def unregister_process(process):
    with _LOCK:
        _PROCESSES.discard(process)


def close_process(process):
    if process.poll() is None:
        try:
            process.terminate()
            process.wait(timeout=3)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
        except OSError:
            if process.poll() is None:
                raise


def cancel_processes():
    _CANCEL.set()
    with _LOCK:
        processes = list(_PROCESSES)
    for process in processes:
        close_process(process)
    return bool(processes)
