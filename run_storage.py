"""Atomic files and provenance for independently archived search runs."""
from contextlib import contextmanager
from dataclasses import asdict, is_dataclass
import hashlib
import json
import os
from pathlib import Path
import tempfile
from datetime import datetime
from uuid import uuid4


def new_run_id():
    return datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid4().hex[:12]


@contextmanager
def atomic_file(filename, newline=None):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline=newline) as stream:
            yield stream
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json(filename, data):
    with atomic_file(filename) as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False,
                  default=lambda item: asdict(item) if is_dataclass(item) else str(item))


def file_hash(filename):
    digest = hashlib.sha256()
    with Path(filename).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
