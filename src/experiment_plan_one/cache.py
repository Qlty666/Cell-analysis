"""Content-addressed local caches; unsigned legacy outputs are never reused."""
from __future__ import annotations

import platform
import sys
from importlib import metadata
from pathlib import Path
from typing import Any

from common.fingerprints import atomic_json, file_hash, fingerprint, read_state


def signature(*, files=(), parameters: Any = None) -> str:
    software_versions: dict[str, str] = {
        "python": sys.version,
        "platform": platform.platform(),
    }
    for name in ("numpy", "pandas", "scipy", "scikit-learn", "scanpy", "anndata"):
        try:
            software_versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            software_versions[name] = "not-installed"
    return fingerprint({
        "files": {str(Path(p).resolve()): file_hash(p) for p in files},
        "parameters": parameters,
        "software": software_versions,
    })


def valid(path: Path, key: str, outputs=()) -> bool:
    state = read_state(path)
    return (state.get("signature") == key and bool(outputs)
            and all(Path(p).is_file() and Path(p).stat().st_size > 0
                    and state.get("outputs", {}).get(str(Path(p).resolve()))
                    == file_hash(p) for p in outputs))


def save(path: Path, key: str, outputs) -> None:
    atomic_json(path, {"signature": key, "outputs": {
        str(Path(p).resolve()): file_hash(p) for p in outputs
    }})
