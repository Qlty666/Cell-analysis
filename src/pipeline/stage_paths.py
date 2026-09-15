"""Workdir layout helpers shared by the full-pipeline stages."""

from __future__ import annotations

import json
import logging
from pathlib import Path

log = logging.getLogger("full_pipeline")


def _integration_dir(workdir: Path) -> Path:
    return workdir / "outputs" / "integration"


def _stage_dir(workdir: Path) -> Path:
    return _integration_dir(workdir) / ".stages"


def _marker(workdir: Path, code: str, name: str) -> Path:
    return _stage_dir(workdir) / f"{code}_{name}.done"


def _read_json(path: Path, default=None) -> dict:
    if not path.exists():
        return default or {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else default or {}
    except Exception:
        return default or {}
