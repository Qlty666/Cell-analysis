"""Shared environment discovery helpers."""

from __future__ import annotations

import shutil
from pathlib import Path

R_BASES = (
    Path(r"C:\Program Files\R"),
    Path(r"C:\Program Files\Microsoft\R Open"),
    Path.home() / "AppData" / "Local" / "Programs" / "R",
)


def find_rscript() -> str | None:
    """Return the first available Rscript path, or None."""
    for name in ("Rscript", "Rscript.exe"):
        found = shutil.which(name)
        if found:
            return found
    for base in R_BASES:
        if not base.exists():
            continue
        candidates = sorted(base.glob("R-*/bin/Rscript.exe"), reverse=True)
        if candidates:
            return str(candidates[0])
    return None


def require_rscript(message: str = "Rscript not found") -> str:
    """Return an Rscript path or raise RuntimeError."""
    found = find_rscript()
    if found is None:
        raise RuntimeError(message)
    return found
