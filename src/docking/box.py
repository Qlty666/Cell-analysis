#!/usr/bin/env python3
"""Detect a docking box from a cocrystal ligand or protein centroid."""

from __future__ import annotations

import math
from pathlib import Path

from .config import ResolvedConfig, save_config
from .utils import DockingError

WATER_NAMES = {"HOH", "WAT", "DOD"}


def detect_box_data(path: Path) -> tuple[list[float], list[float], str]:
    """Return (center, size, mode) for a receptor PDB/PDBQT file."""
    if not path.exists():
        raise DockingError(f"receptor file not found: {path}")
    atoms: list[list[float]] = []
    het: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM"):
            coord = _parse_coord(line)
            if coord:
                atoms.append(coord)
        elif line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname not in WATER_NAMES:
                coord = _parse_coord(line)
                if coord:
                    het.append(coord)
    if het:
        points = het
        mode = "cocrystal_ligand"
    elif atoms:
        points = atoms
        mode = "protein_centroid"
    else:
        raise DockingError("no ATOM/HETATM coordinates found in receptor")

    center = [round(sum(p[i] for p in points) / len(points), 2) for i in range(3)]
    size = _suggest_size(points, mode)
    return center, size, mode


def detect_and_update_config(cfg: ResolvedConfig, log):
    detect_value = cfg.get("receptor", "detect_input")
    detect_path = (
        Path(detect_value) if detect_value else cfg.receptor_input()
    )
    if not detect_path.is_absolute():
        detect_path = cfg.workdir / detect_path
    center, size, mode = detect_box_data(detect_path)
    cfg.data["receptor"]["center"] = center
    cfg.data["receptor"]["size"] = size
    save_config(cfg, cfg.config_path)
    log.info(
        "detect-box: mode=%s center=%s size=%s -> %s",
        mode,
        center,
        size,
        cfg.config_path,
    )
    return center, size, mode


def _parse_coord(line: str):
    try:
        return [
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        ]
    except ValueError:
        return None


def _suggest_size(points: list[list[float]], mode: str) -> list[float]:
    if mode == "protein_centroid":
        return [30.0, 30.0, 30.0]
    size = []
    for axis in range(3):
        span = max(p[axis] for p in points) - min(p[axis] for p in points)
        value = max(18.0, math.ceil(span) + 6.0)
        size.append(round(value * 2) / 2)
    return size
