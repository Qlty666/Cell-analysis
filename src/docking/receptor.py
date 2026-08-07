#!/usr/bin/env python3
"""Receptor preparation stage."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

from .config import ResolvedConfig
from .utils import (
    DockingError,
    ToolNotFoundError,
    find_script,
    find_tool,
    run_command,
    tail,
    tool_command,
)


def prepare_receptor(cfg: ResolvedConfig, log):
    inp = cfg.receptor_input()
    out = cfg.receptor_output()
    if not inp.exists():
        raise DockingError(f"receptor input not found: {inp}")
    out.parent.mkdir(parents=True, exist_ok=True)

    if inp.suffix.lower() == ".pdbqt":
        if inp.resolve() != out.resolve():
            shutil.copyfile(inp, out)
        log.info("receptor PDBQT ready: %s", out)
        return out

    ph = float(cfg.get("ligand", "ph", 7.4))
    script = find_script("mk_prepare_receptor.py")
    if script:
        try:
            _run_meeko_receptor(script, inp, out, ph)
            log.info("prepared receptor with Meeko: %s", out)
            return out
        except DockingError as exc:
            log.warning("Meeko receptor preparation failed (%s); trying fallbacks", exc)

    obabel = find_tool("obabel")
    if obabel:
        cmd = [obabel, str(inp), "-O", str(out), "-xr", "-p", str(ph)]
        _check(run_command(cmd, timeout=600), "obabel")
        log.info("prepared receptor with Open Babel: %s", out)
        return out

    adt = find_script("prepare_receptor4.py")
    if adt:
        cmd = [
            sys.executable,
            adt,
            "-r",
            str(inp),
            "-o",
            str(out),
            "-A",
            "checkhydrogens",
            "-U",
            "nphs_lps_waters",
        ]
        _check(run_command(cmd, timeout=600), "prepare_receptor4.py")
        log.info("prepared receptor with MGLTools: %s", out)
        return out

    raise ToolNotFoundError(
        "no receptor preparation tool found; install Meeko, Open Babel or MGLTools"
    )


def _run_meeko_receptor(script: str, inp: Path, out: Path, ph: float) -> None:
    cmd = tool_command(script) + ["-i", str(inp), "-o", str(out), "--ph", str(ph)]
    try:
        _check(run_command(cmd, timeout=600), "mk_prepare_receptor.py")
    except DockingError:
        cmd = tool_command(script) + ["-i", str(inp), "-o", str(out)]
        _check(run_command(cmd, timeout=600), "mk_prepare_receptor.py")
    if not out.exists() or out.stat().st_size == 0:
        raise DockingError("mk_prepare_receptor.py produced an empty output file")


def _check(result, label: str) -> None:
    if result.returncode != 0:
        raise DockingError(
            f"{label} failed (code {result.returncode}): "
            f"{tail(result.stderr or result.stdout)}"
        )
