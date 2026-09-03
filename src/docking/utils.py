#!/usr/bin/env python3
"""Shared helpers: logging, subprocess, tool discovery and Vina parsing."""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


class DockingError(Exception):
    """Base error for the docking pipeline."""


class ToolNotFoundError(DockingError):
    """Raised when an external tool required by a stage is missing."""


def setup_logging(verbose: bool = False) -> logging.Logger:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    return logging.getLogger("docking")


def run_command(
    cmd: list[str],
    timeout: int = 300,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    stdin_text: str | None = None,
) -> subprocess.CompletedProcess:
    logger = logging.getLogger("docking")
    logger.debug("$ %s", " ".join(str(part) for part in cmd))
    try:
        return subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            env=env,
            input=stdin_text,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
    except FileNotFoundError as exc:
        raise ToolNotFoundError(f"executable not found: {cmd[0]}") from exc
    except subprocess.TimeoutExpired as exc:
        raise DockingError(f"command timed out after {timeout}s: {cmd[0]}") from exc


def tail(text: str, limit: int = 3000) -> str:
    text = (text or "").strip()
    return text[-limit:] if len(text) > limit else text


def find_tool(name: str) -> str | None:
    direct = Path(name)
    if direct.is_absolute() and direct.exists():
        return str(direct)
    hit = shutil.which(name)
    if hit:
        return hit
    variants = [name]
    if not Path(name).suffix:
        variants += [name + ".exe", name + ".bat", name + ".cmd"]
    project_root = Path(__file__).resolve().parent.parent.parent
    candidates: list[Path] = []
    for variant in variants:
        candidates.append(Path(sys.executable).parent / variant)
        candidates.append(Path(sys.executable).parent.parent / "Scripts" / variant)
        candidates.append(project_root / "dock" / "tools" / variant)
        conda_bases: list[Path] = []
        for env_key in ["CONDA_PREFIX", "CONDA_ROOT", "CONDA_HOME", "MAMBA_ROOT_PREFIX"]:
            value = os.environ.get(env_key)
            if value:
                conda_bases.append(Path(value))
        conda_bases += [Path.home() / "anaconda3", Path.home() / "miniconda3"]
        for base in conda_bases:
            if not base.exists():
                continue
            for env in [
                "",
                "envs/vscreen",
                "envs/dock",
                "envs/StudyDnn",
                "envs/base",
            ]:
                for sub in ["Scripts", "bin"]:
                    candidates.append(base / env / sub / variant)
    for cand in candidates:
        if cand.exists():
            return str(cand)
    return None


def find_script(name: str) -> str | None:
    hit = find_tool(name)
    if hit:
        return hit
    stem = name[:-3] if name.lower().endswith(".py") else name
    hit = find_tool(stem)
    if hit:
        return hit
    try:
        import meeko  # type: ignore

        pkg = Path(meeko.__file__).resolve().parent
        for cand in [
            pkg.parent / "Scripts" / name,
            pkg.parent / "bin" / name,
            pkg / name,
        ]:
            if cand.exists():
                return str(cand)
    except Exception:
        pass
    return None


def tool_command(path: str) -> list[str]:
    """Return a command prefix that runs a Python script or native executable."""
    if path.lower().endswith(".py"):
        return [sys.executable, path]
    return [path]


def parse_vina_affinities(text: str) -> list[dict]:
    """Parse the affinity table printed by AutoDock Vina."""
    modes: list[dict] = []
    for line in (text or "").splitlines():
        parts = line.strip().split()
        if not parts or not parts[0].isdigit():
            continue
        try:
            mode = int(parts[0])
            affinity = float(parts[1])
        except (ValueError, IndexError):
            continue
        entry = {"mode": mode, "affinity": affinity}
        if len(parts) >= 4:
            try:
                entry["rmsd_lb"] = float(parts[2])
                entry["rmsd_ub"] = float(parts[3])
            except ValueError:
                pass
        modes.append(entry)
    return modes


def safe_name(value: str | None, fallback: str = "ligand") -> str:
    if value is None:
        value = fallback
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("._")
    return (cleaned or fallback)[:80]


def read_json(path: Path, default=None):
    if not path.exists():
        return default
    try:
        import json

        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return default


def write_json(path: Path, data) -> None:
    import json

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def read_text(path: Path, limit: int | None = None) -> str:
    if not path.exists():
        return ""
    text = path.read_text(encoding="utf-8", errors="replace")
    return text if limit is None else text[-limit:]


def human_seconds(seconds: float) -> str:
    return f"{seconds:.1f}s" if seconds < 60 else f"{seconds / 60:.1f}min"


def now_iso() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
