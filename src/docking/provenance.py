#!/usr/bin/env python3
"""Reproducibility manifests for pipeline runs.

Every scored or exported stage writes a ``run_manifest.json`` that records the
config hash, input file hashes, package versions and the exact parameters used.
This makes intermediate results traceable and re-runnable on another machine.
"""

from __future__ import annotations

import hashlib
import importlib
import platform
import subprocess
import sys
from pathlib import Path

from .config import ResolvedConfig
from .utils import now_iso, write_json


def sha256_file(path: Path, chunk_size: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as fh:
            while True:
                block = fh.read(chunk_size)
                if not block:
                    break
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return "unreadable"


def package_versions() -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for name in [
        "numpy",
        "pandas",
        "scipy",
        "sklearn",
        "torch",
        "joblib",
        "rdkit",
        "matplotlib",
    ]:
        try:
            module = importlib.import_module(name)
            out[name] = getattr(module, "__version__", None)
        except Exception:
            out[name] = None
    return out


def git_revision(root: Path) -> str | None:
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=str(root),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=10,
        )
    except Exception:
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout.strip() or None


def write_run_manifest(
    output_dir: Path,
    cfg: ResolvedConfig,
    stage: str,
    input_files: dict[str, Path | None] | None = None,
    extra: dict | None = None,
) -> Path:
    """Write a JSON manifest describing one reproducible pipeline stage."""
    inputs: dict[str, dict] = {}
    for name, path in (input_files or {}).items():
        if path is None:
            inputs[name] = {"path": None, "sha256": None}
            continue
        path = Path(path)
        inputs[name] = {
            "path": str(path),
            "sha256": sha256_file(path) if path.exists() else "missing",
        }

    manifest = {
        "stage": stage,
        "timestamp": now_iso(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "git_revision": git_revision(cfg.root),
        "config_path": str(cfg.config_path),
        "config_sha256": (
            sha256_file(cfg.config_path)
            if cfg.config_path.exists()
            else "missing"
        ),
        "packages": package_versions(),
        "parameters": extra or {},
        "inputs": inputs,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "run_manifest.json"
    write_json(manifest_path, manifest)
    return manifest_path
