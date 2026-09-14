"""Reproducibility manifest for the integrated scientific pipeline."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from importlib import metadata
from pathlib import Path
from typing import Any


PACKAGES = [
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "matplotlib",
    "rdkit",
    "torch",
    "openpyxl",
    "requests",
    "PoseBusters",
]


def _sha256(path: Path) -> str:
    try:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return "missing"


def _git_state(root: Path) -> dict[str, Any]:
    def run(*arguments: str) -> str:
        try:
            result = subprocess.run(
                ["git", *arguments],
                cwd=str(root),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=10,
            )
        except Exception:
            return ""
        return result.stdout.strip() if result.returncode == 0 else ""

    status = run("status", "--porcelain")
    return {
        "commit": run("rev-parse", "HEAD"),
        "short_commit": run("rev-parse", "--short", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(status),
    }


def _package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for package in PACKAGES:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def _file_record(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {"path": None, "sha256": None, "exists": False}
    path = Path(path)
    return {
        "path": str(path),
        "sha256": _sha256(path) if path.exists() else "missing",
        "exists": path.exists(),
    }


def write_reproducibility_manifest(
    workdir: Path,
    single_cell_root: Path,
    config_path: Path,
    docking_config_path: Path,
    extra_parameters: dict[str, Any] | None = None,
) -> Path:
    """Record code, data, dependency and parameter state for a pipeline run."""
    out_dir = workdir / "outputs" / "integration"
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs = {
        "full_pipeline_config": config_path,
        "docking_config": docking_config_path,
        "single_cell_summary": single_cell_root / "results" / "summary.json",
        "single_cell_complete": single_cell_root
        / "results"
        / "pipeline_complete.json",
        "candidate_universe": out_dir / "candidate_universe.csv",
        "expanded_candidate_universe": out_dir
        / "candidate_universe_evidence_expanded.csv",
        "key_genes": out_dir / "key_genes.csv",
        "gene_evidence": out_dir / "gene_evidence.csv",
        "target_priority": out_dir / "target_priority.csv",
        "integrated_target_priority": out_dir
        / "integrated_target_priority.csv",
        "evidence_hub_database": out_dir
        / "evidence_hub"
        / "evidence.sqlite",
    }
    outputs = {
        "integration_report": out_dir / "integration_report.html",
        "integration_summary": out_dir / "integration_summary.json",
        "run_manifest": out_dir / "run_manifest.json",
        "docking_targets": out_dir / "docking_targets.csv",
        "cadd_targets": out_dir / "cadd_targets.csv",
    }
    manifest = {
        "schema_version": "1.0",
        "stage": "11_report",
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "git": _git_state(Path(__file__).resolve().parents[2]),
        "packages": _package_versions(),
        "parameters": dict(extra_parameters or {}),
        "inputs": {
            name: _file_record(path) for name, path in inputs.items()
        },
        "outputs": {
            name: _file_record(path) for name, path in outputs.items()
        },
    }
    output = out_dir / "reproducibility_manifest.json"
    output.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return output
