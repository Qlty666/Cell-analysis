#!/usr/bin/env python3
"""Orchestrator with per-stage resume markers for the docking pipeline."""

from __future__ import annotations

import logging
import copy
import shutil
from datetime import datetime
from pathlib import Path

from . import analysis, docking, ligands, receptor, redock, report
from .config import ResolvedConfig
from .utils import DockingError
from common.fingerprints import file_hash, fingerprint, read_state, atomic_json

STAGES = [
    ("01", "prepare-receptor", receptor.prepare_receptor),
    ("02", "prepare-ligands", ligands.prepare_ligands),
    ("03", "dock", docking.run_docking),
    ("04", "analyze", analysis.analyze_results),
    ("05", "redock", redock.run_redock),
    ("06", "report", report.generate_report),
]


def run_pipeline(
    cfg: ResolvedConfig,
    force: bool = False,
    start_stage: str | None = None,
    stages: list | None = None,
    logger_name: str = "docking",
    complete_message: str = "pipeline complete",
    stage_error: str = "refusing to remove an unexpected stage directory",
) -> None:
    """Run the configured stages with per-stage resume markers.

    ``stages``/``logger_name``/``complete_message`` let the standalone
    molecular docking board reuse this runner with its own stage list.
    """
    log = logging.getLogger(logger_name)
    cfg = copy.deepcopy(cfg)
    if force:
        cfg.data.setdefault("docking", {})["resume"] = False
        cfg.data.setdefault("redock", {})["resume"] = False
    stage_dir = cfg.stage_dir()
    if force and stage_dir.exists():
        resolved_out = cfg.output_dir.resolve()
        resolved_stage = stage_dir.resolve()
        if (
            resolved_stage.parent != resolved_out
            or resolved_stage.name != ".stages"
        ):
            raise DockingError(stage_error)
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    changed = False
    for code, name, fn in stages or STAGES:
        marker = stage_dir / f"{code}_{name}.done"
        signature = _stage_fingerprint(cfg, code)
        previous = read_state(marker)
        outputs = previous.get("outputs", {})
        valid = bool(outputs) and all(file_hash(p) == h and h != "missing" for p, h in outputs.items())
        if not force and not changed and previous.get("signature") == signature and valid:
            log.info("skip stage %s %s (already done)", code, name)
            continue
        if start_stage and code < start_stage:
            log.info("skip stage %s %s (start at %s)", code, name, start_stage)
            continue
        log.info("=== stage %s %s ===", code, name)
        marker.unlink(missing_ok=True)
        for later in stage_dir.glob("*.done"):
            if later.name[:2] > code:
                later.unlink()
        result = fn(cfg, log)
        paths = _stage_outputs(cfg, code, result)
        if not paths or any(not p.is_file() or p.stat().st_size == 0 for p in paths):
            raise DockingError(f"stage {name} did not produce valid required outputs")
        atomic_json(marker, {"signature": _stage_fingerprint(cfg, code), "outputs": {str(p): file_hash(p) for p in paths}, "completed": datetime.now().isoformat()})
        changed = True
        log.info("stage %s %s complete", code, name)
    log.info("%s", complete_message)


def _stage_fingerprint(cfg, code):
    sections = {"01": ["receptor"], "02": ["ligand"], "03": ["receptor", "docking"], "04": ["analysis"], "05": ["receptor", "docking", "redock"], "06": list(cfg.data)}
    inputs = {"01": [cfg.receptor_input()], "02": [cfg.ligand_input()], "03": [cfg.receptor_output(), cfg.manifest_path()], "04": [cfg.results_path()], "05": [cfg.receptor_output(), cfg.manifest_path(), cfg.analysis_dir() / "data/fig_46_47_ranked_results.csv"], "06": [cfg.results_path(), cfg.analysis_dir() / "summary.json", cfg.redock_dir() / "summary.json"]}.get(code, [])
    if code in {"03", "05"}:
        inputs += [Path(r["pdbqt"]) for r in docking._read_csv(cfg.manifest_path()) if r.get("pdbqt")]
        executable = shutil.which(str(cfg.get("docking", "executable", "vina"))) or str(cfg.get("docking", "executable", "vina"))
        inputs += [Path(executable)] + cfg.receptor_flexible()
    config = {key: {k: v for k, v in cfg.data.get(key, {}).items() if k != "resume"} if isinstance(cfg.data.get(key), dict) else cfg.data.get(key) for key in sections.get(code, [])}
    return fingerprint({"config": config, "inputs": {str(p): file_hash(p) for p in inputs}, "code": {str(p): file_hash(p) for p in Path(__file__).parent.glob("*.py")}})


def _stage_outputs(cfg, code, result):
    if code == "01":
        return [cfg.receptor_output()]
    if code == "02":
        return [cfg.manifest_path()] + [Path(r["pdbqt"]) for r in docking._read_csv(cfg.manifest_path()) if r.get("status") == "ok" and r.get("pdbqt")]
    if code in {"03", "05"}:
        table = cfg.results_path() if code == "03" else cfg.redock_dir() / "data/fig_49_redock_results.csv"
        rows = docking._read_csv(table)
        # Partial failure must remain retryable instead of receiving a done marker.
        if not rows or any(r.get("status") != "ok" for r in rows):
            raise DockingError(f"stage {code} has incomplete docking tasks; rerun to retry")
        return [table] + [Path(r.get("pose_file") or "") for r in rows]
    if code == "04":
        return [cfg.analysis_dir() / "summary.json", cfg.analysis_dir() / "data/fig_46_47_ranked_results.csv"]
    if isinstance(result, (str, Path)):
        return [Path(result)]
    return [p for p in (cfg.reports_dir() / "docking_report.html", cfg.reports_dir() / "molecular_docking_report.html") if p.is_file()]
