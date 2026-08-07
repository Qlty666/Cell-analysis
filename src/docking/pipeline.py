#!/usr/bin/env python3
"""Orchestrator with per-stage resume markers for the docking pipeline."""

from __future__ import annotations

import logging
import shutil
from datetime import datetime

from . import analysis, docking, ligands, receptor
from .config import ResolvedConfig
from .utils import DockingError

STAGES = [
    ("01", "prepare-receptor", receptor.prepare_receptor),
    ("02", "prepare-ligands", ligands.prepare_ligands),
    ("03", "dock", docking.run_docking),
    ("04", "analyze", analysis.analyze_results),
]


def run_pipeline(
    cfg: ResolvedConfig,
    force: bool = False,
    start_stage: str | None = None,
) -> None:
    log = logging.getLogger("docking")
    stage_dir = cfg.stage_dir()
    if force and stage_dir.exists():
        resolved_out = cfg.output_dir.resolve()
        resolved_stage = stage_dir.resolve()
        if (
            resolved_stage.parent != resolved_out
            or resolved_stage.name != ".stages"
        ):
            raise DockingError("refusing to remove an unexpected stage directory")
        shutil.rmtree(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)

    for code, name, fn in STAGES:
        marker = stage_dir / f"{code}_{name}.done"
        if not force and marker.exists():
            log.info("skip stage %s %s (already done)", code, name)
            continue
        if start_stage and code < start_stage:
            log.info("skip stage %s %s (start at %s)", code, name, start_stage)
            continue
        log.info("=== stage %s %s ===", code, name)
        fn(cfg, log)
        marker.write_text(datetime.now().isoformat(timespec="seconds"), encoding="utf-8")
        log.info("stage %s %s complete", code, name)
    log.info("pipeline complete")
