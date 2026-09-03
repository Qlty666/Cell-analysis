"""Stage runner for the standalone molecular docking board."""

from __future__ import annotations

import logging
import shutil
from datetime import datetime

from docking import analysis, docking, ligands, receptor, redock
from docking.config import ResolvedConfig
from docking.utils import DockingError

from .report import generate_report

STAGES = [
    ("01", "prepare-receptor", receptor.prepare_receptor),
    ("02", "prepare-ligands", ligands.prepare_ligands),
    ("03", "dock", docking.run_docking),
    ("04", "analyze", analysis.analyze_results),
    ("05", "redock", redock.run_redock),
    ("06", "report", generate_report),
]


def run_pipeline(
    cfg: ResolvedConfig,
    force: bool = False,
    start_stage: str | None = None,
) -> None:
    """Run receptor/ligand preparation, docking, analysis and reporting."""
    log = logging.getLogger("molecular_docking")
    stage_dir = cfg.stage_dir()
    if force and stage_dir.exists():
        resolved_out = cfg.output_dir.resolve()
        resolved_stage = stage_dir.resolve()
        if (
            resolved_stage.parent != resolved_out
            or resolved_stage.name != ".stages"
        ):
            raise DockingError(
                "refusing to remove an unexpected molecular docking stage dir"
            )
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
    log.info("molecular docking pipeline complete")
