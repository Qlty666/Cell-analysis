"""Stage runner for the standalone molecular docking board."""

from __future__ import annotations

from docking import analysis, docking, ligands, receptor, redock
from docking.config import ResolvedConfig
from docking.pipeline import run_pipeline as _run_pipeline

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
    _run_pipeline(
        cfg,
        force=force,
        start_stage=start_stage,
        stages=STAGES,
        logger_name="molecular_docking",
        complete_message="molecular docking pipeline complete",
        stage_error=(
            "refusing to remove an unexpected molecular docking stage dir"
        ),
    )
