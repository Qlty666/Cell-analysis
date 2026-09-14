#!/usr/bin/env python3
"""End-to-end automation from single-cell analysis to virtual screening.

The module wires the existing pieces together:

1. run the GEO single-cell pipeline and export a sample-level pseudobulk matrix;
2. rank significant DEGs into a compact key-gene table;
3. enrich genes with UniProt/PDB/ChEMBL/STRING/Reactome/Open Targets/KEGG evidence (network optional, cached);
4. build the virtual-knockout inputs and run multidimensional target scoring;
5. for genes with a PDB structure, collect known ligands and run the full
   AutoDock Vina pipeline in an isolated per-target workdir;
6. prepare GROMACS MD inputs, rescore with ML when trained data exists and
   export MD/external tool handoffs for each successful target;
7. run network toxicology and FAERS screening when the user provides inputs;
8. close the loop through cell feedback and write the integrated report.

Every stage writes a marker file so a rerun resumes where it stopped.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent.parent
SRC = APP_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from docking import (  # noqa: E402
    box,
    evidence as evidence_mod,
    handoff,
    md_simulation,
    ml as docking_ml,
    network_toxicology,
    pipeline as docking_pipeline,
    signal_detection,
)
from docking.config import load_config, save_config  # noqa: E402
from docking.knockout import run_knockout  # noqa: E402
from docking.provenance import write_run_manifest  # noqa: E402
from docking.utils import DockingError, ToolNotFoundError, safe_name, write_json  # noqa: E402
from docking.validation import export_validation  # noqa: E402
from common.html_utils import esc as _esc  # noqa: E402
from data.geo_downloader import canonical_accession  # noqa: E402
from evidence import (  # noqa: E402
    EvidenceContext,
    EvidenceHub,
    EvidenceRecord,
    EvidenceTier,
    SQLiteEvidenceStore,
)

from . import cell_feedback, orchestrator  # noqa: E402
# The private helpers remain re-exported for existing callers.
from .differential import (  # noqa: E402
    _betacf,
    _bh_adjust,
    _chi2_contingency,
    _regularized_incomplete_beta,
    _welch_ttest,
    run_differential_abundance,
)
from .errors import IntegrationError, PauseRequested  # noqa: E402
from .integrated_report import generate_integrated_report  # noqa: E402
from .key_targets import (  # noqa: E402
    DEFAULT_GENE_BLACKLIST,
    extract_key_genes,
)
from .qc import (  # noqa: E402
    _as_float,
    _min_affinity_text,
    collect_qc_metrics,
    evaluate_qc_gate,
    write_qc_metrics,
)
from .target_priority import (  # noqa: E402
    build_target_priority,
    write_target_priority_summary,
)
from .structural_quality import build_structural_quality  # noqa: E402
from .omics_qc import run_omics_qc  # noqa: E402
from .stage_paths import (  # noqa: E402
    _integration_dir,
    _marker,
    _read_json,
    _stage_dir,
)

log = logging.getLogger("full_pipeline")

STAGES = [
    ("01", "single_cell", "expression analysis (download, QC, annotation, DEG)"),
    ("02", "key_targets", "build the DEG candidate universe and compact key-gene table"),
    ("03", "evidence", "collect per-target structural/ligand evidence and rank targets with the multi-source evidence hub"),
    ("04", "knockout_inputs", "build pseudobulk expression and knockout inputs"),
    ("05", "knockout", "heuristic perturbation scoring and integrated target prioritization"),
    ("06", "docking", "per-target virtual screening with AutoDock Vina"),
    ("07", "cadd_downstream", "MD preparation, ML rescoring and MD/external handoff"),
    ("08", "network", "compound-disease network toxicology on optional user evidence"),
    ("09", "faers", "FAERS-style disproportionality signal screening on optional event table"),
    ("10", "cell_feedback", "re-score single-cell targets from knockout/docking results"),
    ("11", "report", "integrated HTML report and provenance manifest"),
]

# Required outputs per stage. Paths are relative to the full-pipeline workdir.
STAGE_OUTPUTS = {
    "01": ("results/pipeline_complete.json",),
    "02": (
        "outputs/integration/candidate_universe.csv",
        "outputs/integration/key_genes.csv",
        "outputs/integration/key_genes_summary.json",
    ),
    "03": (
        "outputs/integration/gene_evidence.csv",
        "outputs/integration/candidate_universe_evidence_expanded.csv",
        "outputs/integration/target_priority.csv",
        "outputs/integration/target_priority_summary.json",
    ),
    "04": (
        "data/knockout/expression.csv",
        "data/knockout/metadata.csv",
        "data/knockout/inputs_summary.json",
        "outputs/integration/omics_qc_summary.json",
        "outputs/integration/omics_qc_sample_metrics.csv",
    ),
    "05": (
        "outputs/integration/knockout_summary.json",
        "outputs/integration/integrated_target_priority.csv",
        "outputs/run_001/results/04_knockout/data/fig_52_53_ranked_knockout.csv",
        "outputs/run_001/results/05_validation/data/validation_plan.md",
    ),
    "06": (
        "outputs/integration/docking_summary.json",
        "outputs/integration/docking_targets.csv",
    ),
    "07": (
        "outputs/integration/cadd_downstream_summary.json",
        "outputs/integration/cadd_targets.csv",
        "outputs/integration/structural_quality_summary.json",
        "outputs/integration/structural_quality_targets.csv",
        "outputs/integration/pose_qc_summary.json",
        "outputs/integration/pose_qc_results.csv",
    ),
    "08": ("outputs/integration/network_summary.json",),
    "09": ("outputs/integration/faers_summary.json",),
    "10": ("outputs/integration/cell_feedback/cell_feedback_summary.json",),
    "11": (
        "outputs/integration/integration_report.html",
        "outputs/integration/integration_summary.json",
        "outputs/integration/target_validation_scores.csv",
        "outputs/integration/target_validation_summary.json",
        "outputs/integration/external_validation_summary.json",
        "outputs/integration/reproducibility_manifest.json",
        "outputs/integration/run_manifest.json",
    ),
}

DEFAULT_QC_GATE = {
    "enabled": True,
    "min_cells_after_qc": 0,
    "min_genes": 0,
    "max_doublet_rate": None,
    "min_deg_genes": 0,
    "require_pseudobulk": False,
    "fail_on_missing_metrics": False,
}

DEFAULT_DIFFERENTIAL_ABUNDANCE = {
    "enabled": True,
    "min_cells": 5,
    "fdr": 0.05,
}

DEFAULT_MD_SIMULATION = {
    "enabled": True,
    "mode": "prepare",
    "top_n": 1,
}

DEFAULT_HANDOFF = {
    "enabled": True,
}

DEFAULT_DOCKING_ML = {
    "enabled": True,
    "model": "rf",
    "training_csv": None,
    "label_column": "active",
}

DEFAULT_EVIDENCE = {
    "fetch": True,
    "max_workers": 6,
    "timeout": 90,
    "hub_enabled": True,
    "hub_config": "config/evidence_sources.json",
    "disease_name": "liver cancer",
    "disease_id": "",
    "max_targets": 300,
    "max_records_per_source": 1000,
    "hub_timeout": 120,
    "allow_network": True,
    "strict": False,
    "legacy_pool_size": 50,
    "candidate_expansion_max_targets": 1000,
    "benchmark_positive_targets": "",
    "benchmark_negative_targets": "",
    "benchmark_top_n": 20,
}

DEFAULT_TARGET_PRIORITY = {
    "weights": {
        "expression": 0.25,
        "evidence": 0.45,
        "knockout": 0.20,
        "advanced": 0.10,
    },
    "go_min_score": 0.75,
    "go_min_coverage": 0.50,
    "conditional_min_score": 0.50,
}

DEFAULT_DOCKING_SELECTION = {
    "allow_review": False,
}

DEFAULT_EXTERNAL_VALIDATION = {
    "enabled": False,
    "path": None,
    "target_column": "gene",
    "score_column": "score",
    "label_column": "label",
    "threshold": 0.5,
    "bootstrap": 1000,
}

DEFAULT_NETWORK_TOXICOLOGY = {
    "enabled": True,
    "compound_name": None,
    "disease_name": None,
    "compound_targets_csv": None,
    "target_sources": None,
    "target_sources_dir": None,
    "disease_genes_csv": None,
    "disease_gene_column": None,
    "ppi_network_csv": None,
    "venn": True,
    "output_dir": "outputs/run_001/network_toxicology",
    "cytoscape": "auto",
    "cytoscape_url": "http://127.0.0.1:1234",
    "cytoscape_layout": "cose",
    "cytoscape_save_session": False,
    "max_ppi_edges": 2000,
    "run_enrichment": False,
    "enrichment_timeout": 900,
}

DEFAULT_FAERS = {
    "enabled": True,
    "input_csv": None,
    "drug_column": "drug",
    "event_column": "event",
    "count_column": None,
    "min_count": 3,
    "output_dir": "outputs/run_001/faers",
}

EVIDENCE_COLUMNS = [
    "gene",
    "entrez",
    "uniprot",
    "ensembl",
    "chembl_target_id",
    "known_ligands",
    "chembl_bioactivities",
    "pdb_structures",
    "pdb_ids",
    "off_target_paralogs",
    "safety_concern",
    "string_partners",
    "string_partner_ids",
    "reactome_pathways",
    "reactome_pathway_ids",
    "pharmgkb_annotations",
    "pharmgkb_ids",
    "alphafold_structures",
    "alphafold_ids",
    "opentargets_hits",
    "opentargets_target_ids",
    "kegg_pathways",
    "kegg_pathway_ids",
    "database_sources",
]








def _sha256_file(path: Path) -> str:
    try:
        digest = hashlib.sha256()
        with path.open("rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError as exc:
        log.warning("could not hash %s: %s", path, exc)
        return "missing"


def _json_sorted(value) -> str:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            ensure_ascii=False,
            default=str,
        )
    except Exception as exc:
        log.debug(
            "JSON serialization fell back to str() for %s: %s",
            type(value).__name__,
            exc,
        )
        return str(value)


def _read_stage_marker(workdir: Path, code: str, name: str) -> dict | None:
    marker = _marker(workdir, code, name)
    if not marker.exists():
        return None
    try:
        data = json.loads(marker.read_text(encoding="utf-8"))
        if isinstance(data, dict) and data.get("signature"):
            return data
    except (OSError, ValueError) as exc:
        log.warning("unreadable stage marker %s: %s", marker, exc)
    return None


def _write_stage_marker(
    workdir: Path,
    code: str,
    name: str,
    signature: str,
    note: str = "",
) -> None:
    marker = _marker(workdir, code, name)
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(
        json.dumps(
            {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "stage": f"{code}_{name}",
                "signature": signature,
                "note": note,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )


def _stage_signature(code: str, args, workdir: Path, ctx: dict) -> str:
    """Fingerprint the parameters and inputs that can invalidate a stage."""
    integration = _integration_dir(workdir)
    config_path = Path(str(getattr(args, "config", "") or "")).resolve()
    docking_config = Path(
        str(ctx.get("docking_config") or getattr(args, "docking_config", "") or "")
    ).resolve()
    payload: dict = {"stage": code}

    if code == "01":
        payload.update(
            {
                "single_cell_root": str(ctx.get("single_cell_root") or ""),
                "accession": getattr(args, "accession", None),
                "species": getattr(args, "species", None),
                "skip_scrna": bool(getattr(args, "skip_scrna", False)),
                "skip_download": bool(getattr(args, "skip_download", False)),
                "skip_deps": bool(getattr(args, "skip_deps", False)),
                "qc_gate": _json_sorted(getattr(args, "qc_gate", {})),
                "differential_abundance": _json_sorted(
                    getattr(args, "differential_abundance", {})
                ),
            }
        )
    elif code == "02":
        payload.update(
            {
                "single_cell_root": str(ctx.get("single_cell_root") or ""),
                "top_genes": int(getattr(args, "top_genes", 50) or 50),
                "candidate_universe_size": int(
                    getattr(args, "candidate_universe_size", 1000) or 0
                ),
                "keep_all_genes": bool(getattr(args, "keep_all_genes", False)),
                "gene_blacklist": sorted(
                    getattr(args, "gene_blacklist", DEFAULT_GENE_BLACKLIST) or []
                ),
            }
        )
    elif code == "03":
        key_genes = ctx.get("key_genes_path") or integration / "key_genes.csv"
        candidate_universe = (
            ctx.get("candidate_universe_path")
            or integration / "candidate_universe.csv"
        )
        hub_config = Path(
            str(getattr(args, "evidence_hub_config", "") or "")
        ).resolve()
        payload.update(
            {
                "key_genes_csv": str(key_genes),
                "key_genes_sha256": _sha256_file(Path(str(key_genes))),
                "candidate_universe_csv": str(candidate_universe),
                "candidate_universe_sha256": _sha256_file(
                    Path(str(candidate_universe))
                ),
                "fetch": bool(not getattr(args, "skip_evidence_fetch", False)),
                "max_workers": int(getattr(args, "evidence_workers", 6) or 6),
                "timeout": int(getattr(args, "evidence_timeout", 90) or 90),
                "hub_enabled": bool(
                    getattr(args, "evidence_hub_enabled", False)
                ),
                "hub_config": str(hub_config),
                "hub_config_sha256": _sha256_file(hub_config),
                "disease": getattr(args, "evidence_disease", ""),
                "max_targets": int(
                    getattr(args, "evidence_max_targets", 300) or 0
                ),
                "max_records": int(
                    getattr(args, "evidence_max_records", 1000) or 1000
                ),
                "hub_timeout": int(
                    getattr(args, "evidence_hub_timeout", 120) or 120
                ),
                "allow_network": bool(
                    getattr(args, "evidence_hub_allow_network", True)
                ),
                "strict": bool(
                    getattr(args, "evidence_hub_strict", False)
                ),
                "legacy_pool_size": int(
                    getattr(args, "evidence_legacy_pool_size", 50) or 50
                ),
                "candidate_expansion_max_targets": int(
                    getattr(args, "candidate_expansion_max_targets", 1000)
                    or 0
                ),
                "benchmark_positive_targets": str(
                    getattr(args, "benchmark_positive_targets", "")
                ),
                "benchmark_negative_targets": str(
                    getattr(args, "benchmark_negative_targets", "")
                ),
                "benchmark_top_n": int(
                    getattr(args, "benchmark_top_n", 20) or 20
                ),
                "target_priority_weights": _json_sorted(
                    getattr(args, "target_priority_weights", {})
                ),
            }
        )
    elif code == "04":
        payload.update(
            {
                "single_cell_root": str(ctx.get("single_cell_root") or ""),
                "skip_pseudobulk": bool(getattr(args, "skip_pseudobulk", False)),
            }
        )
    elif code == "05":
        evidence_priority = (
            integration / "evidence_hub" / "target_priority.csv"
        )
        payload.update(
            {
                "docking_config": str(docking_config),
                "docking_config_sha256": _sha256_file(docking_config),
                "case_label": getattr(args, "case_label", None),
                "normal_label": getattr(args, "normal_label", None),
                "ko_top_n": getattr(args, "ko_top_n", None),
                "depmap_csv": getattr(args, "depmap_csv", None),
                "skip_knockout": bool(getattr(args, "skip_knockout", False)),
                "evidence_priority_csv": str(evidence_priority),
                "evidence_priority_sha256": _sha256_file(
                    evidence_priority
                ),
                "target_priority_weights": _json_sorted(
                    getattr(args, "target_priority_weights", {})
                ),
            }
        )
    elif code == "06":
        key_genes = ctx.get("key_genes_path") or integration / "key_genes.csv"
        evidence = ctx.get("evidence_path") or integration / "gene_evidence.csv"
        target_priority = (
            integration / "integrated_target_priority.csv"
        )
        if not target_priority.exists():
            target_priority = integration / "target_priority.csv"
        # The ChEMBL ligand library is fetched at docking time and persisted to
        # outputs/integration/ligands/ (with a per-gene sha256 recorded in
        # docking_summary.json). It is an output, so it cannot be part of the
        # input signature without permanently invalidating the stage; the
        # persisted library is reused on rerun/resume instead.
        payload.update(
            {
                "docking_config": str(docking_config),
                "docking_config_sha256": _sha256_file(docking_config),
                "key_genes_csv": str(key_genes),
                "key_genes_sha256": _sha256_file(Path(str(key_genes))),
                "evidence_csv": str(evidence),
                "evidence_sha256": _sha256_file(Path(str(evidence))),
                "target_priority_csv": str(target_priority),
                "target_priority_sha256": _sha256_file(target_priority),
                "max_targets": int(getattr(args, "docking_targets", 3) or 3),
                "ligand_library": getattr(args, "ligand_library", None),
                "skip_docking": bool(getattr(args, "skip_docking", False)),
                "allow_review": bool(
                    getattr(args, "allow_review_docking", False)
                ),
            }
        )
    elif code == "07":
        payload.update(
            {
                "docking_config": str(docking_config),
                "docking_config_sha256": _sha256_file(docking_config),
                "md": _json_sorted(
                    getattr(args, "md_simulation", DEFAULT_MD_SIMULATION) or {}
                ),
                "handoff": _json_sorted(
                    getattr(args, "handoff", DEFAULT_HANDOFF) or {}
                ),
                "docking_ml": _json_sorted(
                    getattr(args, "docking_ml", DEFAULT_DOCKING_ML) or {}
                ),
                "training_csv": getattr(args, "docking_ml_training_csv", None),
                "training_csv_sha256": _sha256_file(
                    Path(
                        str(getattr(args, "docking_ml_training_csv", "") or "")
                    )
                ),
            }
        )
    elif code == "08":
        payload.update(
            {
                "docking_config": str(docking_config),
                "docking_config_sha256": _sha256_file(docking_config),
                "network_toxicology": _json_sorted(
                    getattr(
                        args,
                        "network_toxicology",
                        DEFAULT_NETWORK_TOXICOLOGY,
                    )
                    or {}
                ),
            }
        )
    elif code == "09":
        payload.update(
            {
                "docking_config": str(docking_config),
                "docking_config_sha256": _sha256_file(docking_config),
                "faers": _json_sorted(
                    getattr(args, "faers", DEFAULT_FAERS) or {}
                ),
            }
        )
    elif code == "10":
        payload.update(
            {
                "workdir": str(workdir),
                "single_cell_root": str(ctx.get("single_cell_root") or ""),
                "top_n": int(getattr(args, "feedback_top_n", 12) or 12),
                "max_features": int(getattr(args, "feedback_max_features", 8) or 8),
                "timeout": int(getattr(args, "feedback_timeout", 3600) or 3600),
                "skip_cell_feedback": bool(
                    getattr(args, "skip_cell_feedback", False)
                ),
            }
        )
    elif code == "11":
        payload.update(
            {
                "workdir": str(workdir),
                "single_cell_root": str(ctx.get("single_cell_root") or ""),
                "docking_config": str(docking_config),
                "docking_config_sha256": _sha256_file(docking_config),
            }
        )

    payload["config_sha256"] = _sha256_file(config_path)
    return hashlib.sha256(
        json.dumps(
            payload,
            sort_keys=True,
            ensure_ascii=False,
            default=str,
        ).encode("utf-8")
    ).hexdigest()




def _run_context_path(workdir: Path) -> Path:
    return _stage_dir(workdir) / "run_context.json"


def _read_run_context(workdir: Path) -> dict:
    return _read_json(_run_context_path(workdir))


def _write_run_context(workdir: Path, single_cell_root: Path) -> None:
    stage_dir = _stage_dir(workdir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    write_json(
        _run_context_path(workdir),
        {"single_cell_root": str(single_cell_root)},
    )


def _recorded_single_cell_root(workdir: Path) -> Path | None:
    root = _read_run_context(workdir).get("single_cell_root")
    if root:
        return Path(str(root)).resolve()
    summary = _read_json(_integration_dir(workdir) / "key_genes_summary.json")
    deg_table = summary.get("deg_table")
    if not deg_table:
        return None
    parts = Path(str(deg_table)).resolve().parts
    for i in range(len(parts) - 2):
        if parts[i].lower() == "results" and parts[i + 1].lower() == "data":
            return Path(*parts[:i])
    return None


def _dataset_mode_from_root(root: Path) -> str:
    summary_path = root / "results" / "summary.json"
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            mode = str(summary.get("dataset_mode", "")).strip()
            if mode:
                return mode
        except (OSError, ValueError) as exc:
            log.warning("could not read dataset mode from %s: %s", summary_path, exc)
    return "single_cell"


def _invalidate_markers_for_changed_root(
    workdir: Path,
    single_cell_root: Path,
) -> bool:
    recorded = _recorded_single_cell_root(workdir)
    if recorded is None or recorded == single_cell_root.resolve():
        return False
    stage_dir = _stage_dir(workdir)
    if not stage_dir.exists():
        return False
    for marker in stage_dir.glob("*.done"):
        marker.unlink(missing_ok=True)
    log.warning(
        "single-cell output changed from %s to %s; resetting stage markers",
        recorded,
        single_cell_root,
    )
    return True


def _clear_downstream_markers(workdir: Path, from_index: int = 1) -> None:
    for code, name, _description in STAGES[from_index:]:
        marker = _marker(workdir, code, name)
        if marker.exists():
            marker.unlink()
            log.warning(
                "stage rerun invalidates stage %s %s marker",
                code,
                name,
            )


def _prune_stale_stage_markers(workdir: Path) -> int:
    """Remove legacy markers that no longer match the current stage list."""
    stage_dir = _stage_dir(workdir)
    if not stage_dir.exists():
        return 0
    current = {f"{code}_{name}" for code, name, _description in STAGES}
    removed = 0
    for marker in stage_dir.glob("*.done"):
        if marker.stem not in current:
            marker.unlink(missing_ok=True)
            removed += 1
    if removed:
        log.warning(
            "removed %s stale stage marker(s) from %s",
            removed,
            stage_dir,
        )
    return removed


def _single_cell_outputs_ready(root: Path) -> bool:
    """True when the single-cell stage produced the files later stages need."""
    return (
        (root / "results" / "pipeline_complete.json").exists()
        and (
            (
                root
                / "results"
                / "data"
                / "05_deg"
                / "fig_09_deg_significant.csv"
            ).exists()
            or (
                root
                / "results"
                / "data"
                / "05_deg"
                / "fig_08_deg_all.csv"
            ).exists()
        )
    )


def _stage_output_paths(
    code: str,
    workdir: Path,
    ctx: dict,
    args,
) -> list[Path]:
    if code == "01":
        root = ctx.get("single_cell_root")
        if root:
            return [Path(str(root)) / "results" / "pipeline_complete.json"]
        return [workdir / "results" / "pipeline_complete.json"]
    rels = list(STAGE_OUTPUTS.get(code, ()))
    if code == "05" and getattr(args, "skip_knockout", False):
        rels = ("outputs/integration/knockout_summary.json",)
    elif code == "06" and getattr(args, "skip_docking", False):
        rels = ("outputs/integration/docking_summary.json",)
    elif code == "10" and getattr(args, "skip_cell_feedback", False):
        rels = ("outputs/integration/cell_feedback/cell_feedback_summary.json",)
    return [workdir / rel for rel in rels]


def _stage_outputs_ready(code: str, workdir: Path, ctx: dict, args) -> bool:
    if code == "01":
        return _single_cell_outputs_ready(ctx.get("single_cell_root"))
    try:
        return all(
            path.exists() and path.stat().st_size > 0
            for path in _stage_output_paths(code, workdir, ctx, args)
        )
    except OSError as exc:
        log.warning(
            "could not stat stage %s outputs in %s: %s", code, workdir, exc
        )
        return False


def _verify_stage_outputs(code: str, workdir: Path, ctx: dict, args) -> None:
    missing = [
        path
        for path in _stage_output_paths(code, workdir, ctx, args)
        if not path.exists() or path.stat().st_size == 0
    ]
    if missing:
        raise IntegrationError(
            f"stage {code} completed but required outputs are missing: "
            + ", ".join(str(p) for p in missing)
        )


def _resolve_path(value, base: Path) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def find_rscript() -> str:
    return orchestrator.find_rscript()


def export_pseudobulk(single_cell_root: Path, out_dir: Path) -> dict:
    """Aggregate single-cell counts by sample with the bundled R helper."""
    try:
        rscript = find_rscript()
    except RuntimeError as exc:
        raise IntegrationError(str(exc)) from exc
    out_dir.mkdir(parents=True, exist_ok=True)
    script = APP_ROOT / "src" / "pipeline" / "export_pseudobulk.R"
    if not script.exists():
        raise IntegrationError(f"pseudobulk export script missing: {script}")
    log.info("exporting pseudobulk matrix (Rscript: %s)", rscript)
    proc = subprocess.run(
        [rscript, str(script), str(single_cell_root), str(out_dir)],
        cwd=str(APP_ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=1800,
    )
    if proc.returncode != 0:
        raise IntegrationError(
            "pseudobulk export failed:\n" + (proc.stderr or proc.stdout)[-3000:]
        )
    required = ["pseudobulk_expression.csv", "pseudobulk_metadata.csv"]
    missing = [name for name in required if not (out_dir / name).exists()]
    if missing:
        raise IntegrationError(f"pseudobulk export missing files: {missing}")
    return {
        "expression_csv": str(out_dir / "pseudobulk_expression.csv"),
        "metadata_csv": str(out_dir / "pseudobulk_metadata.csv"),
    }


def _http_json(url: str, payload: dict | None = None, timeout: int = 90) -> dict:
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    data = None
    if payload is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(payload).encode("utf-8")
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8", errors="replace"))
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"HTTP request failed: {url}: {last_error}")


def _record_fetch_failure(
    failures: list[dict] | None,
    source: str,
    exc: Exception,
) -> None:
    """Log an evidence fetch failure and collect it for the stage summary."""
    log.warning("evidence fetch failed (%s): %s", source, exc)
    if failures is not None:
        failures.append({"source": source, "error": str(exc)})


def _mygene_info(
    gene: str,
    timeout: int,
    failures: list[dict] | None = None,
) -> dict:
    payload = {
        "q": gene,
        "scopes": "symbol",
        "fields": "entrezgene,uniprot,ensembl.gene",
        "species": "human",
        "size": 5,
    }
    url = "https://mygene.info/v3/query"
    try:
        hits = _http_json(
            url,
            payload,
            timeout,
        )
    except Exception as exc:  # noqa: BLE001
        _record_fetch_failure(failures, f"mygene:{gene}:{url}", exc)
        return {}
    for hit in hits or []:
        if str(hit.get("query", "")).upper() != gene:
            continue
        uniprot = ""
        if isinstance(hit.get("uniprot"), dict):
            uniprot = (
                hit["uniprot"].get("Swiss-Prot")
                or hit["uniprot"].get("SWISSPROT")
                or ""
            )
        elif isinstance(hit.get("uniprot"), str):
            uniprot = hit["uniprot"]
        ensembl = ""
        if isinstance(hit.get("ensembl"), dict):
            ensembl = hit["ensembl"].get("gene") or ""
        elif isinstance(hit.get("ensembl"), str):
            ensembl = hit["ensembl"]
        return {
            "entrez": str(hit.get("entrezgene") or ""),
            "uniprot": uniprot,
            "ensembl": ensembl,
        }
    return {}


def _chembl_evidence(
    uniprot: str,
    timeout: int,
    gene: str = "",
    failures: list[dict] | None = None,
) -> tuple[str, int]:
    if not uniprot:
        return "", 0
    chembl_id = ""
    target_url = (
        "https://www.ebi.ac.uk/chembl/api/data/target.json"
        f"?target_components__accession={uniprot}&limit=50"
    )
    try:
        res = _http_json(
            target_url,
            timeout=timeout,
        )
        targets = res.get("targets") or []
        for target in targets:
            if str(target.get("organism", "")).lower() == "homo sapiens":
                chembl_id = target.get("target_chembl_id") or ""
                break
        if not chembl_id and targets:
            chembl_id = targets[0].get("target_chembl_id") or ""
    except Exception as exc:  # noqa: BLE001
        _record_fetch_failure(
            failures, f"chembl_target:{gene or uniprot}:{target_url}", exc
        )
        chembl_id = ""
    if not chembl_id:
        return chembl_id, 0
    activity_url = (
        "https://www.ebi.ac.uk/chembl/api/data/activity.json"
        f"?target_chembl_id={chembl_id}&limit=1"
    )
    try:
        res = _http_json(
            activity_url,
            timeout=timeout,
        )
        meta = res.get("page_meta") or {}
        return chembl_id, int(meta.get("total_count") or 0)
    except Exception as exc:  # noqa: BLE001
        _record_fetch_failure(
            failures, f"chembl_activity:{gene or chembl_id}:{activity_url}", exc
        )
        return chembl_id, 0


def _rcsb_evidence(
    uniprot: str,
    timeout: int,
    gene: str = "",
    failures: list[dict] | None = None,
) -> tuple[int, list[str]]:
    if not uniprot:
        return 0, []
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_accession"
                ),
                "operator": "exact_match",
                "value": uniprot,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 5}},
    }
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    try:
        res = _http_json(
            url,
            query,
            timeout,
        )
    except Exception as exc:  # noqa: BLE001
        _record_fetch_failure(failures, f"rcsb:{gene or uniprot}:{url}", exc)
        return 0, []
    total = int(res.get("total_count") or 0)
    ids = [
        item.get("identifier", "")
        for item in res.get("result_set") or []
        if item.get("identifier")
    ]
    return total, ids[:5]


def _empty_evidence(gene: str) -> dict:
    return {
        "gene": gene,
        "entrez": "",
        "uniprot": "",
        "ensembl": "",
        "chembl_target_id": "",
        "known_ligands": 0,
        "chembl_bioactivities": 0,
        "pdb_structures": 0,
        "pdb_ids": "",
        "off_target_paralogs": 0,
        "safety_concern": 0,
        "string_partners": 0,
        "string_partner_ids": "",
        "reactome_pathways": 0,
        "reactome_pathway_ids": "",
        "pharmgkb_annotations": 0,
        "pharmgkb_ids": "",
        "alphafold_structures": 0,
        "alphafold_ids": "",
        "opentargets_hits": 0,
        "opentargets_target_ids": "",
        "kegg_pathways": 0,
        "kegg_pathway_ids": "",
        "database_sources": "",
    }


def _evidence_for_gene(gene: str, timeout: int = 90) -> dict:
    failures: list[dict] = []
    info = _mygene_info(gene, timeout, failures)
    uniprot = info.get("uniprot") or ""
    ensembl = info.get("ensembl") or ""
    chembl_id, bioactivities = _chembl_evidence(
        uniprot, timeout, gene, failures
    )
    pdb_count, pdb_ids = _rcsb_evidence(uniprot, timeout, gene, failures)
    database = evidence_mod.collect_gene_database_evidence(
        gene,
        max_items=10,
        timeout=timeout,
        uniprot=uniprot,
        ensembl=ensembl,
    )
    row = _empty_evidence(gene)
    row.update(
        {
            "entrez": info.get("entrez") or "",
            "uniprot": uniprot,
            "ensembl": ensembl,
            "chembl_target_id": chembl_id,
            "known_ligands": bioactivities,
            "chembl_bioactivities": bioactivities,
            "pdb_structures": pdb_count,
            "pdb_ids": ",".join(pdb_ids),
        }
    )
    row.update(database)
    log.info(
        "evidence %s: ligands=%s chembl=%s pdb=%s string=%s reactome=%s "
        "opentargets=%s kegg=%s sources=%s",
        gene,
        bioactivities,
        chembl_id,
        pdb_count,
        database.get("string_partners", 0),
        database.get("reactome_pathways", 0),
        database.get("opentargets_hits", 0),
        database.get("kegg_pathways", 0),
        database.get("database_sources", ""),
    )
    row["_fetch_failures"] = failures
    return row


def ensure_gene_evidence(
    genes: list[str],
    workdir: Path,
    fetch: bool = True,
    max_workers: int = 6,
    timeout: int = 90,
) -> pd.DataFrame:
    """Return per-gene evidence, reusing the local cache when possible."""
    out_path = _integration_dir(workdir) / "gene_evidence.csv"
    cache: dict[str, dict] = {}
    if out_path.exists():
        try:
            old = pd.read_csv(out_path)
            for _, row in old.iterrows():
                cache[str(row["gene"])] = row.to_dict()
        except (OSError, ValueError) as exc:
            log.warning("could not reuse cached evidence %s: %s", out_path, exc)
            cache = {}

    failures: list[dict] = []
    missing = [gene for gene in genes if gene not in cache]
    if missing and fetch:
        log.info("fetching evidence for %s genes", len(missing))
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [
                pool.submit(_evidence_for_gene, gene, timeout)
                for gene in missing
            ]
            for gene, future in zip(missing, futures):
                try:
                    row = future.result()
                    failures.extend(row.pop("_fetch_failures", []))
                    cache[str(row["gene"])] = row
                except Exception as exc:  # noqa: BLE001
                    log.warning("evidence fetch failed for %s: %s", gene, exc)
                    failures.append({"source": f"gene:{gene}", "error": str(exc)})

    rows = []
    for gene in genes:
        row = cache.get(gene) or _empty_evidence(gene)
        rows.append({key: row.get(key, _empty_evidence(gene)[key]) for key in EVIDENCE_COLUMNS})
    frame = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, index=False)
    write_json(
        out_path.parent / "evidence_summary.json",
        {
            "genes": len(frame),
            "fetched": len(missing) if fetch else 0,
            "evidence_failures": failures,
        },
    )
    if failures:
        log.warning(
            "gene evidence recorded %s fetch failure(s); see %s",
            len(failures),
            out_path.parent / "evidence_summary.json",
        )
    log.info(
        "gene evidence ready: %s genes, %s with PDB structures",
        len(frame),
        int((pd.to_numeric(frame["pdb_structures"], errors="coerce").fillna(0) > 0).sum()),
    )
    return frame


def build_knockout_inputs(
    single_cell_root: Path,
    workdir: Path,
    skip_pseudobulk: bool = False,
) -> dict:
    """Write expression/metadata/prognosis/druggability inputs for knockout."""
    ko_dir = workdir / "data" / "knockout"
    ko_dir.mkdir(parents=True, exist_ok=True)
    expression_dst = ko_dir / "expression.csv"
    metadata_dst = ko_dir / "metadata.csv"

    if not _knockout_inputs_ready(expression_dst, metadata_dst):
        pseudo_dir = ko_dir / "_pseudobulk"
        if (
            not (pseudo_dir / "pseudobulk_expression.csv").exists()
            or not (pseudo_dir / "pseudobulk_metadata.csv").exists()
            or not _knockout_inputs_ready(
                pseudo_dir / "pseudobulk_expression.csv",
                pseudo_dir / "pseudobulk_metadata.csv",
            )
        ):
            if skip_pseudobulk:
                raise IntegrationError(
                    "pseudobulk files are missing and --skip-pseudobulk was set; "
                    "run the single-cell pipeline first"
                )
            if pseudo_dir.exists():
                shutil.rmtree(pseudo_dir)
            export_pseudobulk(single_cell_root, pseudo_dir)
        if expression_dst.exists():
            expression_dst.unlink()
        if metadata_dst.exists():
            metadata_dst.unlink()
        shutil.copyfile(
            pseudo_dir / "pseudobulk_expression.csv",
            expression_dst,
        )
        shutil.copyfile(
            pseudo_dir / "pseudobulk_metadata.csv",
            metadata_dst,
        )

    expression = pd.read_csv(expression_dst)
    genes = expression.iloc[:, 0].astype(str).tolist()
    metadata = pd.read_csv(metadata_dst)
    if "sample" not in metadata.columns or "condition" not in metadata.columns:
        raise IntegrationError(
            "pseudobulk metadata must contain 'sample' and 'condition' columns"
        )

    evidence_path = _integration_dir(workdir) / "gene_evidence.csv"
    evidence = (
        pd.read_csv(evidence_path)
        if evidence_path.exists()
        else pd.DataFrame([_empty_evidence(gene) for gene in genes])
    )
    evidence = evidence.drop_duplicates("gene", keep="first")

    druggability = evidence[
        [
            "gene",
            "known_ligands",
            "chembl_bioactivities",
            "pdb_structures",
            "off_target_paralogs",
            "safety_concern",
        ]
    ].copy()
    druggability.to_csv(ko_dir / "druggability.csv", index=False)
    off_target = evidence[
        ["gene", "off_target_paralogs", "safety_concern"]
    ].copy()
    off_target.to_csv(ko_dir / "off_target.csv", index=False)

    # No clinical prognosis cohort is available at this stage, so we must not
    # fabricate hazard ratios. Write an empty (header-only) prognosis frame:
    # docking.knockout._prognosis_scores returns (None, None) for an empty CSV
    # and drops the prognosis dimension from target scoring instead of treating
    # every gene as HR=1.0.
    prognosis = pd.DataFrame(columns=["gene", "hr", "p"])
    prognosis.to_csv(ko_dir / "prognosis.csv", index=False)

    summary = {
        "expression_csv": str(expression_dst),
        "metadata_csv": str(metadata_dst),
        "prognosis_csv": str(ko_dir / "prognosis.csv"),
        # Placeholder only: the prognosis dimension is excluded from scoring
        # until a real clinical cohort is supplied.
        "prognosis_source": "not_available",
        "druggability_csv": str(ko_dir / "druggability.csv"),
        "off_target_csv": str(ko_dir / "off_target.csv"),
        "genes": len(genes),
        "samples": int(expression.shape[1] - 1),
        "groups": metadata["condition"].drop_duplicates().tolist(),
    }
    write_json(ko_dir / "inputs_summary.json", summary)
    log.info(
        "knockout inputs ready: %s genes x %s samples -> %s",
        summary["genes"],
        summary["samples"],
        ko_dir,
    )
    return summary


def _knockout_inputs_ready(expression_path: Path, metadata_path: Path) -> bool:
    """Return False when cached knockout inputs are malformed or stale."""
    try:
        if not expression_path.exists() or not metadata_path.exists():
            return False
        expr = pd.read_csv(expression_path)
        meta = pd.read_csv(metadata_path)
        if expr.empty or meta.empty or len(expr.columns) < 2:
            return False
        numeric = expr.drop(columns=[expr.columns[0]]).apply(
            pd.to_numeric, errors="coerce"
        )
        if numeric.dropna(how="all").empty:
            return False
        if not {"sample", "condition"}.issubset(meta.columns):
            return False
        if meta["sample"].isna().any() or meta["sample"].duplicated().any():
            return False
        if meta["condition"].nunique() < 2:
            return False
        return True
    except Exception as exc:
        log.warning(
            "could not validate pseudobulk inputs %s / %s: %s",
            expression_path,
            metadata_path,
            exc,
        )
        return False


def run_knockout_stage(
    workdir: Path,
    docking_config: Path,
    inputs: dict,
    case_label: str | None = None,
    normal_label: str | None = None,
    ko_top_n: int | None = None,
    depmap_csv: str | None = None,
    ppi_network_csv: str | None = None,
) -> dict:
    overrides = {
        "workdir": str(workdir),
        "expression_csv": inputs["expression_csv"],
        "metadata_csv": inputs["metadata_csv"],
        "prognosis_csv": inputs["prognosis_csv"],
        "druggability_csv": inputs["druggability_csv"],
        "off_target_csv": inputs["off_target_csv"],
    }
    if case_label:
        overrides["case_label"] = case_label
    if normal_label:
        overrides["normal_label"] = normal_label
    if ko_top_n:
        overrides["ko_top_n"] = int(ko_top_n)
    if depmap_csv:
        overrides["depmap_csv"] = depmap_csv
    if ppi_network_csv:
        overrides["ppi_network_csv"] = ppi_network_csv

    metadata = pd.read_csv(inputs["metadata_csv"])
    if "cell_type" in metadata.columns:
        overrides["cell_type_column"] = "cell_type"

    cfg = load_config(docking_config, overrides)
    ko_summary = run_knockout(cfg, log)
    val_summary = export_validation(cfg, log)
    result = {"knockout": ko_summary, "validation": val_summary}
    write_json(
        _integration_dir(workdir) / "knockout_summary.json",
        result,
    )
    log.info(
        "knockout + validation complete: %s genes scored, %s candidates",
        ko_summary.get("genes_scored", 0),
        val_summary.get("candidates", 0),
    )
    return result


def _split_pdb_ids(value) -> list[str]:
    text = str(value or "")
    return [
        part.strip().upper()
        for part in re.split(r"[,;|\s]+", text)
        if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", part.strip())
    ]


def _download_pdb(pdb_id: str, target_dir: Path, timeout: int = 90) -> Path | None:
    out = target_dir / "data" / "receptors" / f"{pdb_id}.pdb"
    if out.exists() and out.stat().st_size > 0:
        return out
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "Mozilla/5.0"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            if "\nATOM" not in text and not text.startswith("ATOM"):
                log.warning("PDB %s has no ATOM records", pdb_id)
                return None
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(text, encoding="utf-8")
            return out
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt < 2:
                time.sleep(1.5 * (attempt + 1))
    log.warning("PDB download failed for %s: %s", pdb_id, last_error)
    return None


def _valid_docking_box(center, size) -> bool:
    try:
        return (
            len(center) == 3
            and len(size) == 3
            and all(np.isfinite(center))
            and all(np.isfinite(size))
            and all(float(value) > 0 for value in size)
        )
    except (TypeError, ValueError) as exc:
        log.debug(
            "invalid docking box center=%r size=%r: %s", center, size, exc
        )
        return False


def _ligand_library_dir(workdir: Path) -> Path:
    return _integration_dir(workdir) / "ligands"


def _persisted_ligand_library(workdir: Path, gene: str) -> Path | None:
    """Return the ligand library persisted by an earlier docking run."""
    base = _ligand_library_dir(workdir) / safe_name(str(gene), str(gene))
    for suffix in (".csv", ".smi", ".sdf"):
        candidate = base.with_suffix(suffix)
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    return None


def _persist_ligand_library(
    workdir: Path,
    gene: str,
    library: Path,
) -> tuple[Path, str]:
    """Copy the per-target ligand library into the stage dir and hash it.

    The library may come from a live ChEMBL fetch, which is not part of the
    stage signature. Persisting it here keeps the docked ligand set auditable
    and reusable on resume/rerun.
    """
    dest_dir = _ligand_library_dir(workdir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{safe_name(str(gene), str(gene))}{library.suffix or '.csv'}"
    try:
        if library.resolve() != dest.resolve():
            shutil.copyfile(library, dest)
    except OSError as exc:
        log.warning("could not persist ligand library for %s: %s", gene, exc)
        return library, _sha256_file(library)
    return dest, _sha256_file(dest)


def _prepare_ligand_library(
    target_dir: Path,
    known_ligands_csv: Path,
    fallback_library: str | None,
    pdb_path: Path | None = None,
) -> Path | None:
    if known_ligands_csv.exists():
        try:
            df = pd.read_csv(known_ligands_csv)
            smi_col = next(
                (c for c in df.columns if c.lower() in {"smiles", "canonical_smiles"}),
                None,
            )
            if smi_col is not None:
                df = df[df[smi_col].notna()]
                df["smiles"] = df[smi_col].astype(str).str.strip()
                df = df[df["smiles"] != ""].drop_duplicates("smiles")
                df = df.head(50)
                if not df.empty:
                    id_col = "ligand_id" if "ligand_id" in df.columns else None
                    df["ID"] = [
                        safe_name(str(row.get(id_col, "")), f"ligand_{i + 1}")
                        if id_col
                        else f"ligand_{i + 1}"
                        for i, (_, row) in enumerate(df.iterrows())
                    ]
                    lib = target_dir / "data" / "ligands" / "library.csv"
                    lib.parent.mkdir(parents=True, exist_ok=True)
                    df[["ID", "smiles"]].rename(columns={"smiles": "SMILES"}).to_csv(
                        lib,
                        index=False,
                    )
                    log.info("using %s known ligands for %s", len(df), target_dir.name)
                    return lib
        except Exception as exc:  # noqa: BLE001
            log.warning("known-ligand CSV unusable: %s", exc)

    candidates = []
    if fallback_library:
        candidates.append(Path(fallback_library))
    for name in ["library.smi", "library.sdf", "library.csv"]:
        candidates.append(target_dir.parent.parent / "data" / "ligands" / name)
    for cand in candidates:
        if cand.exists():
            lib = target_dir / "data" / "ligands" / cand.name
            lib.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(cand, lib)
            log.info("using fallback ligand library: %s", cand)
            return lib
    if pdb_path is not None:
        ligands = _extract_cocrystal_ligands(pdb_path)
        if ligands:
            lib = target_dir / "data" / "ligands" / "cocrystal_library.csv"
            lib.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(ligands).rename(columns={"smiles": "SMILES"}).to_csv(
                lib,
                index=False,
            )
            log.info("using %s cocrystal ligands from %s", len(ligands), pdb_path.name)
            return lib
    return None


def _extract_cocrystal_ligands(pdb_path: Path, max_ligands: int = 5) -> list[dict]:
    """Extract non-water HETATM residues as SMILES when DB ligands are missing."""
    try:
        from rdkit import Chem
    except ImportError as exc:
        log.debug(
            "rdkit unavailable; skipping cocrystal ligand extraction for %s: %s",
            pdb_path,
            exc,
        )
        return []
    try:
        lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        log.warning("could not read PDB %s: %s", pdb_path, exc)
        return []
    small_ions = {
        "HOH", "WAT", "DOD", "CL", "NA", "K", "MG", "CA", "ZN",
        "SO4", "PO4", "GOL", "EDO", "ACT", "FMT", "DMS", "PEG",
        "IOD", "BR", "CO", "CU", "FE", "MN", "NI",
    }
    het_atoms: list[str] = []
    conect: list[str] = []
    for line in lines:
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname in small_ions:
                continue
            het_atoms.append(line)
        elif line.startswith("CONECT"):
            conect.append(line)
    if len(het_atoms) < 3:
        return []

    groups: dict[tuple[str, str, str], list[str]] = {}
    for line in het_atoms:
        key = (line[21], line[22:26].strip(), line[17:20].strip())
        groups.setdefault(key, []).append(line)

    ligands: list[dict] = []
    for (chain, resseq, resname), atom_lines in groups.items():
        atom_ids = set()
        for line in atom_lines:
            try:
                atom_ids.add(int(line[6:11]))
            except ValueError as exc:
                log.debug("skipping PDB atom with bad serial %r: %s", line[6:11], exc)
                continue
        block_lines = list(atom_lines)
        for line in conect:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                ids = [int(part) for part in parts[1:] if part.isdigit()]
            except ValueError as exc:
                log.debug("skipping unparseable CONECT record %r: %s", line, exc)
                continue
            if any(atom_id in atom_ids for atom_id in ids):
                block_lines.append(line)
        block = (
            "REMARK generated cocrystal ligand\n"
            + "\n".join(block_lines)
            + "\nEND\n"
        )
        try:
            mol = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=True)
            if mol is None:
                continue
            if mol.GetNumHeavyAtoms() < 6:
                continue
            smiles = Chem.MolToSmiles(mol)
            if not smiles:
                continue
            ligands.append(
                {
                    "id": f"{resname}_{chain}{resseq}",
                    "smiles": smiles,
                }
            )
            if len(ligands) >= max_ligands:
                break
        except Exception as exc:
            log.debug(
                "skipping cocrystal ligand %s_%s%s from %s: %s",
                resname,
                chain,
                resseq,
                pdb_path,
                exc,
            )
            continue
    return ligands


def run_target_docking(
    gene: str,
    workdir: Path,
    docking_config: Path,
    evidence: pd.DataFrame,
    ligand_library: str | None,
    force: bool = False,
) -> dict:
    row = evidence[evidence["gene"].astype(str) == gene]
    base = {
        "gene": gene,
        "status": "skipped",
        "pdb_id": "",
        "uniprot": "",
        "box_mode": "",
        "ligand_count": 0,
        "hits": 0,
        "best_affinity": "",
        "output_dir": "",
        "error": "",
        "ligand_library": "",
        "ligand_library_sha256": "",
    }
    if row.empty:
        base["error"] = "no evidence row"
        return base
    info = row.iloc[0].to_dict()
    uniprot_value = info.get("uniprot")
    uniprot = "" if uniprot_value is None or pd.isna(uniprot_value) else str(uniprot_value)
    pdb_ids = _split_pdb_ids(
        "" if info.get("pdb_ids") is None or pd.isna(info.get("pdb_ids"))
        else info.get("pdb_ids")
    )
    if not pdb_ids:
        base["error"] = "no PDB structure"
        return base

    target_dir = workdir / "work" / safe_name(gene, gene)
    target_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = None
    for pdb_id in pdb_ids:
        pdb_path = _download_pdb(pdb_id, target_dir)
        if pdb_path is not None:
            break
    if pdb_path is None:
        base["error"] = "PDB download failed for all candidates"
        return base

    cfg = load_config(
        docking_config,
        {
            "workdir": str(target_dir),
            "target_name": gene,
            "uniprot": uniprot,
            "pdb": pdb_path.stem,
        },
    )
    chembl_id = info.get("chembl_target_id")
    persisted = _persisted_ligand_library(workdir, gene)
    if persisted is not None:
        # Reuse the library captured by the first run so a rerun/resume cannot
        # silently dock a different ligand set from a live database fetch.
        log.info("reusing persisted ligand library for %s: %s", gene, persisted)
    else:
        if chembl_id and not pd.isna(chembl_id):
            cfg.data.setdefault("evidence", {})["chembl_target_id"] = str(chembl_id)
        try:
            evidence_mod.gather_evidence(cfg, log)
        except Exception as exc:  # noqa: BLE001
            log.warning("evidence collection failed for %s: %s", gene, exc)

    known = target_dir / "evidence" / "known_ligands.csv"
    library = _prepare_ligand_library(
        target_dir,
        persisted if persisted is not None else known,
        ligand_library,
        pdb_path=pdb_path,
    )
    if library is None:
        base["error"] = "no ligand library available"
        return base
    library_path, library_digest = _persist_ligand_library(workdir, gene, library)

    try:
        center, size, mode = box.detect_box_data(pdb_path)
    except Exception as exc:  # noqa: BLE001
        base["status"] = "failed"
        base["error"] = f"docking box detection failed: {exc}"
        base["pdb_id"] = pdb_path.stem
        return base
    base["box_mode"] = mode
    if not _valid_docking_box(center, size):
        base["status"] = "skipped"
        base["error"] = (
            f"invalid docking box (mode={mode}, center={center}, size={size})"
        )
        base["pdb_id"] = pdb_path.stem
        base["uniprot"] = uniprot
        return base
    cfg.data["receptor"]["input"] = str(pdb_path)
    cfg.data["receptor"]["output"] = str(
        target_dir / "data" / "receptors" / f"{pdb_path.stem}.pdbqt"
    )
    cfg.data["receptor"]["center"] = center
    cfg.data["receptor"]["size"] = size
    cfg.data["receptor"]["detect_input"] = None
    cfg.data["ligand"]["input"] = str(library)
    save_config(cfg, target_dir / "config" / "docking_config.json")
    log.info(
        "docking target %s: PDB %s, box mode %s, center %s size %s",
        gene,
        pdb_path.stem,
        mode,
        center,
        size,
    )

    try:
        docking_pipeline.run_pipeline(cfg, force=force)
    except (DockingError, ToolNotFoundError) as exc:
        base["status"] = "failed"
        base["error"] = str(exc)
        base["pdb_id"] = pdb_path.stem
        base["uniprot"] = uniprot
        return base

    report_dir = cfg.analysis_dir()
    summary = _read_json(report_dir / "summary.json")
    ranked = report_dir / "data" / "fig_46_47_ranked_results.csv"
    hits = 0
    best = ""
    if ranked.exists():
        try:
            ranked_df = pd.read_csv(ranked)
            hits = int((ranked_df["affinity"] <= float(cfg.get("analysis", "cutoff", -7.0))).sum())
            if "affinity" in ranked_df.columns:
                best = _min_affinity_text(ranked_df["affinity"])
        except Exception as exc:
            log.warning(
                "could not read ranked docking results %s; using summary: %s",
                ranked,
                exc,
            )
            hits = int(summary.get("hits", 0))
            best = str(summary.get("best_affinity", ""))
    else:
        hits = int(summary.get("hits", 0))
        best = str(summary.get("best_affinity", ""))

    result = {
        "gene": gene,
        "status": "ok",
        "pdb_id": pdb_path.stem,
        "uniprot": uniprot,
        "box_mode": mode,
        "ligand_count": int(summary.get("total_docked", 0)),
        "hits": hits,
        "best_affinity": best,
        "output_dir": str(cfg.output_dir),
        "error": "",
        "ligand_library": str(library_path),
        "ligand_library_sha256": library_digest,
    }
    write_json(
        target_dir / "outputs" / "integration" / "target_summary.json",
        result,
    )
    log.info(
        "docking target %s complete: %s ligands, %s hits, best %s",
        gene,
        result["ligand_count"],
        hits,
        best,
    )
    return result


def run_docking_stage(
    workdir: Path,
    docking_config: Path,
    key_genes_csv: Path,
    evidence_csv: Path,
    max_targets: int,
    ligand_library: str | None,
    force: bool = False,
    priority_csv: Path | None = None,
    allow_review: bool = False,
) -> dict:
    out_dir = _integration_dir(workdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    selection_source = "key_genes"
    selection = pd.DataFrame()
    if priority_csv is not None and Path(priority_csv).exists():
        try:
            priority = pd.read_csv(priority_csv)
            if "gene" in priority.columns and not priority.empty:
                selection = priority
                selection_source = "integrated_target_priority"
        except (OSError, ValueError):
            selection = pd.DataFrame()
    if selection.empty:
        selection = pd.read_csv(key_genes_csv)
    if "decision" in selection.columns and not allow_review:
        eligible = selection[
            selection["decision"].astype(str).str.upper().isin(
                {"GO", "CONDITIONAL_GO"}
            )
        ]
        selection = eligible
    if selection.empty:
        empty = pd.DataFrame(
            columns=[
                "gene",
                "status",
                "pdb_id",
                "uniprot",
                "box_mode",
                "ligand_count",
                "hits",
                "best_affinity",
                "output_dir",
                "error",
                "ligand_library",
                "ligand_library_sha256",
            ]
        )
        empty.to_csv(out_dir / "docking_targets.csv", index=False)
        summary = {
            "status": "skipped",
            "reason": (
                "no GO/CONDITIONAL_GO targets selected; "
                "enable allow_review to dock REVIEW targets"
            ),
            "ok": 0,
            "failed": 0,
            "skipped": 0,
            "selection_source": selection_source,
            "selection_csv": str(priority_csv or key_genes_csv),
        }
        write_json(out_dir / "docking_summary.json", summary)
        return summary
    genes = (
        selection["gene"]
        .astype(str)
        .str.strip()
        .str.upper()
        .head(max_targets)
        .tolist()
    )
    if not genes:
        summary = {"status": "skipped", "reason": "no key genes"}
        write_json(out_dir / "docking_summary.json", summary)
        return summary
    evidence = pd.read_csv(evidence_csv)
    rows = []
    for gene in genes:
        try:
            rows.append(
                run_target_docking(
                    gene,
                    workdir,
                    docking_config,
                    evidence,
                    ligand_library,
                    force=force,
                )
            )
        except Exception as exc:  # noqa: BLE001
            log.error("docking target %s crashed: %s", gene, exc)
            rows.append(
                {
                    "gene": gene,
                    "status": "failed",
                    "pdb_id": "",
                    "uniprot": "",
                    "box_mode": "",
                    "ligand_count": 0,
                    "hits": 0,
                    "best_affinity": "",
                    "output_dir": "",
                    "error": str(exc),
                    "ligand_library": "",
                    "ligand_library_sha256": "",
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "docking_targets.csv", index=False)
    ligand_digests = {
        str(row.get("gene", "")): str(row.get("ligand_library_sha256", ""))
        for _, row in frame.iterrows()
    }
    summary = {
        "status": "completed",
        "targets_requested": len(genes),
        "ok": int((frame["status"] == "ok").sum()),
        "failed": int((frame["status"] == "failed").sum()),
        "skipped": int((frame["status"] == "skipped").sum()),
        "total_hits": int(pd.to_numeric(frame["hits"], errors="coerce").fillna(0).sum()),
        "best_affinity": (
            _min_affinity_text(
                frame.loc[
                    pd.to_numeric(frame["hits"], errors="coerce").fillna(0) > 0,
                    "best_affinity",
                ]
            )
            if len(frame)
            else ""
        ),
        "output_csv": str(out_dir / "docking_targets.csv"),
        "ligand_library_digest": hashlib.sha256(
            json.dumps(ligand_digests, sort_keys=True).encode("utf-8")
        ).hexdigest(),
        "ligand_libraries": ligand_digests,
        "selection_source": selection_source,
        "selection_csv": str(priority_csv or key_genes_csv),
    }
    write_json(out_dir / "docking_summary.json", summary)
    log.info(
        "docking stage complete: %s ok / %s failed / %s skipped",
        summary["ok"],
        summary["failed"],
        summary["skipped"],
    )
    return summary






def _stage_single_cell(args, workdir: Path, ctx: dict) -> None:
    root = ctx["single_cell_root"]
    accession = str(getattr(args, "accession", "") or "").strip().upper()
    dataset_mode = "single_cell"
    if args.skip_scrna:
        if not (root / "results" / "pipeline_complete.json").exists():
            raise IntegrationError(
                f"single-cell outputs not found under {root}; remove --skip-scrna"
            )
        log.info("using existing single-cell outputs: %s", root)
        dataset_mode = _dataset_mode_from_root(root)
    else:
        if accession:
            manifest_path = (
                root / "data" / f"{accession}_manifest.json"
            )
            try:
                manifest = json.loads(
                    manifest_path.read_text(encoding="utf-8")
                )
            except (OSError, ValueError) as exc:
                log.warning("could not read manifest %s: %s", manifest_path, exc)
                manifest = {}
            dataset_mode = (
                "single_cell"
                if manifest.get("mode") == "single_cell"
                else "sample_level"
            )
        code = orchestrator.run_pipeline(
            args.force,
            args.skip_download,
            args.skip_deps,
            args.accession,
            str(root),
            args.species,
            getattr(args, "ml_model", "xgb"),
        )
        if code == 98:
            raise PauseRequested("single-cell pipeline paused; run again to resume")
        if code != 0:
            raise IntegrationError(f"single-cell pipeline exited with code {code}")
        dataset_mode = _dataset_mode_from_root(root)
    ctx["dataset_mode"] = dataset_mode

    metrics = write_qc_metrics(workdir, root, getattr(args, "qc_gate", {}))
    gate = metrics["qc_gate"]
    if gate["status"] == "fail":
        raise IntegrationError("QC gate failed: " + str(gate.get("summary", "")))
    if gate["status"] == "warn":
        log.warning("QC gate warning: %s", gate.get("summary", ""))
    ctx["qc_metrics"] = metrics

    da_config = getattr(args, "differential_abundance", {}) or {}
    if da_config.get("enabled", True) and dataset_mode == "single_cell":
        ctx["differential_abundance"] = run_differential_abundance(
            root,
            _integration_dir(workdir),
            da_config,
        )
    else:
        reason = (
            "differential abundance is not applicable to sample-level datasets"
            if dataset_mode != "single_cell"
            else "differential abundance disabled by configuration"
        )
        summary = {
            "status": "skipped",
            "reason": reason,
            "output_csv": "",
        }
        write_json(
            _integration_dir(workdir) / "differential_abundance_summary.json",
            summary,
        )
        ctx["differential_abundance"] = summary


def _stage_key_targets(args, workdir: Path, ctx: dict) -> None:
    frame = extract_key_genes(
        ctx["single_cell_root"],
        _integration_dir(workdir),
        top_n=args.top_genes,
        keep_all=args.keep_all_genes,
        advanced_priority_csv=getattr(args, "advanced_priority_csv", None),
        universe_size=int(getattr(args, "candidate_universe_size", 0) or 0),
    )
    ctx["key_genes_path"] = _integration_dir(workdir) / "key_genes.csv"
    ctx["key_genes"] = frame
    ctx["candidate_universe_path"] = (
        _integration_dir(workdir) / "candidate_universe.csv"
    )


def _candidate_universe_frame(workdir: Path, ctx: dict) -> pd.DataFrame:
    candidate_path = (
        ctx.get("candidate_universe_expanded_path")
        or _integration_dir(workdir)
        / "candidate_universe_evidence_expanded.csv"
    )
    if not Path(candidate_path).exists():
        candidate_path = (
            ctx.get("candidate_universe_path")
            or _integration_dir(workdir) / "candidate_universe.csv"
        )
    if not Path(candidate_path).exists():
        candidate_path = (
            ctx.get("key_genes_path")
            or _integration_dir(workdir) / "key_genes.csv"
        )
    key_path = (
        ctx.get("key_genes_path")
        or _integration_dir(workdir) / "key_genes.csv"
    )
    source = candidate_path if Path(candidate_path).exists() else key_path
    if not Path(source).exists():
        raise IntegrationError(
            "candidate_universe.csv/key_genes.csv missing; run stage 02 first"
        )
    frame = pd.read_csv(source)
    if "gene" not in frame.columns:
        raise IntegrationError(f"candidate gene table has no gene column: {source}")
    frame["gene"] = frame["gene"].astype(str).str.strip().str.upper()
    return frame[frame["gene"] != ""].drop_duplicates("gene", keep="first")


def _legacy_evidence_records(
    frame: pd.DataFrame,
    source_version: str = "full-pipeline legacy evidence",
) -> list[EvidenceRecord]:
    """Convert legacy per-gene evidence counts into canonical evidence records."""
    records: list[EvidenceRecord] = []
    if frame is None or frame.empty or "gene" not in frame.columns:
        return records

    def number(row: pd.Series, column: str) -> float:
        value = pd.to_numeric(row.get(column), errors="coerce")
        return 0.0 if pd.isna(value) else float(value)

    for _, row in frame.iterrows():
        gene = str(row.get("gene") or "").strip().upper()
        if not gene:
            continue
        chembl = max(
            number(row, "chembl_bioactivities"),
            number(row, "known_ligands"),
        )
        if chembl > 0:
            records.append(
                EvidenceRecord(
                    source="LegacyTargetEvidence",
                    source_record_id=f"{gene}:chembl_bioactivity",
                    evidence_type="direct_bioactivity_aggregate",
                    subject_type="compound",
                    subject_id="ChEMBL aggregate",
                    relation="targets",
                    object_type="target",
                    object_id=gene,
                    target_symbol=gene,
                    tier=EvidenceTier.EXPERIMENTAL,
                    source_version=source_version,
                    source_group="chembl",
                    score=float(
                        np.clip(np.log1p(chembl) / np.log1p(50.0), 0.0, 1.0)
                    ),
                    payload={"bioactivity_count": chembl},
                )
            )
        pdb_count = number(row, "pdb_structures")
        if pdb_count > 0:
            records.append(
                EvidenceRecord(
                    source="LegacyTargetEvidence",
                    source_record_id=f"{gene}:pdb_structure",
                    evidence_type="experimental_structure",
                    subject_type="target",
                    subject_id=gene,
                    relation="has_structure",
                    object_type="structure",
                    object_id=f"PDB:{pdb_count}",
                    target_symbol=gene,
                    tier=EvidenceTier.EXPERIMENTAL,
                    source_version=source_version,
                    source_group="pdb",
                    score=float(min(1.0, 0.65 + 0.1 * pdb_count)),
                    payload={"pdb_structures": pdb_count},
                )
            )
        alphafold_count = number(row, "alphafold_structures")
        if alphafold_count > 0:
            records.append(
                EvidenceRecord(
                    source="LegacyTargetEvidence",
                    source_record_id=f"{gene}:alphafold_structure",
                    evidence_type="predicted_structure",
                    subject_type="target",
                    subject_id=gene,
                    relation="has_structure",
                    object_type="structure",
                    object_id="AlphaFold",
                    target_symbol=gene,
                    tier=EvidenceTier.PREDICTED,
                    source_version=source_version,
                    source_group="alphafold",
                    score=0.5,
                    payload={"alphafold_structures": alphafold_count},
                )
            )
        for source, column, source_group in (
            ("Reactome", "reactome_pathways", "reactome"),
            ("KEGG", "kegg_pathways", "kegg"),
        ):
            count = number(row, column)
            if count <= 0:
                continue
            records.append(
                EvidenceRecord(
                    source=source,
                    source_record_id=f"{gene}:{column}",
                    evidence_type="pathway_membership",
                    subject_type="target",
                    subject_id=gene,
                    relation="participates_in",
                    object_type="pathway",
                    object_id=source,
                    target_symbol=gene,
                    tier=EvidenceTier.CURATED,
                    source_version=source_version,
                    source_group=source_group,
                    score=float(min(1.0, 0.4 + 0.1 * count)),
                    payload={column: count},
                )
            )
    return records


def _evidence_target_origin(evidence: pd.DataFrame) -> dict[str, str]:
    """Assign a conservative origin label to targets found by evidence only."""
    origins: dict[str, str] = {}
    if evidence is None or evidence.empty or "target_symbol" not in evidence.columns:
        return origins
    priority = {
        "GENETIC": 4,
        "DISEASE": 3,
        "DEPENDENCY": 2,
        "EVIDENCE": 1,
    }
    for _, row in evidence.iterrows():
        gene = str(row.get("target_symbol") or "").strip().upper()
        if not gene:
            continue
        evidence_type = str(row.get("evidence_type") or "").lower()
        relation = str(row.get("relation") or "").lower()
        if "genetic" in evidence_type:
            origin = "GENETIC"
        elif relation in {"associated_with", "implicated_in"}:
            origin = "DISEASE"
        elif relation in {"dependent_on", "essential_in"}:
            origin = "DEPENDENCY"
        else:
            origin = "EVIDENCE"
        if priority[origin] > priority.get(origins.get(gene, "EVIDENCE"), 0):
            origins[gene] = origin
    return origins


def _expand_candidate_frame(
    candidate_frame: pd.DataFrame,
    evidence: pd.DataFrame,
    max_extra: int = 1000,
) -> tuple[pd.DataFrame, dict]:
    """Union DEG candidates with targets discovered by disease/genetic evidence."""
    base = candidate_frame.copy()
    base["gene"] = base["gene"].astype(str).str.strip().str.upper()
    base["candidate_origin"] = "DEG"
    if evidence is None or evidence.empty:
        return base, {"expanded_targets": 0, "origins": {"DEG": len(base)}}

    origins = _evidence_target_origin(evidence)
    known = set(base["gene"])
    score_by_gene: dict[str, float] = {}
    if "target_symbol" in evidence.columns:
        score_frame = evidence[["target_symbol", "score"]].copy()
        score_frame["target_symbol"] = (
            score_frame["target_symbol"].astype(str).str.strip().str.upper()
        )
        score_frame["score"] = pd.to_numeric(
            score_frame["score"],
            errors="coerce",
        )
        score_by_gene = (
            score_frame.groupby("target_symbol")["score"]
            .max()
            .fillna(0.0)
            .to_dict()
        )
    discovered = [
        (gene, origin, float(score_by_gene.get(gene, 0.0) or 0.0))
        for gene, origin in origins.items()
        if gene
        and gene not in known
        and re.fullmatch(r"[A-Z][A-Z0-9-]{1,14}", gene)
    ]
    discovered.sort(key=lambda item: (-item[2], item[0]))
    if max_extra > 0:
        discovered = discovered[:max_extra]
    if not discovered:
        return base, {"expanded_targets": 0, "origins": {"DEG": len(base)}}

    rows: list[dict] = []
    for gene, origin, _ in discovered:
        row = {column: np.nan for column in base.columns}
        row["gene"] = gene
        row["candidate_origin"] = origin
        rows.append(row)
    expanded = pd.concat([base, pd.DataFrame(rows)], ignore_index=True)
    expanded["gene"] = expanded["gene"].astype(str).str.strip().str.upper()
    expanded = expanded.drop_duplicates("gene", keep="first")
    origin_counts = {
        str(key): int(value)
        for key, value in expanded["candidate_origin"].value_counts().items()
    }
    return expanded, {
        "expanded_targets": len(discovered),
        "origins": origin_counts,
    }


def _local_source_fingerprints(
    config: dict,
    config_path: Path,
    out_dir: Path,
) -> dict[str, str]:
    fingerprints: dict[str, str] = {}
    for name, options in (config.get("sources") or {}).items():
        options = dict(options or {})
        if str(options.get("connector") or "") != "local":
            continue
        raw_path = str(options.get("path") or "").strip()
        if not raw_path:
            continue
        candidates = (
            Path(raw_path),
            config_path.parent / raw_path,
            out_dir / raw_path,
        )
        path = next(
            (candidate for candidate in candidates if candidate.exists()),
            candidates[0],
        )
        fingerprints[str(name)] = _sha256_file(path)
    return fingerprints


def _run_evidence_hub(
    workdir: Path,
    candidate_frame: pd.DataFrame,
    args,
    legacy_evidence: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """Collect, score and export the multi-source evidence hub."""
    out_dir = _integration_dir(workdir) / "evidence_hub"
    out_dir.mkdir(parents=True, exist_ok=True)
    strict = bool(getattr(args, "evidence_hub_strict", False))
    base_candidates, _ = _expand_candidate_frame(
        candidate_frame,
        pd.DataFrame(),
    )
    enabled = bool(getattr(args, "evidence_hub_enabled", False))
    config_value = str(
        getattr(args, "evidence_hub_config", "") or ""
    ).strip()
    if not enabled:
        summary = {
            "status": "skipped",
            "reason": "multi-source evidence hub disabled",
            "output_dir": str(out_dir),
        }
        write_json(out_dir / "evidence_hub_summary.json", summary)
        return pd.DataFrame(), summary, base_candidates
    if not config_value:
        summary = {
            "status": "skipped",
            "reason": "evidence hub config is empty",
            "output_dir": str(out_dir),
        }
        write_json(out_dir / "evidence_hub_summary.json", summary)
        if strict:
            raise IntegrationError(summary["reason"])
        return pd.DataFrame(), summary, base_candidates

    config_path = _resolve_path(config_value, APP_ROOT)
    if not config_path.exists():
        summary = {
            "status": "failed",
            "reason": f"evidence hub config not found: {config_path}",
            "output_dir": str(out_dir),
        }
        write_json(out_dir / "evidence_hub_summary.json", summary)
        log.warning(summary["reason"])
        if strict:
            raise IntegrationError(summary["reason"])
        return pd.DataFrame(), summary, base_candidates

    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
        if not isinstance(config, dict):
            raise ValueError("evidence hub config must be a JSON object")
    except Exception as exc:  # noqa: BLE001
        summary = {
            "status": "failed",
            "reason": f"invalid evidence hub config: {exc}",
            "output_dir": str(out_dir),
            "config": str(config_path),
        }
        write_json(out_dir / "evidence_hub_summary.json", summary)
        log.warning(summary["reason"])
        if strict:
            raise IntegrationError(summary["reason"]) from exc
        return pd.DataFrame(), summary, base_candidates

    max_targets = int(getattr(args, "evidence_max_targets", 300) or 0)
    genes = candidate_frame["gene"].astype(str).tolist()
    if max_targets > 0:
        genes = genes[:max_targets]
    if not genes:
        summary = {
            "status": "skipped",
            "reason": "candidate universe is empty",
            "output_dir": str(out_dir),
        }
        write_json(out_dir / "evidence_hub_summary.json", summary)
        return pd.DataFrame(), summary, base_candidates

    query_fingerprint = hashlib.sha256(
        json.dumps(
            {
                "config": str(config_path),
                "config_sha256": _sha256_file(config_path),
                "local_source_hashes": _local_source_fingerprints(
                    config,
                    config_path,
                    out_dir,
                ),
                "disease": str(
                    getattr(args, "evidence_disease", "") or "liver cancer"
                ).strip(),
                "targets": genes,
                "allow_network": bool(
                    getattr(args, "evidence_hub_allow_network", True)
                    and not getattr(args, "skip_evidence_fetch", False)
                ),
            },
            sort_keys=True,
            ensure_ascii=False,
        ).encode("utf-8")
    ).hexdigest()
    database = out_dir / "evidence.sqlite"
    previous_summary = _read_json(out_dir / "evidence_hub_summary.json")
    if (
        database.exists()
        and previous_summary.get("query_fingerprint") != query_fingerprint
    ):
        try:
            database.unlink()
        except OSError as exc:
            log.warning(
                "could not remove stale evidence database %s: %s",
                database,
                exc,
            )

    ensembl_ids: dict[str, str] = {}
    for column in ("ensembl", "ensembl_id", "gene_id"):
        if column not in candidate_frame.columns:
            continue
        for gene, value in zip(
            candidate_frame["gene"],
            candidate_frame[column],
        ):
            text = str(value or "").strip()
            if text and text.lower() != "nan" and gene not in ensembl_ids:
                ensembl_ids[str(gene).upper()] = text

    allow_network = bool(
        getattr(args, "evidence_hub_allow_network", True)
        and not getattr(args, "skip_evidence_fetch", False)
    )
    disease_name = str(
        getattr(args, "evidence_disease", "") or "liver cancer"
    ).strip()
    disease_id = str(
        getattr(args, "evidence_disease_id", "") or ""
    ).strip()
    try:
        limits = config.get("limits") or {}
        context = EvidenceContext(
            disease={"name": disease_name, "id": disease_id},
            target_symbols=genes,
            ensembl_ids=ensembl_ids,
            cache_dir=out_dir,
            max_records_per_source=int(
                getattr(
                    args,
                    "evidence_max_records",
                    limits.get("max_records_per_source", 1000),
                )
                or 1000
            ),
            timeout_seconds=int(
                getattr(
                    args,
                    "evidence_hub_timeout",
                    limits.get("timeout_seconds", 120),
                )
                or 120
            ),
            allow_network=allow_network,
            source_options={
                str(name): dict(options or {})
                for name, options in (config.get("sources") or {}).items()
            },
        )
        if strict:
            config["strict"] = True
        with SQLiteEvidenceStore(database) as store:
            hub = EvidenceHub.from_config(
                database,
                config,
                store=store,
            )
            status = hub.collect(context)
            failed_sources = []
            if {"status", "source"}.issubset(status.columns):
                failed_sources = (
                    status.loc[
                        status["status"].astype(str) == "failed",
                        "source",
                    ]
                    .astype(str)
                    .tolist()
                )
            if failed_sources:
                removed = store.delete_source_records(failed_sources)
                if removed:
                    log.warning(
                        "removed %s stale evidence records from failed sources: %s",
                        removed,
                        ", ".join(failed_sources),
                    )
            legacy_records = _legacy_evidence_records(
                legacy_evidence
                if legacy_evidence is not None
                else pd.DataFrame()
            )
            if legacy_records:
                store.upsert_evidence(legacy_records)
            all_evidence = store.records()
            expanded, expansion_summary = _expand_candidate_frame(
                candidate_frame,
                all_evidence,
                max_extra=int(
                    getattr(
                        args,
                        "candidate_expansion_max_targets",
                        1000,
                    )
                    or 0
                ),
            )
            expanded_genes = expanded["gene"].astype(str).tolist()
            context.target_symbols = expanded_genes
            benchmark_positives = [
                value.strip().upper()
                for value in str(
                    getattr(args, "benchmark_positive_targets", "") or ""
                ).replace("\n", ",").split(",")
                if value.strip()
            ]
            benchmark_negatives = [
                value.strip().upper()
                for value in str(
                    getattr(args, "benchmark_negative_targets", "") or ""
                ).replace("\n", ",").split(",")
                if value.strip()
            ]
            scored = hub.score(
                target_symbols=expanded_genes,
                benchmark_positives=benchmark_positives or None,
                benchmark_negatives=benchmark_negatives or None,
                benchmark_top_n=int(
                    getattr(args, "benchmark_top_n", 20) or 20
                ),
            )
            paths = hub.export(
                out_dir,
                context=context,
                source_status=status,
                scored=scored,
            )
            store_summary = store.summary()
        priority = scored.get("priority")
        if not isinstance(priority, pd.DataFrame):
            priority = pd.DataFrame()
        summary = {
            "status": "completed",
            "reason": "",
            "targets_queried": len(genes),
            "candidate_targets_expanded": len(expanded),
            "candidate_expansion": expansion_summary,
            "sources": {
                str(row["source"]): {
                    "status": str(row.get("status", "")),
                    "record_count": int(row.get("record_count", 0) or 0),
                    "error": str(row.get("error", "")),
                }
                for row in status.to_dict(orient="records")
            },
            "source_count": int(store_summary.get("n_sources", 0) or 0),
            "successful_source_count": int(
                sum(
                    str(row.get("status", "")) == "completed"
                    and int(row.get("record_count", 0) or 0) > 0
                    for row in status.to_dict(orient="records")
                )
            ),
            "targets_with_evidence": int(store_summary.get("n_targets", 0) or 0),
            "evidence_records": int(store_summary.get("n_records", 0) or 0),
            "outputs": paths,
            "output_dir": str(out_dir),
            "config": str(config_path),
            "disease": disease_name,
            "allow_network": allow_network,
            "query_fingerprint": query_fingerprint,
            "benchmark": scored.get("benchmark"),
        }
        write_json(out_dir / "evidence_hub_summary.json", summary)
        return priority, summary, expanded
    except Exception as exc:  # noqa: BLE001
        summary = {
            "status": "failed",
            "reason": str(exc),
            "output_dir": str(out_dir),
            "config": str(config_path),
            "disease": disease_name,
            "allow_network": allow_network,
        }
        write_json(out_dir / "evidence_hub_summary.json", summary)
        log.warning("multi-source evidence hub failed: %s", exc)
        if strict:
            raise IntegrationError(
                f"multi-source evidence hub failed: {exc}"
            ) from exc
        return pd.DataFrame(), summary, base_candidates


def _stage_evidence(args, workdir: Path, ctx: dict) -> None:
    candidate_frame = _candidate_universe_frame(workdir, ctx)
    key_path = (
        ctx.get("key_genes_path")
        or _integration_dir(workdir) / "key_genes.csv"
    )
    key_frame = pd.read_csv(key_path)
    key_genes = key_frame["gene"].astype(str).str.strip().str.upper().tolist()

    legacy_pool_size = int(
        getattr(args, "evidence_legacy_pool_size", 50) or 50
    )
    legacy_pool_size = max(
        legacy_pool_size,
        int(getattr(args, "docking_targets", 0) or 0) * 2,
    )
    discovery_genes = candidate_frame["gene"].astype(str).tolist()
    legacy_genes = list(
        dict.fromkeys(key_genes + discovery_genes[:legacy_pool_size])
    )
    legacy = ensure_gene_evidence(
        legacy_genes,
        workdir,
        fetch=not getattr(args, "skip_evidence_fetch", False),
        max_workers=int(getattr(args, "evidence_workers", 6) or 6),
        timeout=int(getattr(args, "evidence_timeout", 90) or 90),
    )

    evidence_priority, hub_summary, expanded = _run_evidence_hub(
        workdir,
        candidate_frame,
        args,
        legacy_evidence=legacy,
    )
    expanded_path = (
        _integration_dir(workdir) / "candidate_universe_evidence_expanded.csv"
    )
    expanded.to_csv(expanded_path, index=False)
    target_priority_path = _integration_dir(workdir) / "target_priority.csv"
    target_priority = build_target_priority(
        expanded_path,
        _integration_dir(workdir) / "evidence_hub" / "target_priority.csv"
        if not evidence_priority.empty
        else None,
        target_priority_path,
        weights=getattr(args, "target_priority_weights", None),
        evidence_queried=(
            str(hub_summary.get("status", "")) == "completed"
        ),
        go_min_score=float(getattr(args, "target_go_min_score", 0.75)),
        go_min_coverage=float(
            getattr(args, "target_go_min_coverage", 0.50)
        ),
        conditional_min_score=float(
            getattr(args, "target_conditional_min_score", 0.50)
        ),
    )
    ctx["candidate_universe"] = candidate_frame
    ctx["candidate_universe_path"] = (
        _integration_dir(workdir) / "candidate_universe.csv"
    )
    ctx["candidate_universe_expanded_path"] = expanded_path
    ctx["candidate_universe_expanded"] = expanded
    ctx["evidence_hub_summary"] = hub_summary
    ctx["evidence_hub_priority"] = evidence_priority
    ctx["target_priority"] = target_priority
    ctx["target_priority_path"] = target_priority_path
    ctx["evidence"] = legacy
    ctx["evidence_path"] = _integration_dir(workdir) / "gene_evidence.csv"
    write_target_priority_summary(
        target_priority,
        _integration_dir(workdir) / "target_priority_summary.json",
        evidence_hub_used=(
            str(hub_summary.get("status", "")) == "completed"
        ),
        evidence_source_count=int(hub_summary.get("source_count", 0) or 0),
        evidence_targets_queried=int(
            hub_summary.get("targets_queried", 0) or 0
        ),
        benchmark=hub_summary.get("benchmark") or {},
    )


def _stage_knockout_inputs(args, workdir: Path, ctx: dict) -> None:
    ctx["knockout_inputs"] = build_knockout_inputs(
        ctx["single_cell_root"],
        workdir,
        skip_pseudobulk=args.skip_pseudobulk,
    )
    ctx["omics_qc"] = run_omics_qc(
        workdir,
        _integration_dir(workdir),
    )


def _stage_knockout(args, workdir: Path, ctx: dict) -> None:
    if args.skip_knockout:
        summary = {"status": "skipped", "reason": "knockout disabled by arguments"}
        write_json(_integration_dir(workdir) / "knockout_summary.json", summary)
        for stale in (
            _integration_dir(workdir) / "integrated_target_priority.csv",
            _integration_dir(workdir)
            / "integrated_target_priority_summary.json",
        ):
            try:
                stale.unlink(missing_ok=True)
            except OSError as exc:
                log.warning("could not remove stale %s: %s", stale, exc)
        ctx["knockout"] = summary
        return
    ctx["knockout"] = run_knockout_stage(
        workdir,
        ctx["docking_config"],
        ctx["knockout_inputs"],
        case_label=args.case_label,
        normal_label=args.normal_label,
        ko_top_n=args.ko_top_n,
        depmap_csv=args.depmap_csv,
        ppi_network_csv=args.ppi_network_csv,
    )
    ko_candidates = [
        workdir
        / "outputs"
        / "run_001"
        / "results"
        / "04_knockout"
        / "data"
        / name
        for name in (
            "fig_52_53_ranked_knockout.csv",
            "fig_52_target_candidates.csv",
        )
    ]
    ko_path = next((path for path in ko_candidates if path.exists()), None)
    candidate_path = (
        ctx.get("candidate_universe_path")
        or _integration_dir(workdir) / "candidate_universe.csv"
    )
    evidence_priority_path = (
        _integration_dir(workdir) / "evidence_hub" / "target_priority.csv"
    )
    integrated_path = (
        _integration_dir(workdir) / "integrated_target_priority.csv"
    )
    if Path(candidate_path).exists():
        integrated = build_target_priority(
            Path(candidate_path),
            evidence_priority_path if evidence_priority_path.exists() else None,
            integrated_path,
            knockout_csv=ko_path,
            weights=getattr(args, "target_priority_weights", None),
            evidence_queried=(
                str(
                    (
                        ctx.get("evidence_hub_summary")
                        or _read_json(
                            _integration_dir(workdir)
                            / "evidence_hub"
                            / "evidence_hub_summary.json"
                        )
                    ).get("status", "")
                )
                == "completed"
            ),
            go_min_score=float(getattr(args, "target_go_min_score", 0.75)),
            go_min_coverage=float(
                getattr(args, "target_go_min_coverage", 0.50)
            ),
            conditional_min_score=float(
                getattr(args, "target_conditional_min_score", 0.50)
            ),
        )
        ctx["integrated_target_priority"] = integrated
        ctx["integrated_target_priority_path"] = integrated_path
        write_target_priority_summary(
            integrated,
            _integration_dir(workdir)
            / "integrated_target_priority_summary.json",
            evidence_hub_used=(
                str(
                    (
                        ctx.get("evidence_hub_summary")
                        or _read_json(
                            _integration_dir(workdir)
                            / "evidence_hub"
                            / "evidence_hub_summary.json"
                        )
                    ).get("status", "")
                )
                == "completed"
            ),
            evidence_source_count=int(
                (ctx.get("evidence_hub_summary") or {}).get(
                    "source_count",
                    0,
                )
                or 0
            ),
            evidence_targets_queried=int(
                (
                    ctx.get("evidence_hub_summary")
                    or _read_json(
                        _integration_dir(workdir)
                        / "evidence_hub"
                        / "evidence_hub_summary.json"
                    )
                ).get("targets_queried", 0)
                or 0
            ),
            benchmark=(
                ctx.get("evidence_hub_summary")
                or _read_json(
                    _integration_dir(workdir)
                    / "evidence_hub"
                    / "evidence_hub_summary.json"
                )
            ).get("benchmark")
            or {},
        )


def _stage_docking(args, workdir: Path, ctx: dict) -> None:
    if args.skip_docking or args.docking_targets <= 0:
        summary = {
            "status": "skipped",
            "reason": "docking disabled by arguments",
            "ok": 0,
            "failed": 0,
            "skipped": 0,
        }
        write_json(_integration_dir(workdir) / "docking_summary.json", summary)
        ctx["docking"] = summary
        return
    if "key_genes_path" not in ctx:
        key_genes_path = _integration_dir(workdir) / "key_genes.csv"
        if not key_genes_path.exists():
            raise IntegrationError(
                "key_genes.csv missing; run stage 02 before docking"
            )
        ctx["key_genes_path"] = key_genes_path
    if "evidence_path" not in ctx:
        evidence_path = _integration_dir(workdir) / "gene_evidence.csv"
        if not evidence_path.exists():
            raise IntegrationError(
                "gene_evidence.csv missing; run stage 03 before docking"
            )
        ctx["evidence_path"] = evidence_path
    integrated_path = (
        _integration_dir(workdir) / "integrated_target_priority.csv"
    )
    knockout_summary = _read_json(
        _integration_dir(workdir) / "knockout_summary.json"
    )
    integrated_ready = (
        not getattr(args, "skip_knockout", False)
        and bool(knockout_summary.get("knockout"))
        and integrated_path.exists()
        and integrated_path.stat().st_size > 0
    )
    priority_candidates = [integrated_path] if integrated_ready else []
    priority_candidates.append(_integration_dir(workdir) / "target_priority.csv")
    priority_path = None
    for candidate in priority_candidates:
        if candidate.exists() and candidate.stat().st_size > 0:
            priority_path = candidate
            break
    ctx["docking"] = run_docking_stage(
        workdir,
        ctx["docking_config"],
        ctx["key_genes_path"],
        ctx["evidence_path"],
        max_targets=args.docking_targets,
        ligand_library=args.ligand_library,
        force=args.force,
        priority_csv=priority_path,
        allow_review=bool(
            getattr(args, "allow_review_docking", False)
        ),
    )


def _target_dir(workdir: Path, gene: str) -> Path:
    return workdir / "work" / safe_name(str(gene), str(gene))


def _load_target_cfg(workdir: Path, gene: str, docking_config: Path):
    target = _target_dir(workdir, gene)
    config_path = target / "config" / "docking_config.json"
    if config_path.exists():
        return load_config(config_path)
    return load_config(
        docking_config,
        {"workdir": str(target), "target_name": str(gene)},
    )


def _run_target_cadd(
    gene: str,
    workdir: Path,
    docking_config: Path,
    md_settings: dict,
    handoff_settings: dict,
    docking_ml_settings: dict,
    training_csv: str | None,
) -> dict:
    """Run MD preparation, handoff and optional ML rescoring for one target."""
    rec = {
        "gene": str(gene),
        "md_status": "skipped",
        "md_mode": md_settings.get("mode", "prepare"),
        "md_requested": 0,
        "md_completed": 0,
        "md_prepared": 0,
        "md_failed": 0,
        "handoff_status": "skipped",
        "ml_status": "skipped",
        "ml_scored": 0,
        "error": "",
    }
    target = _target_dir(workdir, gene)
    try:
        cfg = _load_target_cfg(workdir, gene, docking_config)
    except Exception as exc:  # noqa: BLE001
        rec["handoff_status"] = "failed"
        rec["error"] = f"target config could not be loaded: {exc}"
        return rec

    if md_settings.get("enabled", True):
        try:
            mode = str(md_settings.get("mode") or "prepare")
            cfg.data.setdefault("md_simulation", {})["mode"] = mode
            if md_settings.get("top_n"):
                cfg.data["md_simulation"]["top_n"] = int(md_settings["top_n"])
            md_summary = md_simulation.run_md_simulation(cfg, log, mode=mode)
            rec["md_mode"] = mode
            rec["md_requested"] = int(md_summary.get("requested", 0))
            rec["md_completed"] = int(md_summary.get("completed", 0))
            rec["md_prepared"] = int(md_summary.get("prepared", 0))
            requested = rec["md_requested"]
            rec["md_failed"] = int(md_summary.get("failed", 0))
            failed = rec["md_failed"]
            if requested == 0:
                rec["md_status"] = "skipped"
                rec["error"] = "no top docking poses selected for MD"
            elif failed == requested:
                rec["md_status"] = "failed"
            elif failed > 0:
                rec["md_status"] = "partial"
            elif rec["md_completed"]:
                rec["md_status"] = "completed"
            else:
                rec["md_status"] = "prepared"
        except Exception as exc:  # noqa: BLE001
            rec["md_status"] = "failed"
            rec["md_failed"] = max(1, int(rec.get("md_failed", 0)))
            rec["error"] = f"MD simulation failed: {exc}"
            log.error("target %s MD failed: %s", gene, exc)

    if handoff_settings.get("enabled", True):
        try:
            handoff.export_md(cfg, log)
            handoff.export_external(cfg, log)
            rec["handoff_status"] = "completed"
        except Exception as exc:  # noqa: BLE001
            rec["handoff_status"] = "failed"
            if rec["error"]:
                rec["error"] += f"; handoff failed: {exc}"
            else:
                rec["error"] = f"handoff failed: {exc}"
            log.error("target %s handoff failed: %s", gene, exc)

    if docking_ml_settings.get("enabled", True):
        try:
            model_file = cfg.ml_dir() / "data" / "ml_model_info.json"
            training_value = str(training_csv or "").strip()
            train_path = (
                Path(training_value).expanduser() if training_value else None
            )
            if train_path is not None and train_path.is_file():
                docking_ml.train_ml(
                    cfg,
                    log,
                    model_type=str(
                        docking_ml_settings.get("model")
                        or cfg.get("ml", "model", "rf")
                        or "rf"
                    ),
                    label_column=str(
                        docking_ml_settings.get("label_column")
                        or cfg.get("ml", "label_column", "active")
                        or "active"
                    ),
                    training_csv=str(train_path),
                )
            if not model_file.exists():
                rec["ml_status"] = "skipped"
                rec["ml_scored"] = 0
            else:
                summary = docking_ml.predict_ml(cfg, log)
                rec["ml_status"] = "completed"
                rec["ml_scored"] = int(summary.get("scored", 0))
        except Exception as exc:  # noqa: BLE001
            rec["ml_status"] = "failed"
            if rec["error"]:
                rec["error"] += f"; ML rescoring failed: {exc}"
            else:
                rec["error"] = f"ML rescoring failed: {exc}"
            log.error("target %s ML rescoring failed: %s", gene, exc)

    write_json(
        target / "outputs" / "integration" / "cadd_target_summary.json",
        rec,
    )
    return rec


def _stage_cadd_downstream(args, workdir: Path, ctx: dict) -> None:
    docking_csv = _integration_dir(workdir) / "docking_targets.csv"
    rows: list[dict] = []
    if docking_csv.exists():
        docking = pd.read_csv(docking_csv)
        ok_genes = docking.loc[
            docking.get("status", pd.Series(dtype=str)).astype(str) == "ok",
            "gene",
        ].astype(str).tolist()
        for gene in ok_genes:
            rows.append(
                _run_target_cadd(
                    gene,
                    workdir,
                    ctx["docking_config"],
                    getattr(args, "md_simulation", DEFAULT_MD_SIMULATION),
                    getattr(args, "handoff", DEFAULT_HANDOFF),
                    getattr(args, "docking_ml", DEFAULT_DOCKING_ML),
                    getattr(args, "docking_ml_training_csv", None),
                )
            )
    if not rows:
        summary = {
            "status": "skipped",
            "reason": "no successful docking targets",
            "targets": 0,
        }
        _integration_dir(workdir).mkdir(parents=True, exist_ok=True)
        empty = pd.DataFrame(
            columns=[
                "gene",
                "md_status",
                "md_mode",
                "md_requested",
                "md_completed",
                "md_prepared",
                "md_failed",
                "handoff_status",
                "ml_status",
                "ml_scored",
                "error",
            ]
        )
        empty.to_csv(
            _integration_dir(workdir) / "cadd_targets.csv",
            index=False,
        )
    else:
        frame = pd.DataFrame(rows)
        frame.to_csv(
            _integration_dir(workdir) / "cadd_targets.csv",
            index=False,
        )
        md_failed = int((frame["md_status"] == "failed").sum())
        md_completed = int(frame["md_completed"].sum())
        md_prepared = int(frame["md_prepared"].sum())
        handoff_ok = int((frame["handoff_status"] == "completed").sum())
        handoff_failed = int((frame["handoff_status"] == "failed").sum())
        ml_failed = int((frame["ml_status"] == "failed").sum())
        ml_scored = int(frame["ml_scored"].sum())
        md_enabled = bool(
            (
                getattr(args, "md_simulation", DEFAULT_MD_SIMULATION) or {}
            ).get("enabled", True)
        )
        handoff_enabled = bool(
            (getattr(args, "handoff", DEFAULT_HANDOFF) or {}).get(
                "enabled",
                True,
            )
        )
        docking_ml_enabled = bool(
            (getattr(args, "docking_ml", DEFAULT_DOCKING_ML) or {}).get(
                "enabled",
                True,
            )
        )
        all_md_failed = (
            md_enabled
            and md_failed
            and not (md_completed or md_prepared)
        )
        all_handoff_failed = (
            handoff_enabled
            and handoff_failed == len(rows)
            and handoff_ok == 0
        )
        all_ml_failed = (
            docking_ml_enabled
            and ml_failed == len(rows)
            and len(rows) > 0
        )
        if all_md_failed or all_handoff_failed or all_ml_failed:
            status = "failed"
            details = []
            if all_md_failed:
                details.append("MD failed for all target(s)")
            if all_handoff_failed:
                details.append("handoff failed for all target(s)")
            if all_ml_failed:
                details.append("ML rescoring failed for all target(s)")
            reason = "; ".join(details)
        elif md_failed or handoff_failed or ml_failed:
            status = "partial"
            reason = (
                f"MD failed for {md_failed}, handoff failed for "
                f"{handoff_failed}, ML rescoring failed for {ml_failed} "
                f"of {len(rows)} target(s)"
            )
        else:
            status = "completed"
            reason = ""
        summary = {
            "status": status,
            "reason": reason,
            "targets": len(rows),
            "md_failed": md_failed,
            "md_completed": md_completed,
            "md_prepared": md_prepared,
            "handoff_ok": handoff_ok,
            "ml_scored": ml_scored,
            "output_csv": str(
                _integration_dir(workdir) / "cadd_targets.csv"
            ),
        }
    structural_quality = build_structural_quality(
        workdir,
        _integration_dir(workdir),
        protein_rmsd_std_cutoff_nm=float(
            getattr(
                args,
                "md_rmsd_stable_std_nm",
                0.15,
            )
            or 0.15
        ),
        ligand_rmsd_std_cutoff_nm=float(
            getattr(
                args,
                "md_rmsd_stable_std_nm",
                0.15,
            )
            or 0.15
        ),
    )
    summary["structural_quality"] = structural_quality
    write_json(
        _integration_dir(workdir) / "cadd_downstream_summary.json",
        summary,
    )
    ctx["cadd_downstream"] = summary
    if summary.get("status") == "failed":
        raise IntegrationError(summary.get("reason") or "CADD downstream failed")


def _resolved_section_paths(section: dict, base: Path) -> dict:
    section = dict(section or {})
    for key in (
        "compound_targets_csv",
        "disease_genes_csv",
        "ppi_network_csv",
        "input_csv",
    ):
        value = section.get(key)
        if value:
            section[key] = str(_resolve_path(value, base))
    targets = section.get("target_sources")
    if isinstance(targets, dict):
        section["target_sources"] = {
            str(source): str(_resolve_path(path, base))
            for source, path in targets.items()
            if path
        }
    output_dir = section.get("output_dir")
    if output_dir:
        section["output_dir"] = str(_resolve_path(output_dir, base))
    return section


def _stage_network(args, workdir: Path, ctx: dict) -> None:
    section = _resolved_section_paths(
        getattr(args, "network_toxicology", DEFAULT_NETWORK_TOXICOLOGY),
        workdir,
    )
    summary_path = _integration_dir(workdir) / "network_summary.json"
    if not section.get("enabled", True):
        write_json(
            summary_path,
            {"status": "skipped", "reason": "network toxicology disabled"},
        )
        return

    disease_csv = section.get("disease_genes_csv")
    if not disease_csv:
        key_genes = _integration_dir(workdir) / "key_genes.csv"
        if key_genes.exists():
            disease_csv = str(key_genes)
            section["disease_genes_csv"] = disease_csv
    if not section.get("compound_targets_csv") and not section.get(
        "target_sources"
    ):
        write_json(
            summary_path,
            {
                "status": "skipped",
                "reason": (
                    "no compound-target evidence provided "
                    "(network_toxicology.compound_targets_csv or "
                    "target_sources)"
                ),
            },
        )
        return

    cfg = load_config(
        ctx["docking_config"],
        {"workdir": str(workdir)},
    )
    cfg.data["network_toxicology"] = section
    try:
        result = network_toxicology.run_network_toxicology(cfg, log)
    except Exception as exc:  # noqa: BLE001
        summary = {
            "status": "failed",
            "reason": str(exc),
            "output_dir": section.get("output_dir", ""),
        }
        write_json(summary_path, summary)
        log.error("network toxicology failed: %s", exc)
        raise IntegrationError(
            "network toxicology stage failed: " + str(exc)
        ) from exc
    summary = {"status": "completed", **result}
    summary["output_dir"] = section.get("output_dir", "")
    write_json(summary_path, summary)
    ctx["network"] = summary


def _stage_faers(args, workdir: Path, ctx: dict) -> None:
    section = _resolved_section_paths(
        getattr(args, "faers", DEFAULT_FAERS),
        workdir,
    )
    summary_path = _integration_dir(workdir) / "faers_summary.json"
    if not section.get("enabled", True):
        write_json(
            summary_path,
            {"status": "skipped", "reason": "FAERS screening disabled"},
        )
        return
    if not section.get("input_csv"):
        write_json(
            summary_path,
            {
                "status": "skipped",
                "reason": (
                    "no FAERS event table provided (faers.input_csv)"
                ),
            },
        )
        return

    cfg = load_config(
        ctx["docking_config"],
        {"workdir": str(workdir)},
    )
    cfg.data["faers"] = section
    try:
        result = signal_detection.run_faers(cfg, log)
    except Exception as exc:  # noqa: BLE001
        summary = {
            "status": "failed",
            "reason": str(exc),
            "output_dir": section.get("output_dir", ""),
        }
        write_json(summary_path, summary)
        log.error("FAERS signal detection failed: %s", exc)
        raise IntegrationError(
            "FAERS signal detection stage failed: " + str(exc)
        ) from exc
    summary = {"status": "completed", **result}
    summary["output_dir"] = section.get("output_dir", "")
    write_json(summary_path, summary)
    ctx["faers"] = summary


def _stage_cell_feedback(args, workdir: Path, ctx: dict) -> None:
    if ctx.get("dataset_mode") is None:
        ctx["dataset_mode"] = _dataset_mode_from_root(ctx["single_cell_root"])
    sample_level_mode = ctx.get("dataset_mode") != "single_cell"
    if args.skip_cell_feedback or sample_level_mode:
        summary = {
            "status": "skipped",
            "reason": (
                "cell-level feedback is not applicable to sample-level datasets"
                if sample_level_mode
                else "cell feedback disabled by arguments"
            ),
        }
        write_json(
            _integration_dir(workdir) / "cell_feedback" / "cell_feedback_summary.json",
            summary,
        )
        ctx["cell_feedback"] = summary
        return
    ctx["cell_feedback"] = cell_feedback.run_cell_feedback(
        workdir,
        ctx["single_cell_root"],
        top_n=args.feedback_top_n,
        max_features=args.feedback_max_features,
        timeout_seconds=args.feedback_timeout,
        species=_resolve_feedback_species(args, ctx),
    )


def _resolve_feedback_species(args, ctx: dict) -> str:
    species = str(getattr(args, "species", "auto") or "auto").strip().lower()
    if species in ("hs", "mm"):
        return species
    root = Path(ctx["single_cell_root"])
    accession = str(getattr(args, "accession", "") or "").strip().upper()
    if accession:
        manifest_path = root / "data" / f"{accession}_manifest.json"
        if manifest_path.exists():
            try:
                organism = str(
                    json.loads(
                        manifest_path.read_text(encoding="utf-8")
                    ).get("organism", "")
                ).lower()
            except (OSError, ValueError) as exc:
                log.warning("could not read manifest %s: %s", manifest_path, exc)
                organism = ""
            if organism in ("hs", "mm"):
                return organism
            if "mus musculus" in organism or "mouse" in organism:
                return "mm"
    return "hs"


def _stage_report(args, workdir: Path, ctx: dict) -> None:
    ctx["report"] = generate_integrated_report(
        workdir,
        ctx["single_cell_root"],
        ctx["docking_config"],
        ctx,
    )


def _dry_run_stages(args, workdir: Path, ctx: dict) -> int:
    """Show what the pipeline would run without executing anything."""
    print(f"full pipeline dry run (workdir: {workdir})")
    for code, name, description in STAGES:
        if args.start_stage and code < args.start_stage:
            print(f"{code} {name:<20} SKIP  {description} (before --start-stage)")
            continue
        signature = _stage_signature(code, args, workdir, ctx)
        outputs_ready = _stage_outputs_ready(code, workdir, ctx, args)
        outdated, reason = _stage_outdated(
            workdir,
            code,
            name,
            signature,
            outputs_ready,
        )
        if args.force or outdated:
            print(f"{code} {name:<20} RUN   {description} ({reason})")
        else:
            print(f"{code} {name:<20} DONE  {description} (marker up to date)")
    return 0


def _stage_outdated(
    workdir: Path,
    code: str,
    name: str,
    signature: str,
    outputs_ready: bool,
) -> tuple[bool, str]:
    marker = _read_stage_marker(workdir, code, name)
    if marker is None:
        return True, "no marker"
    if not outputs_ready:
        return True, "required outputs missing"
    if marker.get("signature") != signature:
        return True, "inputs or parameters changed"
    return False, ""


def run_full_pipeline(args) -> int:
    """Run the integrated pipeline with provenance-aware resume markers."""
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    cfg_path = Path(args.docking_config).resolve()
    ctx = {
        "single_cell_root": Path(args.output).resolve(),
        "workdir": workdir,
        "full_config": Path(args.config).resolve(),
        "docking_config": cfg_path,
        "external_validation_path": getattr(
            args,
            "external_validation_path",
            None,
        ),
        "external_validation_target_column": getattr(
            args,
            "external_validation_target_column",
            "gene",
        ),
        "external_validation_score_column": getattr(
            args,
            "external_validation_score_column",
            "score",
        ),
        "external_validation_label_column": getattr(
            args,
            "external_validation_label_column",
            "label",
        ),
        "external_validation_threshold": getattr(
            args,
            "external_validation_threshold",
            0.5,
        ),
        "external_validation_bootstrap": getattr(
            args,
            "external_validation_bootstrap",
            1000,
        ),
    }
    if not args.force:
        _invalidate_markers_for_changed_root(workdir, ctx["single_cell_root"])
    if args.dry_run:
        return _dry_run_stages(args, workdir, ctx)
    _prune_stale_stage_markers(workdir)

    stage_fns = {
        "01": _stage_single_cell,
        "02": _stage_key_targets,
        "03": _stage_evidence,
        "04": _stage_knockout_inputs,
        "05": _stage_knockout,
        "06": _stage_docking,
        "07": _stage_cadd_downstream,
        "08": _stage_network,
        "09": _stage_faers,
        "10": _stage_cell_feedback,
        "11": _stage_report,
    }

    for index, (code, name, _description) in enumerate(STAGES):
        signature = _stage_signature(code, args, workdir, ctx)
        marker = _marker(workdir, code, name)
        if args.start_stage and code < args.start_stage:
            _write_stage_marker(
                workdir,
                code,
                name,
                signature,
                note="skipped by start-stage",
            )
            log.info("stage %s %s skipped by --start-stage", code, name)
            continue

        outputs_ready = _stage_outputs_ready(code, workdir, ctx, args)
        if not args.force:
            outdated, reason = _stage_outdated(
                workdir,
                code,
                name,
                signature,
                outputs_ready,
            )
            if not outdated:
                log.info("skip stage %s %s (already done)", code, name)
                continue
            if reason == "no marker":
                log.info("stage %s %s not run yet; executing", code, name)
            else:
                log.warning(
                    "stage %s %s outdated: %s; rerunning",
                    code,
                    name,
                    reason,
                )
            marker.unlink(missing_ok=True)

        log.info("=== stage %s %s ===", code, name)
        try:
            stage_fns[code](args, workdir, ctx)
        except PauseRequested as exc:
            log.info(str(exc))
            return 98
        _verify_stage_outputs(code, workdir, ctx, args)
        _write_stage_marker(workdir, code, name, signature)
        _clear_downstream_markers(workdir, from_index=index + 1)
        if code == "01":
            _write_run_context(workdir, ctx["single_cell_root"])
        log.info("stage %s %s complete", code, name)
    _write_run_context(workdir, ctx["single_cell_root"])
    log.info("full pipeline complete: %s", _integration_dir(workdir))
    return 0


def load_full_config(path: Path) -> dict:
    defaults = {
        "accession": "GSE125449",
        "single_cell_output": "",
        "workdir": "",
        "species": "auto",
        "top_genes": 50,
        "candidate_universe_size": 0,
        "docking_targets": 3,
        "ml_model": "xgb",
        "keep_all_genes": False,
        "case_label": None,
        "normal_label": None,
        "ligand_library": None,
        "ko_top_n": None,
        "depmap_csv": None,
        "ppi_network_csv": None,
        "md_simulation": DEFAULT_MD_SIMULATION,
        "handoff": DEFAULT_HANDOFF,
        "docking_ml": DEFAULT_DOCKING_ML,
        "docking_ml_training_csv": None,
        "network_toxicology": DEFAULT_NETWORK_TOXICOLOGY,
        "faers": DEFAULT_FAERS,
        "cell_feedback": {
            "enabled": True,
            "top_n": 12,
            "max_features": 8,
            "timeout_seconds": 3600,
        },
        "evidence": DEFAULT_EVIDENCE,
        "target_priority": DEFAULT_TARGET_PRIORITY,
        "docking_selection": DEFAULT_DOCKING_SELECTION,
        "external_validation": DEFAULT_EXTERNAL_VALIDATION,
        "qc_gate": DEFAULT_QC_GATE,
        "differential_abundance": DEFAULT_DIFFERENTIAL_ABUNDANCE,
        "gene_blacklist": DEFAULT_GENE_BLACKLIST,
    }
    raw = {}
    if path.exists():
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise IntegrationError(
                f"config file is not valid JSON: {path}: {exc}"
            ) from exc
        if not isinstance(raw, dict):
            raise IntegrationError(
                f"config file must contain a JSON object: {path}"
            )
    config = dict(defaults)
    config.update(raw or {})
    config["cell_feedback"] = dict(defaults["cell_feedback"])
    config["cell_feedback"].update((raw.get("cell_feedback") or {}))
    config["evidence"] = dict(defaults["evidence"])
    config["evidence"].update((raw.get("evidence") or {}))
    target_priority = dict(defaults["target_priority"])
    target_priority_weights = dict(
        defaults["target_priority"].get("weights") or {}
    )
    target_priority_weights.update(
        ((raw.get("target_priority") or {}).get("weights") or {})
    )
    target_priority.update(raw.get("target_priority") or {})
    target_priority["weights"] = target_priority_weights
    config["target_priority"] = target_priority
    docking_selection = dict(defaults["docking_selection"])
    docking_selection.update((raw.get("docking_selection") or {}))
    config["docking_selection"] = docking_selection
    external_validation = dict(defaults["external_validation"])
    external_validation.update((raw.get("external_validation") or {}))
    config["external_validation"] = external_validation
    config["qc_gate"] = dict(defaults["qc_gate"])
    config["qc_gate"].update((raw.get("qc_gate") or {}))
    config["differential_abundance"] = dict(
        defaults["differential_abundance"]
    )
    config["differential_abundance"].update(
        (raw.get("differential_abundance") or {})
    )
    config["gene_blacklist"] = (
        raw.get("gene_blacklist") or defaults["gene_blacklist"]
    )
    for key, default in [
        ("md_simulation", DEFAULT_MD_SIMULATION),
        ("handoff", DEFAULT_HANDOFF),
        ("docking_ml", DEFAULT_DOCKING_ML),
        ("network_toxicology", DEFAULT_NETWORK_TOXICOLOGY),
        ("faers", DEFAULT_FAERS),
    ]:
        merged = dict(default)
        merged.update((raw.get(key) or {}))
        config[key] = merged
    return config


def _apply_defaults(args, config: dict) -> None:
    if args.accession is None:
        args.accession = config.get("accession", "GSE125449")
    if args.output is None:
        args.output = config.get("single_cell_output", "../liver_cancer")
    if args.workdir is None:
        args.workdir = config.get("workdir", "")
    if args.species is None:
        args.species = config.get("species", "auto")
    if args.top_genes is None:
        args.top_genes = int(config.get("top_genes", 50))
    if getattr(args, "candidate_universe_size", None) is None:
        args.candidate_universe_size = int(
            config.get("candidate_universe_size", 0)
        )
    if getattr(args, "ml_model", None) is None:
        args.ml_model = config.get("ml_model", "xgb")
    if args.docking_targets is None:
        args.docking_targets = int(config.get("docking_targets", 3))
    if args.ligand_library is None:
        args.ligand_library = config.get("ligand_library")
    if args.case_label is None:
        args.case_label = config.get("case_label")
    if args.normal_label is None:
        args.normal_label = config.get("normal_label")
    if args.ko_top_n is None:
        args.ko_top_n = config.get("ko_top_n")
    if args.depmap_csv is None:
        args.depmap_csv = config.get("depmap_csv")
    if args.ppi_network_csv is None:
        args.ppi_network_csv = config.get("ppi_network_csv")
    if getattr(args, "advanced_priority_csv", None) is None:
        args.advanced_priority_csv = config.get("advanced_priority_csv")

    if getattr(args, "skip_md", None) is None:
        args.skip_md = not bool(
            config.get("md_simulation", {}).get("enabled", True)
        )
    if getattr(args, "md_mode", None) is None:
        args.md_mode = config.get("md_simulation", {}).get("mode", "prepare")
    if getattr(args, "md_top_n", None) is None:
        args.md_top_n = config.get("md_simulation", {}).get("top_n", 1)
    args.md_simulation = {
        "enabled": not bool(args.skip_md),
        "mode": str(args.md_mode or "prepare"),
        "top_n": int(args.md_top_n or 1),
    }
    if getattr(args, "skip_handoff", None) is None:
        args.skip_handoff = not bool(
            config.get("handoff", {}).get("enabled", True)
        )
    args.handoff = {"enabled": not bool(args.skip_handoff)}
    if getattr(args, "skip_docking_ml", None) is None:
        args.skip_docking_ml = not bool(
            config.get("docking_ml", {}).get("enabled", True)
        )
    if getattr(args, "docking_ml_model", None) is None:
        args.docking_ml_model = config.get("docking_ml", {}).get("model", "rf")
    if getattr(args, "docking_ml_label_column", None) is None:
        args.docking_ml_label_column = config.get("docking_ml", {}).get(
            "label_column",
            "active",
        )
    if getattr(args, "docking_ml_training_csv", None) is None:
        args.docking_ml_training_csv = config.get(
            "docking_ml_training_csv"
        ) or config.get("docking_ml", {}).get("training_csv")
    args.docking_ml = {
        "enabled": not bool(args.skip_docking_ml),
        "model": str(args.docking_ml_model or "rf"),
        "label_column": str(args.docking_ml_label_column or "active"),
    }

    network_section = dict(DEFAULT_NETWORK_TOXICOLOGY)
    network_section.update(config.get("network_toxicology") or {})
    if getattr(args, "skip_network", None) is None:
        args.skip_network = not bool(network_section.get("enabled", True))
    network_section["enabled"] = not bool(args.skip_network)
    for attr, key in [
        ("network_compound_targets_csv", "compound_targets_csv"),
        ("network_target_sources_dir", "target_sources_dir"),
        ("network_disease_genes_csv", "disease_genes_csv"),
        ("network_disease_gene_column", "disease_gene_column"),
        ("network_ppi_network_csv", "ppi_network_csv"),
        ("network_output_dir", "output_dir"),
        ("network_cytoscape", "cytoscape"),
        ("network_cytoscape_url", "cytoscape_url"),
        ("network_cytoscape_layout", "cytoscape_layout"),
        ("network_max_ppi_edges", "max_ppi_edges"),
        ("network_run_enrichment", "run_enrichment"),
        ("network_enrichment_timeout", "enrichment_timeout"),
    ]:
        value = getattr(args, attr, None)
        if value is None:
            value = network_section.get(key)
        if value is not None:
            network_section[key] = value
    if getattr(args, "network_cytoscape_session", None) is not None:
        network_section["cytoscape_save_session"] = bool(
            args.network_cytoscape_session
        )
    args.network_toxicology = network_section

    faers_section = dict(DEFAULT_FAERS)
    faers_section.update(config.get("faers") or {})
    if getattr(args, "skip_faers", None) is None:
        args.skip_faers = not bool(faers_section.get("enabled", True))
    faers_section["enabled"] = not bool(args.skip_faers)
    for attr, key in [
        ("faers_input", "input_csv"),
        ("faers_drug_column", "drug_column"),
        ("faers_event_column", "event_column"),
        ("faers_count_column", "count_column"),
        ("faers_min_count", "min_count"),
    ]:
        value = getattr(args, attr, None)
        if value is None:
            value = faers_section.get(key)
        if value is not None:
            faers_section[key] = value
    args.faers = faers_section

    if args.skip_cell_feedback is None:
        args.skip_cell_feedback = not bool(
            config.get("cell_feedback", {}).get("enabled", True)
        )
    if args.feedback_top_n is None:
        args.feedback_top_n = int(
            config.get("cell_feedback", {}).get("top_n", 12)
        )
    if args.feedback_max_features is None:
        args.feedback_max_features = int(
            config.get("cell_feedback", {}).get("max_features", 8)
        )
    if args.feedback_timeout is None:
        args.feedback_timeout = int(
            config.get("cell_feedback", {}).get("timeout_seconds", 3600)
        )
    if args.keep_all_genes is None:
        args.keep_all_genes = bool(config.get("keep_all_genes", False))
    if args.skip_evidence_fetch is None:
        args.skip_evidence_fetch = not bool(
            config.get("evidence", {}).get("fetch", True)
        )
    evidence_section = dict(DEFAULT_EVIDENCE)
    evidence_section.update(config.get("evidence") or {})
    if getattr(args, "skip_evidence_hub", False):
        args.evidence_hub_enabled = False
    if getattr(args, "evidence_hub_offline", False):
        args.evidence_hub_allow_network = False
    evidence_defaults = {
        "evidence_hub_enabled": "hub_enabled",
        "evidence_hub_config": "hub_config",
        "evidence_disease": "disease_name",
        "evidence_disease_id": "disease_id",
        "evidence_max_targets": "max_targets",
        "evidence_max_records": "max_records_per_source",
        "evidence_hub_timeout": "hub_timeout",
        "evidence_hub_allow_network": "allow_network",
        "evidence_hub_strict": "strict",
        "evidence_legacy_pool_size": "legacy_pool_size",
        "candidate_expansion_max_targets": "candidate_expansion_max_targets",
        "benchmark_positive_targets": "benchmark_positive_targets",
        "benchmark_negative_targets": "benchmark_negative_targets",
        "benchmark_top_n": "benchmark_top_n",
    }
    for attr, key in evidence_defaults.items():
        if getattr(args, attr, None) is None:
            setattr(args, attr, evidence_section.get(key))
    target_priority = dict(DEFAULT_TARGET_PRIORITY)
    target_priority.update(config.get("target_priority") or {})
    target_weights = dict(DEFAULT_TARGET_PRIORITY["weights"])
    target_weights.update(
        (config.get("target_priority") or {}).get("weights") or {}
    )
    args.target_priority_weights = target_weights
    args.target_go_min_score = float(
        target_priority.get("go_min_score", 0.75)
    )
    args.target_go_min_coverage = float(
        target_priority.get("go_min_coverage", 0.50)
    )
    args.target_conditional_min_score = float(
        target_priority.get("conditional_min_score", 0.50)
    )
    if getattr(args, "allow_review_docking", None) is None:
        args.allow_review_docking = bool(
            (config.get("docking_selection") or {}).get(
                "allow_review",
                False,
            )
        )
    external_validation = dict(DEFAULT_EXTERNAL_VALIDATION)
    external_validation.update(config.get("external_validation") or {})
    for attr in (
        "external_validation_path",
        "external_validation_target_column",
        "external_validation_score_column",
        "external_validation_label_column",
    ):
        if getattr(args, attr, None) is None:
            key = attr.replace("external_validation_", "")
            setattr(args, attr, external_validation.get(key))
    if getattr(args, "external_validation_threshold", None) is None:
        args.external_validation_threshold = float(
            external_validation.get("threshold", 0.5)
        )
    if getattr(args, "external_validation_bootstrap", None) is None:
        args.external_validation_bootstrap = int(
            external_validation.get("bootstrap", 1000)
        )
    if args.evidence_workers is None:
        args.evidence_workers = int(config.get("evidence", {}).get("max_workers", 6))
    if args.evidence_timeout is None:
        args.evidence_timeout = int(config.get("evidence", {}).get("timeout", 90))
    qc_gate = dict(DEFAULT_QC_GATE)
    qc_gate.update(config.get("qc_gate") or {})
    if getattr(args, "skip_qc_gate", False):
        qc_gate["enabled"] = False
    args.qc_gate = qc_gate
    differential_abundance = dict(DEFAULT_DIFFERENTIAL_ABUNDANCE)
    differential_abundance.update(config.get("differential_abundance") or {})
    if getattr(args, "skip_differential_abundance", False):
        differential_abundance["enabled"] = False
    args.differential_abundance = differential_abundance
    if not str(args.output or "").strip():
        raise IntegrationError(
            "single-cell output is required; provide --output or set "
            "single_cell_output in config"
        )
    if not str(args.workdir or "").strip():
        raise IntegrationError(
            "workdir is required; provide --workdir or set workdir in config"
        )
    args.accession = canonical_accession(
        str(args.accession or "").strip().upper()
    )
    if not re.fullmatch(
        r"(?:GSE\d+|E-[A-Z0-9]+-\d+|S-BSST\d+)",
        args.accession,
    ):
        raise IntegrationError(
            "dataset accession must look like GSE125449, E-MTAB-1234, "
            "or S-BSST123"
        )
    if args.workdir:
        args.workdir = str(_resolve_path(args.workdir, APP_ROOT))
    if args.output:
        args.output = str(_resolve_path(args.output, APP_ROOT))
    if args.ligand_library:
        args.ligand_library = str(_resolve_path(args.ligand_library, Path.cwd()))
    if args.depmap_csv:
        args.depmap_csv = str(_resolve_path(args.depmap_csv, Path.cwd()))
    if args.ppi_network_csv:
        args.ppi_network_csv = str(
            _resolve_path(args.ppi_network_csv, Path.cwd())
        )
    if args.advanced_priority_csv:
        args.advanced_priority_csv = str(
            _resolve_path(args.advanced_priority_csv, Path.cwd())
        )
    if getattr(args, "evidence_hub_config", None):
        args.evidence_hub_config = str(
            _resolve_path(args.evidence_hub_config, APP_ROOT)
        )
    if getattr(args, "external_validation_path", None):
        args.external_validation_path = str(
            _resolve_path(args.external_validation_path, Path.cwd())
        )
    for attr in [
        "docking_ml_training_csv",
        "network_compound_targets_csv",
        "network_target_sources_dir",
        "network_disease_genes_csv",
        "network_ppi_network_csv",
        "faers_input",
    ]:
        value = getattr(args, attr, None)
        if value:
            setattr(args, attr, str(_resolve_path(value, Path.cwd())))
    section_updates = {
        "network_toxicology": {
            "compound_targets_csv": "network_compound_targets_csv",
            "target_sources_dir": "network_target_sources_dir",
            "disease_genes_csv": "network_disease_genes_csv",
            "ppi_network_csv": "network_ppi_network_csv",
        },
        "faers": {"input_csv": "faers_input"},
    }
    for section_name, mapping in section_updates.items():
        section = getattr(args, section_name, {})
        if not isinstance(section, dict):
            continue
        for section_key, attr_name in mapping.items():
            value = getattr(args, attr_name, None)
            if value is not None:
                section[section_key] = value
    config_path_keys = {
        "network_toxicology": [
            "compound_targets_csv",
            "target_sources_dir",
            "disease_genes_csv",
            "ppi_network_csv",
        ],
        "faers": ["input_csv"],
    }
    for section_name, keys in config_path_keys.items():
        section = getattr(args, section_name, {})
        if not isinstance(section, dict):
            continue
        for key in keys:
            value = section.get(key)
            if value and not Path(str(value)).is_absolute():
                section[key] = str(_resolve_path(value, Path.cwd()))
        targets = section.get("target_sources")
        if isinstance(targets, dict):
            section["target_sources"] = {
                str(source): str(_resolve_path(path, Path.cwd()))
                for source, path in targets.items()
                if path
            }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_full_pipeline",
        description="Automated scRNA-seq -> key targets -> virtual screening -> knockout",
    )
    parser.add_argument(
        "--config",
        default=str(APP_ROOT / "config" / "full_pipeline_config.json"),
        help="full pipeline config JSON",
    )
    parser.add_argument(
        "--accession",
        help=(
            "dataset accession: GSE125449, E-MTAB-1234 or S-BSST123 "
            "(default from config)"
        ),
    )
    parser.add_argument("--output", help="single-cell output root")
    parser.add_argument("--workdir", help="docking/integration work root")
    parser.add_argument(
        "--docking-config",
        default=str(APP_ROOT / "config" / "docking_config.json"),
        help="base docking config JSON",
    )
    parser.add_argument("--species", choices=["hs", "mm", "auto"])
    parser.add_argument("--top-genes", type=int, help="number of key genes to keep")
    parser.add_argument(
        "--candidate-universe-size",
        type=int,
        default=None,
        help="maximum DEG candidates kept before evidence ranking; 0 keeps all",
    )
    parser.add_argument(
        "--ml-model",
        choices=["xgb", "rf", "gbm", "mlp", "lasso_svm"],
        default=None,
        help="single-cell ML model used by the integrated pipeline",
    )
    parser.add_argument("--docking-targets", type=int, help="max genes to dock")
    parser.add_argument(
        "--allow-review-docking",
        action="store_true",
        default=None,
        help="allow REVIEW targets to enter docking; default is GO/CONDITIONAL_GO only",
    )
    parser.add_argument("--external-validation-path", default=None)
    parser.add_argument("--external-validation-target-column", default=None)
    parser.add_argument("--external-validation-score-column", default=None)
    parser.add_argument("--external-validation-label-column", default=None)
    parser.add_argument(
        "--external-validation-threshold",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--external-validation-bootstrap",
        type=int,
        default=None,
    )
    parser.add_argument("--ligand-library", help="ligand library file (.smi/.sdf/.csv)")
    parser.add_argument("--case-label", help="case group label for knockout")
    parser.add_argument("--normal-label", help="normal group label for knockout")
    parser.add_argument("--ko-top-n", type=int, help="top N knockout report genes")
    parser.add_argument("--depmap-csv", help="DepMap gene effect CSV")
    parser.add_argument(
        "--ppi-network-csv",
        help="STRING-style PPI edge table for knockout PPI hub scoring",
    )
    parser.add_argument(
        "--advanced-priority-csv",
        default=None,
        help=(
            "integrated_priority.csv from the advanced analysis; "
            "used to order key genes"
        ),
    )
    parser.add_argument("--skip-md", action="store_true", default=None)
    parser.add_argument(
        "--md-mode",
        choices=["prepare", "auto"],
        default=None,
        help="GROMACS MD mode used by the integrated pipeline (default: prepare)",
    )
    parser.add_argument("--md-top-n", type=int, default=None)
    parser.add_argument("--skip-handoff", action="store_true", default=None)
    parser.add_argument("--skip-docking-ml", action="store_true", default=None)
    parser.add_argument(
        "--docking-ml-model",
        choices=["rf", "gbm", "mlp", "lasso_svm", "torch"],
        default=None,
        help="docking rescoring model when a labeled training CSV is provided",
    )
    parser.add_argument("--docking-ml-label-column", default=None)
    parser.add_argument("--docking-ml-training-csv", default=None)
    parser.add_argument("--skip-network", action="store_true", default=None)
    parser.add_argument("--network-compound-targets-csv", default=None)
    parser.add_argument("--network-target-sources-dir", default=None)
    parser.add_argument("--network-disease-genes-csv", default=None)
    parser.add_argument("--network-disease-gene-column", default=None)
    parser.add_argument("--network-ppi-network-csv", default=None)
    parser.add_argument("--network-output-dir", default=None)
    parser.add_argument(
        "--network-cytoscape",
        choices=["auto", "on", "off"],
        default=None,
        help="Cytoscape live export mode for the network stage (default: auto)",
    )
    parser.add_argument("--network-cytoscape-url", default=None)
    parser.add_argument("--network-cytoscape-layout", default=None)
    parser.add_argument(
        "--network-cytoscape-session",
        action="store_true",
        default=None,
        help="also save a .cys session when Cytoscape live export runs",
    )
    parser.add_argument("--network-max-ppi-edges", type=int, default=None)
    parser.add_argument(
        "--network-run-enrichment",
        action="store_true",
        default=None,
        help="run GO/KEGG enrichment on the network overlap genes",
    )
    parser.add_argument(
        "--network-enrichment-timeout",
        type=int,
        default=None,
        help="timeout in seconds for network GO/KEGG enrichment",
    )
    parser.add_argument("--skip-faers", action="store_true", default=None)
    parser.add_argument("--faers-input", default=None)
    parser.add_argument("--faers-drug-column", default=None)
    parser.add_argument("--faers-event-column", default=None)
    parser.add_argument("--faers-count-column", default=None)
    parser.add_argument("--faers-min-count", type=int, default=None)
    parser.add_argument("--evidence-workers", type=int, default=None)
    parser.add_argument("--evidence-timeout", type=int, default=None)
    parser.add_argument(
        "--skip-evidence-hub",
        action="store_true",
        default=None,
        help="skip the multi-source evidence hub while keeping legacy target evidence",
    )
    parser.add_argument("--evidence-hub-config", default=None)
    parser.add_argument("--evidence-disease", default=None)
    parser.add_argument("--evidence-disease-id", default=None)
    parser.add_argument("--evidence-max-targets", type=int, default=None)
    parser.add_argument("--evidence-max-records", type=int, default=None)
    parser.add_argument("--evidence-hub-timeout", type=int, default=None)
    parser.add_argument(
        "--evidence-hub-strict",
        action="store_true",
        default=None,
        help="fail the evidence stage when a configured source fails",
    )
    parser.add_argument(
        "--evidence-hub-offline",
        action="store_true",
        default=None,
        help="use only cached or local evidence sources",
    )
    parser.add_argument("--evidence-legacy-pool-size", type=int, default=None)
    parser.add_argument(
        "--candidate-expansion-max-targets",
        type=int,
        default=None,
        help="maximum non-DEG disease/genetic targets added to the candidate universe",
    )
    parser.add_argument(
        "--benchmark-positive",
        dest="benchmark_positive_targets",
        default=None,
        help="comma-separated positive control target symbols",
    )
    parser.add_argument(
        "--benchmark-negative",
        dest="benchmark_negative_targets",
        default=None,
        help="comma-separated negative control target symbols",
    )
    parser.add_argument("--benchmark-top-n", type=int, default=None)
    parser.add_argument("--skip-scrna", action="store_true")
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-deps", action="store_true")
    parser.add_argument("--skip-evidence-fetch", action="store_true", default=None)
    parser.add_argument("--skip-pseudobulk", action="store_true")
    parser.add_argument("--skip-knockout", action="store_true")
    parser.add_argument("--skip-docking", action="store_true")
    parser.add_argument("--skip-cell-feedback", action="store_true", default=None)
    parser.add_argument("--skip-qc-gate", action="store_true", default=None)
    parser.add_argument(
        "--skip-differential-abundance",
        action="store_true",
        default=None,
    )
    parser.add_argument("--feedback-top-n", type=int, default=None)
    parser.add_argument("--feedback-max-features", type=int, default=None)
    parser.add_argument("--feedback-timeout", type=int, default=None)
    parser.add_argument("--keep-all-genes", action="store_true", default=None)
    parser.add_argument("--force", action="store_true", help="rerun stages from scratch")
    parser.add_argument("--start-stage", default=None, help="stage code to start from")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show which stages would run without executing them",
    )
    parser.add_argument("--list-stages", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.list_stages:
        for code, name, description in STAGES:
            print(f"{code}  {name:<20} {description}")
        return 0

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    try:
        config = load_full_config(Path(args.config))
        _apply_defaults(args, config)
        return run_full_pipeline(args)
    except (IntegrationError, DockingError, ToolNotFoundError) as exc:
        log.error("ERROR: %s", exc)
        return 1
    except (OSError, ValueError) as exc:
        log.error("ERROR: invalid configuration or JSON input: %s", exc)
        if getattr(args, "verbose", False):
            import traceback

            traceback.print_exc()
        return 1
    except KeyboardInterrupt:
        log.info("interrupted")
        return 130


if __name__ == "__main__":
    sys.exit(main())
