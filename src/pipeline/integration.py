#!/usr/bin/env python3
"""End-to-end automation from single-cell analysis to virtual screening.

The module wires the existing pieces together:

1. run the GEO single-cell pipeline and export a sample-level pseudobulk matrix;
2. rank significant DEGs into a compact key-gene table;
3. enrich genes with UniProt/PDB/ChEMBL/STRING/Reactome/Open Targets/KEGG evidence (network optional, cached);
4. build the virtual-knockout inputs and run multidimensional target scoring;
5. for genes with a PDB structure, collect known ligands and run the full
   AutoDock Vina pipeline in an isolated per-target workdir;
6. export the wet-lab validation plan and an integrated HTML report.

Every stage writes a marker file so a rerun resumes where it stopped.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import logging
import math
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

from docking import box, evidence as evidence_mod, pipeline as docking_pipeline  # noqa: E402
from docking.config import load_config, save_config  # noqa: E402
from docking.knockout import run_knockout  # noqa: E402
from docking.provenance import write_run_manifest  # noqa: E402
from docking.utils import DockingError, ToolNotFoundError, safe_name, write_json  # noqa: E402
from docking.validation import export_validation  # noqa: E402
from data.geo_downloader import canonical_accession  # noqa: E402

from . import cell_feedback, orchestrator  # noqa: E402

log = logging.getLogger("full_pipeline")

STAGES = [
    ("01", "single_cell", "expression analysis (download, QC, annotation, DEG)"),
    ("02", "key_targets", "extract and rank key genes/proteins from DEGs"),
    ("03", "evidence", "enrich genes with UniProt/PDB/ChEMBL/STRING/Reactome/Open Targets/KEGG evidence"),
    ("04", "knockout_inputs", "build pseudobulk expression and knockout inputs"),
    ("05", "knockout", "virtual knockout and multidimensional target scoring"),
    ("06", "docking", "per-target virtual screening with AutoDock Vina"),
    ("07", "cell_feedback", "re-score single-cell targets from knockout/docking results"),
    ("08", "report", "integrated HTML report and provenance manifest"),
]

# Required outputs per stage. Paths are relative to the full-pipeline workdir.
STAGE_OUTPUTS = {
    "01": ("results/pipeline_complete.json",),
    "02": (
        "outputs/integration/key_genes.csv",
        "outputs/integration/key_genes_summary.json",
    ),
    "03": ("outputs/integration/gene_evidence.csv",),
    "04": (
        "data/knockout/expression.csv",
        "data/knockout/metadata.csv",
        "data/knockout/inputs_summary.json",
    ),
    "05": (
        "outputs/integration/knockout_summary.json",
        "outputs/run_001/results/04_knockout/data/fig_52_53_ranked_knockout.csv",
        "outputs/run_001/results/05_validation/data/validation_plan.md",
    ),
    "06": (
        "outputs/integration/docking_summary.json",
        "outputs/integration/docking_targets.csv",
    ),
    "07": ("outputs/integration/cell_feedback/cell_feedback_summary.json",),
    "08": (
        "outputs/integration/integration_report.html",
        "outputs/integration/integration_summary.json",
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

DEFAULT_GENE_BLACKLIST = [
    r"^RPL",
    r"^RPS",
    r"^MRPL",
    r"^MRPS",
    r"^MT-",
    r"^MTRNR",
    r"^SNORD",
    r"^SCGB",
    r"^IGH",
    r"^IGK",
    r"^IGL",
    r"^TRA",
    r"^TRB",
    r"^TRG",
    r"^HLA-D",
    r"^LINC",
    r"^RP[0-9]",
    r"^AC[0-9]",
    r"^AL[0-9]",
]

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


class IntegrationError(RuntimeError):
    """Raised when the integrated pipeline cannot continue."""


class PauseRequested(Exception):
    """Raised when the underlying single-cell pipeline pauses."""


def _integration_dir(workdir: Path) -> Path:
    return workdir / "outputs" / "integration"


def _stage_dir(workdir: Path) -> Path:
    return _integration_dir(workdir) / ".stages"


def _marker(workdir: Path, code: str, name: str) -> Path:
    return _stage_dir(workdir) / f"{code}_{name}.done"


def _sha256_file(path: Path) -> str:
    try:
        digest = hashlib.sha256()
        with path.open("rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return "missing"


def _json_sorted(value) -> str:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            ensure_ascii=False,
            default=str,
        )
    except Exception:
        return str(value)


def _read_stage_marker(workdir: Path, code: str, name: str) -> dict | None:
    marker = _marker(workdir, code, name)
    if not marker.exists():
        return None
    try:
        data = json.loads(marker.read_text(encoding="utf-8"))
        if isinstance(data, dict) and data.get("signature"):
            return data
    except Exception:
        pass
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
                "keep_all_genes": bool(getattr(args, "keep_all_genes", False)),
                "gene_blacklist": sorted(
                    getattr(args, "gene_blacklist", DEFAULT_GENE_BLACKLIST) or []
                ),
            }
        )
    elif code == "03":
        key_genes = ctx.get("key_genes_path") or integration / "key_genes.csv"
        payload.update(
            {
                "key_genes_csv": str(key_genes),
                "key_genes_sha256": _sha256_file(Path(str(key_genes))),
                "fetch": bool(not getattr(args, "skip_evidence_fetch", False)),
                "max_workers": int(getattr(args, "evidence_workers", 6) or 6),
                "timeout": int(getattr(args, "evidence_timeout", 90) or 90),
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
        payload.update(
            {
                "docking_config": str(docking_config),
                "docking_config_sha256": _sha256_file(docking_config),
                "case_label": getattr(args, "case_label", None),
                "normal_label": getattr(args, "normal_label", None),
                "ko_top_n": getattr(args, "ko_top_n", None),
                "depmap_csv": getattr(args, "depmap_csv", None),
                "skip_knockout": bool(getattr(args, "skip_knockout", False)),
            }
        )
    elif code == "06":
        key_genes = ctx.get("key_genes_path") or integration / "key_genes.csv"
        evidence = ctx.get("evidence_path") or integration / "gene_evidence.csv"
        payload.update(
            {
                "docking_config": str(docking_config),
                "docking_config_sha256": _sha256_file(docking_config),
                "key_genes_csv": str(key_genes),
                "key_genes_sha256": _sha256_file(Path(str(key_genes))),
                "evidence_csv": str(evidence),
                "evidence_sha256": _sha256_file(Path(str(evidence))),
                "max_targets": int(getattr(args, "docking_targets", 3) or 3),
                "ligand_library": getattr(args, "ligand_library", None),
                "skip_docking": bool(getattr(args, "skip_docking", False)),
            }
        )
    elif code == "07":
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
    elif code == "08":
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


def _read_json(path: Path, default=None) -> dict:
    if not path.exists():
        return default or {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else default or {}
    except Exception:
        return default or {}


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
        except (OSError, ValueError):
            pass
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
    elif code == "07" and getattr(args, "skip_cell_feedback", False):
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
    except OSError:
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


def _as_float(value) -> float | None:
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def collect_qc_metrics(single_cell_root: Path, workdir: Path) -> dict:
    """Aggregate single-cell and downstream metrics for the QC gate."""
    data_dir = single_cell_root / "results" / "data"
    metrics: dict = {
        "single_cell": _read_json(single_cell_root / "results" / "summary.json"),
        "doublet_rate": None,
        "doublet_rate_by_sample": [],
        "pseudobulk_used": None,
        "pseudobulk_warning": "",
        "knockout": {},
        "docking": {},
    }

    doublet_csv = data_dir / "08_publication" / "fig_29_doublet_rate_by_sample.csv"
    if not doublet_csv.exists():
        doublet_csv = data_dir / "08_publication" / "fig_29_doublet_rate_sample.csv"
    if doublet_csv.exists():
        try:
            frame = pd.read_csv(doublet_csv)
            rates = pd.to_numeric(frame.get("doublet_rate"), errors="coerce")
            metrics["doublet_rate_by_sample"] = frame.to_dict(orient="records")
            if rates.notna().any():
                metrics["doublet_rate"] = float(rates.mean())
        except Exception:
            pass
    else:
        doublet_tbl = data_dir / "02_doublets" / "fig_02_doublet_results.csv"
        if doublet_tbl.exists():
            try:
                frame = pd.read_csv(doublet_tbl)
                if "doublet_call" in frame.columns and "sample" in frame.columns:
                    grouped = (
                        frame.groupby("sample")["doublet_call"]
                        .agg(
                            n_cells="count",
                            n_doublets=lambda s: int(
                                (s.astype(str).str.lower() == "doublet").sum()
                            ),
                        )
                        .reset_index()
                    )
                    grouped["doublet_rate"] = (
                        grouped["n_doublets"] / grouped["n_cells"]
                    )
                    metrics["doublet_rate_by_sample"] = grouped.to_dict(
                        orient="records"
                    )
                    rates = pd.to_numeric(
                        grouped["doublet_rate"], errors="coerce"
                    )
                    if rates.notna().any():
                        metrics["doublet_rate"] = float(rates.mean())
            except Exception:
                pass

    warn_path = data_dir / "pseudobulk_warning.txt"
    if warn_path.exists():
        try:
            metrics["pseudobulk_warning"] = (
                warn_path.read_text(encoding="utf-8", errors="replace").strip()
            )
            metrics["pseudobulk_used"] = False
        except OSError:
            pass

    integration = _integration_dir(workdir)
    ko_summary = _read_json(integration / "knockout_summary.json")
    metrics["knockout"] = ko_summary.get("knockout") or {}
    metrics["docking"] = _read_json(integration / "docking_summary.json")
    return metrics


def evaluate_qc_gate(metrics: dict, config: dict) -> dict:
    """Turn collected metrics into an explicit pass/warn/fail gate."""
    if not config.get("enabled", True):
        return {
            "status": "skipped",
            "checks": [],
            "summary": "QC gate disabled by configuration",
        }

    sc = metrics.get("single_cell") or {}
    checks: list[dict] = []
    fail_on_missing = bool(config.get("fail_on_missing_metrics", False))

    def check(
        name: str,
        value,
        threshold,
        message: str,
        missing_message: str,
        below_is_bad: bool = True,
    ) -> None:
        number = _as_float(threshold)
        if number is None or number == 0:
            return
        current = _as_float(value)
        if current is None:
            checks.append(
                {
                    "name": name,
                    "level": "fail" if fail_on_missing else "warn",
                    "ok": False,
                    "message": missing_message,
                }
            )
            return
        ok = current >= number if below_is_bad else current <= number
        checks.append(
            {
                "name": name,
                "level": "fail" if not ok else "pass",
                "ok": ok,
                "message": message.format(value=current, threshold=number),
            }
        )

    check(
        "min_cells_after_qc",
        sc.get("n_cells_after_qc"),
        config.get("min_cells_after_qc"),
        "cells after QC = {value:.0f} (min {threshold:.0f})",
        "n_cells_after_qc metric is missing",
    )
    check(
        "min_genes",
        sc.get("n_genes"),
        config.get("min_genes"),
        "genes = {value:.0f} (min {threshold:.0f})",
        "n_genes metric is missing",
    )
    check(
        "min_deg_genes",
        sc.get("deg_total"),
        config.get("min_deg_genes"),
        "DEGs = {value:.0f} (min {threshold:.0f})",
        "deg_total metric is missing",
    )
    check(
        "max_doublet_rate",
        metrics.get("doublet_rate"),
        config.get("max_doublet_rate"),
        "doublet rate = {value:.3f} (max {threshold:.3f})",
        "doublet rate metric is missing",
        below_is_bad=False,
    )

    if config.get("require_pseudobulk") and metrics.get("pseudobulk_used") is False:
        checks.append(
            {
                "name": "require_pseudobulk",
                "level": "fail",
                "ok": False,
                "message": metrics.get("pseudobulk_warning")
                or "DE fell back to cells-as-replicates",
            }
        )

    failed = [c for c in checks if c["level"] == "fail" and not c["ok"]]
    warned = [c for c in checks if c["level"] == "warn" and not c["ok"]]
    status = "fail" if failed else ("warn" if warned else "pass")
    return {
        "status": status,
        "checks": checks,
        "summary": (
            f"{len(failed)} failed, {len(warned)} warned"
            if checks
            else "no QC thresholds configured"
        ),
    }


def write_qc_metrics(
    workdir: Path,
    single_cell_root: Path,
    qc_config: dict,
) -> dict:
    """Collect and evaluate QC metrics, writing qc_metrics.json."""
    metrics = collect_qc_metrics(single_cell_root, workdir)
    metrics["qc_gate"] = evaluate_qc_gate(metrics, qc_config)
    out_dir = _integration_dir(workdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_json(out_dir / "qc_metrics.json", metrics)
    return metrics


def _chi2_contingency(table) -> tuple[float, float]:
    """Pearson chi-square (Yates-corrected 2x2) with p from the erfc formula."""
    matrix = np.asarray(table, dtype=float)
    if matrix.shape != (2, 2):
        return float("nan"), float("nan")
    total = float(matrix.sum())
    if total == 0 or np.any(matrix < 0):
        return float("nan"), float("nan")
    row_totals = matrix.sum(axis=1)
    col_totals = matrix.sum(axis=0)
    expected = np.outer(row_totals, col_totals) / total
    with np.errstate(divide="ignore", invalid="ignore"):
        corrected = float(
            np.nansum((np.abs(matrix - expected) - 0.5) ** 2 / expected)
        )
    chi2 = max(0.0, corrected)
    if chi2 == 0:
        return 0.0, 1.0
    p_value = math.erfc(math.sqrt(chi2 / 2.0))
    return chi2, min(max(p_value, 0.0), 1.0)


def _bh_adjust(p_values) -> list[float]:
    """Benjamini-Hochberg FDR adjustment without extra dependencies."""
    values = np.asarray(p_values, dtype=float)
    n = len(values)
    if n == 0:
        return []
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out = np.empty_like(adjusted)
    out[order] = adjusted
    return [float(v) for v in out]


def run_differential_abundance(
    single_cell_root: Path,
    out_dir: Path,
    config: dict | None = None,
) -> dict:
    """Pair DE with a cell-type composition shift test (skill-guided)."""
    cfg = config or {}
    out_dir.mkdir(parents=True, exist_ok=True)
    ann_path = (
        single_cell_root
        / "results"
        / "data"
        / "04_annotation"
        / "fig_05_16_17_cell_annotations.csv"
    )
    if not ann_path.exists():
        return {
            "status": "skipped",
            "reason": "cell annotation CSV not found",
            "output_csv": "",
        }
    try:
        ann = pd.read_csv(ann_path)
    except Exception:
        return {
            "status": "skipped",
            "reason": "cell annotation CSV unreadable",
            "output_csv": "",
        }
    if not {"celltype_annot", "condition"}.issubset(ann.columns):
        return {
            "status": "skipped",
            "reason": "annotation CSV lacks celltype_annot/condition columns",
            "output_csv": "",
        }
    conditions = [str(x) for x in ann["condition"].dropna().unique()]
    if len(conditions) < 2:
        return {
            "status": "skipped",
            "reason": "fewer than two conditions for composition test",
            "output_csv": "",
        }
    cond0, cond1 = conditions[0], conditions[1]
    counts = pd.crosstab(ann["celltype_annot"], ann["condition"])
    if cond0 not in counts.columns or cond1 not in counts.columns:
        return {
            "status": "skipped",
            "reason": "condition columns missing from crosstab",
            "output_csv": "",
        }

    total0 = int(counts[cond0].sum())
    total1 = int(counts[cond1].sum())
    min_cells = int(cfg.get("min_cells", 5))
    rows: list[dict] = []
    for celltype, row in counts.iterrows():
        a = int(row[cond0])
        b = int(row[cond1])
        if a + b < min_cells or (total0 + total1) == 0:
            continue
        table = np.array(
            [[a, total0 - a], [b, total1 - b]],
            dtype=float,
        )
        chi2, p_value = _chi2_contingency(table)
        fraction0 = a / total0 if total0 else 0.0
        fraction1 = b / total1 if total1 else 0.0
        rows.append(
            {
                "celltype": str(celltype),
                f"{cond0}_cells": a,
                f"{cond1}_cells": b,
                f"{cond0}_fraction": round(fraction0, 6),
                f"{cond1}_fraction": round(fraction1, 6),
                "n_cells": a + b,
                "chi2": None if math.isnan(chi2) else round(chi2, 6),
                "p_value": None if math.isnan(p_value) else p_value,
                "direction": (
                    "enriched_in_" + cond0
                    if fraction0 > fraction1
                    else "enriched_in_" + cond1
                ),
            }
        )

    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["p_adjust"] = _bh_adjust(frame["p_value"].fillna(1.0))
        frame["significant"] = frame["p_adjust"] < float(cfg.get("fdr", 0.05))
        frame = frame.sort_values(
            ["p_adjust", "p_value"],
            na_position="last",
        ).reset_index(drop=True)

    csv_path = out_dir / "differential_abundance.csv"
    frame.to_csv(csv_path, index=False)
    summary = {
        "status": "completed",
        "conditions": [cond0, cond1],
        "celltypes_tested": int(len(frame)),
        "significant_celltypes": int(
            frame["significant"].sum() if not frame.empty else 0
        ),
        "output_csv": str(csv_path),
    }
    write_json(out_dir / "differential_abundance_summary.json", summary)
    log.info(
        "differential abundance: %s cell types tested, %s significant",
        summary["celltypes_tested"],
        summary["significant_celltypes"],
    )
    return summary


def extract_key_genes(
    single_cell_root: Path,
    out_dir: Path,
    top_n: int = 50,
    keep_all: bool = False,
    blacklist_patterns: list[str] | None = None,
) -> pd.DataFrame:
    """Rank significant DEGs into a compact key-gene table."""
    data_dir = single_cell_root / "results" / "data"
    significant_path = data_dir / "05_deg" / "fig_09_deg_significant.csv"
    all_path = data_dir / "05_deg" / "fig_08_deg_all.csv"
    deg_path = significant_path if significant_path.exists() else all_path
    if not deg_path.exists():
        raise IntegrationError(f"DEG table not found under {data_dir}")

    frame = pd.read_csv(deg_path)
    if frame.empty and deg_path == significant_path and all_path.exists():
        log.warning(
            "significant DEG table is empty (%s); "
            "falling back to the full DEG table",
            significant_path,
        )
        deg_path = all_path
        frame = pd.read_csv(all_path)
    if frame.empty:
        raise IntegrationError(f"DEG table is empty: {deg_path}")

    rename = {
        "avg_log2FC": "avg_log2fc",
        "log2FoldChange": "avg_log2fc",
        "pvalue": "p_val",
        "padj": "p_val_adj",
    }
    frame = frame.rename(columns=rename)
    # Some DEG exports contain both avg_log2FC and log2FoldChange. After the
    # rename both become avg_log2fc, which turns column access into a DataFrame.
    frame = frame.loc[:, ~frame.columns.duplicated()]
    if "gene" not in frame.columns:
        raise IntegrationError(f"DEG table has no gene column: {deg_path}")

    if "significant" in frame.columns:
        flag = (
            frame["significant"]
            .astype(str)
            .str.strip()
            .str.upper()
            .isin(["TRUE", "1", "YES"])
        )
        frame = frame[flag]
    if "direction" in frame.columns and not keep_all:
        frame = frame[frame["direction"].astype(str).str.strip().isin(["Up", "Down"])]

    if "avg_log2fc" not in frame.columns:
        raise IntegrationError(f"DEG table has no log2FC column: {deg_path}")
    if "p_val_adj" not in frame.columns:
        frame["p_val_adj"] = np.nan

    frame = frame.copy()
    frame["gene"] = frame["gene"].astype(str)
    frame["avg_log2fc"] = pd.to_numeric(frame["avg_log2fc"], errors="coerce")
    frame["p_val_adj"] = pd.to_numeric(frame["p_val_adj"], errors="coerce")
    frame["abs_log2fc"] = frame["avg_log2fc"].abs()
    frame = frame.sort_values(
        ["p_val_adj", "abs_log2fc"],
        ascending=[True, False],
        na_position="last",
    ).reset_index(drop=True)
    frame["deg_rank"] = np.arange(1, len(frame) + 1)

    total_before = len(frame)
    if not keep_all:
        patterns = blacklist_patterns or DEFAULT_GENE_BLACKLIST
        if patterns:
            expr = re.compile("|".join(patterns), re.IGNORECASE)
            frame = frame[~frame["gene"].str.match(expr)]
    frame = frame.reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)

    ml_path = data_dir / "07_ml" / "fig_24_ml_feature_importance.csv"
    if ml_path.exists():
        try:
            ml = pd.read_csv(ml_path, index_col=0)
            ml.index = ml.index.astype(str)
            ml_values = ml.iloc[:, 0].astype(float)
            frame["ml_importance"] = frame["gene"].map(ml_values).fillna(0.0)
        except Exception:
            frame["ml_importance"] = 0.0
    else:
        frame["ml_importance"] = 0.0

    out_cols = [
        "rank",
        "gene",
        "direction",
        "avg_log2fc",
        "p_val_adj",
        "pct.1",
        "pct.2",
        "deg_rank",
        "ml_importance",
    ]
    for col in out_cols:
        if col not in frame.columns:
            frame[col] = ""
    frame = frame.head(top_n)[out_cols].reset_index(drop=True)
    frame["source"] = "DEG"

    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "key_genes.csv"
    frame.to_csv(csv_path, index=False)
    summary = {
        "deg_table": str(deg_path),
        "deg_total": int(total_before),
        "after_blacklist": int(len(frame)),
        "top_n": int(top_n),
        "keep_all": bool(keep_all),
        "output_csv": str(csv_path),
    }
    write_json(out_dir / "key_genes_summary.json", summary)
    log.info(
        "key targets: %s genes kept from %s DEGs -> %s",
        len(frame),
        total_before,
        csv_path,
    )
    return frame


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


def _mygene_info(gene: str, timeout: int) -> dict:
    payload = {
        "q": gene,
        "scopes": "symbol",
        "fields": "entrezgene,uniprot,ensembl.gene",
        "species": "human",
        "size": 5,
    }
    try:
        hits = _http_json(
            "https://mygene.info/v3/query",
            payload,
            timeout,
        )
    except Exception:
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


def _chembl_evidence(uniprot: str, timeout: int) -> tuple[str, int]:
    if not uniprot:
        return "", 0
    chembl_id = ""
    try:
        res = _http_json(
            "https://www.ebi.ac.uk/chembl/api/data/target.json"
            f"?target_components__accession={uniprot}&limit=50",
            timeout=timeout,
        )
        targets = res.get("targets") or []
        for target in targets:
            if str(target.get("organism", "")).lower() == "homo sapiens":
                chembl_id = target.get("target_chembl_id") or ""
                break
        if not chembl_id and targets:
            chembl_id = targets[0].get("target_chembl_id") or ""
    except Exception:
        chembl_id = ""
    if not chembl_id:
        return chembl_id, 0
    try:
        res = _http_json(
            "https://www.ebi.ac.uk/chembl/api/data/activity.json"
            f"?target_chembl_id={chembl_id}&limit=1",
            timeout=timeout,
        )
        meta = res.get("page_meta") or {}
        return chembl_id, int(meta.get("total_count") or 0)
    except Exception:
        return chembl_id, 0


def _rcsb_evidence(uniprot: str, timeout: int) -> tuple[int, list[str]]:
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
    try:
        res = _http_json(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            query,
            timeout,
        )
    except Exception:
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
    info = _mygene_info(gene, timeout)
    uniprot = info.get("uniprot") or ""
    ensembl = info.get("ensembl") or ""
    chembl_id, bioactivities = _chembl_evidence(uniprot, timeout)
    pdb_count, pdb_ids = _rcsb_evidence(uniprot, timeout)
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
        except Exception:
            cache = {}

    missing = [gene for gene in genes if gene not in cache]
    if missing and fetch:
        log.info("fetching evidence for %s genes", len(missing))
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [
                pool.submit(_evidence_for_gene, gene, timeout)
                for gene in missing
            ]
            for future in futures:
                try:
                    row = future.result()
                    cache[str(row["gene"])] = row
                except Exception as exc:  # noqa: BLE001
                    log.warning("evidence fetch failed: %s", exc)

    rows = []
    for gene in genes:
        row = cache.get(gene) or _empty_evidence(gene)
        rows.append({key: row.get(key, _empty_evidence(gene)[key]) for key in EVIDENCE_COLUMNS})
    frame = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, index=False)
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

    prognosis = pd.DataFrame({"gene": genes, "hr": 1.0, "p": 1.0})
    prognosis.to_csv(ko_dir / "prognosis.csv", index=False)

    summary = {
        "expression_csv": str(expression_dst),
        "metadata_csv": str(metadata_dst),
        "prognosis_csv": str(ko_dir / "prognosis.csv"),
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
    except Exception:
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
    except (TypeError, ValueError):
        return False


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
    except ImportError:
        return []
    try:
        lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
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
            except ValueError:
                continue
        block_lines = list(atom_lines)
        for line in conect:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                ids = [int(part) for part in parts[1:] if part.isdigit()]
            except ValueError:
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
        except Exception:
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
    if chembl_id and not pd.isna(chembl_id):
        cfg.data.setdefault("evidence", {})["chembl_target_id"] = str(chembl_id)
    try:
        evidence_mod.gather_evidence(cfg, log)
    except Exception as exc:  # noqa: BLE001
        log.warning("evidence collection failed for %s: %s", gene, exc)

    known = target_dir / "evidence" / "known_ligands.csv"
    library = _prepare_ligand_library(
        target_dir,
        known,
        ligand_library,
        pdb_path=pdb_path,
    )
    if library is None:
        base["error"] = "no ligand library available"
        return base

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
                best = str(ranked_df["affinity"].min())
        except Exception:
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
) -> dict:
    out_dir = _integration_dir(workdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = (
        pd.read_csv(key_genes_csv)["gene"]
        .astype(str)
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
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "docking_targets.csv", index=False)
    summary = {
        "status": "completed",
        "targets_requested": len(genes),
        "ok": int((frame["status"] == "ok").sum()),
        "failed": int((frame["status"] == "failed").sum()),
        "skipped": int((frame["status"] == "skipped").sum()),
        "total_hits": int(pd.to_numeric(frame["hits"], errors="coerce").fillna(0).sum()),
        "best_affinity": (
            str(
                frame.loc[
                    pd.to_numeric(frame["hits"], errors="coerce").fillna(0) > 0,
                    "best_affinity",
                ].min()
            )
            if len(frame)
            else ""
        ),
        "output_csv": str(out_dir / "docking_targets.csv"),
    }
    write_json(out_dir / "docking_summary.json", summary)
    log.info(
        "docking stage complete: %s ok / %s failed / %s skipped",
        summary["ok"],
        summary["failed"],
        summary["skipped"],
    )
    return summary


def _esc(value) -> str:
    import html

    return html.escape(str(value if value is not None else ""))


def _render_table(frame: pd.DataFrame, columns: list[str]) -> str:
    if frame is None or frame.empty:
        return '<p class="muted">No data.</p>'
    head = "".join(f"<th>{_esc(c)}</th>" for c in columns)
    body = ""
    for _, row in frame.head(20).iterrows():
        cells = "".join(f"<td>{_esc(row.get(c, ''))}</td>" for c in columns)
        body += f"<tr>{cells}</tr>"
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def generate_integrated_report(
    workdir: Path,
    single_cell_root: Path,
    docking_config: Path,
    ctx: dict,
) -> Path:
    out_dir = _integration_dir(workdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sc_summary = _read_json(single_cell_root / "results" / "summary.json")
    dataset_mode = str(sc_summary.get("dataset_mode", "single_cell"))
    sample_label = "samples" if dataset_mode != "single_cell" else "cells"
    key_genes = pd.read_csv(out_dir / "key_genes.csv") if (out_dir / "key_genes.csv").exists() else pd.DataFrame()
    ko_summary = _read_json(out_dir / "knockout_summary.json")
    ko_top = pd.DataFrame()
    ko_ranked = (
        workdir
        / "outputs"
        / "run_001"
        / "results"
        / "04_knockout"
        / "data"
        / "fig_52_53_ranked_knockout.csv"
    )
    if not ko_ranked.exists():
        ko_ranked = (
            workdir
            / "outputs"
            / "run_001"
            / "results"
            / "04_knockout"
            / "data"
            / "fig_52_53_ranked_knockout.csv"
        )
    if ko_ranked.exists():
        ko_top = pd.read_csv(ko_ranked)
    docking_summary = _read_json(out_dir / "docking_summary.json")
    docking = pd.read_csv(out_dir / "docking_targets.csv") if (out_dir / "docking_targets.csv").exists() else pd.DataFrame()
    evidence = pd.read_csv(out_dir / "gene_evidence.csv") if (out_dir / "gene_evidence.csv").exists() else pd.DataFrame()
    feedback_summary = _read_json(out_dir / "cell_feedback" / "cell_feedback_summary.json")
    feedback_targets = (
        pd.read_csv(out_dir / "cell_feedback" / "data" / "feedback_targets.csv")
        if (out_dir / "cell_feedback" / "data" / "feedback_targets.csv").exists()
        else pd.DataFrame()
    )
    feedback_enrichment = (
        pd.read_csv(out_dir / "cell_feedback" / "data" / "celltype_enrichment.csv")
        if (out_dir / "cell_feedback" / "data" / "celltype_enrichment.csv").exists()
        else pd.DataFrame()
    )
    feedback_deg = (
        pd.read_csv(out_dir / "cell_feedback" / "data" / "feedback_deg.csv")
        if (out_dir / "cell_feedback" / "data" / "feedback_deg.csv").exists()
        else pd.DataFrame()
    )
    feedback_go = (
        pd.read_csv(out_dir / "cell_feedback" / "data" / "feedback_enrichment_go.csv")
        if (out_dir / "cell_feedback" / "data" / "feedback_enrichment_go.csv").exists()
        else pd.DataFrame()
    )
    feedback_kegg = (
        pd.read_csv(
            out_dir / "cell_feedback" / "data" / "feedback_enrichment_kegg.csv"
        )
        if (
            out_dir
            / "cell_feedback"
            / "data"
            / "feedback_enrichment_kegg.csv"
        ).exists()
        else pd.DataFrame()
    )
    qc_metrics = _read_json(out_dir / "qc_metrics.json")
    differential_abundance = (
        pd.read_csv(out_dir / "differential_abundance.csv")
        if (out_dir / "differential_abundance.csv").exists()
        else pd.DataFrame()
    )
    differential_abundance_summary = _read_json(
        out_dir / "differential_abundance_summary.json"
    )
    qc_gate = qc_metrics.get("qc_gate") or {}
    qc_gate_frame = pd.DataFrame(
        qc_gate.get("checks") or [],
        columns=["name", "level", "ok", "message"],
    )

    sc_html = _render_table(
        pd.DataFrame(
            [
                {
                    "accession": sc_summary.get("dataset", ""),
                    sample_label: sc_summary.get("n_cells_after_doublet_removal", ""),
                    "genes": sc_summary.get("n_genes", ""),
                    "deg_up": sc_summary.get("deg_up", ""),
                    "deg_down": sc_summary.get("deg_down", ""),
                }
            ]
        ),
        ["accession", sample_label, "genes", "deg_up", "deg_down"],
    )
    ko_cols = [
        c
        for c in [
            "rank",
            "gene",
            "target_class",
            "target_score",
            "knockout_score",
            "druggability_score",
        ]
        if c in ko_top.columns
    ]
    dock_cols = [
        c
        for c in [
            "gene",
            "status",
            "pdb_id",
            "ligand_count",
            "hits",
            "best_affinity",
        ]
        if c in docking.columns
    ]
    ev_cols = [
        c
        for c in [
            "gene",
            "uniprot",
            "known_ligands",
            "pdb_structures",
            "pdb_ids",
            "string_partners",
            "reactome_pathways",
            "pharmgkb_annotations",
            "alphafold_structures",
            "opentargets_hits",
            "kegg_pathways",
            "database_sources",
        ]
        if c in evidence.columns
    ]
    feedback_cols = [
        c
        for c in [
            "gene",
            "source",
            "feedback_score",
            "target_score",
            "knockout_score",
            "docking_hits",
            "cell_detection_rate",
            "celltype_specificity",
            "cell_support_score",
            "top_celltype",
        ]
        if c in feedback_targets.columns
    ]
    feedback_enrichment_cols = [
        c
        for c in ["celltype", "n_cells", "module_mean", "module_diff", "p_adjust"]
        if c in feedback_enrichment.columns
    ]
    feedback_deg_cols = [
        c
        for c in [
            "gene",
            "avg_log2FC",
            "pct.1",
            "pct.2",
            "p_val_adj",
            "direction",
            "significant",
        ]
        if c in feedback_deg.columns
    ]
    feedback_go_cols = [
        c
        for c in [
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "pvalue",
            "p.adjust",
            "Count",
            "geneID",
        ]
        if c in feedback_go.columns
    ]
    if not feedback_go_cols and "note" in feedback_go.columns:
        feedback_go_cols = ["note"]
    feedback_kegg_cols = [
        c
        for c in [
            "ID",
            "Description",
            "GeneRatio",
            "BgRatio",
            "pvalue",
            "p.adjust",
            "Count",
            "geneID",
        ]
        if c in feedback_kegg.columns
    ]
    if not feedback_kegg_cols and "note" in feedback_kegg.columns:
        feedback_kegg_cols = ["note"]
    differential_abundance_cols = [
        c
        for c in [
            "celltype",
            "n_cells",
            "chi2",
            "p_value",
            "p_adjust",
            "significant",
            "direction",
        ]
        if c in differential_abundance.columns
    ]

    def rel(path):
        try:
            return os.path.relpath(path, out_dir)
        except ValueError:
            return str(path)

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Integrated Discovery Pipeline Report</title>
<style>
body {{ font-family: "Segoe UI", Arial, sans-serif; margin: 24px; color: #1f2933; background: #f5f7fa; }}
h1 {{ font-size: 24px; }}
h2 {{ font-size: 18px; margin-top: 22px; }}
.card {{ background: #fff; border: 1px solid #e4e7eb; border-radius: 8px; padding: 16px; margin-bottom: 14px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
th, td {{ border: 1px solid #e5e7eb; padding: 6px 7px; text-align: left; }}
th {{ background: #eef2f7; }}
.muted {{ color: #6b7280; }}
a {{ color: #1d4ed8; }}
</style>
</head>
<body>
<h1>Integrated Discovery Pipeline Report</h1>
<div class="card">
  <p><b>Analysis output:</b> {rel(single_cell_root)}</p>
  <p><b>Integration output:</b> {rel(out_dir)}</p>
  <p><b>Docking summary:</b> {_esc(docking_summary)}</p>
</div>
<div class="card">
  <h2>Expression analysis summary</h2>
  {sc_html}
</div>
<div class="card">
  <h2>QC gate ({_esc(qc_gate.get("status", "skipped"))})</h2>
  <p class="muted">{_esc(qc_gate.get("summary", ""))}</p>
  {_render_table(qc_gate_frame, ["name", "level", "ok", "message"])}
</div>
<div class="card">
  <h2>Differential abundance (cell type composition)</h2>
  <p class="muted">{_esc(differential_abundance_summary.get("reason", ""))}</p>
  {_render_table(differential_abundance, differential_abundance_cols)}
</div>
<div class="card">
  <h2>Key genes (top 20)</h2>
  {_render_table(key_genes, ["rank", "gene", "direction", "avg_log2fc", "p_val_adj"])}
</div>
<div class="card">
  <h2>Virtual knockout targets (top 20)</h2>
  {_render_table(ko_top, ko_cols)}
</div>
<div class="card">
  <h2>Docking per target</h2>
  {_render_table(docking, dock_cols)}
</div>
<div class="card">
  <h2>Cell feedback targets</h2>
  {_render_table(feedback_targets, feedback_cols)}
</div>
<div class="card">
  <h2>Cell type enrichment</h2>
  {_render_table(feedback_enrichment, feedback_enrichment_cols)}
</div>
<div class="card">
  <h2>Cell feedback differential expression</h2>
  {_render_table(feedback_deg, feedback_deg_cols)}
</div>
<div class="card">
  <h2>Cell feedback GO enrichment (top 5 network)</h2>
  {_render_table(feedback_go, feedback_go_cols)}
</div>
<div class="card">
  <h2>Cell feedback KEGG enrichment (top 5 network)</h2>
  {_render_table(feedback_kegg, feedback_kegg_cols)}
</div>
<div class="card">
  <h2>Gene evidence</h2>
  {_render_table(evidence, ev_cols)}
</div>
<div class="card">
  <h2>Outputs</h2>
  <ul>
    <li><a href="{rel(out_dir / 'key_genes.csv')}">key_genes.csv</a></li>
    <li><a href="{rel(out_dir / 'differential_abundance.csv') if (out_dir / 'differential_abundance.csv').exists() else '#'}">differential_abundance.csv</a></li>
    <li><a href="{rel(out_dir / 'qc_metrics.json')}">qc_metrics.json</a></li>
    <li><a href="{rel(ko_ranked) if ko_ranked.exists() else '#'}">fig_52_53_ranked_knockout.csv</a></li>
    <li><a href="{rel(out_dir / 'docking_targets.csv') if (out_dir / 'docking_targets.csv').exists() else '#'}">docking_targets.csv</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'data' / 'feedback_targets.csv') if (out_dir / 'cell_feedback' / 'data' / 'feedback_targets.csv').exists() else '#'}">cell_feedback_targets.csv</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'data' / 'feedback_deg.csv') if (out_dir / 'cell_feedback' / 'data' / 'feedback_deg.csv').exists() else '#'}">feedback_deg.csv</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'data' / 'feedback_enrichment_go.csv') if (out_dir / 'cell_feedback' / 'data' / 'feedback_enrichment_go.csv').exists() else '#'}">feedback_enrichment_go.csv</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'data' / 'feedback_enrichment_kegg.csv') if (out_dir / 'cell_feedback' / 'data' / 'feedback_enrichment_kegg.csv').exists() else '#'}">feedback_enrichment_kegg.csv</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'figures' / 'fig_54_feedback_module_umap.png') if (out_dir / 'cell_feedback' / 'figures' / 'fig_54_feedback_module_umap.png').exists() else '#'}">fig_54_feedback_module_umap.png</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'figures' / 'fig_59_feedback_targets_volcano.png') if (out_dir / 'cell_feedback' / 'figures' / 'fig_59_feedback_targets_volcano.png').exists() else '#'}">fig_59_feedback_targets_volcano.png</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'figures' / 'fig_60_feedback_condition_violin.png') if (out_dir / 'cell_feedback' / 'figures' / 'fig_60_feedback_condition_violin.png').exists() else '#'}">fig_60_feedback_condition_violin.png</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'figures' / 'fig_61_feedback_go_network.png') if (out_dir / 'cell_feedback' / 'figures' / 'fig_61_feedback_go_network.png').exists() else '#'}">fig_61_feedback_go_network.png</a></li>
    <li><a href="{rel(out_dir / 'cell_feedback' / 'figures' / 'fig_62_feedback_kegg_network.png') if (out_dir / 'cell_feedback' / 'figures' / 'fig_62_feedback_kegg_network.png').exists() else '#'}">fig_62_feedback_kegg_network.png</a></li>
  </ul>
</div>
</body>
</html>
"""
    report_path = out_dir / "integration_report.html"
    report_path.write_text(html_text, encoding="utf-8")

    summary = {
        "single_cell": sc_summary,
        "qc_gate": qc_gate,
        "differential_abundance": differential_abundance_summary,
        "key_genes": len(key_genes),
        "knockout": {
            "genes_scored": (ko_summary.get("knockout") or {}).get("genes_scored", 0),
            "validation_candidates": (ko_summary.get("validation") or {}).get("candidates", 0),
        },
        "docking": docking_summary,
        "cell_feedback": {
            "status": feedback_summary.get("status", "skipped"),
            "genes_matched": feedback_summary.get("genes_matched", 0),
            "deg_genes": len(feedback_deg),
            "go_terms": len(feedback_go) if "ID" in feedback_go.columns else 0,
            "kegg_terms": len(feedback_kegg) if "ID" in feedback_kegg.columns else 0,
            "go_top5": (
                feedback_go["Description"].head(5).tolist()
                if "Description" in feedback_go.columns
                else []
            ),
            "kegg_top5": (
                feedback_kegg["Description"].head(5).tolist()
                if "Description" in feedback_kegg.columns
                else []
            ),
            "n_celltypes": feedback_summary.get("n_celltypes", 0),
            "top_celltypes": feedback_summary.get("top_celltypes", []),
            "figures": feedback_summary.get("figures", []),
        },
        "evidence_genes": len(evidence),
        "evidence_database_sources": (
            ",".join(
                sorted(
                    {
                        source
                        for value in evidence.get("database_sources", [])
                        for source in str(value).split(",")
                        if source
                    }
                )
            )
            if "database_sources" in evidence.columns
            else ""
        ),
        "report_html": str(report_path),
        "finished_at": datetime.now().isoformat(timespec="seconds"),
    }
    write_json(out_dir / "integration_summary.json", summary)
    cfg = load_config(docking_config, {"workdir": str(workdir)})
    write_run_manifest(
        out_dir,
        cfg,
        "full-pipeline",
        {
            "key_genes_csv": out_dir / "key_genes.csv",
            "gene_evidence_csv": out_dir / "gene_evidence.csv",
            "integration_summary_json": out_dir / "integration_summary.json",
        },
        summary,
    )
    log.info("integrated report generated: %s", report_path)
    return report_path


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
            except (OSError, ValueError):
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
    )
    ctx["key_genes_path"] = _integration_dir(workdir) / "key_genes.csv"
    ctx["key_genes"] = frame


def _stage_evidence(args, workdir: Path, ctx: dict) -> None:
    genes = pd.read_csv(ctx["key_genes_path"])["gene"].astype(str).tolist()
    frame = ensure_gene_evidence(
        genes,
        workdir,
        fetch=not args.skip_evidence_fetch,
        max_workers=args.evidence_workers,
        timeout=args.evidence_timeout,
    )
    ctx["evidence"] = frame
    ctx["evidence_path"] = _integration_dir(workdir) / "gene_evidence.csv"


def _stage_knockout_inputs(args, workdir: Path, ctx: dict) -> None:
    ctx["knockout_inputs"] = build_knockout_inputs(
        ctx["single_cell_root"],
        workdir,
        skip_pseudobulk=args.skip_pseudobulk,
    )


def _stage_knockout(args, workdir: Path, ctx: dict) -> None:
    if args.skip_knockout:
        summary = {"status": "skipped", "reason": "knockout disabled by arguments"}
        write_json(_integration_dir(workdir) / "knockout_summary.json", summary)
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
    ctx["docking"] = run_docking_stage(
        workdir,
        ctx["docking_config"],
        ctx["key_genes_path"],
        ctx["evidence_path"],
        max_targets=args.docking_targets,
        ligand_library=args.ligand_library,
        force=args.force,
    )


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
            except (OSError, ValueError):
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
        "docking_config": cfg_path,
    }
    if not args.force:
        _invalidate_markers_for_changed_root(workdir, ctx["single_cell_root"])
    if args.dry_run:
        return _dry_run_stages(args, workdir, ctx)

    stage_fns = {
        "01": _stage_single_cell,
        "02": _stage_key_targets,
        "03": _stage_evidence,
        "04": _stage_knockout_inputs,
        "05": _stage_knockout,
        "06": _stage_docking,
        "07": _stage_cell_feedback,
        "08": _stage_report,
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
        "docking_targets": 3,
        "ml_model": "xgb",
        "keep_all_genes": False,
        "case_label": None,
        "normal_label": None,
        "ligand_library": None,
        "ko_top_n": None,
        "depmap_csv": None,
        "ppi_network_csv": None,
        "cell_feedback": {
            "enabled": True,
            "top_n": 12,
            "max_features": 8,
            "timeout_seconds": 3600,
        },
        "evidence": {"fetch": True, "max_workers": 6, "timeout": 90},
        "qc_gate": DEFAULT_QC_GATE,
        "differential_abundance": DEFAULT_DIFFERENTIAL_ABUNDANCE,
        "gene_blacklist": DEFAULT_GENE_BLACKLIST,
    }
    raw = {}
    if path.exists():
        raw = json.loads(path.read_text(encoding="utf-8"))
    config = dict(defaults)
    config.update(raw or {})
    config["cell_feedback"] = dict(defaults["cell_feedback"])
    config["cell_feedback"].update((raw.get("cell_feedback") or {}))
    config["evidence"] = dict(defaults["evidence"])
    config["evidence"].update((raw.get("evidence") or {}))
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
        "--ml-model",
        choices=["xgb", "rf", "gbm", "mlp", "lasso_svm"],
        default=None,
        help="single-cell ML model used by the integrated pipeline",
    )
    parser.add_argument("--docking-targets", type=int, help="max genes to dock")
    parser.add_argument("--ligand-library", help="ligand library file (.smi/.sdf/.csv)")
    parser.add_argument("--case-label", help="case group label for knockout")
    parser.add_argument("--normal-label", help="normal group label for knockout")
    parser.add_argument("--ko-top-n", type=int, help="top N knockout report genes")
    parser.add_argument("--depmap-csv", help="DepMap gene effect CSV")
    parser.add_argument(
        "--ppi-network-csv",
        help="STRING-style PPI edge table for knockout PPI hub scoring",
    )
    parser.add_argument("--evidence-workers", type=int, default=None)
    parser.add_argument("--evidence-timeout", type=int, default=None)
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
    except KeyboardInterrupt:
        log.info("interrupted")
        return 130


if __name__ == "__main__":
    sys.exit(main())
