"""Integrated HTML report for the full pipeline."""

from __future__ import annotations

import logging
import os
from datetime import datetime
from pathlib import Path

import pandas as pd

from common.html_utils import esc as _esc  # noqa: E402
from docking.config import load_config  # noqa: E402
from docking.provenance import write_run_manifest  # noqa: E402
from docking.utils import write_json  # noqa: E402

from .stage_paths import _integration_dir, _read_json  # noqa: E402
from .reproducibility import write_reproducibility_manifest  # noqa: E402
from .target_validation import build_target_validation  # noqa: E402
from .external_validation import run_external_validation  # noqa: E402
from .decision_report import build_decision_report  # noqa: E402

log = logging.getLogger("full_pipeline")


def _render_table(frame: pd.DataFrame, columns: list[str]) -> str:
    if frame is None or frame.empty:
        return '<p class="muted">No data.</p>'
    head = "".join(f"<th>{_esc(c)}</th>" for c in columns)
    body = ""
    for _, row in frame.head(20).iterrows():
        cells = "".join(f"<td>{_esc(row.get(c, ''))}</td>" for c in columns)
        body += f"<tr>{cells}</tr>"
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _publication_readiness(
    target_priority_summary: dict,
    evidence_hub_summary: dict,
    cadd_summary: dict,
    docking_config: Path,
    out_dir: Path,
    structural_quality: dict | None = None,
    external_validation_summary: dict | None = None,
    omics_qc_summary: dict | None = None,
) -> dict:
    """Summarize which publication-strength gates have actually been met."""
    try:
        cfg = load_config(docking_config)
        positive_control = bool(
            cfg.get("docking", "positive_control_pdbqt", None)
        )
        replicates = int(cfg.get("docking", "replicates", 1) or 1)
    except Exception:
        positive_control = False
        replicates = 1
    structural_gates = (structural_quality or {}).get("gates") or {}
    positive_control = bool(
        structural_gates.get("positive_control", positive_control)
    )
    replicate_consensus = bool(
        structural_gates.get(
            "replicate_consensus",
            replicates >= 3,
        )
    )
    md_completed = bool(
        structural_gates.get(
            "md_completed",
            int(cadd_summary.get("md_completed", 0) or 0) > 0,
        )
    )
    workdir = out_dir.parents[1]
    benchmark = (
        target_priority_summary.get("benchmark")
        or evidence_hub_summary.get("benchmark")
        or {}
    )
    benchmark_auroc = pd.to_numeric(
        pd.Series([benchmark.get("auroc")]),
        errors="coerce",
    ).iloc[0]
    benchmark_recall = pd.to_numeric(
        pd.Series([benchmark.get("recall_at_n")]),
        errors="coerce",
    ).iloc[0]
    benchmark_independent = bool(
        benchmark.get("benchmark_independent_verified", False)
    )
    benchmark_source = str(benchmark.get("benchmark_source") or "").strip()
    benchmark_version = str(
        benchmark.get("benchmark_version") or ""
    ).strip()
    benchmark_provenance = bool(
        (benchmark_source and benchmark_version)
        or benchmark.get("benchmark_evidence_mode") == "source_holdout"
    )
    benchmark_passed = bool(
        benchmark_independent
        and benchmark_provenance
        and (
            (pd.notna(benchmark_auroc) and float(benchmark_auroc) >= 0.70)
            or (
                pd.notna(benchmark_recall)
                and float(benchmark_recall) >= 0.50
            )
        )
    )
    external_validation_paths = (
        out_dir / "external_validation" / "validation_summary.json",
        out_dir
        / "external_validation"
        / "ml_external_validation.csv",
        workdir
        / "outputs"
        / "run_001"
        / "results"
        / "05_validation"
        / "data"
        / "ml_external_validation.csv",
        workdir
        / "outputs"
        / "run_001"
        / "results"
        / "05_validation"
        / "data"
        / "validation_external.csv",
    )
    mechanistic_path = (
        workdir
        / "outputs"
        / "run_001"
        / "results"
        / "04_knockout"
        / "in_silico"
        / "insilico_summary.json"
    )
    external_validation_completed = (
        str((external_validation_summary or {}).get("status", ""))
        == "completed"
    ) or any(path.exists() for path in external_validation_paths)
    external_auroc = pd.to_numeric(
        pd.Series(
            [(external_validation_summary or {}).get("auroc")]
        ),
        errors="coerce",
    ).iloc[0]
    checks = {
        "omics_qc": bool(
            str((omics_qc_summary or {}).get("status", "")) == "completed"
            and (omics_qc_summary or {}).get("gate_passed", False)
        ),
        "multi_source_evidence": (
            str(evidence_hub_summary.get("status", "")) == "completed"
            and int(
                evidence_hub_summary.get(
                    "successful_source_count",
                    evidence_hub_summary.get("source_count", 0),
                )
                or 0
            )
            >= 3
        ),
        "candidate_universe_not_top50_only": int(
            target_priority_summary.get("targets", 0) or 0
        ) > 50,
        "evidence_target_coverage": (
            int(target_priority_summary.get("targets", 0) or 0) > 0
            and int(
                target_priority_summary.get("evidence_targets_queried", 0)
                or 0
            )
            / max(1, int(target_priority_summary.get("targets", 0) or 0))
            >= 0.50
        ),
        "docking_positive_control": positive_control,
        "docking_replicates": replicate_consensus,
        "md_completed": md_completed,
        "md_rmsd_stability": bool(
            structural_gates.get("md_rmsd_stability", False)
        ),
        "benchmark": benchmark_passed,
        "benchmark_independent": benchmark_independent,
        "benchmark_provenance": benchmark_provenance,
        "external_validation": bool(
            external_validation_completed
            and (external_validation_summary or {}).get(
                "score_provenance_valid",
                False,
            )
            and float(
                (external_validation_summary or {}).get(
                    "target_match_rate",
                    1.0,
                )
                or 0.0
            )
            >= 0.80
            and pd.notna(external_auroc)
            and float(external_auroc) >= 0.70
        ),
        "mechanistic_perturbation": mechanistic_path.exists(),
        "mmpbsa_available": bool(
            structural_gates.get("mmpbsa_available", False)
        ),
        "posebusters_valid": bool(
            structural_gates.get("posebusters_valid", False)
        ),
    }
    required = (
        "omics_qc",
        "multi_source_evidence",
        "candidate_universe_not_top50_only",
        "evidence_target_coverage",
        "docking_positive_control",
        "docking_replicates",
        "md_completed",
        "md_rmsd_stability",
        "benchmark",
        "benchmark_independent",
        "benchmark_provenance",
        "external_validation",
    )
    passed = sum(bool(checks[name]) for name in required)
    if (
        passed == len(required)
        and checks["mechanistic_perturbation"]
        and checks["mmpbsa_available"]
        and checks["posebusters_valid"]
    ):
        level = "publication_grade"
    elif passed == len(required):
        level = "paper_supporting"
    else:
        level = "exploratory"
    missing = [name for name in required if not checks[name]]
    if not checks["mechanistic_perturbation"]:
        missing.append("mechanistic_perturbation")
    return {
        "level": level,
        "passed": passed,
        "required": len(required),
        "checks": checks,
        "missing": missing,
    }


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
    expanded_candidate_path = (
        out_dir / "candidate_universe_evidence_expanded.csv"
    )
    candidate_universe = (
        pd.read_csv(expanded_candidate_path)
        if expanded_candidate_path.exists()
        else (
            pd.read_csv(out_dir / "candidate_universe.csv")
            if (out_dir / "candidate_universe.csv").exists()
            else pd.DataFrame()
        )
    )
    target_priority = (
        pd.read_csv(out_dir / "integrated_target_priority.csv")
        if (out_dir / "integrated_target_priority.csv").exists()
        else (
            pd.read_csv(out_dir / "target_priority.csv")
            if (out_dir / "target_priority.csv").exists()
            else pd.DataFrame()
        )
    )
    target_priority_summary = _read_json(
        out_dir / "integrated_target_priority_summary.json"
    ) or _read_json(out_dir / "target_priority_summary.json")
    target_validation, target_validation_summary = build_target_validation(
        workdir,
        out_dir,
    )
    target_decisions, target_decision_summary = build_decision_report(out_dir)
    external_validation = run_external_validation(
        Path(str(ctx.get("external_validation_path")))
        if ctx.get("external_validation_path")
        else None,
        out_dir,
        pipeline_scores_path=(
            out_dir / "integrated_target_priority.csv"
            if (out_dir / "integrated_target_priority.csv").exists()
            else out_dir / "target_priority.csv"
        ),
        allow_external_score=bool(
            ctx.get("external_validation_allow_external_score", False)
        ),
        target_column=str(
            ctx.get("external_validation_target_column") or "gene"
        ),
        score_column=str(
            ctx.get("external_validation_score_column") or "score"
        ),
        label_column=str(
            ctx.get("external_validation_label_column") or "label"
        ),
        threshold=float(
            ctx.get("external_validation_threshold") or 0.5
        ),
        bootstrap=int(
            ctx.get("external_validation_bootstrap") or 1000
        ),
    )
    evidence_hub_summary = _read_json(
        out_dir / "evidence_hub" / "evidence_hub_summary.json"
    )
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
            / "fig_52_target_candidates.csv"
        )
    if ko_ranked.exists():
        ko_top = pd.read_csv(ko_ranked)
    docking_summary = _read_json(out_dir / "docking_summary.json")
    docking = pd.read_csv(out_dir / "docking_targets.csv") if (out_dir / "docking_targets.csv").exists() else pd.DataFrame()
    evidence = pd.read_csv(out_dir / "gene_evidence.csv") if (out_dir / "gene_evidence.csv").exists() else pd.DataFrame()
    evidence_summary = _read_json(out_dir / "evidence_summary.json")
    cadd_summary = _read_json(out_dir / "cadd_downstream_summary.json")
    cadd_targets = (
        pd.read_csv(out_dir / "cadd_targets.csv")
        if (out_dir / "cadd_targets.csv").exists()
        else pd.DataFrame()
    )
    structural_quality = _read_json(
        out_dir / "structural_quality_summary.json"
    )
    network_summary = _read_json(out_dir / "network_summary.json")
    network_overlap = pd.DataFrame()
    network_overlap_csv = (
        (network_summary.get("outputs") or {}).get("overlap_csv")
        or network_summary.get("overlap_csv")
        or ""
    )
    if network_overlap_csv and Path(str(network_overlap_csv)).exists():
        try:
            network_overlap = pd.read_csv(network_overlap_csv)
        except Exception:
            network_overlap = pd.DataFrame()
    faers_summary = _read_json(out_dir / "faers_summary.json")
    faers_signals = pd.DataFrame()
    faers_csv = faers_summary.get("output_csv") or ""
    if faers_csv and Path(str(faers_csv)).exists():
        try:
            faers_signals = pd.read_csv(faers_csv)
        except Exception:
            faers_signals = pd.DataFrame()
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
    omics_qc_summary = _read_json(out_dir / "omics_qc_summary.json")
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
    readiness = _publication_readiness(
        target_priority_summary,
        evidence_hub_summary,
        cadd_summary,
        docking_config,
        out_dir,
        structural_quality,
        external_validation,
        omics_qc_summary,
    )
    readiness_frame = pd.DataFrame(
        [
            {
                "check": name,
                "passed": bool(value),
            }
            for name, value in (readiness.get("checks") or {}).items()
        ]
    )
    benchmark = (
        target_priority_summary.get("benchmark")
        or evidence_hub_summary.get("benchmark")
        or {}
    )
    benchmark_frame = pd.DataFrame([benchmark]) if benchmark else pd.DataFrame()
    benchmark_cols = [
        column
        for column in (
            "recall_at_n",
            "precision_at_n",
            "auroc",
            "auroc_ci_low",
            "auroc_ci_high",
            "auprc",
            "auprc_ci_low",
            "auprc_ci_high",
            "enrichment_factor",
            "permutation_p_value",
            "n_positive",
            "n_negative",
            "benchmark_source",
            "benchmark_version",
            "benchmark_independent",
            "benchmark_independent_verified",
            "benchmark_evidence_mode",
            "source_leakage",
            "source_overlap",
        )
        if column in benchmark_frame.columns
    ]
    external_validation_frame = (
        pd.DataFrame([external_validation])
        if external_validation
        else pd.DataFrame()
    )
    external_validation_cols = [
        column
        for column in (
            "status",
            "n_samples",
            "n_positive",
            "n_negative",
            "threshold",
            "auroc",
            "auroc_ci_low",
            "auroc_ci_high",
            "auprc",
            "auprc_ci_low",
            "auprc_ci_high",
            "precision",
            "recall",
            "f1",
            "specificity",
            "score_origin",
            "score_provenance_valid",
            "target_match_rate",
        )
        if column in external_validation_frame.columns
    ]

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
    priority_cols = [
        c
        for c in [
            "priority_rank",
            "gene",
            "candidate_origin",
            "decision",
            "integrated_score",
            "evidence_score",
            "genetic_causal",
            "clinical_precedent",
            "safety_risk",
            "coverage_ratio",
            "knockout_score",
            "evidence_status",
            "missing_categories",
        ]
        if c in target_priority.columns
    ]
    validation_cols = [
        c
        for c in [
            "validation_rank",
            "gene",
            "decision",
            "adjusted_score",
            "evidence_completeness",
            "disease_association",
            "druggability",
            "chemical_matter",
            "clinical_precedent",
            "structural_data",
            "clinical_evidence_quality",
            "clinical_evidence_cap_applied",
            "safety_penalty",
        ]
        if c in target_validation.columns
    ]
    decision_cols = [
        column
        for column in (
            "gene",
            "base_decision",
            "target_decision",
            "final_decision",
            "composite_decision",
            "target_action",
            "platform_action",
            "composite_action",
            "action",
            "adjusted_score",
            "safety_risk",
            "target_flags",
            "platform_conflicts",
            "conflict_flags",
        )
        if column in target_decisions.columns
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
    cadd_cols = [
        c
        for c in [
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
        if c in cadd_targets.columns
    ]
    network_cols = [
        c
        for c in [
            "gene",
            "n_sources",
            "sources",
            "ppi_degree",
            "ppi_hub_score",
        ]
        if c in network_overlap.columns
    ]
    faers_cols = [
        c
        for c in [
            "drug",
            "event",
            "a",
            "ror",
            "ror_lower",
            "prr",
            "ic",
            "ebgm",
            "signal",
        ]
        if c in faers_signals.columns
    ]

    def rel(path):
        try:
            return os.path.relpath(path, out_dir)
        except ValueError:
            return str(path)

    def display_paths(value, base=out_dir):
        if isinstance(value, dict):
            return {
                key: display_paths(item, base=base)
                for key, item in value.items()
            }
        if (
            isinstance(value, str)
            and value
            and "\\" in value
            and Path(value).is_absolute()
        ):
            try:
                return os.path.relpath(value, base).replace("\\", "/")
            except ValueError:
                return value
        return value

    docking_display = display_paths(docking_summary)

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
  <p><b>Docking summary:</b> {_esc(docking_display)}</p>
</div>
<div class="card">
  <h2>Publication readiness ({_esc(readiness.get("level", "exploratory"))})</h2>
  <p class="muted">
    Passed {_esc(readiness.get("passed", 0))} of
    {_esc(readiness.get("required", 0))} required gates.
    Missing: {_esc(", ".join(readiness.get("missing") or []) or "none")}.
  </p>
  {_render_table(readiness_frame, ["check", "passed"])}
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
  <h2>Pseudobulk omics QC ({_esc(omics_qc_summary.get("status", "skipped"))})</h2>
  <p class="muted">
    genes: {_esc(omics_qc_summary.get("n_genes", "NA"))};
    samples: {_esc(omics_qc_summary.get("n_samples", "NA"))};
    missing fraction: {_esc(omics_qc_summary.get("missing_fraction", "NA"))};
    groups: {_esc(omics_qc_summary.get("group_counts", {}))}
  </p>
</div>
<div class="card">
  <h2>Final target decisions</h2>
  <p class="muted">
    Actions: {_esc(target_decision_summary.get("actions", {}))}
  </p>
  {_render_table(target_decisions, decision_cols)}
</div>
<div class="card">
  <h2>External validation</h2>
  <p class="muted">{_esc(external_validation.get("reason", ""))}</p>
  {_render_table(external_validation_frame, external_validation_cols)}
</div>
<div class="card">
  <h2>Five-axis target validation (top 20)</h2>
  <p class="muted">
    GO: {_esc((target_validation_summary.get("decisions") or {}).get("GO", 0))};
    median adjusted score:
    {_esc(target_validation_summary.get("median_adjusted_score", "NA"))}
  </p>
  {_render_table(target_validation, validation_cols)}
</div>
<div class="card">
  <h2>Evidence-aware target priority (top 20)</h2>
  <p class="muted">
    Candidate universe: {_esc(len(candidate_universe))};
    evidence hub: {_esc(evidence_hub_summary.get("status", "skipped"))};
    GO: {_esc((target_priority_summary.get("decisions") or {}).get("GO", 0))};
    median coverage:
    {_esc(target_priority_summary.get("median_coverage_ratio", "NA"))}
  </p>
  {_render_table(target_priority, priority_cols)}
  {_render_table(benchmark_frame, benchmark_cols)}
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
  <h2>CADD downstream ({_esc(cadd_summary.get("status", "skipped"))})</h2>
  <p class="muted">{_esc(cadd_summary)}</p>
  {_render_table(cadd_targets, cadd_cols)}
</div>
<div class="card">
  <h2>Network toxicology ({_esc(network_summary.get("status", "skipped"))})</h2>
  <p class="muted">{_esc(network_summary.get("reason", ""))}</p>
  {_render_table(network_overlap, network_cols)}
</div>
<div class="card">
  <h2>FAERS signals ({_esc(faers_summary.get("status", "skipped"))})</h2>
  <p class="muted">{_esc(faers_summary.get("reason", ""))}</p>
  {_render_table(faers_signals, faers_cols)}
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
    <li><a href="{rel(out_dir / 'candidate_universe.csv') if (out_dir / 'candidate_universe.csv').exists() else '#'}">candidate_universe.csv</a></li>
    <li><a href="{rel(out_dir / 'candidate_universe_evidence_expanded.csv') if (out_dir / 'candidate_universe_evidence_expanded.csv').exists() else '#'}">candidate_universe_evidence_expanded.csv</a></li>
    <li><a href="{rel(out_dir / 'target_priority.csv') if (out_dir / 'target_priority.csv').exists() else '#'}">target_priority.csv</a></li>
    <li><a href="{rel(out_dir / 'integrated_target_priority.csv') if (out_dir / 'integrated_target_priority.csv').exists() else '#'}">integrated_target_priority.csv</a></li>
    <li><a href="{rel(out_dir / 'target_validation_scores.csv') if (out_dir / 'target_validation_scores.csv').exists() else '#'}">target_validation_scores.csv</a></li>
    <li><a href="{rel(out_dir / 'external_validation_summary.json') if (out_dir / 'external_validation_summary.json').exists() else '#'}">external_validation_summary.json</a></li>
    <li><a href="{rel(out_dir / 'target_decision_report.csv') if (out_dir / 'target_decision_report.csv').exists() else '#'}">target_decision_report.csv</a></li>
    <li><a href="{rel(out_dir / 'evidence_hub' / 'evidence_coverage.csv') if (out_dir / 'evidence_hub' / 'evidence_coverage.csv').exists() else '#'}">evidence_coverage.csv</a></li>
    <li><a href="{rel(out_dir / 'evidence_hub' / 'source_ablation.csv') if (out_dir / 'evidence_hub' / 'source_ablation.csv').exists() else '#'}">source_ablation.csv</a></li>
    <li><a href="{rel(out_dir / 'differential_abundance.csv') if (out_dir / 'differential_abundance.csv').exists() else '#'}">differential_abundance.csv</a></li>
    <li><a href="{rel(out_dir / 'qc_metrics.json')}">qc_metrics.json</a></li>
    <li><a href="{rel(ko_ranked) if ko_ranked.exists() else '#'}">fig_52_53_ranked_knockout.csv</a></li>
    <li><a href="{rel(out_dir / 'docking_targets.csv') if (out_dir / 'docking_targets.csv').exists() else '#'}">docking_targets.csv</a></li>
    <li><a href="{rel(out_dir / 'cadd_targets.csv') if (out_dir / 'cadd_targets.csv').exists() else '#'}">cadd_targets.csv</a></li>
    <li><a href="{rel(out_dir / 'network_summary.json') if (out_dir / 'network_summary.json').exists() else '#'}">network_summary.json</a></li>
    <li><a href="{rel(out_dir / 'faers_summary.json') if (out_dir / 'faers_summary.json').exists() else '#'}">faers_summary.json</a></li>
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
        "omics_qc": omics_qc_summary,
        "differential_abundance": differential_abundance_summary,
        "key_genes": len(key_genes),
        "candidate_universe": len(candidate_universe),
        "candidate_origins": target_priority_summary.get(
            "candidate_origins",
            {},
        ),
        "target_priority": {
            **target_priority_summary,
            "top_targets": target_priority.head(20).to_dict(orient="records"),
        },
        "target_validation": target_validation_summary,
        "external_validation": external_validation,
        "target_decisions": target_decision_summary,
        "publication_readiness": readiness,
        "knockout": {
            "genes_scored": (ko_summary.get("knockout") or {}).get("genes_scored", 0),
            "validation_candidates": (ko_summary.get("validation") or {}).get("candidates", 0),
        },
        "docking": display_paths(docking_summary),
        "cadd_downstream": display_paths(cadd_summary),
        "network": display_paths(network_summary),
        "faers": display_paths(faers_summary),
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
        "evidence_hub": evidence_hub_summary,
        "evidence_failures": evidence_summary.get("evidence_failures", []),
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
        "report_html": os.path.relpath(
            report_path, out_dir
        ).replace("\\", "/"),
        "finished_at": datetime.now().isoformat(timespec="seconds"),
    }
    write_json(out_dir / "integration_summary.json", summary)
    config_path = Path(
        ctx.get("full_config")
        or (Path(__file__).resolve().parents[2] / "config" / "full_pipeline_config.json")
    )
    reproducibility_path = write_reproducibility_manifest(
        workdir,
        single_cell_root,
        config_path,
        docking_config,
        {
            "targets": len(target_priority),
            "target_validation": target_validation_summary.get(
                "decisions",
                {},
            ),
            "external_validation": external_validation,
            "candidate_origins": target_priority_summary.get(
                "candidate_origins",
                {},
            ),
            "benchmark": benchmark,
            "readiness": readiness,
            "omics_qc_gate_passed": omics_qc_summary.get(
                "gate_passed",
                False,
            ),
        },
    )
    summary["reproducibility_manifest"] = os.path.relpath(
        reproducibility_path,
        out_dir,
    ).replace("\\", "/")
    write_json(out_dir / "integration_summary.json", summary)
    cfg = load_config(docking_config, {"workdir": str(workdir)})
    write_run_manifest(
        out_dir,
        cfg,
        "full-pipeline",
        {
            "key_genes_csv": out_dir / "key_genes.csv",
            "candidate_universe_csv": out_dir / "candidate_universe.csv",
            "candidate_universe_evidence_expanded_csv": (
                out_dir / "candidate_universe_evidence_expanded.csv"
            ),
            "target_priority_csv": out_dir / "target_priority.csv",
            "integrated_target_priority_csv": (
                out_dir / "integrated_target_priority.csv"
            ),
            "target_validation_scores_csv": (
                out_dir / "target_validation_scores.csv"
            ),
            "target_validation_summary_json": (
                out_dir / "target_validation_summary.json"
            ),
            "external_validation_summary_json": (
                out_dir / "external_validation_summary.json"
            ),
            "target_decision_report_csv": (
                out_dir / "target_decision_report.csv"
            ),
            "target_decision_report_json": (
                out_dir / "target_decision_report.json"
            ),
            "gene_evidence_csv": out_dir / "gene_evidence.csv",
            "omics_qc_summary_json": out_dir / "omics_qc_summary.json",
            "evidence_hub_summary_json": (
                out_dir / "evidence_hub" / "evidence_hub_summary.json"
            ),
            "integration_summary_json": out_dir / "integration_summary.json",
            "reproducibility_manifest_json": reproducibility_path,
        },
        summary,
    )
    log.info("integrated report generated: %s", report_path)
    return report_path
