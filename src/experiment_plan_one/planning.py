"""Frozen study-plan, cohort-manifest and provenance helpers."""

from __future__ import annotations

import hashlib
import json
import shutil
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from .common import ensure_dir, sha256_file, write_json

PLAN_DIR = "00_plan"
GOVERNANCE_DIR = "governance"
COHORT_COLUMNS = [
    "accession",
    "sample_id",
    "library_id",
    "donor_id",
    "species",
    "platform",
    "tissue",
    "condition",
    "original_diagnosis",
    "paired_patient",
    "batch",
    "exposure_status",
    "included",
    "inclusion_reason",
    "data_version",
    "source",
]

GOVERNANCE_FILES = {
    "script_freeze": "experiment_plan_one_freeze.json",
    "delivery_tiers": "experiment_plan_one_delivery_tiers.json",
    "issue_register": "experiment_plan_one_issue_register.csv",
    "method_changes": "experiment_plan_one_method_changes.csv",
    "verification": "experiment_plan_one_verification.json",
}

ISSUE_REGISTER_COLUMNS = {
    "issue_id",
    "category",
    "scope",
    "problem_statement",
    "evidence",
    "responsibility_role",
    "disposition",
    "status",
    "delivery_tier",
    "in_scope_now",
    "last_reviewed",
}

METHOD_CHANGE_COLUMNS = {
    "change_id",
    "panel_scope",
    "original_method",
    "implemented_method",
    "reason",
    "status",
    "equivalent",
    "requires_protocol_revision",
    "reported_name",
    "responsibility_role",
    "evidence",
}

ISSUE_STATUSES = {
    "resolved_verified",
    "recorded_substitution",
    "limited",
    "blocked",
    "not_run",
}

METHOD_CHANGE_STATUSES = {
    "recorded_substitution",
    "blocked",
    "not_run",
    "equivalent_verified",
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _read_json(path: Path, default: Any) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def governance_sources(root: Path | None = None) -> dict[str, Path]:
    base = (root or _repo_root()) / "config"
    return {
        key: base / filename
        for key, filename in GOVERNANCE_FILES.items()
    }


def analysis_code_files(root: Path | None = None) -> list[Path]:
    base = root or _repo_root()
    files = list((base / "src" / "experiment_plan_one").glob("*.py"))
    files.extend((base / "src" / "experiment_plan_one" / "R").glob("*.R"))
    files.extend(
        [
            base / "src" / "docking" / "insilico.py",
            base / "src" / "docking" / "insilico_enrichment.R",
            base / "src" / "docking" / "md_simulation.py",
        ]
    )
    return sorted(path for path in files if path.is_file())


def analysis_code_fingerprint(root: Path | None = None) -> str:
    base = root or _repo_root()
    hashes = {
        str(path.relative_to(base)): sha256_file(path)
        for path in analysis_code_files(base)
    }
    return hashlib.sha256(
        json.dumps(hashes, sort_keys=True).encode("utf-8")
    ).hexdigest()


def load_governance(root: Path | None = None) -> dict[str, Any]:
    """Load and validate the freeze policy issue register and method log."""
    sources = governance_sources(root)
    missing = [str(path) for path in sources.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"missing governance files: {missing}")
    freeze = _read_json(sources["script_freeze"], {})
    tiers = _read_json(sources["delivery_tiers"], {})
    verification = _read_json(sources["verification"], {})
    issues = pd.read_csv(sources["issue_register"]).fillna("")
    changes = pd.read_csv(sources["method_changes"]).fillna("")
    if not ISSUE_REGISTER_COLUMNS.issubset(issues.columns):
        raise ValueError("issue register is missing required columns")
    if not METHOD_CHANGE_COLUMNS.issubset(changes.columns):
        raise ValueError("method change log is missing required columns")
    invalid_issue_statuses = set(issues["status"].astype(str)) - ISSUE_STATUSES
    if invalid_issue_statuses:
        raise ValueError(
            f"invalid issue statuses: {sorted(invalid_issue_statuses)}"
        )
    invalid_method_statuses = (
        set(changes["status"].astype(str)) - METHOD_CHANGE_STATUSES
    )
    if invalid_method_statuses:
        raise ValueError(
            f"invalid method-change statuses: {sorted(invalid_method_statuses)}"
        )
    return {
        "sources": sources,
        "freeze": freeze,
        "delivery_tiers": tiers,
        "verification": verification,
        "issue_register": issues,
        "method_changes": changes,
    }


def write_governance_artifacts(
    output_root: Path,
    root: Path | None = None,
) -> dict[str, Path]:
    """Copy freeze governance inputs into the run output and summarize them."""
    governance = load_governance(root)
    out_dir = ensure_dir(Path(output_root) / PLAN_DIR / GOVERNANCE_DIR)
    copied: dict[str, Path] = {}
    for key, source in governance["sources"].items():
        destination = out_dir / source.name
        shutil.copy2(source, destination)
        copied[key] = destination
    issues = governance["issue_register"]
    changes = governance["method_changes"]
    current_fingerprint = analysis_code_fingerprint(root)
    expected_fingerprint = str(
        governance["verification"].get("analysis_code_fingerprint") or ""
    )
    verification_status = str(
        governance["verification"].get("status") or ""
    )
    if (
        str(governance["freeze"].get("state") or "") == "frozen"
        and verification_status != "passed"
    ):
        raise RuntimeError(
            "script freeze verification is not passed; update the "
            "verification record before running the frozen pipeline"
        )
    if expected_fingerprint and current_fingerprint != expected_fingerprint:
        raise RuntimeError(
            "analysis code changed after freeze verification; rerun tests and "
            "update the frozen analysis_code_fingerprint"
        )
    write_json(
        out_dir / "issue_register_summary.json",
        {
            "status_counts": (
                issues["status"].value_counts().to_dict()
                if not issues.empty
                else {}
            ),
            "category_counts": (
                issues["category"].value_counts().to_dict()
                if not issues.empty
                else {}
            ),
            "in_scope_now": int(
                issues["in_scope_now"].astype(str).str.lower().eq("yes").sum()
            )
            if not issues.empty
            else 0,
            "method_change_count": int(len(changes)),
            "analysis_code_fingerprint": current_fingerprint,
            "analysis_code_fingerprint_match": (
                True
                if not expected_fingerprint
                else current_fingerprint == expected_fingerprint
            ),
            "non_equivalent_method_changes": int(
                changes["equivalent"]
                .astype(str)
                .str.lower()
                .ne("true")
                .sum()
            )
            if not changes.empty
            else 0,
        },
    )
    copied["issue_register_summary"] = out_dir / "issue_register_summary.json"
    return copied


def _full_rerun_status(output_root: Path) -> dict[str, Any]:
    manifest = _read_json(
        Path(output_root) / "10_reports" / "analysis_manifest.json",
        {},
    )
    stages = manifest.get("stages") or {}
    governance = load_governance()
    required = list(governance["freeze"].get("full_rerun_stages") or [])
    statuses = {
        stage: str((stages.get(stage) or {}).get("status") or "not_run")
        for stage in required
    }
    incomplete = {
        stage: status
        for stage, status in statuses.items()
        if status not in {"completed", "valid_negative", "valid_positive"}
    }
    return {
        "status": "completed" if required and not incomplete else "not_run",
        "stages": statuses,
        "incomplete_stages": incomplete,
    }


def audit_delivery_readiness(
    output_root: Path,
    *,
    coverage_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Evaluate the three delivery tiers without claiming unrun science."""
    governance = load_governance()
    freeze = governance["freeze"]
    verification = governance["verification"]
    issues = governance["issue_register"]
    changes = governance["method_changes"]
    coverage_path = (
        Path(output_root)
        / "10_reports"
        / "plan_coverage"
        / "plan_coverage.csv"
    )
    coverage = (
        pd.read_csv(coverage_path).fillna("")
        if coverage_path.exists()
        else pd.DataFrame()
    )
    panel_statuses = (
        coverage["current_status"].astype(str).value_counts().to_dict()
        if not coverage.empty and "current_status" in coverage
        else {}
    )
    unresolved = {
        status: int(count)
        for status, count in panel_statuses.items()
        if status
        in {
            "missing",
            "missing_audit",
            "blocked",
            "not_run",
            "prepared_not_run",
            "needs_revision",
        }
    }
    full_rerun = _full_rerun_status(output_root)
    non_equivalent_changes = changes[
        changes["equivalent"].astype(str).str.lower().ne("true")
    ]
    blocking_method_changes = changes[
        changes["status"].astype(str).eq("blocked")
        | changes["reported_name"].astype(str).str.strip().eq("")
    ]
    unresolved_issues = issues[
        issues["status"].astype(str).isin({"blocked", "not_run"})
    ]
    tier1_met = (
        str(freeze.get("state") or "") == "frozen"
        and bool(freeze.get("analysis_modules_locked"))
        and str(verification.get("status") or "") == "passed"
        and (
            not verification.get("analysis_code_fingerprint")
            or str(verification.get("analysis_code_fingerprint"))
            == analysis_code_fingerprint()
        )
    )
    tier2_met = (
        tier1_met
        and full_rerun["status"] == "completed"
        and not unresolved
        and blocking_method_changes.empty
        and (
            not coverage_summary
            or not coverage_summary.get("missing_or_unmet_panels")
        )
    )
    experimental_manifest = _read_json(
        Path(output_root) / PLAN_DIR / "experimental_validation_manifest.json",
        {},
    )
    experimental_records = experimental_manifest.get("records") or []
    tier3_met = bool(experimental_records) and all(
        str(record.get("status") or "") == "completed_verified"
        for record in experimental_records
        if isinstance(record, dict)
    )
    tier_status = {
        "T1": {
            "name": "reproducible_computational_pipeline",
            "status": "met" if tier1_met else "not_met",
        },
        "T2": {
            "name": "complete_42_panel_evidence",
            "status": "met" if tier2_met else "not_run",
            "full_rerun": full_rerun,
            "unresolved_panel_statuses": unresolved,
            "non_equivalent_method_changes": (
                non_equivalent_changes[
                    [
                        "change_id",
                        "panel_scope",
                        "original_method",
                        "implemented_method",
                        "status",
                    ]
                ].to_dict("records")
                if not non_equivalent_changes.empty
                else []
            ),
            "blocking_method_changes": (
                blocking_method_changes[
                    [
                        "change_id",
                        "panel_scope",
                        "original_method",
                        "implemented_method",
                        "status",
                    ]
                ].to_dict("records")
                if not blocking_method_changes.empty
                else []
            ),
        },
        "T3": {
            "name": "wet_lab_causal_validation",
            "status": "met" if tier3_met else "blocked",
            "experimental_manifest": str(
                Path(output_root)
                / PLAN_DIR
                / "experimental_validation_manifest.json"
            ),
            "reason": (
                ""
                if tier3_met
                else "real exposure intervention rescue or binding validation has not been supplied"
            ),
        },
    }
    return {
        "freeze_id": freeze.get("freeze_id"),
        "script_freeze_state": freeze.get("state"),
        "publication_grade": bool(tier1_met and tier2_met and tier3_met),
        "publication_grade_note": (
            "Publication-grade status requires a completed full rerun "
            "resolved panels and wet-lab support for causal claims."
        ),
        "tiers": tier_status,
        "unresolved_issue_counts": (
            unresolved_issues["status"].value_counts().to_dict()
            if not unresolved_issues.empty
            else {}
        ),
        "policy": governance["delivery_tiers"].get("publication_grade_policy")
        or {},
    }


def dataset_registry(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Return the validated endpoint/data-role registry used by the pipeline."""
    datasets = config.get("datasets") or {}
    defaults: dict[str, dict[str, Any]] = {
        "GSE89632": {
            "species": "Homo sapiens",
            "platform": "GPL14951",
            "tissue": "liver",
            "role": "development",
            "endpoint": "healthy_control_vs_NAFLD",
            "validation_status": "development_cohort",
            "biological_unit": "patient/sample annotation",
            "metadata_note": "Preserve HC, SS and NASH labels; SS/NASH are combined only in the explicit HC-vs-NAFLD contrast.",
        },
        "GSE135251": {
            "species": "Homo sapiens",
            "platform": "RNA-seq",
            "tissue": "liver",
            "role": "same_endpoint_external_candidate",
            "endpoint": "healthy_control_vs_NAFLD",
            "validation_status": "previously_exposed_development_support",
            "biological_unit": "patient",
            "metadata_note": "Do not describe as completely unseen if repository history used this cohort during model or preprocessing development.",
        },
        "GSE49541": {
            "species": "Homo sapiens",
            "platform": "GPL570",
            "tissue": "liver",
            "role": "separate_endpoint_fibrosis_association",
            "endpoint": "advanced_vs_mild_fibrosis",
            "validation_status": "not_same_endpoint_validation",
            "biological_unit": "patient/sample annotation",
            "metadata_note": "Mild and advanced groups are ordinal categories; do not call the 0-1/3-4 labels a continuous stage model.",
        },
        "GSE164441": {
            "species": "Homo sapiens",
            "platform": "RNA-seq",
            "tissue": "liver",
            "role": "HCC_extension_analysis",
            "endpoint": "tumor_vs_patient_matched_adjacent_non_tumor",
            "validation_status": "different_endpoint",
            "biological_unit": "patient",
            "metadata_note": "Adjacent non-tumor is not healthy liver; patient pairing is mandatory.",
        },
        "GSE270583": {
            "species": "Mus musculus",
            "platform": "10x scRNA-seq",
            "tissue": "liver",
            "role": "descriptive_cell_localization",
            "endpoint": "NCD_vs_HFD_descriptive",
            "validation_status": "insufficient_biological_replication_for_disease_inference",
            "biological_unit": "library/animal when verifiable",
            "metadata_note": "Four libraries are NCD x1, HFD x2 and JQF x1; JQF is not a normal control.",
        },
        "GSE202379": {
            "species": "Homo sapiens",
            "platform": "snRNA-seq",
            "tissue": "liver",
            "role": "human_patient_level_localization_validation",
            "endpoint": "healthy_vs_MASLD_MASH_spectrum",
            "validation_status": "patient_level_validation_pending_metadata_reaudit",
            "biological_unit": "donor",
            "metadata_note": "Multiple libraries may map to one donor; use the corrected 2025-06-06 SOFT metadata and never count GSMs as independent patients.",
        },
    }
    rows: list[dict[str, Any]] = []
    for accession, base in defaults.items():
        configured = dict(datasets.get(accession) or {})
        rows.append(
            {
                "accession": accession,
                **base,
                **configured,
                "configured_role": configured.get("role", base["role"]),
                "condition_column": configured.get(
                    "condition_column", "condition"
                ),
                "exposure_status": configured.get("exposure_status", "unknown"),
            }
        )
    return rows


def write_analysis_plan(
    output_root: Path,
    config: dict[str, Any],
) -> dict[str, Path]:
    """Write the frozen analysis contract and initial cohort manifest."""
    plan_dir = ensure_dir(Path(output_root) / PLAN_DIR)
    registry = dataset_registry(config)
    plan = {
        "schema_version": 1,
        "run_id": str(config.get("run_id") or "assigned_at_runtime"),
        "status": "frozen_for_execution",
        "primary_question": "healthy_control_vs_NAFLD_expression_classification",
        "secondary_questions": [
            "SS_vs_NASH_within_GSE89632",
            "fibrosis_severity_association_in_GSE49541",
            "HCC_tumor_vs_patient_matched_adjacent_in_GSE164441",
            "descriptive_mouse_cell_localization",
            "human_donor_level_cell_localization",
        ],
        "primary_development_cohort": "GSE89632",
        "external_validation_policy": {
            "same_endpoint_candidate": "GSE135251",
            "previously_exposed": True,
            "different_endpoints_are_not_external_validation": [
                "GSE49541",
                "GSE164441",
            ],
            "auc_threshold_is_not_a_gate": True,
            "calibration_p_value_alone_is_not_sufficient": True,
        },
        "statistical_units": {
            "bulk": "patient/sample as annotated",
            "single_cell_expression": "donor x cell_type pseudobulk",
            "cell_composition": "donor/sample first, then group summary",
            "communication": "independent donor condition permutation",
        },
        "reporting_rules": {
            "preserve_negative_results": True,
            "mark_methods_by_actual_implementation": True,
            "do_not_claim_wet_lab_validation_from_computation": True,
            "unknown_exposure_remains_unknown": True,
        },
        "datasets": registry,
    }
    plan_path = plan_dir / "analysis_plan.yaml"
    try:
        import yaml

        plan_path.write_text(
            yaml.safe_dump(plan, sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )
    except Exception:
        plan_path.write_text(
            json.dumps(plan, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )

    cohort = pd.DataFrame(
        [
            {
                "accession": row["accession"],
                "sample_id": "",
                "library_id": "",
                "donor_id": "",
                "species": row["species"],
                "platform": row["platform"],
                "tissue": row["tissue"],
                "condition": row["endpoint"],
                "original_diagnosis": "",
                "paired_patient": "yes" if row["accession"] == "GSE164441" else "unknown",
                "batch": "",
                "exposure_status": "unknown",
                "included": "planned",
                "inclusion_reason": row["role"],
                "data_version": "to_be_frozen_after_download",
                "source": "analysis_plan",
            }
            for row in registry
        ],
        columns=COHORT_COLUMNS,
    )
    cohort_path = plan_dir / "cohort_manifest.tsv"
    cohort.to_csv(cohort_path, sep="\t", index=False)

    inclusion = pd.DataFrame(
        [
            {
                "accession": row["accession"],
                "role": row["role"],
                "endpoint": row["endpoint"],
                "biological_unit": row["biological_unit"],
                "validation_status": row["validation_status"],
                "exposure_status": row.get("exposure_status", "unknown"),
                "decision": "retain_with_documented_scope",
                "reason": row["metadata_note"],
            }
            for row in registry
        ]
    )
    inclusion_path = plan_dir / "sample_inclusion_exclusion.tsv"
    inclusion.to_csv(inclusion_path, sep="\t", index=False)

    corrections = pd.DataFrame(
        [
            {
                "accession": "GSE202379",
                "record_id": "GSM6112262;GSM6112263",
                "official_revision_date": "2025-06-06",
                "correction": "sample metadata exchanged; data files were not exchanged",
                "pipeline_action": "parse current GEO SOFT and map every GSM to donor_id before inference",
                "status": "required_audit_at_runtime",
            },
            {
                "accession": "all",
                "record_id": "",
                "official_revision_date": "",
                "correction": "historical NAFLD/NASH labels retained unless diagnostic evidence supports MASLD/MASH mapping",
                "pipeline_action": "store raw diagnosis and normalized analysis condition separately",
                "status": "policy",
            },
        ]
    )
    corrections_path = plan_dir / "metadata_corrections.tsv"
    corrections.to_csv(corrections_path, sep="\t", index=False)

    checksums_path = plan_dir / "data_checksums.json"
    write_json(
        checksums_path,
        {
            "status": "planned",
            "files": {},
            "note": "Updated by data and processing stages after files are available.",
        },
    )
    return {
        "analysis_plan": plan_path,
        "cohort_manifest": cohort_path,
        "sample_inclusion_exclusion": inclusion_path,
        "metadata_corrections": corrections_path,
        "data_checksums": checksums_path,
    }


def update_data_checksums(output_root: Path, files: Iterable[Path]) -> Path:
    path = ensure_dir(Path(output_root) / PLAN_DIR) / "data_checksums.json"
    existing: dict[str, Any] = {}
    if path.exists():
        try:
            existing = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            existing = {}
    records = dict(existing.get("files") or {})
    for file_path in files:
        resolved = Path(file_path).resolve()
        if not resolved.is_file():
            continue
        records[str(resolved)] = {
            "sha256": sha256_file(resolved),
            "size_bytes": int(resolved.stat().st_size),
        }
    write_json(
        path,
        {
            "status": "updated",
            "files": records,
            "note": "Content hashes are frozen inputs, not biological validation.",
        },
    )
    return path


def append_cohort_rows(output_root: Path, rows: Iterable[dict[str, Any]]) -> Path:
    """Merge sample/library/donor rows into the cohort manifest."""
    path = ensure_dir(Path(output_root) / PLAN_DIR) / "cohort_manifest.tsv"
    existing = (
        pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if path.exists()
        else pd.DataFrame(columns=COHORT_COLUMNS)
    )
    incoming = pd.DataFrame(list(rows)).fillna("")
    for column in COHORT_COLUMNS:
        if column not in incoming:
            incoming[column] = ""
    incoming = incoming[COHORT_COLUMNS]
    if existing.empty:
        combined = incoming
    else:
        for column in COHORT_COLUMNS:
            if column not in existing:
                existing[column] = ""
        existing = existing[COHORT_COLUMNS]
        keys = ["accession", "sample_id", "library_id"]
        incoming_keys = set(map(tuple, incoming[keys].astype(str).to_numpy()))
        keep = [
            tuple(row)
            not in incoming_keys
            for row in existing[keys].astype(str).to_numpy()
        ]
        combined = pd.concat([existing.loc[keep], incoming], ignore_index=True)
    combined = combined.drop_duplicates(
        subset=["accession", "sample_id", "library_id"],
        keep="last",
    )
    combined.to_csv(path, sep="\t", index=False)
    return path
