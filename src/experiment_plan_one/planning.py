"""Frozen study-plan, cohort-manifest and provenance helpers."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from .common import ensure_dir, sha256_file, write_json

PLAN_DIR = "00_plan"
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
