"""Frozen study-plan, cohort-manifest and provenance helpers."""

from __future__ import annotations

import hashlib
import json
import re
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
SOURCE_AUTHORIZATIONS_FILE = (
    "experiment_plan_one_source_authorizations.json"
)
REQUIRED_FINAL_ACCESSIONS = {
    "GSE89632",
    "GSE49541",
    "GSE164441",
    "GSE135251",
    "GSE270583",
    "GSE202379",
}
PATIENT_LEVEL_ACCESSIONS = {"GSE164441", "GSE202379"}
ANIMAL_LEVEL_ACCESSIONS = {"GSE270583"}

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
    "retired_original_restored",
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
            base / "src" / "docking" / "R" / "sctenifoldknk.R",
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
    cohort_freeze = _read_json(
        Path(output_root) / PLAN_DIR / "cohort_freeze.json",
        {},
    )
    retired_changes = changes[
        changes["status"].astype(str).eq("retired_original_restored")
    ]
    active_changes = changes[
        ~changes["status"].astype(str).eq("retired_original_restored")
    ]
    non_equivalent_changes = active_changes[
        active_changes["equivalent"].astype(str).str.lower().ne("true")
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
        and str(cohort_freeze.get("status") or "") == "frozen_complete"
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
            "cohort_freeze_status": cohort_freeze.get("status"),
            "cohort_status": cohort_freeze.get("cohort_status"),
            "donor_mapping_status": cohort_freeze.get(
                "donor_mapping_status"
            ),
            "authorization_status": cohort_freeze.get(
                "authorization_status"
            ),
            "cohort_validation_errors": cohort_freeze.get(
                "cohort_validation_errors"
            )
            or [],
            "source_authorization_errors": cohort_freeze.get(
                "source_authorization_errors"
            )
            or [],
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
            "retired_original_method_changes": (
                retired_changes[
                    [
                        "change_id",
                        "panel_scope",
                        "original_method",
                        "implemented_method",
                        "status",
                    ]
                ].to_dict("records")
                if not retired_changes.empty
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


def _source_authorizations_config_path(
    root: Path | None = None,
) -> Path:
    return (root or _repo_root()) / "config" / SOURCE_AUTHORIZATIONS_FILE


def _resolve_authorization_file(
    value: Any,
    *,
    root: Path,
) -> Path | None:
    text = str(value or "").strip()
    if not text:
        return None
    path = Path(text).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def freeze_source_authorizations(
    output_root: Path,
    *,
    root: Path | None = None,
) -> tuple[Path, dict[str, Any]]:
    """Freeze licensed/open data provenance and optional authorization files."""
    repo_root = root or _repo_root()
    config_path = _source_authorizations_config_path(repo_root)
    errors = (
        []
        if config_path.is_file()
        else [f"source authorization config is missing: {config_path}"]
    )
    payload: dict[str, Any] = _read_json(config_path, {})
    sources = payload.get("sources") if isinstance(payload, dict) else []
    frozen_sources: list[dict[str, Any]] = []
    for raw in sources or []:
        if not isinstance(raw, dict):
            errors.append("source authorization entry is not an object")
            continue
        source = dict(raw)
        source_id = str(source.get("source_id") or "").strip()
        access_type = str(source.get("access_type") or "licensed").lower()
        local_file = _resolve_authorization_file(
            source.get("local_file"),
            root=repo_root,
        )
        authorization_file = _resolve_authorization_file(
            source.get("authorization_file"),
            root=repo_root,
        )
        authorization_reference = str(
            source.get("authorization_reference") or ""
        ).strip()
        status = str(source.get("status") or "").strip()
        if access_type == "open":
            status = "available_open_access"
        elif local_file is None:
            status = status or "not_provided"
            if status != "not_required":
                errors.append(
                    f"{source_id}: licensed local export is not provided"
                )
        elif not local_file.is_file():
            status = "blocked_missing_file"
            errors.append(
                f"{source_id}: local export does not exist: {local_file}"
            )
        elif not authorization_reference:
            status = "blocked_missing_authorization_reference"
            errors.append(
                f"{source_id}: authorization reference is required"
            )
        else:
            status = "available_with_authorization"
        if (
            authorization_file is not None
            and not authorization_file.is_file()
        ):
            status = "blocked_missing_authorization_file"
            errors.append(
                f"{source_id}: authorization file does not exist: "
                f"{authorization_file}"
            )
        record = {
            **source,
            "source_id": source_id,
            "status": status,
            "local_file": str(local_file) if local_file else "",
            "local_file_sha256": (
                sha256_file(local_file)
                if local_file and local_file.is_file()
                else ""
            ),
            "authorization_file": (
                str(authorization_file) if authorization_file else ""
            ),
            "authorization_file_sha256": (
                sha256_file(authorization_file)
                if authorization_file and authorization_file.is_file()
                else ""
            ),
        }
        frozen_sources.append(record)
    status_counts: dict[str, int] = {}
    for record in frozen_sources:
        status = str(record.get("status") or "unknown")
        status_counts[status] = status_counts.get(status, 0) + 1
    frozen = {
        "schema_version": int(payload.get("schema_version", 1) or 1),
        "status": "complete" if not errors else "blocked_incomplete",
        "config_file": str(config_path),
        "config_sha256": (
            sha256_file(config_path) if config_path.is_file() else ""
        ),
        "status_counts": status_counts,
        "errors": errors,
        "sources": frozen_sources,
    }
    output_path = (
        ensure_dir(Path(output_root) / PLAN_DIR)
        / "source_authorizations.frozen.json"
    )
    write_json(output_path, frozen)
    return output_path, frozen


def build_donor_mapping(
    cohort: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    """Build and validate the library-to-donor mapping used for inference."""
    frame = cohort.copy().fillna("")
    required = {"accession", "sample_id", "library_id", "donor_id", "condition"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(
            f"cohort manifest is missing donor-mapping columns: {sorted(missing)}"
        )
    rows: list[dict[str, Any]] = []
    errors: list[str] = []
    for row in frame.to_dict("records"):
        accession = str(row.get("accession") or "")
        donor_id = str(row.get("donor_id") or "").strip()
        sample_id = str(row.get("sample_id") or "")
        library_id = str(row.get("library_id") or sample_id)
        if accession == "GSE164441" and donor_id:
            donor_id = re.sub(r"[NT]$", "", donor_id, flags=re.IGNORECASE)
        if accession in ANIMAL_LEVEL_ACCESSIONS:
            mapping_status = (
                "verified_library_to_animal"
                if donor_id
                else "unresolved_library_level_only"
            )
            evidence = (
                "current GEO SOFT/RAW animal identifier"
                if donor_id
                else "four GEO libraries only; independent animal IDs not established"
            )
        elif accession in PATIENT_LEVEL_ACCESSIONS:
            mapping_status = (
                "verified_library_to_donor"
                if donor_id
                else "unresolved_patient_mapping"
            )
            evidence = (
                "corrected GEO SOFT donor/patient field"
                if donor_id
                else "donor/patient identifier unavailable"
            )
        elif donor_id:
            mapping_status = "donor_available_not_primary_inference_unit"
            evidence = "metadata donor/patient field"
        else:
            mapping_status = "sample_level_only"
            evidence = "sample-level endpoint analysis; donor mapping not required"
        rows.append(
            {
                "accession": accession,
                "sample_id": sample_id,
                "library_id": library_id,
                "donor_id": donor_id,
                "condition": str(row.get("condition") or ""),
                "species": str(row.get("species") or ""),
                "tissue": str(row.get("tissue") or ""),
                "biological_unit": (
                    "animal"
                    if accession in ANIMAL_LEVEL_ACCESSIONS
                    else "donor"
                    if accession in PATIENT_LEVEL_ACCESSIONS
                    else "sample"
                ),
                "mapping_status": mapping_status,
                "mapping_evidence": evidence,
                "included": str(row.get("included") or ""),
                "source": str(row.get("source") or ""),
            }
        )
    mapping = pd.DataFrame(rows)
    if mapping.empty:
        return mapping, errors
    key_columns = ["accession", "sample_id", "library_id"]
    duplicate_mask = mapping.duplicated(key_columns, keep=False)
    if duplicate_mask.any():
        errors.append("duplicate accession/sample/library rows in donor mapping")
    mapping["donor_conflict"] = False
    mapped = mapping.loc[mapping["donor_id"].astype(str).str.strip().ne("")]
    if not mapped.empty:
        conflicting: set[tuple[str, str]] = set()
        for (accession, donor_id), group in mapped.groupby(
            ["accession", "donor_id"],
            observed=True,
        ):
            conditions = set(
                group["condition"].astype(str).str.strip()
            )
            conditions.discard("")
            allowed_pair = (
                str(accession) == "GSE164441"
                and conditions.issubset({"tumor", "adjacent_normal"})
            )
            if len(conditions) > 1 and not allowed_pair:
                conflicting.add((str(accession), str(donor_id)))
        if conflicting:
            mapping["donor_conflict"] = [
                (str(row.accession), str(row.donor_id)) in conflicting
                for row in mapping.itertuples(index=False)
            ]
            errors.append(
                "one or more donors map to conflicting conditions"
            )
    return mapping, errors


def freeze_cohort_artifacts(
    output_root: Path,
    *,
    root: Path | None = None,
    final: bool = False,
    force: bool = False,
) -> dict[str, Path | dict[str, Any]]:
    """Freeze the cohort manifest, donor mapping and source authorizations."""
    root = root or _repo_root()
    plan_dir = ensure_dir(Path(output_root) / PLAN_DIR)
    cohort_path = plan_dir / "cohort_manifest.tsv"
    cohort_errors: list[str] = []
    if not cohort_path.is_file():
        if final:
            raise FileNotFoundError(
                f"cohort manifest not found: {cohort_path}"
            )
        cohort = pd.DataFrame(columns=COHORT_COLUMNS)
        cohort_errors.append(f"cohort manifest not found: {cohort_path}")
    else:
        cohort = pd.read_csv(cohort_path, sep="\t", dtype=str).fillna("")
    missing_columns = set(COHORT_COLUMNS) - set(cohort.columns)
    if missing_columns:
        cohort_errors.append(
            f"cohort manifest missing columns: {sorted(missing_columns)}"
        )
    mapping, mapping_errors = build_donor_mapping(cohort)
    cohort_errors.extend(mapping_errors)
    mapping_path = plan_dir / "donor_mapping.tsv"
    mapping.to_csv(mapping_path, sep="\t", index=False)
    authorization_path, authorization = freeze_source_authorizations(
        output_root,
        root=root,
    )
    if final:
        if "included" in cohort and cohort["included"].astype(str).eq(
            "planned"
        ).any():
            cohort_errors.append(
                "cohort manifest still contains planned placeholder rows"
            )
        present = set(cohort["accession"].astype(str)) if "accession" in cohort else set()
        missing_accessions = sorted(REQUIRED_FINAL_ACCESSIONS - present)
        if missing_accessions:
            cohort_errors.append(
                f"cohort manifest lacks required accessions: {missing_accessions}"
            )
        if not mapping.empty:
            unresolved_patient = mapping[
                mapping["accession"].astype(str).isin(
                    PATIENT_LEVEL_ACCESSIONS
                )
                & mapping["donor_id"].astype(str).str.strip().eq("")
            ]
            if not unresolved_patient.empty:
                cohort_errors.append(
                    "patient-level cohorts still contain unresolved donor mappings"
                )
    authorization_errors = [
        str(value) for value in authorization.get("errors", [])
    ]
    frozen_manifest_path = plan_dir / "cohort_manifest.frozen.tsv"
    existing_frozen_hash = (
        sha256_file(frozen_manifest_path)
        if frozen_manifest_path.is_file()
        else ""
    )
    current_manifest_hash = (
        sha256_file(cohort_path) if cohort_path.is_file() else ""
    )
    if (
        existing_frozen_hash
        and existing_frozen_hash != current_manifest_hash
        and not force
    ):
        cohort_errors.append(
            "cohort manifest changed after an existing frozen snapshot"
        )
    cohort_complete = bool(final and not cohort_errors)
    donor_mapping_complete = bool(
        cohort_complete
        and not mapping_errors
        and (
            mapping.empty
            or not (
                mapping["accession"].astype(str).isin(
                    PATIENT_LEVEL_ACCESSIONS
                )
                & mapping["donor_id"].astype(str).str.strip().eq("")
            ).any()
        )
    )
    donor_mapping_limitations: list[str] = []
    if not mapping.empty:
        animal_unresolved = mapping[
            mapping["accession"].astype(str).isin(ANIMAL_LEVEL_ACCESSIONS)
            & mapping["donor_id"].astype(str).str.strip().eq("")
        ]
        if not animal_unresolved.empty:
            donor_mapping_limitations.append(
                "animal-level libraries have unresolved independent IDs; "
                "disease inference remains descriptive"
            )
    authorization_complete = (
        str(authorization.get("status") or "") == "complete"
    )
    if final and cohort_complete and authorization_complete:
        status = "frozen_complete"
    elif final and cohort_complete:
        status = "partially_frozen_authorization_blocked"
    elif final:
        status = "blocked_incomplete"
    elif not cohort_errors and not authorization_errors:
        status = "draft_valid"
    else:
        status = "draft_incomplete"
    validation_errors = [
        *cohort_errors,
        *authorization_errors,
    ]
    freeze_payload = {
        "schema_version": 1,
        "status": status,
        "final": bool(final),
        "cohort_status": (
            "frozen_complete" if cohort_complete else "blocked_incomplete"
        ),
        "donor_mapping_status": (
            "complete" if donor_mapping_complete else "incomplete"
        ),
        "donor_mapping_limitations": donor_mapping_limitations,
        "authorization_status": (
            "complete" if authorization_complete else "blocked_incomplete"
        ),
        "cohort_manifest": str(cohort_path),
        "cohort_manifest_sha256": current_manifest_hash,
        "cohort_manifest_frozen": (
            str(frozen_manifest_path) if cohort_complete else ""
        ),
        "cohort_manifest_frozen_sha256": (
            current_manifest_hash if cohort_complete else ""
        ),
        "donor_mapping": str(mapping_path),
        "donor_mapping_sha256": sha256_file(mapping_path),
        "source_authorizations": str(authorization_path),
        "source_authorization_status": authorization.get("status"),
        "source_authorization_errors": authorization.get("errors") or [],
        "rows": int(len(cohort)),
        "sample_level_rows": (
            int(cohort["sample_id"].astype(str).str.strip().ne("").sum())
            if "sample_id" in cohort
            else 0
        ),
        "accessions": (
            cohort["accession"].value_counts().to_dict()
            if "accession" in cohort
            else {}
        ),
        "donor_mapping_status_counts": (
            mapping["mapping_status"].value_counts().to_dict()
            if not mapping.empty
            else {}
        ),
        "cohort_validation_errors": cohort_errors,
        "authorization_validation_errors": authorization_errors,
        "validation_errors": validation_errors,
    }
    freeze_path = plan_dir / (
        "cohort_freeze.json" if final else "cohort_freeze.draft.json"
    )
    write_json(freeze_path, freeze_payload)
    if cohort_complete:
        frozen_manifest_path.write_text(
            cohort_path.read_text(encoding="utf-8"),
            encoding="utf-8",
        )
        freeze_payload["cohort_manifest_frozen_sha256"] = sha256_file(
            frozen_manifest_path
        )
        write_json(freeze_path, freeze_payload)
    return {
        "cohort_freeze": freeze_path,
        "donor_mapping": mapping_path,
        "source_authorizations": authorization_path,
        "validation": freeze_payload,
    }


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
    cohort_freeze = freeze_cohort_artifacts(
        output_root,
        final=False,
    )
    return {
        "analysis_plan": plan_path,
        "cohort_manifest": cohort_path,
        "sample_inclusion_exclusion": inclusion_path,
        "metadata_corrections": corrections_path,
        "data_checksums": checksums_path,
        "cohort_freeze": Path(cohort_freeze["cohort_freeze"]),
        "donor_mapping": Path(cohort_freeze["donor_mapping"]),
        "source_authorizations": Path(
            cohort_freeze["source_authorizations"]
        ),
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
        sample_level_accessions = set(
            incoming.loc[
                incoming["sample_id"].astype(str).str.strip().ne("")
                | incoming["library_id"].astype(str).str.strip().ne(""),
                "accession",
            ].astype(str)
        )
        placeholder_mask = (
            existing["accession"].astype(str).isin(sample_level_accessions)
            & existing["sample_id"].astype(str).str.strip().eq("")
            & existing["library_id"].astype(str).str.strip().eq("")
        )
        existing = existing.loc[~placeholder_mask].copy()
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
