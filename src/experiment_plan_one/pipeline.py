"""End-to-end orchestration for the experiment-plan-one workflow."""

from __future__ import annotations

import hashlib
import html
import json
import platform
import re
import shutil
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from evidence import (
    EvidenceContext,
    EvidenceHub,
    EvidenceRecord,
    EvidenceTier,
    SQLiteEvidenceStore,
)

from . import __version__
from .bulk import (
    GSE164441_COUNT_URL,
    SERIES_MATRIX_URLS,
    candidate_heatmap,
    differential_expression_limma,
    differential_summary_table,
    prepare_bulk_data,
    validation_boxplots,
)
from .classify import classify_experiment_plan_results
from .coexpression import run_coexpression_analysis
from .common import (
    FIGURE_PALETTE,
    LOG,
    STATUS_COLORS,
    configure_logging,
    download_file,
    ensure_dir,
    extract_tar,
    read_json,
    save_figure,
    sha256_file,
    slug,
    write_json,
)
from .coverage import audit_plan_coverage, audit_plan_environment
from .docking_md import prepare_or_run_md, run_docking_for_targets
from .enrichment import run_go_kegg
from .figure_audit import audit_figures
from .md_figures import generate_plan_md_figures
from .ml import run_ml_validation
from .planning import (
    append_cohort_rows,
    audit_delivery_readiness,
    dataset_registry,
    load_governance,
    update_data_checksums,
    write_analysis_plan,
    write_governance_artifacts,
)
from .ppi import run_ppi_analysis
from .single_cell import (
    run_human_single_cell,
    run_mouse_single_cell,
)
from .targets import (
    clinvar_disease_targets,
    collect_compound_targets,
    combine_disease_sources,
    compound_properties,
    gwas_catalog_disease_targets,
    load_compound_target_file,
    load_disease_source_file,
    make_venn_figure,
    make_workflow_figure,
    open_targets_disease_targets,
    rdkit_descriptors,
    write_compound_figures,
)

STAGES = (
    "plan",
    "data",
    "targets",
    "disease",
    "evidence",
    "ppi",
    "bulk",
    "ml",
    "mouse",
    "human",
    "docking",
    "md",
    "classify",
    "figure_audit",
    "report",
)

ARCHIVES = {
    "GSE270583_RAW.tar": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/"
            "GSE270583/suppl/GSE270583_RAW.tar"
        ),
        "size": 910_407_680,
        "extract_to": "GSE270583",
    },
    "GSE202379_RAW.tar": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/"
            "GSE202379/suppl/GSE202379_RAW.tar"
        ),
        "size": 595_322_880,
        "extract_to": "GSE202379",
    },
    "GSE135251_RAW.tar": {
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/"
            "GSE135251/suppl/GSE135251_RAW.tar"
        ),
        "size": 45_854_720,
        "extract_to": "GSE135251",
    },
}

SOFT_URLS = {
    "GSE270583_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE270nnn/"
        "GSE270583/soft/GSE270583_family.soft.gz"
    ),
    "GSE202379_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/"
        "GSE202379/soft/GSE202379_family.soft.gz"
    ),
    "GSE135251_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/"
        "GSE135251/soft/GSE135251_family.soft.gz"
    ),
}


@dataclass
class PipelineContext:
    root: Path
    config: dict[str, Any]
    output_root: Path
    raw_dir: Path
    processed_dir: Path
    state_dir: Path
    log_dir: Path
    force: bool = False

    def dir(self, name: str) -> Path:
        return ensure_dir(self.output_root / name)

    def fingerprint_file(self, path: Path) -> str:
        path = Path(path)
        if not path.exists() or not path.is_file():
            return "missing"
        return f"sha256:{sha256_file(path)}"

    def stage_signature(self, stage: str, previous_stages: list[str]) -> str:
        payload: dict[str, Any] = {
            "stage": stage,
            "pipeline_version": __version__,
            "config": self.config,
            "previous": {},
            "inputs": {},
            "code": {
                str(path.relative_to(self.root)): self.fingerprint_file(path)
                for path in sorted((self.root / "src" / "experiment_plan_one").rglob("*"))
                if path.is_file() and path.suffix in {".py", ".R"}
            },
        }
        for previous in previous_stages:
            state_path = self.state_dir / f"{previous}.json"
            payload["previous"][previous] = self.fingerprint_file(state_path)
        for path in STAGE_INPUT_PATHS.get(stage, []):
            resolved = path
            if not resolved.is_absolute():
                resolved = self.output_root / resolved
            if resolved.is_dir():
                files = sorted(
                    item
                    for item in resolved.rglob("*")
                    if item.is_file()
                )
                payload["inputs"][str(resolved)] = {
                    str(file.relative_to(resolved)): self.fingerprint_file(file)
                    for file in files
                }
            else:
                payload["inputs"][str(resolved)] = self.fingerprint_file(resolved)
        if stage == "evidence":
            evidence_config = str(
                (self.config.get("evidence") or {}).get("config_file")
                or "config/evidence_sources.json"
            )
            evidence_path = Path(evidence_config).expanduser()
            if not evidence_path.is_absolute():
                evidence_path = self.root / evidence_path
            payload["inputs"]["evidence_config"] = self.fingerprint_file(
                evidence_path.resolve()
            )
        if stage == "plan":
            for filename in (
                "experiment_plan_one_freeze.json",
                "experiment_plan_one_delivery_tiers.json",
                "experiment_plan_one_issue_register.csv",
                "experiment_plan_one_method_changes.csv",
                "experiment_plan_one_verification.json",
            ):
                payload["inputs"][filename] = self.fingerprint_file(
                    self.root / "config" / filename
                )
        encoded = json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            default=str,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()


STAGE_INPUT_PATHS: dict[str, list[Path]] = {
    "plan": [],
    "data": [],
    "targets": [Path("00_data/cache")],
    "disease": [],
    "evidence": [
        Path("01_compound_characterization/compound_properties.csv"),
        Path("01_compound_characterization/compound_targets.csv"),
        Path("02_disease_targets/disease_targets.csv"),
    ],
    "ppi": [
        Path("03_intersection_ppi/compound_disease_overlap.csv"),
        Path("02b_evidence/target_priority.csv"),
    ],
    "bulk": [Path("00_data/raw"), Path("00_data/processed")],
    "ml": [Path("00_data/processed"), Path("04_bulk_training")],
    "mouse": [
        Path("00_data/raw/extracted/GSE270583"),
        Path("05_machine_learning/ml_core_genes.json"),
    ],
    "human": [
        Path("00_data/raw/extracted/GSE202379"),
        Path("05_machine_learning/ml_core_genes.json"),
    ],
    "docking": [
        Path("05_machine_learning/ml_core_genes.json"),
        Path("02b_evidence/target_priority.csv"),
    ],
    "md": [Path("08_docking")],
    "classify": [Path("10_reports")],
    "figure_audit": [Path("按方案分类")],
    "report": [Path("10_reports")],
}

STAGE_REQUIRED_OUTPUTS: dict[str, list[Path]] = {
    "plan": [
        Path("00_plan/analysis_plan.yaml"),
        Path("00_plan/cohort_manifest.tsv"),
        Path("00_plan/sample_inclusion_exclusion.tsv"),
        Path("00_plan/metadata_corrections.tsv"),
        Path("00_plan/data_checksums.json"),
        Path("00_plan/governance/experiment_plan_one_freeze.json"),
        Path("00_plan/governance/experiment_plan_one_delivery_tiers.json"),
        Path("00_plan/governance/experiment_plan_one_issue_register.csv"),
        Path("00_plan/governance/experiment_plan_one_method_changes.csv"),
        Path("00_plan/governance/experiment_plan_one_verification.json"),
    ],
    "data": [
        Path("00_data/raw/dataset_inventory.json"),
        Path("00_data/raw/dataset_receipts.json"),
        *[Path("00_data/raw") / filename for filename in ARCHIVES],
        *[Path("00_data/raw") / filename for filename in SOFT_URLS],
    ],
    "targets": [
        Path("01_compound_characterization/compound_properties.csv"),
        Path("01_compound_characterization/compound_targets.csv"),
    ],
    "disease": [Path("02_disease_targets/disease_targets.csv")],
    "evidence": [
        Path("02b_evidence/evidence_summary.json"),
        Path("02b_evidence/target_priority.csv"),
    ],
    "ppi": [Path("03_intersection_ppi/ppi_summary.json")],
    "bulk": [
        Path("04_bulk_training/bulk_summary.json"),
        Path("04_bulk_training/bulk_source_manifest.json"),
    ],
    "ml": [Path("05_machine_learning/ml_summary.json")],
    "mouse": [Path("06_single_cell_mouse/mouse_single_cell_summary.json")],
    "human": [Path("07_single_cell_human/human_single_cell_summary.json")],
    "docking": [Path("08_docking/docking_summary.json")],
    "md": [
        Path("09_md_mmpbsa/md_stage_summary.json"),
        Path("09_md_mmpbsa/md_plan_figures.json"),
    ],
    "classify": [Path("按方案分类/分类汇总.json")],
    "figure_audit": [
        Path("10_reports/figure_quality_audit/figure_quality_audit.json"),
        Path("10_reports/plan_coverage/plan_coverage.json"),
        Path("10_reports/plan_coverage/environment_audit.json"),
        Path("10_reports/governance/delivery_readiness.json"),
    ],
    "report": [Path("10_reports/experiment_plan_one_report.html")],
}


def default_config() -> dict[str, Any]:
    return {
        "completion_target_percent": 90.0,
        "run_id": None,
        "known_reference_targets": [
            "TNF",
            "IL1B",
            "IL6",
            "TP53",
            "MAPK14",
            "PTGS2",
        ],
        "novelty_exclusions": ["TNF", "IL1B", "IL6", "TP53", "MAPK14", "PTGS2"],
        "compound": {
            "name": "6PPD-Q",
            "pubchem_cid": "154926030",
            "target_databases": [
                "ChEMBL",
                "STITCH",
                "SwissTargetPrediction",
                "SEA",
            ],
            "target_prediction_file": None,
            "target_prediction_gene_column": None,
            "target_prediction_score_column": None,
            "target_prediction_files": [],
            "sea": {
                "p_value_threshold": 0.05,
                "min_zscore": 0.0,
                "max_targets": 100,
                "timeout_seconds": 900,
                "poll_interval_seconds": 5,
            },
        },
        "disease": {
            "name": "NAFLD",
            "open_targets_terms": [
                "nonalcoholic fatty liver disease",
                "metabolic dysfunction-associated steatotic liver disease",
            ],
            "open_targets_min_score": 0.01,
            "gene_cards_file": None,
            "omim_file": None,
            "ttd_file": None,
            "gene_column": None,
            "use_open_evidence_sources": True,
            "source_thresholds": {
                "GeneCards": "relevance_score_above_median_verify_licensed_export",
                "OMIM": "gene_map_or_morbid_map_entry",
                "TTD": "approved_or_clinical_marker_entry",
            },
        },
        "datasets": {
            "GSE89632": {
                "role": "training",
                "endpoint": "healthy_control_vs_NAFLD",
                "platform": "GPL14951",
                "condition_column": "condition",
            },
            "GSE49541": {
                "role": "separate_endpoint_fibrosis_association",
                "endpoint": "advanced_vs_mild_fibrosis",
                "validation_status": "not_same_endpoint_validation",
                "platform": "GPL570",
                "condition_column": "condition",
            },
            "GSE164441": {
                "role": "hcc_extension_analysis",
                "endpoint": "tumor_vs_patient_matched_adjacent_non_tumor",
                "validation_status": "different_endpoint",
                "condition_column": "condition",
                "pair_column": "patient",
            },
            "GSE135251": {
                "role": "same_endpoint_external_candidate",
                "endpoint": "healthy_vs_NAFLD",
                "validation_status": "previously_exposed_development_support",
                "previously_used_for_model_development": True,
                "condition_column": "condition",
                "expected_condition_counts": {
                    "control": 10,
                    "NAFLD": 206,
                },
            },
            "GSE270583": {
                "role": "mouse_single_cell_discovery",
                "species": "mm",
                "validation_status": "insufficient_biological_replication",
            },
            "GSE202379": {
                "role": "human_single_cell_validation",
                "species": "hs",
                "biological_unit": "donor",
                "metadata_revision": "2025-06-06",
            },
        },
        "ppi": {
            "required_score": 700,
            "top_n": 20,
            "add_nodes": 50,
            "module_method": "mcl",
            "mcl_inflation": 2.0,
        },
        "ml": {"cv_folds": 5, "cv_repeats": 5, "seed": 42},
        "coexpression": {
            "enabled": True,
            "max_genes": 2000,
            "min_module_size": 20,
            "distance_threshold": 0.5,
            "kme_threshold": 0.5,
            "top_hubs_per_module": 25,
        },
        "evidence": {
            "enabled": True,
            "config_file": "config/evidence_sources.json",
            "target_scope": "intersection",
            "max_target_symbols": 500,
            "top_n": 50,
            "benchmark_positive": [],
            "benchmark_negative": [],
            "benchmark_source": "",
            "benchmark_version": "",
            "benchmark_independent": False,
            "benchmark_exclude_sources": [],
        },
        "single_cell": {
            "mouse": {
                "max_cells": 5000,
                "cellchat_backend": "r_cellchat",
                "cellchat_max_cells": 5000,
                "cellchat_permutations": 100,
                "cellchat_seed": 123,
                "inference_unit": "library_or_verified_animal",
            },
            "human": {
                "max_cells_per_sample": 0,
                "seed": 42,
                "inference_unit": "donor",
            },
            "insilico": {
                "engine": "sctenifoldknk_r",
                "scientific_role": "predicted_perturbation_response",
                "perturbation_scope": "pre_specified_primary_target",
                "scTenifold_n_networks": 5,
                "scTenifold_n_cells": 500,
                "scTenifold_n_comp": 3,
                "scTenifold_q": 0.95,
                "scTenifold_K": 3,
                "scTenifold_ma_dim": 2,
                "scTenifold_jobs": 1,
                "scTenifold_timeout_seconds": 86400,
            },
        },
        "docking": {
            "targets": 5,
            "exhaustiveness": 100,
            "num_modes": 20,
            "replicates": 3,
            "cpu": 4,
            "timeout_seconds": 3600,
        },
        "md": {
            "run": False,
            "gpu": True,
            "cpu": 4,
            "protein_forcefield": "amber14sb",
            "timeout_seconds": 172800,
        },
    }


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    output = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(output.get(key), dict):
            output[key] = merge_config(output[key], value)
        else:
            output[key] = value
    return output


def load_config(path: Path | None) -> dict[str, Any]:
    config = default_config()
    if path is None:
        return config
    if not path.exists():
        raise FileNotFoundError(f"config not found: {path}")
    override = json.loads(path.read_text(encoding="utf-8"))
    return merge_config(config, override)


class ExperimentPlanOne:
    """Stage runner with a small JSON state file per completed stage."""

    def __init__(
        self,
        output_root: Path,
        config: dict[str, Any],
        *,
        force: bool = False,
        verbose: bool = False,
    ) -> None:
        self.context = PipelineContext(
            root=Path(__file__).resolve().parents[2],
            config=config,
            output_root=output_root.resolve(),
            raw_dir=(output_root / "00_data" / "raw").resolve(),
            processed_dir=(output_root / "00_data" / "processed").resolve(),
            state_dir=(output_root / "12_reports" / ".stages").resolve(),
            log_dir=(output_root / "logs").resolve(),
            force=force,
        )
        for path in (
            self.context.raw_dir,
            self.context.processed_dir,
            self.context.state_dir,
            self.context.log_dir,
        ):
            ensure_dir(path)
        configure_logging(self.context.log_dir / "experiment_plan_one.log", verbose)
        identity_path = self.context.output_root / "run_identity.json"
        existing_identity = read_json(identity_path, {})
        requested_run_id = str(self.context.config.get("run_id") or "").strip()
        if (
            requested_run_id
            and existing_identity.get("run_id")
            and str(existing_identity["run_id"]) != requested_run_id
        ):
            raise RuntimeError(
                "output_root already belongs to a different run_id; use a "
                "new output_root to avoid mixing runs"
            )
        self.run_id = (
            requested_run_id
            or str(existing_identity.get("run_id") or "")
            or f"run_{time.strftime('%Y%m%d_%H%M%S')}"
        )
        write_json(
            identity_path,
            {
                "run_id": self.run_id,
                "output_root": str(self.context.output_root),
                "created_at": existing_identity.get("created_at")
                or time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "config_fingerprint": hashlib.sha256(
                    json.dumps(
                        self.context.config,
                        ensure_ascii=False,
                        sort_keys=True,
                        default=str,
                    ).encode("utf-8")
                ).hexdigest(),
            },
        )
        self.results: dict[str, Any] = {}

    def run(self, stages: list[str]) -> dict[str, Any]:
        stage_functions: dict[str, Callable[[], dict[str, Any]]] = {
            "plan": self.stage_plan,
            "data": self.stage_data,
            "targets": self.stage_targets,
            "disease": self.stage_disease,
            "evidence": self.stage_evidence,
            "ppi": self.stage_ppi,
            "bulk": self.stage_bulk,
            "ml": self.stage_ml,
            "mouse": self.stage_mouse,
            "human": self.stage_human,
            "docking": self.stage_docking,
            "md": self.stage_md,
            "classify": self.stage_classify,
            "figure_audit": self.stage_figure_audit,
            "report": self.stage_report,
        }
        for stage in stages:
            if stage not in stage_functions:
                raise ValueError(f"unknown stage: {stage}; expected one of {STAGES}")
        for stage in stages:
            state_path = self.context.state_dir / f"{stage}.json"
            previous_stages = list(STAGES[: STAGES.index(stage)])
            signature = self.context.stage_signature(stage, previous_stages)
            missing_outputs = self._missing_stage_outputs(stage)
            successful_statuses = (
                {"completed", "prepared"} if stage == "md" else {"completed"}
            )
            if state_path.exists() and not self.context.force:
                state = read_json(state_path, {})
                same_signature = state.get("signature") == signature
                output_ready = not missing_outputs
                if (
                    state.get("status") in successful_statuses
                    and same_signature
                    and output_ready
                ):
                    LOG.info(
                        "stage %s already complete with matching signature; reusing %s",
                        stage,
                        state_path,
                    )
                    self.results[stage] = state
                    continue
                LOG.warning(
                    "stage %s will rerun: signature_match=%s missing_outputs=%s",
                    stage,
                    same_signature,
                    missing_outputs or "none",
                )
            started = time.time()
            LOG.info("starting stage %s", stage)
            result = stage_functions[stage]()
            result = _serializable(result)
            missing_outputs = self._missing_stage_outputs(stage)
            result_status = (
                str(result.get("status") or "completed")
                if isinstance(result, dict)
                else "completed"
            )
            if result_status == "failed" or missing_outputs:
                result_status = "failed"
            state = {
                "stage": stage,
                "status": result_status,
                "signature": signature,
                "elapsed_seconds": round(time.time() - started, 3),
                "result": result,
            }
            write_json(state_path, state)
            self.results[stage] = state
            if result_status not in successful_statuses:
                reason = (
                    f"missing outputs: {missing_outputs}"
                    if missing_outputs
                    else str(result.get("reason") or "stage returned failed")
                    if isinstance(result, dict)
                    else "stage returned failed"
                )
                raise RuntimeError(f"stage {stage} failed: {reason}")
            if result_status == "prepared":
                LOG.info(
                    "stage %s prepared successfully in %.1f s",
                    stage,
                    state["elapsed_seconds"],
                )
            else:
                LOG.info(
                    "completed stage %s in %.1f s",
                    stage,
                    state["elapsed_seconds"],
                )
        self._write_manifest()
        self._write_delivery_readiness()
        return self.results

    def _missing_stage_outputs(self, stage: str) -> list[str]:
        missing: list[str] = []
        for relative in STAGE_REQUIRED_OUTPUTS.get(stage, []):
            path = self.context.output_root / relative
            if not path.exists() or not path.is_file() or path.stat().st_size == 0:
                missing.append(relative.as_posix())
        return missing

    def stage_plan(self) -> dict[str, Any]:
        outputs = write_analysis_plan(
            self.context.output_root,
            {**self.context.config, "run_id": self.run_id},
        )
        outputs.update(
            write_governance_artifacts(self.context.output_root)
        )
        update_data_checksums(
            self.context.output_root,
            [
                self.context.root / "config" / "experiment_plan_one.json",
                *[
                    path
                    for key, path in outputs.items()
                    if key
                    in {
                        "script_freeze",
                        "delivery_tiers",
                        "issue_register",
                        "method_changes",
                        "verification",
                    }
                ],
            ],
        )
        return {
            "status": "completed",
            "outputs": {key: str(value) for key, value in outputs.items()},
            "note": (
                "The plan freezes endpoint roles and known exposure; it does "
                "not assert that external validation or wet-lab work is complete."
            ),
        }

    def stage_data(self) -> dict[str, Any]:
        outputs: dict[str, str] = {}
        receipts: dict[str, dict[str, Any]] = {}
        for filename, spec in ARCHIVES.items():
            archive = self.context.raw_dir / filename
            download_file(
                str(spec["url"]),
                archive,
                expected_size=int(spec["size"]),
                timeout=300,
            )
            extraction = self.context.raw_dir / "extracted" / str(spec["extract_to"])
            extract_tar(archive, extraction)
            outputs[filename] = str(archive)
            outputs[f"{spec['extract_to']}_extracted"] = str(extraction)
            receipts[filename] = {
                "url": spec["url"],
                "expected_size": int(spec["size"]),
                "size": archive.stat().st_size,
                "sha256": sha256_file(archive),
            }
        for filename, url in SOFT_URLS.items():
            destination = self.context.raw_dir / filename
            download_file(url, destination, timeout=180)
            outputs[filename] = str(destination)
            receipts[filename] = {
                "url": url,
                "size": destination.stat().st_size,
                "sha256": sha256_file(destination),
            }
        # Bulk series matrices and platforms are handled by their dedicated stage.
        write_json(self.context.raw_dir / "dataset_inventory.json", outputs)
        write_json(self.context.raw_dir / "dataset_receipts.json", receipts)
        update_data_checksums(
            self.context.output_root,
            [
                self.context.raw_dir / filename
                for filename in (*ARCHIVES.keys(), *SOFT_URLS.keys())
            ],
        )
        return {"outputs": outputs}

    def stage_targets(self) -> dict[str, Any]:
        out_dir = self.context.dir("01_compound_characterization")
        retrieved_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        cid = str(self.context.config["compound"]["pubchem_cid"])
        properties = compound_properties(cid)
        descriptors = rdkit_descriptors(properties["canonical_smiles"])
        compound = {**properties, **descriptors}
        pd.DataFrame([compound]).to_csv(out_dir / "compound_properties.csv", index=False)
        figures = write_compound_figures(
            properties["canonical_smiles"],
            properties,
            descriptors,
            out_dir,
        )
        figures["workflow"] = make_workflow_figure(out_dir / "fig1a_workflow.png")
        targets, statuses = collect_compound_targets(
            properties,
            out_dir,
            cache_dir=self.context.output_root / "00_data" / "cache",
            enabled_sources=list(
                self.context.config["compound"].get("target_databases") or []
            ),
            sea_config=self.context.config["compound"].get("sea") or {},
        )
        source_dir = out_dir / "sources"
        extra_prediction, extra_status = load_compound_target_file(
            "LocalPrediction",
            self.context.config["compound"].get("target_prediction_file"),
            gene_column=self.context.config["compound"].get(
                "target_prediction_gene_column"
            ),
            score_column=self.context.config["compound"].get(
                "target_prediction_score_column"
            ),
        )
        statuses["LocalPrediction"] = extra_status
        local_prediction_written = not extra_prediction.empty
        if not extra_prediction.empty:
            extra_prediction.to_csv(
                source_dir / "LocalPrediction.csv",
                index=False,
            )
        for index, spec in enumerate(
            self.context.config["compound"].get("target_prediction_files") or [],
            start=1,
        ):
            if not isinstance(spec, dict):
                LOG.warning(
                    "ignoring target_prediction_files[%s]: expected an object",
                    index,
                )
                continue
            name = str(spec.get("name") or f"LocalPrediction{index}").strip()
            frame, status = load_compound_target_file(
                name,
                spec.get("path"),
                gene_column=spec.get("gene_column"),
                score_column=spec.get("score_column"),
            )
            statuses[name] = status
            if not frame.empty:
                frame.to_csv(
                    source_dir / f"{slug(name)}.csv",
                    index=False,
                )
                local_prediction_written = True
        for name, status in statuses.items():
            if isinstance(status, dict):
                status.setdefault("retrieved_at", retrieved_at)
                status.setdefault("species", "Homo sapiens")
                status.setdefault("evidence_type", "compound-target association")
                if name == "SEA":
                    status.setdefault(
                        "threshold",
                        self.context.config["compound"].get("sea") or {},
                    )
        pd.DataFrame(
            [
                {
                    "source": name,
                    **(
                        dict(status)
                        if isinstance(status, dict)
                        else {"status": str(status)}
                    ),
                }
                for name, status in statuses.items()
            ]
        ).to_csv(out_dir / "target_evidence_provenance.csv", index=False)
        if local_prediction_written:
            source_frames: dict[str, pd.DataFrame] = {}
            for path in sorted(source_dir.glob("*.csv")):
                frame = pd.read_csv(path)
                if frame.empty or "gene" not in frame.columns:
                    continue
                score_column = next(
                    (
                        column
                        for column in ("score", "probability", "best_score")
                        if column in frame.columns
                    ),
                    None,
                )
                if score_column is None:
                    frame["score"] = 1.0
                    score_column = "score"
                source_frames[path.stem] = frame[
                    ["gene", score_column]
                ].rename(columns={score_column: "score"})
            targets = combine_disease_sources(source_frames)
            targets.to_csv(out_dir / "compound_targets.csv", index=False)
        source_sets: dict[str, set[str]] = {}
        for path in sorted(source_dir.glob("*.csv")):
            if not path.exists():
                continue
            frame = pd.read_csv(path)
            if frame.empty or "gene" not in frame.columns:
                continue
            genes = set(frame["gene"].dropna().astype(str).str.upper())
            if genes:
                source_sets[path.stem] = genes
        if len(source_sets) >= 2:
            (out_dir / "fig1d_compound_target_source_counts.png").unlink(
                missing_ok=True
            )
            make_venn_figure(
                source_sets,
                out_dir / "fig1d_compound_target_source_venn.png",
                title="6PPD-Q compound-target sources",
            )
        else:
            (out_dir / "fig1d_compound_target_source_venn.png").unlink(
                missing_ok=True
            )
            _plot_target_source_counts(
                source_sets,
                out_dir / "fig1d_compound_target_source_counts.png",
            )
        write_json(
            out_dir / "target_collection_summary.json",
            {
                "compound": properties,
                "compound_descriptors": descriptors,
            "source_status": statuses,
            "provenance": str(out_dir / "target_evidence_provenance.csv"),
            "target_count": len(targets),
            },
        )
        return {
            "compound_properties": str(out_dir / "compound_properties.csv"),
            "compound_targets": str(out_dir / "compound_targets.csv"),
            "figures": {key: str(value) for key, value in figures.items()},
            "n_targets": len(targets),
            "source_status": statuses,
        }

    def stage_disease(self) -> dict[str, Any]:
        out_dir = self.context.dir("02_disease_targets")
        retrieved_at = time.strftime("%Y-%m-%dT%H:%M:%S%z")
        sources: dict[str, pd.DataFrame] = {}
        statuses: dict[str, Any] = {}
        for term in self.context.config["disease"].get("open_targets_terms") or [
            self.context.config["disease"]["name"]
        ]:
            try:
                frame, status = open_targets_disease_targets(
                    str(term),
                    page_size=500,
                    min_score=float(
                        self.context.config["disease"].get("open_targets_min_score", 0.0)
                    ),
                )
            except Exception as exc:  # noqa: BLE001
                frame, status = pd.DataFrame(), {"status": "failed", "reason": str(exc)}
            statuses[f"OpenTargets:{term}"] = status
            if not frame.empty:
                frame = frame.copy()
                source_name = f"OpenTargets:{status.get('disease_id') or slug(str(term))}"
                frame["source"] = source_name
                sources[source_name] = frame
        for name, key in (
            ("GeneCards", "gene_cards_file"),
            ("OMIM", "omim_file"),
            ("TTD", "ttd_file"),
        ):
            frame, status = load_disease_source_file(
                name,
                self.context.config["disease"].get(key),
                gene_column=self.context.config["disease"].get("gene_column"),
            )
            statuses[name] = status
            if not frame.empty:
                sources[name] = frame
        if (
            bool(
                self.context.config["disease"].get(
                    "use_open_evidence_sources",
                    True,
                )
            )
            and bool(
                (self.context.config.get("evidence") or {}).get(
                    "allow_network",
                    True,
                )
            )
        ):
            terms = [
                str(term)
                for term in (
                    self.context.config["disease"].get("open_targets_terms")
                    or [self.context.config["disease"]["name"]]
                )
                if str(term).strip()
            ]
            gwas_frame, gwas_status = gwas_catalog_disease_targets(terms)
            statuses["GWAS_Catalog"] = gwas_status
            if not gwas_frame.empty:
                sources["GWAS_Catalog"] = gwas_frame
            seed_genes = list(
                dict.fromkeys(
                    str(value)
                    for frame in sources.values()
                    for value in frame["gene"].astype(str).tolist()
                    if str(value).strip()
                )
            )
            clinvar_frame, clinvar_status = clinvar_disease_targets(
                seed_genes[:120]
            )
            statuses["ClinVar"] = clinvar_status
            if not clinvar_frame.empty:
                sources["ClinVar"] = clinvar_frame
        for name, status in statuses.items():
            if isinstance(status, dict):
                status.setdefault("retrieved_at", retrieved_at)
                status.setdefault("species", "Homo sapiens")
                status.setdefault("disease", self.context.config["disease"]["name"])
                status.setdefault(
                    "threshold",
                    (
                        self.context.config["disease"].get(
                            "source_thresholds",
                            {},
                        )
                        or {}
                    ).get(name, "source-specific exported fields"),
                )
        pd.DataFrame(
            [
                {
                    "source": name,
                    **(
                        dict(status)
                        if isinstance(status, dict)
                        else {"status": str(status)}
                    ),
                }
                for name, status in statuses.items()
            ]
        ).to_csv(out_dir / "disease_target_provenance.csv", index=False)
        if not sources:
            raise RuntimeError(
                "no disease-target source could be loaded; provide local "
                "GeneCards/OMIM/TTD tables or check Open Targets access"
            )
        combined = combine_disease_sources(sources)
        combined.to_csv(out_dir / "disease_targets.csv", index=False)
        for name, frame in sources.items():
            frame.to_csv(out_dir / f"source_{slug(name)}.csv", index=False)
        if len(sources) >= 2:
            (out_dir / "fig1e_disease_target_source_status.png").unlink(
                missing_ok=True
            )
            make_venn_figure(
                {name: set(frame["gene"]) for name, frame in sources.items()},
                out_dir / "fig1e_disease_target_venn.png",
                title="NAFLD disease target sources",
            )
        else:
            (out_dir / "fig1e_disease_target_venn.png").unlink(missing_ok=True)
            _plot_source_coverage(statuses, out_dir / "fig1e_disease_target_source_status.png")
        write_json(
            out_dir / "disease_target_summary.json",
            {
                "sources": statuses,
                "counts": {name: len(frame) for name, frame in sources.items()},
                "combined_targets": len(combined),
                "note": (
                    "GeneCards, OMIM and TTD bulk downloads require licensed or "
                    "credentialed access unless local exports are supplied."
                ),
                "provenance": str(out_dir / "disease_target_provenance.csv"),
            },
        )
        return {
            "disease_targets": str(out_dir / "disease_targets.csv"),
            "sources": statuses,
            "n_targets": len(combined),
        }

    def stage_evidence(self) -> dict[str, Any]:
        """Collect independent evidence sources and produce target priorities."""
        out_dir = self.context.dir("02b_evidence")
        config_section = dict(self.context.config.get("evidence") or {})
        if not config_section.get("enabled", True):
            summary = {"status": "skipped", "reason": "evidence stage disabled"}
            pd.DataFrame(
                columns=[
                    "rank",
                    "target_symbol",
                    "priority_score",
                    "coverage_ratio",
                ]
            ).to_csv(out_dir / "target_priority.csv", index=False)
            self._write_priority_markdown(
                pd.DataFrame(),
                out_dir / "target_priority.md",
            )
            write_json(out_dir / "evidence_summary.json", summary)
            return summary

        config_file = str(
            config_section.get("config_file")
            or "config/evidence_sources.json"
        )
        source_config_path = Path(config_file).expanduser()
        if not source_config_path.is_absolute():
            source_config_path = (self.context.root / source_config_path).resolve()
        evidence_config = (
            json.loads(source_config_path.read_text(encoding="utf-8"))
            if source_config_path.exists()
            else {}
        )
        overrides = {
            key: value
            for key, value in config_section.items()
            if key
            not in {
                "enabled",
                "config_file",
                "target_scope",
                "top_n",
                "benchmark_positive",
                "benchmark_negative",
                "benchmark_source",
                "benchmark_version",
                "benchmark_independent",
                "benchmark_exclude_sources",
            }
        }
        evidence_config = merge_config(evidence_config, overrides)

        compound_dir = (
            self.context.output_root / "01_compound_characterization"
        )
        disease_dir = self.context.output_root / "02_disease_targets"
        compound_frame = pd.read_csv(compound_dir / "compound_targets.csv")
        disease_frame = pd.read_csv(disease_dir / "disease_targets.csv")
        compound_genes = {
            str(value).upper()
            for value in compound_frame.get("gene", pd.Series(dtype=str)).dropna()
        }
        disease_genes = {
            str(value).upper()
            for value in disease_frame.get("gene", pd.Series(dtype=str)).dropna()
        }
        scope = str(config_section.get("target_scope") or "intersection").lower()
        if scope == "all":
            target_symbols = sorted(compound_genes | disease_genes)
        elif scope == "compound":
            target_symbols = sorted(compound_genes)
        elif scope == "disease":
            target_symbols = sorted(disease_genes)
        else:
            target_symbols = sorted(compound_genes & disease_genes)
            if not target_symbols:
                target_symbols = sorted(compound_genes | disease_genes)
        max_targets = max(1, int(config_section.get("max_target_symbols", 500)))
        target_symbols = target_symbols[:max_targets]

        ensembl_map: dict[str, str] = {}
        if "ensembl_id" in disease_frame.columns:
            for row in disease_frame[["gene", "ensembl_id"]].dropna().to_dict(
                "records"
            ):
                gene = str(row.get("gene") or "").upper().strip()
                ensembl = str(row.get("ensembl_id") or "").strip()
                if gene and ensembl:
                    ensembl_map[gene] = ensembl
        compound_row = pd.read_csv(
            compound_dir / "compound_properties.csv"
        ).iloc[0].to_dict()
        disease_row = {"name": self.context.config["disease"]["name"]}
        context = EvidenceContext(
            compound=compound_row,
            disease=disease_row,
            target_symbols=target_symbols,
            ensembl_ids=ensembl_map,
            cache_dir=out_dir,
            max_records_per_source=int(
                (evidence_config.get("limits") or {}).get(
                    "max_records_per_source",
                    1000,
                )
            ),
            timeout_seconds=int(
                (evidence_config.get("limits") or {}).get(
                    "timeout_seconds",
                    120,
                )
            ),
            allow_network=bool(config_section.get("allow_network", True)),
            source_options={
                str(name): dict(options or {})
                for name, options in (evidence_config.get("sources") or {}).items()
            },
        )
        database = self.context.output_root / "00_data" / "evidence" / "target_evidence.sqlite"
        with SQLiteEvidenceStore(database) as store:
            inherited_started = time.time()
            inherited = self._inherited_target_evidence(
                compound_row,
                str(disease_row["name"]),
            )
            store.upsert_evidence(inherited)
            if inherited:
                store.record_source_run(
                    source="ExperimentPlanInherited",
                    source_version="current stage outputs",
                    status="completed",
                    started_at=time.strftime(
                        "%Y-%m-%dT%H:%M:%S%z",
                        time.localtime(inherited_started),
                    ),
                    completed_at=time.strftime(
                        "%Y-%m-%dT%H:%M:%S%z",
                        time.localtime(time.time()),
                    ),
                    record_count=len(inherited),
                    query_hash="local-stage-outputs",
                    metadata={"mode": "local import"},
                )
            hub = EvidenceHub.from_config(
                database,
                evidence_config,
                store=store,
            )
            status = hub.collect(context)
            scored = hub.score(
                benchmark_positives=config_section.get("benchmark_positive") or [],
                benchmark_negatives=config_section.get("benchmark_negative") or [],
                benchmark_top_n=int(config_section.get("top_n", 50)),
                benchmark_source=str(
                    config_section.get("benchmark_source") or ""
                ),
                benchmark_version=str(
                    config_section.get("benchmark_version") or ""
                ),
                benchmark_independent=bool(
                    config_section.get("benchmark_independent", False)
                ),
                benchmark_exclude_sources=(
                    config_section.get("benchmark_exclude_sources") or []
                ),
            )
            paths = hub.export(
                out_dir,
                context=context,
                source_status=status,
                scored=scored,
            )
        priority = scored["priority"]
        summary = {
            "status": "completed",
            "database": str(database),
            "target_scope": scope,
            "requested_targets": len(target_symbols),
            "evidence_records": len(scored["evidence"]),
            "prioritized_targets": len(priority),
            "collector_status": status.to_dict(orient="records"),
            "benchmark": scored.get("benchmark"),
            **paths,
        }
        write_json(out_dir / "evidence_summary.json", summary)
        self._write_priority_markdown(priority, out_dir / "target_priority.md")
        return summary

    def _inherited_target_evidence(
        self,
        compound: dict[str, Any],
        disease_name: str,
    ) -> list[EvidenceRecord]:
        records: list[EvidenceRecord] = []
        compound_id = str(
            compound.get("pubchem_cid")
            or compound.get("inchi_key")
            or "compound"
        )
        compound_sources = {
            "SwissTargetPrediction": (EvidenceTier.PREDICTED, "predicted_target"),
            "STITCH": (EvidenceTier.PREDICTED, "predicted_target"),
            "ChEMBL_similarity": (EvidenceTier.PREDICTED, "ligand_similarity"),
        }
        source_dir = (
            self.context.output_root
            / "01_compound_characterization"
            / "sources"
        )
        for source_name, (tier, evidence_type) in compound_sources.items():
            path = source_dir / f"{source_name}.csv"
            if not path.exists():
                continue
            frame = pd.read_csv(path)
            if "gene" not in frame.columns:
                continue
            score_column = "probability" if "probability" in frame.columns else None
            for index, row in frame.iterrows():
                symbol = str(row.get("gene") or "").upper().strip()
                if not symbol:
                    continue
                records.append(
                    EvidenceRecord(
                        source=source_name,
                        source_record_id=f"pipeline:{index}:{symbol}",
                        evidence_type=evidence_type,
                        subject_type="compound",
                        subject_id=compound_id,
                        relation="targets",
                        object_type="target",
                        object_id=symbol,
                        target_symbol=symbol,
                        tier=tier,
                        source_version="experiment-plan-one inherited output",
                        source_group=source_name.lower(),
                        score=(
                            float(row.get(score_column))
                            if score_column
                            and pd.notna(row.get(score_column))
                            else None
                        ),
                        species="Homo sapiens",
                        payload={"inherited_from": str(path)},
                    )
                )
        disease_dir = self.context.output_root / "02_disease_targets"
        for path in sorted(disease_dir.glob("source_*.csv")):
            source_name = path.stem.removeprefix("source_")
            if source_name.lower().startswith("opentargets"):
                continue
            frame = pd.read_csv(path)
            if "gene" not in frame.columns:
                continue
            for index, row in frame.iterrows():
                symbol = str(row.get("gene") or "").upper().strip()
                if not symbol:
                    continue
                records.append(
                    EvidenceRecord(
                        source=source_name,
                        source_record_id=f"pipeline:{index}:{symbol}",
                        evidence_type="disease_gene_association",
                        subject_type="disease",
                        subject_id=disease_name,
                        relation="associated_with",
                        object_type="target",
                        object_id=symbol,
                        target_symbol=symbol,
                        tier=EvidenceTier.CURATED,
                        source_version="local or licensed snapshot",
                        source_group=source_name.lower(),
                        score=(
                            float(row.get("score"))
                            if "score" in frame.columns
                            and pd.notna(row.get("score"))
                            else 1.0
                        ),
                        species="Homo sapiens",
                        payload={"inherited_from": str(path)},
                    )
                )
        return records

    @staticmethod
    def _write_priority_markdown(priority: pd.DataFrame, output: Path) -> None:
        lines = [
            "# Target evidence priorities",
            "",
            "Scores use only available evidence categories. Missing database "
            "records are not treated as negative evidence.",
            "",
            "| Rank | Target | Priority | Coverage | Sources | Missing categories |",
            "|---:|---|---:|---:|---:|---|",
        ]
        for row in priority.head(50).to_dict("records"):
            lines.append(
                "| {rank} | {target} | {score:.3f} | {coverage:.2f} | "
                "{sources} | {missing} |".format(
                    rank=int(row.get("rank") or 0),
                    target=row.get("target_symbol") or "",
                    score=float(row.get("priority_score") or 0.0),
                    coverage=float(row.get("coverage_ratio") or 0.0),
                    sources=int(row.get("source_group_count") or 0),
                    missing=row.get("missing_categories") or "",
                )
            )
        if len(priority) == 0:
            lines.append("| - | No target with usable evidence | 0 | 0 | 0 | - |")
        output.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def stage_ppi(self) -> dict[str, Any]:
        out_dir = self.context.dir("03_intersection_ppi")
        compound = pd.read_csv(
            self.context.output_root / "01_compound_characterization" / "compound_targets.csv"
        )
        disease = pd.read_csv(
            self.context.output_root / "02_disease_targets" / "disease_targets.csv"
        )
        compound_genes = set(compound["gene"].astype(str).str.upper())
        disease_genes = set(disease["gene"].astype(str).str.upper())
        overlap = sorted(compound_genes & disease_genes)
        overlap_frame = pd.DataFrame(
            {
                "gene": overlap,
                "compound_source": [
                    ";".join(
                        compound.loc[compound["gene"] == gene, "sources"].astype(str)
                    )
                    for gene in overlap
                ],
                "disease_source": [
                    ";".join(
                        disease.loc[disease["gene"] == gene, "sources"].astype(str)
                    )
                    for gene in overlap
                ],
            }
        )
        excluded = {
            str(value).upper()
            for value in (
                self.context.config.get("known_reference_targets")
                or self.context.config.get("novelty_exclusions", [])
            )
        }
        overlap_frame["previously_reported_in_plan"] = overlap_frame["gene"].isin(
            excluded
        )
        overlap_frame["novelty_status"] = np.where(
            overlap_frame["previously_reported_in_plan"],
            "known_positive_reference_not_removed",
            "not_in_legacy_reference_list",
        )
        overlap_frame.to_csv(out_dir / "compound_disease_overlap.csv", index=False)
        overlap_frame.to_csv(out_dir / "target_novelty_assessment.csv", index=False)
        make_venn_figure(
            {"Compound targets": compound_genes, "Disease targets": disease_genes},
            out_dir / "fig1f_compound_disease_venn.png",
            title="6PPD-Q targets and NAFLD targets",
        )
        # Preserve a small/empty intersection; never replace it with a union.
        analysis_genes = overlap
        ppi = run_ppi_analysis(
            analysis_genes,
            out_dir,
            required_score=int(self.context.config["ppi"]["required_score"]),
            top_n=int(self.context.config["ppi"]["top_n"]),
            add_nodes=int(self.context.config["ppi"].get("add_nodes", 50)),
            force=self.context.force,
            module_method=str(
                self.context.config["ppi"].get("module_method", "mcl")
            ),
            mcl_inflation=float(
                self.context.config["ppi"].get("mcl_inflation", 2.0)
            ),
        )
        enrichment = run_go_kegg(
            analysis_genes,
            out_dir / "enrichment",
            universe=sorted(disease_genes),  # conditional background: disease-associated candidate pool
            species="hs",
            force=self.context.force,
        )
        return {
            "overlap": str(out_dir / "compound_disease_overlap.csv"),
            "novelty_assessment": str(out_dir / "target_novelty_assessment.csv"),
            "n_overlap": len(overlap),
            "ppi": {key: _serializable(value) for key, value in ppi.items()},
            "enrichment": enrichment,
        }

    def stage_bulk(self) -> dict[str, Any]:
        out_dir = self.context.dir("04_bulk_training")
        paths = prepare_bulk_data(
            self.context.raw_dir,
            self.context.processed_dir,
            force=self.context.force,
        )
        gse89632_expr = Path(paths["GSE89632_expression"])
        gse89632_meta = Path(paths["GSE89632_metadata"])
        gse49541_expr = Path(paths["GSE49541_expression"])
        gse49541_meta = Path(paths["GSE49541_metadata"])
        gse164441_expr = Path(paths["GSE164441_expression"])
        gse164441_meta = Path(paths["GSE164441_metadata"])
        gse135251_expr = Path(paths["GSE135251_expression"])
        gse135251_meta = Path(paths["GSE135251_metadata"])
        registry = {
            row["accession"]: row for row in dataset_registry(self.context.config)
        }
        cohort_rows: list[dict[str, Any]] = []
        for accession, metadata_path in (
            ("GSE89632", gse89632_meta),
            ("GSE49541", gse49541_meta),
            ("GSE164441", gse164441_meta),
            ("GSE135251", gse135251_meta),
        ):
            metadata = pd.read_csv(metadata_path, index_col=0, dtype=str).fillna("")
            descriptor = registry.get(accession, {})
            donor_column = next(
                (
                    column
                    for column in ("patient", "patient_id", "donor_id", "donor")
                    if column in metadata.columns
                ),
                "",
            )
            diagnosis_column = next(
                (
                    column
                    for column in (
                        "diagnosis",
                        "disease",
                        "disease_status",
                        "original_diagnosis",
                    )
                    if column in metadata.columns
                ),
                "",
            )
            for sample_id, row in metadata.iterrows():
                condition = str(row.get("condition") or "")
                donor_id = str(row.get(donor_column) or "") if donor_column else ""
                cohort_rows.append(
                    {
                        "accession": accession,
                        "sample_id": str(sample_id),
                        "library_id": str(sample_id),
                        "donor_id": donor_id or str(sample_id),
                        "species": descriptor.get("species", ""),
                        "platform": descriptor.get("platform", ""),
                        "tissue": descriptor.get("tissue", "liver"),
                        "condition": condition,
                        "original_diagnosis": (
                            str(row.get(diagnosis_column) or "")
                            if diagnosis_column
                            else condition
                        ),
                        "paired_patient": "yes" if donor_column else "no",
                        "batch": str(row.get("batch") or row.get("platform_id") or ""),
                        "exposure_status": "unknown",
                        "included": "yes",
                        "inclusion_reason": descriptor.get(
                            "validation_status", descriptor.get("role", "")
                        ),
                        "data_version": "frozen_processed_matrix",
                        "source": str(metadata_path),
                    }
                )
        append_cohort_rows(self.context.output_root, cohort_rows)
        update_data_checksums(
            self.context.output_root,
            [
                *(
                    Path(path)
                    for path in paths.values()
                    if Path(path).is_file()
                ),
                *(
                    source_path
                    for source_path in (
                        self.context.raw_dir / name
                        for name in (
                            "GSE89632_series_matrix.txt.gz",
                            "GSE49541_series_matrix.txt.gz",
                            "GSE164441_count.txt.gz",
                            "GSE135251_RAW.tar",
                            "GSE135251_family.soft.gz",
                        )
                    )
                    if source_path.is_file()
                ),
            ],
        )

        bulk_sources = {
            "GSE89632_series_matrix": {
                "url": SERIES_MATRIX_URLS["GSE89632"],
                "path": str(self.context.raw_dir / "GSE89632_series_matrix.txt.gz"),
            },
            "GSE49541_series_matrix": {
                "url": SERIES_MATRIX_URLS["GSE49541"],
                "path": str(self.context.raw_dir / "GSE49541_series_matrix.txt.gz"),
            },
            "GSE164441_count": {
                "url": GSE164441_COUNT_URL,
                "path": str(self.context.raw_dir / "GSE164441_count.txt.gz"),
            },
            "GSE135251_raw": {
                "url": ARCHIVES["GSE135251_RAW.tar"]["url"],
                "path": str(self.context.raw_dir / "GSE135251_RAW.tar"),
            },
            "GSE135251_soft": {
                "url": SOFT_URLS["GSE135251_family.soft.gz"],
                "path": str(
                    self.context.raw_dir / "GSE135251_family.soft.gz"
                ),
            },
        }
        for record in bulk_sources.values():
            source_path = Path(record["path"])
            record.update({"size": source_path.stat().st_size, "sha256": sha256_file(source_path)})
        write_json(out_dir / "bulk_source_manifest.json", bulk_sources)

        hc_vs_nafld = differential_expression_limma(
            gse89632_expr,
            gse89632_meta,
            out_dir / "gse89632_HC_vs_NAFLD_limma.csv",
            condition_column="condition",
            group1=["HC"],
            group2=["SS", "NASH"],
            comparison="HC_vs_NAFLD",
        )
        ss_vs_nash = differential_expression_limma(
            gse89632_expr,
            gse89632_meta,
            out_dir / "gse89632_SS_vs_NASH_limma.csv",
            condition_column="condition",
            group1=["SS"],
            group2=["NASH"],
            comparison="SS_vs_NASH",
        )
        fibrosis = differential_expression_limma(
            gse49541_expr,
            gse49541_meta,
            out_dir / "gse49541_mild_vs_advanced_limma.csv",
            condition_column="condition",
            group1=["mild_fibrosis"],
            group2=["advanced_fibrosis"],
            comparison="mild_vs_advanced_fibrosis",
        )
        tumor = differential_expression_limma(
            Path(paths["GSE164441_counts"]),
            gse164441_meta,
            out_dir / "gse164441_adjacent_vs_tumor_limma.csv",
            condition_column="condition",
            group1=["adjacent_normal"],
            group2=["tumor"],
            comparison="adjacent_normal_vs_tumor",
            input_type="counts",
            pair_column="patient",
        )
        gse135251_disease = differential_expression_limma(
            Path(paths["GSE135251_counts"]),
            gse135251_meta,
            out_dir / "gse135251_control_vs_NAFLD_limma.csv",
            condition_column="condition",
            group1=["control"],
            group2=["NAFLD"],
            comparison="control_vs_NAFLD",
            input_type="counts",
        )
        for name, frame in {
            "gse89632_HC_vs_NAFLD": hc_vs_nafld,
            "gse89632_SS_vs_NASH": ss_vs_nash,
            "gse49541_mild_vs_advanced": fibrosis,
            "gse164441_adjacent_vs_tumor": tumor,
            "gse135251_control_vs_NAFLD": gse135251_disease,
        }.items():
            differential_summary_table(frame).to_csv(
                out_dir / f"{name}_with_fdr.csv",
                index=False,
            )

        coexpression_config = dict(
            self.context.config.get("coexpression") or {}
        )
        if bool(coexpression_config.get("enabled", True)):
            try:
                coexpression = run_coexpression_analysis(
                    gse89632_expr,
                    gse89632_meta,
                    out_dir / "coexpression",
                    max_genes=int(coexpression_config.get("max_genes", 2000)),
                    min_module_size=int(
                        coexpression_config.get("min_module_size", 20)
                    ),
                    distance_threshold=float(
                        coexpression_config.get("distance_threshold", 0.5)
                    ),
                    kme_threshold=float(
                        coexpression_config.get("kme_threshold", 0.5)
                    ),
                    top_hubs_per_module=int(
                        coexpression_config.get("top_hubs_per_module", 25)
                    ),
                )
            except Exception as exc:  # noqa: BLE001
                LOG.warning("co-expression screening failed: %s", exc)
                coexpression = {
                    "status": "failed",
                    "reason": str(exc),
                }
                write_json(
                    out_dir
                    / "coexpression"
                    / "coexpression_summary.json",
                    coexpression,
                )
        else:
            coexpression = {
                "status": "skipped",
                "reason": "coexpression.enabled is false",
            }

        candidates = self._core_candidate_genes()
        ordered = [gene for gene in candidates if gene in set(hc_vs_nafld["gene"])]
        top_candidates = ordered[:30] or hc_vs_nafld["gene"].head(30).tolist()
        heatmap = candidate_heatmap(
            gse89632_expr,
            gse89632_meta,
            top_candidates,
            out_dir / "fig2f_candidate_gene_heatmap.png",
            group_order=["HC", "SS", "NASH"],
        )
        validation_boxplots(
            gse49541_expr,
            gse49541_meta,
            top_candidates[:6],
            out_dir,
            comparisons=("mild_fibrosis", "advanced_fibrosis"),
            prefix="fig2g_gse49541",
        )
        validation_boxplots(
            gse164441_expr,
            gse164441_meta,
            top_candidates[:6],
            out_dir,
            comparisons=("adjacent_normal", "tumor"),
            prefix="fig2h_gse164441",
            pair_column="patient",
        )
        write_json(
            out_dir / "bulk_summary.json",
            {
                "datasets": {
                    "GSE89632": {
                        "samples": int(pd.read_csv(gse89632_meta, index_col=0).shape[0]),
                        "genes": int(pd.read_csv(gse89632_expr, index_col=0).shape[0]),
                        "condition_counts": (
                            pd.read_csv(gse89632_meta, index_col=0)["condition"]
                            .astype(str)
                            .value_counts()
                            .to_dict()
                        ),
                        "endpoint": registry["GSE89632"]["endpoint"],
                    },
                    "GSE49541": {
                        "samples": int(pd.read_csv(gse49541_meta, index_col=0).shape[0]),
                        "genes": int(pd.read_csv(gse49541_expr, index_col=0).shape[0]),
                        "condition_counts": (
                            pd.read_csv(gse49541_meta, index_col=0)["condition"]
                            .astype(str)
                            .value_counts()
                            .to_dict()
                        ),
                        "endpoint": registry["GSE49541"]["endpoint"],
                        "validation_status": registry["GSE49541"][
                            "validation_status"
                        ],
                    },
                    "GSE164441": {
                        "samples": int(pd.read_csv(gse164441_meta, index_col=0).shape[0]),
                        "genes": int(pd.read_csv(gse164441_expr, index_col=0).shape[0]),
                        "condition_counts": (
                            pd.read_csv(gse164441_meta, index_col=0)["condition"]
                            .astype(str)
                            .value_counts()
                            .to_dict()
                        ),
                        "endpoint": registry["GSE164441"]["endpoint"],
                        "validation_status": registry["GSE164441"][
                            "validation_status"
                        ],
                        "paired_complete": int(
                            pd.read_csv(paths["GSE164441_pair_audit"])[
                                "paired_complete"
                            ].sum()
                        ),
                        "pair_audit": str(paths["GSE164441_pair_audit"]),
                    },
                    "GSE135251": {
                        "samples": int(pd.read_csv(gse135251_meta, index_col=0).shape[0]),
                        "genes": int(pd.read_csv(gse135251_expr, index_col=0).shape[0]),
                        "condition_counts": (
                            pd.read_csv(gse135251_meta, index_col=0)["condition"]
                            .astype(str)
                            .value_counts()
                            .to_dict()
                        ),
                        "endpoint": registry["GSE135251"]["endpoint"],
                        "validation_status": registry["GSE135251"][
                            "validation_status"
                        ],
                        "expected_condition_counts": registry["GSE135251"].get(
                            "expected_condition_counts",
                            {},
                        ),
                        "expected_condition_counts_match": (
                            {
                                str(key): int(value)
                                for key, value in (
                                    pd.read_csv(gse135251_meta, index_col=0)[
                                        "condition"
                                    ]
                                    .astype(str)
                                    .value_counts()
                                    .to_dict()
                                ).items()
                            }
                            == {
                                str(key): int(value)
                                for key, value in (
                                    registry["GSE135251"].get(
                                        "expected_condition_counts",
                                        {},
                                    )
                                ).items()
                            }
                        ),
                    },
                },
                "candidate_genes": top_candidates,
                "gse164441_caveat": (
                    "GSE164441 is NAFLD-associated HCC tumor versus adjacent "
                    "non-tumor tissue, not healthy liver versus NAFLD. Its AUC "
                    "is reported as a different clinical endpoint."
                ),
                "coexpression": coexpression,
            },
        )
        return {
            "processed": {key: str(value) for key, value in paths.items()},
            "deg": {
                "hc_vs_nafld": str(out_dir / "gse89632_HC_vs_NAFLD_limma.csv"),
                "ss_vs_nash": str(out_dir / "gse89632_SS_vs_NASH_limma.csv"),
                "fibrosis": str(out_dir / "gse49541_mild_vs_advanced_limma.csv"),
                "tumor": str(out_dir / "gse164441_adjacent_vs_tumor_limma.csv"),
            },
            "heatmap": _serializable(heatmap),
            "candidate_genes": top_candidates,
            "coexpression": coexpression,
        }

    def stage_ml(self) -> dict[str, Any]:
        out_dir = self.context.dir("05_machine_learning")
        paths = prepare_bulk_data(
            self.context.raw_dir,
            self.context.processed_dir,
            force=False,
        )
        ppi = pd.read_csv(
            self.context.output_root / "03_intersection_ppi" / "ppi_hub_metrics.csv"
        )
        overlap = pd.read_csv(
            self.context.output_root / "03_intersection_ppi" / "compound_disease_overlap.csv"
        )
        disease = pd.read_csv(
            self.context.output_root / "02_disease_targets" / "disease_targets.csv"
        )
        candidate_sets = {
            "compound_disease_intersection": overlap["gene"].astype(str).tolist(),
            "ppi_top20": ppi.head(20)["gene"].astype(str).tolist(),
            "ppi_plus_disease_top30": list(
                dict.fromkeys(
                    ppi.head(15)["gene"].astype(str).tolist()
                    + disease.head(15)["gene"].astype(str).tolist()
                )
            ),
        }
        # Full-cohort, label-associated modules are exploratory only. They must
        # not enter CV as preselected static features.
        validation = {
            "GSE49541_fibrosis": {
                "expression": str(paths["GSE49541_expression"]),
                "metadata": str(paths["GSE49541_metadata"]),
                "condition_map": {"mild_fibrosis": 0, "advanced_fibrosis": 1},
                "comparison": ("advanced fibrosis", "mild fibrosis"),
                "title": "GSE49541: advanced vs mild fibrosis",
                "evaluate_auc_target": False,
                "endpoint_role": "different_endpoint_fibrosis",
                "endpoint_warning": (
                    "Fibrosis stage is not the same endpoint as healthy "
                    "versus NAFLD; do not interpret this AUC as the primary "
                    "NAFLD endpoint."
                ),
            },
            "GSE164441_tumor": {
                "expression": str(paths["GSE164441_expression"]),
                "metadata": str(paths["GSE164441_metadata"]),
                "condition_map": {"adjacent_normal": 0, "tumor": 1},
                "comparison": ("tumor", "adjacent normal"),
                "title": "GSE164441: tumor vs adjacent non-tumor",
                "evaluate_auc_target": False,
                "endpoint_role": "different_endpoint_hcc",
                "endpoint_warning": (
                    "Tumor versus adjacent normal is an HCC endpoint, not a "
                    "healthy versus NAFLD endpoint."
                ),
            },
            "GSE135251_NAFLD": {
                "expression": str(paths["GSE135251_expression"]),
                "metadata": str(paths["GSE135251_metadata"]),
                "condition_map": {"control": 0, "NAFLD": 1},
                "comparison": ("NAFLD", "control"),
                "title": "GSE135251: NAFLD vs healthy control",
                "evaluate_auc_target": False,
                "endpoint_role": "supplementary_true_nafld",
                "endpoint_warning": (
                    "The cohort is the closest available same-endpoint "
                    "candidate, but repository history records prior model "
                    "or preprocessing exposure; call it reused validation "
                    "evidence rather than a fully unseen validation cohort."
                ),
            },
        }
        result = run_ml_validation(
            Path(paths["GSE89632_expression"]),
            Path(paths["GSE89632_metadata"]),
            candidate_sets,
            validation,
            out_dir,
            seed=int(self.context.config["ml"]["seed"]),
            cv_folds=int(self.context.config["ml"]["cv_folds"]),
            cv_repeats=int(self.context.config["ml"].get("cv_repeats", 5)),
        )
        aliases = {
            "GSE49541_fibrosis_roc.png": "fig3c_gse49541_cross_endpoint_roc.png",
            "GSE135251_NAFLD_roc.png": "fig3d_gse135251_same_endpoint_roc.png",
            "GSE164441_tumor_roc.png": "fig3d_supplementary_gse164441_hcc_roc.png",
        }
        for source_name, target_name in aliases.items():
            source = out_dir / source_name
            if source.exists():
                shutil.copy2(source, out_dir / target_name)
                for suffix in (".pdf", ".svg"):
                    vector = source.with_suffix(suffix)
                    if vector.exists():
                        shutil.copy2(vector, (out_dir / target_name).with_suffix(suffix))
        core_genes = list(result["core_genes"])[:3]
        if core_genes:
            validation_boxplots(
                Path(paths["GSE49541_expression"]),
                Path(paths["GSE49541_metadata"]),
                core_genes,
                out_dir,
                comparisons=("mild_fibrosis", "advanced_fibrosis"),
                prefix="fig3h_gse49541_core_genes",
            )
        else:
            write_json(
                out_dir / "fig3h_status.json",
                {
                    "status": "missing",
                    "reason": (
                        "SHAP did not produce a valid model explanation; no "
                        "surrogate importance was renamed as SHAP."
                    ),
                },
            )
        write_json(
            out_dir / "ml_core_genes.json",
            {
                "core_genes": core_genes,
                "source": "SHAP ranking of explicitly identified uncalibrated explanatory model",
            },
        )
        return {**_serializable(result), "core_genes": core_genes}

    def stage_mouse(self) -> dict[str, Any]:
        out_dir = self.context.dir("06_single_cell_mouse")
        extracted = self.context.raw_dir / "extracted" / "GSE270583"
        soft = self.context.raw_dir / "GSE270583_family.soft.gz"
        core_genes = self._ml_core_genes()
        result = run_mouse_single_cell(
            extracted,
            soft,
            core_genes,
            out_dir,
            cellchat_backend=str(
                self.context.config["single_cell"]["mouse"].get(
                    "cellchat_backend",
                    "r_cellchat",
                )
            ),
            cellchat_max_cells=int(
                self.context.config["single_cell"]["mouse"].get(
                    "cellchat_max_cells",
                    self.context.config["single_cell"]["mouse"].get(
                        "max_cells",
                        5000,
                    ),
                )
            ),
            cellchat_permutations=int(
                self.context.config["single_cell"]["mouse"].get(
                    "cellchat_permutations",
                    100,
                )
            ),
            cellchat_seed=int(
                self.context.config["single_cell"]["mouse"].get(
                    "cellchat_seed",
                    123,
                )
            ),
            allow_network=bool(
                (self.context.config.get("evidence") or {}).get(
                    "allow_network",
                    True,
                )
            ),
        )
        append_cohort_rows(
            self.context.output_root,
            result.get("cohort_rows") or [],
        )
        update_data_checksums(
            self.context.output_root,
            [Path(result["h5ad"])],
        )
        knockout = self._run_insilico_knockout(out_dir, core_genes)
        return {
            "h5ad": str(result["h5ad"]),
            "core_genes": result["core_genes"],
            "cellchat": {
                "interactions_csv": str(out_dir / "cellchat_interactions.csv"),
                "pathways_csv": str(out_dir / "cellchat_pathways.csv"),
                "n_interactions": len(result["cellchat"]["interactions"]),
                "n_pathways": len(result["cellchat"]["pathways"]),
            },
            "knockout": knockout,
        }

    def stage_human(self) -> dict[str, Any]:
        out_dir = self.context.dir("07_single_cell_human")
        extracted = self.context.raw_dir / "extracted" / "GSE202379"
        soft = self.context.raw_dir / "GSE202379_family.soft.gz"
        core_genes = self._ml_core_genes()
        result = run_human_single_cell(
            extracted,
            soft,
            core_genes,
            out_dir,
            max_cells_per_sample=int(
                self.context.config["single_cell"]["human"]["max_cells_per_sample"]
            ),
            seed=int(self.context.config["single_cell"]["human"]["seed"]),
        )
        append_cohort_rows(
            self.context.output_root,
            result.get("cohort_rows") or [],
        )
        update_data_checksums(
            self.context.output_root,
            [Path(result["h5ad"])],
        )
        return {
            "h5ad": str(result["h5ad"]),
            "core_genes": result["core_genes"],
            "expression": str(result["expression"]),
        }

    def stage_docking(self) -> dict[str, Any]:
        out_dir = self.context.dir("08_docking")
        ppi = pd.read_csv(
            self.context.output_root / "03_intersection_ppi" / "ppi_hub_metrics.csv"
        )
        core = self._ml_core_genes()
        ordered = [
            gene
            for gene in dict.fromkeys(core + ppi["gene"].astype(str).tolist())
            if not re.match(r"^(MIR|LET|LINC|SNOR|SCARNA|RMRP)", gene, flags=re.IGNORECASE)
        ]
        limit = int(self.context.config["docking"]["targets"])
        selected = ordered[:limit]
        target_rows = []
        for gene in selected:
            try:
                from .docking_md import uniprot_for_gene

                target_rows.append(uniprot_for_gene(gene))
            except Exception as exc:  # noqa: BLE001
                LOG.warning("UniProt mapping failed for docking target %s: %s", gene, exc)
                target_rows.append({"gene": gene})
        compound_path = (
            self.context.output_root
            / "01_compound_characterization"
            / "compound_properties.csv"
        )
        compound = pd.read_csv(compound_path).iloc[0].to_dict()
        result = run_docking_for_targets(
            target_rows,
            compound,
            out_dir,
            exhaustiveness=int(self.context.config["docking"]["exhaustiveness"]),
            num_modes=int(self.context.config["docking"].get("num_modes", 20)),
            replicates=int(self.context.config["docking"].get("replicates", 3)),
            cpu=int(self.context.config["docking"]["cpu"]),
            timeout=int(self.context.config["docking"]["timeout_seconds"]),
        )
        update_data_checksums(
            self.context.output_root,
            [
                out_dir / "docking_scores.csv",
                out_dir / "docking_run_manifest.json",
                out_dir / "fig5b_status.json",
            ],
        )
        return _serializable(result)

    def stage_md(self) -> dict[str, Any]:
        out_dir = self.context.dir("09_md_mmpbsa")
        docking_dir = self.context.output_root / "08_docking"
        result = prepare_or_run_md(
            docking_dir,
            run=bool(self.context.config["md"]["run"]),
            gpu=bool(self.context.config["md"]["gpu"]),
            cpu=int(self.context.config["md"]["cpu"]),
            protein_forcefield=str(
                self.context.config["md"].get(
                    "protein_forcefield",
                    "amber14sb",
                )
            ),
            timeout=int(self.context.config["md"]["timeout_seconds"]),
        )
        write_json(out_dir / "md_stage_summary.json", result)
        nested = Path(str(result.get("run_dir") or ""))
        if str(nested) and nested.is_dir():
            for path in nested.rglob("*"):
                if path.is_file() and path.suffix.lower() in {
                    ".mdp",
                    ".csv",
                    ".json",
                    ".html",
                    ".pdb",
                }:
                    relative = path.relative_to(nested)
                    destination = out_dir / "prepared_inputs" / relative
                    ensure_dir(destination.parent)
                    shutil.copy2(path, destination)
        md_figures = generate_plan_md_figures(docking_dir, out_dir)
        update_data_checksums(
            self.context.output_root,
            [
                docking_dir / "md_run_manifest.json",
                out_dir / "md_plan_figures.json",
            ],
        )
        _plot_md_status(
            result,
            out_dir / "fig5d_md_preparation_status.png",
        )
        if not self.context.config["md"]["run"]:
            (out_dir / "SIMULATION_STATUS.txt").write_text(
                "GROMACS inputs were prepared for 100 ns. The 100 ns production "
                "run was not started in this invocation; execute the generated "
                "run script or rerun the stage with md.run=true.\n",
                encoding="utf-8",
            )
            for panel in ("d_rmsd", "e_ligand_rmsd", "f_rmsf", "g_rg", "h_mmpbsa"):
                (out_dir / f"fig5{panel}_NOT_RUN.txt").write_text(
                    "The 100 ns production trajectory was not started. "
                    "Inputs and run configuration are prepared and no numerical "
                    "trajectory result is claimed.\n",
                    encoding="utf-8",
                )
                _plot_not_run_panel(
                    panel,
                    out_dir / f"fig5{panel}_NOT_RUN.png",
                )
        return _serializable({**result, "plan_figures": md_figures})

    def stage_report(self) -> dict[str, Any]:
        out_dir = self.context.dir("10_reports")
        inventory = self._inventory()
        inventory.to_csv(out_dir / "result_inventory.csv", index=False)
        figure_index = _figure_index(self.context.output_root)
        figure_index.to_csv(out_dir / "figure_index.csv", index=False)
        status = self._stage_statuses()
        coverage = read_json(
            out_dir / "plan_coverage" / "plan_coverage.json",
            {},
        )
        delivery = audit_delivery_readiness(
            self.context.output_root,
            coverage_summary=coverage,
        )
        governance_dir = ensure_dir(out_dir / "governance")
        write_json(governance_dir / "delivery_readiness.json", delivery)
        report = _render_report(self.context, inventory, status, figure_index)
        report_path = out_dir / "experiment_plan_one_report.html"
        report_path.write_text(report, encoding="utf-8")
        markdown = _render_markdown_report(
            self.context,
            inventory,
            status,
            figure_index,
        )
        markdown_path = out_dir / "experiment_plan_one_report.md"
        markdown_path.write_text(markdown, encoding="utf-8")
        summary_path = self.context.output_root / "RESULTS_SUMMARY.md"
        summary_path.write_text(
            _render_results_summary(self.context, status),
            encoding="utf-8",
        )
        write_json(
            out_dir / "analysis_status.json",
            {
                "version": __version__,
                "run_id": self.run_id,
                "output_root": str(self.context.output_root),
                "stages": status,
                "files": len(inventory),
                "plan_coverage": read_json(
                    out_dir / "plan_coverage" / "plan_coverage.json",
                    {},
                ),
                "delivery_readiness": delivery,
                "completion_target_percent": float(
                    self.context.config.get("completion_target_percent", 90.0)
                ),
                "projected_completion_target_met": bool(
                    float(
                        read_json(
                            out_dir
                            / "plan_coverage"
                            / "plan_coverage.json",
                            {},
                        ).get(
                            "projected_completion_percent_with_prerequisites",
                            0.0,
                        )
                        or 0.0
                    )
                    >= float(
                        self.context.config.get(
                            "completion_target_percent",
                            90.0,
                        )
                    )
                ),
                "notes": [
                    "GSE164441 is tumor versus adjacent non-tumor, not a healthy-versus-NAFLD cohort.",
                    "GSE135251 is a same-endpoint candidate but prior analysis exposure is recorded; it is not described as fully unseen.",
                    "Communication analysis uses R CellChat when configured; cell-level probabilities are never treated as independent biological replication.",
                    "Predicted perturbation results use R scTenifoldKnk in the experiment-plan configuration and remain network predictions, not wet-lab knockouts.",
                    "100 ns MD production is prepared but started only when md.run=true.",
                    "Publication-grade status remains false until a full rerun resolves all panels and causal claims have wet-lab evidence.",
                ],
            },
        )
        return {
            "report_html": str(report_path),
            "report_markdown": str(markdown_path),
            "inventory": str(out_dir / "result_inventory.csv"),
            "figure_index": str(out_dir / "figure_index.csv"),
            "summary": str(summary_path),
            "n_files": len(inventory),
        }

    def stage_classify(self) -> dict[str, Any]:
        result = classify_experiment_plan_results(self.context.output_root)
        return _serializable(result)

    def stage_figure_audit(self) -> dict[str, Any]:
        result = audit_figures(self.context.output_root)
        coverage = audit_plan_coverage(self.context.output_root)
        delivery = audit_delivery_readiness(
            self.context.output_root,
            coverage_summary=coverage,
        )
        governance_dir = ensure_dir(
            self.context.output_root / "10_reports" / "governance"
        )
        write_json(governance_dir / "delivery_readiness.json", delivery)
        pd.DataFrame(
            [
                {
                    "tier": tier,
                    **(
                        dict(payload)
                        if isinstance(payload, dict)
                        else {"status": str(payload)}
                    ),
                }
                for tier, payload in (delivery.get("tiers") or {}).items()
            ]
        ).to_csv(governance_dir / "delivery_readiness.csv", index=False)
        return _serializable(
            {
                "figure_audit": result,
                "plan_coverage": coverage,
                "delivery_readiness": delivery,
            }
        )

    def _core_candidate_genes(self) -> list[str]:
        path = (
            self.context.output_root
            / "03_intersection_ppi"
            / "compound_disease_overlap.csv"
        )
        if not path.exists():
            return []
        return pd.read_csv(path)["gene"].astype(str).tolist()

    def _ml_core_genes(self) -> list[str]:
        path = self.context.output_root / "05_machine_learning" / "ml_core_genes.json"
        payload = read_json(path, {}) if path.exists() else {}
        core = list(payload.get("core_genes") or [])
        core = [
            gene
            for gene in dict.fromkeys(core)
            if not re.match(r"^(MIR|LET|LINC|SNOR|SCARNA|RMRP)", gene, flags=re.IGNORECASE)
        ]
        return core[:3]

    def _run_insilico_knockout(
        self,
        mouse_dir: Path,
        core_genes: list[str],
    ) -> dict[str, Any]:
        if not core_genes:
            return {"status": "skipped", "reason": "no core genes"}
        insilico_config = dict(
            (self.context.config.get("single_cell") or {}).get("insilico")
            or {}
        )
        engine = str(insilico_config.get("engine") or "local_grn")
        scientific_role = str(
            insilico_config.get("scientific_role")
            or "predicted_perturbation_response"
        )
        mouse_config = dict(
            (self.context.config.get("single_cell") or {}).get("mouse")
            or {}
        )
        insilico_max_cells = int(
            insilico_config.get(
                "max_cells",
                mouse_config.get("max_cells", 5000),
            )
        )
        insilico_max_genes = int(
            insilico_config.get(
                "max_genes",
                mouse_config.get("max_genes", 1800),
            )
        )
        try:
            import anndata as ad

            data = ad.read_h5ad(mouse_dir / "mouse_liver_processed.h5ad")
            knockout_dir = ensure_dir(mouse_dir / "virtual_knockout")
            gene = next(
                (value for value in core_genes if value in data.var_names),
                core_genes[0],
            )
            previous_config_path = knockout_dir / "insilico_config.json"
            previous_report = (
                knockout_dir
                / "outputs"
                / "run_001"
                / "results"
                / "04_knockout"
                / "in_silico"
                / "in_silico_knockout_report.html"
            )
            input_signature = hashlib.sha256(
                json.dumps(
                    {
                        "gene": str(gene).upper(),
                        "h5ad": self.context.fingerprint_file(
                            mouse_dir / "mouse_liver_processed.h5ad"
                        ),
                        "engine": engine,
                        "scientific_role": scientific_role,
                        "max_cells": insilico_max_cells,
                        "max_genes": insilico_max_genes,
                        "seed": 123,
                    },
                    sort_keys=True,
                ).encode("utf-8")
            ).hexdigest()
            if previous_config_path.exists() and previous_report.exists():
                previous = read_json(previous_config_path, {})
                previous_gene = (
                    (previous.get("insilico_knockout") or {}).get("ko_gene")
                )
                previous_signature = (
                    previous.get("insilico_knockout") or {}
                ).get("_input_signature")
                if (
                    str(previous_gene).upper() == str(gene).upper()
                    and previous_signature == input_signature
                ):
                    _export_knockout_outputs(knockout_dir, mouse_dir)
                    return {
                        "status": "completed",
                        "gene": gene,
                        "engine": f"cached {engine} output",
                        "scientific_role": scientific_role,
                        "report": str(previous_report),
                        "data": {
                            "go": str(mouse_dir / "virtual_knockout_go_enrichment.csv"),
                            "kegg": str(mouse_dir / "virtual_knockout_kegg_enrichment.csv"),
                            "target_changes": str(mouse_dir / "virtual_knockout_target_changes.csv"),
                        },
                    }
            selected_cells = _balanced_cell_sample(
                data,
                max_cells=insilico_max_cells,
            )
            data = data[selected_cells].copy()
            counts = data.layers.get("counts")
            if counts is None:
                raise RuntimeError("mouse AnnData has no raw count layer")
            expr = pd.DataFrame.sparse.from_spmatrix(
                counts.T,
                index=data.var_names.astype(str),
                columns=data.obs_names.astype(str),
            )
            expr.to_csv(
                knockout_dir / "expression.csv.gz",
                index=True,
                compression="gzip",
            )
            metadata = pd.DataFrame(
                {
                    "cell": data.obs_names.astype(str),
                    "cell_type": data.obs["cell_type"].astype(str).to_numpy(),
                    "condition": data.obs["condition"].astype(str).to_numpy(),
                }
            )
            metadata.to_csv(knockout_dir / "metadata.csv", index=False)
            embedding = pd.DataFrame(
                data.obsm["X_umap"],
                columns=["umap_1", "umap_2"],
            )
            embedding.insert(0, "cell", data.obs_names.astype(str))
            embedding.to_csv(knockout_dir / "embedding.csv", index=False)
            config = {
                "workdir": str(knockout_dir),
                "output_dir": "outputs/run_001",
                "receptor": {},
                "ligand": {},
                "docking": {},
                "analysis": {},
                "knockout": {
                    "expression_csv": str(knockout_dir / "expression.csv.gz"),
                    "metadata_csv": str(knockout_dir / "metadata.csv"),
                    "cell_type_column": "cell_type",
                    "group_column": "condition",
                    "case_label": "HFD",
                    "normal_label": "NCD",
                },
                "insilico_knockout": {
                    "enabled": True,
                    "ko_gene": gene,
                    "engine": engine,
                    "scientific_role": scientific_role,
                    "raw_count_input": True,
                    "species": "mm",
                    "embedding_csv": str(knockout_dir / "embedding.csv"),
                    "min_cells": 100,
                    "max_cells": insilico_max_cells,
                    "max_genes": insilico_max_genes,
                    "seed": 123,
                    "knn_impute_neighbors": 66,
                    "n_propagation": 3,
                    "network_edges_per_regulator": 300,
                    "target_top_n": 15,
                    "bar_top_n": 10,
                    "run_enrichment": True,
                    "enrichment_genes": 150,
                    "enrichment_timeout": 900,
                    "figures": True,
                    "_input_signature": input_signature,
                    "_provenance_output_dir": str(knockout_dir),
                },
            }
            config_path = knockout_dir / "insilico_config.json"
            write_json(config_path, config)
            from docking.config import ResolvedConfig
            from docking.insilico import run_insilico_knockout
            from docking.utils import setup_logging

            cfg = ResolvedConfig(config, config_path)
            result = run_insilico_knockout(cfg, setup_logging(str(knockout_dir / "knockout.log")), gene)
            _export_knockout_outputs(knockout_dir, mouse_dir)
            return {"status": "completed", "gene": gene, "result": _serializable(result)}
        except Exception as exc:  # noqa: BLE001
            LOG.warning("virtual knockout failed: %s", exc)
            reason = str(exc)
            return {
                "status": (
                    "blocked"
                    if "blocked" in reason.lower()
                    else "failed"
                ),
                "reason": reason,
            }

    def _inventory(self) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        for path in sorted(self.context.output_root.rglob("*")):
            if not path.is_file() or path.parts[-2:-1] == (".stages",):
                continue
            relative = path.relative_to(self.context.output_root)
            rows.append(
                {
                    "category": relative.parts[0] if relative.parts else "",
                    "file": str(relative),
                    "absolute_path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "modified": time.strftime(
                        "%Y-%m-%d %H:%M:%S",
                        time.localtime(path.stat().st_mtime),
                    ),
                }
            )
        return pd.DataFrame(rows)

    def _stage_statuses(self) -> dict[str, Any]:
        status = {}
        for stage in STAGES:
            path = self.context.state_dir / f"{stage}.json"
            status[stage] = read_json(path, {"status": "not_run"})
        return status

    def _write_manifest(self) -> None:
        out_dir = self.context.dir("10_reports")
        software = {
            "python": sys.version,
            "platform": platform.platform(),
            "packages": _package_versions(
                [
                    "numpy",
                    "pandas",
                    "scipy",
                    "scikit-learn",
                    "scanpy",
                    "anndata",
                    "matplotlib",
                    "rdkit",
                    "plip",
                    "scTenifoldpy",
                    "openbabel",
                    "xgboost",
                    "lightgbm",
                    "shap",
                ]
            ),
        }
        write_json(
            out_dir / "analysis_manifest.json",
            {
                "version": __version__,
                "run_id": self.run_id,
                "script_freeze": read_json(
                    self.context.root
                    / "config"
                    / "experiment_plan_one_freeze.json",
                    {},
                ),
                "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                "output_root": str(self.context.output_root),
                "config": self.context.config,
                "software": software,
                "external_tools": audit_plan_environment(),
                "method_decisions": load_governance()[
                    "method_changes"
                ].to_dict("records"),
                "reproducibility": {
                    "git_commit": _git_commit(self.context.root),
                    "random_seeds": {
                        "ml": self.context.config.get("ml", {}).get("seed"),
                        "single_cell_mouse": self.context.config.get("single_cell", {}).get("mouse", {}).get("cellchat_seed"),
                        "single_cell_human": self.context.config.get("single_cell", {}).get("human", {}).get("seed"),
                    },
                    "dataset_sources": {
                        name: {"role": spec.get("role"), "platform": spec.get("platform"), "species": spec.get("species")}
                        for name, spec in (self.context.config.get("datasets") or {}).items()
                    },
                    "known_reference_targets": (
                        self.context.config.get("known_reference_targets")
                        or self.context.config.get("novelty_exclusions", [])
                    ),
                    "known_reference_policy": (
                        "retained as positive references; not removed from "
                        "candidate, model, or docking pools"
                    ),
                },
                "stages": {
                    name: read_json(self.context.state_dir / f"{name}.json", {"status": "not_run"})
                    for name in STAGES
                },
            },
        )

    def _write_delivery_readiness(self) -> None:
        coverage = read_json(
            self.context.output_root
            / "10_reports"
            / "plan_coverage"
            / "plan_coverage.json",
            {},
        )
        delivery = audit_delivery_readiness(
            self.context.output_root,
            coverage_summary=coverage,
        )
        governance_dir = ensure_dir(
            self.context.output_root / "10_reports" / "governance"
        )
        write_json(governance_dir / "delivery_readiness.json", delivery)
        pd.DataFrame(
            [
                {
                    "tier": tier,
                    **(
                        dict(payload)
                        if isinstance(payload, dict)
                        else {"status": str(payload)}
                    ),
                }
                for tier, payload in (delivery.get("tiers") or {}).items()
            ]
        ).to_csv(governance_dir / "delivery_readiness.csv", index=False)


def _balanced_cell_sample(data: Any, max_cells: int, seed: int = 123) -> np.ndarray:
    import anndata as ad  # noqa: F401

    rng = np.random.default_rng(seed)
    groups = data.obs.groupby(["condition", "cell_type"], observed=True).indices
    selected: list[int] = []
    per_group = max(10, max_cells // max(len(groups), 1))
    for indices in groups.values():
        indices = np.asarray(indices)
        if len(indices) > per_group:
            indices = rng.choice(indices, size=per_group, replace=False)
        selected.extend(int(value) for value in indices)
    if len(selected) > max_cells:
        selected = rng.choice(selected, size=max_cells, replace=False).tolist()
    return np.asarray(sorted(set(int(value) for value in selected)), dtype=int)


def _plot_source_coverage(statuses: dict[str, Any], output: Path) -> None:
    names = list(statuses)
    fig, ax = plt.subplots(figsize=(6.6, max(3.2, 0.48 * len(names) + 0.7)))
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, len(names) - 0.4)
    ax.axis("off")
    for row_index, name in enumerate(reversed(names)):
        status = str(statuses[name].get("status") or "unknown")
        count = statuses[name].get("count")
        color = STATUS_COLORS.get(status, "#7A7A7A")
        y = row_index
        ax.scatter([0.05], [y], s=52, color=color, zorder=3)
        ax.text(
            0.11,
            y,
            str(name),
            ha="left",
            va="center",
            fontsize=8.2,
        )
        status_label = status.replace("_", " ")
        if count is not None:
            status_label += f" ({int(count)} targets)"
        ax.text(
            0.98,
            y,
            status_label,
            ha="right",
            va="center",
            fontsize=7.8,
            color="#3f4b55",
        )
        ax.plot([0.03, 0.97], [y - 0.36, y - 0.36], color="#e1e6ea", lw=0.6)
    ax.set_title("Disease target source availability", fontweight="bold")
    save_figure(fig, output)


def _plot_target_source_counts(sources: dict[str, set[str]], output: Path) -> None:
    pairs = sorted(
        ((name, len(genes)) for name, genes in sources.items()),
        key=lambda item: (item[1], item[0]),
    )
    names = [pair[0] for pair in pairs]
    counts = [pair[1] for pair in pairs]
    if not names:
        names = ["No source loaded"]
        counts = [0]
    fig, ax = plt.subplots(figsize=(6.2, max(3.0, 0.55 * len(names) + 1.0)))
    colors = [
        FIGURE_PALETTE[index % len(FIGURE_PALETTE)]
        for index in range(len(names))
    ]
    bars = ax.barh(names, counts, color=colors, edgecolor="white", linewidth=0.6)
    for bar, count in zip(bars, counts):
        ax.text(
            float(bar.get_width()) + max(counts, default=0) * 0.02,
            bar.get_y() + bar.get_height() / 2,
            str(count),
            ha="left",
            va="center",
            fontsize=7.5,
            fontweight="bold",
        )
    ax.set_xlim(0, max(max(counts, default=0) * 1.2, 1))
    ax.set_xlabel("Unique target genes")
    ax.set_title("Compound-target source coverage", fontweight="bold")
    ax.text(
        0,
        -0.24,
        "Only mapped human target genes are counted; unavailable sources are "
        "listed separately in source-status files.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
        color="#4a5560",
    )
    save_figure(fig, output)


def _export_knockout_outputs(knockout_root: Path, mouse_dir: Path) -> None:
    """Expose stable Figure 4g/h filenames from the existing knockout wrapper."""
    nested = (
        knockout_root
        / "outputs"
        / "run_001"
        / "results"
        / "04_knockout"
        / "in_silico"
    )
    figures = nested / "figures"
    data = nested / "data"
    mappings = {
        figures / "fig_63_ko_target_expression_table.png": mouse_dir
        / "fig4g_virtual_knockout_top10.png",
        figures / "fig_65_ko_regulatory_network.png": mouse_dir
        / "fig4g_virtual_knockout_network.png",
        figures / "fig_66_ko_shift_umap.png": mouse_dir
        / "fig4g_virtual_knockout_shift.png",
        data / "fig_63_ko_target_top15.csv": mouse_dir
        / "virtual_knockout_top15.csv",
        data / "insilico_target_changes.csv": mouse_dir
        / "virtual_knockout_target_changes.csv",
        data / "insilico_go_enrichment.csv": mouse_dir
        / "virtual_knockout_go_enrichment.csv",
        data / "insilico_kegg_enrichment.csv": mouse_dir
        / "virtual_knockout_kegg_enrichment.csv",
    }
    for source, destination in mappings.items():
        if source.exists():
            shutil.copy2(source, destination)
    top15_path = mouse_dir / "virtual_knockout_top15.csv"
    if top15_path.exists():
        top15 = pd.read_csv(top15_path)
        required = {"gene", "wt_mean", "ko_mean"}
        if required.issubset(top15.columns):
            top15 = top15.copy()
            top15["_effect"] = (
                pd.to_numeric(top15["delta"], errors="coerce").abs()
                if "delta" in top15.columns
                else (
                    pd.to_numeric(top15["wt_mean"], errors="coerce")
                    - pd.to_numeric(top15["ko_mean"], errors="coerce")
                ).abs()
            )
            top15 = top15.sort_values("_effect", ascending=False).head(10)
            genes = top15["gene"].astype(str).tolist()
            wt = pd.to_numeric(top15["wt_mean"], errors="coerce").to_numpy()
            ko = pd.to_numeric(top15["ko_mean"], errors="coerce").to_numpy()
            fig, ax = plt.subplots(figsize=(7.2, 4.8))
            y = np.arange(len(genes))
            height = 0.36
            ax.barh(y - height / 2, wt, height=height, color="#466b8a", label="WT")
            ax.barh(
                y + height / 2,
                ko,
                height=height,
                color="#c05b4d",
                label="Predicted perturbation",
            )
            ax.set_yticks(y)
            ax.set_yticklabels(genes)
            ax.invert_yaxis()
            ax.set_xlabel("Mean expression")
            ax.set_title(
                "Top predicted perturbation-response changes",
                fontweight="bold",
            )
            ax.legend(frameon=False)
            ax.grid(axis="x", alpha=0.2)
            save_figure(fig, mouse_dir / "fig4g_virtual_knockout_top10.png")
    go = data / "insilico_go_enrichment.csv"
    kegg = data / "insilico_kegg_enrichment.csv"
    go_text = go.read_text(encoding="utf-8", errors="replace") if go.exists() else ""
    kegg_text = kegg.read_text(encoding="utf-8", errors="replace") if kegg.exists() else ""
    if "no significant" in go_text.lower() or "no significant" in kegg_text.lower():
        fig, ax = plt.subplots(figsize=(7.2, 3.6))
        ax.axis("off")
        ax.text(
            0.5,
            0.62,
            "No significant GO/KEGG terms",
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )
        ax.text(
            0.5,
            0.36,
            "The predicted perturbation-response target set did not reach the configured "
            "FDR threshold. This is retained as an explicit negative result.",
            ha="center",
            va="center",
            fontsize=9,
        )
        save_figure(fig, mouse_dir / "fig4h_virtual_knockout_enrichment.png")


def _plot_md_status(result: dict[str, Any], output: Path) -> None:
    status = str(result.get("status") or "unknown")
    mode = str(result.get("mode") or "prepare")
    target = str(result.get("target") or "")
    simulation_ns = result.get("simulation_ns", 100.0)
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.axis("off")
    ax.text(
        0.5,
        0.78,
        "GROMACS MD status",
        ha="center",
        va="center",
        fontsize=15,
        fontweight="bold",
    )
    rows = [
        ("Target", target or "not available"),
        ("Mode", mode),
        ("Configured production", f"{simulation_ns:g} ns"),
        ("Temperature", f"{result.get('temperature_k', 310.0):g} K"),
        ("Pressure", f"{result.get('pressure_bar', 1.0):g} bar"),
        ("GPU requested", "yes" if result.get("gpu") else "no"),
        ("Status", status),
    ]
    y = 0.61
    for label, value in rows:
        ax.text(0.08, y, label, ha="left", va="center", fontsize=10)
        ax.text(0.92, y, value, ha="right", va="center", fontsize=10, fontweight="bold")
        y -= 0.09
    ax.text(
        0.5,
        0.03,
        "A status of 'prepared' means no trajectory metric was fabricated.",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#5d6670",
    )
    save_figure(fig, output)


def _plot_not_run_panel(panel: str, output: Path) -> None:
    labels = {
        "d_rmsd": "Protein backbone RMSD",
        "e_ligand_rmsd": "Ligand RMSD",
        "f_rmsf": "Residue RMSF",
        "g_rg": "Radius of gyration",
        "h_mmpbsa": "MM-PBSA energy decomposition",
    }
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.axis("off")
    ax.text(
        0.5,
        0.63,
        labels.get(panel, panel),
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.41,
        "Not run / not available",
        ha="center",
        va="center",
        fontsize=12,
        color="#b24c3c",
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.22,
        "The 100 ns GROMACS production run was not started. "
        "No numerical trajectory metric is plotted.",
        ha="center",
        va="center",
        fontsize=8.5,
        color="#5d6670",
    )
    save_figure(fig, output)


def _package_versions(packages: list[str]) -> dict[str, str]:
    from importlib import metadata

    output: dict[str, str] = {}
    for package in packages:
        try:
            output[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            output[package] = "not_installed"
    return output


def _git_commit(root: Path) -> str | None:
    """Return the source revision when the workflow is run from a checkout."""
    try:
        import subprocess

        result = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    commit = result.stdout.strip()
    return commit or None


def _serializable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, pd.DataFrame):
        return value.to_dict(orient="records")
    if isinstance(value, dict):
        return {str(key): _serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_serializable(item) for item in value]
    return value


def _render_report(
    context: PipelineContext,
    inventory: pd.DataFrame,
    status: dict[str, Any],
    figure_index: pd.DataFrame,
) -> str:
    rows = "".join(
        f"<tr><td>{html.escape(str(row.category))}</td>"
        f"<td>{html.escape(str(row.file))}</td>"
        f"<td>{int(row.size_bytes):,}</td></tr>"
        for row in inventory.itertuples(index=False)
    )
    stage_rows = "".join(
        f"<tr><td>{html.escape(name)}</td><td>{html.escape(str(value.get('status')))}</td>"
        f"<td>{value.get('elapsed_seconds', '')}</td></tr>"
        for name, value in status.items()
    )
    figure_rows = "".join(
        f"<tr><td>{html.escape(str(row.figure))}</td>"
        f"<td>{html.escape(str(row.panel))}</td>"
        f"<td>{html.escape(str(row.purpose))}</td>"
        f"<td>{html.escape(str(row.status))}</td>"
        f"<td>{html.escape(str(row.file))}</td></tr>"
        for row in figure_index.itertuples(index=False)
    )
    return f"""<!doctype html>
<html lang="zh">
<head>
<meta charset="utf-8">
<title>6PPD-Q / NAFLD experiment plan one results</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 36px; color: #25313d; }}
h1 {{ font-size: 25px; }}
table {{ border-collapse: collapse; width: 100%; margin: 12px 0 28px; }}
th, td {{ border-bottom: 1px solid #d8dee4; padding: 7px 8px; text-align: left; }}
th {{ background: #edf2f6; }}
code {{ background: #f1f4f6; padding: 2px 4px; }}
.note {{ background: #fff8e6; border-left: 4px solid #d99b34; padding: 12px; }}
</style>
</head>
<body>
<h1>6PPD-Q / NAFLD: experiment plan one results</h1>
<p>Output root: <code>{html.escape(str(context.output_root))}</code></p>
<div class="note">
<b>Method notes.</b> GSE164441 in the source plan is NAFLD-associated HCC
versus adjacent non-tumor tissue, not a healthy-versus-NAFLD cohort.
SwissTargetPrediction is queried directly. ChEMBL is represented by an
explicit similarity-based target screen when no exact ChEMBL molecule exists.
SEA is queried through its public API when enabled; cached or local target
tables are accepted when a service is unavailable. A WGCNA-style local
co-expression screen is supplementary to the five-plan figures and is not
presented as the original R WGCNA implementation.
GeneCards, OMIM and TTD bulk exports require licensed access unless local
files are supplied. The CellChat panel runs the original R implementation when
the configured backend is <code>r_cellchat</code>; cell-level probabilities
remain descriptive without biological-unit replication. The 100 ns MD
production run is prepared but not started unless <code>md.run=true</code>.
</div>
<h2>Stages</h2>
<table><thead><tr><th>Stage</th><th>Status</th><th>Elapsed seconds</th></tr></thead>
<tbody>{stage_rows}</tbody></table>
<h2>Figure index</h2>
<table><thead><tr><th>Figure</th><th>Panel</th><th>Purpose</th><th>Status</th><th>File</th></tr></thead>
<tbody>{figure_rows}</tbody></table>
<h2>Result inventory</h2>
<table><thead><tr><th>Category</th><th>File</th><th>Bytes</th></tr></thead>
<tbody>{rows}</tbody></table>
</body>
</html>
"""


def _render_markdown_report(
    context: PipelineContext,
    inventory: pd.DataFrame,
    status: dict[str, Any],
    figure_index: pd.DataFrame,
) -> str:
    lines = [
        "# 6PPD-Q / NAFLD experiment plan one",
        "",
        f"- Output root: `{context.output_root}`",
        f"- Result files: {len(inventory)}",
        "",
        "## Method notes",
        "",
        "- The source plan labels GSE164441 as a NAFLD validation set; GEO "
        "describes it as NAFLD-associated HCC tumor versus adjacent non-tumor "
        "tissue. Its model AUC is reported for that different endpoint.",
        "- ChEMBL has no exact molecule record for every query; an explicit "
        "similarity-based ChEMBL target source is used and labelled.",
        "- SEA is queried when enabled, with cache and explicit failure states; "
        "local prediction tables can supply PharmMapper or other sources.",
        "- The supplementary co-expression module screen is a local "
        "WGCNA-style approximation, not the original R WGCNA implementation.",
        "- GeneCards, OMIM and TTD require licensed/credentialed bulk access "
        "unless local source tables are supplied.",
        "- The communication panel runs the original R CellChat implementation "
        "when the configured backend is `r_cellchat`; cell-level probabilities "
        "are descriptive unless biological-unit replication is sufficient.",
        "- GROMACS inputs are prepared for 100 ns; production is not started "
        "unless `md.run` is true.",
        "",
        "## Stages",
        "",
        "| Stage | Status | Elapsed seconds |",
        "|---|---|---:|",
    ]
    for name, value in status.items():
        lines.append(
            f"| {name} | {value.get('status')} | "
            f"{value.get('elapsed_seconds', '')} |"
        )
    lines += [
        "",
        "## Figure index",
        "",
        "| Figure | Panel | Purpose | Status | File |",
        "|---|---|---|---|---|",
    ]
    for row in figure_index.itertuples(index=False):
        lines.append(
            f"| {row.figure} | {row.panel} | {row.purpose} | "
            f"{row.status} | `{row.file}` |"
        )
    lines += ["", "## Results", "", "| Category | File | Bytes |", "|---|---|---:|"]
    for row in inventory.itertuples(index=False):
        lines.append(f"| {row.category} | `{row.file}` | {int(row.size_bytes):,} |")
    return "\n".join(lines) + "\n"


def _render_results_summary(
    context: PipelineContext,
    status: dict[str, Any],
) -> str:
    lines = [
        "# 6PPD-Q / NAFLD analysis results",
        "",
        f"Output root: `{context.output_root}`",
        "",
        "## Key Results",
        "",
    ]
    ml_path = context.output_root / "05_machine_learning" / "ml_summary.json"
    if ml_path.exists():
        ml = read_json(ml_path, {})
        lines.extend(
            [
                f"- Best cross-validated classifier: {ml.get('best_model')} "
                f"({ml.get('best_feature_set')}), CV AUC "
                f"{float(ml.get('cv_auc', float('nan'))):.3f}; "
                f"{ml.get('cv_folds', 'NA')}-fold x "
                f"{ml.get('cv_repeats', 'NA')} repeats.",
                f"- SHAP core genes: {', '.join(ml.get('core_genes') or [])}.",
                "",
                "| External cohort | Endpoint | n | AUC | Endpoint target evaluation |",
                "|---|---|---:|---:|---|",
            ]
        )
        for row in ml.get("external_validation") or []:
            if row.get("evaluate_auc_target") is False:
                evaluation = "not evaluated (different endpoint)"
            elif row.get("target_met") is True:
                evaluation = "met"
            else:
                evaluation = "not met"
            lines.append(
                f"| {row.get('dataset')} | {row.get('comparison')} | "
                f"{row.get('n')} | {float(row.get('auc', float('nan'))):.3f} | "
                f"{evaluation} |"
            )
    mouse_path = (
        context.output_root / "06_single_cell_mouse" / "mouse_single_cell_summary.json"
    )
    human_path = (
        context.output_root / "07_single_cell_human" / "human_single_cell_summary.json"
    )
    if mouse_path.exists():
        mouse = read_json(mouse_path, {})
        lines += [
            "",
            f"- Mouse scRNA-seq: {mouse.get('n_cells', 0):,} cells, "
            f"{mouse.get('n_genes', 0):,} genes, "
            f"{len(mouse.get('cell_types') or {})} major cell types.",
        ]
    if human_path.exists():
        human = read_json(human_path, {})
        lines += [
            f"- Human snRNA-seq: {human.get('n_cells', 0):,} nuclei, "
            f"{human.get('samples', 0)} GEO samples.",
        ]
    coexpression_path = (
        context.output_root
        / "04_bulk_training"
        / "coexpression"
        / "coexpression_summary.json"
    )
    if coexpression_path.exists():
        coexpression = read_json(coexpression_path, {})
        if coexpression.get("status") == "completed":
            lines += [
                "",
                (
                    "- WGCNA-style co-expression screen: "
                    f"{coexpression.get('n_modules', 0)} modules, "
                    f"{coexpression.get('n_hubs', 0)} disease-module hubs; "
                    f"soft power {coexpression.get('selected_soft_power', 'NA')}."
                ),
                (
                    "- Co-expression method note: local signed-correlation "
                    "module screening, not the original R WGCNA TOM/dynamic "
                    "tree-cut implementation."
                ),
            ]
    docking_path = context.output_root / "08_docking" / "docking_scores.csv"
    if docking_path.exists():
        docking = pd.read_csv(docking_path).sort_values("best_affinity_kcal_mol")
        lines += [
            "",
            "| Target | Structure | Vina affinity (kcal/mol) |",
            "|---|---|---:|",
        ]
        for row in docking.itertuples(index=False):
            pdb_id = getattr(row, "pdb_id", "")
            structure_label = (
                f"{row.structure_source} {pdb_id}"
                if isinstance(pdb_id, str) and pdb_id and pdb_id.lower() != "nan"
                else str(row.structure_source)
            )
            lines.append(
                f"| {row.gene} | {structure_label} | "
                f"{float(row.best_affinity_kcal_mol):.3f} |"
            )
    audit_path = (
        context.output_root
        / "10_reports"
        / "figure_quality_audit"
        / "figure_quality_audit.json"
    )
    if audit_path.exists():
        audit = read_json(audit_path, {})
        overlap_counts = audit.get("label_overlap_counts") or {}
        lines += [
            "",
            "## Figure Audit",
            "",
            f"- Overall verdict: **{audit.get('overall_verdict', 'unknown')}**",
            (
                f"- Available panels: {audit.get('available_panels', 0)}/"
                f"{audit.get('panel_entries', 42)}"
            ),
            f"- Median DPI: {audit.get('median_dpi_x', 'NA')}",
            f"- Major issues: {audit.get('major_issue_count', 0)}",
            (
                "- Resolution/vector issues: "
                f"{audit.get('resolution_or_vector_issue_count', 0)}"
            ),
            (
                "- Label overlap: "
                f"severe={overlap_counts.get('严重', 0)}, "
                f"minor={overlap_counts.get('轻微', 0)}, "
                f"none={overlap_counts.get('无', 0)}"
            ),
            (
                "- Duplicate SVG text anchors: "
                f"{audit.get('duplicate_text_anchor_panel_count', 0)}"
            ),
            (
                "- Detailed review: "
                "`10_reports/figure_quality_audit/figure_quality_audit.md`"
            ),
        ]
    coverage_path = (
        context.output_root
        / "10_reports"
        / "plan_coverage"
        / "plan_coverage.json"
    )
    if coverage_path.exists():
        coverage = read_json(coverage_path, {})
        lines += [
            "",
            "## Plan Coverage",
            "",
            (
                "- Code/implementation coverage: "
                f"{coverage.get('implementation_completion_percent', 'NA')}%"
            ),
            (
                "- Current executed-result coverage: "
                f"{coverage.get('current_result_completion_percent', 'NA')}%"
            ),
            (
                "- Projected coverage with prerequisites: "
                f"{coverage.get('projected_completion_percent_with_prerequisites', 'NA')}%"
            ),
            (
                "- Performance targets met: "
                f"{coverage.get('performance_targets_met', 'NA')}/"
                f"{coverage.get('performance_targets_evaluable', 'NA')} evaluable; "
                f"{coverage.get('performance_targets_not_evaluated', 'NA')} "
                "not evaluated because the available cohort measures a "
                "different endpoint"
            ),
            "- Detailed panel matrix: `10_reports/plan_coverage/plan_coverage.md`",
        ]
    lines += [
        "",
        "## Important Limits",
        "",
        "- GSE49541 and GSE164441 measure different clinical endpoints from "
        "healthy-versus-NAFLD; their AUC values are not interchangeable.",
        "- GSE135251 is included as a supplementary NAFLD-versus-control "
        "validation cohort.",
        "- GeneCards, OMIM and TTD are not available without licensed or "
        "credentialed bulk access; local files can be supplied in the config.",
          "- The communication panel uses R CellChat and retains cell-level "
          "probabilities as descriptive evidence; group-level significance "
          "still requires independent biological units.",
        "- The Figure 5d-h panels are generated from GROMACS/MM-PBSA outputs "
        "only when the 100 ns production run is actually executed; otherwise "
        "they are explicitly marked not run.",
        "",
        "## Main Folders",
        "",
        "- `01_compound_characterization`",
        "- `02_disease_targets`",
        "- `03_intersection_ppi`",
        "- `04_bulk_training`",
        "- `05_machine_learning`",
        "- `06_single_cell_mouse`",
        "- `07_single_cell_human`",
        "- `08_docking`",
        "- `09_md_mmpbsa`",
        "- `10_reports`",
        "- `按方案分类`（按 Figure 1-5 和 Panel 重组）",
        "",
        "Full file inventory and stage status: `10_reports/experiment_plan_one_report.md`.",
    ]
    return "\n".join(lines) + "\n"


def _figure_index(root: Path) -> pd.DataFrame:
    specs = [
        ("Figure 1", "a", "研究流程", "01_compound_characterization/fig1a_workflow.png"),
        ("Figure 1", "b", "6PPD-Q 2D 结构", "01_compound_characterization/fig1b_compound_2d.png"),
        ("Figure 1", "c", "3D 结构与理化性质", "01_compound_characterization/fig1c_*.png"),
        ("Figure 1", "d", "化合物靶点来源", "01_compound_characterization/fig1d_*.png"),
        ("Figure 1", "e", "疾病靶点来源", "02_disease_targets/fig1e_*.png"),
        ("Figure 1", "f", "化合物-疾病靶点交集", "03_intersection_ppi/fig1f_*.png"),
        ("Figure 1", "g", "KEGG 富集", "03_intersection_ppi/enrichment/fig1g_*.png"),
        ("Figure 1", "h", "GO 富集", "03_intersection_ppi/enrichment/fig1h_*.png"),
        ("Figure 2", "a", "STRING PPI 网络", "03_intersection_ppi/fig2a_*.png"),
        ("Figure 2", "b", "PPI 模块网络", "03_intersection_ppi/fig2b_*.png"),
        ("Figure 2", "c", "Degree 排名", "03_intersection_ppi/fig2c_*.png"),
        ("Figure 2", "d", "Betweenness 排名", "03_intersection_ppi/fig2d_*.png"),
        ("Figure 2", "e", "MCC/Degree 共识 Venn", "03_intersection_ppi/fig2e_*.png"),
        ("Figure 2", "f", "候选基因热图", "04_bulk_training/fig2f_*.png"),
        ("Figure 2", "g", "GSE49541 箱线图", "04_bulk_training/fig2g_*.png"),
        ("Figure 2", "h", "GSE164441 箱线图", "04_bulk_training/fig2h_*.png"),
        ("Figure 3", "a", "11 种算法 AUC 热图", "05_machine_learning/fig3a_*.png"),
        ("Figure 3", "b", "最优模型训练 ROC", "05_machine_learning/fig3b_*.png"),
        ("Figure 3", "c", "GSE49541 外部 ROC", "05_machine_learning/fig3c_*.png"),
        ("Figure 3", "d", "GSE164441 外部 ROC", "05_machine_learning/fig3d_*.png"),
        ("Figure 3", "e", "校准曲线", "05_machine_learning/fig3e_*.png"),
        ("Figure 3", "f", "SHAP 重要性", "05_machine_learning/fig3f_*.png"),
        ("Figure 3", "g", "SHAP 蜂群图", "05_machine_learning/fig3g_*.png"),
        ("Figure 3", "h", "核心基因与纤维化关联", "05_machine_learning/fig3h_*.png"),
        ("Figure 4", "a", "小鼠 UMAP", "06_single_cell_mouse/fig4a_*.png"),
        ("Figure 4", "b", "细胞类型 marker", "06_single_cell_mouse/fig4b_*.png"),
        ("Figure 4", "c", "核心基因小提琴图", "06_single_cell_mouse/fig4c_*.png"),
        ("Figure 4", "d", "核心基因 UMAP", "06_single_cell_mouse/fig4d_*.png"),
        ("Figure 4", "e", "细胞组成", "06_single_cell_mouse/fig4e_*.png"),
        ("Figure 4", "f", "细胞通讯", "06_single_cell_mouse/fig4f_*.png"),
        ("Figure 4", "g", "虚拟敲除", "06_single_cell_mouse/fig4g_*.png"),
        ("Figure 4", "h", "虚拟敲除富集", "06_single_cell_mouse/fig4h_*.png"),
        ("Figure 4", "i", "人类 UMAP", "07_single_cell_human/fig4i_*.png"),
        ("Figure 4", "j", "人类核心基因验证", "07_single_cell_human/fig4j_*.png"),
        ("Figure 5", "a", "3D 对接构象", "08_docking/fig5a_*.png"),
        ("Figure 5", "b", "相互作用平面图", "08_docking/fig5b_*.png"),
        ("Figure 5", "c", "结合能热图", "08_docking/fig5c_*.png"),
        ("Figure 5", "d", "蛋白主链 RMSD", "09_md_mmpbsa/fig5d_*.png"),
        ("Figure 5", "e", "配体 RMSD", "09_md_mmpbsa/fig5e_*.png"),
        ("Figure 5", "f", "关键残基 RMSF", "09_md_mmpbsa/fig5f_*.png"),
        ("Figure 5", "g", "回旋半径 Rg", "09_md_mmpbsa/fig5g_*.png"),
        ("Figure 5", "h", "MM-PBSA", "09_md_mmpbsa/fig5h_*.png"),
    ]
    rows: list[dict[str, str]] = []
    for figure, panel, purpose, pattern in specs:
        matches = sorted(root.glob(pattern))
        status = "available"
        if not matches:
            status = "missing"
        elif all("NOT_RUN" in path.name for path in matches):
            status = "prepared_not_run"
        relative = ";".join(
            str(path.relative_to(root)) for path in matches
        )
        rows.append(
            {
                "figure": figure,
                "panel": panel,
                "purpose": purpose,
                "status": status,
                "file": relative,
            }
        )
    return pd.DataFrame(rows)
