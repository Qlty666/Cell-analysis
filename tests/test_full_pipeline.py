#!/usr/bin/env python3
"""Tests for the integrated scRNA -> targets -> knockout pipeline."""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from pipeline.integration import (  # noqa: E402
    STAGES,
    _clear_downstream_markers,
    _dataset_mode_from_root,
    _download_pdb,
    _extract_cocrystal_ligands,
    _invalidate_markers_for_changed_root,
    _resolve_feedback_species,
    _stage_cell_feedback,
    _stage_outdated,
    _stage_output_paths,
    _stage_outputs_ready,
    _stage_signature,
    _valid_docking_box,
    _write_stage_marker,
    _write_run_context,
    build_knockout_inputs,
    collect_qc_metrics,
    evaluate_qc_gate,
    extract_key_genes,
    main,
    run_differential_abundance,
    write_qc_metrics,
)
from pipeline.cell_feedback import (  # noqa: E402
    build_feedback_manifest,
    run_cell_feedback,
)
from pipeline import orchestrator  # noqa: E402


def _write_deg(root: Path) -> None:
    rows = [
        ["GENE1", 0.0, 3.0, "Up", 0.5, 0.1, True],
        ["GENE2", 0.0, -2.5, "Down", 0.6, 0.2, True],
        ["RPLP0", 0.0, 4.0, "Up", 0.9, 0.3, True],
        ["MT-ND1", 0.0, 3.5, "Up", 0.8, 0.2, True],
        ["GENE3", 1e-12, 1.5, "Up", 0.4, 0.2, True],
        ["GENE4", 1e-6, 0.5, "NS", 0.3, 0.2, False],
    ]
    frame = pd.DataFrame(
        rows,
        columns=[
            "gene",
            "p_val_adj",
            "avg_log2FC",
            "direction",
            "pct.1",
            "pct.2",
            "significant",
        ],
    )
    (root / "results" / "data" / "05_deg").mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        root / "results" / "data" / "05_deg" / "fig_09_deg_significant.csv",
        index=False,
    )


def _write_pseudobulk(workdir: Path, n_genes: int = 24) -> None:
    pseudo = workdir / "data" / "knockout" / "_pseudobulk"
    pseudo.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    signature = [
        "MKI67",
        "PCNA",
        "TOP2A",
        "AURKA",
        "CDC20",
        "CDK1",
        "CCNB1",
        "BUB1",
        "BIRC5",
        "CENPA",
    ]
    genes = signature + [f"GENE{i:02d}" for i in range(n_genes)]
    genes = list(dict.fromkeys(genes))[:n_genes]
    rows = {"gene": genes}
    for sample in ["T1", "T2", "N1", "N2"]:
        base = 4.0 if sample.startswith("T") else 1.0
        rows[sample] = base + rng.normal(0, 0.2, len(genes))
    expr = pd.DataFrame(rows)
    expr.loc[expr["gene"].isin(signature), ["T1", "T2"]] += 3.0
    expr.to_csv(pseudo / "pseudobulk_expression.csv", index=False)
    pd.DataFrame(
        {
            "sample": ["T1", "T2", "N1", "N2"],
            "condition": ["Tumor", "Tumor", "Normal", "Normal"],
            "cell_type": ["Hepatocyte", "T_NK", "Hepatocyte", "T_NK"],
        }
    ).to_csv(pseudo / "pseudobulk_metadata.csv", index=False)


def _pipeline_args(
    config_path: Path,
    docking_config: Path,
    output: Path,
    workdir: Path,
    **overrides,
) -> argparse.Namespace:
    values = {
        "config": str(config_path),
        "docking_config": str(docking_config),
        "accession": "GSE999999",
        "species": "hs",
        "skip_scrna": True,
        "skip_download": False,
        "skip_deps": True,
        "qc_gate": {"enabled": True},
        "differential_abundance": {"enabled": True},
        "top_genes": 10,
        "keep_all_genes": False,
        "gene_blacklist": [],
        "skip_evidence_fetch": True,
        "evidence_workers": 2,
        "evidence_timeout": 30,
        "skip_pseudobulk": False,
        "case_label": None,
        "normal_label": None,
        "ko_top_n": None,
        "depmap_csv": None,
        "skip_knockout": True,
        "docking_targets": 0,
        "ligand_library": None,
        "skip_docking": True,
        "feedback_top_n": 5,
        "feedback_max_features": 3,
        "feedback_timeout": 60,
        "skip_cell_feedback": True,
        "dry_run": False,
        "force": False,
        "start_stage": None,
        "output": str(output),
        "workdir": str(workdir),
    }
    values.update(overrides)
    return argparse.Namespace(**values)


class TestExtractKeyGenes(unittest.TestCase):
    def test_rank_and_blacklist(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            root.mkdir(parents=True)
            _write_deg(root)
            out = root / "integration_out"
            frame = extract_key_genes(root, out, top_n=10)
            self.assertNotIn("RPLP0", frame["gene"].tolist())
            self.assertNotIn("MT-ND1", frame["gene"].tolist())
            self.assertNotIn("GENE4", frame["gene"].tolist())
            self.assertEqual(frame.iloc[0]["gene"], "GENE1")
            self.assertAlmostEqual(float(frame.iloc[0]["avg_log2fc"]), 3.0)
            self.assertTrue((out / "key_genes.csv").exists())
            self.assertTrue((out / "key_genes_summary.json").exists())

    def test_duplicate_log2fc_columns_are_deduped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            data_dir = root / "results" / "data"
            (data_dir / "05_deg").mkdir(parents=True)
            frame = pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2"],
                    "avg_log2FC": [3.0, -2.0],
                    "log2FoldChange": [2.9, -1.9],
                    "p_val_adj": [0.01, 0.02],
                    "pct.1": [0.5, 0.6],
                    "pct.2": [0.2, 0.3],
                    "direction": ["Up", "Down"],
                    "significant": [True, True],
                }
            )
            frame.to_csv(
                data_dir / "05_deg" / "fig_09_deg_significant.csv",
                index=False,
            )
            out = root / "integration_out"
            result = extract_key_genes(root, out, top_n=10)
            self.assertEqual(result.iloc[0]["gene"], "GENE1")
            self.assertAlmostEqual(float(result.iloc[0]["avg_log2fc"]), 3.0)

    def test_empty_significant_table_falls_back_to_all_deg_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            data_dir = root / "results" / "data"
            (data_dir / "05_deg").mkdir(parents=True)
            pd.DataFrame(
                columns=[
                    "gene",
                    "p_val_adj",
                    "avg_log2FC",
                    "direction",
                    "pct.1",
                    "pct.2",
                    "significant",
                ]
            ).to_csv(
                data_dir / "05_deg" / "fig_09_deg_significant.csv",
                index=False,
            )
            all_frame = pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2"],
                    "p_val_adj": [0.01, 0.02],
                    "avg_log2FC": [2.0, -1.5],
                    "direction": ["Up", "Down"],
                    "pct.1": [0.5, 0.6],
                    "pct.2": [0.2, 0.3],
                    "significant": [True, True],
                }
            )
            all_frame.to_csv(
                data_dir / "05_deg" / "fig_08_deg_all.csv",
                index=False,
            )
            out = root / "integration_out"
            result = extract_key_genes(root, out, top_n=10)
            self.assertEqual(result.iloc[0]["gene"], "GENE1")
            self.assertGreater(len(result), 0)


class TestCocrystalLigandFallback(unittest.TestCase):
    def test_extract_hetatm_ligand(self):
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp) / "lig.pdb"
            pdb.write_text(
                "HETATM    1  C1  LIG A 101       1.000   2.000   3.000  1.00 20.00           C\n"
                "HETATM    2  N1  LIG A 101       2.000   2.000   3.000  1.00 20.00           N\n"
                "HETATM    3  C2  LIG A 101       3.000   2.000   3.000  1.00 20.00           C\n"
                "HETATM    4  O1  LIG A 101       4.000   2.000   3.000  1.00 20.00           O\n"
                "HETATM    5  C3  LIG A 101       1.000   1.000   3.000  1.00 20.00           C\n"
                "HETATM    6  C4  LIG A 101       1.000   3.000   3.000  1.00 20.00           C\n"
                "CONECT    1    2    5    6\n"
                "CONECT    2    1    3\n"
                "CONECT    3    2    4\n"
                "CONECT    4    3\n"
                "END\n",
                encoding="utf-8",
            )
            ligands = _extract_cocrystal_ligands(pdb)
            self.assertEqual(len(ligands), 1)
            self.assertTrue(ligands[0]["smiles"])


class TestBuildKnockoutInputs(unittest.TestCase):
    def test_inputs_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            root.mkdir(parents=True)
            _write_deg(root)
            workdir = Path(tmp) / "work"
            _write_pseudobulk(workdir)
            evidence = pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2"],
                    "uniprot": ["P00001", ""],
                    "known_ligands": [5, 0],
                    "chembl_bioactivities": [5, 0],
                    "pdb_structures": [2, 0],
                    "pdb_ids": ["1ABC,2DEF", ""],
                    "off_target_paralogs": [1, 0],
                    "safety_concern": [0, 0],
                    "entrez": ["1", "2"],
                    "ensembl": ["ENSG1", "ENSG2"],
                    "chembl_target_id": ["CHEMBL1", ""],
                }
            )
            (workdir / "outputs" / "integration").mkdir(parents=True, exist_ok=True)
            evidence.to_csv(
                workdir / "outputs" / "integration" / "gene_evidence.csv",
                index=False,
            )
            summary = build_knockout_inputs(root, workdir)
            ko_dir = workdir / "data" / "knockout"
            self.assertTrue((ko_dir / "expression.csv").exists())
            self.assertTrue((ko_dir / "metadata.csv").exists())
            self.assertTrue((ko_dir / "prognosis.csv").exists())
            self.assertTrue((ko_dir / "druggability.csv").exists())
            self.assertTrue((ko_dir / "off_target.csv").exists())
            self.assertEqual(summary["samples"], 4)
            druggable = pd.read_csv(ko_dir / "druggability.csv")
            self.assertEqual(
                int(druggable.loc[druggable["gene"] == "GENE1", "known_ligands"].iloc[0]),
                5,
            )


class TestCellFeedbackManifest(unittest.TestCase):
    def _write_screen_results(self, workdir: Path) -> None:
        ko_dir = (
            workdir
            / "outputs"
            / "run_001"
            / "results"
            / "04_knockout"
            / "data"
        )
        ko_dir.mkdir(parents=True)
        pd.DataFrame(
            {
                "rank": [1, 2],
                "gene": ["GENE1", "GENE2"],
                "target_score": [0.9, 0.7],
                "knockout_score": [0.8, 0.6],
                "target_class": ["core_driver", "biomarker"],
            }
        ).to_csv(ko_dir / "fig_52_53_ranked_knockout.csv", index=False)
        integration = workdir / "outputs" / "integration"
        integration.mkdir(parents=True)
        pd.DataFrame(
            {
                "gene": ["GENE1", "GENE2"],
                "status": ["ok", "failed"],
                "hits": [5, 0],
                "best_affinity": ["-8.5", ""],
                "pdb_id": ["1ABC", ""],
            }
        ).to_csv(integration / "docking_targets.csv", index=False)

    def test_manifest_merges_knockout_and_docking(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            self._write_screen_results(workdir)
            frame = build_feedback_manifest(workdir, top_n=10)
            self.assertEqual(frame.iloc[0]["gene"], "GENE1")
            self.assertIn("knockout", frame.iloc[0]["source"])
            self.assertIn("docking", frame.iloc[0]["source"])
            self.assertEqual(int(frame.iloc[0]["docking_hits"]), 5)
            self.assertGreater(frame.iloc[0]["feedback_score"], 0.9)
            self.assertEqual(int(frame.iloc[1]["docking_hits"]), 0)

    def test_run_cell_feedback_skips_without_screen_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            root = Path(tmp) / "single_cell"
            summary = run_cell_feedback(workdir, root, top_n=10)
            self.assertEqual(summary["status"], "skipped")
            self.assertTrue(
                (
                    workdir
                    / "outputs"
                    / "integration"
                    / "cell_feedback"
                    / "cell_feedback_summary.json"
                ).exists()
            )

    def test_feedback_species_resolves_from_dataset_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            data_dir = root / "data"
            data_dir.mkdir(parents=True)
            (data_dir / "GSE123456_manifest.json").write_text(
                json.dumps({"organism": "mm"}),
                encoding="utf-8",
            )
            ctx = {"single_cell_root": root}
            args = argparse.Namespace(
                species="auto",
                accession="GSE123456",
            )
            self.assertEqual(_resolve_feedback_species(args, ctx), "mm")

            args.species = "hs"
            self.assertEqual(_resolve_feedback_species(args, ctx), "hs")


class TestFullPipelineMarkers(unittest.TestCase):
    def _write_markers(self, workdir: Path) -> None:
        stage_dir = workdir / "outputs" / "integration" / ".stages"
        stage_dir.mkdir(parents=True, exist_ok=True)
        for code, name, _description in STAGES:
            (stage_dir / f"{code}_{name}.done").write_text(
                "done",
                encoding="utf-8",
            )

    def test_changed_single_cell_root_resets_stage_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            workdir = base / "work"
            self._write_markers(workdir)
            _write_run_context(workdir, base / "old_root")
            changed = _invalidate_markers_for_changed_root(
                workdir,
                base / "new_root",
            )
            self.assertTrue(changed)
            stage_dir = workdir / "outputs" / "integration" / ".stages"
            self.assertEqual(list(stage_dir.glob("*.done")), [])

    def test_changed_root_detected_from_key_genes_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            workdir = base / "work"
            self._write_markers(workdir)
            integration = workdir / "outputs" / "integration"
            old_root = base / "old_root"
            (integration / "key_genes_summary.json").write_text(
                json.dumps(
                    {
                        "deg_table": str(
            (
                old_root
                / "results"
                / "data"
                / "05_deg"
                / "fig_09_deg_significant.csv"
            )
                        )
                    }
                ),
                encoding="utf-8",
            )
            changed = _invalidate_markers_for_changed_root(
                workdir,
                base / "new_root",
            )
            self.assertTrue(changed)
            stage_dir = workdir / "outputs" / "integration" / ".stages"
            self.assertEqual(list(stage_dir.glob("*.done")), [])

    def test_clear_downstream_markers_keeps_stage01(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            self._write_markers(workdir)
            _clear_downstream_markers(workdir)
            stage_dir = workdir / "outputs" / "integration" / ".stages"
            remaining = {path.stem for path in stage_dir.glob("*.done")}
            self.assertEqual(remaining, {"01_single_cell"})


class TestFullPipeline(unittest.TestCase):
    def test_bulk_manifest_is_supported_by_pipeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            (root / "data").mkdir(parents=True)
            (root / "data" / "GSE299321_manifest.json").write_text(
                json.dumps(
                    {
                        "accession": "GSE299321",
                        "mode": "bulk",
                        "organism": "hs",
                        "files": {
                            "matrix": ["GSM9037276_GC119559.txt.gz"],
                            "barcodes": [],
                            "genes": [],
                            "metadata": [],
                            "series_matrices": ["GSE299321_series_matrix.txt.gz"],
                        },
                    }
                ),
                encoding="utf-8",
            )
            with (
                mock.patch.object(
                    orchestrator,
                    "run_r_pipeline",
                    return_value=(0, Path(tmp) / "pipeline_r.log"),
                ) as mock_r,
                mock.patch.object(orchestrator, "verify_outputs"),
                mock.patch.object(orchestrator, "run_ml_analysis"),
                mock.patch.object(orchestrator, "generate_report"),
            ):
                orchestrator.run_pipeline(
                    force=False,
                    skip_download=True,
                    skip_deps=True,
                    accession="GSE299321",
                    output_root=str(root),
                    species="auto",
                )
            self.assertEqual(
                mock_r.call_count,
                1,
            )

    def test_biostudies_accession_is_accepted_by_pipeline(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            (root / "data").mkdir(parents=True)
            (root / "data" / "E-MTAB-1234_manifest.json").write_text(
                json.dumps(
                    {
                        "accession": "E-MTAB-1234",
                        "mode": "bulk",
                        "organism": "hs",
                        "files": {
                            "matrix": ["counts.txt"],
                            "barcodes": [],
                            "genes": [],
                            "metadata": [],
                            "series_matrices": [],
                        },
                    }
                ),
                encoding="utf-8",
            )
            with (
                mock.patch.object(
                    orchestrator,
                    "run_r_pipeline",
                    return_value=(0, Path(tmp) / "pipeline_r.log"),
                ) as mock_r,
                mock.patch.object(orchestrator, "verify_outputs"),
                mock.patch.object(orchestrator, "run_ml_analysis"),
                mock.patch.object(orchestrator, "generate_report"),
            ):
                orchestrator.run_pipeline(
                    force=False,
                    skip_download=True,
                    skip_deps=True,
                    accession="E-MTAB-1234",
                    output_root=str(root),
                    species="auto",
                )
            self.assertEqual(mock_r.call_args.args[0], "E-MTAB-1234")

    def test_egood_accession_is_canonicalized_to_gse(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            (root / "data").mkdir(parents=True)
            (root / "data" / "GSE1_manifest.json").write_text(
                json.dumps(
                    {
                        "accession": "GSE1",
                        "mode": "bulk",
                        "organism": "hs",
                        "files": {
                            "matrix": ["counts.txt"],
                            "barcodes": [],
                            "genes": [],
                            "metadata": [],
                            "series_matrices": [],
                        },
                    }
                ),
                encoding="utf-8",
            )
            with (
                mock.patch.object(
                    orchestrator,
                    "run_r_pipeline",
                    return_value=(0, Path(tmp) / "pipeline_r.log"),
                ) as mock_r,
                mock.patch.object(orchestrator, "verify_outputs"),
                mock.patch.object(orchestrator, "run_ml_analysis"),
                mock.patch.object(orchestrator, "generate_report"),
            ):
                orchestrator.run_pipeline(
                    force=False,
                    skip_download=True,
                    skip_deps=True,
                    accession="E-GEOD-1",
                    output_root=str(root),
                    species="auto",
                )
            self.assertEqual(mock_r.call_args.args[0], "GSE1")

    def test_dataset_mode_from_bulk_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            (root / "results").mkdir(parents=True)
            (root / "results" / "summary.json").write_text(
                json.dumps({"dataset_mode": "sample_level"}),
                encoding="utf-8",
            )
            self.assertEqual(_dataset_mode_from_root(root), "sample_level")

    def test_bulk_stage_cell_feedback_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            root = Path(tmp) / "single_cell"
            args = argparse.Namespace(skip_cell_feedback=False)
            ctx = {"single_cell_root": root, "dataset_mode": "sample_level"}
            _stage_cell_feedback(args, workdir, ctx)
            summary_path = (
                workdir
                / "outputs"
                / "integration"
                / "cell_feedback"
                / "cell_feedback_summary.json"
            )
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["status"], "skipped")
            self.assertIn("sample-level", summary["reason"])

    def test_invalid_accession_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            root.mkdir(parents=True)
            workdir = Path(tmp) / "work"
            cfg = Path(tmp) / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps(
                    {
                        "accession": "../../evil",
                        "single_cell_output": str(root),
                        "workdir": str(workdir),
                        "top_genes": 10,
                        "docking_targets": 2,
                        "species": "hs",
                    }
                ),
                encoding="utf-8",
            )
            code = main(
                [
                    "--config",
                    str(cfg),
                    "--skip-scrna",
                    "--skip-docking",
                    "--skip-evidence-fetch",
                    "--skip-cell-feedback",
                    "--docking-config",
                    str(APP_ROOT / "config" / "docking_config.json"),
                ]
            )
            self.assertEqual(code, 1)

    def test_single_cell_output_required_when_config_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            cfg = Path(tmp) / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps(
                    {
                        "accession": "GSE999999",
                        "single_cell_output": "",
                        "workdir": str(workdir),
                        "top_genes": 10,
                        "docking_targets": 2,
                        "species": "hs",
                    }
                ),
                encoding="utf-8",
            )
            code = main(
                [
                    "--config",
                    str(cfg),
                    "--skip-scrna",
                    "--skip-docking",
                    "--skip-evidence-fetch",
                    "--skip-cell-feedback",
                    "--docking-config",
                    str(APP_ROOT / "config" / "docking_config.json"),
                ]
            )
            self.assertEqual(code, 1)

    def test_workdir_required_when_config_empty(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            root.mkdir(parents=True)
            cfg = Path(tmp) / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps(
                    {
                        "accession": "GSE999999",
                        "single_cell_output": str(root),
                        "workdir": "",
                        "top_genes": 10,
                        "docking_targets": 2,
                        "species": "hs",
                    }
                ),
                encoding="utf-8",
            )
            code = main(
                [
                    "--config",
                    str(cfg),
                    "--skip-scrna",
                    "--skip-docking",
                    "--skip-evidence-fetch",
                    "--docking-config",
                    str(APP_ROOT / "config" / "docking_config.json"),
                ]
            )
            self.assertEqual(code, 1)

    def test_end_to_end_without_docking(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            root.mkdir(parents=True)
            _write_deg(root)
            (root / "results").mkdir(parents=True, exist_ok=True)
            (root / "results" / "pipeline_complete.json").write_text(
                json.dumps({"ok": True}),
                encoding="utf-8",
            )
            (root / "results" / "summary.json").write_text(
                json.dumps({"dataset": "GSE999999", "n_genes": 100}),
                encoding="utf-8",
            )
            workdir = Path(tmp) / "work"
            _write_pseudobulk(workdir)
            cfg = Path(tmp) / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps(
                    {
                        "accession": "GSE999999",
                        "single_cell_output": str(root),
                        "workdir": str(workdir),
                        "top_genes": 10,
                        "docking_targets": 2,
                        "species": "hs",
                    }
                ),
                encoding="utf-8",
            )
            code = main(
                [
                    "--config",
                    str(cfg),
                    "--skip-scrna",
                    "--skip-docking",
                    "--skip-evidence-fetch",
                    "--skip-cell-feedback",
                    "--docking-config",
                    str(APP_ROOT / "config" / "docking_config.json"),
                ]
            )
            self.assertEqual(code, 0)
            out = workdir / "outputs" / "integration"
            self.assertTrue((out / "key_genes.csv").exists())
            self.assertTrue((out / "gene_evidence.csv").exists())
            self.assertTrue(
                (
                    workdir
                    / "outputs"
                    / "run_001"
                    / "results"
                    / "04_knockout"
                    / "data"
                    / "fig_52_53_ranked_knockout.csv"
                ).exists()
            )
            self.assertTrue(
                (
                    workdir
                    / "outputs"
                    / "run_001"
                    / "results"
                    / "05_validation"
                    / "data"
                    / "validation_plan.md"
                ).exists()
            )
            self.assertTrue((out / "integration_report.html").exists())
            self.assertTrue((out / "integration_summary.json").exists())
            self.assertTrue((out / ".stages" / "08_report.done").exists())
            evidence = pd.read_csv(out / "gene_evidence.csv")
            for column in [
                "string_partners",
                "reactome_pathways",
                "kegg_pathways",
                "database_sources",
            ]:
                self.assertIn(column, evidence.columns)

    def test_stale_stage01_marker_reruns_when_outputs_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            root.mkdir(parents=True)
            workdir = Path(tmp) / "work"
            stage_dir = workdir / "outputs" / "integration" / ".stages"
            stage_dir.mkdir(parents=True)
            (stage_dir / "01_single_cell.done").write_text(
                "stale", encoding="utf-8"
            )
            cfg = Path(tmp) / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps(
                    {
                        "accession": "GSE999999",
                        "single_cell_output": str(root),
                        "workdir": str(workdir),
                        "top_genes": 10,
                        "docking_targets": 2,
                        "species": "hs",
                    }
                ),
                encoding="utf-8",
            )
            code = main(
                [
                    "--config",
                    str(cfg),
                    "--skip-scrna",
                    "--skip-docking",
                    "--skip-evidence-fetch",
                    "--docking-config",
                    str(APP_ROOT / "config" / "docking_config.json"),
                ]
            )
            self.assertEqual(code, 1)
            self.assertFalse((stage_dir / "01_single_cell.done").exists())


class TestStageProvenance(unittest.TestCase):
    def test_signature_change_invalidates_stage(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            cfg = base / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps({"accession": "GSE999999"}),
                encoding="utf-8",
            )
            dock = base / "docking_config.json"
            dock.write_text("{}", encoding="utf-8")
            workdir = base / "work"
            output = base / "single_cell"
            output.mkdir()
            args = _pipeline_args(cfg, dock, output, workdir)
            ctx = {
                "single_cell_root": output,
                "workdir": workdir,
                "docking_config": dock,
            }
            signature = _stage_signature("02", args, workdir, ctx)
            _write_stage_marker(workdir, "02", "key_targets", signature)
            outdated, _ = _stage_outdated(
                workdir,
                "02",
                "key_targets",
                signature,
                True,
            )
            self.assertFalse(outdated)

            args.top_genes = 25
            new_signature = _stage_signature("02", args, workdir, ctx)
            self.assertNotEqual(signature, new_signature)
            outdated, reason = _stage_outdated(
                workdir,
                "02",
                "key_targets",
                new_signature,
                True,
            )
            self.assertTrue(outdated)
            self.assertIn("parameters", reason)

    def test_legacy_plain_marker_is_outdated(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            cfg = base / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps({"accession": "GSE999999"}),
                encoding="utf-8",
            )
            dock = base / "docking_config.json"
            dock.write_text("{}", encoding="utf-8")
            workdir = base / "work"
            output = base / "single_cell"
            output.mkdir()
            args = _pipeline_args(cfg, dock, output, workdir)
            marker_dir = (
                workdir / "outputs" / "integration" / ".stages"
            )
            marker_dir.mkdir(parents=True, exist_ok=True)
            (marker_dir / "03_evidence.done").write_text(
                "old plain marker",
                encoding="utf-8",
            )
            signature = _stage_signature("03", args, workdir, {})
            outdated, reason = _stage_outdated(
                workdir,
                "03",
                "evidence",
                signature,
                True,
            )
            self.assertTrue(outdated)
            self.assertIn("no marker", reason)


class TestQcGate(unittest.TestCase):
    def _write_summary(self, root: Path) -> None:
        (root / "results").mkdir(parents=True, exist_ok=True)
        (root / "results" / "summary.json").write_text(
            json.dumps(
                {
                    "n_cells_after_qc": 50,
                    "n_genes": 500,
                    "deg_total": 10,
                }
            ),
            encoding="utf-8",
        )

    def test_fail_when_thresholds_violated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            workdir = Path(tmp) / "work"
            self._write_summary(root)
            config = {
                "enabled": True,
                "min_cells_after_qc": 100,
                "min_genes": 1000,
                "min_deg_genes": 5,
                "max_doublet_rate": 0.3,
                "require_pseudobulk": False,
                "fail_on_missing_metrics": False,
            }
            metrics = collect_qc_metrics(root, workdir)
            gate = evaluate_qc_gate(metrics, config)
            self.assertEqual(gate["status"], "fail")
            written = write_qc_metrics(workdir, root, config)
            self.assertEqual(written["qc_gate"]["status"], "fail")
            self.assertTrue(
                (
                    workdir / "outputs" / "integration" / "qc_metrics.json"
                ).exists()
            )

    def test_pass_when_metrics_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            workdir = Path(tmp) / "work"
            self._write_summary(root)
            config = {
                "enabled": True,
                "min_cells_after_qc": 10,
                "min_genes": 100,
                "min_deg_genes": 1,
                "max_doublet_rate": None,
                "require_pseudobulk": False,
                "fail_on_missing_metrics": False,
            }
            gate = evaluate_qc_gate(collect_qc_metrics(root, workdir), config)
            self.assertEqual(gate["status"], "pass")


class TestDifferentialAbundance(unittest.TestCase):
    def test_composition_shift_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "single_cell"
            ann_dir = root / "results" / "data" / "04_annotation"
            ann_dir.mkdir(parents=True)
            rows = []
            for i in range(40):
                rows.append(
                    {
                        "cell": f"t{i}",
                        "celltype_annot": "T_cell",
                        "condition": "Tumor",
                    }
                )
            for i in range(10):
                rows.append(
                    {
                        "cell": f"h{i}",
                        "celltype_annot": "Hepatocyte",
                        "condition": "Tumor",
                    }
                )
            for i in range(10):
                rows.append(
                    {
                        "cell": f"n{i}",
                        "celltype_annot": "T_cell",
                        "condition": "Normal",
                    }
                )
            for i in range(40):
                rows.append(
                    {
                        "cell": f"nh{i}",
                        "celltype_annot": "Hepatocyte",
                        "condition": "Normal",
                    }
                )
            pd.DataFrame(rows).to_csv(
                ann_dir / "fig_05_16_17_cell_annotations.csv",
                index=False,
            )
            out = Path(tmp) / "integration"
            summary = run_differential_abundance(
                root,
                out,
                {"min_cells": 5, "fdr": 0.05},
            )
            self.assertEqual(summary["status"], "completed")
            self.assertGreaterEqual(summary["celltypes_tested"], 2)
            self.assertGreaterEqual(summary["significant_celltypes"], 1)
            frame = pd.read_csv(out / "differential_abundance.csv")
            self.assertLessEqual(frame["p_adjust"].max(), 1.0)
            self.assertGreaterEqual(frame["p_adjust"].min(), 0.0)


class TestDockingRobustness(unittest.TestCase):
    def test_valid_docking_box(self):
        self.assertTrue(_valid_docking_box([1.0, 2.0, 3.0], [20, 20, 20]))
        self.assertFalse(_valid_docking_box([0, 0, 0], [0, 0, 0]))
        self.assertFalse(_valid_docking_box([float("nan"), 0, 0], [20, 20, 20]))

    def test_pdb_download_retries_and_writes(self):
        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def read(self):
                return (
                    b"ATOM      1  N   MET A   1       1.000   2.000   3.000"
                    b"  1.00 20.00           N\nEND\n"
                )

        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch("pipeline.integration.time.sleep", return_value=None):
                with mock.patch(
                    "pipeline.integration.urllib.request.urlopen",
                    side_effect=[OSError("boom"), FakeResponse()],
                ) as urlopen:
                    out = _download_pdb("1ABC", Path(tmp), timeout=5)
            self.assertIsNotNone(out)
            self.assertEqual(urlopen.call_count, 2)
            self.assertTrue(out.read_text(encoding="utf-8").startswith("ATOM"))


class TestDryRun(unittest.TestCase):
    def test_dry_run_returns_zero_without_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            root = base / "single_cell"
            workdir = base / "work"
            cfg = base / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps(
                    {
                        "accession": "GSE999999",
                        "single_cell_output": str(root),
                        "workdir": str(workdir),
                        "species": "hs",
                    }
                ),
                encoding="utf-8",
            )
            dock = base / "docking_config.json"
            dock.write_text("{}", encoding="utf-8")
            code = main(
                [
                    "--config",
                    str(cfg),
                    "--docking-config",
                    str(dock),
                    "--skip-scrna",
                    "--dry-run",
                ]
            )
            self.assertEqual(code, 0)


class TestStageOutputVerification(unittest.TestCase):
    def test_output_paths_are_checked(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            cfg = base / "full_pipeline_config.json"
            cfg.write_text(
                json.dumps({"accession": "GSE999999"}),
                encoding="utf-8",
            )
            dock = base / "docking_config.json"
            dock.write_text("{}", encoding="utf-8")
            workdir = base / "work"
            output = base / "single_cell"
            output.mkdir()
            args = _pipeline_args(cfg, dock, output, workdir)
            ctx = {"single_cell_root": output}
            paths = _stage_output_paths("03", workdir, ctx, args)
            self.assertEqual(len(paths), 1)
            self.assertFalse(_stage_outputs_ready("03", workdir, ctx, args))
            paths[0].parent.mkdir(parents=True, exist_ok=True)
            paths[0].write_text("evidence", encoding="utf-8")
            self.assertTrue(_stage_outputs_ready("03", workdir, ctx, args))


if __name__ == "__main__":
    unittest.main()
