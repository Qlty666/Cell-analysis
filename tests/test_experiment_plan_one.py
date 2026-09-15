#!/usr/bin/env python3
"""Focused tests for the experiment-plan-one analysis package."""

from __future__ import annotations

import re
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import anndata as ad

from experiment_plan_one.common import bh_fdr, split_gene_symbol
from experiment_plan_one.classify import classify_experiment_plan_results
from experiment_plan_one.coverage import audit_plan_coverage
from experiment_plan_one.md_figures import generate_plan_md_figures
from experiment_plan_one.ml import _nested_cv_evaluation, _pipeline, _model_zoo
from experiment_plan_one.pipeline import (
    ExperimentPlanOne,
    default_config,
    merge_config,
)
from experiment_plan_one.single_cell import _cellchat_like_analysis
from experiment_plan_one.targets import (
    _parse_swiss_target_table,
    load_compound_target_file,
    make_venn_figure,
)


class TestExperimentPlanOne(unittest.TestCase):
    def test_split_gene_symbol(self):
        self.assertEqual(split_gene_symbol("EGFR"), "EGFR")
        self.assertEqual(split_gene_symbol("EGFR /// ERBB1"), "EGFR")
        self.assertEqual(split_gene_symbol("N/A"), "")

    def test_bh_fdr(self):
        adjusted = bh_fdr([0.01, 0.02, 0.5])
        self.assertEqual(len(adjusted), 3)
        self.assertTrue((adjusted >= 0).all())
        self.assertTrue((adjusted <= 1).all())
        self.assertLessEqual(adjusted[0], adjusted[-1])

    def test_swiss_target_parser_uses_genecards_anchor(self):
        page = """
        <table id="resultTable"><tr>
        <th>Target</th><th>Gene</th><th>UniProt</th><th>ChEMBL</th>
        <th>Class</th><th>Probability*</th></tr>
        <tr>
        <td>Epidermal growth factor receptor</td>
        <td><a href="https://www.genecards.org/cgi-bin/carddisp.pl?gene=EGFR">EGFR</a></td>
        <td><a href="https://www.uniprot.org/uniprot/P00533">P00533</a></td>
        <td>CHEMBL203</td><td>Kinase</td>
        <td><span>0.2234</span></td>
        </tr></table>
        """
        frame = _parse_swiss_target_table(page)
        self.assertEqual(frame.iloc[0]["gene"], "EGFR")
        self.assertAlmostEqual(frame.iloc[0]["probability"], 0.2234)

    def test_venn_figure_creation(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            output = Path(tmp) / "venn.png"
            make_venn_figure(
                {
                    "A": {"A1", "A2", "A3", "SHARED"},
                    "B": {"SHARED", "B1", "B2", "B3", "B4"},
                },
                output,
                title="test",
            )
            self.assertTrue(output.exists())
            self.assertGreater(output.stat().st_size, 100)
            svg_path = output.with_suffix(".svg")
            self.assertTrue(svg_path.exists())
            svg = svg_path.read_text(encoding="utf-8")
            number_positions = {
                value: x
                for x, value in re.findall(
                    r'<text[^>]* x="([^"]+)"[^>]*>(\d+)</text>',
                    svg,
                )
            }
            self.assertEqual(set(number_positions), {"1", "3", "4"})
            self.assertEqual(len(set(number_positions.values())), 3)

    def test_config_merge(self):
        config = merge_config(default_config(), {"ml": {"seed": 7}})
        self.assertEqual(config["ml"]["seed"], 7)
        self.assertEqual(config["ml"]["cv_folds"], 5)

    def test_stage_signature_changes_with_inputs(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            runner = ExperimentPlanOne(root, default_config())
            source = root / "00_data" / "processed" / "input.csv"
            source.parent.mkdir(parents=True, exist_ok=True)
            source.write_text("gene,value\nEGFR,1\n", encoding="utf-8")
            first = runner.context.stage_signature("ml", ["bulk"])
            source.write_text("gene,value\nEGFR,2\n", encoding="utf-8")
            second = runner.context.stage_signature("ml", ["bulk"])
            self.assertNotEqual(first, second)

    def test_missing_stage_output_detection(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            runner = ExperimentPlanOne(Path(tmp), default_config())
            self.assertIn(
                "02b_evidence/evidence_summary.json",
                runner._missing_stage_outputs("evidence"),
            )

    def test_prepared_md_stage_is_a_successful_terminal_state(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            md_dir = root / "09_md_mmpbsa"
            md_dir.mkdir(parents=True)
            (md_dir / "md_stage_summary.json").write_text("{}", encoding="utf-8")
            (md_dir / "md_plan_figures.json").write_text("{}", encoding="utf-8")
            runner = ExperimentPlanOne(root, default_config())
            with mock.patch.object(
                ExperimentPlanOne,
                "stage_md",
                return_value={"status": "prepared", "mode": "prepare"},
            ):
                results = runner.run(["md"])
            self.assertEqual(results["md"]["status"], "prepared")
            state_path = root / "12_reports" / ".stages" / "md.json"
            state = json.loads(state_path.read_text(encoding="utf-8"))
            self.assertEqual(state["status"], "prepared")

    def test_evidence_stage_uses_local_connector(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            compound_dir = root / "01_compound_characterization"
            disease_dir = root / "02_disease_targets"
            compound_dir.mkdir(parents=True)
            disease_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {
                        "pubchem_cid": "1",
                        "canonical_smiles": "CCO",
                        "isomeric_smiles": "CCO",
                        "inchi_key": "TEST",
                    }
                ]
            ).to_csv(compound_dir / "compound_properties.csv", index=False)
            pd.DataFrame(
                [{"gene": "EGFR", "source_count": 1}]
            ).to_csv(compound_dir / "compound_targets.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "gene": "EGFR",
                        "ensembl_id": "ENSG00000146648",
                        "source_count": 1,
                    }
                ]
            ).to_csv(disease_dir / "disease_targets.csv", index=False)
            local_source = root / "local.csv"
            pd.DataFrame(
                {"gene": ["EGFR"], "score": [0.9]}
            ).to_csv(local_source, index=False)
            local_config = root / "evidence.json"
            local_config.write_text(
                json.dumps(
                    {
                        "strict": True,
                        "sources": {
                            "Local": {
                                "enabled": True,
                                "path": str(local_source),
                                "evidence_type": "direct_bioactivity",
                                "relation": "targets",
                                "tier": "EXPERIMENTAL",
                                "subject_type": "compound",
                                "subject_value": "CID1",
                                "target_column": "gene",
                                "score_column": "score",
                            }
                        },
                        "scoring": {},
                    }
                ),
                encoding="utf-8",
            )
            config = default_config()
            config["evidence"]["config_file"] = str(local_config)
            config["evidence"]["allow_network"] = False
            runner = ExperimentPlanOne(root, config)
            result = runner.stage_evidence()
            self.assertEqual(result["status"], "completed")
            self.assertTrue(
                (root / "02b_evidence" / "target_priority.csv").exists()
            )

    def test_nested_model_selection_produces_outer_predictions(self):
        rng = np.random.default_rng(7)
        X = pd.DataFrame(
            rng.normal(size=(30, 4)),
            columns=["G1", "G2", "G3", "G4"],
        )
        y = pd.Series([0] * 15 + [1] * 15)
        zoo = _model_zoo(seed=7)
        candidates = {
            ("all", "Ridge"): _pipeline(zoo["Ridge"]),
            ("all", "RandomForest"): _pipeline(zoo["RandomForest"]),
        }
        selected, metrics, probabilities = _nested_cv_evaluation(
            candidates,
            {"all": X},
            y,
            seed=7,
            outer_folds=3,
            inner_folds=2,
        )
        self.assertEqual(len(selected), 3)
        self.assertEqual(len(probabilities), len(y))
        self.assertGreaterEqual(metrics["auc"], 0.0)

    def test_optional_boosting_models_can_be_absent(self):
        with mock.patch(
            "experiment_plan_one.ml._xgboost",
            side_effect=RuntimeError("xgboost unavailable"),
        ), mock.patch(
            "experiment_plan_one.ml._lightgbm",
            side_effect=RuntimeError("lightgbm unavailable"),
        ):
            models = _model_zoo(seed=7)
        self.assertNotIn("XGBoost", models)
        self.assertNotIn("LightGBM", models)
        self.assertIn("RandomForest", models)

    def test_classify_results_creates_plan_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "01_compound_characterization"
            source.mkdir(parents=True)
            (source / "fig1a_workflow.png").write_bytes(b"png")
            result = classify_experiment_plan_results(root)
            classified = Path(result["classified_root"])
            self.assertTrue(
                (
                    classified
                    / "Figure1_化合物表征_靶点预测与通路富集"
                    / "Fig1a_研究全局流程图.png"
                ).exists()
            )
            self.assertTrue((classified / "分类清单.csv").exists())

    def test_compound_target_file_loader(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "predictions.csv"
            pd.DataFrame(
                {
                    "gene": ["EGFR", "GPAT3"],
                    "probability": [0.9, 0.4],
                }
            ).to_csv(path, index=False)
            frame, status = load_compound_target_file(
                "LocalPrediction",
                path,
                score_column="probability",
            )
            self.assertEqual(status["status"], "completed")
            self.assertEqual(set(frame["gene"]), {"EGFR", "GPAT3"})
            self.assertEqual(set(frame["source"]), {"LocalPrediction"})

    def test_md_figures_parse_real_xvg_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            docking = root / "08_docking"
            run = docking / "targets" / "GENE1" / "outputs" / "md" / "06_md" / "lig1"
            run.mkdir(parents=True)
            pd.DataFrame(
                {
                    "id": ["lig1"],
                    "mmpbsa_status": ["completed"],
                    "mmpbsa_delta_total_kj_mol": [-42.0],
                }
            ).to_csv(docking / "md_simulation_results.csv", index=False)
            for name, values in {
                "rmsd_protein.xvg": [0.1, 0.11, 0.12],
                "rmsd_ligand.xvg": [0.05, 0.06, 0.07],
                "gyrate_protein.xvg": [2.1, 2.11, 2.12],
            }.items():
                run.joinpath(name).write_text(
                    "\n".join(
                        f"{index * 10} {value}"
                        for index, value in enumerate(values)
                    ),
                    encoding="utf-8",
                )
            run.joinpath("rmsf_protein_residue.xvg").write_text(
                "1 0.1\n2 0.2\n",
                encoding="utf-8",
            )
            result = generate_plan_md_figures(docking, root / "09_md_mmpbsa")
            self.assertTrue(result["panels"]["d_rmsd"])
            self.assertTrue(result["panels"]["e_ligand_rmsd"])
            self.assertTrue(result["panels"]["f_rmsf"])
            self.assertTrue(result["panels"]["g_rg"])
            self.assertTrue(result["panels"]["h_mmpbsa"])

    def test_cellchat_permutation_outputs_fdr(self):
        rng = np.random.default_rng(3)
        data = ad.AnnData(
            X=rng.poisson(1, size=(40, 4)).astype(float),
            obs=pd.DataFrame(
                {
                    "cell_type": ["A"] * 20 + ["B"] * 20,
                    "condition": ["NCD", "HFD"] * 20,
                }
            ),
        )
        data.var_names = ["TNF", "TNFRSF1A", "IL6", "IL6R"]
        result = _cellchat_like_analysis(
            data,
            n_permutations=10,
            seed=1,
        )
        self.assertFalse(result["interactions"].empty)
        self.assertIn("p_value", result["interactions"].columns)
        self.assertIn("fdr", result["interactions"].columns)
        self.assertIn("significant", result["interactions"].columns)

    def test_plan_coverage_separates_implementation_from_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            audit_dir = root / "10_reports" / "figure_quality_audit"
            audit_dir.mkdir(parents=True)
            pd.DataFrame(
                {
                    "figure": [
                        "Figure1_化合物表征_靶点预测与通路富集",
                    ],
                    "panel": ["a"],
                    "status": ["available"],
                    "audit_verdict": ["可用"],
                }
            ).to_csv(audit_dir / "figure_quality_audit.csv", index=False)
            summary = audit_plan_coverage(root)
            self.assertGreaterEqual(
                summary["implementation_completion_percent"],
                0.0,
            )
            self.assertIn("current_result_completion_percent", summary)
            self.assertIn("environment", summary)


if __name__ == "__main__":
    unittest.main()
