#!/usr/bin/env python3
"""Focused tests for the experiment-plan-one analysis package."""

from __future__ import annotations

import re
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import anndata as ad

from experiment_plan_one.common import bh_fdr, split_gene_symbol
from experiment_plan_one.bulk import (
    _parse_paired_count_sample,
    differential_expression_limma,
)
from experiment_plan_one import cache as plan_cache
from experiment_plan_one.classify import classify_experiment_plan_results
from experiment_plan_one.coexpression import run_coexpression_analysis
from experiment_plan_one.coverage import audit_plan_coverage
from experiment_plan_one.figure_audit import _dynamic_result_reviews, _summary
from experiment_plan_one.md_figures import generate_plan_md_figures
from experiment_plan_one.ml import (
    _model_zoo,
    _nested_cv_evaluation,
    _pipeline,
    run_ml_validation,
)
from experiment_plan_one.pipeline import (
    ExperimentPlanOne,
    default_config,
    merge_config,
)
from experiment_plan_one.planning import (
    audit_delivery_readiness,
    load_governance,
    write_analysis_plan,
    write_governance_artifacts,
)
from experiment_plan_one.ppi import run_ppi_analysis
from experiment_plan_one.single_cell import (
    _cellchat_like_analysis,
    _donor_expression,
    _patient_id,
    map_human_to_mouse_homologs,
)
from experiment_plan_one.targets import (
    _parse_swiss_target_table,
    _sea_result_frame,
    collect_compound_targets,
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

    def test_sea_result_frame_keeps_only_significant_human_targets(self):
        rows = [
            {
                "Target_Chembl_ID": "CHEMBL210",
                "Target_Species": "Homo",
                "Z_score": 3.2,
                "p_value": 0.003,
                "Target_Name": "Beta-2 adrenergic receptor",
                "MaxTc": 0.36,
                "Uniprot_Accession": "P07550",
            },
            {
                "Target_Chembl_ID": "CHEMBL3440",
                "Target_Species": "Mus",
                "Z_score": 4.0,
                "p_value": 0.001,
                "Target_Name": "Beta-1 adrenergic receptor",
                "MaxTc": 0.34,
                "Uniprot_Accession": "P34971",
            },
            {
                "Target_Chembl_ID": "CHEMBL999",
                "Target_Species": "Homo",
                "Z_score": 0.2,
                "p_value": 0.5,
                "Target_Name": "Non-significant target",
                "MaxTc": 0.1,
                "Uniprot_Accession": "P00000",
            },
        ]
        with mock.patch(
            "experiment_plan_one.targets._chembl_target_gene",
            side_effect=lambda target_id: "ADRB2"
            if target_id == "CHEMBL210"
            else "",
        ), mock.patch(
            "experiment_plan_one.targets._map_uniprot_symbols",
            return_value={},
        ):
            frame = _sea_result_frame(
                rows,
                p_value_threshold=0.05,
                min_zscore=0.0,
                max_targets=10,
            )
        self.assertEqual(frame["gene"].tolist(), ["ADRB2"])
        self.assertEqual(frame["source"].tolist(), ["SEA"])
        self.assertAlmostEqual(frame.iloc[0]["p_value"], 0.003)

    def test_collect_compound_targets_can_run_sea_as_independent_source(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            output = Path(tmp) / "targets"
            sea = pd.DataFrame(
                {
                    "gene": ["ADRB2"],
                    "source": ["SEA"],
                    "probability": [2.5],
                }
            )
            with mock.patch(
                "experiment_plan_one.targets._fetch_sea_targets",
                return_value=(
                    sea,
                    {
                        "status": "completed",
                        "task_id": "test",
                        "targets": 1,
                    },
                ),
            ):
                combined, statuses = collect_compound_targets(
                    {
                        "canonical_smiles": "CCO",
                        "iupac_name": "test",
                    },
                    output,
                    enabled_sources=["SEA"],
                )
            self.assertEqual(combined["gene"].tolist(), ["ADRB2"])
            self.assertEqual(statuses["SEA"]["status"], "completed")
            self.assertEqual(
                statuses["SwissTargetPrediction"]["status"],
                "disabled",
            )

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

    def test_analysis_plan_freezes_endpoint_roles_and_manifest(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            outputs = write_analysis_plan(root, default_config())
            for path in outputs.values():
                self.assertTrue(path.exists(), path)
            plan_text = outputs["analysis_plan"].read_text(encoding="utf-8")
            manifest = pd.read_csv(outputs["cohort_manifest"], sep="\t")
            self.assertIn("same_endpoint_external_candidate", plan_text)
            self.assertIn("different_endpoints_are_not_external_validation", plan_text)
            self.assertIn("GSE164441", set(manifest["accession"]))

    def test_governance_freeze_register_and_method_changes_are_complete(self):
        governance = load_governance()
        self.assertEqual(governance["freeze"]["state"], "frozen")
        self.assertTrue(governance["freeze"]["analysis_modules_locked"])
        self.assertEqual(
            set(governance["freeze"]["allowed_change_categories"]),
            {
                "code_defect",
                "data_alignment",
                "status_determination",
                "traceability",
            },
        )
        self.assertEqual(len(governance["issue_register"]), 26)
        self.assertGreaterEqual(
            governance["issue_register"]["status"].isin(
                {"blocked", "not_run"}
            ).sum(),
            8,
        )
        self.assertFalse(
            governance["method_changes"]["equivalent"]
            .astype(str)
            .str.lower()
            .eq("true")
            .any()
        )

    def test_governance_artifacts_are_written_to_run(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            paths = write_governance_artifacts(root)
            for key in (
                "script_freeze",
                "delivery_tiers",
                "issue_register",
                "method_changes",
                "verification",
                "issue_register_summary",
            ):
                self.assertTrue(paths[key].exists(), key)

    def test_delivery_readiness_never_claims_publication_without_rerun(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            readiness = audit_delivery_readiness(Path(tmp))
            self.assertFalse(readiness["publication_grade"])
            self.assertEqual(readiness["tiers"]["T1"]["status"], "met")
            self.assertEqual(readiness["tiers"]["T2"]["status"], "not_run")
            self.assertEqual(readiness["tiers"]["T3"]["status"], "blocked")

    def test_paired_count_sample_parser_accepts_explicit_n_t_labels(self):
        self.assertEqual(
            _parse_paired_count_sample("651N_count"),
            ("651", "adjacent_normal"),
        )
        self.assertEqual(
            _parse_paired_count_sample("651T_count"),
            ("651", "tumor"),
        )
        with self.assertRaises(ValueError):
            _parse_paired_count_sample("651_count")

    @unittest.skipUnless(
        shutil.which("Rscript") or shutil.which("Rscript.exe"),
        "Rscript is unavailable",
    )
    def test_limma_combined_disease_contrast_and_gene_alignment(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            genes = ["G1", "G2", "G3", "G4", "G5", "G6"]
            samples = ["S1", "S2", "S3", "S4", "S5", "S6"]
            counts = pd.DataFrame(
                np.column_stack(
                    [
                        np.full(6, 100.0),
                        np.full(6, 100.0),
                        np.full(6, 100.0),
                        np.full(6, 100.0),
                        np.full(6, 100.0),
                        np.full(6, 100.0),
                    ]
                ),
                index=genes,
                columns=samples,
            )
            counts.loc["G1", :] = [100, 100, 100, 100, 100, 100]
            counts.loc["G2", :] = [100, 100, 500, 500, 50, 50]
            counts.loc["G3", :] = [100, 100, 50, 50, 500, 500]
            counts.loc["G4", :] = [100, 100, 50, 50, 50, 50]
            counts.loc["G5", :] = [100, 100, 100, 100, 100, 100]
            counts.loc["G6", :] = [500, 500, 200, 200, 200, 200]
            counts.to_csv(root / "counts.csv")
            pd.DataFrame(
                {
                    "condition": [
                        "HC",
                        "HC",
                        "SS",
                        "SS",
                        "NASH",
                        "NASH",
                    ]
                },
                index=samples,
            ).to_csv(root / "metadata.csv")
            result = differential_expression_limma(
                root / "counts.csv",
                root / "metadata.csv",
                root / "result.csv",
                condition_column="condition",
                group1=["HC"],
                group2=["SS", "NASH"],
                comparison="HC_vs_NAFLD",
                input_type="counts",
                timeout=180,
            )
            by_gene = result.set_index("gene")
            self.assertGreater(by_gene.loc["G2", "delta"], 0)
            self.assertLess(by_gene.loc["G4", "delta"], 0)
            self.assertGreater(by_gene.loc["G2", "group2_mean"], 0)
            self.assertFalse(
                result[
                    ["gene", "group1_mean", "group2_mean", "delta"]
                ].isna().any().any()
            )

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

    def test_cache_requires_same_signature_and_content(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            output = root / "output.csv"
            state = root / "state.json"
            output.write_text("gene,value\nEGFR,1\n", encoding="utf-8")
            key = plan_cache.signature(
                files=[],
                parameters={"genes": ["EGFR"], "software": "test"},
            )
            plan_cache.save(state, key, [output])
            self.assertTrue(plan_cache.valid(state, key, [output]))
            changed_key = plan_cache.signature(
                files=[],
                parameters={"genes": ["TP53"], "software": "test"},
            )
            self.assertFalse(plan_cache.valid(state, changed_key, [output]))
            output.write_text("gene,value\nEGFR,2\n", encoding="utf-8")
            self.assertFalse(plan_cache.valid(state, key, [output]))

    def test_output_root_cannot_mix_run_ids(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            ExperimentPlanOne(root, {**default_config(), "run_id": "run-A"})
            with self.assertRaisesRegex(RuntimeError, "different run_id"):
                ExperimentPlanOne(
                    root,
                    {**default_config(), "run_id": "run-B"},
                )

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

    def test_static_coexpression_selection_is_rejected_as_cv_input(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            samples = [f"S{index}" for index in range(10)]
            pd.DataFrame(
                np.arange(30, dtype=float).reshape(3, 10),
                index=["G1", "G2", "G3"],
                columns=samples,
            ).to_csv(root / "expression.csv")
            pd.DataFrame(
                {"condition": ["HC", "NASH"] * 5},
                index=samples,
            ).to_csv(root / "metadata.csv")
            with self.assertRaisesRegex(
                ValueError,
                "coexpression",
            ):
                run_ml_validation(
                    root / "expression.csv",
                    root / "metadata.csv",
                    {"coexpression_disease_module": ["G1", "G2"]},
                    {},
                    root / "output",
                    cv_folds=5,
                )

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

    def test_empty_ppi_intersection_is_valid_negative_not_fabricated_network(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            output = Path(tmp) / "ppi"
            result = run_ppi_analysis([], output)
            self.assertEqual(result["status"], "valid_negative")
            summary = json.loads(
                (output / "ppi_summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(summary["status"], "valid_negative")
            self.assertEqual(summary["network_edges"], 0)

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
            (docking / "md_run_manifest.json").write_text(
                json.dumps(
                    {
                        "target_id": "GENE1",
                        "ligand_id": "lig1",
                        "run_id": "GENE1|lig1|md|run_001",
                        "run_dir": str(run),
                        "status": "completed",
                        "production_ns": 100.0,
                    }
                ),
                encoding="utf-8",
            )
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
                        f"{index * 50} {value}"
                        for index, value in enumerate(values)
                    ),
                    encoding="utf-8",
                )
            run.joinpath("rmsf_protein_residue.xvg").write_text(
                "1 0.1\n2 0.2\n",
                encoding="utf-8",
            )
            run.joinpath("FINAL_DECOMP_MMPBSA.dat").write_text(
                "#Residue TOTAL STD\n"
                "ALA1 -5.0 0.2\n"
                "GLY2 -4.0 0.1\n",
                encoding="utf-8",
            )
            result = generate_plan_md_figures(docking, root / "09_md_mmpbsa")
            self.assertEqual(result["status"], "completed")
            self.assertTrue(result["panels"]["d_rmsd"])
            self.assertTrue(result["panels"]["e_ligand_rmsd"])
            self.assertTrue(result["panels"]["f_rmsf"])
            self.assertTrue(result["panels"]["g_rg"])
            self.assertTrue(result["panels"]["h_mmpbsa"])

    def test_md_total_energy_alone_does_not_complete_mmpbsa_panel(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            docking = root / "08_docking"
            run = docking / "targets" / "GENE1" / "outputs" / "md" / "06_md" / "lig1"
            run.mkdir(parents=True)
            (docking / "md_run_manifest.json").write_text(
                json.dumps(
                    {
                        "target_id": "GENE1",
                        "run_id": "GENE1|lig1|md|run_001",
                        "run_dir": str(run),
                        "status": "completed",
                        "production_ns": 100.0,
                    }
                ),
                encoding="utf-8",
            )
            for name in (
                "rmsd_protein.xvg",
                "rmsd_ligand.xvg",
                "gyrate_protein.xvg",
            ):
                run.joinpath(name).write_text(
                    "0 0.1\n50 0.11\n100 0.12\n",
                    encoding="utf-8",
                )
            run.joinpath("rmsf_protein_residue.xvg").write_text(
                "1 0.1\n2 0.2\n",
                encoding="utf-8",
            )
            pd.DataFrame(
                {
                    "id": ["lig1"],
                    "mmpbsa_status": ["completed"],
                    "mmpbsa_delta_total_kj_mol": [-42.0],
                }
            ).to_csv(docking / "md_simulation_results.csv", index=False)
            result = generate_plan_md_figures(docking, root / "09_md_mmpbsa")
            self.assertFalse(result["panels"]["h_mmpbsa"])
            self.assertEqual(
                result["mmpbsa"]["status"],
                "incomplete_total_energy_only",
            )

    def test_cellchat_permutation_outputs_fdr(self):
        rng = np.random.default_rng(3)
        data = ad.AnnData(
            X=rng.poisson(1, size=(40, 4)).astype(float),
            obs=pd.DataFrame(
                {
                    "cell_type": ["A"] * 20 + ["B"] * 20,
                    "condition": (["NCD"] * 10 + ["HFD"] * 10) * 2,
                    "donor_id": (["D1"] * 5 + ["D2"] * 5 + ["D3"] * 5 + ["D4"] * 5)
                    * 2,
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
        self.assertEqual(result["status"]["status"], "completed")
        self.assertEqual(result["status"]["n_NCD_units"], 2)
        self.assertEqual(result["status"]["n_HFD_units"], 2)

    def test_cellchat_does_not_infer_significance_with_one_donor_per_group(self):
        data = ad.AnnData(
            X=np.ones((8, 4), dtype=float),
            obs=pd.DataFrame(
                {
                    "cell_type": ["A"] * 4 + ["B"] * 4,
                    "condition": ["NCD"] * 2 + ["HFD"] * 2 + ["NCD"] * 2 + ["HFD"] * 2,
                    "donor_id": ["N1"] * 2 + ["H1"] * 2 + ["N1"] * 2 + ["H1"] * 2,
                }
            ),
        )
        data.var_names = ["TNF", "TNFRSF1A", "IL6", "IL6R"]
        result = _cellchat_like_analysis(data, n_permutations=10, seed=1)
        self.assertEqual(result["status"]["status"], "descriptive_only")
        self.assertFalse(result["interactions"]["significant"].any())
        self.assertTrue(result["interactions"]["fdr"].isna().all())

    def test_cellchat_library_ids_are_not_treated_as_animal_replicates(self):
        data = ad.AnnData(
            X=np.ones((8, 4), dtype=float),
            obs=pd.DataFrame(
                {
                    "cell_type": ["A"] * 4 + ["B"] * 4,
                    "condition": ["NCD"] * 2 + ["HFD"] * 2 + ["NCD"] * 2 + ["HFD"] * 2,
                    "sample_id": ["L1"] * 2 + ["L2"] * 2 + ["L1"] * 2 + ["L2"] * 2,
                }
            ),
        )
        data.var_names = ["TNF", "TNFRSF1A", "IL6", "IL6R"]
        result = _cellchat_like_analysis(data, n_permutations=10, seed=1)
        self.assertEqual(result["status"]["status"], "descriptive_only")
        self.assertFalse(result["status"]["verified_biological_units"])
        self.assertIn("library IDs", result["status"]["reason"])

    def test_human_donor_expression_merges_technical_libraries(self):
        frame = pd.DataFrame(
            {
                "patient_id": ["P1", "P1", "P2"],
                "gene": ["G1", "G1", "G1"],
                "cell_type": ["A", "A", "A"],
                "condition": ["MASH", "MASH", "Healthy"],
                "expression": [1.0, 3.0, 2.0],
            }
        )
        merged = _donor_expression(frame, "patient_id")
        self.assertEqual(len(merged), 2)
        self.assertEqual(
            float(merged.loc[merged["patient_id"].eq("P1"), "expression"].iloc[0]),
            2.0,
        )
        self.assertEqual(
            _patient_id({"title": "PCL17-end stage-SITTE11"}),
            "PCL17",
        )

    def test_mouse_homolog_mapping_is_cached_with_source_version(self):
        def fake_get_json(url, **_kwargs):
            if "/query?" in url:
                return {
                    "hits": [
                        {
                            "symbol": "GPAT3",
                            "homologene": {
                                "id": 13099,
                                "genes": [
                                    [9606, 84803],
                                    [10090, 231510],
                                ],
                            },
                        }
                    ]
                }
            if "/gene/" in url:
                return {"symbol": "Gpat3"}
            raise AssertionError(url)

        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            output = Path(tmp) / "homologs.tsv"
            with mock.patch(
                "experiment_plan_one.single_cell.get_json",
                side_effect=fake_get_json,
            ):
                mapping = map_human_to_mouse_homologs(["GPAT3"], output)
            self.assertEqual(mapping, {"GPAT3": "GPAT3"})
            frame = pd.read_csv(output, sep="\t")
            self.assertEqual(frame.iloc[0]["mapping_status"], "one_to_one")
            self.assertEqual(frame.iloc[0]["source_version"], "HomoloGene:13099")

    def test_coexpression_module_analysis_outputs_hubs(self):
        rng = np.random.default_rng(3)
        sample_count = 24
        genes = [f"G{index:03d}" for index in range(90)]
        signals = [
            rng.normal(size=sample_count),
            rng.normal(size=sample_count),
            rng.normal(size=sample_count),
        ]
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            values = np.vstack(
                [
                    signal + rng.normal(scale=0.15, size=sample_count)
                    for signal in signals
                    for _ in range(30)
                ]
            )
            pd.DataFrame(
                values,
                index=genes,
                columns=[f"S{index}" for index in range(sample_count)],
            ).to_csv(root / "expression.csv")
            pd.DataFrame(
                {
                    "condition": ["HC"] * 12 + ["NASH"] * 12,
                },
                index=[f"S{index}" for index in range(sample_count)],
            ).to_csv(root / "metadata.csv")
            result = run_coexpression_analysis(
                root / "expression.csv",
                root / "metadata.csv",
                root / "coexpression",
                max_genes=90,
                min_module_size=10,
            )
            self.assertEqual(result["status"], "completed")
            self.assertGreaterEqual(result["n_modules"], 2)
            self.assertTrue(
                (
                    root
                    / "coexpression"
                    / "coexpression_hubs.csv"
                ).exists()
            )

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
            self.assertEqual(summary["performance_targets_evaluable"], 0)
            self.assertEqual(summary["performance_targets_unknown"], 1)
            self.assertEqual(summary["performance_targets_not_evaluated"], 2)

    def test_plan_coverage_marks_missing_tools_as_blocked_not_met(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            audit_dir = root / "10_reports" / "figure_quality_audit"
            audit_dir.mkdir(parents=True)
            pd.DataFrame(
                {
                    "figure": [
                        "Figure5_分子对接与分子动力学模拟",
                    ],
                    "panel": ["d"],
                    "status": ["prepared_not_run"],
                    "audit_verdict": ["不满足方案"],
                }
            ).to_csv(audit_dir / "figure_quality_audit.csv", index=False)
            summary = audit_plan_coverage(root)
            self.assertGreater(summary["blocked_panels"], 0)
            self.assertGreater(summary["not_run_panels"], 0)
            self.assertIn("blocked", summary["evidence_state_vocabulary"])
            self.assertIn("not_run", summary["evidence_state_vocabulary"])

    def test_ml_figure_reviews_follow_current_metrics(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp)
            ml_dir = root / "05_machine_learning"
            ml_dir.mkdir(parents=True)
            (ml_dir / "ml_summary.json").write_text(
                json.dumps(
                    {
                        "external_validation": [
                            {
                                "dataset": "GSE49541_fibrosis",
                                "auc": 0.309,
                                "target_auc": None,
                                "target_met": None,
                                "evaluate_auc_target": False,
                            },
                            {
                                "dataset": "GSE164441_tumor",
                                "auc": 0.86,
                                "target_auc": None,
                                "target_met": None,
                                "evaluate_auc_target": False,
                            },
                            {
                                "dataset": "GSE135251_NAFLD",
                                "auc": 0.909,
                                "target_auc": 0.8,
                                "target_met": True,
                            },
                        ],
                        "calibration": {
                            "hosmer_lemeshow_p": 0.1398,
                            "brier": 0.0806,
                            "calibration_slope": 1.5871,
                            "calibration_intercept": -0.3084,
                            "target_met": True,
                        },
                    }
                ),
                encoding="utf-8",
            )
            reviews = _dynamic_result_reviews(root)
            figure = "Figure3_机器学习模型构建与SHAP核心特征"
            self.assertEqual(reviews[(figure, "c")].verdict, "需限定解释")
            self.assertIn("0.309", reviews[(figure, "c")].notes)
            self.assertEqual(reviews[(figure, "d")].verdict, "需限定解释")
            self.assertIn("0.909", reviews[(figure, "d")].notes)
            self.assertEqual(reviews[(figure, "e")].verdict, "需限定解释")
            self.assertIn("0.140", reviews[(figure, "e")].notes)

    def test_figure_summary_counts_unique_panels_not_file_aliases(self):
        base = {
            "figure": "Figure1",
            "panel": "a",
            "content": "panel a",
            "status": "available",
            "automatic_quality_score": 95,
            "dpi_x": 600,
            "dpi_y": 600,
            "physical_width_in": 5.0,
            "has_pdf": True,
            "has_svg": True,
            "recommended_output": "主图候选",
            "label_overlap_severity": "无",
            "label_overlap_notes": "",
            "duplicate_svg_text_anchor_groups": 0,
            "audit_verdict": "可用",
            "notes": "",
        }
        audit = pd.DataFrame(
            [
                base,
                {**base, "status": "missing"},
                {
                    **base,
                    "panel": "b",
                    "status": "missing",
                    "automatic_quality_score": None,
                },
            ]
        )
        summary = _summary(audit)
        self.assertEqual(summary["panel_entries"], 2)
        self.assertEqual(summary["image_entries"], 3)
        self.assertEqual(summary["available_panels"], 1)
        self.assertEqual(summary["missing_panels"], 1)


if __name__ == "__main__":
    unittest.main()
