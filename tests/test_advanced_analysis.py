#!/usr/bin/env python3
"""Tests for the multi-cohort advanced bulk analysis helpers."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from analysis.advanced_analysis import (  # noqa: E402
    _evaluate_models,
    _load_optional_evidence,
    differential_expression,
    integrate_priorities,
    normalize_expression,
    run,
)


class TestAdvancedBulkHelpers(unittest.TestCase):
    def test_count_detection_and_log_cpm(self):
        frame = pd.DataFrame(
            {
                "S1": [10, 20, 30],
                "S2": [20, 40, 60],
            },
            index=["A", "B", "C"],
        )
        normalized, kind = normalize_expression(frame)
        self.assertEqual(kind, "counts")
        self.assertTrue(np.isfinite(normalized.to_numpy()).all())
        self.assertGreater(float(normalized.iloc[0, 0]), 0.0)

    def test_differential_expression_and_priority_integration(self):
        rng = np.random.default_rng(11)
        expression = pd.DataFrame(
            rng.normal(size=(6, 10)),
            index=[f"G{i}" for i in range(6)],
            columns=[f"S{i}" for i in range(10)],
        )
        expression.loc["G0", "S5":] += 3.0
        labels = pd.Series(
            ["Control"] * 5 + ["Disease"] * 5,
            index=expression.columns,
        )
        deg = differential_expression(
            expression,
            labels,
            "Disease",
            "Control",
        )
        self.assertIn("padj", deg.columns)
        self.assertEqual(deg.iloc[0]["gene"], "G0")
        importance = pd.DataFrame(
            {"importance": [0.9, 0.1]},
            index=["G0", "G1"],
        )
        integrated = integrate_priorities(
            deg,
            None,
            importance,
            None,
            {"weights": {"deg_score": 1.0, "ml_score": 1.0}},
        )
        self.assertIn("priority_score", integrated.columns)
        self.assertEqual(integrated.iloc[0]["gene"], "G0")

    def test_run_without_optional_r(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(5)
            genes = [f"G{i}" for i in range(12)]
            samples = [f"S{i}" for i in range(12)]
            matrix = rng.integers(20, 100, size=(12, 12)).astype(float)
            matrix[:3, 6:] += 100
            pd.DataFrame(matrix, index=genes, columns=samples).rename_axis(
                "gene"
            ).to_csv(root / "expression.csv")
            pd.DataFrame(
                {
                    "sample": samples,
                    "condition": ["Control"] * 6 + ["Disease"] * 6,
                }
            ).to_csv(root / "metadata.csv", index=False)
            args = Namespace(
                config=None,
                expression=str(root / "expression.csv"),
                metadata=str(root / "metadata.csv"),
                cohort_name="discovery",
                condition_column="condition",
                gene_column=None,
                case_label="Disease",
                control_label="Control",
                skip_r=True,
                skip_ml=False,
                skip_wgcna=True,
                skip_immune=True,
                skip_survival=True,
                feature_cap=8,
                cv_folds=2,
                cv_repeats=1,
                seed=7,
                output=str(root / "out"),
            )
            self.assertEqual(run(args), 0)
            self.assertTrue((root / "out" / "ml_model_comparison.csv").exists())
            self.assertTrue((root / "out" / "integrated_priority.csv").exists())

    def test_case_label_sort_order_does_not_invert_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(12)
            genes = [f"G{i}" for i in range(8)]
            samples = [f"S{i}" for i in range(24)]
            matrix = rng.normal(size=(8, 24))
            matrix[:4, 12:] += 3.0
            expression = pd.DataFrame(
                matrix,
                index=genes,
                columns=samples,
            )
            labels = pd.Series(
                ["A"] * 12 + ["B"] * 12,
                index=samples,
            )
            comparison, _, _ = _evaluate_models(
                expression,
                labels,
                genes,
                {
                    "case_label": "A",
                    "cv_folds": 3,
                    "cv_repeats": 1,
                    "seed": 7,
                    "feature_cap": 8,
                },
                root,
            )
            self.assertGreater(
                float(comparison.iloc[0]["cv_auc_mean"]),
                0.8,
            )
            summary = json.loads(
                (root / "ml_summary.json").read_text(encoding="utf-8")
            )
            self.assertGreater(float(summary["cv_auc"]), 0.8)
            self.assertEqual(summary["positive_class"], "A")

    def test_explicit_evidence_score_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = root / "ppi_hub_scores.csv"
            pd.DataFrame(
                {
                    "gene": ["A", "B"],
                    "ppi_degree": [99, 1],
                    "ppi_hub_score": [0.2, 0.9],
                }
            ).to_csv(evidence, index=False)
            frame = _load_optional_evidence(
                {
                    "network_hub_score": {
                        "path": str(evidence),
                        "score_column": "ppi_hub_score",
                    }
                },
                root,
                ["A", "B"],
            )
            self.assertAlmostEqual(frame.loc["A", "network_hub_score"], 0.2)
            self.assertAlmostEqual(frame.loc["B", "network_hub_score"], 0.9)


if __name__ == "__main__":
    unittest.main()
