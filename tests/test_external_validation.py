"""Tests for external target-ranking validation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.external_validation import run_external_validation


class TestExternalValidation(unittest.TestCase):
    def test_computes_discrimination_metrics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "validation.csv"
            pd.DataFrame(
                {
                    "gene": ["A", "B", "C", "D"],
                    "label": [1, 1, 0, 0],
                }
            ).to_csv(table, index=False)
            priority = root / "priority.csv"
            pd.DataFrame(
                {
                    "gene": ["A", "B", "C", "D"],
                    "integrated_score": [0.95, 0.80, 0.30, 0.10],
                }
            ).to_csv(priority, index=False)
            summary = run_external_validation(
                table,
                root / "out",
                pipeline_scores_path=priority,
                bootstrap=50,
            )
            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["auroc"], 1.0)
            self.assertEqual(summary["auprc"], 1.0)
            self.assertEqual(summary["score_origin"], "pipeline")
            self.assertTrue(summary["score_provenance_valid"])
            self.assertTrue(
                (root / "out" / "external_validation_predictions.csv").exists()
            )

    def test_missing_table_is_reported_as_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary = run_external_validation(
                Path(tmp) / "missing.csv",
                Path(tmp) / "out",
            )
            self.assertEqual(summary["status"], "skipped")

    def test_external_score_requires_explicit_opt_in(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "validation.csv"
            pd.DataFrame(
                {
                    "gene": ["A", "B"],
                    "score": [0.9, 0.1],
                    "label": [1, 0],
                }
            ).to_csv(table, index=False)
            summary = run_external_validation(table, root / "out")
            self.assertEqual(summary["status"], "skipped")
            self.assertIn("allow_external_score", summary["reason"])

    def test_explicit_external_score_is_not_pipeline_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "validation.csv"
            pd.DataFrame(
                {
                    "gene": ["A", "B"],
                    "score": [0.9, 0.1],
                    "label": [1, 0],
                }
            ).to_csv(table, index=False)
            summary = run_external_validation(
                table,
                root / "out",
                allow_external_score=True,
                bootstrap=20,
            )
            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["score_origin"], "external_allowed")
            self.assertFalse(summary["score_provenance_valid"])

    def test_reports_unmatched_pipeline_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            table = root / "validation.csv"
            pd.DataFrame(
                {
                    "gene": ["A", "B", "UNKNOWN"],
                    "label": [1, 0, 1],
                }
            ).to_csv(table, index=False)
            priority = root / "priority.csv"
            pd.DataFrame(
                {
                    "gene": ["A", "B"],
                    "integrated_score": [0.9, 0.1],
                }
            ).to_csv(priority, index=False)
            summary = run_external_validation(
                table,
                root / "out",
                pipeline_scores_path=priority,
                bootstrap=20,
            )
            self.assertEqual(summary["target_match_count"], 2)
            self.assertAlmostEqual(summary["target_match_rate"], 2 / 3)
            self.assertIn("UNKNOWN", summary["unmatched_targets"])
            self.assertIn("f1", summary)
            self.assertIn("specificity", summary)


if __name__ == "__main__":
    unittest.main()
