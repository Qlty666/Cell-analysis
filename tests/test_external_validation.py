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
                    "score": [0.95, 0.80, 0.30, 0.10],
                    "label": [1, 1, 0, 0],
                }
            ).to_csv(table, index=False)
            summary = run_external_validation(
                table,
                root / "out",
                bootstrap=50,
            )
            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["auroc"], 1.0)
            self.assertEqual(summary["auprc"], 1.0)
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


if __name__ == "__main__":
    unittest.main()
