"""Tests for pseudobulk omics quality control."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.omics_qc import run_omics_qc


class TestOmicsQc(unittest.TestCase):
    def test_reports_matrix_and_metadata_gates(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            data = workdir / "data" / "knockout"
            data.mkdir(parents=True)
            pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2", "GENE3"],
                    "T1": [5.0, 4.0, 0.0],
                    "T2": [5.5, 4.5, 0.0],
                    "N1": [2.0, 3.0, 0.0],
                    "N2": [2.5, 3.5, 0.0],
                }
            ).to_csv(data / "expression.csv", index=False)
            pd.DataFrame(
                {
                    "sample": ["T1", "T2", "N1", "N2"],
                    "condition": ["Tumor", "Tumor", "Normal", "Normal"],
                    "cell_type": ["Hepatocyte"] * 4,
                }
            ).to_csv(data / "metadata.csv", index=False)
            summary = run_omics_qc(workdir, workdir / "outputs" / "integration")
            self.assertEqual(summary["n_genes"], 3)
            self.assertEqual(summary["n_samples"], 4)
            self.assertTrue(summary["gate_passed"])
            self.assertTrue(
                (
                    workdir
                    / "outputs"
                    / "integration"
                    / "omics_qc_summary.json"
                ).exists()
            )


if __name__ == "__main__":
    unittest.main()
