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
            self.assertTrue(
                (
                    workdir
                    / "outputs"
                    / "integration"
                    / "omics_qc_pca.csv"
                ).exists()
            )
            self.assertTrue(
                (
                    workdir
                    / "outputs"
                    / "integration"
                    / "omics_qc_correlation.csv"
                ).exists()
            )

    def test_zero_fraction_uses_only_finite_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            data = workdir / "data" / "knockout"
            data.mkdir(parents=True)
            pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2"],
                    "T1": [0.0, float("nan")],
                    "T2": [0.0, float("nan")],
                    "N1": [5.0, float("nan")],
                    "N2": [6.0, float("nan")],
                }
            ).to_csv(data / "expression.csv", index=False)
            pd.DataFrame(
                {
                    "sample": ["T1", "T2", "N1", "N2"],
                    "condition": ["Tumor", "Tumor", "Normal", "Normal"],
                }
            ).to_csv(data / "metadata.csv", index=False)
            summary = run_omics_qc(
                workdir,
                workdir / "outputs" / "integration",
            )
            self.assertEqual(summary["zero_fraction"], 0.5)
            sample_metrics = pd.read_csv(
                workdir
                / "outputs"
                / "integration"
                / "omics_qc_sample_metrics.csv"
            )
            self.assertTrue(
                sample_metrics["zero_fraction"].notna().all()
            )

    def test_reports_batch_confounding(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            data = workdir / "data" / "knockout"
            data.mkdir(parents=True)
            pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2", "GENE3"],
                    "T1": [8.0, 7.0, 6.0],
                    "T2": [8.5, 7.5, 6.5],
                    "N1": [2.0, 3.0, 4.0],
                    "N2": [2.5, 3.5, 4.5],
                }
            ).to_csv(data / "expression.csv", index=False)
            pd.DataFrame(
                {
                    "sample": ["T1", "T2", "N1", "N2"],
                    "condition": ["Tumor", "Tumor", "Normal", "Normal"],
                    "batch": ["B1", "B1", "B2", "B2"],
                }
            ).to_csv(data / "metadata.csv", index=False)
            summary = run_omics_qc(
                workdir,
                workdir / "outputs" / "integration",
            )
            self.assertTrue(summary["has_batch_metadata"])
            self.assertTrue(summary["batch_condition_confounding"])
            self.assertFalse(
                summary["gates"]["batch_confounding_acceptable"]
            )


if __name__ == "__main__":
    unittest.main()
