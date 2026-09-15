#!/usr/bin/env python3
"""Tests for local MR harmonisation and sensitivity analysis."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from analysis.mr_coloc import _ivw, harmonise, run  # noqa: E402


class TestMRColoc(unittest.TestCase):
    def test_harmonise_flips_outcome_alleles(self):
        exposure = pd.DataFrame(
            {
                "snp": ["rs1", "rs2"],
                "beta": [0.2, 0.3],
                "se": [0.02, 0.03],
                "effect_allele": ["A", "C"],
                "other_allele": ["G", "T"],
                "pval": [1e-8, 1e-8],
            }
        )
        outcome = pd.DataFrame(
            {
                "snp": ["rs1", "rs2"],
                "beta": [0.1, 0.15],
                "se": [0.02, 0.03],
                "effect_allele": ["G", "T"],
                "other_allele": ["A", "C"],
                "pval": [0.01, 0.01],
            }
        )
        merged = harmonise(exposure, outcome)
        self.assertEqual(len(merged), 2)
        self.assertTrue((merged["allele_status"] == "flipped").all())
        self.assertLess(float(merged["ratio"].mean()), 0.0)

    def test_palindromic_snp_without_eaf_is_removed(self):
        exposure = pd.DataFrame(
            {
                "snp": ["rs1"],
                "beta": [0.2],
                "se": [0.02],
                "effect_allele": ["A"],
                "other_allele": ["T"],
                "pval": [1e-8],
            }
        )
        outcome = exposure.copy()
        merged = harmonise(exposure, outcome)
        self.assertTrue(merged.empty)
        self.assertEqual(merged.attrs["missing_eaf_palindromic_snps"], 1)

    def test_distance_clumping_removes_correlated_neighbor(self):
        exposure = pd.DataFrame(
            {
                "snp": ["rs1", "rs2", "rs3"],
                "beta": [0.2, 0.2, 0.2],
                "se": [0.02, 0.02, 0.02],
                "effect_allele": ["A", "A", "A"],
                "other_allele": ["G", "G", "G"],
                "pval": [1e-8, 2e-8, 1e-9],
                "chromosome": ["1", "1", "2"],
                "position": [1000, 2000, 1000],
            }
        )
        # Exercise the public run path with a 10 kb distance.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exposure.to_csv(root / "exposure.csv", index=False)
            pd.DataFrame(
                {
                    **{
                        column: exposure[column]
                        for column in exposure.columns
                    },
                    "pval": [0.2, 0.2, 0.2],
                }
            ).to_csv(root / "outcome.csv", index=False)
            config = {
                "_config_dir": str(root),
                "exposure": {
                    "file": "exposure.csv",
                    "snp": "snp",
                    "beta": "beta",
                    "se": "se",
                    "effect_allele": "effect_allele",
                    "other_allele": "other_allele",
                    "pval": "pval",
                    "chromosome": "chromosome",
                    "position": "position",
                },
                "outcome": {
                    "file": "outcome.csv",
                    "snp": "snp",
                    "beta": "beta",
                    "se": "se",
                    "effect_allele": "effect_allele",
                    "other_allele": "other_allele",
                    "pval": "pval",
                },
                "p_threshold": 1.0,
                "clump": {
                    "enabled": True,
                    "distance_kb": 10,
                },
                "coloc": {"enabled": False},
            }
            summary = run(config, root / "out", skip_r=True)
            self.assertEqual(summary["n_harmonised_snps"], 2)
            self.assertEqual(
                summary["clumping"]["method"],
                "distance_pruning_no_ld_reference",
            )
            self.assertEqual(summary["clumping"]["removed"], 1)

    def test_ivw_and_full_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rng = np.random.default_rng(4)
            snps = [f"rs{i}" for i in range(20)]
            bx = rng.normal(0.2, 0.05, len(snps))
            by = bx * 0.7 + rng.normal(0, 0.02, len(snps))
            exposure = pd.DataFrame(
                {
                    "SNP": snps,
                    "BETA": bx,
                    "SE": rng.uniform(0.02, 0.04, len(snps)),
                    "EA": "A",
                    "OA": "G",
                    "P": 1e-8,
                }
            )
            outcome = pd.DataFrame(
                {
                    "SNP": snps,
                    "BETA": by,
                    "SE": rng.uniform(0.02, 0.04, len(snps)),
                    "EA": "A",
                    "OA": "G",
                    "P": 0.2,
                }
            )
            exposure.to_csv(root / "exposure.csv", index=False)
            outcome.to_csv(root / "outcome.csv", index=False)
            config = {
                "_config_dir": str(root),
                "exposure": {
                    "file": "exposure.csv",
                    "snp": "SNP",
                    "beta": "BETA",
                    "se": "SE",
                    "effect_allele": "EA",
                    "other_allele": "OA",
                    "pval": "P",
                },
                "outcome": {
                    "file": "outcome.csv",
                    "snp": "SNP",
                    "beta": "BETA",
                    "se": "SE",
                    "effect_allele": "EA",
                    "other_allele": "OA",
                    "pval": "P",
                },
                "p_threshold": 1.0,
                "coloc": {"enabled": False},
            }
            summary = run(config, root / "out", skip_r=True)
            self.assertEqual(summary["status"], "completed")
            self.assertAlmostEqual(
                summary["primary"]["estimate"],
                0.7,
                delta=0.15,
            )
            self.assertTrue((root / "out" / "mr_methods.csv").exists())

    def test_ivw_estimator(self):
        frame = pd.DataFrame(
            {
                "beta_exposure": [0.1, 0.2],
                "beta_outcome": [0.05, 0.1],
                "se_outcome": [0.1, 0.1],
                "ratio": [0.5, 0.5],
            }
        )
        result = _ivw(frame)
        self.assertAlmostEqual(result["estimate"], 0.5, places=6)


if __name__ == "__main__":
    unittest.main()
