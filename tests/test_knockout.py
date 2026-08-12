#!/usr/bin/env python3
"""Unit tests for the virtual knockout scoring module."""

from __future__ import annotations

import logging
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from docking.config import load_config  # noqa: E402
from docking.knockout import run_knockout  # noqa: E402

DEFAULT_CONFIG = APP_ROOT / "config" / "docking_config.json"
LOG = logging.getLogger("test_knockout")


class TestKnockoutConfig(unittest.TestCase):
    def test_defaults_and_overrides(self):
        cfg = load_config(DEFAULT_CONFIG)
        self.assertTrue(cfg.get("knockout", "enabled", False))
        self.assertEqual(cfg.get("knockout", "top_n", None), 50)

        cfg = load_config(
            DEFAULT_CONFIG,
            {
                "expression_csv": "data/ko/expression.csv",
                "metadata_csv": "data/ko/metadata.csv",
                "depmap_csv": "data/ko/depmap.csv",
                "case_label": "Tumor",
                "normal_label": "Normal",
                "ko_top_n": 10,
            },
        )
        self.assertEqual(
            cfg.get("knockout", "expression_csv"),
            "data/ko/expression.csv",
        )
        self.assertEqual(cfg.get("knockout", "top_n"), 10)
        self.assertEqual(cfg.get("knockout", "case_label"), "Tumor")


class TestRunKnockout(unittest.TestCase):
    def _make_inputs(self, workdir: Path, n_genes: int = 40):
        data_dir = workdir / "data" / "knockout"
        data_dir.mkdir(parents=True, exist_ok=True)
        rng = np.random.default_rng(7)
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
        genes = genes[:n_genes]
        rows = {"gene": genes}
        for sample in ["T1", "T2", "N1", "N2"]:
            base = 3.0 if sample.startswith("T") else 1.0
            rows[sample] = base + rng.normal(0, 0.15, len(genes))
        expr = pd.DataFrame(rows)
        expr.loc[expr["gene"].isin(signature), ["T1", "T2"]] += 2.5
        expr.to_csv(data_dir / "expression.csv", index=False)

        pd.DataFrame(
            {
                "sample": ["T1", "T2", "N1", "N2"],
                "condition": ["Tumor", "Tumor", "Normal", "Normal"],
            }
        ).to_csv(data_dir / "metadata.csv", index=False)

        pd.DataFrame(
            {
                "gene": genes,
                "effect": [
                    -1.1 if gene in signature else 0.3 for gene in genes
                ],
            }
        ).to_csv(data_dir / "depmap.csv", index=False)

    def test_wide_matrix_with_depmap(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            self._make_inputs(workdir)
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "expression_csv": "data/knockout/expression.csv",
                    "metadata_csv": "data/knockout/metadata.csv",
                    "depmap_csv": "data/knockout/depmap.csv",
                    "case_label": "Tumor",
                    "normal_label": "Normal",
                    "ko_top_n": 10,
                },
            )
            summary = run_knockout(cfg, LOG)
            ko_dir = cfg.knockout_dir()
            ranked = ko_dir / "data" / "fig_52_53_ranked_knockout.csv"
            self.assertTrue(ranked.exists())
            self.assertTrue((ko_dir / "summary.json").exists())
            self.assertTrue((ko_dir / "data" / "knockout_report.md").exists())

            frame = pd.read_csv(ranked)
            self.assertEqual(len(frame), 40)
            self.assertGreaterEqual(frame["knockout_score"].min(), 0.0)
            self.assertLessEqual(frame["knockout_score"].max(), 1.0)
            self.assertIn("TOP2A", frame["gene"].head(10).tolist())
            self.assertEqual(summary["case_samples"], 2)
            self.assertEqual(summary["normal_samples"], 2)
            self.assertEqual(summary["depmap_lines"], 40)
            self.assertTrue(summary["figures"])

    def test_long_format_without_depmap(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            data_dir = workdir / "data" / "knockout"
            data_dir.mkdir(parents=True)
            rows = []
            genes = ["ALB", "GPC3", "AFP", "CDK1", "TOP2A", "MKI67", "X1", "X2"]
            for i, gene in enumerate(genes):
                for sample, group in [
                    ("S1", "Tumor"),
                    ("S2", "Tumor"),
                    ("S3", "Normal"),
                    ("S4", "Normal"),
                ]:
                    rows.append(
                        {
                            "gene": gene,
                            "sample": sample,
                            "value": 1.0 + (i % 4) * 0.5,
                            "condition": group,
                        }
                    )
            pd.DataFrame(rows).to_csv(
                data_dir / "expression_long.csv",
                index=False,
            )
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "expression_csv": "data/knockout/expression_long.csv",
                    "case_label": "Tumor",
                    "normal_label": "Normal",
                },
            )
            summary = run_knockout(cfg, LOG)
            ranked = cfg.knockout_dir() / "data" / "fig_52_53_ranked_knockout.csv"
            frame = pd.read_csv(ranked)
            self.assertEqual(len(frame), len(genes))
            self.assertFalse(summary["depmap_included"])
            self.assertIn("log2fc", frame.columns)


if __name__ == "__main__":
    unittest.main()
