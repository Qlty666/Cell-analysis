#!/usr/bin/env python3
"""Tests for the single-cell in-silico knockout module."""

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
from docking.insilico import run_insilico_knockout  # noqa: E402

DEFAULT_CONFIG = APP_ROOT / "config" / "docking_config.json"
LOG = logging.getLogger("test_insilico_knockout")


class TestInSilicoKnockout(unittest.TestCase):
    def _write_inputs(self, workdir: Path) -> Path:
        data_dir = workdir / "data" / "knockout"
        data_dir.mkdir(parents=True, exist_ok=True)
        rng = np.random.default_rng(42)
        n_a = 160
        n_b = 160
        cells = [f"A{i}" for i in range(n_a)] + [f"B{i}" for i in range(n_b)]
        driver_a = rng.normal(2.0, 0.4, n_a)
        driver_b = rng.normal(0.3, 0.3, n_b)
        driver = np.concatenate([driver_a, driver_b])
        target_genes = [f"T{i:02d}" for i in range(12)]
        genes = ["TF1", "TF2", *target_genes, "N1", "N2", "N3"]
        base = rng.normal(0.0, 0.25, size=(len(genes), len(cells)))
        frame = pd.DataFrame(base, index=genes, columns=cells)
        frame.loc["TF1"] = driver
        for gene in target_genes:
            frame.loc[gene] = driver * (0.6 + 0.2 * rng.random()) + rng.normal(
                0, 0.12, len(cells)
            )
        frame = np.log1p(np.clip(frame, 0, None))
        frame = frame.reset_index().rename(columns={"index": "gene"})
        frame.to_csv(data_dir / "expression.csv", index=False)

        pd.DataFrame(
            {
                "cell": cells,
                "cell_type": ["Erythroid"] * n_a + ["Myeloid"] * n_b,
                "umap_1": np.concatenate(
                    [rng.normal(-3, 1.0, n_a), rng.normal(3, 1.0, n_b)]
                ),
                "umap_2": np.concatenate(
                    [rng.normal(1, 1.0, n_a), rng.normal(-2, 1.0, n_b)]
                ),
            }
        ).to_csv(data_dir / "metadata.csv", index=False)
        pd.DataFrame({"regulator": ["TF1", "TF2"]}).to_csv(
            data_dir / "regulators.csv", index=False
        )
        return data_dir

    def test_run_generates_report_and_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            data_dir = self._write_inputs(workdir)
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "expression_csv": "data/knockout/expression.csv",
                    "metadata_csv": "data/knockout/metadata.csv",
                    "insilico_gene": "TF1",
                    "insilico_regulators_csv": "data/knockout/regulators.csv",
                },
            )
            cfg.data["insilico_knockout"].update(
                {
                    "enabled": True,
                    "run_enrichment": False,
                    "max_genes": 50,
                    "min_cells": 10,
                    "n_propagation": 3,
                    "network_edges_per_regulator": 20,
                }
            )
            summary = run_insilico_knockout(cfg, LOG, ko_gene="TF1")
            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["ko_gene"], "TF1")
            in_silico = cfg.knockout_dir() / "in_silico"
            self.assertTrue((in_silico / "insilico_summary.json").exists())
            self.assertTrue((in_silico / "in_silico_knockout_report.html").exists())
            self.assertTrue(
                (in_silico / "data" / "insilico_target_changes.csv").exists()
            )
            self.assertTrue(
                (in_silico / "data" / "fig_63_ko_target_top15.csv").exists()
            )
            for figure in summary["figures"]:
                self.assertTrue((in_silico / "figures" / figure).exists())
            changes = pd.read_csv(
                in_silico / "data" / "insilico_target_changes.csv"
            )
            self.assertIn("TF1", changes["gene"].head(5).tolist())
            row = changes[changes["gene"] == "TF1"].iloc[0]
            self.assertLess(row["ko_mean"], row["wt_mean"])
            self.assertTrue((data_dir / "expression.csv").exists())


if __name__ == "__main__":
    unittest.main()
