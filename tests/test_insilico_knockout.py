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
from docking.insilico import _looks_like_raw_counts  # noqa: E402
from docking.insilico import _merge_scTenifold  # noqa: E402
from docking.insilico import _scTenifold_available  # noqa: E402
from docking.insilico import _plot_enrichment_bubble  # noqa: E402
from docking.insilico import _plot_umap_shift  # noqa: E402

DEFAULT_CONFIG = APP_ROOT / "config" / "docking_config.json"
LOG = logging.getLogger("test_insilico_knockout")


class TestInSilicoKnockout(unittest.TestCase):
    def test_go_bubble_splits_into_bp_cc_mf(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows = []
            for ont, base_p in (("BP", 1e-5), ("CC", 2e-4), ("MF", 3e-4)):
                for idx in range(8):
                    p = base_p * (idx + 1)
                    rows.append(
                        {
                            "ONTOLOGY": ont,
                            "Description": f"{ont} term {idx}",
                            "GeneRatio": f"{idx + 2}/100",
                            "Count": idx + 2,
                            "pvalue": p,
                            "p.adjust": p,
                        }
                    )
            csv_path = Path(tmp) / "go.csv"
            out_path = Path(tmp) / "fig_67_ko_go_enrichment.png"
            pd.DataFrame(rows).to_csv(csv_path, index=False)
            ok = _plot_enrichment_bubble(
                csv_path,
                out_path,
                "GO Enrichment",
                group_col="ONTOLOGY",
            )
            self.assertTrue(ok)
            self.assertGreater(out_path.stat().st_size, 10_000)

    def test_dense_umap_arrows_use_grid_without_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            n = 1200
            rng = np.random.default_rng(5)
            embedding = pd.DataFrame(
                {
                    "umap_1": rng.normal(0, 4, n),
                    "umap_2": rng.normal(0, 3, n),
                },
                index=[f"C{i}" for i in range(n)],
            )
            shift = pd.DataFrame(
                {
                    "shift_1": rng.normal(0, 0.5, n),
                    "shift_2": rng.normal(0, 0.5, n),
                },
                index=embedding.index,
            )
            cell_types = pd.Series(
                rng.choice(["A", "B", "C"], n), index=embedding.index
            )
            out = Path(tmp) / "fig_66_dense.png"
            _plot_umap_shift(embedding, shift, cell_types, "TF1", out)
            self.assertTrue(out.exists())
            self.assertGreater(out.stat().st_size, 10_000)

    def _write_inputs(self, workdir: Path, raw: bool = False) -> Path:
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
        if not raw:
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

    def test_looks_like_raw_counts_detects_count_matrix(self):
        counts = pd.DataFrame(
            {
                "C1": [0, 4, 100, 5],
                "C2": [1, 0, 60, 3],
                "C3": [2, 7, 120, 0],
            },
            index=["A", "B", "C", "D"],
        )
        self.assertTrue(_looks_like_raw_counts(counts))
        logged = np.log1p(counts)
        self.assertFalse(_looks_like_raw_counts(logged))

    def test_merge_scTenifold_adds_combined_rank(self):
        changes = pd.DataFrame(
            {
                "gene": ["A", "B", "C", "D"],
                "wt_mean": [1.0, 1.0, 1.0, 1.0],
                "ko_mean": [0.7, 0.8, 0.9, 1.0],
                "delta": [-0.3, -0.2, -0.1, 0.0],
                "abs_delta": [0.3, 0.2, 0.1, 0.0],
            }
        )
        diff = pd.DataFrame(
            {
                "gene": ["D", "A", "B", "C"],
                "sctenifold_distance": [0.5, 0.2, 0.1, 0.05],
                "sctenifold_padj": [0.001, 0.01, 0.02, 0.1],
            }
        )
        merged, kept = _merge_scTenifold(changes, diff)
        self.assertIsNotNone(kept)
        self.assertIn("sctenifold_score", merged.columns)
        self.assertIn("combined_impact", merged.columns)
        self.assertEqual(merged["gene"].iloc[0], "A")

    @unittest.skipUnless(
        _scTenifold_available(),
        "scTenifoldpy not installed",
    )
    def test_run_uses_real_scTenifold_engine_on_raw_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            data_dir = self._write_inputs(workdir, raw=True)
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "expression_csv": "data/knockout/expression.csv",
                    "metadata_csv": "data/knockout/metadata.csv",
                    "insilico_gene": "TF1",
                    "insilico_engine": "triple",
                    "insilico_raw_count_input": True,
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
                    "scTenifold_min_lib_size": 1,
                    "scTenifold_min_percent": 0.01,
                    "scTenifold_min_exp_sum": 1,
                    "scTenifold_remove_outlier_cells": False,
                    "scTenifold_n_networks": 1,
                    "scTenifold_n_cells": 50,
                    "scTenifold_q": 0.5,
                    "scTenifold_K": 2,
                }
            )
            summary = run_insilico_knockout(cfg, LOG, ko_gene="TF1")
            self.assertEqual(summary["status"], "completed")
            self.assertIn("scTenifoldKnk", summary["engine"])
            self.assertIsNotNone(summary["scTenifoldKnk"])
            self.assertEqual(summary["drugreflector"]["status"], "skipped")
            in_silico = cfg.knockout_dir() / "in_silico"
            sc_csv = in_silico / "data" / "insilico_scTenifold_results.csv"
            self.assertTrue(sc_csv.exists())
            changes = pd.read_csv(
                in_silico / "data" / "insilico_target_changes.csv"
            )
            self.assertIn("sctenifold_padj", changes.columns)
            self.assertIn("combined_impact", changes.columns)
            self.assertTrue((data_dir / "expression.csv").exists())

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
