#!/usr/bin/env python3
"""Tests for compound-disease network toxicology and PPI hub scoring."""

from __future__ import annotations

import logging
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from docking.config import load_config  # noqa: E402
from docking.knockout import run_knockout  # noqa: E402
from docking.network_toxicology import (  # noqa: E402
    overlap_analysis,
    ppi_hub_scores,
    run_network_toxicology,
)

DEFAULT_CONFIG = APP_ROOT / "config" / "docking_config.json"
LOG = logging.getLogger("test_network_toxicology")


class TestNetworkToxicology(unittest.TestCase):
    def _write_inputs(self, workdir: Path):
        data_dir = workdir / "data" / "network"
        data_dir.mkdir(parents=True)
        pd.DataFrame(
            {
                "gene": ["ALB", "GPC3", "MMP9", "AKT1", "SRC"],
                "source": ["CTD", "CTD", "ChEMBL", "STITCH", "CTD"],
            }
        ).to_csv(data_dir / "targets.csv", index=False)
        pd.DataFrame(
            {"gene": ["ALB", "GPC3", "MMP9", "EGFR", "TP53"]}
        ).to_csv(data_dir / "disease.csv", index=False)
        pd.DataFrame(
            {
                "protein1": ["ALB", "GPC3", "MMP9", "SRC"],
                "protein2": ["GPC3", "MMP9", "AKT1", "EGFR"],
            }
        ).to_csv(data_dir / "ppi.tsv", sep="\t", index=False)
        return data_dir

    def test_overlap_and_full_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            data_dir = self._write_inputs(workdir)
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "compound_name": "Bisphenol A",
                    "disease_name": "Liver Cancer",
                    "compound_targets_csv": str(data_dir / "targets.csv"),
                    "disease_genes_csv": str(data_dir / "disease.csv"),
                    "ppi_network_csv": str(data_dir / "ppi.tsv"),
                    "network_output_dir": "outputs/network",
                },
            )
            summary = run_network_toxicology(cfg, LOG)
            out_dir = workdir / "outputs" / "network"
            overlap = pd.read_csv(out_dir / "data" / "compound_disease_overlap.csv")
            self.assertEqual(set(overlap["gene"]), {"ALB", "GPC3", "MMP9"})
            self.assertTrue((out_dir / "figures" / "compound_disease_venn.png").exists())
            self.assertTrue((out_dir / "data" / "ppi_hub_scores.csv").exists())
            self.assertTrue((out_dir / "data" / "ctpd_nodes.csv").exists())
            self.assertTrue((out_dir / "data" / "ctpd_edges.csv").exists())
            self.assertTrue((out_dir / "data" / "ctpd_network.html").exists())
            self.assertEqual(summary["overlap_genes"], 3)
            self.assertTrue(summary["ppi_hub_scored"])

    def test_overlap_analysis_counts_sources(self):
        targets = {
            "ctd": pd.DataFrame({"gene": ["ALB", "MMP9"]}),
            "chembl": pd.DataFrame({"gene": ["ALB", "GPC3"]}),
        }
        overlap = overlap_analysis(targets, {"ALB", "MMP9", "GPC3"})
        self.assertEqual(len(overlap), 3)
        alb = overlap[overlap["gene"] == "ALB"].iloc[0]
        self.assertEqual(alb["n_sources"], 2)

    def test_ppi_hub_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "edges.tsv"
            pd.DataFrame(
                {
                    "protein1": ["A", "A", "B", "C"],
                    "protein2": ["B", "C", "D", "D"],
                }
            ).to_csv(path, sep="\t", index=False)
            frame = ppi_hub_scores(path, genes=["A", "B", "C", "D", "E"])
            self.assertEqual(frame.loc[frame["gene"] == "A", "ppi_degree"].iloc[0], 2)
            self.assertGreaterEqual(
                frame.loc[frame["gene"] == "A", "ppi_hub_score"].iloc[0],
                frame.loc[frame["gene"] == "E", "ppi_hub_score"].iloc[0],
            )

    def test_knockout_includes_ppi_hub(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            data_dir = workdir / "data" / "knockout"
            data_dir.mkdir(parents=True)
            genes = ["ALB", "GPC3", "MMP9", "AKT1", "TP53", "EGFR", "SRC", "CDK1"]
            pd.DataFrame(
                {
                    "gene": genes,
                    **{
                        sample: [1.0 + (i % 3) * 0.3 for i in range(len(genes))]
                        for sample in ["T1", "T2", "N1", "N2"]
                    },
                }
            ).to_csv(data_dir / "expression.csv", index=False)
            pd.DataFrame(
                {
                    "sample": ["T1", "T2", "N1", "N2"],
                    "condition": ["Tumor", "Tumor", "Normal", "Normal"],
                }
            ).to_csv(data_dir / "metadata.csv", index=False)
            pd.DataFrame(
                {
                    "protein1": ["ALB", "GPC3", "MMP9", "AKT1"],
                    "protein2": ["GPC3", "MMP9", "AKT1", "TP53"],
                }
            ).to_csv(data_dir / "ppi.tsv", sep="\t", index=False)
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "expression_csv": "data/knockout/expression.csv",
                    "metadata_csv": "data/knockout/metadata.csv",
                    "ppi_network_csv": "data/knockout/ppi.tsv",
                    "case_label": "Tumor",
                    "normal_label": "Normal",
                },
            )
            summary = run_knockout(cfg, LOG)
            self.assertTrue(summary["ppi_hub_included"])
            frame = pd.read_csv(
                cfg.knockout_dir() / "data" / "fig_52_53_ranked_knockout.csv"
            )
            self.assertIn("ppi_hub_score", frame.columns)
            self.assertIn("ppi_degree", frame.columns)


if __name__ == "__main__":
    unittest.main()
