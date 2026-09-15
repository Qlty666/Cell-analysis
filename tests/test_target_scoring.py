#!/usr/bin/env python3
"""Tests for multi-dimensional target scoring and validation handoff."""

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
from docking.validation import export_validation  # noqa: E402

DEFAULT_CONFIG = APP_ROOT / "config" / "docking_config.json"
LOG = logging.getLogger("test_target_scoring")

SIGNATURE = [
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


def _make_inputs(workdir: Path):
    data_dir = workdir / "data" / "knockout"
    data_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    genes = SIGNATURE + [
        "ALB",
        "GPC3",
        "AFP",
        "CD8A",
        "CD68",
        "COL1A1",
        "VIM",
        "SNAI1",
        "ZEB1",
        "TP53",
        "MDM2",
        "BAX",
        "PIK3CA",
        "PTEN",
        "GENE01",
        "GENE02",
        "GENE03",
        "GENE04",
        "GENE05",
        "GENE06",
        "GENE07",
        "GENE08",
        "GENE09",
        "GENE10",
    ]
    rows = {"gene": genes}
    for i, sample in enumerate(["T1", "T2", "T3", "T4", "N1", "N2", "N3", "N4"]):
        base = 3.0 if sample.startswith("T") else 1.0
        values = base + rng.normal(0, 0.2, len(genes))
        values[: len(SIGNATURE)] += 2.0 if sample.startswith("T") else 0.0
        values[genes.index("ALB")] += 2.0 if i < 2 else 0.0
        values[genes.index("CD8A")] += 2.0 if i in (2, 3) else 0.0
        values[genes.index("COL1A1")] += 2.0 if i in (4, 5) else 0.0
        rows[sample] = values
    expr = pd.DataFrame(rows)
    expr.to_csv(data_dir / "expression.csv", index=False)

    pd.DataFrame(
        {
            "sample": ["T1", "T2", "T3", "T4", "N1", "N2", "N3", "N4"],
            "condition": ["Tumor"] * 4 + ["Normal"] * 4,
            "cell_type": [
                "Hepatocyte",
                "Hepatocyte",
                "T_Cell",
                "T_Cell",
                "Fibroblast",
                "Fibroblast",
                "Macrophage",
                "Macrophage",
            ],
        }
    ).to_csv(data_dir / "metadata.csv", index=False)

    pd.DataFrame(
        {
            "gene": genes,
            "effect": [
                -1.1 if gene in SIGNATURE else 0.2 for gene in genes
            ],
        }
    ).to_csv(data_dir / "depmap.csv", index=False)

    pd.DataFrame(
        {
            "gene": genes,
            "hr": [1.8 if gene in SIGNATURE else 0.8 for gene in genes],
            "p": [0.001 if gene in SIGNATURE else 0.2 for gene in genes],
        }
    ).to_csv(data_dir / "prognosis.csv", index=False)

    pd.DataFrame(
        {
            "gene": genes,
            "known_ligands": [3 if gene in SIGNATURE else 0 for gene in genes],
            "pdb_structures": [2 if gene in SIGNATURE else 0 for gene in genes],
            "chembl_bioactivities": [
                10 if gene in SIGNATURE else 0 for gene in genes
            ],
        }
    ).to_csv(data_dir / "druggability.csv", index=False)

    pd.DataFrame(
        {
            "gene": genes,
            "off_target_paralogs": [2 if gene == "GENE01" else 0 for gene in genes],
            "safety_concern": [1 if gene == "GENE01" else 0 for gene in genes],
        }
    ).to_csv(data_dir / "off_target.csv", index=False)


class TestTargetScoring(unittest.TestCase):
    def _run(self, workdir: Path, extra_overrides: dict | None = None):
        overrides = {
            "workdir": str(workdir),
            "expression_csv": "data/knockout/expression.csv",
            "metadata_csv": "data/knockout/metadata.csv",
            "depmap_csv": "data/knockout/depmap.csv",
            "prognosis_csv": "data/knockout/prognosis.csv",
            "druggability_csv": "data/knockout/druggability.csv",
            "off_target_csv": "data/knockout/off_target.csv",
            "case_label": "Tumor",
            "normal_label": "Normal",
            "ko_top_n": 10,
        }
        overrides.update(extra_overrides or {})
        cfg = load_config(DEFAULT_CONFIG, overrides)
        summary = run_knockout(cfg, LOG)
        return cfg, summary

    def test_multidimensional_scores_and_classification(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            _make_inputs(workdir)
            cfg, summary = self._run(workdir)
            ko_dir = cfg.knockout_dir()
            frame = pd.read_csv(ko_dir / "data" / "fig_52_53_ranked_knockout.csv")

            for col in [
                "reversal_score",
                "pathway_score",
                "specificity_score",
                "prognosis_score",
                "druggability_score",
                "target_score",
                "target_class",
                "off_target_paralogs",
                "safety_concern",
            ]:
                self.assertIn(col, frame.columns)
            self.assertGreaterEqual(frame["target_score"].min(), 0.0)
            self.assertLessEqual(frame["target_score"].max(), 1.0)
            self.assertTrue(summary["multidimensional_scoring"])
            self.assertGreaterEqual(
                summary["target_class_counts"].get("core_driver", 0), 1
            )
            self.assertTrue(
                (ko_dir / "data" / "fig_52_target_candidates.csv").exists()
            )
            self.assertTrue((ko_dir / "data" / "target_report.md").exists())
            self.assertTrue((ko_dir / "run_manifest.json").exists())

            manifest = __import__("json").loads(
                (ko_dir / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["stage"], "virtual-knockout")
            self.assertIn("config_sha256", manifest)
            self.assertIn("inputs", manifest)

    def test_off_target_penalty_lowers_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            _make_inputs(workdir)
            cfg, _ = self._run(workdir, {"ko_top_n": 50})
            frame = pd.read_csv(
                cfg.knockout_dir() / "data" / "fig_52_53_ranked_knockout.csv"
            )
            flagged = frame.loc[frame["gene"] == "GENE01"].iloc[0]
            self.assertEqual(flagged["safety_concern"], 1)
            self.assertEqual(flagged["off_target_paralogs"], 2)

    def test_export_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            _make_inputs(workdir)
            cfg, _ = self._run(workdir)
            summary = export_validation(cfg, LOG)
            val_dir = cfg.validation_dir()
            self.assertTrue((val_dir / "data" / "validation_candidates.csv").exists())
            self.assertTrue((val_dir / "data" / "validation_plan.md").exists())
            self.assertTrue((val_dir / "data" / "summary.json").exists())
            self.assertTrue((val_dir / "data" / "run_manifest.json").exists())
            self.assertEqual(summary["candidates"], 10)
            candidates = pd.read_csv(val_dir / "data" / "validation_candidates.csv")
            self.assertIn("phase_1_cell_line_assay", candidates.columns)
            self.assertIn("target_class", candidates.columns)


if __name__ == "__main__":
    unittest.main()
