"""Tests for docking and MD structural quality gates."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.structural_quality import build_structural_quality


class TestStructuralQuality(unittest.TestCase):
    def test_collects_docking_and_md_gates(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            integration = workdir / "outputs" / "integration"
            integration.mkdir(parents=True)
            pd.DataFrame(
                {
                    "gene": ["GENE1"],
                    "status": ["ok"],
                    "best_affinity": [-8.0],
                    "hits": [2],
                }
            ).to_csv(integration / "docking_targets.csv", index=False)
            (integration / "docking_summary.json").write_text(
                json.dumps({"positive_control": {"passed": True}}),
                encoding="utf-8",
            )
            target = (
                workdir
                / "work"
                / "GENE1"
                / "outputs"
                / "run_001"
                / "results"
            )
            (target / "01_analysis").mkdir(parents=True)
            (target / "06_md").mkdir(parents=True)
            (target / "01_analysis" / "summary.json").write_text(
                json.dumps({"replicate_consensus": True}),
                encoding="utf-8",
            )
            (target / "06_md" / "md_simulation_summary.json").write_text(
                json.dumps(
                    {
                        "mode": "auto",
                        "completed": 1,
                        "prepared": 0,
                        "total_ns": 100.0,
                    }
                ),
                encoding="utf-8",
            )
            pd.DataFrame(
                {
                    "id": ["lig1"],
                    "rmsd_protein_tail_std_nm": [0.08],
                    "rmsd_ligand_tail_std_nm": [0.10],
                    "mmpbsa_delta_g": [-30.0],
                }
            ).to_csv(
                target / "06_md" / "md_simulation_results.csv",
                index=False,
            )
            summary = build_structural_quality(
                workdir,
                integration,
            )
            self.assertTrue(summary["all_structural_gates"])
            self.assertTrue(summary["gates"]["md_rmsd_stability"])
            self.assertTrue(summary["gates"]["mmpbsa_available"])


if __name__ == "__main__":
    unittest.main()
