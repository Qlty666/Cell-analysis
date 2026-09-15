"""Tests for docking and MD structural quality gates."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from pipeline.structural_quality import build_structural_quality


class TestStructuralQuality(unittest.TestCase):
    def _build_case(
        self,
        tmp: str,
        *,
        protein_std: float,
        ligand_std: float,
        mmpbsa_values: list[float],
    ):
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
                "rmsd_protein_tail_std_nm": [protein_std],
                "rmsd_ligand_tail_std_nm": [ligand_std],
                "mmpbsa_delta_g": mmpbsa_values,
            }
        ).to_csv(
            target / "06_md" / "md_simulation_results.csv",
            index=False,
        )
        return workdir, integration, target

    def test_collects_docking_and_md_gates(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir, integration, _ = self._build_case(
                tmp,
                protein_std=0.08,
                ligand_std=0.10,
                mmpbsa_values=[-30.0],
            )
            summary = build_structural_quality(
                workdir,
                integration,
            )
            self.assertTrue(summary["all_structural_gates"])
            self.assertTrue(summary["gates"]["md_rmsd_stability"])
            self.assertTrue(summary["gates"]["mmpbsa_available"])
            target_frame = pd.read_csv(
                integration / "structural_quality_targets.csv"
            )
            self.assertTrue(bool(target_frame.loc[0, "md_rmsd_stable"]))

    def test_protein_stability_does_not_mask_unstable_ligand(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir, integration, _ = self._build_case(
                tmp,
                protein_std=0.05,
                ligand_std=0.50,
                mmpbsa_values=[-30.0],
            )
            summary = build_structural_quality(workdir, integration)
            self.assertFalse(summary["gates"]["md_rmsd_stability"])

    def test_ligand_stability_does_not_mask_unstable_protein(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir, integration, _ = self._build_case(
                tmp,
                protein_std=0.50,
                ligand_std=0.05,
                mmpbsa_values=[-30.0],
            )
            summary = build_structural_quality(workdir, integration)
            self.assertFalse(summary["gates"]["md_rmsd_stability"])

    def test_all_zero_or_missing_mmpbsa_is_not_available(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir, integration, _ = self._build_case(
                tmp,
                protein_std=0.08,
                ligand_std=0.10,
                mmpbsa_values=[0.0],
            )
            summary = build_structural_quality(workdir, integration)
            self.assertFalse(summary["gates"]["mmpbsa_available"])


if __name__ == "__main__":
    unittest.main()
