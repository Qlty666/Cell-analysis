"""Tests for optional PoseBusters pose validation."""

from __future__ import annotations

import tempfile
import sys
import types
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from pipeline.pose_qc import run_pose_qc


class TestPoseQc(unittest.TestCase):
    def test_missing_sdf_is_reported_as_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            summary = run_pose_qc(root, root / "out")
            self.assertEqual(summary["status"], "skipped")
            self.assertFalse(summary["gate_passed"])
            self.assertTrue((root / "out" / "pose_qc_results.csv").exists())

    def test_pdbqt_pose_is_converted_before_posebusters_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "work" / "GENE1" / "outputs" / "docked"
            target.mkdir(parents=True)
            pdbqt = target / "pose1.pdbqt"
            pdbqt.write_text("REMARK test\n", encoding="utf-8")
            converted = root / "out" / "poses_sdf" / "GENE1" / "pose1.sdf"
            receptor_dir = root / "work" / "GENE1" / "data" / "receptors"
            receptor_dir.mkdir(parents=True)
            receptor = receptor_dir / "receptor.pdb"
            receptor.write_text("ATOM\n", encoding="utf-8")

            def fake_convert(path, sdf_dir, gene):
                converted.parent.mkdir(parents=True, exist_ok=True)
                converted.write_text("$$$$\n", encoding="utf-8")
                return converted, ""

            with (
                mock.patch(
                    "pipeline.pose_qc._convert_pdbqt_to_sdf",
                    side_effect=fake_convert,
                ),
                mock.patch(
                    "pipeline.pose_qc.importlib.util.find_spec",
                    return_value=None,
                ),
            ):
                summary = run_pose_qc(root, root / "out")
            self.assertEqual(summary["status"], "unavailable")
            self.assertEqual(summary["sdf_files"], 1)
            self.assertEqual(summary["pose_format"], "pdbqt_converted")

    def test_converted_pdbqt_keeps_original_receptor_for_posebusters(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / "work" / "GENE1" / "outputs" / "docked"
            target.mkdir(parents=True)
            pdbqt = target / "pose1.pdbqt"
            pdbqt.write_text("REMARK test\n", encoding="utf-8")
            receptor_dir = root / "work" / "GENE1" / "data" / "receptors"
            receptor_dir.mkdir(parents=True)
            receptor = receptor_dir / "receptor.pdb"
            receptor.write_text("ATOM\n", encoding="utf-8")
            converted = root / "out" / "poses_sdf" / "GENE1" / "pose1.sdf"

            def fake_convert(path, sdf_dir, gene):
                converted.parent.mkdir(parents=True, exist_ok=True)
                converted.write_text("$$$$\n", encoding="utf-8")
                return converted, ""

            calls = []

            class FakePoseBusters:
                def __init__(self, config):
                    self.config = config

                def bust(self, *, mol_pred, mol_cond):
                    calls.append((mol_pred, mol_cond))
                    return pd.DataFrame(
                        {
                            "mol_pred_loaded": [True],
                            "sanity": [True],
                            "rmsd": [0.1],
                        }
                    )

            fake_module = types.ModuleType("posebusters")
            fake_module.PoseBusters = FakePoseBusters
            with (
                mock.patch(
                    "pipeline.pose_qc._convert_pdbqt_to_sdf",
                    side_effect=fake_convert,
                ),
                mock.patch(
                    "pipeline.pose_qc.importlib.util.find_spec",
                    return_value=object(),
                ),
                mock.patch.dict(sys.modules, {"posebusters": fake_module}),
            ):
                summary = run_pose_qc(root, root / "out")

            self.assertEqual(summary["status"], "completed")
            self.assertTrue(summary["gate_passed"])
            self.assertEqual(summary["pose_format"], "pdbqt_converted")
            self.assertEqual(summary["no_receptor"], 0)
            self.assertEqual(summary["receptors_resolved"], 1)
            self.assertEqual(len(calls), 1)
            self.assertEqual(Path(calls[0][1]), receptor)


if __name__ == "__main__":
    unittest.main()
