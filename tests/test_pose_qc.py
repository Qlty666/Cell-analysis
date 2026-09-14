"""Tests for optional PoseBusters pose validation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pipeline.pose_qc import run_pose_qc


class TestPoseQc(unittest.TestCase):
    def test_missing_sdf_is_reported_as_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            summary = run_pose_qc(root, root / "out")
            self.assertEqual(summary["status"], "skipped")
            self.assertFalse(summary["gate_passed"])
            self.assertTrue((root / "out" / "pose_qc_results.csv").exists())


if __name__ == "__main__":
    unittest.main()
