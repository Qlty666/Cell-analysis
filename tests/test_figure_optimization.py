#!/usr/bin/env python3
"""Tests for non-destructive experiment-plan figure optimization."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from PIL import Image

from optimize_experiment_plan_figures import optimize


class TestFigureOptimization(unittest.TestCase):
    def test_optimizer_writes_rgb_manifest_and_not_run_panels(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp) / "results"
            output = Path(tmp) / "optimized"
            root.mkdir()
            manifest = optimize(root, output)
            self.assertTrue(
                (output / "figure_optimization_manifest.json").exists()
            )
            self.assertTrue(manifest["files"])
            panel = output / "09_md_mmpbsa" / "fig5d_rmsd.png"
            self.assertTrue(panel.exists())
            with Image.open(panel) as image:
                self.assertEqual(image.mode, "RGB")
                self.assertFalse(image.info.get("transparency"))

    def test_optimizer_does_not_modify_source_tables(self):
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmp:
            root = Path(tmp) / "results"
            output = Path(tmp) / "optimized"
            source = root / "01_compound_characterization" / "source.csv"
            source.parent.mkdir(parents=True)
            source.write_text("gene,score\nEGFR,0.9\n", encoding="utf-8")
            before = source.read_bytes()
            optimize(root, output)
            self.assertEqual(source.read_bytes(), before)
            self.assertIn(
                "note",
                json.loads(
                    (
                        output / "figure_optimization_manifest.json"
                    ).read_text(encoding="utf-8")
                ),
            )


if __name__ == "__main__":
    unittest.main()
