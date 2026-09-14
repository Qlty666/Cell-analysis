"""Tests for the reproducibility manifest."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from pipeline.reproducibility import write_reproducibility_manifest


class TestReproducibilityManifest(unittest.TestCase):
    def test_manifest_records_git_packages_and_hashes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            workdir = root / "work"
            single_cell = root / "single_cell"
            (single_cell / "results").mkdir(parents=True)
            (single_cell / "results" / "summary.json").write_text(
                "{}",
                encoding="utf-8",
            )
            config = root / "config.json"
            config.write_text("{}", encoding="utf-8")
            docking = root / "docking.json"
            docking.write_text("{}", encoding="utf-8")
            manifest_path = write_reproducibility_manifest(
                workdir,
                single_cell,
                config,
                docking,
                {"seed": 42},
            )
            manifest = json.loads(
                manifest_path.read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["schema_version"], "1.0")
            self.assertIn("commit", manifest["git"])
            self.assertEqual(
                manifest["inputs"]["full_pipeline_config"]["sha256"],
                manifest["inputs"]["docking_config"]["sha256"],
            )
            self.assertIn("pandas", manifest["packages"])


if __name__ == "__main__":
    unittest.main()
