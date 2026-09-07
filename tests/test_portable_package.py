#!/usr/bin/env python3
"""Tests for the clean source package used on new computers."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "launchers") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "launchers"))

from package_portable import _excluded, build_package, source_files  # noqa: E402


class TestPortablePackage(unittest.TestCase):
    def test_source_files_are_clean(self):
        files = source_files()
        self.assertTrue(files)
        relative = [path.relative_to(APP_ROOT).as_posix() for path in files]
        for prefix in (
            "data_cache/",
            "logs/",
            "results/",
            "dock/tools/",
            "dock/outputs/",
            "portable/",
        ):
            self.assertFalse(
                any(
                    name.startswith(prefix) and not name.endswith(".gitkeep")
                    for name in relative
                ),
                prefix,
            )
        self.assertIn("src/docking/network_toxicology.py", relative)
        self.assertFalse(
            any(name.startswith("scripts/_test_h5ad_reader.R") for name in relative)
        )
        for path in files:
            if path.suffix.lower() not in {
                ".py",
                ".r",
                ".md",
                ".json",
                ".yml",
                ".yaml",
                ".bat",
                ".sh",
            }:
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            self.assertNotIn("D:\\AAA Liver cancer", text, str(path))
            self.assertNotIn("C:\\Users\\20338", text, str(path))

    def test_build_package_dry_run_does_not_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "portable" / "package.zip"
            planned = build_package(output, dry_run=True)
            self.assertEqual(planned, output.resolve())
            self.assertFalse(output.exists())

    def test_excluded_ignores_generated_directories(self):
        self.assertTrue(_excluded(APP_ROOT / "data_cache" / "index.json"))
        self.assertTrue(_excluded(APP_ROOT / "logs" / "run.log"))
        self.assertTrue(_excluded(APP_ROOT / "src" / "__pycache__" / "x.pyc"))


if __name__ == "__main__":
    unittest.main()
