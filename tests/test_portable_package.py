#!/usr/bin/env python3
"""Tests for the clean source package used on new computers."""

from __future__ import annotations

import sys
import tempfile
import unittest
import zipfile
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
        home = Path.home()
        # Personal absolute paths must never be committed into the portable
        # source. ``Path.home()`` is resolved at runtime, so no username or
        # machine-specific directory is hardcoded here. When the checkout
        # itself lives inside the home directory, the absolute home prefix is
        # not a reliable leak marker, so it is skipped.
        home_inside_repo = home == APP_ROOT or APP_ROOT in home.parents
        home_tokens = (
            ()
            if home_inside_repo
            else tuple(
                dict.fromkeys(
                    (
                        str(home),
                        str(home).replace("\\", "/"),
                        f"Users\\{home.name}",
                        f"Users/{home.name}",
                        f"/home/{home.name}",
                    )
                )
            )
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
            repo_tokens = (
                str(APP_ROOT),
                str(APP_ROOT).replace("\\", "/"),
            )
            for token in repo_tokens:
                self.assertNotIn(token, text, str(path))
            for token in home_tokens:
                self.assertNotIn(token, text, str(path))

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

    def test_excluded_rejects_secret_files(self):
        for name in (
            ".env",
            ".env.local",
            "server.pem",
            "client.p12",
            "private.key",
            "bundle.pfx",
        ):
            with self.subTest(name=name):
                self.assertTrue(_excluded(APP_ROOT / name))
        self.assertTrue(
            _excluded(APP_ROOT / "config" / ".env.production")
        )

    def test_excluded_keeps_normal_sources(self):
        self.assertFalse(
            _excluded(APP_ROOT / "src" / "docking" / "network_toxicology.py")
        )
        self.assertFalse(_excluded(APP_ROOT / "requirements.txt"))

    def test_build_package_refuses_existing_output_without_force(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "package.zip"
            output.write_bytes(b"stale archive")
            with self.assertRaises(SystemExit):
                build_package(output)
            self.assertEqual(output.read_bytes(), b"stale archive")
            build_package(output, force=True)
            self.assertTrue(zipfile.is_zipfile(output))
            self.assertGreater(output.stat().st_size, len(b"stale archive"))


if __name__ == "__main__":
    unittest.main()
