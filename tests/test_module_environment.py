#!/usr/bin/env python3
"""Tests for the per-module environment installer."""

from __future__ import annotations

import io
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "launchers") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "launchers"))

from install_environment import (  # noqa: E402
    MODULES,
    _download,
    check_module,
    install_module,
    main,
    md5_file,
    resolve_module,
    sha256_file,
    verify_checksum,
)


class TestChecksums(unittest.TestCase):
    PAYLOAD = b"liver-cancer"
    SHA256 = "e46aa0f5f3e14b088ae65ff4dc8158ca1ff4dffbf3a73ea6b9778c2b1cc276cc"
    MD5 = "0c16659c807c8c6df298926073219cf3"

    def test_digest_helpers(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "payload.bin"
            path.write_bytes(self.PAYLOAD)
            self.assertEqual(sha256_file(path), self.SHA256)
            self.assertEqual(md5_file(path), self.MD5)

    def test_verify_checksum_accepts_matching_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "payload.bin"
            path.write_bytes(self.PAYLOAD)
            self.assertTrue(verify_checksum(path, self.SHA256))
            self.assertTrue(verify_checksum(path, self.MD5, "md5"))

    def test_verify_checksum_rejects_missing_file_and_bad_digest(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "payload.bin"
            self.assertFalse(verify_checksum(path, self.SHA256))
            path.write_bytes(self.PAYLOAD)
            self.assertFalse(verify_checksum(path, "0" * 64))
            self.assertFalse(verify_checksum(path, ""))
            self.assertFalse(verify_checksum(path, "not-a-digest", "md5"))

    def test_verify_checksum_rejects_unknown_algorithm(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "payload.bin"
            path.write_bytes(self.PAYLOAD)
            with self.assertRaises(ValueError):
                verify_checksum(path, self.SHA256, "crc32")

    def test_download_failure_removes_partial_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "partial.bin"
            dest.write_bytes(b"stale")
            missing = (Path(tmp) / "missing.bin").as_uri()
            self.assertFalse(_download(missing, dest))
            self.assertFalse(dest.exists())


class TestModuleEnvironment(unittest.TestCase):
    def test_all_module_entrypoints_and_requirements_exist(self):
        for name, spec in MODULES.items():
            with self.subTest(module=name):
                for relative in spec["entrypoints"]:
                    if " " in relative:
                        continue
                    self.assertTrue(
                        (APP_ROOT / relative).is_file(),
                        f"{relative} should exist",
                    )
                for relative in spec["pip_requirements"]:
                    self.assertTrue(
                        (APP_ROOT / relative).is_file(),
                        f"{relative} should exist",
                    )

    def test_module_aliases_resolve(self):
        self.assertEqual(resolve_module("single"), "expression")
        self.assertEqual(resolve_module("pipeline"), "expression")
        self.assertEqual(resolve_module("dock"), "docking")
        self.assertEqual(resolve_module("all"), "full")

    def test_unknown_module_raises(self):
        with self.assertRaises(ValueError):
            resolve_module("not-a-module")

    def test_list_prints_all_modules(self):
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["list"])
        self.assertEqual(code, 0)
        for name in MODULES:
            self.assertIn(name, out.getvalue())

    def test_web_check_dry_run(self):
        self.assertEqual(check_module("web", dry_run=True), 0)

    def test_expression_install_dry_run(self):
        self.assertEqual(
            install_module(
                "expression",
                auto_install_r=False,
                dry_run=True,
            ),
            0,
        )

    def test_full_install_dry_run(self):
        self.assertEqual(
            install_module(
                "full",
                auto_install_r=False,
                dry_run=True,
            ),
            0,
        )


if __name__ == "__main__":
    unittest.main()
