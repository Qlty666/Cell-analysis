#!/usr/bin/env python3
"""Security hardening tests."""

from __future__ import annotations

import io
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))
if str(APP_ROOT / "web") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "web"))

from data.geo_downloader import _extract_archive, normalize_accession  # noqa: E402
from web_ui import start_full_job  # noqa: E402


class TestAccessionValidation(unittest.TestCase):
    def test_normalize_accession_accepts_valid_gse(self):
        self.assertEqual(normalize_accession("gse125449"), "GSE125449")

    def test_normalize_accession_rejects_path_traversal(self):
        for value in ["../../evil", "GSE1/../x", "GSE", "GSE12x"]:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    normalize_accession(value)

    def test_start_full_job_rejects_invalid_accession(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            with self.assertRaises(ValueError):
                start_full_job(
                    {
                        "output": str(base / "out"),
                        "workdir": str(base / "work"),
                        "accession": "../../evil",
                    }
                )


class TestSafeArchiveExtraction(unittest.TestCase):
    def test_tar_member_escaping_destination_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            archive = base / "bad.tar"
            dest = base / "dest"
            dest.mkdir()
            with tarfile.open(archive, "w") as tf:
                data = b"evil"
                info = tarfile.TarInfo("../escape.txt")
                info.size = len(data)
                tf.addfile(info, io.BytesIO(data))
            with self.assertRaises(RuntimeError):
                _extract_archive(archive, dest)
            self.assertFalse((base / "escape.txt").exists())


if __name__ == "__main__":
    unittest.main()
