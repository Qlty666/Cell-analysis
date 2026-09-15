#!/usr/bin/env python3
"""Tests for the DOCX/PDF report export status."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from report import export_report  # noqa: E402


def _fpdf_available() -> bool:
    try:
        import fpdf  # noqa: F401
    except ImportError:
        return False
    return True


class TestExportReport(unittest.TestCase):
    @staticmethod
    def _make_root(tmp: str) -> Path:
        root = Path(tmp)
        results = root / "results"
        results.mkdir(parents=True)
        (results / "summary.json").write_text(
            json.dumps(
                {
                    "dataset": "GSE_TEST",
                    "title": "test",
                    "dataset_mode": "single_cell",
                    "n_cells_raw": 10,
                    "n_cells_after_qc": 8,
                    "n_cells_after_doublet_removal": 8,
                    "n_genes": 100,
                    "n_clusters": 3,
                    "deg_up": 2,
                    "deg_down": 1,
                }
            ),
            encoding="utf-8",
        )
        return root

    def test_missing_summary_returns_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.object(
                sys, "argv", ["export_report.py", str(Path(tmp) / "missing")]
            ):
                self.assertEqual(export_report.main(), 1)

    def test_pdf_is_not_claimed_when_fpdf_is_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self._make_root(tmp)
            with mock.patch.dict(sys.modules, {"fpdf": None}):
                with mock.patch.object(
                    sys, "argv", ["export_report.py", str(root)]
                ):
                    self.assertEqual(export_report.main(), 0)
            status = (root / "results" / "export_status.txt").read_text(
                encoding="utf-8"
            )
            self.assertIn("DOCX export completed", status)
            self.assertIn("PDF export skipped", status)
            self.assertNotIn("PDF export completed", status)
            self.assertFalse((root / "results" / "result_report.pdf").exists())

    @unittest.skipUnless(_fpdf_available(), "fpdf2 not installed")
    def test_pdf_is_claimed_only_after_it_is_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self._make_root(tmp)
            with mock.patch.object(sys, "argv", ["export_report.py", str(root)]):
                self.assertEqual(export_report.main(), 0)
            pdf = root / "results" / "result_report.pdf"
            self.assertTrue(pdf.is_file())
            self.assertGreater(pdf.stat().st_size, 0)
            status = (root / "results" / "export_status.txt").read_text(
                encoding="utf-8"
            )
            self.assertIn("PDF export completed", status)


if __name__ == "__main__":
    unittest.main()
