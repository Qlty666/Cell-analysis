#!/usr/bin/env python3
"""Tests for the generic GEO supplementary-file downloader."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from data.geo_downloader import _select_files  # noqa: E402


class TestSelectFiles(unittest.TestCase):
    def test_recognizes_rna_seq_raw_count_table(self):
        names = [
            "GSE342934_RNA-Seq.txt.gz",
            "GSE342934_all.fpkm_anno.txt.gz",
            "GSE342934_series_matrix.txt.gz",
        ]
        selected = _select_files(names)
        self.assertEqual(selected["matrix"], ["GSE342934_RNA-Seq.txt.gz"])

    def test_skips_normalized_expression_tables(self):
        names = [
            "GSE12345_raw_counts.tsv.gz",
            "GSE12345_normalized_tpm.tsv.gz",
            "GSE12345_fpkm_anno.txt.gz",
        ]
        selected = _select_files(names)
        self.assertEqual(selected["matrix"], ["GSE12345_raw_counts.tsv.gz"])


if __name__ == "__main__":
    unittest.main()
