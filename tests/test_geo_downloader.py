#!/usr/bin/env python3
"""Tests for the generic GEO supplementary-file downloader."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from data.geo_downloader import (  # noqa: E402
    BULK_COUNT_TABLE_RE,
    _select_files,
)


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

    def test_detects_per_sample_bulk_count_tables(self):
        names = [
            "GSM9037276_GC119559.txt.gz",
            "GSM9037277_GC119560.txt.gz",
            "GSE299321_series_matrix.txt.gz",
        ]
        selected = _select_files(names)
        self.assertEqual(
            selected["matrix"],
            ["GSM9037276_GC119559.txt.gz", "GSM9037277_GC119560.txt.gz"],
        )
        self.assertTrue(selected["bulk"])

    def test_10x_files_are_not_bulk(self):
        names = [
            "GSE125449_Set1_barcodes.tsv.gz",
            "GSE125449_Set1_genes.tsv.gz",
            "GSE125449_Set1_matrix.mtx.gz",
        ]
        selected = _select_files(names)
        self.assertEqual(selected["matrix"], ["GSE125449_Set1_matrix.mtx.gz"])
        self.assertFalse(selected["bulk"])

    def test_bulk_count_table_regex(self):
        self.assertTrue(BULK_COUNT_TABLE_RE.search("GSM9037276_GC119559.txt.gz"))
        self.assertTrue(
            BULK_COUNT_TABLE_RE.search("_extracted/GSM9037276_counts.tsv.gz")
        )
        self.assertFalse(BULK_COUNT_TABLE_RE.search("GSE125449_Set1_matrix.mtx.gz"))
        self.assertFalse(BULK_COUNT_TABLE_RE.search("GSE299321_RAW.tar"))

if __name__ == "__main__":
    unittest.main()
