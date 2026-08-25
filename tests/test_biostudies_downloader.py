#!/usr/bin/env python3
"""Tests for the EBI BioStudies/ArrayExpress downloader."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from data import biostudies_downloader as bsd  # noqa: E402


class TestBioStudiesAccession(unittest.TestCase):
    def test_normalize_accession(self):
        self.assertEqual(bsd.normalize_accession("e-mtab-1234"), "E-MTAB-1234")
        self.assertEqual(bsd.normalize_accession("s-bsst1"), "S-BSST1")

    def test_normalize_accession_rejects_geo(self):
        with self.assertRaises(ValueError):
            bsd.normalize_accession("GSE1")


class TestSelectProcessedFiles(unittest.TestCase):
    def test_selects_count_matrix_and_skips_raw(self):
        items = [
            {"path": "GSM1_counts.txt", "Section": "processed-data"},
            {"path": "GSM2_sample_table.txt", "Section": "processed-data"},
            {"path": "genes.tsv", "Section": "processed-data"},
            {"path": "sample1.fastq.gz", "Section": "raw-data"},
        ]
        selected = bsd._select_processed_files(items)
        self.assertEqual(
            selected["matrix"],
            ["GSM1_counts.txt", "GSM2_sample_table.txt"],
        )
        self.assertEqual(selected["genes"], ["genes.tsv"])

    def test_rejects_traversal_path(self):
        with self.assertRaises(RuntimeError):
            bsd._select_processed_files(
                [{"path": "../escape.txt", "Section": "processed-data"}]
            )


class TestEnsureBioStudiesDataset(unittest.TestCase):
    def test_writes_manifest(self):
        def fake_json(url):
            if url.endswith("/files"):
                return {
                    "items": [
                        {
                            "path": "GSM1_counts.txt",
                            "Section": "processed-data",
                            "Size": 100,
                        },
                        {
                            "path": "genes.tsv",
                            "Section": "processed-data",
                            "Size": 100,
                        },
                        {
                            "path": "sample1.fastq.gz",
                            "Section": "raw-data",
                            "Size": 100,
                        },
                    ]
                }
            return {
                "attributes": [{"name": "RootPath", "value": "E-MTAB-1"}],
                "section": {
                    "attributes": [
                        {"name": "Title", "value": "HCC counts"},
                        {"name": "Organism", "value": "Homo sapiens"},
                        {"name": "Description", "value": "bulk RNA-seq"},
                    ]
                },
            }

        def fake_download(url, out, log):
            out.parent.mkdir(parents=True, exist_ok=True)
            if out.name.lower().endswith(".txt"):
                out.write_bytes(b"gene\tcell1\nA\t1\nB\t2\n")
            else:
                out.write_bytes(b"A\nB\n")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "root"
            logs: list[str] = []
            with (
                mock.patch.object(bsd, "_http_get_json", side_effect=fake_json),
                mock.patch.object(bsd.gd, "_download", side_effect=fake_download),
            ):
                manifest = bsd.ensure_biostudies_dataset(
                    "E-MTAB-1",
                    root,
                    logs.append,
                )
            self.assertEqual(manifest["accession"], "E-MTAB-1")
            self.assertEqual(manifest["organism"], "hs")
            self.assertIn("GSM1_counts.txt", manifest["files"]["matrix"])
            self.assertEqual(manifest["files"]["genes"], ["genes.tsv"])
            manifest_path = root / "data" / "E-MTAB-1_manifest.json"
            self.assertTrue(manifest_path.exists())
            self.assertTrue((root / "data" / "raw" / "E-MTAB-1" / "GSM1_counts.txt").exists())

    def test_no_matrix_raises(self):
        def fake_json(url):
            if url.endswith("/files"):
                return {
                    "items": [
                        {
                            "path": "sample1.fastq.gz",
                            "Section": "raw-data",
                            "Size": 100,
                        }
                    ]
                }
            return {
                "attributes": [{"name": "RootPath", "value": "E-MTAB-2"}],
                "section": {
                    "attributes": [
                        {"name": "Title", "value": "raw only"},
                        {"name": "Organism", "value": "Homo sapiens"},
                        {"name": "Description", "value": "raw"},
                    ]
                },
            }

        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.object(bsd, "_http_get_json", side_effect=fake_json):
                with self.assertRaisesRegex(RuntimeError, "No processed count matrix"):
                    bsd.ensure_biostudies_dataset(
                        "E-MTAB-2",
                        Path(tmp) / "root",
                        lambda _msg: None,
                    )


if __name__ == "__main__":
    unittest.main()
