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
    def _ensure(
        self,
        tmp: str,
        organism: str = "Homo sapiens",
        description: str = "bulk RNA-seq",
        accession: str = "E-MTAB-1",
    ) -> dict:
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
                "attributes": [{"name": "RootPath", "value": accession}],
                "section": {
                    "attributes": [
                        {"name": "Title", "value": "HCC counts"},
                        {"name": "Organism", "value": organism},
                        {"name": "Description", "value": description},
                    ]
                },
            }

        def fake_download(url, out, log):
            out.parent.mkdir(parents=True, exist_ok=True)
            if out.name.lower().endswith(".txt"):
                out.write_bytes(b"gene\tcell1\nA\t1\nB\t2\n")
            else:
                out.write_bytes(b"A\nB\n")

        with (
            mock.patch.object(bsd, "_http_get_json", side_effect=fake_json),
            mock.patch.object(bsd.gd, "_download", side_effect=fake_download),
            mock.patch.object(bsd, "CACHE_ROOT", Path(tmp) / "cache"),
        ):
            return bsd.ensure_biostudies_dataset(
                accession,
                Path(tmp) / "root",
                lambda _msg: None,
            )

    def test_writes_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = self._ensure(tmp)
            self.assertEqual(manifest["accession"], "E-MTAB-1")
            self.assertEqual(manifest["organism"], "hs")
            self.assertIn("GSM1_counts.txt", manifest["files"]["matrix"])
            self.assertEqual(manifest["files"]["genes"], ["genes.tsv"])
            manifest_path = Path(tmp) / "root" / "data" / "E-MTAB-1_manifest.json"
            self.assertTrue(manifest_path.exists())
            self.assertTrue(
                (
                    Path(tmp)
                    / "root"
                    / "data"
                    / "raw"
                    / "E-MTAB-1"
                    / "GSM1_counts.txt"
                ).exists()
            )

    def test_maps_mouse_and_zebrafish_organisms(self):
        with tempfile.TemporaryDirectory() as tmp:
            mouse = self._ensure(tmp, organism="Mus musculus")
            self.assertEqual(mouse["organism"], "mm")
        with tempfile.TemporaryDirectory() as tmp:
            zebrafish = self._ensure(tmp, organism="Danio rerio")
            self.assertEqual(zebrafish["organism"], "dr")

    def test_unknown_organism_is_recorded_not_defaulted_to_human(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = self._ensure(
                tmp,
                organism="Glycine max",
                description="soybean seed development",
            )
            self.assertEqual(manifest["organism"], "unknown")

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
