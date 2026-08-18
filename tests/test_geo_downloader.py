#!/usr/bin/env python3
"""Tests for the generic GEO supplementary-file downloader."""

from __future__ import annotations

import gzip
import io
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest import mock

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from data.geo_downloader import (  # noqa: E402
    BULK_COUNT_TABLE_RE,
    _convert_downloaded,
    _download,
    _expand_archive_files,
    _matrix_header_looks_single_cell,
    _refresh_manifest_mode,
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

    def test_archive_names_are_not_selected_as_matrix(self):
        names = [
            "GSE254513_CG_raw_feature_bc_matrix.tar.gz",
            "GSE254513_CG_Seurat_metadata.csv.gz",
            "GSE254513_series_matrix.txt.gz",
        ]
        selected = _select_files(names)
        self.assertEqual(selected["matrix"], [])
        self.assertEqual(selected["metadata"], ["GSE254513_CG_Seurat_metadata.csv.gz"])
        self.assertFalse(selected["bulk"])

    def test_filtered_feature_bc_matrix_rds_not_bulk(self):
        names = [
            "_extracted/GSM9947997_Organoid_pool_DMSO_filtered_feature_bc_matrix.rds",
            "_extracted/GSM9947998_Organoid_pool_JTE607_filtered_feature_bc_matrix.rds",
            "_extracted/GSM9947999_LACO1_FPKM_quan.tab.gz",
        ]
        selected = _select_files(names)
        self.assertEqual(selected["matrix"], names[:2])
        self.assertFalse(selected["bulk"])

    def test_rds_count_matrix_not_bulk(self):
        names = ["GSM9037276_counts.rds"]
        selected = _select_files(names)
        self.assertEqual(selected["matrix"], names)
        self.assertFalse(selected["bulk"])

    def test_recognizes_h5ad_and_loom_as_single_cell_matrices(self):
        names = [
            "_extracted/GSM9440169_E12_B.h5ad",
            "_extracted/GSM9440170_E12_F1.loom",
        ]
        selected = _select_files(names)
        self.assertEqual(selected["matrix"], names)
        self.assertFalse(selected["bulk"])

    def test_bulk_count_table_regex(self):
        self.assertTrue(BULK_COUNT_TABLE_RE.search("GSM9037276_GC119559.txt.gz"))
        self.assertTrue(
            BULK_COUNT_TABLE_RE.search("_extracted/GSM9037276_counts.tsv.gz")
        )
        self.assertFalse(BULK_COUNT_TABLE_RE.search("GSE125449_Set1_matrix.mtx.gz"))
        self.assertFalse(BULK_COUNT_TABLE_RE.search("GSE299321_RAW.tar"))


class TestRefreshManifestMode(unittest.TestCase):
    def test_reclassifies_seurat_rds_from_old_bulk_manifest(self):
        manifest = {
            "accession": "GSE343226",
            "mode": "bulk",
            "organism": "hs",
            "files": {
                "matrix": [
                    "_extracted/GSM9947997_Organoid_pool_DMSO_filtered_feature_bc_matrix.rds",
                    "_extracted/GSM9947998_Organoid_pool_JTE607_filtered_feature_bc_matrix.rds",
                ],
                "barcodes": [],
                "genes": [],
                "metadata": [],
                "series_matrices": ["GSE343226_series_matrix.txt.gz"],
            },
        }
        updated = _refresh_manifest_mode(manifest)
        self.assertEqual(updated["mode"], "single_cell")

    def test_keeps_per_sample_count_tables_bulk(self):
        manifest = {
            "accession": "GSE299321",
            "mode": "bulk",
            "organism": "hs",
            "files": {
                "matrix": ["GSM9037276_GC119559.txt.gz"],
                "barcodes": [],
                "genes": [],
                "metadata": [],
                "series_matrices": ["GSE299321_series_matrix.txt.gz"],
            },
        }
        updated = _refresh_manifest_mode(manifest)
        self.assertEqual(updated["mode"], "bulk")

    def test_single_cell_hint_overrides_bulk_file_pattern(self):
        manifest = {
            "accession": "GSE165816",
            "mode": "bulk",
            "single_cell_hint": True,
            "files": {
                "matrix": ["_extracted/GSM5050521_G1counts.csv.gz"],
                "barcodes": [],
                "genes": [],
                "metadata": [],
                "series_matrices": [],
            },
        }
        updated = _refresh_manifest_mode(manifest)
        self.assertEqual(updated["mode"], "single_cell")

    def test_barcode_header_csv_is_reclassified_as_single_cell(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            matrix = base / "GSM5050521_G1counts.csv.gz"
            with gzip.open(matrix, "wt", encoding="utf-8") as fh:
                fh.write(
                    "AAACCTGAGGCTCTTA,AAACCTGCAACTTGAC,AAACCTGCACGAAATA\n"
                    "AL627309.1,0,0,1\n"
                    "AL669831.5,1,0,0\n"
                )
            manifest = {
                "accession": "GSE165816",
                "mode": "bulk",
                "files": {
                    "matrix": [matrix.name],
                    "barcodes": [],
                    "genes": [],
                    "metadata": [],
                    "series_matrices": [],
                },
            }
            self.assertTrue(_matrix_header_looks_single_cell(matrix))
            updated = _refresh_manifest_mode(manifest, base)
            self.assertEqual(updated["mode"], "single_cell")


class TestExpandArchiveFiles(unittest.TestCase):
    def test_extracts_feature_bc_matrix_archive_into_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            archive = base / "GSE254513_CG_raw_feature_bc_matrix.tar.gz"
            members = {
                "matrices/S1/raw_feature_bc_matrix/barcodes.tsv": "AAAC-1\n",
                "matrices/S1/raw_feature_bc_matrix/features.tsv": "GENE1\n",
                "matrices/S1/raw_feature_bc_matrix/matrix.mtx": (
                    "%%MatrixMarket matrix coordinate integer general\n"
                    "1 1 1\n"
                    "1 1 1\n"
                ),
            }
            with tarfile.open(archive, "w:gz") as tf:
                for rel, data in members.items():
                    info = tarfile.TarInfo(rel)
                    info.size = len(data)
                    tf.addfile(info, io.BytesIO(data.encode("utf-8")))

            files = {
                "matrix": [archive.name],
                "barcodes": [],
                "genes": [],
                "metadata": [],
                "series_matrices": [],
            }
            logs = []
            expanded = _expand_archive_files(files, base, logs.append)

            prefix = "_extracted/matrices/S1/raw_feature_bc_matrix"
            self.assertEqual(
                expanded["matrix"],
                [f"{prefix}/matrix.mtx"],
            )
            self.assertEqual(
                expanded["barcodes"],
                [f"{prefix}/barcodes.tsv"],
            )
            self.assertEqual(
                expanded["genes"],
                [f"{prefix}/features.tsv"],
            )
            self.assertTrue((base / prefix / "matrix.mtx").exists())
            self.assertEqual(logs, ["expanding archive " + archive.name])


class TestDownloadRetry(unittest.TestCase):
    def test_h5ad_conversion_preserves_extracted_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp)
            (raw / "_extracted").mkdir()
            (raw / "_extracted" / "GSM9440169_E12_B.h5ad").write_bytes(b"x")
            converted = {
                "matrix": "GSM9440169_E12_B.matrix.mtx.gz",
                "barcodes": "GSM9440169_E12_B.barcodes.tsv.gz",
                "genes": "GSM9440169_E12_B.genes.tsv.gz",
            }
            with mock.patch(
                "data.h5_converter.convert_h5ad",
                return_value=converted,
            ):
                downloaded = _convert_downloaded(
                    {
                        "matrix": ["_extracted/GSM9440169_E12_B.h5ad"],
                        "barcodes": [],
                        "genes": [],
                        "metadata": [],
                    },
                    raw,
                )
            self.assertEqual(
                downloaded["matrix"],
                ["_extracted/GSM9440169_E12_B.matrix.mtx.gz"],
            )
            self.assertEqual(
                downloaded["barcodes"],
                ["_extracted/GSM9440169_E12_B.barcodes.tsv.gz"],
            )
            self.assertEqual(
                downloaded["genes"],
                ["_extracted/GSM9440169_E12_B.genes.tsv.gz"],
            )

    def test_skips_existing_nonempty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "GSE123456_RAW.tar"
            out.write_bytes(b"complete")
            log_lines = []
            with mock.patch(
                "data.geo_downloader.subprocess.run",
                side_effect=AssertionError("download should be skipped"),
            ):
                _download(
                    "https://example.test/GSE123456_RAW.tar",
                    out,
                    log_lines.append,
                )
            self.assertEqual(out.read_bytes(), b"complete")
            self.assertTrue(any("already exists" in line for line in log_lines))

    def test_retries_transient_curl_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "GSE123456_RAW.tar"
            calls = {"count": 0}

            def fake_run(cmd, **kwargs):
                calls["count"] += 1
                if calls["count"] < 3:
                    raise subprocess.CalledProcessError(56, cmd)
                out.write_bytes(b"downloaded")
                return subprocess.CompletedProcess(cmd, 0)

            log_lines = []
            with (
                mock.patch("data.geo_downloader._curl", return_value="curl"),
                mock.patch(
                    "data.geo_downloader.subprocess.run",
                    side_effect=fake_run,
                ),
                mock.patch("data.geo_downloader.time.sleep"),
            ):
                _download(
                    "https://example.test/GSE123456_RAW.tar",
                    out,
                    log_lines.append,
                )

            self.assertEqual(calls["count"], 3)
            self.assertEqual(out.read_bytes(), b"downloaded")
            self.assertEqual(len(log_lines), 5)


if __name__ == "__main__":
    unittest.main()
