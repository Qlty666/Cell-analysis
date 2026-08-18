#!/usr/bin/env python3
"""Unit tests for h5ad/loom to 10x MTX conversion."""

from __future__ import annotations

import gzip
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from data.h5_converter import convert_h5ad, convert_loom  # noqa: E402


def _read_gzip(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return fh.read()


def _write_sparse_h5ad(
    path: Path,
    with_names: bool = True,
    with_counts_layer: bool = False,
) -> None:
    with h5py.File(path, "w") as f:
        x = f.create_group("X")
        x.create_dataset("data", data=np.array([1.0, 2.0], dtype=np.float32))
        x.create_dataset("indices", data=np.array([0, 1], dtype=np.int64))
        x.create_dataset("indptr", data=np.array([0, 1, 2], dtype=np.int64))
        x.attrs["shape"] = np.array([2, 2], dtype=np.int64)
        if with_counts_layer:
            layers = f.create_group("layers")
            counts = layers.create_group("counts")
            counts.create_dataset(
                "data", data=np.array([5.0, 6.0], dtype=np.float32)
            )
            counts.create_dataset("indices", data=np.array([0, 1], dtype=np.int64))
            counts.create_dataset("indptr", data=np.array([0, 1, 2], dtype=np.int64))
            counts.attrs["shape"] = np.array([2, 2], dtype=np.int64)
        obs = f.create_group("obs")
        var = f.create_group("var")
        if with_names:
            obs.create_dataset(
                "_index",
                data=np.array([b"AAAC-1", b"TTTG-1"], dtype="S6"),
            )
            var.create_dataset(
                "_index",
                data=np.array([b"TP53", b"EGFR"], dtype="S4"),
            )


class TestH5adConversion(unittest.TestCase):
    def test_decodes_names_and_writes_barcodes_genes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.h5ad"
            _write_sparse_h5ad(path, with_names=True)
            result = convert_h5ad(path)
            barcodes = _read_gzip(path.with_name(result["barcodes"]))
            genes = _read_gzip(path.with_name(result["genes"]))
            self.assertEqual(barcodes.splitlines(), ["AAAC-1", "TTTG-1"])
            self.assertEqual(genes.splitlines(), ["TP53", "EGFR"])

    def test_falls_back_to_generated_names_when_obs_var_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.h5ad"
            _write_sparse_h5ad(path, with_names=False)
            result = convert_h5ad(path)
            barcodes = _read_gzip(path.with_name(result["barcodes"]))
            genes = _read_gzip(path.with_name(result["genes"]))
            self.assertEqual(barcodes.splitlines(), ["Cell1", "Cell2"])
            self.assertEqual(genes.splitlines(), ["Gene1", "Gene2"])

    def test_prefers_counts_layer_over_normalized_x(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.h5ad"
            _write_sparse_h5ad(
                path,
                with_names=True,
                with_counts_layer=True,
            )
            result = convert_h5ad(path)
            matrix_text = _read_gzip(path.with_name(result["matrix"]))
            self.assertIn(" 5", matrix_text)
            self.assertIn(" 6", matrix_text)


class TestLoomConversion(unittest.TestCase):
    def test_decodes_names_and_falls_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.loom"
            with h5py.File(path, "w") as f:
                f.create_dataset(
                    "matrix",
                    data=np.array([[1.0, 2.0], [3.0, 4.0]]),
                )
                row = f.create_group("row_attrs")
                row.create_dataset(
                    "Gene",
                    data=np.array([b"GENE1", b"GENE2"], dtype="S5"),
                )
                col = f.create_group("col_attrs")
                col.create_dataset(
                    "CellID",
                    data=np.array([b"CELL-1", b"CELL-2"], dtype="S6"),
                )
            result = convert_loom(path)
            barcodes = _read_gzip(path.with_name(result["barcodes"]))
            genes = _read_gzip(path.with_name(result["genes"]))
            self.assertEqual(barcodes.splitlines(), ["CELL-1", "CELL-2"])
            self.assertEqual(genes.splitlines(), ["GENE1", "GENE2"])


if __name__ == "__main__":
    unittest.main()
