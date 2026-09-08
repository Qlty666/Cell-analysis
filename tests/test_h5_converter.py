#!/usr/bin/env python3
"""Unit tests for h5ad/loom to 10x MTX conversion."""

from __future__ import annotations

import gzip
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from data.h5_converter import (  # noqa: E402
    convert_h5ad,
    convert_h5ad_gz,
    convert_loom,
)


def _read_gzip(path: Path) -> str:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return fh.read()


def _mtx_entries(path: Path) -> dict[tuple[int, int], int]:
    """Parse a MatrixMarket file into {(row, col): value} (1-based)."""
    entries: dict[tuple[int, int], int] = {}
    for line in _read_gzip(path).splitlines()[2:]:
        row, col, value = line.split()
        entries[(int(row), int(col))] = int(value)
    return entries


def _create_sparse_group(
    parent,
    name: str,
    data,
    indices,
    indptr,
    shape,
    encoding_type: str | None = None,
):
    group = parent.create_group(name)
    group.create_dataset("data", data=np.array(data, dtype=np.float32))
    group.create_dataset("indices", data=np.array(indices, dtype=np.int64))
    group.create_dataset("indptr", data=np.array(indptr, dtype=np.int64))
    group.attrs["shape"] = np.array(shape, dtype=np.int64)
    if encoding_type:
        group.attrs["encoding-type"] = encoding_type
    return group


def _write_sparse_h5ad(
    path: Path,
    with_names: bool = True,
    with_counts_layer: bool = False,
    with_attrs_index: bool = False,
    x_data=(1.0, 2.0),
    with_raw_counts: bool = False,
) -> None:
    with h5py.File(path, "w") as f:
        _create_sparse_group(f, "X", x_data, [0, 1], [0, 1, 2], (2, 2))
        if with_counts_layer:
            _create_sparse_group(
                f.create_group("layers"),
                "counts",
                [5.0, 6.0],
                [0, 1],
                [0, 1, 2],
                (2, 2),
            )
        if with_raw_counts:
            _create_sparse_group(
                f.create_group("raw"),
                "X",
                [7.0, 8.0],
                [0, 1],
                [0, 1, 2],
                (2, 2),
            )
        obs = f.create_group("obs")
        var = f.create_group("var")
        if with_attrs_index:
            obs.create_dataset(
                "sample_id",
                data=np.array([b"AAAC-1", b"TTTG-1"], dtype="S6"),
            )
            obs.attrs["_index"] = "sample_id"
            var.create_dataset(
                "gene_id",
                data=np.array([b"TP53", b"EGFR"], dtype="S4"),
            )
            var.attrs["_index"] = "gene_id"
        elif with_names:
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

    def test_fallback_names_are_warned_and_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.h5ad"
            _write_sparse_h5ad(path, with_names=False)
            with self.assertWarns(RuntimeWarning):
                result = convert_h5ad(path)
            self.assertEqual(len(result["warnings"]), 2)
            self.assertTrue(
                any("Cell1" in message for message in result["warnings"]),
                result["warnings"],
            )

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

    def test_uses_obs_var_attrs_index_when_underscore_index_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.h5ad"
            _write_sparse_h5ad(path, with_attrs_index=True)
            result = convert_h5ad(path)
            barcodes = _read_gzip(path.with_name(result["barcodes"]))
            genes = _read_gzip(path.with_name(result["genes"]))
            self.assertEqual(barcodes.splitlines(), ["AAAC-1", "TTTG-1"])
            self.assertEqual(genes.splitlines(), ["TP53", "EGFR"])

    def test_reads_csc_encoded_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "csc.h5ad"
            with h5py.File(path, "w") as f:
                _create_sparse_group(
                    f,
                    "X",
                    data=[1.0, 3.0, 2.0],
                    indices=[0, 2, 1],
                    indptr=[0, 2, 3],
                    shape=(3, 2),
                    encoding_type="csc_matrix",
                )
            result = convert_h5ad(path)
            self.assertEqual(
                _mtx_entries(path.parent / result["matrix"]),
                {(1, 1): 1, (3, 1): 3, (2, 2): 2},
            )

    def test_reads_dense_x_dataset(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "dense.h5ad"
            with h5py.File(path, "w") as f:
                f.create_dataset(
                    "X",
                    data=np.array([[1, 2], [3, 4]], dtype=np.int32),
                )
            result = convert_h5ad(path)
            self.assertEqual(
                _mtx_entries(path.parent / result["matrix"]),
                {(1, 1): 1, (1, 2): 2, (2, 1): 3, (2, 2): 4},
            )

    def test_rejects_normalized_x_without_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "normalized.h5ad"
            _write_sparse_h5ad(path, with_names=True, x_data=(0.5, 1.5))
            with self.assertRaisesRegex(ValueError, "integer count matrix"):
                convert_h5ad(path)

    def test_prefers_raw_counts_over_normalized_x(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "raw.h5ad"
            _write_sparse_h5ad(
                path,
                with_names=True,
                x_data=(0.5, 1.5),
                with_raw_counts=True,
            )
            result = convert_h5ad(path)
            self.assertEqual(
                _mtx_entries(path.parent / result["matrix"]),
                {(1, 1): 7, (2, 2): 8},
            )

    def test_gz_conversion_leaves_outputs_next_to_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            plain = Path(tmp) / "sample.h5ad"
            _write_sparse_h5ad(plain, with_names=True)
            gz_path = Path(tmp) / "sample.h5ad.gz"
            with plain.open("rb") as src, gzip.open(gz_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
            result = convert_h5ad_gz(gz_path)
            for key in ("matrix", "barcodes", "genes"):
                out = gz_path.parent / result[key]
                self.assertTrue(out.exists(), f"{key} output missing: {out}")
            self.assertEqual(
                _read_gzip(gz_path.parent / result["barcodes"]).splitlines(),
                ["AAAC-1", "TTTG-1"],
            )
            self.assertEqual(
                _mtx_entries(gz_path.parent / result["matrix"]),
                {(1, 1): 1, (2, 2): 2},
            )


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

    def test_loom_fallback_names_are_warned(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bare.loom"
            with h5py.File(path, "w") as f:
                f.create_dataset("matrix", data=np.array([[1.0, 2.0]]))
            with self.assertWarns(RuntimeWarning):
                result = convert_loom(path)
            self.assertTrue(result["warnings"])

    def test_loom_rejects_non_integer_matrix(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "normalized.loom"
            with h5py.File(path, "w") as f:
                f.create_dataset("matrix", data=np.array([[0.5, 1.5]]))
            with self.assertRaisesRegex(ValueError, "integer counts"):
                convert_loom(path)


if __name__ == "__main__":
    unittest.main()
