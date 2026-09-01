"""Convert h5ad and loom matrices to 10x MTX files."""

import gzip
import shutil
import tempfile
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp


def _write_mtx(path: Path, matrix: sp.spmatrix, genes, cells) -> None:
    matrix = matrix.tocoo()
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{matrix.shape[0]} {matrix.shape[1]} {matrix.nnz}\n")
        for i, j, v in zip(matrix.row, matrix.col, matrix.data):
            try:
                value = int(round(float(v)))
            except (TypeError, ValueError, OverflowError):
                value = 0
            fh.write(f"{i + 1} {j + 1} {value}\n")
    with gzip.open(path.with_name(path.name.replace("matrix.mtx", "barcodes.tsv")), "wt", encoding="utf-8") as fh:
        fh.write("\n".join(cells) + "\n")
    with gzip.open(path.with_name(path.name.replace("matrix.mtx", "genes.tsv")), "wt", encoding="utf-8") as fh:
        fh.write("\n".join(genes) + "\n")


def _decode_value(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8", "replace")
    return str(value)


def _read_text(ds) -> list[str]:
    if ds is None:
        return []
    return [
        text
        for value in np.asarray(ds).reshape(-1)
        if value is not None and (text := _decode_value(value))
    ]


def _find_dataset(group, keys: list[str]):
    if not isinstance(group, h5py.Group):
        return None
    index_name = group.attrs.get("_index")
    if isinstance(index_name, str) and index_name in group:
        item = group[index_name]
        if isinstance(item, h5py.Dataset):
            return item
    for key in keys:
        if key in group:
            item = group[key]
            if isinstance(item, h5py.Dataset):
                return item
    return None


def _fallback_names(prefix: str, count: int) -> list[str]:
    return [f"{prefix}{i + 1}" for i in range(count)]


def _load_count_matrix(x) -> sp.spmatrix:
    if isinstance(x, h5py.Group):
        data = x["data"][:]
        indices = x["indices"][:]
        indptr = x["indptr"][:]
        shape = tuple(int(v) for v in x.attrs["shape"])
        return sp.csr_matrix((data, indices, indptr), shape=shape)
    return sp.csr_matrix(np.asarray(x))


def convert_h5ad(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        x = None
        layers = f.get("layers")
        if isinstance(layers, h5py.Group):
            for name in ("counts", "raw_counts", "count"):
                if name in layers:
                    x = layers[name]
                    break
        if x is None:
            x = f["X"]
        matrix = _load_count_matrix(x)

        obs = f.get("obs")
        var = f.get("var")
        cells = _read_text(
            _find_dataset(
                obs,
                ["_index", "index", "cell_id", "barcode", "cell_names", "obs_names"],
            )
        )
        if not cells:
            cells = _fallback_names("Cell", matrix.shape[1])
        genes = _read_text(
            _find_dataset(
                var,
                ["_index", "index", "gene_names", "gene_ids", "feature_name", "name", "genes"],
            )
        )
        if not genes:
            genes = _fallback_names("Gene", matrix.shape[0])

    prefix = path.with_name(path.stem + ".matrix.mtx.gz")
    _write_mtx(prefix, matrix, genes, cells)
    return {
        "matrix": prefix.name,
        "barcodes": prefix.name.replace("matrix.mtx", "barcodes.tsv"),
        "genes": prefix.name.replace("matrix.mtx", "genes.tsv"),
    }


def convert_h5ad_gz(path: Path) -> dict:
    """Convert a gzipped .h5ad file by decompressing to a temporary file."""
    with tempfile.TemporaryDirectory(prefix="h5ad_gz_") as tmp:
        plain = Path(tmp) / path.stem
        with gzip.open(path, "rb") as src, plain.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        return convert_h5ad(plain)


def convert_loom(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        matrix = sp.csr_matrix(np.asarray(f["matrix"][:]))
        row = f.get("row_attrs")
        col = f.get("col_attrs")
        genes = _read_text(_find_dataset(row, ["Gene", "gene_names"]))
        cells = _read_text(_find_dataset(col, ["CellID", "cell_names"]))
        if not genes:
            genes = _fallback_names("Gene", matrix.shape[0])
        if not cells:
            cells = _fallback_names("Cell", matrix.shape[1])

    prefix = path.with_name(path.stem + ".matrix.mtx.gz")
    _write_mtx(prefix, matrix, genes, cells)
    return {
        "matrix": prefix.name,
        "barcodes": prefix.name.replace("matrix.mtx", "barcodes.tsv"),
        "genes": prefix.name.replace("matrix.mtx", "genes.tsv"),
    }
