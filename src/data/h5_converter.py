"""Convert h5ad and loom matrices to 10x MTX files."""

import gzip
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
            fh.write(f"{i + 1} {j + 1} {int(v)}\n")
    with gzip.open(path.with_name(path.name.replace("matrix.mtx", "barcodes.tsv")), "wt", encoding="utf-8") as fh:
        fh.write("\n".join(cells) + "\n")
    with gzip.open(path.with_name(path.name.replace("matrix.mtx", "genes.tsv")), "wt", encoding="utf-8") as fh:
        fh.write("\n".join(genes) + "\n")


def _read_text(ds) -> list[str]:
    return [str(x) for x in np.asarray(ds).reshape(-1)]


def convert_h5ad(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        x = f["X"]
        if isinstance(x, h5py.Group):
            data = x["data"][:]
            indices = x["indices"][:]
            indptr = x["indptr"][:]
            shape = tuple(int(v) for v in x.attrs["shape"])
            matrix = sp.csr_matrix((data, indices, indptr), shape=shape)
        else:
            matrix = sp.csr_matrix(np.asarray(x))

        obs = f["obs"]
        var = f["var"]
        cells = _read_text(obs["_index"]) if "_index" in obs else []
        if var and "_index" in var:
            genes = _read_text(var["_index"])
        else:
            genes = [f"Gene{i}" for i in range(matrix.shape[0])]

    prefix = path.with_name(path.stem + ".matrix.mtx.gz")
    _write_mtx(prefix, matrix, genes, cells)
    return {
        "matrix": prefix.name,
        "barcodes": prefix.name.replace("matrix.mtx", "barcodes.tsv"),
        "genes": prefix.name.replace("matrix.mtx", "genes.tsv"),
    }


def convert_loom(path: Path) -> dict:
    with h5py.File(path, "r") as f:
        matrix = sp.csr_matrix(np.asarray(f["matrix"][:]))
        row = f["row_attrs"]
        col = f["col_attrs"]
        genes = _read_text(row["Gene"]) if "Gene" in row else _read_text(row["gene_names"])
        cells = _read_text(col["CellID"]) if "CellID" in col else _read_text(col["cell_names"])

    prefix = path.with_name(path.stem + ".matrix.mtx.gz")
    _write_mtx(prefix, matrix, genes, cells)
    return {
        "matrix": prefix.name,
        "barcodes": prefix.name.replace("matrix.mtx", "barcodes.tsv"),
        "genes": prefix.name.replace("matrix.mtx", "genes.tsv"),
    }
