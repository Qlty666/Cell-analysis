"""Convert h5ad and loom matrices to 10x MTX files."""

import gzip
import shutil
import tempfile
import warnings
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

# Layers (in preference order) that are expected to hold raw integer counts.
COUNTS_LAYER_NAMES = ("counts", "raw_counts", "count")

# AnnData stores sparse matrices with this attribute. Files written by older
# tools omit it; CSR is the AnnData default in that case.
SPARSE_ENCODING_TYPES = {
    "csr_matrix": sp.csr_matrix,
    "csc_matrix": sp.csc_matrix,
}


def _as_integer_counts(matrix, source: str) -> sp.spmatrix:
    """Return ``matrix`` with integer count values or raise ValueError.

    MTX "integer general" output cannot represent log-normalized or scaled
    values, so reject them instead of silently truncating every value to 0/1.
    """
    matrix = sp.csr_matrix(matrix)
    data = np.asarray(matrix.data, dtype=float)
    if data.size == 0:
        return matrix
    if not np.all(np.isfinite(data)):
        raise ValueError(
            f"{source} contains non-finite values; expected raw integer counts"
        )
    rounded = np.rint(data)
    if (rounded < 0).any() or not np.allclose(data, rounded, rtol=0.0, atol=1e-6):
        raise ValueError(
            f"{source} does not contain non-negative integer counts; it looks "
            "normalized or scaled. Provide a raw counts layer (layers['counts'] "
            "or raw/X) and retry."
        )
    matrix = matrix.copy()
    matrix.data = rounded
    return matrix


def _is_integer_counts(matrix) -> bool:
    """Return True when every stored value is a non-negative integer."""
    data = np.asarray(sp.csr_matrix(matrix).data, dtype=float)
    if data.size == 0:
        return True
    if not np.all(np.isfinite(data)):
        return False
    if (data < 0).any():
        return False
    return bool(np.allclose(data, np.rint(data), rtol=0.0, atol=1e-6))


def _write_mtx(path: Path, matrix: sp.spmatrix, genes, cells) -> None:
    matrix = _as_integer_counts(matrix, path.name).tocoo()
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write("%%MatrixMarket matrix coordinate integer general\n")
        fh.write(f"{matrix.shape[0]} {matrix.shape[1]} {matrix.nnz}\n")
        for i, j, v in zip(matrix.row, matrix.col, matrix.data):
            fh.write(f"{i + 1} {j + 1} {int(v)}\n")
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


def _fallback_names(prefix: str, count: int, warnings_out: list[str] | None = None) -> list[str]:
    message = (
        f"missing {prefix.lower()} names in the input file; generated "
        f"placeholder names {prefix}1..{prefix}{count}"
    )
    warnings.warn(message, RuntimeWarning, stacklevel=2)
    if warnings_out is not None:
        warnings_out.append(message)
    return [f"{prefix}{i + 1}" for i in range(count)]


def _load_matrix(x) -> sp.spmatrix:
    """Load an h5ad/loom matrix node honoring AnnData's ``encoding-type``.

    Sparse groups carry a ``csc_matrix`` or ``csr_matrix`` encoding-type; a
    CSC group must not be read with the CSR ``indptr`` interpretation. Plain
    datasets are dense arrays.
    """
    if isinstance(x, h5py.Group):
        encoding = x.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode("utf-8", "replace")
        encoding = str(encoding).strip().lower()
        constructor = SPARSE_ENCODING_TYPES.get(encoding)
        if constructor is None:
            if encoding:
                raise ValueError(f"unsupported sparse encoding-type: {encoding!r}")
            constructor = sp.csr_matrix
        data = x["data"][:]
        indices = x["indices"][:]
        indptr = x["indptr"][:]
        shape = tuple(int(v) for v in x.attrs["shape"])
        return constructor((data, indices, indptr), shape=shape)
    return sp.csr_matrix(np.asarray(x))


def _count_matrix_candidates(f) -> list[tuple[str, object]]:
    """Return ``(label, node)`` count-matrix candidates in preference order."""
    candidates: list[tuple[str, object]] = []
    layers = f.get("layers")
    if isinstance(layers, h5py.Group):
        for name in COUNTS_LAYER_NAMES:
            if name in layers:
                candidates.append((f"layers/{name}", layers[name]))
    raw = f.get("raw")
    if isinstance(raw, h5py.Group):
        for name in ("X",) + COUNTS_LAYER_NAMES:
            if name in raw:
                candidates.append((f"raw/{name}", raw[name]))
    x = f.get("X")
    if x is not None:
        candidates.append(("X", x))
    return candidates


def _select_count_matrix(f) -> tuple[sp.spmatrix, str]:
    """Pick the first candidate layer holding non-negative integer counts.

    ``X`` frequently stores log-normalized values, so it is only accepted when
    its values are integers; ``layers['counts']`` and ``raw/X`` are preferred.
    """
    checked: list[str] = []
    for label, node in _count_matrix_candidates(f):
        matrix = _load_matrix(node)
        if _is_integer_counts(matrix):
            return matrix, label
        checked.append(label)
    raise ValueError(
        "no integer count matrix found in h5ad file; checked "
        f"{', '.join(checked) if checked else 'no count candidates'}. "
        "The available matrices look normalized/scaled; supply a file with "
        "layers['counts'] or raw/X counts."
    )


def convert_h5ad(path: Path) -> dict:
    warnings_out: list[str] = []
    with h5py.File(path, "r") as f:
        matrix, _layer = _select_count_matrix(f)

        obs = f.get("obs")
        var = f.get("var")
        cells = _read_text(
            _find_dataset(
                obs,
                ["_index", "index", "cell_id", "barcode", "cell_names", "obs_names"],
            )
        )
        if not cells:
            cells = _fallback_names("Cell", matrix.shape[1], warnings_out)
        genes = _read_text(
            _find_dataset(
                var,
                ["_index", "index", "gene_names", "gene_ids", "feature_name", "name", "genes"],
            )
        )
        if not genes:
            genes = _fallback_names("Gene", matrix.shape[0], warnings_out)

    prefix = path.with_name(path.stem + ".matrix.mtx.gz")
    _write_mtx(prefix, matrix, genes, cells)
    return {
        "matrix": prefix.name,
        "barcodes": prefix.name.replace("matrix.mtx", "barcodes.tsv"),
        "genes": prefix.name.replace("matrix.mtx", "genes.tsv"),
        "warnings": warnings_out,
    }


def convert_h5ad_gz(path: Path) -> dict:
    """Convert a gzipped .h5ad file by decompressing to a temporary file.

    ``convert_h5ad`` writes its outputs next to its input, so the converted
    files are moved out of the temporary directory before it is removed.
    """
    path = Path(path)
    with tempfile.TemporaryDirectory(prefix="h5ad_gz_") as tmp:
        plain = Path(tmp) / path.stem
        with gzip.open(path, "rb") as src, plain.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        result = convert_h5ad(plain)
        moved: dict = {"warnings": result.get("warnings", [])}
        for key in ("matrix", "barcodes", "genes"):
            name = result[key]
            shutil.move(str(Path(tmp) / name), str(path.parent / name))
            moved[key] = name
    return moved


def convert_loom(path: Path) -> dict:
    warnings_out: list[str] = []
    with h5py.File(path, "r") as f:
        matrix = sp.csr_matrix(np.asarray(f["matrix"][:]))
        row = f.get("row_attrs")
        col = f.get("col_attrs")
        genes = _read_text(_find_dataset(row, ["Gene", "gene_names"]))
        cells = _read_text(_find_dataset(col, ["CellID", "cell_names"]))
        if not genes:
            genes = _fallback_names("Gene", matrix.shape[0], warnings_out)
        if not cells:
            cells = _fallback_names("Cell", matrix.shape[1], warnings_out)

    prefix = path.with_name(path.stem + ".matrix.mtx.gz")
    _write_mtx(prefix, matrix, genes, cells)
    return {
        "matrix": prefix.name,
        "barcodes": prefix.name.replace("matrix.mtx", "barcodes.tsv"),
        "genes": prefix.name.replace("matrix.mtx", "genes.tsv"),
        "warnings": warnings_out,
    }
