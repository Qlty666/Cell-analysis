"""Convert GEO series-matrix expression tables to gene x sample TSV."""

from __future__ import annotations

import csv
import gzip
from pathlib import Path

import pandas as pd


def _open_text(path: Path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _parse_fields(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t"))


def series_matrix_to_tsv(source: str | Path, destination: str | Path) -> Path:
    """Extract the ``!series_matrix_table_begin`` block into a TSV file."""
    source = Path(source)
    destination = Path(destination)
    with _open_text(source) as handle:
        lines = handle.read().splitlines()
    begin = next(
        (
            index
            for index, line in enumerate(lines)
            if line.strip().lower() == "!series_matrix_table_begin"
        ),
        None,
    )
    if begin is None:
        begin = next(
            (
                index
                for index, line in enumerate(lines)
                if line.lstrip("\"'").upper().startswith("ID_REF")
            ),
            None,
        )
        if begin is not None:
            begin -= 1
    if begin is None:
        raise ValueError(f"series matrix table not found: {source}")
    header_index = begin + 1
    if header_index >= len(lines):
        raise ValueError(f"series matrix table is empty: {source}")
    end = len(lines)
    for index in range(header_index + 1, len(lines)):
        if lines[index].strip().lower() == "!series_matrix_table_end":
            end = index
            break
    header = _parse_fields(lines[header_index])
    if len(header) < 2:
        raise ValueError(f"series matrix has fewer than two columns: {source}")
    rows = []
    for line in lines[header_index + 1 : end]:
        if not line.strip() or line.startswith("!"):
            continue
        fields = _parse_fields(line)
        if len(fields) < len(header):
            fields.extend([""] * (len(header) - len(fields)))
        rows.append(fields[: len(header)])
    if not rows:
        raise ValueError(f"series matrix has no expression rows: {source}")
    frame = pd.DataFrame(rows, columns=header)
    first = str(frame.columns[0])
    frame = frame.rename(columns={first: "gene"})
    for column in frame.columns[1:]:
        frame[column] = pd.to_numeric(
            frame[column].replace(
                {"": pd.NA, "null": pd.NA, "NA": pd.NA, "N/A": pd.NA}
            ),
            errors="coerce",
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if str(destination).lower().endswith(".gz") else None
    frame.to_csv(
        destination,
        sep="\t",
        index=False,
        compression=compression,
    )
    return destination
