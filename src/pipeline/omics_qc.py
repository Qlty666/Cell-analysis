"""Quality-control summary for the pseudobulk expression matrix."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from docking.utils import write_json


def _load_expression(path: Path) -> tuple[pd.DataFrame, list[str]]:
    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError(f"expression matrix is empty: {path}")
    if {"gene", "sample", "value"}.issubset(frame.columns):
        samples = sorted(frame["sample"].astype(str).unique())
        matrix = (
            frame.pivot_table(
                index="gene",
                columns="sample",
                values="value",
                aggfunc="mean",
            )
            .sort_index()
        )
        matrix.index = matrix.index.astype(str)
        return matrix, samples
    if "gene" not in frame.columns:
        raise ValueError("expression matrix requires a gene column")
    matrix = frame.set_index("gene")
    matrix.index = matrix.index.astype(str)
    numeric = matrix.apply(pd.to_numeric, errors="coerce")
    return numeric, [str(column) for column in numeric.columns]


def run_omics_qc(workdir: Path, out_dir: Path) -> dict:
    """Assess matrix integrity and sample metadata before downstream ranking."""
    expression_path = workdir / "data" / "knockout" / "expression.csv"
    metadata_path = workdir / "data" / "knockout" / "metadata.csv"
    if not expression_path.exists():
        summary = {
            "status": "skipped",
            "reason": "pseudobulk expression matrix not found",
            "gate_passed": False,
        }
        write_json(out_dir / "omics_qc_summary.json", summary)
        return summary
    matrix, samples = _load_expression(expression_path)
    duplicate_genes = int(matrix.index.duplicated().sum())
    matrix = matrix[~matrix.index.duplicated(keep="first")]
    values = matrix.to_numpy(dtype=float)
    finite = np.isfinite(values)
    total_values = int(values.size)
    missing_fraction = (
        float((~finite).sum() / total_values) if total_values else 1.0
    )
    zero_fraction = (
        float((values[finite] == 0).sum() / finite.sum())
        if finite.any()
        else 0.0
    )
    negative_values = int((values[finite] < 0).sum())

    metadata_match = False
    group_counts: dict[str, int] = {}
    has_cell_type = False
    has_batch = False
    metadata = pd.DataFrame()
    if metadata_path.exists():
        metadata = pd.read_csv(metadata_path)
        sample_column = next(
            (
                column
                for column in ("sample", "sample_id", "barcode", "cell")
                if column in metadata.columns
            ),
            None,
        )
        if sample_column:
            metadata_samples = set(metadata[sample_column].astype(str))
            metadata_match = set(samples).issubset(metadata_samples)
        group_column = next(
            (
                column
                for column in ("condition", "group", "label", "status")
                if column in metadata.columns
            ),
            None,
        )
        if group_column:
            group_counts = {
                str(key): int(value)
                for key, value in metadata[group_column]
                .astype(str)
                .value_counts()
                .items()
            }
        has_cell_type = "cell_type" in metadata.columns
        has_batch = any(
            column in metadata.columns
            for column in ("batch", "dataset", "study", "source")
        )

    gates = {
        "sample_metadata_match": bool(metadata_match),
        "no_duplicate_genes": duplicate_genes == 0,
        "acceptable_missingness": missing_fraction <= 0.20,
        "nonnegative_expression": negative_values == 0,
        "at_least_two_groups": len(group_counts) >= 2,
        "each_group_at_least_two_samples": bool(
            group_counts and min(group_counts.values()) >= 2
        ),
    }
    summary = {
        "status": "completed",
        "expression_path": str(expression_path),
        "metadata_path": str(metadata_path) if metadata_path.exists() else "",
        "n_genes": int(matrix.shape[0]),
        "n_samples": int(matrix.shape[1]),
        "duplicate_genes": duplicate_genes,
        "missing_fraction": missing_fraction,
        "zero_fraction": zero_fraction,
        "negative_values": negative_values,
        "sample_metadata_match": metadata_match,
        "group_counts": group_counts,
        "has_cell_type": has_cell_type,
        "has_batch_metadata": has_batch,
        "gates": gates,
        "gate_passed": all(gates.values()),
    }
    sample_metrics = pd.DataFrame(
        {
            "sample": matrix.columns.astype(str),
            "missing_fraction": np.isnan(values).mean(axis=0),
            "mean_expression": np.nanmean(values, axis=0),
            "median_expression": np.nanmedian(values, axis=0),
            "zero_fraction": np.nanmean(values == 0, axis=0),
        }
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_metrics.to_csv(out_dir / "omics_qc_sample_metrics.csv", index=False)
    write_json(out_dir / "omics_qc_summary.json", summary)
    return summary
