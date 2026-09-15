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


def _finite_matrix(matrix: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    values = matrix.to_numpy(dtype=float)
    return values, np.isfinite(values)


def _pca_scores(
    matrix: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return sample PCA coordinates, explained variance and retained genes."""
    samples_by_gene = matrix.T
    numeric = samples_by_gene.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.replace([np.inf, -np.inf], np.nan)
    keep = numeric.notna().sum(axis=0) >= 2
    numeric = numeric.loc[:, keep]
    if numeric.shape[1] == 0 or numeric.shape[0] < 2:
        return (
            np.empty((numeric.shape[0], 0), dtype=float),
            np.array([], dtype=float),
            [],
        )
    medians = numeric.median(axis=0, skipna=True)
    imputed = numeric.fillna(medians)
    variances = imputed.var(axis=0, skipna=True)
    retained = variances[variances > 0].index.tolist()
    if not retained:
        return (
            np.empty((numeric.shape[0], 0), dtype=float),
            np.array([], dtype=float),
            [],
        )
    x = imputed[retained].to_numpy(dtype=float)
    x = (x - np.nanmean(x, axis=0)) / np.where(
        np.nanstd(x, axis=0) == 0,
        1.0,
        np.nanstd(x, axis=0),
    )
    try:
        from sklearn.decomposition import PCA

        components = min(2, x.shape[0] - 1, x.shape[1])
        pca = PCA(n_components=components, random_state=42)
        scores = pca.fit_transform(x)
        if scores.shape[1] < 2:
            scores = np.column_stack(
                [scores[:, 0], np.zeros(scores.shape[0])]
            )
        return scores, pca.explained_variance_ratio_, retained
    except Exception:
        centered = x - np.nanmean(x, axis=0)
        _, singular_values, vt = np.linalg.svd(
            centered,
            full_matrices=False,
        )
        scores = centered @ vt[:2].T
        if scores.shape[1] < 2:
            scores = np.column_stack(
                [scores[:, 0], np.zeros(scores.shape[0])]
            )
        variance = singular_values**2
        explained = (
            variance[:2] / variance.sum()
            if variance.sum() > 0
            else np.zeros(2)
        )
        return scores, explained, retained


def _detect_pca_outliers(
    pca_scores: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    n_samples = pca_scores.shape[0]
    if n_samples == 0 or pca_scores.shape[1] == 0:
        return np.zeros(n_samples, dtype=bool), np.zeros(n_samples, dtype=bool)
    isolation = np.zeros(n_samples, dtype=bool)
    if n_samples >= 5:
        try:
            from sklearn.ensemble import IsolationForest

            model = IsolationForest(
                n_estimators=200,
                random_state=42,
                contamination="auto",
            )
            isolation = model.fit_predict(pca_scores[:, :2]) == -1
        except Exception:
            isolation = np.zeros(n_samples, dtype=bool)
    robust = np.zeros(n_samples, dtype=bool)
    if n_samples >= 4:
        for column in range(min(2, pca_scores.shape[1])):
            values = pca_scores[:, column]
            median = float(np.median(values))
            mad = float(np.median(np.abs(values - median)))
            scale = max(1.4826 * mad, 1e-12)
            robust |= np.abs(values - median) / scale > 3.5
    return isolation, robust


def _correlation_values(
    matrix: pd.DataFrame,
) -> tuple[pd.DataFrame, float | None, float | None, list[dict[str, str]]]:
    values = matrix.to_numpy(dtype=float)
    if values.shape[1] < 2:
        return pd.DataFrame(), None, None, []
    # Correlations are calculated on finite observations only for each pair.
    rows: list[dict[str, float | str]] = []
    pair_values_by_sample: dict[str, list[float]] = {
        str(sample): [] for sample in matrix.columns
    }
    for left in range(values.shape[1]):
        for right in range(left + 1, values.shape[1]):
            pair = values[:, [left, right]]
            finite = np.isfinite(pair).all(axis=1)
            if finite.sum() < 2:
                correlation = np.nan
            else:
                x = pair[finite, 0]
                y = pair[finite, 1]
                if np.std(x) == 0 or np.std(y) == 0:
                    correlation = np.nan
                else:
                    correlation = float(np.corrcoef(x, y)[0, 1])
            left_name = str(matrix.columns[left])
            right_name = str(matrix.columns[right])
            rows.append(
                {
                    "sample_a": left_name,
                    "sample_b": right_name,
                    "correlation": correlation,
                }
            )
            if np.isfinite(correlation):
                pair_values_by_sample[left_name].append(correlation)
                pair_values_by_sample[right_name].append(correlation)
    correlation_frame = pd.DataFrame(rows)
    valid = correlation_frame["correlation"].dropna()
    correlation_median = float(valid.median()) if not valid.empty else None
    correlation_min = float(valid.min()) if not valid.empty else None
    low_pairs = (
        correlation_frame.loc[
            correlation_frame["correlation"] < 0.50,
            ["sample_a", "sample_b"],
        ]
        .astype(str)
        .to_dict("records")
        if "correlation" in correlation_frame
        else []
    )
    return correlation_frame, correlation_median, correlation_min, low_pairs


def _batch_metrics(
    metadata: pd.DataFrame,
    samples: list[str],
    pca_scores: np.ndarray,
) -> dict:
    batch_column = next(
        (
            column
            for column in ("batch", "dataset", "study", "source")
            if column in metadata.columns
        ),
        None,
    )
    condition_column = next(
        (
            column
            for column in ("condition", "group", "label", "status")
            if column in metadata.columns
        ),
        None,
    )
    sample_column = next(
        (
            column
            for column in ("sample", "sample_id", "barcode", "cell")
            if column in metadata.columns
        ),
        None,
    )
    result = {
        "batch_column": batch_column or "",
        "batch_silhouette_score": None,
        "batch_condition_confounding": False,
        "batch_condition_association_p": None,
        "batch_counts": {},
    }
    if (
        batch_column is None
        or sample_column is None
        or pca_scores.shape[0] < 2
    ):
        return result
    batch_by_sample = (
        metadata.assign(
            __sample=metadata[sample_column].astype(str)
        )
        .drop_duplicates("__sample")
        .set_index("__sample")[batch_column]
        .astype(str)
    )
    labels = [batch_by_sample.get(str(sample), "") for sample in samples]
    result["batch_counts"] = {
        str(key): int(value)
        for key, value in pd.Series(labels).value_counts().items()
        if str(key)
    }
    if len(set(labels)) >= 2 and all(labels):
        try:
            from sklearn.metrics import silhouette_score

            result["batch_silhouette_score"] = float(
                silhouette_score(pca_scores[:, :2], labels)
            )
        except Exception:
            result["batch_silhouette_score"] = None
    if condition_column is None or len(set(labels)) < 2:
        return result
    condition_by_sample = (
        metadata.assign(
            __sample=metadata[sample_column].astype(str)
        )
        .drop_duplicates("__sample")
        .set_index("__sample")[condition_column]
        .astype(str)
    )
    conditions = [
        condition_by_sample.get(str(sample), "")
        for sample in samples
    ]
    cross_tab = pd.crosstab(
        pd.Series(labels, name="batch"),
        pd.Series(conditions, name="condition"),
    )
    if cross_tab.shape[0] > 1 and cross_tab.shape[1] > 1:
        row_totals = cross_tab.sum(axis=1).replace(0, np.nan)
        row_shares = cross_tab.div(row_totals, axis=0).fillna(0.0)
        result["batch_condition_confounding"] = bool(
            (row_shares.max(axis=1) >= 0.80).any()
        )
        try:
            from scipy.stats import chi2_contingency

            _, p_value, _, _ = chi2_contingency(cross_tab)
            result["batch_condition_association_p"] = float(p_value)
        except Exception:
            result["batch_condition_association_p"] = None
    return result


def _metadata_context(
    metadata_path: Path,
    samples: list[str],
) -> tuple[pd.DataFrame, bool, dict[str, int], bool, bool]:
    metadata = pd.DataFrame()
    metadata_match = False
    group_counts: dict[str, int] = {}
    has_cell_type = False
    has_batch = False
    if not metadata_path.exists():
        return (
            metadata,
            metadata_match,
            group_counts,
            has_cell_type,
            has_batch,
        )
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
    return (
        metadata,
        metadata_match,
        group_counts,
        has_cell_type,
        has_batch,
    )


def run_omics_qc(workdir: Path, out_dir: Path) -> dict:
    """Assess matrix integrity, sample outliers and batch metadata."""
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
    values, finite = _finite_matrix(matrix)
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

    (
        metadata,
        metadata_match,
        group_counts,
        has_cell_type,
        has_batch,
    ) = _metadata_context(metadata_path, samples)
    pca_scores, explained_variance, retained_genes = _pca_scores(matrix)
    isolation_outliers, robust_outliers = _detect_pca_outliers(pca_scores)
    outliers = isolation_outliers | robust_outliers
    (
        correlation_frame,
        correlation_median,
        correlation_min,
        low_correlation_pairs,
    ) = _correlation_values(matrix)
    batch_metrics = _batch_metrics(metadata, samples, pca_scores)

    finite_count = finite.sum(axis=0)
    sample_missing = (
        1.0 - finite_count / values.shape[0]
        if values.shape[0]
        else np.zeros(values.shape[1], dtype=float)
    )
    sample_zero = np.asarray(
        [
            float((values[finite[:, index], index] == 0).mean())
            if finite[:, index].any()
            else 0.0
            for index in range(values.shape[1])
        ]
    )
    sample_mean = np.asarray(
        [
            float(np.nanmean(values[:, index]))
            if finite[:, index].any()
            else np.nan
            for index in range(values.shape[1])
        ]
    )
    sample_median = np.asarray(
        [
            float(np.nanmedian(values[:, index]))
            if finite[:, index].any()
            else np.nan
            for index in range(values.shape[1])
        ]
    )
    correlation_by_sample = {}
    if not correlation_frame.empty:
        for row in correlation_frame.to_dict("records"):
            for key in ("sample_a", "sample_b"):
                correlation_by_sample.setdefault(str(row[key]), []).append(
                    row["correlation"]
                )
    sample_correlation_values = []
    for sample in matrix.columns:
        correlation_values = [
            float(value)
            for value in correlation_by_sample.get(str(sample), [])
            if pd.notna(value) and np.isfinite(value)
        ]
        sample_correlation_values.append(
            float(np.mean(correlation_values))
            if correlation_values
            else np.nan
        )
    sample_correlation = np.asarray(sample_correlation_values)
    outlier_fraction = (
        float(outliers.sum() / len(outliers)) if len(outliers) else 0.0
    )
    low_correlation_fraction = (
        len(low_correlation_pairs)
        / max(1, len(correlation_frame))
        if not correlation_frame.empty
        else 0.0
    )
    median_cv = None
    if len(sample_median) >= 3 and np.isfinite(sample_median).any():
        finite_medians = sample_median[np.isfinite(sample_median)]
        if len(finite_medians) >= 3 and abs(float(np.mean(finite_medians))) > 0:
            median_cv = float(
                np.std(finite_medians, ddof=1)
                / abs(float(np.mean(finite_medians)))
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
        "sample_outliers_acceptable": outlier_fraction <= 0.20,
        "sample_correlation_acceptable": bool(
            correlation_median is None or correlation_median >= 0.50
        )
        and low_correlation_fraction <= 0.20,
        "batch_confounding_acceptable": not bool(
            batch_metrics.get("batch_condition_confounding", False)
        ),
        "sample_distribution_acceptable": bool(
            median_cv is None or median_cv <= 0.50
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
        "pca_retained_genes": len(retained_genes),
        "pca_explained_variance": [
            float(value) for value in explained_variance
        ],
        "sample_outlier_count": int(outliers.sum()),
        "sample_outlier_fraction": outlier_fraction,
        "outlier_samples": [
            str(sample)
            for sample, outlier in zip(matrix.columns, outliers)
            if outlier
        ],
        "correlation_median": correlation_median,
        "correlation_min": correlation_min,
        "low_correlation_pairs": low_correlation_pairs,
        "sample_median_cv": median_cv,
        "batch_silhouette_score": batch_metrics.get(
            "batch_silhouette_score"
        ),
        "batch_condition_confounding": batch_metrics.get(
            "batch_condition_confounding",
            False,
        ),
        "batch_condition_association_p": batch_metrics.get(
            "batch_condition_association_p"
        ),
        "batch_counts": batch_metrics.get("batch_counts", {}),
        "gates": gates,
        "gate_passed": all(gates.values()),
    }
    sample_metrics = pd.DataFrame(
        {
            "sample": matrix.columns.astype(str),
            "missing_fraction": sample_missing,
            "zero_fraction": sample_zero,
            "mean_expression": sample_mean,
            "median_expression": sample_median,
            "mean_correlation": sample_correlation,
            "outlier": outliers,
            "pca_outlier": isolation_outliers,
            "robust_outlier": robust_outliers,
        }
    )
    pca_frame = pd.DataFrame(
        {
            "sample": matrix.columns.astype(str),
            "PC1": (
                pca_scores[:, 0]
                if pca_scores.shape[1] > 0
                else np.zeros(len(matrix.columns))
            ),
            "PC2": (
                pca_scores[:, 1]
                if pca_scores.shape[1] > 1
                else np.zeros(len(matrix.columns))
            ),
            "outlier": outliers,
        }
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_metrics.to_csv(out_dir / "omics_qc_sample_metrics.csv", index=False)
    pca_frame.to_csv(out_dir / "omics_qc_pca.csv", index=False)
    correlation_frame.to_csv(
        out_dir / "omics_qc_correlation.csv",
        index=False,
    )
    write_json(out_dir / "omics_qc_summary.json", summary)
    return summary
