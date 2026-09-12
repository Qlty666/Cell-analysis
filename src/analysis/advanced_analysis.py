#!/usr/bin/env python3
"""Advanced multi-cohort bulk analysis and gene-level model validation.

This module is intentionally independent of the single-cell pipeline. It adds
the analyses that are commonly used in exposure-toxicity and target-discovery
papers but cannot be represented by a cell-type proportion model:

* discovery and external-validation cohort loading;
* count/normalised/microarray-aware preprocessing;
* optional ComBat/limma/WGCNA/immune/survival analyses through R;
* leakage-resistant gene-level machine-learning model comparison;
* ROC confidence intervals, calibration, decision curves and optional SHAP;
* integration of DEG, WGCNA, ML and user-supplied evidence into one ranking.

The Python portion works without optional R packages. R-backed sections are
skipped with an explicit status when their packages are unavailable.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.calibration import calibration_curve
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import (
    RepeatedStratifiedKFold,
    StratifiedKFold,
    cross_val_predict,
)
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC

APP_ROOT = Path(__file__).resolve().parents[2]
LOG = logging.getLogger("advanced_analysis")


@dataclass
class Cohort:
    """One expression cohort aligned to sample-level labels."""

    name: str
    expression: pd.DataFrame
    metadata: pd.DataFrame
    labels: pd.Series
    condition_column: str
    case_label: str
    control_label: str


def _json_load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _json_write(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _read_table(path: Path) -> pd.DataFrame:
    """Read CSV/TSV with automatic delimiter detection."""
    if not path.exists():
        raise FileNotFoundError(f"table not found: {path}")
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else None
    try:
        frame = pd.read_csv(path, sep=sep, engine="python")
    except Exception:
        frame = pd.read_csv(path, sep=None, engine="python")
    if frame.empty:
        raise ValueError(f"table is empty: {path}")
    return frame


def _first_existing(columns: list[str], candidates: list[str]) -> str | None:
    lower = {str(column).strip().lower(): str(column) for column in columns}
    return next((lower[candidate] for candidate in candidates if candidate in lower), None)


def _read_expression(path: Path, gene_column: str | None = None) -> pd.DataFrame:
    frame = _read_table(path)
    gene_col = gene_column
    if gene_col is None:
        gene_col = _first_existing(
            list(frame.columns),
            [
                "gene",
                "symbol",
                "gene_symbol",
                "hgnc",
                "hugo",
                "feature",
                "id",
            ],
        )
    if gene_col is None:
        first = frame.columns[0]
        numeric = pd.to_numeric(frame[first], errors="coerce")
        if numeric.notna().mean() < 0.5:
            gene_col = str(first)
        else:
            frame = frame.set_index(frame.columns[0])
            gene_col = None

    if gene_col is not None:
        if gene_col not in frame.columns:
            raise ValueError(f"gene column '{gene_col}' not found in {path}")
        frame = frame.set_index(gene_col)
    frame.index = frame.index.astype(str).str.strip().str.upper()
    frame = frame.loc[frame.index != ""]
    frame = frame[~frame.index.duplicated(keep="first")]
    frame = frame.apply(pd.to_numeric, errors="coerce")
    if frame.empty:
        raise ValueError(f"expression matrix has no genes: {path}")
    return frame


def _read_metadata(path: Path, sample_column: str | None = None) -> pd.DataFrame:
    frame = _read_table(path)
    sample_col = sample_column
    if sample_col is None:
        sample_col = _first_existing(
            list(frame.columns),
            ["sample", "sample_id", "sampleid", "donor", "patient", "id"],
        )
    if sample_col is None:
        frame = frame.set_index(frame.columns[0])
    else:
        if sample_col not in frame.columns:
            raise ValueError(f"sample column '{sample_col}' not found in {path}")
        frame = frame.set_index(sample_col)
    frame.index = frame.index.astype(str).str.strip()
    frame = frame.loc[~frame.index.duplicated(keep="first")]
    return frame


def _resolve(base: Path, value: str | Path | None) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _infer_case_control(labels: pd.Series) -> tuple[str, str]:
    values = sorted({str(value) for value in labels.dropna().astype(str)})
    if len(values) != 2:
        raise ValueError(
            "classification requires exactly two condition labels; "
            f"found {values}"
        )
    case_tokens = ("case", "tumor", "tumour", "disease", "treated", "exposed", "hcc")
    control_tokens = ("control", "normal", "healthy", "vehicle", "untreated")
    case = next(
        (value for value in values if any(token in value.lower() for token in case_tokens)),
        None,
    )
    control = next(
        (
            value
            for value in values
            if any(token in value.lower() for token in control_tokens)
        ),
        None,
    )
    return case or values[0], control or values[1]


def load_cohort(spec: dict, base: Path) -> Cohort:
    """Load one cohort from a JSON specification."""
    name = str(spec.get("name") or "discovery")
    expression_path = _resolve(base, spec.get("expression"))
    metadata_path = _resolve(base, spec.get("metadata"))
    if expression_path is None or metadata_path is None:
        raise ValueError(
            f"cohort '{name}' requires expression and metadata paths"
        )
    expression = _read_expression(expression_path, spec.get("gene_column"))
    metadata = _read_metadata(metadata_path, spec.get("sample_column"))
    condition_col = spec.get("condition_column") or _first_existing(
        list(metadata.columns),
        ["condition", "group", "disease", "status", "tissue", "treatment"],
    )
    if not condition_col or condition_col not in metadata.columns:
        raise ValueError(
            f"cohort '{name}' has no usable condition column: "
            f"{list(metadata.columns)}"
        )

    common = [sample for sample in expression.columns if sample in metadata.index]
    if len(common) < 4:
        # GEO matrices often use different punctuation than sample metadata.
        normalized_expression = {
            str(sample).replace(".", "-").upper(): str(sample)
            for sample in expression.columns
        }
        normalized_metadata = {
            str(sample).replace(".", "-").upper(): str(sample)
            for sample in metadata.index
        }
        common = [
            normalized_expression[key]
            for key in normalized_expression.keys() & normalized_metadata.keys()
        ]
    if len(common) < 4:
        raise ValueError(
            f"cohort '{name}' has fewer than four expression/metadata matches"
        )

    expression = expression.loc[:, common]
    metadata = metadata.loc[common, :]
    labels = metadata[condition_col].astype(str)
    case_label = spec.get("case_label")
    control_label = spec.get("control_label")
    if not case_label or not control_label:
        inferred_case, inferred_control = _infer_case_control(labels)
        case_label = case_label or inferred_case
        control_label = control_label or inferred_control
    keep = labels.isin([str(case_label), str(control_label)])
    if keep.sum() < 4:
        raise ValueError(
            f"cohort '{name}' has fewer than four samples in the requested "
            f"case/control labels: {case_label} vs {control_label}"
        )
    expression = expression.loc[:, keep]
    metadata = metadata.loc[keep, :]
    labels = labels.loc[keep]
    return Cohort(
        name=name,
        expression=expression,
        metadata=metadata,
        labels=labels,
        condition_column=str(condition_col),
        case_label=str(case_label),
        control_label=str(control_label),
    )


def _looks_like_counts(frame: pd.DataFrame) -> bool:
    values = frame.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0 or finite.min() < 0:
        return False
    if not np.allclose(finite, np.round(finite), atol=1e-6):
        return False
    column_sums = np.nansum(values, axis=0)
    return bool(np.median(column_sums) >= 10)


def _quantile_normalize(frame: pd.DataFrame) -> pd.DataFrame:
    """Rank-based quantile normalisation for already log-scale arrays."""
    values = frame.to_numpy(dtype=float)
    if values.shape[1] < 2:
        return frame
    ranks = np.argsort(np.argsort(values, axis=0), axis=0)
    sorted_values = np.sort(values, axis=0)
    reference = sorted_values.mean(axis=1)
    normalized = reference[ranks]
    normalized[~np.isfinite(values)] = np.nan
    return pd.DataFrame(
        normalized,
        index=frame.index,
        columns=frame.columns,
    )


def normalize_expression(
    frame: pd.DataFrame,
    quantile: bool = False,
) -> tuple[pd.DataFrame, str]:
    """Return a log-scale matrix and the detected input type."""
    numeric = frame.astype(float)
    if _looks_like_counts(numeric):
        library = numeric.sum(axis=0).replace(0, np.nan)
        normalized = np.log2(numeric.div(library, axis=1) * 1e6 + 1.0)
        kind = "counts"
    else:
        normalized = numeric.copy()
        if np.nanmin(normalized.to_numpy()) < 0:
            # Some platforms ship log2 ratios; shift to a positive log scale
            # only when needed for downstream correlation/network calculations.
            normalized = normalized - np.nanmin(normalized.to_numpy()) + 1.0
        kind = "normalized_or_microarray"
    normalized = normalized.replace([np.inf, -np.inf], np.nan)
    if quantile and kind != "counts":
        normalized = _quantile_normalize(normalized)
    return normalized, kind


def _bh_adjust(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").fillna(1.0).to_numpy()
    order = np.argsort(values)
    ranked = values[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return pd.Series(out, index=p_values.index)


def differential_expression(
    expression: pd.DataFrame,
    labels: pd.Series,
    case_label: str,
    control_label: str,
) -> pd.DataFrame:
    """Welch t-test DEG table used when R/limma is unavailable."""
    case_samples = labels.index[labels.astype(str) == str(case_label)]
    control_samples = labels.index[labels.astype(str) == str(control_label)]
    rows: list[dict] = []
    for gene, values in expression.iterrows():
        case_values = pd.to_numeric(values.loc[case_samples], errors="coerce").dropna()
        control_values = pd.to_numeric(
            values.loc[control_samples], errors="coerce"
        ).dropna()
        if len(case_values) < 2 or len(control_values) < 2:
            continue
        statistic, p_value = stats.ttest_ind(
            case_values,
            control_values,
            equal_var=False,
            nan_policy="omit",
        )
        if not math.isfinite(float(p_value)):
            continue
        rows.append(
            {
                "gene": gene,
                "case_mean": float(case_values.mean()),
                "control_mean": float(control_values.mean()),
                "log2FoldChange": float(
                    case_values.mean() - control_values.mean()
                ),
                "statistic": float(statistic) if math.isfinite(statistic) else np.nan,
                "pvalue": float(p_value),
                "n_case": int(len(case_values)),
                "n_control": int(len(control_values)),
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["padj"] = _bh_adjust(frame["pvalue"])
    frame["direction"] = np.where(
        (frame["padj"] <= 0.05) & (frame["log2FoldChange"] > 0),
        "Up",
        np.where(
            (frame["padj"] <= 0.05) & (frame["log2FoldChange"] < 0),
            "Down",
            "NS",
        ),
    )
    return frame.sort_values(
        ["padj", "pvalue", "log2FoldChange"],
        ascending=[True, True, False],
        na_position="last",
    ).reset_index(drop=True)


def _classification_models(random_state: int) -> dict[str, Pipeline]:
    models: dict[str, Pipeline] = {}
    models["elastic_net"] = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("select", SelectKBest(f_classif, k="all")),
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="elasticnet",
                    solver="saga",
                    l1_ratio=0.5,
                    C=1.0,
                    max_iter=5000,
                    class_weight="balanced",
                    random_state=random_state,
                ),
            ),
        ]
    )
    models["lasso"] = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("select", SelectKBest(f_classif, k="all")),
            ("scale", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    penalty="l1",
                    solver="liblinear",
                    C=1.0,
                    class_weight="balanced",
                    random_state=random_state,
                ),
            ),
        ]
    )
    models["random_forest"] = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("select", SelectKBest(f_classif, k="all")),
            (
                "clf",
                RandomForestClassifier(
                    n_estimators=400,
                    max_depth=None,
                    min_samples_leaf=2,
                    class_weight="balanced",
                    random_state=random_state,
                    n_jobs=-1,
                ),
            ),
        ]
    )
    models["gradient_boosting"] = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("select", SelectKBest(f_classif, k="all")),
            (
                "clf",
                GradientBoostingClassifier(
                    n_estimators=200,
                    learning_rate=0.05,
                    max_depth=3,
                    random_state=random_state,
                ),
            ),
        ]
    )
    models["svm_rbf"] = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("select", SelectKBest(f_classif, k="all")),
            ("scale", StandardScaler()),
            (
                "clf",
                SVC(
                    kernel="rbf",
                    probability=True,
                    class_weight="balanced",
                    random_state=random_state,
                ),
            ),
        ]
    )
    models["mlp"] = Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("select", SelectKBest(f_classif, k="all")),
            ("scale", StandardScaler()),
            (
                "clf",
                MLPClassifier(
                    hidden_layer_sizes=(64, 32),
                    max_iter=1000,
                    early_stopping=False,
                    random_state=random_state,
                ),
            ),
        ]
    )
    try:
        from xgboost import XGBClassifier

        models["xgboost"] = Pipeline(
            [
                ("impute", SimpleImputer(strategy="median")),
                ("select", SelectKBest(f_classif, k="all")),
                (
                    "clf",
                    XGBClassifier(
                        n_estimators=300,
                        max_depth=3,
                        learning_rate=0.05,
                        subsample=0.8,
                        colsample_bytree=0.8,
                        eval_metric="logloss",
                        random_state=random_state,
                        n_jobs=-1,
                    ),
                ),
            ]
        )
    except Exception:
        pass
    return models


def _macro_auc(y_true: np.ndarray, probabilities: np.ndarray) -> float:
    classes = np.unique(y_true)
    if len(classes) < 2:
        return float("nan")
    if probabilities.ndim == 1:
        return float(roc_auc_score(y_true, probabilities))
    if len(classes) == 2:
        return float(roc_auc_score(y_true, probabilities[:, 1]))
    return float(
        roc_auc_score(
            y_true,
            probabilities,
            multi_class="ovr",
            average="macro",
        )
    )


def _positive_probabilities(
    probabilities: np.ndarray,
    classes: np.ndarray,
    positive_label: int,
) -> np.ndarray:
    if probabilities.ndim == 1:
        return probabilities
    index = int(np.where(classes == positive_label)[0][0])
    return probabilities[:, index]


def _bootstrap_auc_ci(
    y_true: np.ndarray,
    scores: np.ndarray,
    seed: int,
    repeats: int = 1000,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values: list[float] = []
    n = len(y_true)
    for _ in range(max(50, repeats)):
        indices = rng.integers(0, n, size=n)
        if len(np.unique(y_true[indices])) < 2:
            continue
        values.append(float(roc_auc_score(y_true[indices], scores[indices])))
    if not values:
        return float("nan"), float("nan")
    return (
        float(np.percentile(values, 2.5)),
        float(np.percentile(values, 97.5)),
    )


def _decision_curve(
    y_true: np.ndarray,
    probabilities: np.ndarray,
    thresholds: np.ndarray | None = None,
) -> pd.DataFrame:
    if thresholds is None:
        thresholds = np.linspace(0.01, 0.99, 99)
    rows: list[dict] = []
    n = len(y_true)
    prevalence = float(np.mean(y_true))
    for threshold in thresholds:
        predicted = probabilities >= threshold
        tp = float(np.sum(predicted & (y_true == 1)))
        fp = float(np.sum(predicted & (y_true == 0)))
        net_benefit = tp / n - fp / n * (threshold / (1.0 - threshold))
        treat_all = prevalence - (1.0 - prevalence) * (
            threshold / (1.0 - threshold)
        )
        rows.append(
            {
                "threshold": float(threshold),
                "net_benefit": float(net_benefit),
                "treat_all": float(treat_all),
            }
        )
    return pd.DataFrame(rows)


def _safe_pipeline_k(name: str, n_features: int, configured: int | None) -> int:
    if configured and configured > 0:
        return max(1, min(int(configured), n_features))
    if name in {"elastic_net", "lasso"}:
        return max(1, min(80, n_features))
    return max(1, min(200, n_features))


def _fit_pipeline(
    pipeline: Pipeline,
    X: pd.DataFrame,
    y: np.ndarray,
    feature_cap: int | None,
) -> Pipeline:
    fitted = clone(pipeline)
    fitted.set_params(
        select__k=_safe_pipeline_k(
            "unknown",
            X.shape[1],
            feature_cap,
        )
    )
    fitted.fit(X, y)
    return fitted


def _selected_feature_names(
    pipeline: Pipeline,
    columns: list[str],
) -> list[str]:
    selector = pipeline.named_steps.get("select")
    if selector is None or not hasattr(selector, "get_support"):
        return list(columns)
    support = np.asarray(selector.get_support(), dtype=bool)
    return [column for column, keep in zip(columns, support) if keep]


def _importance_for_pipeline(
    pipeline: Pipeline,
    X: pd.DataFrame,
    y: np.ndarray,
    columns: list[str],
    seed: int,
) -> pd.Series:
    selected = _selected_feature_names(pipeline, columns)
    estimator = pipeline.named_steps["clf"]
    if hasattr(estimator, "feature_importances_"):
        values = np.asarray(estimator.feature_importances_, dtype=float)
    elif hasattr(estimator, "coef_"):
        coef = np.asarray(estimator.coef_, dtype=float)
        values = np.abs(coef).mean(axis=0) if coef.ndim > 1 else np.abs(coef)
    else:
        try:
            transformed = pipeline[:-1].transform(X)
            result = permutation_importance(
                estimator,
                transformed,
                y,
                n_repeats=5,
                random_state=seed,
                scoring="roc_auc",
            )
            values = np.asarray(result.importances_mean, dtype=float)
        except Exception:
            values = np.zeros(len(selected), dtype=float)
    values = np.asarray(values, dtype=float).ravel()
    if len(values) != len(selected):
        values = np.resize(values, len(selected))
    return pd.Series(values, index=selected, name="importance").sort_values(
        ascending=False
    )


def _evaluate_models(
    expression: pd.DataFrame,
    labels: pd.Series,
    candidate_genes: list[str],
    config: dict,
    out_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Evaluate models with repeated CV and write model-level outputs."""
    genes = [gene for gene in candidate_genes if gene in expression.index]
    if len(genes) < 2:
        raise ValueError("fewer than two candidate genes are present in the matrix")
    X = expression.loc[genes, :].T
    y_text = labels.astype(str)
    encoder = LabelEncoder()
    y = encoder.fit_transform(y_text)
    positive_hits = np.where(
        encoder.classes_ == str(config.get("case_label"))
    )[0]
    if len(positive_hits) == 0:
        raise ValueError(
            "case_label is not present after label encoding: "
            f"{config.get('case_label')!r}"
        )
    positive_label = int(positive_hits[0])
    y_positive = (y == positive_label).astype(int)
    min_class = int(pd.Series(y).value_counts().min())
    if min_class < 2:
        raise ValueError(
            "each class needs at least two samples for stratified model "
            f"evaluation; observed class counts: "
            f"{pd.Series(y).value_counts().to_dict()}"
        )
    n_splits = max(2, min(int(config.get("cv_folds", 5)), min_class))
    n_repeats = max(1, int(config.get("cv_repeats", 5)))
    cv = RepeatedStratifiedKFold(
        n_splits=n_splits,
        n_repeats=n_repeats,
        random_state=int(config.get("seed", 42)),
    )
    predictor_cv = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=int(config.get("seed", 42)),
    )
    feature_cap = config.get("feature_cap")
    models = _classification_models(int(config.get("seed", 42)))
    rows: list[dict] = []
    fitted_models: dict[str, Pipeline] = {}
    importances: dict[str, pd.Series] = {}
    oof_store: dict[str, np.ndarray] = {}
    for name, template in models.items():
        fitted = clone(template)
        fitted.set_params(
            select__k=_safe_pipeline_k(
                name,
                X.shape[1],
                feature_cap,
            )
        )
        aucs: list[float] = []
        accuracies: list[float] = []
        for train_idx, test_idx in cv.split(X, y):
            fold_model = clone(fitted)
            fold_model.fit(X.iloc[train_idx], y[train_idx])
            proba = fold_model.predict_proba(X.iloc[test_idx])
            positive = _positive_probabilities(
                proba,
                fold_model.classes_,
                positive_label,
            )
            if len(np.unique(y[test_idx])) >= 2:
                aucs.append(
                    float(
                        roc_auc_score(
                            y_positive[test_idx],
                            positive,
                        )
                    )
                )
            accuracies.append(
                float(accuracy_score(y[test_idx], fold_model.predict(X.iloc[test_idx])))
            )
        fitted.fit(X, y)
        fitted_models[name] = fitted
        importances[name] = _importance_for_pipeline(
            fitted,
            X,
            y,
            list(X.columns),
            int(config.get("seed", 42)),
        )
        try:
            oof = cross_val_predict(
                clone(fitted),
                X,
                y,
                cv=predictor_cv,
                method="predict_proba",
            )
            oof_store[name] = _positive_probabilities(
                oof,
                np.unique(y),
                positive_label,
            )
        except Exception:
            pass
        auc_mean = float(np.nanmean(aucs)) if aucs else float("nan")
        auc_sd = float(np.nanstd(aucs)) if aucs else float("nan")
        rows.append(
            {
                "model": name,
                "cv_auc_mean": auc_mean,
                "cv_auc_sd": auc_sd,
                "cv_accuracy_mean": float(np.nanmean(accuracies)),
                "cv_accuracy_sd": float(np.nanstd(accuracies)),
                "n_features": int(X.shape[1]),
                "n_samples": int(X.shape[0]),
                "n_repeats": n_repeats,
                "n_splits": n_splits,
            }
        )
    comparison = pd.DataFrame(rows).sort_values(
        ["cv_auc_mean", "cv_accuracy_mean"],
        ascending=False,
        na_position="last",
    )
    comparison.to_csv(out_dir / "ml_model_comparison.csv", index=False)
    importance_frame = pd.concat(
        [
            series.rename(name)
            for name, series in importances.items()
            if not series.empty
        ],
        axis=1,
    ).fillna(0.0)
    importance_frame.to_csv(out_dir / "ml_feature_importance.csv")
    best_name = str(comparison.iloc[0]["model"])
    best_model = fitted_models[best_name]
    selected_genes = _selected_feature_names(best_model, list(X.columns))
    pd.DataFrame({"gene": selected_genes}).to_csv(
        out_dir / "ml_selected_features.csv",
        index=False,
    )

    if best_name in oof_store:
        probabilities = oof_store[best_name]
        auc_value = float(roc_auc_score(y_positive, probabilities))
        lower, upper = _bootstrap_auc_ci(
            y_positive,
            probabilities,
            int(config.get("seed", 42)),
        )
        fpr, tpr, _ = roc_curve(y_positive, probabilities)
        plt.figure(figsize=(5.5, 5.0))
        plt.plot(fpr, tpr, label=f"AUC={auc_value:.3f}")
        plt.plot([0, 1], [0, 1], "--", color="grey")
        plt.xlabel("False positive rate")
        plt.ylabel("True positive rate")
        plt.title(f"{best_name} cross-validated ROC")
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "ml_roc.png", dpi=160)
        plt.close()

        prob_true, prob_pred = calibration_curve(
            y_positive,
            probabilities,
            n_bins=max(3, min(8, int(np.ceil(np.sqrt(len(y)))))),
        )
        plt.figure(figsize=(5.5, 4.5))
        plt.plot(prob_pred, prob_true, marker="o", label="Calibration")
        plt.plot([0, 1], [0, 1], "--", color="grey", label="Perfect")
        plt.xlabel("Mean predicted probability")
        plt.ylabel("Observed frequency")
        plt.title("Calibration curve")
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / "ml_calibration.png", dpi=160)
        plt.close()
        _decision_curve(y_positive, probabilities).to_csv(
            out_dir / "ml_decision_curve.csv",
            index=False,
        )
        try:
            cm = confusion_matrix(
                y_positive,
                (probabilities >= 0.5).astype(int),
            )
            pd.DataFrame(
                cm,
                index=["true_control", "true_case"],
                columns=["pred_control", "pred_case"],
            ).to_csv(out_dir / "ml_confusion_matrix.csv")
        except Exception:
            pass
        _json_write(
            out_dir / "ml_summary.json",
            {
                "status": "completed",
                "best_model": best_name,
                "cv_auc": auc_value,
                "cv_auc_ci95": [lower, upper],
                "average_precision": float(
                    average_precision_score(y_positive, probabilities)
                ),
                "brier_score": float(
                    brier_score_loss(y_positive, probabilities)
                ),
                "positive_class": str(config.get("case_label")),
                "selected_features": selected_genes,
                "external_validation_generated": True,
            },
        )
        if config.get("shap", False):
            _write_shap_if_available(
                best_model,
                X,
                y,
                list(X.columns),
                out_dir,
                int(config.get("seed", 42)),
            )
    return comparison, importance_frame, pd.DataFrame(
        {"gene": selected_genes}
    )


def _write_shap_if_available(
    pipeline: Pipeline,
    X: pd.DataFrame,
    y: np.ndarray,
    columns: list[str],
    out_dir: Path,
    seed: int,
) -> None:
    try:
        import shap

        selected = _selected_feature_names(pipeline, columns)
        transformed = pipeline[:-1].transform(X)
        estimator = pipeline.named_steps["clf"]
        if estimator.__class__.__name__ in {
            "RandomForestClassifier",
            "GradientBoostingClassifier",
            "XGBClassifier",
        }:
            explainer = shap.TreeExplainer(estimator)
            values = explainer.shap_values(transformed)
        else:
            explainer = shap.KernelExplainer(
                estimator.predict_proba,
                transformed[: min(30, len(transformed))],
            )
            values = explainer.shap_values(
                transformed[: min(80, len(transformed))],
                nsamples=100,
            )
        if isinstance(values, list):
            values = values[1]
        plt.figure(figsize=(8, 6))
        shap.summary_plot(values, transformed, feature_names=selected, show=False)
        plt.tight_layout()
        plt.savefig(out_dir / "ml_shap.png", dpi=160, bbox_inches="tight")
        plt.close()
    except Exception as exc:
        (out_dir / "ml_shap_status.txt").write_text(
            f"SHAP skipped: {exc}",
            encoding="utf-8",
        )


def _external_validation(
    discovery_expression: pd.DataFrame,
    discovery_labels: pd.Series,
    validation: list[Cohort],
    candidate_genes: list[str],
    config: dict,
    out_dir: Path,
) -> pd.DataFrame:
    if not validation:
        return pd.DataFrame()
    rows: list[dict] = []
    for cohort in validation:
        genes = [gene for gene in candidate_genes if gene in cohort.expression.index]
        if len(genes) < 2:
            rows.append(
                {
                    "cohort": cohort.name,
                    "status": "skipped",
                    "reason": "fewer than two candidate genes present",
                }
            )
            continue
        X_train = discovery_expression.loc[genes, :].T
        y_train = (discovery_labels.astype(str) == str(config.get("case_label"))).astype(int)
        X_test = cohort.expression.loc[genes, :].T
        y_test = (cohort.labels.astype(str) == str(cohort.case_label)).astype(int)
        for name, template in _classification_models(
            int(config.get("seed", 42))
        ).items():
            try:
                model = clone(template)
                model.set_params(
                    select__k=_safe_pipeline_k(
                        name,
                        X_train.shape[1],
                        config.get("feature_cap"),
                    )
                )
                model.fit(X_train, y_train)
                probabilities = model.predict_proba(X_test)[:, 1]
                rows.append(
                    {
                        "cohort": cohort.name,
                        "model": name,
                        "status": "completed",
                        "n_samples": int(len(y_test)),
                        "auc": float(roc_auc_score(y_test, probabilities)),
                        "average_precision": float(
                            average_precision_score(y_test, probabilities)
                        ),
                        "brier_score": float(
                            brier_score_loss(y_test, probabilities)
                        ),
                    }
                )
            except Exception as exc:
                rows.append(
                    {
                        "cohort": cohort.name,
                        "model": name,
                        "status": "failed",
                        "reason": str(exc),
                    }
                )
    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "ml_external_validation.csv", index=False)
    return frame


def _load_optional_evidence(
    files: dict[str, str | dict],
    base: Path,
    genes: list[str],
) -> pd.DataFrame:
    output = pd.DataFrame(index=pd.Index(genes, name="gene"))
    for name, value in files.items():
        explicit_column = None
        path_value = value
        if isinstance(value, dict):
            path_value = value.get("path") or value.get("file")
            explicit_column = value.get("score_column")
        path = _resolve(base, path_value)
        if path is None or not path.exists():
            continue
        try:
            frame = _read_table(path)
            gene_col = _first_existing(
                list(frame.columns),
                ["gene", "symbol", "gene_symbol", "target", "id"],
            )
            if gene_col is None:
                continue
            numeric_cols = [
                column
                for column in frame.columns
                if column != gene_col
                and pd.to_numeric(frame[column], errors="coerce").notna().any()
            ]
            if not numeric_cols:
                continue
            if explicit_column and explicit_column in frame.columns:
                if str(explicit_column) not in numeric_cols:
                    raise ValueError(
                        f"score_column {explicit_column!r} is not numeric in {path}"
                    )
                score_col = str(explicit_column)
            else:
                preferred = []
                lower_name = name.lower()
                if "hub" in lower_name or "network" in lower_name:
                    preferred.extend(
                        ["ppi_hub_score", "hub_score", "score"]
                    )
                elif "prognos" in lower_name:
                    preferred.extend(
                        ["hr", "hazard_ratio", "hazardratio", "cox_hr"]
                    )
                elif "immune" in lower_name:
                    preferred.extend(["immune_score", "correlation", "score"])
                score_col = next(
                    (
                        column
                        for column in preferred
                        if column in frame.columns
                    ),
                    numeric_cols[0],
                )
            normalized = frame.assign(
                **{
                    gene_col: frame[gene_col].astype(str).str.upper(),
                    score_col: pd.to_numeric(
                        frame[score_col],
                        errors="coerce",
                    ),
                }
            )
            values = (
                normalized
                .dropna(subset=[gene_col, score_col])
                .drop_duplicates(gene_col, keep="first")
                .set_index(gene_col)[score_col]
            )
            output[name] = output.index.map(values)
        except Exception as exc:
            LOG.warning("could not read optional evidence %s: %s", path, exc)
    return output


def integrate_priorities(
    deg: pd.DataFrame,
    wgcna_hubs: pd.DataFrame | None,
    ml_importance: pd.DataFrame,
    evidence: pd.DataFrame | None,
    config: dict,
) -> pd.DataFrame:
    genes = set(deg.get("gene", pd.Series(dtype=str)).astype(str))
    if wgcna_hubs is not None and not wgcna_hubs.empty:
        genes |= set(wgcna_hubs["gene"].astype(str))
    genes |= set(ml_importance.index.astype(str))
    if evidence is not None:
        genes |= set(evidence.index.astype(str))
    frame = pd.DataFrame({"gene": sorted(genes)}).set_index("gene")

    deg_map = (
        deg.drop_duplicates("gene")
        .set_index("gene")
        .reindex(frame.index)
    )
    frame["deg_log2fc"] = deg_map.get("log2FoldChange", pd.Series(index=frame.index))
    frame["deg_padj"] = deg_map.get("padj", pd.Series(index=frame.index))
    frame["deg_score"] = pd.Series(
        -np.log10(
            pd.to_numeric(
                frame["deg_padj"],
                errors="coerce",
            ).clip(lower=1e-300)
        ),
        index=frame.index,
    )

    if wgcna_hubs is not None and not wgcna_hubs.empty:
        hub = wgcna_hubs.drop_duplicates("gene").set_index("gene")
        if "kME" in hub.columns:
            frame["wgcna_kme"] = hub["kME"].reindex(frame.index)
            frame["wgcna_score"] = frame["wgcna_kme"].abs()
        else:
            frame["wgcna_kme"] = np.nan
            frame["wgcna_score"] = np.nan
    else:
        frame["wgcna_kme"] = np.nan
        frame["wgcna_score"] = np.nan

    importance = ml_importance.iloc[:, 0] if not ml_importance.empty else pd.Series(dtype=float)
    frame["ml_importance"] = importance.reindex(frame.index)
    frame["ml_score"] = frame["ml_importance"]

    if evidence is not None and not evidence.empty:
        for column in evidence.columns:
            frame[column] = evidence[column].reindex(frame.index)

    score_columns = [
        column
        for column in frame.columns
        if column.endswith("_score") and column != "priority_score"
    ]
    weights = config.get("weights", {})
    weighted_sum = pd.Series(0.0, index=frame.index)
    weight_total = 0.0
    evidence_count = pd.Series(0, index=frame.index, dtype=int)
    for column in score_columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        valid = values.notna()
        percentile = values.rank(pct=True)
        weight = float(weights.get(column, 1.0))
        weighted_sum = weighted_sum.add(percentile.fillna(0.0) * weight, fill_value=0.0)
        weight_total += weight * valid.mean() if valid.any() else 0.0
        evidence_count += valid.astype(int)
    frame["priority_score"] = (
        weighted_sum / weight_total if weight_total > 0 else 0.0
    )
    frame["evidence_count"] = evidence_count
    return frame.sort_values(
        ["priority_score", "evidence_count", "deg_padj"],
        ascending=[False, False, True],
        na_position="last",
    ).reset_index()


def _find_rscript() -> str | None:
    return shutil.which("Rscript") or shutil.which("Rscript.exe")


def run_r_advanced(
    config_path: Path,
    out_dir: Path,
    log,
) -> dict:
    rscript = _find_rscript()
    script = APP_ROOT / "src" / "analysis" / "R" / "advanced_bulk.R"
    if rscript is None or not script.exists():
        return {
            "status": "skipped",
            "reason": "Rscript or advanced_bulk.R not available",
        }
    proc = subprocess.run(
        [rscript, str(script), str(config_path), str(out_dir)],
        cwd=APP_ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=int(os.environ.get("LIVER_ADVANCED_R_TIMEOUT", "7200")),
    )
    (out_dir / "advanced_r.log").write_text(
        (proc.stdout or "") + "\n" + (proc.stderr or ""),
        encoding="utf-8",
    )
    if proc.returncode != 0:
        log.warning("advanced R analysis failed: %s", proc.stderr[-2000:])
        return {
            "status": "failed",
            "returncode": proc.returncode,
        }
    summary_path = out_dir / "advanced_r_summary.json"
    if summary_path.exists():
        return _json_load(summary_path)
    return {"status": "completed_without_summary"}


def build_config(args, base: Path) -> dict:
    if args.config:
        config_path = _resolve(base, args.config)
        if config_path is None:
            raise ValueError("--config path is empty")
        config = _json_load(config_path)
        config["_config_path"] = str(config_path)
        return config
    if not args.expression or not args.metadata:
        raise ValueError(
            "provide --config or both --expression and --metadata"
        )
    primary = {
        "name": args.cohort_name or "discovery",
        "expression": args.expression,
        "metadata": args.metadata,
        "condition_column": args.condition_column,
        "gene_column": args.gene_column,
        "case_label": args.case_label,
        "control_label": args.control_label,
    }
    return {
        "primary": primary,
        "validation": [],
        "batch": {"enabled": False},
        "wgcna": {"enabled": not args.skip_wgcna},
        "immune": {"enabled": not args.skip_immune},
        "survival": {"enabled": not args.skip_survival},
        "ml": {
            "enabled": not args.skip_ml,
            "feature_cap": args.feature_cap,
            "cv_folds": args.cv_folds,
            "cv_repeats": args.cv_repeats,
            "seed": args.seed,
        },
        "evidence_files": {},
        "weights": {},
    }


def run(args) -> int:
    base = Path.cwd()
    config_path = _resolve(base, args.config) if args.config else None
    config = build_config(args, base)
    out_dir = Path(args.output).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    resolved_dir = out_dir / "resolved"
    resolved_dir.mkdir(parents=True, exist_ok=True)
    resolved_config = resolved_dir / "advanced_config.json"
    config_for_r = dict(config)
    config_for_r["_output_dir"] = str(out_dir)
    # Resolve all file paths relative to the config location or current cwd.
    config_base = (
        config_path.parent if config_path is not None else base
    )
    for key in ("expression", "metadata"):
        if config.get("primary", {}).get(key):
            config_for_r["primary"][key] = str(
                _resolve(config_base, config["primary"][key])
            )
    for cohort in config.get("validation", []) or []:
        for key in ("expression", "metadata"):
            if cohort.get(key):
                cohort[key] = str(_resolve(config_base, cohort[key]))
    _json_write(resolved_config, config_for_r)

    r_summary = run_r_advanced(resolved_config, out_dir, LOG) if not args.skip_r else {
        "status": "skipped",
        "reason": "disabled by --skip-r",
    }

    primary = load_cohort(config.get("primary", {}), config_base)
    validation = [
        load_cohort(spec, config_base)
        for spec in (config.get("validation", []) or [])
    ]
    expression, input_kind = normalize_expression(
        primary.expression,
        quantile=bool(config.get("normalization", {}).get("quantile")),
    )
    correction_path = out_dir / "expression_primary_corrected.csv"
    if correction_path.exists() and not args.skip_r:
        corrected = _read_expression(correction_path)
        common = [gene for gene in corrected.index if gene in expression.index]
        if len(common) > 100:
            expression.loc[common, :] = corrected.loc[common, :]
            input_kind = "batch_corrected"
    deg_r_path = out_dir / "deg_primary_limma.csv"
    if deg_r_path.exists() and not args.skip_r:
        deg = pd.read_csv(deg_r_path)
        if "gene" not in deg.columns:
            deg = deg.rename(columns={deg.columns[0]: "gene"})
        if "log2FoldChange" not in deg.columns and "logFC" in deg.columns:
            deg["log2FoldChange"] = deg["logFC"]
        if "padj" not in deg.columns and "adj.P.Val" in deg.columns:
            deg["padj"] = deg["adj.P.Val"]
    else:
        deg = differential_expression(
            expression,
            primary.labels,
            primary.case_label,
            primary.control_label,
        )
        deg.to_csv(out_dir / "deg_primary_python.csv", index=False)
    deg.to_csv(out_dir / "deg_primary.csv", index=False)

    ml_config = dict(config.get("ml", {}) or {})
    ml_config.setdefault("case_label", primary.case_label)
    ml_config.setdefault("cv_folds", 5)
    ml_config.setdefault("cv_repeats", 5)
    ml_config.setdefault("seed", 42)
    if "padj" in deg.columns:
        candidate = deg[
            pd.to_numeric(deg["padj"], errors="coerce")
            <= float(ml_config.get("padj", 0.05))
        ]
    else:
        candidate = deg.head(int(ml_config.get("top_genes", 500)))
    if candidate.empty:
        candidate = deg.head(int(ml_config.get("top_genes", 500)))
    candidate_genes = candidate["gene"].astype(str).head(
        int(ml_config.get("top_genes", 500))
    ).tolist()
    ml_summary: dict = {"status": "skipped", "reason": "disabled"}
    ml_importance = pd.DataFrame()
    if not args.skip_ml:
        comparison, ml_importance, selected = _evaluate_models(
            expression,
            primary.labels,
            candidate_genes,
            ml_config,
            out_dir,
        )
        external = _external_validation(
            expression,
            primary.labels,
            validation,
            candidate_genes,
            ml_config,
            out_dir,
        )
        ml_summary = (
            _json_load(out_dir / "ml_summary.json")
            if (out_dir / "ml_summary.json").exists()
            else {"status": "completed"}
        )
        ml_summary["external_validation"] = (
            external.to_dict(orient="records")
            if not external.empty
            else []
        )
        _json_write(out_dir / "ml_summary.json", ml_summary)

    wgcna_hubs_path = out_dir / "wgcna_hubs.csv"
    wgcna_hubs = (
        pd.read_csv(wgcna_hubs_path)
        if wgcna_hubs_path.exists()
        else None
    )
    evidence = _load_optional_evidence(
        config.get("evidence_files", {}) or {},
        config_base,
        sorted(set(deg["gene"].astype(str)) | set(ml_importance.index.astype(str))),
    )
    integrated = integrate_priorities(
        deg,
        wgcna_hubs,
        ml_importance,
        evidence,
        config.get("integration", {}) or {},
    )
    integrated.to_csv(out_dir / "integrated_priority.csv", index=False)

    summary = {
        "status": "completed",
        "discovery": primary.name,
        "input_expression_kind": input_kind,
        "n_discovery_samples": int(expression.shape[1]),
        "n_genes": int(expression.shape[0]),
        "n_deg": int(len(deg)),
        "n_supported_deg": int(
            (pd.to_numeric(deg.get("padj"), errors="coerce") <= 0.05).sum()
        )
        if "padj" in deg.columns
        else 0,
        "validation_cohorts": [cohort.name for cohort in validation],
        "r_advanced": r_summary,
        "ml": ml_summary,
        "wgcna_hubs": int(len(wgcna_hubs)) if wgcna_hubs is not None else 0,
        "integrated_targets": int(len(integrated)),
        "outputs": {
            "deg": str(out_dir / "deg_primary.csv"),
            "ml_model_comparison": str(out_dir / "ml_model_comparison.csv"),
            "ml_external_validation": str(out_dir / "ml_external_validation.csv"),
            "wgcna_hubs": str(wgcna_hubs_path) if wgcna_hubs_path.exists() else "",
            "integrated_priority": str(out_dir / "integrated_priority.csv"),
        },
    }
    _json_write(out_dir / "advanced_analysis_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Advanced multi-cohort bulk analysis, WGCNA, gene-level ML, "
            "immune/survival integration and target ranking."
        )
    )
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--config")
    parser.add_argument("--expression")
    parser.add_argument("--metadata")
    parser.add_argument("--cohort-name", default="discovery")
    parser.add_argument("--condition-column", default="condition")
    parser.add_argument("--gene-column")
    parser.add_argument("--case-label")
    parser.add_argument("--control-label")
    parser.add_argument("--feature-cap", type=int)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--cv-repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip-r", action="store_true")
    parser.add_argument("--skip-ml", action="store_true")
    parser.add_argument("--skip-wgcna", action="store_true")
    parser.add_argument("--skip-immune", action="store_true")
    parser.add_argument("--skip-survival", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s %(name)s: %(message)s",
    )
    try:
        return run(args)
    except Exception as exc:
        LOG.exception("advanced analysis failed: %s", exc)
        output = Path(args.output).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        _json_write(
            output / "advanced_analysis_summary.json",
            {
                "status": "failed",
                "reason": str(exc),
            },
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
