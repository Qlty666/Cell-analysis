"""Leakage-aware model comparison, external validation and SHAP analysis."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV, calibration_curve
from sklearn.ensemble import (
    ExtraTreesClassifier,
    GradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from .common import ensure_dir, read_expression, save_figure, write_json

LOG = logging.getLogger("experiment_plan_one.ml")
warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")


def _model_zoo(seed: int = 42) -> dict[str, Any]:
    models: dict[str, Any] = {
        "LASSO": LogisticRegression(
            penalty="l1",
            C=0.2,
            solver="liblinear",
            class_weight="balanced",
            max_iter=5000,
            random_state=seed,
        ),
        "Ridge": LogisticRegression(
            penalty="l2",
            C=0.5,
            solver="liblinear",
            class_weight="balanced",
            max_iter=5000,
            random_state=seed,
        ),
        "ElasticNet": LogisticRegression(
            penalty="elasticnet",
            C=0.3,
            l1_ratio=0.5,
            solver="saga",
            class_weight="balanced",
            max_iter=10000,
            random_state=seed,
        ),
        "SVM-linear": SVC(
            kernel="linear",
            C=1.0,
            probability=True,
            class_weight="balanced",
            random_state=seed,
        ),
        "SVM-RBF": SVC(
            kernel="rbf",
            C=2.0,
            gamma="scale",
            probability=True,
            class_weight="balanced",
            random_state=seed,
        ),
        "RandomForest": RandomForestClassifier(
            n_estimators=500,
            max_features="sqrt",
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=seed,
            n_jobs=-1,
        ),
        "ExtraTrees": ExtraTreesClassifier(
            n_estimators=500,
            max_features="sqrt",
            min_samples_leaf=2,
            class_weight="balanced",
            random_state=seed,
            n_jobs=-1,
        ),
        "GradientBoosting": GradientBoostingClassifier(
            n_estimators=150,
            learning_rate=0.04,
            max_depth=2,
            random_state=seed,
        ),
        "MLP": MLPClassifier(
            hidden_layer_sizes=(32, 16),
            alpha=0.001,
            learning_rate_init=0.002,
            max_iter=1500,
            early_stopping=True,
            random_state=seed,
        ),
    }
    for name, factory in {
        "XGBoost": lambda: _xgboost(seed),
        "LightGBM": lambda: _lightgbm(seed),
    }.items():
        try:
            models[name] = factory()
        except RuntimeError as exc:
            LOG.warning("optional model %s unavailable: %s", name, exc)
    return models


def _xgboost(seed: int) -> Any:
    try:
        from xgboost import XGBClassifier
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"xgboost unavailable: {exc}") from exc
    return XGBClassifier(
        n_estimators=350,
        max_depth=3,
        learning_rate=0.03,
        subsample=0.85,
        colsample_bytree=0.8,
        reg_alpha=0.05,
        reg_lambda=1.0,
        eval_metric="logloss",
        random_state=seed,
        n_jobs=2,
    )


def _lightgbm(seed: int) -> Any:
    try:
        from lightgbm import LGBMClassifier
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"lightgbm unavailable: {exc}") from exc
    return LGBMClassifier(
        n_estimators=350,
        num_leaves=15,
        learning_rate=0.03,
        subsample=0.85,
        colsample_bytree=0.8,
        reg_alpha=0.05,
        reg_lambda=0.2,
        random_state=seed,
        n_jobs=2,
        verbose=-1,
    )


def _pipeline(
    estimator: Any,
    *,
    selector: SelectKBest | None = None,
) -> Pipeline:
    steps: list[tuple[str, Any]] = [
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
    ]
    if selector is not None:
        steps.append(("selector", selector))
    steps.append(("model", estimator))
    return Pipeline(steps)


def _fit_calibrated(
    estimator: Pipeline,
    X: pd.DataFrame,
    y: pd.Series,
    *,
    seed: int,
) -> Any:
    """Fit a probability calibrator using training data only."""
    min_class = int(y.value_counts().min())
    folds = min(3, min_class)
    if folds < 2:
        return clone(estimator).fit(X, y)
    calibrated = CalibratedClassifierCV(
        estimator=clone(estimator),
        method="sigmoid",
        cv=folds,
    )
    calibrated.fit(X, y)
    return calibrated


def _nested_cv_evaluation(
    candidate_models: dict[tuple[str, str], Pipeline],
    X_by_feature: dict[str, pd.DataFrame],
    y: pd.Series,
    *,
    seed: int,
    outer_folds: int,
    inner_folds: int,
) -> tuple[pd.DataFrame, dict[str, Any], np.ndarray]:
    """Select feature set and algorithm inside each outer training split."""
    outer = StratifiedKFold(
        n_splits=outer_folds,
        shuffle=True,
        random_state=seed,
    )
    probabilities = np.full(len(y), np.nan, dtype=float)
    selected_rows: list[dict[str, Any]] = []
    for fold_number, (train_index, test_index) in enumerate(
        outer.split(next(iter(X_by_feature.values())), y),
        start=1,
    ):
        candidate_scores: dict[tuple[str, str], float] = {}
        for key, estimator in candidate_models.items():
            X = X_by_feature[key[0]].iloc[train_index]
            y_train = y.iloc[train_index]
            min_class = int(y_train.value_counts().min())
            folds = max(2, min(inner_folds, min_class))
            inner = StratifiedKFold(
                n_splits=folds,
                shuffle=True,
                random_state=seed + fold_number,
            )
            scores: list[float] = []
            for inner_train, inner_test in inner.split(X, y_train):
                try:
                    fitted = clone(estimator).fit(
                        X.iloc[inner_train],
                        y_train.iloc[inner_train],
                    )
                    probability = fitted.predict_proba(
                        X.iloc[inner_test]
                    )[:, 1]
                    scores.append(
                        float(
                            roc_auc_score(
                                y_train.iloc[inner_test],
                                probability,
                            )
                        )
                    )
                except Exception:  # noqa: BLE001
                    continue
            if scores:
                candidate_scores[key] = float(np.mean(scores))
        if not candidate_scores:
            raise RuntimeError(
                f"no nested-CV candidate completed in outer fold {fold_number}"
            )
        best_key = max(
            candidate_scores,
            key=lambda key: (candidate_scores[key], key[1], key[0]),
        )
        X_train = X_by_feature[best_key[0]].iloc[train_index]
        X_test = X_by_feature[best_key[0]].iloc[test_index]
        y_train = y.iloc[train_index]
        model = _fit_calibrated(
            candidate_models[best_key],
            X_train,
            y_train,
            seed=seed,
        )
        probabilities[test_index] = model.predict_proba(X_test)[:, 1]
        selected_rows.append(
            {
                "outer_fold": fold_number,
                "feature_set": best_key[0],
                "model": best_key[1],
                "inner_auc": candidate_scores[best_key],
                "n_train": int(len(train_index)),
                "n_test": int(len(test_index)),
            }
        )
    if np.isnan(probabilities).any():
        raise RuntimeError("nested cross-validation left samples without predictions")
    selected = pd.DataFrame(selected_rows)
    metrics = {
        "auc": float(roc_auc_score(y, probabilities)),
        "average_precision": float(
            average_precision_score(y, probabilities)
        ),
        "brier": float(brier_score_loss(y, probabilities)),
        "selected_models": selected["model"].value_counts().to_dict(),
        "selected_feature_sets": selected["feature_set"].value_counts().to_dict(),
    }
    return selected, metrics, probabilities


def run_ml_validation(
    training_expression_path: Path,
    training_metadata_path: Path,
    candidate_sets: dict[str, list[str]],
    validation_datasets: dict[str, dict[str, Any]],
    output_dir: Path,
    *,
    seed: int = 42,
    cv_folds: int = 5,
    cv_repeats: int = 5,
) -> dict[str, Any]:
    """Compare models and validate the best model without outcome leakage."""
    ensure_dir(output_dir)
    expression = read_expression(training_expression_path)
    metadata = pd.read_csv(training_metadata_path, index_col=0)
    sample_ids = [sample for sample in expression.columns if sample in metadata.index]
    condition = metadata.loc[sample_ids, "condition"].astype(str)
    mapping = {
        "HC": 0,
        "SS": 1,
        "NASH": 1,
        "control": 0,
        "normal": 0,
        "healthy_control": 0,
        "NAFLD": 1,
        "case": 1,
        "disease": 1,
    }
    y = condition.map(mapping)
    valid = y.notna()
    sample_ids = y[valid].index.tolist()
    y = y.loc[sample_ids].astype(int)
    if y.nunique() != 2:
        raise ValueError("training labels must contain exactly two disease classes")

    # Rank harmonization removes platform-specific mean/variance shifts
    # without using outcome labels. Each gene is converted to its within-
    # cohort percentile before model fitting or external prediction.
    all_matrix = (
        expression[sample_ids]
        .T.astype(float)
        .rank(axis=0, method="average", pct=True)
        .fillna(0.5)
    )
    feature_sets: dict[str, tuple[list[str], int | None, str]] = {}
    for name, genes in candidate_sets.items():
        usable = sorted(set(gene.upper() for gene in genes) & set(expression.index))
        if len(usable) >= 2:
            feature_sets[name] = (usable, None, "static gene set")
    feature_sets["DEG_top50"] = (list(expression.index), 50, "fold-specific univariate selection")
    if not feature_sets:
        raise ValueError("none of the candidate gene sets are present in GSE89632")

    rows: list[dict[str, Any]] = []
    fitted_models: dict[tuple[str, str], tuple[Pipeline, pd.Series]] = {}
    cv_predictions: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    candidate_models: dict[tuple[str, str], Pipeline] = {}
    X_by_feature = {
        name: all_matrix[features]
        for name, (features, _, _) in feature_sets.items()
    }
    model_zoo = _model_zoo(seed)
    splitter = RepeatedStratifiedKFold(
        n_splits=cv_folds,
        n_repeats=max(1, int(cv_repeats)),
        random_state=seed,
    )

    for feature_name, (features, dynamic_k, feature_description) in feature_sets.items():
        X = all_matrix[features]
        for model_name, estimator in model_zoo.items():
            selector = (
                SelectKBest(score_func=f_classif, k=min(dynamic_k, len(features)))
                if dynamic_k
                else None
            )
            model = _pipeline(estimator, selector=selector)
            candidate_models[(feature_name, model_name)] = model
            fold_auc: list[float] = []
            fold_ap: list[float] = []
            fold_brier: list[float] = []
            prediction_sum = np.zeros(len(y), dtype=float)
            prediction_count = np.zeros(len(y), dtype=int)
            try:
                for train_index, test_index in splitter.split(X, y):
                    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
                    y_train = y.iloc[train_index]
                    fitted = clone(model).fit(X_train, y_train)
                    probability = fitted.predict_proba(X_test)[:, 1]
                    prediction_sum[test_index] += probability
                    prediction_count[test_index] += 1
                    fold_auc.append(float(roc_auc_score(y.iloc[test_index], probability)))
                    fold_ap.append(
                        float(average_precision_score(y.iloc[test_index], probability))
                    )
                    fold_brier.append(
                        float(brier_score_loss(y.iloc[test_index], probability))
                    )
                fitted = clone(model).fit(X, y)
                prediction = np.divide(
                    prediction_sum,
                    np.maximum(prediction_count, 1),
                )
            except Exception as exc:  # noqa: BLE001
                LOG.warning(
                    "model failed: %s x %s: %s",
                    model_name,
                    feature_name,
                    exc,
                )
                rows.append(
                    {
                            "feature_set": feature_name,
                            "feature_description": feature_description,
                            "model": model_name,
                            "n_features": (
                                min(int(dynamic_k), len(features))
                                if dynamic_k
                                else len(features)
                            ),
                            "n_input_features": len(features),
                            "status": "failed",
                        "error": str(exc),
                    }
                )
                continue
            rows.append(
                {
                    "feature_set": feature_name,
                    "feature_description": feature_description,
                    "model": model_name,
                    "n_features": (
                        min(int(dynamic_k), len(features))
                        if dynamic_k
                        else len(features)
                    ),
                    "n_input_features": len(features),
                    "status": "ok",
                    "cv_auc_mean": float(np.mean(fold_auc)),
                    "cv_auc_sd": float(np.std(fold_auc, ddof=1)),
                    "cv_auc_folds": ";".join(f"{value:.6f}" for value in fold_auc),
                    "cv_auc_ci_low": float(np.percentile(fold_auc, 2.5)),
                    "cv_auc_ci_high": float(np.percentile(fold_auc, 97.5)),
                    "cv_repeats": int(max(1, cv_repeats)),
                    "cv_average_precision": float(np.mean(fold_ap)),
                    "cv_brier": float(np.mean(fold_brier)),
                }
            )
            fitted_models[(feature_name, model_name)] = (fitted, y)
            cv_predictions[(feature_name, model_name)] = {
                "y": y.to_numpy(),
                "probability": prediction,
            }

    performance = pd.DataFrame(rows)
    performance.to_csv(output_dir / "model_performance.csv", index=False)
    valid_performance = performance[performance["status"] == "ok"].sort_values(
        ["cv_auc_mean", "cv_average_precision"], ascending=False
    )
    if valid_performance.empty:
        raise RuntimeError("all ML models failed")
    nested_selection, nested_metrics, nested_probability = _nested_cv_evaluation(
        candidate_models,
        X_by_feature,
        y,
        seed=seed,
        outer_folds=cv_folds,
        inner_folds=max(2, min(3, cv_folds)),
    )
    nested_selection.to_csv(
        output_dir / "nested_model_selection.csv",
        index=False,
    )
    selected_model = max(
        nested_metrics["selected_models"],
        key=lambda value: (
            nested_metrics["selected_models"][value],
            value,
        ),
    )
    selected_feature_set = max(
        nested_metrics["selected_feature_sets"],
        key=lambda value: (
            nested_metrics["selected_feature_sets"][value],
            value,
        ),
    )
    best_key = (selected_feature_set, selected_model)
    best_model_raw, _ = fitted_models[best_key]
    best_model = _fit_calibrated(
        best_model_raw,
        X_by_feature[best_key[0]],
        y,
        seed=seed,
    )
    cv = cv_predictions[best_key]
    best = valid_performance[
        (valid_performance["feature_set"] == best_key[0])
        & (valid_performance["model"] == best_key[1])
    ].iloc[0]

    _plot_auc_heatmap(performance, output_dir / "fig3a_multimodel_auc_heatmap.png")
    _plot_roc(
        cv["y"],
        cv["probability"],
        output_dir / "fig3b_best_model_roc_training.png",
        title=f"{best_key[0]} / {best_key[1]}",
    )
    calibration = _plot_calibration(
        y.to_numpy(),
        nested_probability,
        output_dir / "fig3e_calibration_curve.png",
    )
    decision_curve = _plot_decision_curve(
        y.to_numpy(),
        nested_probability,
        output_dir / "fig3e_decision_curve.png",
    )
    calibration["decision_curve"] = decision_curve
    calibration = {**calibration, **nested_metrics}
    cv_frame = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "observed": cv["y"],
            "probability": cv["probability"],
        }
    )
    cv_frame.to_csv(output_dir / "best_model_cv_predictions.csv", index=False)

    validation_rows: list[dict[str, Any]] = []
    for dataset_name, spec in validation_datasets.items():
        result = validate_external_model(
            best_model,
            training_expression_path=training_expression_path,
            validation_expression_path=Path(spec["expression"]),
            validation_metadata_path=Path(spec["metadata"]),
            condition_map=spec.get("condition_map") or {},
            comparison=spec.get("comparison") or ("case", "control"),
            output_dir=output_dir,
            dataset_name=dataset_name,
            title=spec.get("title") or dataset_name,
            endpoint_warning=str(spec.get("endpoint_warning") or ""),
            evaluate_auc_target=bool(spec.get("evaluate_auc_target", True)),
            endpoint_role=str(spec.get("endpoint_role") or ""),
        )
        validation_rows.append(result)
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(output_dir / "external_validation_metrics.csv", index=False)

    shap = _shap_analysis(
        best_model_raw,
        all_matrix,
        y,
        output_dir,
        feature_key=best_key,
    )
    core_genes = shap["top_genes"][:3]
    write_json(
        output_dir / "ml_summary.json",
        {
            "best_feature_set": best_key[0],
            "best_model": best_key[1],
            "cv_auc": float(best["cv_auc_mean"]),
            "cv_auc_sd": float(best["cv_auc_sd"]),
            "cv_auc_ci_low": float(best["cv_auc_ci_low"]),
            "cv_auc_ci_high": float(best["cv_auc_ci_high"]),
            "cv_folds": int(cv_folds),
            "cv_repeats": int(max(1, cv_repeats)),
            "nested_cv": nested_metrics,
            "calibration_method": "sigmoid_calibrated_outer_fold",
            "cv_auc_target": 0.85,
            "cv_auc_target_met": bool(float(best["cv_auc_mean"]) >= 0.85),
            "calibration": calibration,
            "external_validation": validation.to_dict(orient="records"),
            "core_genes": core_genes,
            "shap_status": shap["status"],
        },
    )
    return {
        "performance": output_dir / "model_performance.csv",
        "external_validation": output_dir / "external_validation_metrics.csv",
        "best_model": best_key,
        "cv_auc": float(best["cv_auc_mean"]),
        "core_genes": core_genes,
        "shap_values": shap.get("values_path"),
    }


def validate_external_model(
    fitted_model: Pipeline,
    *,
    training_expression_path: Path,
    validation_expression_path: Path,
    validation_metadata_path: Path,
    condition_map: dict[str, int],
    comparison: tuple[str, str],
    output_dir: Path,
    dataset_name: str,
    title: str,
    endpoint_warning: str = "",
    evaluate_auc_target: bool = True,
    endpoint_role: str = "",
) -> dict[str, Any]:
    training = read_expression(training_expression_path)
    validation = read_expression(validation_expression_path)
    ranked_training = training.T.astype(float).rank(axis=0, method="average", pct=True)
    ranked_validation = validation.T.astype(float).rank(axis=0, method="average", pct=True)
    metadata = pd.read_csv(validation_metadata_path, index_col=0)
    model_columns = _estimator_feature_names(fitted_model)
    if not model_columns:
        imputer = fitted_model.named_steps.get("imputer")
        model_columns = list(getattr(imputer, "feature_names_in_", []))
    if not model_columns:
        raise ValueError(f"could not determine fitted feature names for {dataset_name}")
    matrix = ranked_validation.reindex(columns=model_columns)
    for column in model_columns:
        if column not in ranked_validation.columns:
            matrix[column] = (
                float(ranked_training[column].median())
                if column in ranked_training.columns
                else 0.5
            )
        matrix[column] = pd.to_numeric(matrix[column], errors="coerce").fillna(
            float(ranked_training[column].median())
            if column in ranked_training.columns
            else 0.5
        )
    condition = metadata.loc[metadata.index.intersection(matrix.index), "condition"].astype(str)
    mapped = condition.map(condition_map)
    valid = mapped.notna()
    if valid.sum() < 4:
        raise ValueError(f"too few mapped validation samples for {dataset_name}")
    X = matrix.loc[valid[valid].index]
    y = mapped[valid].astype(int)
    probabilities = fitted_model.predict_proba(X)[:, 1]
    auc = float(roc_auc_score(y, probabilities))
    rng = np.random.default_rng(42)
    auc_boot: list[float] = []
    ap_boot: list[float] = []
    for _ in range(1000):
        indices = rng.integers(0, len(y), size=len(y))
        sampled = y.to_numpy()[indices]
        if len(set(sampled)) < 2:
            continue
        sampled_probability = probabilities[indices]
        auc_boot.append(
            float(roc_auc_score(sampled, sampled_probability))
        )
        ap_boot.append(
            float(
                average_precision_score(
                    sampled,
                    sampled_probability,
                )
            )
        )
    calibration = _calibration_parameters(
        y.to_numpy(),
        probabilities,
    )
    calibration_curve_frame = _external_calibration_curve(
        y.to_numpy(),
        probabilities,
        output_dir / f"{dataset_name}_calibration.png",
        title=title,
    )
    fpr, tpr, thresholds = roc_curve(y, probabilities)
    prediction = pd.DataFrame(
        {
            "sample_id": X.index,
            "condition": condition.loc[X.index].to_numpy(),
            "label": y.to_numpy(),
            "probability": probabilities,
        }
    )
    prediction.to_csv(output_dir / f"{dataset_name}_predictions.csv", index=False)
    pd.DataFrame(
        {
            "fpr": fpr,
            "tpr": tpr,
            "threshold": thresholds,
        }
    ).to_csv(output_dir / f"{dataset_name}_roc_curve.csv", index=False)
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    ax.plot(fpr, tpr, color="#b24c3c", lw=2, label=f"AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#8a949e", lw=1)
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_title(title, fontweight="bold")
    ax.legend(frameon=False)
    save_figure(fig, output_dir / f"{dataset_name}_roc.png")
    return {
        "dataset": dataset_name,
        "comparison": f"{comparison[0]} vs {comparison[1]}",
        "n": int(len(y)),
        "n_positive": int(y.sum()),
        "auc": auc,
        "average_precision": float(average_precision_score(y, probabilities)),
        "auroc_ci_low": (
            float(np.percentile(auc_boot, 2.5)) if auc_boot else None
        ),
        "auroc_ci_high": (
            float(np.percentile(auc_boot, 97.5)) if auc_boot else None
        ),
        "auprc_ci_low": (
            float(np.percentile(ap_boot, 2.5)) if ap_boot else None
        ),
        "auprc_ci_high": (
            float(np.percentile(ap_boot, 97.5)) if ap_boot else None
        ),
        "brier": float(brier_score_loss(y, probabilities)),
        "calibration_slope": calibration["slope"],
        "calibration_intercept": calibration["intercept"],
        "calibration_bins": calibration_curve_frame.to_dict(orient="records"),
        "target_auc": 0.8 if evaluate_auc_target else None,
        "target_met": (
            bool(auc >= 0.8)
            if evaluate_auc_target
            else None
        ),
        "evaluate_auc_target": bool(evaluate_auc_target),
        "endpoint_role": endpoint_role,
        "endpoint_warning": endpoint_warning,
    }


def _estimator_feature_names(estimator: Any) -> list[str]:
    """Return fitted feature names from plain or calibrated estimators."""
    names = list(getattr(estimator, "feature_names_in_", []))
    if names:
        return names
    if hasattr(estimator, "named_steps"):
        for step in estimator.named_steps.values():
            names = _estimator_feature_names(step)
            if names:
                return names
    calibrated = getattr(estimator, "calibrated_classifiers_", None)
    if calibrated:
        for wrapper in calibrated:
            names = _estimator_feature_names(getattr(wrapper, "estimator", None))
            if names:
                return names
    return []


def _plot_auc_heatmap(performance: pd.DataFrame, output: Path) -> None:
    frame = performance[performance["status"] == "ok"]
    if frame.empty:
        return
    matrix = frame.pivot(index="feature_set", columns="model", values="cv_auc_mean")
    fig, ax = plt.subplots(figsize=(7.2, max(3.8, matrix.shape[0] * 0.7)))
    image = ax.imshow(matrix.to_numpy(), cmap="YlGnBu", vmin=0.5, vmax=1.0)
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=5.8)
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=6)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix.iloc[row, column]
            if pd.notna(value):
                ax.text(
                    column,
                    row,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=6.5,
                    color="white" if value >= 0.82 else "#1f2933",
                )
    ax.set_title("Repeated cross-validated AUC", fontweight="bold")
    fig.colorbar(image, ax=ax, label="AUC", shrink=0.7)
    save_figure(fig, output)


def _plot_roc(y: np.ndarray, probability: np.ndarray, output: Path, title: str) -> None:
    fpr, tpr, _ = roc_curve(y, probability)
    auc = roc_auc_score(y, probability)
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    ax.plot(fpr, tpr, color="#b24c3c", lw=2.2, label=f"AUC={auc:.3f}")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#8a949e", lw=1)
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_title(title, fontweight="bold")
    ax.legend(frameon=False)
    save_figure(fig, output)


def _plot_calibration(
    y: np.ndarray,
    probability: np.ndarray,
    output: Path,
) -> dict[str, Any]:
    observed, predicted = calibration_curve(
        y,
        probability,
        n_bins=min(8, max(3, len(y) // 7)),
        strategy="quantile",
    )
    hl_p = _hosmer_lemeshow(y, probability)
    calibration = _calibration_parameters(y, probability)
    fig, ax = plt.subplots(figsize=(4.8, 4.6))
    ax.plot(predicted, observed, marker="o", color="#2f6bb3", lw=1.8, label="Observed")
    ax.plot([0, 1], [0, 1], linestyle="--", color="#8a949e", lw=1)
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Observed proportion")
    ax.set_title("Calibration", fontweight="bold")
    ax.text(
        0.03,
        0.97,
        (
            f"Hosmer-Lemeshow p={hl_p:.3g}\n"
            f"Brier={brier_score_loss(y, probability):.3f}\n"
            f"slope={calibration['slope']:.2f}; intercept={calibration['intercept']:.2f}"
        ),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    ax.legend(frameon=False)
    save_figure(fig, output)
    return {
        "hosmer_lemeshow_p": hl_p,
        "brier": float(brier_score_loss(y, probability)),
        "calibration_slope": calibration["slope"],
        "calibration_intercept": calibration["intercept"],
        "target_met": bool(hl_p is not None and np.isfinite(hl_p) and hl_p > 0.05),
    }


def _plot_decision_curve(
    y: np.ndarray,
    probability: np.ndarray,
    output: Path,
) -> dict[str, Any]:
    thresholds = np.linspace(0.01, 0.99, 99)
    prevalence = float(np.mean(y))
    rows: list[dict[str, float]] = []
    for threshold in thresholds:
        selected = probability >= threshold
        true_positive = float(np.sum(selected & (y == 1)))
        false_positive = float(np.sum(selected & (y == 0)))
        odds = threshold / (1.0 - threshold)
        model_net_benefit = true_positive / len(y) - false_positive / len(y) * odds
        all_net_benefit = prevalence - (1.0 - prevalence) * odds
        rows.append(
            {
                "threshold": float(threshold),
                "model": float(model_net_benefit),
                "treat_all": float(all_net_benefit),
                "treat_none": 0.0,
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(output.with_suffix(".csv"), index=False)
    fig, ax = plt.subplots(figsize=(4.9, 4.2))
    ax.plot(frame["threshold"], frame["model"], color="#2f6bb3", lw=2, label="Model")
    ax.plot(
        frame["threshold"],
        frame["treat_all"],
        color="#8a949e",
        lw=1.5,
        linestyle="--",
        label="Treat all",
    )
    ax.axhline(0, color="#c7d0d8", lw=1)
    ax.set_xlabel("Threshold probability")
    ax.set_ylabel("Net benefit")
    ax.set_title("Decision curve analysis", fontweight="bold")
    ax.legend(frameon=False)
    save_figure(fig, output)
    positive = frame.loc[frame["model"] > 0, "threshold"]
    return {
        "max_model_net_benefit": float(frame["model"].max()),
        "max_treat_all_net_benefit": float(frame["treat_all"].max()),
        "positive_net_benefit_threshold_low": (
            float(positive.min()) if not positive.empty else None
        ),
        "positive_net_benefit_threshold_high": (
            float(positive.max()) if not positive.empty else None
        ),
    }


def _external_calibration_curve(
    y: np.ndarray,
    probability: np.ndarray,
    output: Path,
    *,
    title: str,
) -> pd.DataFrame:
    observed, predicted = calibration_curve(
        y,
        probability,
        n_bins=min(6, max(2, len(y) // 6)),
        strategy="quantile",
    )
    frame = pd.DataFrame(
        {
            "predicted_probability": predicted,
            "observed_probability": observed,
        }
    )
    frame.to_csv(output.with_suffix(".csv"), index=False)
    fig, ax = plt.subplots(figsize=(4.8, 4.4))
    ax.plot(predicted, observed, marker="o", color="#2f6bb3", lw=1.8)
    ax.plot([0, 1], [0, 1], linestyle="--", color="#8a949e", lw=1)
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Observed proportion")
    ax.set_title(title, fontweight="bold")
    ax.text(
        0.03,
        0.97,
        f"Brier={brier_score_loss(y, probability):.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    save_figure(fig, output)
    return frame


def _calibration_parameters(
    y: np.ndarray,
    probability: np.ndarray,
) -> dict[str, float]:
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    logits = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    try:
        model = LogisticRegression(
            penalty=None,
            solver="lbfgs",
            max_iter=2000,
        )
        model.fit(logits, y)
        return {
            "slope": float(model.coef_[0, 0]),
            "intercept": float(model.intercept_[0]),
        }
    except Exception:  # noqa: BLE001
        return {"slope": float("nan"), "intercept": float("nan")}


def _hosmer_lemeshow(y: np.ndarray, probability: np.ndarray, bins: int = 8) -> float:
    if len(y) < 4:
        return float("nan")
    bins = min(bins, max(2, len(y) // 3))
    try:
        groups = pd.qcut(probability, q=bins, duplicates="drop")
    except ValueError:
        return float("nan")
    statistic = 0.0
    degrees = 0
    for _, indices in pd.Series(np.arange(len(y))).groupby(groups, observed=False):
        idx = indices.to_numpy()
        if len(idx) == 0:
            continue
        observed = float(y[idx].sum())
        expected = float(probability[idx].sum())
        n = len(idx)
        if expected <= 0 or expected >= n:
            continue
        statistic += (observed - expected) ** 2 / (expected * (1 - expected / n))
        degrees += 1
    degrees = max(1, degrees - 2)
    return float(1.0 - stats.chi2.cdf(statistic, degrees))


def _shap_analysis(
    model: Pipeline,
    X: pd.DataFrame,
    y: pd.Series,
    output_dir: Path,
    *,
    feature_key: tuple[str, str],
) -> dict[str, Any]:
    try:
        import shap

        transformed, feature_names = _transform_for_explanation(model, X)
        estimator = model.named_steps["model"]
        if estimator.__class__.__name__ in {
            "RandomForestClassifier",
            "ExtraTreesClassifier",
            "GradientBoostingClassifier",
            "XGBClassifier",
            "LGBMClassifier",
        }:
            explainer = shap.TreeExplainer(estimator)
            values = explainer.shap_values(transformed)
        elif estimator.__class__.__name__ in {
            "LogisticRegression",
            "LinearSVC",
        }:
            explainer = shap.LinearExplainer(estimator, transformed)
            values = explainer.shap_values(transformed)
        else:
            background = shap.sample(transformed, min(30, transformed.shape[0]))
            explainer = shap.KernelExplainer(
                lambda matrix: estimator.predict_proba(matrix)[:, 1],
                background,
            )
            values = explainer.shap_values(
                transformed[: min(80, transformed.shape[0])],
                nsamples=100,
            )
        if isinstance(values, list):
            values = values[1] if len(values) > 1 else values[0]
        values = np.asarray(values)
        if values.ndim == 3:
            values = values[:, :, -1]
        if values.shape[0] != transformed.shape[0]:
            transformed = transformed[: values.shape[0]]
        mean_abs = np.mean(np.abs(values), axis=0)
        importance = (
            pd.DataFrame({"gene": feature_names, "mean_abs_shap": mean_abs})
            .sort_values("mean_abs_shap", ascending=False)
            .reset_index(drop=True)
        )
        importance.to_csv(output_dir / "fig3f_shap_importance.csv", index=False)
        values_path = output_dir / "shap_values.csv"
        pd.DataFrame(values, columns=feature_names).to_csv(values_path, index=False)
        _shap_bar(importance, output_dir / "fig3f_shap_importance.png")
        _shap_beeswarm(
            values,
            transformed,
            feature_names,
            output_dir / "fig3g_shap_beeswarm.png",
        )
        return {
            "status": "completed",
            "top_genes": importance["gene"].tolist(),
            "values_path": values_path,
            "importance": output_dir / "fig3f_shap_importance.csv",
        }
    except Exception as exc:  # noqa: BLE001
        LOG.warning("SHAP failed; using model importance fallback: %s", exc)
        importance = _fallback_importance(model, X)
        importance.to_csv(output_dir / "fig3f_shap_importance.csv", index=False)
        _shap_bar(importance, output_dir / "fig3f_shap_importance.png")
        return {
            "status": f"fallback: {exc}",
            "top_genes": importance["gene"].tolist(),
            "values_path": None,
            "importance": output_dir / "fig3f_shap_importance.csv",
        }


def _transform_for_explanation(model: Pipeline, X: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    transformed = X.to_numpy(dtype=float)
    feature_names = list(X.columns)
    for name, step in model.steps[:-1]:
        transformed = step.transform(transformed)
        if name == "selector" and hasattr(step, "get_support"):
            feature_names = list(np.asarray(feature_names)[step.get_support()])
    return np.asarray(transformed, dtype=float), feature_names


def _fallback_importance(model: Pipeline, X: pd.DataFrame) -> pd.DataFrame:
    estimator = model.named_steps["model"]
    _, feature_names = _transform_for_explanation(model, X)
    if hasattr(estimator, "feature_importances_"):
        values = np.asarray(estimator.feature_importances_, dtype=float)
    elif hasattr(estimator, "coef_"):
        values = np.abs(np.asarray(estimator.coef_, dtype=float)).ravel()
    else:
        values = np.ones(len(feature_names), dtype=float)
    return pd.DataFrame({"gene": feature_names, "mean_abs_shap": values}).sort_values(
        "mean_abs_shap", ascending=False
    )


def _shap_bar(importance: pd.DataFrame, output: Path) -> None:
    top = importance.head(20).sort_values("mean_abs_shap", ascending=True)
    fig, ax = plt.subplots(figsize=(6.4, max(4, len(top) * 0.32)))
    ax.barh(top["gene"], top["mean_abs_shap"], color="#b24c3c")
    ax.set_xlabel("Mean |SHAP value|")
    ax.set_ylabel("")
    ax.set_title("Global feature importance", fontweight="bold")
    save_figure(fig, output)


def _shap_beeswarm(
    values: np.ndarray,
    transformed: np.ndarray,
    feature_names: list[str],
    output: Path,
) -> None:
    try:
        import shap

        shap.summary_plot(
            values,
            transformed,
            feature_names=feature_names,
            max_display=20,
            show=False,
            plot_size=(7.2, 6.2),
        )
        fig = plt.gcf()
        save_figure(fig, output, tight=False)
    except Exception as exc:  # noqa: BLE001
        LOG.warning("SHAP beeswarm rendering failed: %s", exc)
