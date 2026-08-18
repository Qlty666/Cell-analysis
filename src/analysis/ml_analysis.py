#!/usr/bin/env python3
"""ML disease classification, feature importance, and SHAP explainability."""

import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.ensemble import (
    GradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.feature_selection import RFE, SelectFromModel
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    auc,
    average_precision_score,
    classification_report,
    confusion_matrix,
    precision_recall_curve,
    roc_curve,
)
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_predict,
    cross_val_score,
)
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.svm import SVC


def model_name() -> str:
    name = os.environ.get("LIVER_ML_MODEL", "xgb").strip().lower()
    if name not in {"xgb", "rf", "gbm", "mlp", "lasso_svm"}:
        name = "xgb"
    return name


def build_model(name: str, random_state: int = 42):
    if name == "rf":
        return RandomForestClassifier(
            n_estimators=300,
            max_depth=8,
            random_state=random_state,
            n_jobs=-1,
        )
    if name == "gbm":
        return GradientBoostingClassifier(
            n_estimators=200,
            max_depth=3,
            learning_rate=0.05,
            random_state=random_state,
        )
    if name == "mlp":
        return MLPClassifier(
            hidden_layer_sizes=(64, 32),
            max_iter=500,
            early_stopping=True,
            random_state=random_state,
        )
    if name == "lasso_svm":
        return SVC(kernel="linear", probability=True, random_state=random_state)
    try:
        from xgboost import XGBClassifier

        return XGBClassifier(
            n_estimators=200,
            max_depth=4,
            use_label_encoder=False,
            eval_metric="logloss",
            random_state=random_state,
            verbosity=0,
        )
    except Exception:
        return RandomForestClassifier(
            n_estimators=300,
            max_depth=8,
            random_state=random_state,
            n_jobs=-1,
        )


def _importance_values(model, X, y, columns):
    if hasattr(model, "feature_importances_"):
        values = np.asarray(model.feature_importances_, dtype=float)
    elif hasattr(model, "coef_"):
        coef = np.asarray(model.coef_, dtype=float)
        values = np.abs(coef).mean(axis=0) if coef.ndim > 1 else np.abs(coef)
    else:
        try:
            from sklearn.inspection import permutation_importance

            result = permutation_importance(
                model,
                X,
                y,
                n_repeats=3,
                random_state=42,
                scoring="accuracy",
            )
            values = result.importances_mean
        except Exception:
            values = np.zeros(len(columns), dtype=float)
    series = pd.Series(
        values,
        index=columns,
        name="importance",
    ).sort_values(ascending=False)
    return series


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    data_dir = root / "results" / "data" / "07_ml"
    fig_dir = root / "results" / "figures" / "07_ml"
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    status = {
        "status": "skipped",
        "reason": "",
        "model": "",
        "accuracy": None,
        "cv_scores": [],
    }

    try:
        sc_data_dir = root / "results" / "data"
        ann_path = sc_data_dir / "04_annotation" / "fig_05_16_17_cell_annotations.csv"
        qc_path = sc_data_dir / "01_qc" / "fig_01_qc_metrics.csv"
        if not ann_path.exists() or not qc_path.exists():
            status["reason"] = (
                "missing fig_05_16_17_cell_annotations.csv "
                "or fig_01_qc_metrics.csv"
            )
            (data_dir / "ml_model_summary.json").write_text(
                json.dumps(status, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            return 0

        ann = pd.read_csv(ann_path)
        qc = pd.read_csv(qc_path)

        props = pd.crosstab(ann["sample"], ann["celltype_annot"], normalize="index")
        meta = ann[["sample", "condition"]].drop_duplicates()
        qc_means = qc.groupby("sample")[
            ["nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"]
        ].mean()
        features = props.join(meta.set_index("sample"), how="inner")
        features = features.join(qc_means, how="left")
        features = features.dropna(subset=["condition"])

        y = features["condition"]
        le = LabelEncoder()
        y_enc = le.fit_transform(y)
        X = features.drop(columns=["condition"]).fillna(0)
        if len(set(y_enc)) < 2 or X.shape[0] < 4:
            status["reason"] = "insufficient samples or groups for classification"
            (data_dir / "ml_model_summary.json").write_text(
                json.dumps(status, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            return 0

        chosen_model = model_name()
        fit_X = X
        selected_columns = list(X.columns)
        feature_selection = "none"
        if chosen_model == "lasso_svm":
            lasso_selector = SelectFromModel(
                LogisticRegression(
                    penalty="l1",
                    solver="liblinear",
                    C=1.0,
                    max_iter=1000,
                    random_state=42,
                )
            )
            lasso_selector.fit(X, y_enc)
            fit_X = lasso_selector.transform(X)
            selected_columns = X.columns[lasso_selector.get_support()].tolist()
            if fit_X.shape[1] < 1:
                selected_columns = list(X.columns[:5])
                fit_X = X[selected_columns].to_numpy()
            if fit_X.shape[1] >= 2:
                n_features = max(2, min(15, fit_X.shape[1]))
                rfe = RFE(
                    SVC(kernel="linear", probability=True, random_state=42),
                    n_features_to_select=n_features,
                )
                rfe.fit(fit_X, y_enc)
                keep = rfe.get_support()
                fit_X = fit_X[:, keep]
                selected_columns = np.asarray(selected_columns)[keep].tolist()
            feature_selection = "lasso_svm"

        model = build_model(chosen_model)
        min_class = min(pd.Series(y_enc).value_counts())
        n_splits = max(2, min(5, min_class))
        cv = StratifiedKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=42,
        )
        scores = cross_val_score(
            model,
            fit_X,
            y_enc,
            cv=cv,
            scoring="accuracy",
        )
        y_pred = cross_val_predict(model, fit_X, y_enc, cv=cv)
        try:
            y_proba = cross_val_predict(
                model,
                fit_X,
                y_enc,
                cv=cv,
                method="predict_proba",
            )
        except Exception:
            model.fit(fit_X, y_enc)
            y_proba = model.predict_proba(fit_X)
        model.fit(fit_X, y_enc)

        proba = model.predict_proba(fit_X)
        sample_ids = (
            features["sample"]
            if "sample" in features.columns
            else features.index
        )
        confidence = pd.DataFrame({
            "sample": sample_ids,
            "condition": features["condition"],
            "predicted_condition": le.inverse_transform(model.predict(fit_X)),
            "confidence": proba.max(axis=1),
        })
        confidence.to_csv(
            data_dir / "fig_43_44_45_ml_classification_results.csv",
            index=False,
        )

        pd.DataFrame({"feature": selected_columns}).to_csv(
            data_dir / "fig_24_ml_selected_features.csv",
            index=False,
        )
        importance = _importance_values(
            model,
            fit_X,
            y_enc,
            selected_columns,
        )
        importance.to_csv(data_dir / "fig_24_ml_feature_importance.csv")

        fig, ax = plt.subplots(figsize=(8, 6))
        top = importance.head(15).iloc[::-1]
        ax.barh(top.index, top.values, color="#4C72B0")
        ax.set_title("ML feature importance")
        ax.set_xlabel("Importance")
        fig.tight_layout()
        fig.savefig(fig_dir / "fig_24_ml_feature_importance.png", dpi=150)
        plt.close(fig)

        try:
            import shap

            if chosen_model in {"xgb", "rf", "gbm"}:
                explainer = shap.TreeExplainer(model)
                shap_values = explainer.shap_values(fit_X)
            else:
                background = fit_X[: min(20, fit_X.shape[0])]
                explainer = shap.KernelExplainer(
                    model.predict_proba,
                    background,
                )
                shap_values = explainer.shap_values(
                    fit_X[: min(50, fit_X.shape[0])]
                )
            if isinstance(shap_values, list):
                shap_values = shap_values[1]
            shap.summary_plot(
                shap_values,
                fit_X[: shap_values.shape[0]],
                show=False,
                max_display=15,
            )
            plt.savefig(fig_dir / "fig_25_ml_shap.png", dpi=150, bbox_inches="tight")
            plt.close()
        except Exception:
            (data_dir / "ml_shap_status.txt").write_text(
                "SHAP plot skipped",
                encoding="utf-8",
            )

        report_dict = classification_report(
            y_enc,
            y_pred,
            output_dict=True,
            zero_division=0,
        )
        report_rows = []
        for label, metrics in report_dict.items():
            if isinstance(metrics, dict):
                row = {"class": label}
                row.update(metrics)
                report_rows.append(row)
        pd.DataFrame(report_rows).to_csv(
            data_dir / "fig_43_44_45_ml_classification_report.csv",
            index=False,
        )

        cm = confusion_matrix(y_enc, y_pred)
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(cm, cmap="Blues")
        ax.set_xticks(range(len(le.classes_)))
        ax.set_yticks(range(len(le.classes_)))
        ax.set_xticklabels(le.classes_, rotation=45, ha="right")
        ax.set_yticklabels(le.classes_)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        ax.set_title("Cross-validated confusion matrix")
        for i in range(len(le.classes_)):
            for j in range(len(le.classes_)):
                color = (
                    "white"
                    if cm[i, j] > cm.max() / 2
                    else "black"
                )
                ax.text(j, i, str(cm[i, j]), ha="center", va="center", color=color)
        fig.colorbar(im, ax=ax, fraction=0.046)
        fig.tight_layout()
        fig.savefig(fig_dir / "fig_43_ml_confusion_matrix.png", dpi=150)
        plt.close(fig)

        n_classes = len(le.classes_)
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        if n_classes == 2:
            fpr, tpr, _ = roc_curve(y_enc, y_proba[:, 1])
            roc_auc = auc(fpr, tpr)
            precision, recall, _ = precision_recall_curve(y_enc, y_proba[:, 1])
            ap_score = average_precision_score(y_enc, y_proba[:, 1])
            axes[0].plot(fpr, tpr, label=f"AUC={roc_auc:.3f}")
            axes[1].plot(recall, precision, label=f"AP={ap_score:.3f}")
        else:
            roc_auc = 0.0
            ap_score = 0.0
            for i in range(n_classes):
                y_bin = (y_enc == i).astype(int)
                fpr, tpr, _ = roc_curve(y_bin, y_proba[:, i])
                class_auc = auc(fpr, tpr)
                roc_auc += class_auc
                precision, recall, _ = precision_recall_curve(
                    y_bin, y_proba[:, i]
                )
                class_ap = average_precision_score(y_bin, y_proba[:, i])
                ap_score += class_ap
                axes[0].plot(
                    fpr,
                    tpr,
                    label=f"{le.classes_[i]} (AUC={class_auc:.3f})",
                )
                axes[1].plot(
                    recall,
                    precision,
                    label=f"{le.classes_[i]} (AP={class_ap:.3f})",
                )
            roc_auc /= n_classes
            ap_score /= n_classes
        axes[0].plot([0, 1], [0, 1], "--", color="grey")
        axes[0].set_xlabel("False positive rate")
        axes[0].set_ylabel("True positive rate")
        axes[0].set_title("ROC curve")
        axes[0].legend(fontsize=8)
        axes[1].set_xlabel("Recall")
        axes[1].set_ylabel("Precision")
        axes[1].set_title("Precision-recall curve")
        axes[1].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(fig_dir / "fig_44_ml_roc_pr.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.boxplot(scores, patch_artist=True)
        ax.scatter([1] * len(scores), scores, color="#4C72B0", zorder=3)
        ax.set_xticks([1])
        ax.set_xticklabels(["CV accuracy"])
        ax.set_ylabel("Accuracy")
        ax.set_title(f"Cross-validation scores (n={len(scores)})")
        fig.tight_layout()
        fig.savefig(fig_dir / "fig_45_ml_cv_scores.png", dpi=150)
        plt.close(fig)

        if n_classes == 2:
            try:
                n_bins = max(3, min(6, int(np.ceil(np.sqrt(len(y_enc))))))
                prob_true, prob_pred = calibration_curve(
                    y_enc,
                    y_proba[:, 1],
                    n_bins=n_bins,
                )
                fig, ax = plt.subplots(figsize=(6, 4.5))
                ax.plot(
                    prob_pred,
                    prob_true,
                    marker="o",
                    color="#4C72B0",
                    label="Calibration",
                )
                ax.plot(
                    [0, 1],
                    [0, 1],
                    "--",
                    color="grey",
                    label="Perfect",
                )
                ax.set_xlabel("Mean predicted probability")
                ax.set_ylabel("Observed frequency")
                ax.set_title("Calibration curve")
                ax.legend(fontsize=8)
                fig.tight_layout()
                fig.savefig(
                    fig_dir / "fig_45_ml_calibration_curve.png",
                    dpi=150,
                )
                plt.close(fig)
            except Exception:
                (data_dir / "ml_calibration_status.txt").write_text(
                    "Calibration curve skipped",
                    encoding="utf-8",
                )

        status = {
            "status": "completed",
            "model": chosen_model,
            "feature_selection": feature_selection,
            "selected_features": len(selected_columns),
            "accuracy": float(scores.mean()),
            "cv_mean": float(scores.mean()),
            "cv_sd": float(scores.std()),
            "auc": float(roc_auc),
            "average_precision": float(ap_score),
            "cv_scores": [float(s) for s in scores],
            "confusion_matrix": [
                [int(v) for v in row] for row in cm
            ],
            "n_samples": int(X.shape[0]),
            "features": int(X.shape[1]),
            "classes": le.classes_.tolist(),
        }
        (data_dir / "ml_model_summary.json").write_text(
            json.dumps(status, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    except Exception as exc:
        status["status"] = "failed"
        status["reason"] = str(exc)
        (data_dir / "ml_model_summary.json").write_text(
            json.dumps(status, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
