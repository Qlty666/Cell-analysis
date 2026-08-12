#!/usr/bin/env python3
"""ML disease classification, feature importance, and SHAP explainability."""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
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
from sklearn.preprocessing import LabelEncoder


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

        use_xgb = True
        try:
            from xgboost import XGBClassifier

            model = XGBClassifier(
                n_estimators=200,
                max_depth=4,
                use_label_encoder=False,
                eval_metric="logloss",
                random_state=42,
                verbosity=0,
            )
        except Exception:
            use_xgb = False
            from sklearn.ensemble import RandomForestClassifier

            model = RandomForestClassifier(
                n_estimators=300,
                max_depth=8,
                random_state=42,
                n_jobs=-1,
            )

        min_class = min(pd.Series(y_enc).value_counts())
        n_splits = max(2, min(5, min_class))
        cv = StratifiedKFold(
            n_splits=n_splits,
            shuffle=True,
            random_state=42,
        )
        scores = cross_val_score(model, X, y_enc, cv=cv, scoring="accuracy")
        y_pred = cross_val_predict(model, X, y_enc, cv=cv)
        try:
            y_proba = cross_val_predict(
                model, X, y_enc, cv=cv, method="predict_proba"
            )
        except Exception:
            model.fit(X, y_enc)
            y_proba = model.predict_proba(X)
        model.fit(X, y_enc)

        proba = model.predict_proba(X)
        sample_ids = (
            features["sample"]
            if "sample" in features.columns
            else features.index
        )
        confidence = pd.DataFrame({
            "sample": sample_ids,
            "condition": features["condition"],
            "predicted_condition": le.inverse_transform(model.predict(X)),
            "confidence": proba.max(axis=1),
        })
        confidence.to_csv(
            data_dir / "fig_43_44_45_ml_classification_results.csv",
            index=False,
        )

        importance = pd.Series(
            model.feature_importances_,
            index=X.columns,
            name="importance",
        ).sort_values(ascending=False)
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

            explainer = shap.TreeExplainer(model)
            shap_values = explainer.shap_values(X)
            if isinstance(shap_values, list):
                shap_values = shap_values[1]
            shap.summary_plot(
                shap_values,
                X,
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

        status = {
            "status": "completed",
            "model": "XGBoost" if use_xgb else "RandomForest",
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
