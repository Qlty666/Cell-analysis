#!/usr/bin/env python3
"""ML disease classification, feature importance, and SHAP explainability."""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    data_dir = root / "results" / "data"
    fig_dir = root / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    status = {
        "status": "skipped",
        "reason": "",
        "model": "",
        "accuracy": None,
        "cv_scores": [],
    }

    try:
        ann_path = data_dir / "cell_annotations.csv"
        qc_path = data_dir / "qc_metrics.csv"
        if not ann_path.exists() or not qc_path.exists():
            status["reason"] = "missing cell_annotations.csv or qc_metrics.csv"
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
        confidence.to_csv(data_dir / "ml_classification_results.csv", index=False)

        importance = pd.Series(
            model.feature_importances_,
            index=X.columns,
            name="importance",
        ).sort_values(ascending=False)
        importance.to_csv(data_dir / "ml_feature_importance.csv")

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

        status = {
            "status": "completed",
            "model": "XGBoost" if use_xgb else "RandomForest",
            "accuracy": float(scores.mean()),
            "cv_scores": [float(s) for s in scores],
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
