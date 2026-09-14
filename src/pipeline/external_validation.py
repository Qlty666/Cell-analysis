"""External validation metrics for prioritized targets."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    precision_score,
    recall_score,
    roc_auc_score,
)

from docking.utils import write_json


def _binary_label(value) -> int | None:
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "positive", "case", "hcc"}:
        return 1
    if text in {"0", "false", "no", "negative", "control", "normal"}:
        return 0
    try:
        number = float(text)
    except ValueError:
        return None
    return int(number > 0.5)


def run_external_validation(
    path: Path | None,
    out_dir: Path,
    *,
    target_column: str = "gene",
    score_column: str = "score",
    label_column: str = "label",
    threshold: float = 0.5,
    bootstrap: int = 1000,
) -> dict:
    """Evaluate a held-out target table without mixing it into training."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if path is None or not Path(path).exists():
        summary = {
            "status": "skipped",
            "reason": "external validation table not provided",
            "path": str(path) if path else "",
        }
        write_json(out_dir / "external_validation_summary.json", summary)
        return summary
    frame = pd.read_csv(path)
    missing = [
        column
        for column in (target_column, score_column, label_column)
        if column not in frame.columns
    ]
    if missing:
        raise ValueError(
            "external validation table missing columns: "
            + ", ".join(missing)
        )
    frame = frame[[target_column, score_column, label_column]].copy()
    frame.columns = ["target", "score", "label"]
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    frame["label"] = frame["label"].map(_binary_label)
    frame = frame.dropna(subset=["target", "score", "label"])
    if frame["label"].nunique() < 2:
        summary = {
            "status": "skipped",
            "reason": "external validation requires positive and negative labels",
            "path": str(path),
        }
        write_json(out_dir / "external_validation_summary.json", summary)
        return summary
    labels = frame["label"].astype(int).to_numpy()
    scores = frame["score"].astype(float).to_numpy()
    auroc = float(roc_auc_score(labels, scores))
    auprc = float(average_precision_score(labels, scores))
    predictions = (scores >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(
        labels,
        predictions,
        labels=[0, 1],
    ).ravel()
    rng = np.random.default_rng(42)
    boot_auroc = []
    boot_auprc = []
    for _ in range(max(1, int(bootstrap))):
        indices = rng.integers(0, len(labels), size=len(labels))
        sampled_labels = labels[indices]
        sampled_scores = scores[indices]
        if len(set(sampled_labels)) < 2:
            continue
        boot_auroc.append(roc_auc_score(sampled_labels, sampled_scores))
        boot_auprc.append(
            average_precision_score(sampled_labels, sampled_scores)
        )
    summary = {
        "status": "completed",
        "path": str(path),
        "n_samples": int(len(frame)),
        "n_positive": int(labels.sum()),
        "n_negative": int(len(labels) - labels.sum()),
        "threshold": float(threshold),
        "auroc": auroc,
        "auroc_ci_low": (
            float(np.percentile(boot_auroc, 2.5)) if boot_auroc else None
        ),
        "auroc_ci_high": (
            float(np.percentile(boot_auroc, 97.5)) if boot_auroc else None
        ),
        "auprc": auprc,
        "auprc_ci_low": (
            float(np.percentile(boot_auprc, 2.5)) if boot_auprc else None
        ),
        "auprc_ci_high": (
            float(np.percentile(boot_auprc, 97.5)) if boot_auprc else None
        ),
        "precision": float(precision_score(labels, predictions, zero_division=0)),
        "recall": float(recall_score(labels, predictions, zero_division=0)),
        "confusion_matrix": {
            "tn": int(tn),
            "fp": int(fp),
            "fn": int(fn),
            "tp": int(tp),
        },
    }
    frame["prediction"] = predictions
    frame.to_csv(out_dir / "external_validation_predictions.csv", index=False)
    write_json(out_dir / "external_validation_summary.json", summary)
    return summary
