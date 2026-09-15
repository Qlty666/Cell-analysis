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
from docking.provenance import sha256_file


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
    pipeline_scores_path: Path | None = None,
    allow_external_score: bool = False,
    target_column: str = "gene",
    score_column: str | None = "score",
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
        for column in (target_column, label_column)
        if column not in frame.columns
    ]
    if missing:
        raise ValueError(
            "external validation table missing columns: "
            + ", ".join(missing)
        )
    target_values = frame[target_column].astype(str).str.strip().str.upper()
    labels = frame[label_column]
    score_origin = "pipeline"
    score_provenance_note = "scores loaded from the pipeline ranking table"
    pipeline_path = Path(pipeline_scores_path) if pipeline_scores_path else None
    if (
        score_column
        and score_column in frame.columns
        and allow_external_score
    ):
        scores = pd.to_numeric(frame[score_column], errors="coerce")
        frame = frame[[target_column, label_column]].copy()
        frame.columns = ["target", "label"]
        frame["score"] = scores.to_numpy()
        score_origin = "external_allowed"
        score_provenance_note = (
            "external score column explicitly allowed; this is not an "
            "independent validation of the pipeline ranking"
        )
    elif pipeline_path is not None and pipeline_path.exists():
        ranked = pd.read_csv(pipeline_path)
        ranked_gene = next(
            (
                column
                for column in ("gene", "target_symbol")
                if column in ranked.columns
            ),
            None,
        )
        ranked_score = next(
            (
                column
                for column in (
                    "integrated_score",
                    "priority_score",
                    "target_score",
                    "score",
                )
                if column in ranked.columns
            ),
            None,
        )
        if ranked_gene is None or ranked_score is None:
            raise ValueError(
                "pipeline score table requires a gene and score column"
            )
        ranked = ranked[[ranked_gene, ranked_score]].rename(
            columns={ranked_gene: "target", ranked_score: "score"}
        )
        ranked["target"] = (
            ranked["target"].astype(str).str.strip().str.upper()
        )
        ranked["score"] = pd.to_numeric(ranked["score"], errors="coerce")
        joined = pd.DataFrame(
            {
                "target": target_values,
                "label": labels,
            }
        ).merge(
            ranked.drop_duplicates("target", keep="first"),
            on="target",
            how="left",
        )
        frame = joined.copy()
    elif score_column and score_column in frame.columns:
        summary = {
            "status": "skipped",
            "reason": (
                "external score column is present but allow_external_score=false; "
                "provide a pipeline score table instead"
            ),
            "path": str(path),
        }
        write_json(out_dir / "external_validation_summary.json", summary)
        return summary
    else:
        summary = {
            "status": "skipped",
            "reason": "no pipeline score table or external score column provided",
            "path": str(path),
        }
        write_json(out_dir / "external_validation_summary.json", summary)
        return summary
    frame["score"] = pd.to_numeric(frame["score"], errors="coerce")
    frame["label"] = frame["label"].map(_binary_label)
    target_count = int(len(frame))
    unmatched_targets = sorted(
        {
            str(value)
            for value in frame.loc[frame["score"].isna(), "target"].tolist()
            if str(value).strip()
        }
    )
    frame = frame.dropna(subset=["target", "score", "label"])
    match_rate = len(frame) / target_count if target_count else 0.0
    if frame.empty:
        summary = {
            "status": "skipped",
            "reason": "no external validation targets matched pipeline scores",
            "path": str(path),
            "score_origin": score_origin,
            "score_provenance_valid": False,
            "score_provenance_note": score_provenance_note,
            "target_match_count": 0,
            "target_match_rate": 0.0,
            "unmatched_targets": unmatched_targets,
        }
        write_json(out_dir / "external_validation_summary.json", summary)
        return summary
    if frame["label"].nunique() < 2:
        summary = {
            "status": "skipped",
            "reason": "external validation requires positive and negative labels",
            "path": str(path),
            "score_origin": score_origin,
            "score_provenance_valid": bool(score_origin == "pipeline"),
            "score_provenance_note": score_provenance_note,
            "target_match_count": int(len(frame)),
            "target_match_rate": float(match_rate),
            "unmatched_targets": unmatched_targets,
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
    precision_value = float(
        precision_score(labels, predictions, zero_division=0)
    )
    recall_value = float(
        recall_score(labels, predictions, zero_division=0)
    )
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
        "score_origin": score_origin,
        "score_provenance_valid": bool(score_origin == "pipeline"),
        "score_provenance_note": score_provenance_note,
        "pipeline_scores_path": (
            str(pipeline_path) if pipeline_path else ""
        ),
        "pipeline_scores_sha256": (
            sha256_file(pipeline_path)
            if pipeline_path is not None and pipeline_path.exists()
            else None
        ),
        "target_match_count": int(len(frame)),
        "target_match_rate": float(match_rate),
        "unmatched_targets": unmatched_targets,
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
        "precision": precision_value,
        "recall": recall_value,
        "f1": float(
            (
                2
                * precision_value
                * recall_value
            )
            / max(
                1e-12,
                precision_value + recall_value,
            )
        ),
        "specificity": float(tn / max(1, tn + fp)),
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
