#!/usr/bin/env python3
"""ML/DL relevance ranking for GEO dataset search results."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent
if str(APP_ROOT) not in sys.path:
    sys.path.insert(0, str(APP_ROOT))

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import (
    GradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.neural_network import MLPClassifier

import search_datasets  # noqa: E402


def _value(row, name: str):
    if isinstance(row, dict):
        return row.get(name, "")
    return getattr(row, name, "")


def _text(disease: str, direction: str, row: dict) -> str:
    return " ".join(
        [
            disease,
            direction,
            str(_value(row, "accession")),
            str(_value(row, "title")),
            str(_value(row, "summary")),
            str(_value(row, "organism")),
        ]
    )


def load_samples(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = [
        "disease",
        "research_direction",
        "title",
        "summary",
        "label",
    ]
    missing = [name for name in required if name not in frame.columns]
    if missing:
        raise ValueError(f"training CSV missing columns: {missing}")
    frame["label"] = pd.to_numeric(frame["label"], errors="coerce")
    frame = frame.dropna(subset=["label"])
    frame["label"] = frame["label"].astype(int)
    frame["title"] = frame["title"].fillna("").astype(str)
    frame["summary"] = frame["summary"].fillna("").astype(str)
    frame["accession"] = frame.get("accession", pd.Series("", index=frame.index)).fillna("").astype(str)
    frame["organism"] = frame.get("organism", pd.Series("", index=frame.index)).fillna("").astype(str)
    return frame


def _build_model(model_type: str, seed: int):
    model_type = model_type.lower()
    if model_type == "lr":
        return LogisticRegression(max_iter=1000, random_state=seed)
    if model_type == "rf":
        return RandomForestClassifier(
            n_estimators=300,
            random_state=seed,
            n_jobs=-1,
        )
    if model_type == "gbm":
        return GradientBoostingClassifier(random_state=seed)
    if model_type == "mlp":
        return MLPClassifier(
            hidden_layer_sizes=(64, 32),
            activation="relu",
            max_iter=600,
            early_stopping=True,
            random_state=seed,
        )
    raise ValueError(f"unknown model type: {model_type}")


def _fit(frame: pd.DataFrame, model_type: str, seed: int):
    texts = [
        _text(row.disease, row.research_direction, row)
        for row in frame.itertuples(index=False)
    ]
    vectorizer = TfidfVectorizer(
        ngram_range=(1, 2),
        max_features=8000,
        sublinear_tf=True,
        stop_words="english",
    )
    x = vectorizer.fit_transform(texts)
    y = frame["label"].to_numpy()
    model = _build_model(model_type, seed)
    model.fit(x, y)
    return vectorizer, model


def train(
    samples_csv: Path,
    output_model: Path,
    model_type: str = "mlp",
    seed: int = 42,
) -> Path:
    frame = load_samples(samples_csv)
    vectorizer, model = _fit(frame, model_type, seed)
    payload = {
        "vectorizer": vectorizer,
        "model": model,
        "model_type": model_type,
        "classes": sorted(frame["label"].unique().tolist()),
    }
    output_model.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(payload, output_model)
    return output_model


def load_model(path: Path) -> dict:
    return joblib.load(path)


def lexical_scores(
    rows: list[dict],
    disease: str,
    direction: str,
) -> list[float]:
    disease_terms = [disease] + search_datasets._mapping_values(
        search_datasets.DISEASE_SYNONYMS,
        disease,
    )
    direction_terms = [direction] + search_datasets._mapping_values(
        search_datasets.DIRECTION_SYNONYMS,
        direction,
    )
    scores: list[float] = []
    for row in rows:
        text = " ".join(
            [
                str(row.get("accession", "")),
                str(row.get("title", "")),
                str(row.get("summary", "")),
            ]
        ).lower()
        disease_hits = sum(term.lower() in text for term in disease_terms)
        direction_hits = sum(term.lower() in text for term in direction_terms)
        disease_score = disease_hits / len(disease_terms) if disease_terms else 0.0
        direction_score = (
            direction_hits / len(direction_terms)
            if direction_terms
            else 0.0
        )
        scores.append(float(np.sqrt(disease_score * direction_score)))
    return scores


def rerank(
    rows: list[dict],
    disease: str,
    direction: str,
    model: dict | None = None,
    ml_weight: float = 0.6,
) -> list[dict]:
    if not rows:
        return []
    lexical = lexical_scores(rows, disease, direction)
    scores = [float(score) for score in lexical]
    if model is not None:
        frame = pd.DataFrame(rows)
        texts = [
            _text(disease, direction, row)
            for row in frame.itertuples(index=False)
        ]
        x = model["vectorizer"].transform(texts)
        proba = model["model"].predict_proba(x)
        positive = (
            1
            if model["model"].classes_[1] == 1
            else 0
        )
        ml_scores = proba[:, positive]
        scores = [
            ml_weight * float(ml_score)
            + (1.0 - ml_weight) * float(lex_score)
            for ml_score, lex_score in zip(ml_scores, lexical)
        ]
    ranked = [
        dict(row, relevance_score=score)
        for row, score in zip(rows, scores)
    ]
    return sorted(
        ranked,
        key=lambda item: item["relevance_score"],
        reverse=True,
    )


def evaluate(
    samples_csv: Path,
    model_type: str = "mlp",
    seed: int = 42,
    cv: int = 5,
) -> dict:
    frame = load_samples(samples_csv)
    if len(frame) < 4 or frame["label"].nunique() < 2:
        raise ValueError("not enough labeled samples for evaluation")
    texts = [
        _text(row.disease, row.research_direction, row)
        for row in frame.itertuples(index=False)
    ]
    y = frame["label"].to_numpy()
    vectorizer = TfidfVectorizer(
        ngram_range=(1, 2),
        max_features=8000,
        sublinear_tf=True,
        stop_words="english",
    )
    x = vectorizer.fit_transform(texts)
    y_true: list[int] = []
    y_pred: list[int] = []
    y_score: list[float] = []
    splitter = StratifiedKFold(
        n_splits=max(2, min(cv, min(np.bincount(y)))),
        shuffle=True,
        random_state=seed,
    )
    for train_idx, test_idx in splitter.split(x, y):
        model = _build_model(model_type, seed)
        model.fit(x[train_idx], y[train_idx])
        pred = model.predict(x[test_idx])
        proba = model.predict_proba(x[test_idx])
        positive = 1 if model.classes_[1] == 1 else 0
        y_true.extend(y[test_idx].tolist())
        y_pred.extend(pred.tolist())
        y_score.extend(proba[:, positive].tolist())
    accuracy = float(accuracy_score(y_true, y_pred))
    auc = (
        float(roc_auc_score(y_true, y_score))
        if len(set(y_true)) > 1
        else None
    )
    return {
        "model_type": model_type,
        "cv_folds": splitter.get_n_splits(),
        "samples": len(frame),
        "positive": int((y == 1).sum()),
        "negative": int((y == 0).sum()),
        "accuracy": round(accuracy, 4),
        "roc_auc": round(auc, 4) if auc is not None else None,
        "report": classification_report(
            y_true,
            y_pred,
            zero_division=0,
            output_dict=True,
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Train or evaluate ML/DL relevance ranking for GEO dataset search."
        )
    )
    parser.add_argument(
        "--train",
        action="store_true",
        help="train a model from labeled samples",
    )
    parser.add_argument(
        "--eval",
        action="store_true",
        help="cross-validate a model on labeled samples",
    )
    parser.add_argument("--samples", required=True)
    parser.add_argument(
        "--output",
        default=str(
            APP_ROOT / "data_cache" / "dataset_search" / "relevance_model.joblib"
        ),
    )
    parser.add_argument("--model-type", default="mlp")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--cv", type=int, default=5)
    args = parser.parse_args()

    if args.train:
        path = train(
            Path(args.samples),
            Path(args.output),
            model_type=args.model_type,
            seed=args.seed,
        )
        print(f"Model saved: {path}")
    if args.eval:
        result = evaluate(
            Path(args.samples),
            model_type=args.model_type,
            seed=args.seed,
            cv=args.cv,
        )
        out = Path(args.output).with_suffix(".eval.json")
        out.write_text(
            json.dumps(result, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        print(json.dumps(result, ensure_ascii=False, indent=2))
        print(f"Evaluation JSON: {out}")
    if not args.train and not args.eval:
        parser.error("provide --train or --eval")
    return 0


if __name__ == "__main__":
    sys.exit(main())
