"""Nested, optionally group-disjoint classification validation."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold


def splits(X, y, folds=5, seed=42, groups=None):
    y = np.asarray(y)
    if groups is None:
        count = int(pd.Series(y).value_counts().min())
        splitter = StratifiedKFold(n_splits=min(folds, count), shuffle=True, random_state=seed) if count >= 2 else None
    else:
        groups = np.asarray(groups)
        if pd.isna(groups).any():
            raise ValueError("CV group identifiers must not be missing")
        count = min(len(np.unique(groups[y == label])) for label in np.unique(y))
        splitter = StratifiedGroupKFold(n_splits=min(folds, count), shuffle=True, random_state=seed) if count >= 2 else None
    if splitter is None:
        raise ValueError("each class requires at least two independent samples/groups")
    result = list(splitter.split(X, y, groups))
    if any(len(np.unique(y[a])) < 2 or len(np.unique(y[b])) < 2 for a, b in result):
        raise ValueError("group-disjoint folds cannot contain both classes; adjust study design/folds")
    return result


def nested_predictions(models, X, y, positive_label, folds=5, seed=42, groups=None):
    """Outer predictions never use their labels for feature/model selection."""
    y = np.asarray(y)
    groups = None if groups is None else np.asarray(groups)
    probabilities = np.full(len(y), np.nan)
    records = []
    for fold, (train, test) in enumerate(splits(X, y, folds, seed, groups)):
        train_groups = None if groups is None else groups[train]
        try:
            inner = splits(X.iloc[train], y[train], min(3, folds), seed + fold, train_groups)
        except ValueError:
            # Do not use the outer test labels to choose a model in tiny cohorts.
            inner = []
        best_name, best_score = next(iter(models)), -np.inf
        for name, template in models.items():
            values = []
            for a, b in inner:
                model = clone(template).fit(X.iloc[train[a]], y[train[a]])
                p = model.predict_proba(X.iloc[train[b]])[:, list(model.classes_).index(positive_label)]
                values.append(roc_auc_score(y[train[b]] == positive_label, p))
            score = float(np.mean(values)) if values else -np.inf
            if score > best_score:
                best_name, best_score = name, score
        model = clone(models[best_name]).fit(X.iloc[train], y[train])
        probabilities[test] = model.predict_proba(X.iloc[test])[:, list(model.classes_).index(positive_label)]
        selector = model.named_steps.get("select")
        genes = list(X.columns[np.asarray(selector.get_support())]) if selector is not None else list(X.columns)
        records.append({"fold": fold, "model": best_name, "selection": "inner_cv" if inner else "predefined_first_model_insufficient_inner_samples", "train_samples": list(map(str, X.index[train])), "test_samples": list(map(str, X.index[test])), "features": genes})
    return probabilities, records
