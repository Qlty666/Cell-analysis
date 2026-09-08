#!/usr/bin/env python3
"""Tests for the single-cell ML model expansion."""

from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path

import numpy as np

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from analysis.ml_analysis import (  # noqa: E402
    _LassoRfeSelector,
    _importance_values,
    build_model,
)


class TestModelBuilders(unittest.TestCase):
    def test_build_model_supported_types(self):
        from sklearn.ensemble import (
            GradientBoostingClassifier,
            RandomForestClassifier,
        )
        from sklearn.neural_network import MLPClassifier
        from sklearn.svm import SVC

        self.assertIsInstance(build_model("rf"), RandomForestClassifier)
        self.assertIsInstance(build_model("gbm"), GradientBoostingClassifier)
        self.assertIsInstance(build_model("mlp"), MLPClassifier)
        self.assertIsInstance(build_model("lasso_svm"), SVC)
        self.assertTrue(hasattr(build_model("xgb"), "fit"))

    def test_importance_values_uses_coefficients(self):
        model = types.SimpleNamespace(
            coef_=np.array([[1.0, -2.0, 0.5]])
        )
        importance = _importance_values(
            model,
            None,
            None,
            ["A", "B", "C"],
        )
        self.assertEqual(list(importance.index), ["B", "A", "C"])


class TestLassoRfeSelector(unittest.TestCase):
    """Selection must be a transformer so CV fits it per fold (no leakage)."""

    @staticmethod
    def _dataset(seed: int = 0, informative: bool = True):
        rng = np.random.default_rng(seed)
        n_samples, n_features = 40, 8
        X = rng.normal(size=(n_samples, n_features))
        y = (X[:, 0] + X[:, 1] > 0).astype(int) if informative else (
            rng.integers(0, 2, size=n_samples)
        )
        return X, y

    def test_transform_returns_subset_of_columns(self):
        X, y = self._dataset()
        columns = [f"f{index}" for index in range(X.shape[1])]
        selector = _LassoRfeSelector(random_state=42)
        selector.fit(X, y)
        transformed = selector.transform(X)
        self.assertEqual(transformed.shape[0], X.shape[0])
        self.assertGreaterEqual(transformed.shape[1], 1)
        self.assertLessEqual(transformed.shape[1], X.shape[1])
        selected = selector.selected_columns(columns)
        self.assertEqual(len(selected), transformed.shape[1])
        for name in selected:
            self.assertIn(name, columns)

    def test_selector_is_cv_compatible(self):
        from sklearn.model_selection import cross_val_score
        from sklearn.pipeline import Pipeline

        X, y = self._dataset(seed=1)
        pipeline = Pipeline(
            [
                ("select", _LassoRfeSelector(random_state=42)),
                ("clf", build_model("lasso_svm")),
            ]
        )
        scores = cross_val_score(pipeline, X, y, cv=3)
        self.assertEqual(len(scores), 3)
        self.assertTrue(np.all(scores >= 0.0))

    def test_selector_handles_all_zero_support(self):
        """A degenerate fit must still return at least one feature."""
        X = np.zeros((10, 5))
        y = np.array([0, 1] * 5)
        selector = _LassoRfeSelector(random_state=42)
        selector.fit(X, y)
        self.assertGreaterEqual(selector.transform(X).shape[1], 1)


if __name__ == "__main__":
    unittest.main()
