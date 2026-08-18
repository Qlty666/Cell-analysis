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

from analysis.ml_analysis import _importance_values, build_model  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
