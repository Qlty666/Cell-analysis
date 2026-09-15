#!/usr/bin/env python3
"""Tests for FAERS-style disproportionality signal detection."""

from __future__ import annotations

import logging
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from docking.config import load_config  # noqa: E402
from docking.signal_detection import detect_signals, run_faers  # noqa: E402

DEFAULT_CONFIG = APP_ROOT / "config" / "docking_config.json"
LOG = logging.getLogger("test_signal_detection")


class TestSignalDetection(unittest.TestCase):
    def test_known_signal(self):
        df = pd.DataFrame(
            {
                "drug": ["A", "A", "B", "B"],
                "event": ["X", "Y", "X", "Y"],
                "count": [10, 2, 2, 20],
            }
        )
        signals = detect_signals(df, min_count=3)
        row = signals[
            (signals["drug"] == "A") & (signals["event"] == "X")
        ].iloc[0]
        self.assertEqual(row["a"], 10)
        self.assertEqual(row["b"], 2)
        self.assertEqual(row["c"], 2)
        self.assertEqual(row["d"], 20)
        self.assertGreater(row["ror"], 10)
        self.assertGreater(row["ror_lower"], 1.0)
        self.assertTrue(bool(row["signal"]))

    def test_aggregates_duplicate_rows(self):
        df = pd.DataFrame(
            {
                "drug": ["A", "A", "B"],
                "event": ["X", "X", "Y"],
            }
        )
        signals = detect_signals(df)
        row = signals[
            (signals["drug"] == "A") & (signals["event"] == "X")
        ].iloc[0]
        self.assertEqual(row["a"], 2)

    def test_run_faers(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp) / "work"
            workdir.mkdir()
            data_dir = workdir / "data" / "faers"
            data_dir.mkdir(parents=True)
            pd.DataFrame(
                {
                    "drug": ["A", "A", "B", "B"],
                    "event": ["X", "Y", "X", "Y"],
                    "count": [10, 2, 2, 20],
                }
            ).to_csv(data_dir / "events.csv", index=False)
            cfg = load_config(
                DEFAULT_CONFIG,
                {
                    "workdir": str(workdir),
                    "faers_input": "data/faers/events.csv",
                    "faers_min_count": 3,
                },
            )
            summary = run_faers(cfg, LOG)
            out_dir = workdir / "outputs" / "run_001" / "faers"
            self.assertTrue((out_dir / "data" / "faers_signals.csv").exists())
            self.assertTrue((out_dir / "faers_summary.json").exists())
            self.assertGreaterEqual(summary["signals"], 1)


if __name__ == "__main__":
    unittest.main()
