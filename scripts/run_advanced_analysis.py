#!/usr/bin/env python3
"""CLI entry point for advanced bulk analysis."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from analysis.advanced_analysis import main  # noqa: E402


if __name__ == "__main__":
    sys.exit(main())
