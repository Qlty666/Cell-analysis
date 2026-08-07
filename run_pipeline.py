#!/usr/bin/env python3
"""CLI entry point for the liver cancer single-cell pipeline."""

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pipeline.orchestrator import main  # noqa: E402


if __name__ == "__main__":
    sys.exit(main())
