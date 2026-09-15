#!/usr/bin/env python3
"""Suite entry point that keeps liverbio importable as a normal package."""

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from liverbio_suite.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
