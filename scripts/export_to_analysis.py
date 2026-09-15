#!/usr/bin/env python3
"""CLI entry point for exporting runs into a Codex analysis workspace."""

import sys
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from liverbio_suite.analysis_export import main  # noqa: E402


if __name__ == "__main__":
    sys.exit(main())
