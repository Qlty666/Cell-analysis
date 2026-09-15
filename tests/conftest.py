#!/usr/bin/env python3
"""Shared pytest bootstrap for the project test suite.

Adds the repository root and the main source directories to ``sys.path`` once,
so tests can import project modules without repeating path surgery. Individual
test modules keep their own ``sys.path`` tweaks; this file is additive only.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

for _relative in ("", "src", "web", "launchers", "scripts"):
    _path = ROOT / _relative if _relative else ROOT
    _entry = str(_path)
    if _path.is_dir() and _entry not in sys.path:
        sys.path.insert(0, _entry)

del _relative, _path, _entry
