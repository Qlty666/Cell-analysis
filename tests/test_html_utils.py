#!/usr/bin/env python3
"""Tests for the shared HTML escaping helper."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from common.html_utils import esc  # noqa: E402


class TestHtmlEsc(unittest.TestCase):
    def test_escapes_markup(self):
        self.assertEqual(
            esc("<script>alert('x')</script>"),
            "&lt;script&gt;alert(&#x27;x&#x27;)&lt;/script&gt;",
        )
        self.assertEqual(esc("a & b"), "a &amp; b")
        self.assertEqual(esc('"quoted"'), "&quot;quoted&quot;")

    def test_none_renders_empty(self):
        self.assertEqual(esc(None), "")

    def test_numbers_and_zero_are_preserved(self):
        self.assertEqual(esc(0), "0")
        self.assertEqual(esc(1.5), "1.5")


if __name__ == "__main__":
    unittest.main()
