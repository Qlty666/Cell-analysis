#!/usr/bin/env python3
"""Shared HTML escaping helper for report writers."""

from __future__ import annotations

import html


def esc(value) -> str:
    """Return an HTML-escaped string for safe interpolation into markup."""
    return html.escape(str(value if value is not None else ""))
