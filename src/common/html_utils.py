#!/usr/bin/env python3
"""Shared HTML escaping helpers for report generation."""

from __future__ import annotations

import html


def esc(value) -> str:
    """HTML-escape a value for safe embedding in generated reports.

    ``None`` renders as an empty string so missing values do not leak the
    literal text ``None`` into the report.
    """
    return html.escape(str(value if value is not None else ""))
