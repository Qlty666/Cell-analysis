"""Shared checks for generated validation outputs."""

from __future__ import annotations

import json
from pathlib import Path

GSEA_KEGG_STATUS_REL = (
    Path("results")
    / "data"
    / "06_enrichment"
    / "fig_21_gsea_kegg_status.json"
)


def gsea_kegg_problem(root: Path) -> str | None:
    """Return a validation problem when KEGG GSEA did not complete."""
    path = root / GSEA_KEGG_STATUS_REL
    if not path.is_file():
        return str(GSEA_KEGG_STATUS_REL)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return f"{GSEA_KEGG_STATUS_REL}: invalid JSON"
    if not isinstance(payload, dict):
        return f"{GSEA_KEGG_STATUS_REL}: invalid payload"
    status = str(payload.get("status") or "").strip().lower()
    if status not in {"ok", "skipped"}:
        reason = str(payload.get("reason") or status or "missing status")
        return f"GSEA KEGG {reason}"
    return None
