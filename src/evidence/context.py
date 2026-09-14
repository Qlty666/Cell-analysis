"""Inputs shared by evidence-source connectors."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(slots=True)
class EvidenceContext:
    compound: dict[str, Any] = field(default_factory=dict)
    disease: dict[str, Any] = field(default_factory=dict)
    target_symbols: list[str] = field(default_factory=list)
    ensembl_ids: dict[str, str] = field(default_factory=dict)
    cache_dir: Path | None = None
    max_records_per_source: int = 1000
    timeout_seconds: int = 120
    allow_network: bool = True
    source_options: dict[str, dict[str, Any]] = field(default_factory=dict)

    def query_hash(self, source: str) -> str:
        payload = {
            "source": source,
            "compound": self.compound,
            "disease": self.disease,
            "target_symbols": sorted(
                str(value).upper() for value in self.target_symbols if value
            ),
            "ensembl_ids": {
                str(key).upper(): str(value)
                for key, value in sorted(self.ensembl_ids.items())
                if value
            },
        }
        encoded = json.dumps(
            payload,
            ensure_ascii=False,
            sort_keys=True,
            default=str,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()
