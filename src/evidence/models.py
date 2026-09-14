"""Canonical evidence models used across database connectors."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import IntEnum
from typing import Any


class EvidenceTier(IntEnum):
    """Conservative ordering from computational prediction to direct assay."""

    PREDICTED = 1
    TEXT_MINED = 2
    CONTEXT = 3
    GENETIC = 4
    CURATED = 5
    EXPERIMENTAL = 6

    @classmethod
    def parse(cls, value: "EvidenceTier | int | str") -> "EvidenceTier":
        if isinstance(value, cls):
            return value
        if isinstance(value, int):
            return cls(value)
        key = str(value).strip().upper()
        if key.isdigit():
            return cls(int(key))
        return cls[key]


@dataclass(slots=True)
class EvidenceRecord:
    """One auditable evidence assertion.

    ``subject`` and ``object`` make the record usable for compound-target,
    target-disease, expression-context, dependency and pathway evidence
    without inventing one table per relationship type.
    """

    source: str
    source_record_id: str
    evidence_type: str
    subject_type: str
    subject_id: str
    relation: str
    object_type: str
    object_id: str
    target_symbol: str = ""
    tier: EvidenceTier = EvidenceTier.CURATED
    source_version: str = ""
    source_group: str = ""
    score: float | None = None
    effect_size: float | None = None
    p_value: float | None = None
    sample_size: int | None = None
    species: str = ""
    tissue: str = ""
    cell_type: str = ""
    assay: str = ""
    direction: str = "unknown"
    url: str = ""
    license: str = ""
    retrieved_at: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat()
    )
    payload: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        required = {
            "source": self.source,
            "source_record_id": self.source_record_id,
            "evidence_type": self.evidence_type,
            "subject_type": self.subject_type,
            "subject_id": self.subject_id,
            "relation": self.relation,
            "object_type": self.object_type,
            "object_id": self.object_id,
        }
        missing = [name for name, value in required.items() if not str(value).strip()]
        if missing:
            raise ValueError(f"evidence record missing fields: {', '.join(missing)}")

    @property
    def key(self) -> tuple[str, str, str, str, str, str]:
        return (
            self.source,
            self.source_record_id,
            self.evidence_type,
            self.subject_id,
            self.relation,
            self.object_id,
        )

    def fingerprint(self) -> str:
        payload = asdict(self)
        payload.pop("retrieved_at", None)
        payload["tier"] = int(self.tier)
        encoded = json.dumps(
            payload,
            sort_keys=True,
            ensure_ascii=False,
            default=str,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def to_row(self) -> dict[str, Any]:
        self.validate()
        row = asdict(self)
        row["tier"] = int(self.tier)
        row["payload_json"] = json.dumps(
            row.pop("payload"),
            ensure_ascii=False,
            sort_keys=True,
            default=str,
        )
        row["fingerprint"] = self.fingerprint()
        return row
