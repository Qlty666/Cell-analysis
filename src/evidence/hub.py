"""Evidence collection orchestration and reproducible export."""

from __future__ import annotations

import json
import platform
import sys
import time
from dataclasses import asdict
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd

from .connectors import BaseConnector, ConnectorError, connector_from_config
from .context import EvidenceContext
from .models import EvidenceTier
from .scoring import benchmark_ranking, score_targets
from .store import SQLiteEvidenceStore


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _version(package: str) -> str:
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return "not-installed"


class EvidenceHub:
    """Collect source outputs into one local database and score their targets."""

    def __init__(
        self,
        store: SQLiteEvidenceStore,
        connectors: Iterable[BaseConnector],
        *,
        scoring: Mapping | None = None,
        strict: bool = False,
    ) -> None:
        self.store = store
        self.connectors = list(connectors)
        self.scoring = dict(scoring or {})
        self.strict = bool(strict)

    @classmethod
    def from_config(
        cls,
        database: str | Path,
        config: Mapping[str, Any],
        *,
        store: SQLiteEvidenceStore | None = None,
    ) -> "EvidenceHub":
        source_config = config.get("sources") or {}
        connectors: list[BaseConnector] = []
        for name, raw_options in source_config.items():
            options = dict(raw_options or {})
            if isinstance(raw_options, bool):
                options = {}
                enabled = raw_options
            else:
                enabled = bool(options.get("enabled", True))
            if not enabled:
                continue
            connectors.append(connector_from_config(name, options))
        evidence_store = store or SQLiteEvidenceStore(database)
        return cls(
            evidence_store,
            connectors,
            scoring=config.get("scoring") or {},
            strict=bool(config.get("strict", False)),
        )

    def collect(self, context: EvidenceContext) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        for connector in self.connectors:
            started_at = _utc_now()
            started = time.time()
            status = "completed"
            error = ""
            records = []
            if not context.allow_network and connector.requires_network:
                status = "skipped"
                error = "network disabled"
            else:
                try:
                    records = connector.collect(context)
                    inserted = self.store.upsert_evidence(records)
                except Exception as exc:  # noqa: BLE001
                    if self.strict:
                        raise
                    status = "failed"
                    error = str(exc)
                    inserted = 0
                    records = []
            if not context.allow_network and connector.requires_network:
                inserted = 0
            completed_at = _utc_now()
            if status == "completed" and inserted == 0:
                status = "empty"
            self.store.record_source_run(
                source=connector.name,
                source_version=connector.source_version,
                status=status,
                started_at=started_at,
                completed_at=completed_at,
                record_count=inserted,
                query_hash=context.query_hash(connector.name),
                error=error,
                metadata={
                    "source_group": connector.source_group,
                    "tier": connector.tier.name,
                    "duration_seconds": round(time.time() - started, 3),
                },
            )
            rows.append(
                {
                    "source": connector.name,
                    "source_group": connector.source_group,
                    "tier": connector.tier.name,
                    "status": status,
                    "record_count": inserted,
                    "error": error,
                    "duration_seconds": round(time.time() - started, 3),
                }
            )
        return pd.DataFrame(rows)

    def score(
        self,
        *,
        target_symbols: Iterable[str] | None = None,
        benchmark_positives: Iterable[str] | None = None,
        benchmark_negatives: Iterable[str] | None = None,
        benchmark_top_n: int = 20,
        benchmark_source: str = "",
        benchmark_version: str = "",
        benchmark_independent: bool = False,
        benchmark_exclude_sources: Iterable[str] | None = None,
    ) -> dict[str, Any]:
        evidence = self.store.records(target_symbols=target_symbols)
        priority, matrix, ablation = score_targets(
            evidence,
            config=self.scoring,
        )
        result: dict[str, Any] = {
            "evidence": evidence,
            "priority": priority,
            "matrix": matrix,
            "ablation": ablation,
        }
        if benchmark_positives:
            source_column = (
                "source_group"
                if "source_group" in evidence.columns
                else "source"
                if "source" in evidence.columns
                else None
            )
            source_groups_by_target: dict[str, set[str]] = {}
            if source_column is not None and "target_symbol" in evidence.columns:
                for target, group in zip(
                    evidence["target_symbol"],
                    evidence[source_column],
                ):
                    target_key = str(target).upper().strip()
                    group_key = str(group).strip().lower()
                    if target_key and group_key:
                        source_groups_by_target.setdefault(
                            target_key,
                            set(),
                        ).add(group_key)
            benchmark = benchmark_ranking(
                priority,
                positive_targets=benchmark_positives,
                negative_targets=benchmark_negatives,
                top_n=benchmark_top_n,
                benchmark_source=benchmark_source,
                benchmark_version=benchmark_version,
                benchmark_independent=benchmark_independent,
                exclude_sources=benchmark_exclude_sources,
                source_groups_by_target=source_groups_by_target,
            )
            excluded_sources = {
                str(value).strip().lower()
                for value in (benchmark_exclude_sources or [])
                if str(value).strip()
            }
            if excluded_sources and source_column is not None:
                held_out = evidence[
                    ~evidence[source_column]
                    .fillna("")
                    .astype(str)
                    .str.lower()
                    .isin(excluded_sources)
                ]
                holdout_priority, _, _ = score_targets(
                    held_out,
                    config=self.scoring,
                )
                holdout = benchmark_ranking(
                    holdout_priority,
                    positive_targets=benchmark_positives,
                    negative_targets=benchmark_negatives,
                    top_n=benchmark_top_n,
                    benchmark_source=benchmark_source,
                    benchmark_version=benchmark_version,
                    benchmark_independent=True,
                    exclude_sources=excluded_sources,
                    source_groups_by_target=source_groups_by_target,
                )
                benchmark["holdout"] = holdout
                benchmark["benchmark_evidence_mode"] = "source_holdout"
                benchmark["benchmark_independent_verified"] = bool(
                    holdout.get("benchmark_independent_verified", False)
                )
            result["benchmark"] = benchmark
        return result

    def export(
        self,
        output_dir: str | Path,
        *,
        context: EvidenceContext,
        source_status: pd.DataFrame | None = None,
        scored: Mapping[str, Any] | None = None,
    ) -> dict[str, str]:
        output = Path(output_dir).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        store_paths = self.store.export_csv(output)
        scored = dict(
            scored
            or self.score(target_symbols=context.target_symbols)
        )
        outputs = dict(store_paths)
        outputs.update(
            {
                "priority": str(output / "target_priority.csv"),
                "matrix": str(output / "target_evidence_matrix.csv"),
                "source_ablation": str(output / "source_ablation.csv"),
                "source_status": str(output / "collector_status.csv"),
                "run_manifest": str(output / "evidence_run_manifest.json"),
            }
        )
        scored["priority"].to_csv(outputs["priority"], index=False)
        scored["matrix"].to_csv(outputs["matrix"], index=False)
        scored["ablation"].to_csv(outputs["source_ablation"], index=False)
        if source_status is None:
            status_frame = self.store.source_runs()
        else:
            status_frame = source_status
        status_frame.to_csv(outputs["source_status"], index=False)
        manifest = {
            "created_at": _utc_now(),
            "database": str(self.store.path),
            "context": {
                "compound": context.compound,
                "disease": context.disease,
                "target_symbols": sorted(
                    str(value).upper()
                    for value in context.target_symbols
                    if value
                ),
                "allow_network": context.allow_network,
                "max_records_per_source": context.max_records_per_source,
            },
            "connectors": [
                {
                    "name": connector.name,
                    "source_version": connector.source_version,
                    "source_group": connector.source_group,
                    "tier": connector.tier.name,
                }
                for connector in self.connectors
            ],
            "scoring": self.scoring
            or {"strategy": "default source-aware coverage-adjusted scoring"},
            "software": {
                "python": sys.version,
                "platform": platform.platform(),
                "numpy": _version("numpy"),
                "pandas": _version("pandas"),
            },
            "store_summary": self.store.summary(),
            "benchmark": scored.get("benchmark"),
        }
        Path(outputs["run_manifest"]).write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        return outputs
