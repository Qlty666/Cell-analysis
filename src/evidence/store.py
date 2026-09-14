"""SQLite-backed evidence store with source-run provenance."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Iterable, Mapping

import pandas as pd

from .models import EvidenceRecord, EvidenceTier


EVIDENCE_COLUMNS = [
    "source",
    "source_record_id",
    "evidence_type",
    "subject_type",
    "subject_id",
    "relation",
    "object_type",
    "object_id",
    "target_symbol",
    "tier",
    "source_version",
    "source_group",
    "score",
    "effect_size",
    "p_value",
    "sample_size",
    "species",
    "tissue",
    "cell_type",
    "assay",
    "direction",
    "url",
    "license",
    "retrieved_at",
    "payload_json",
    "fingerprint",
]


class SQLiteEvidenceStore:
    """Local evidence database suitable for reproducible desktop workflows."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path).expanduser().resolve()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._connection = sqlite3.connect(self.path)
        self._connection.row_factory = sqlite3.Row
        self._create_schema()

    def __enter__(self) -> "SQLiteEvidenceStore":
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        self.close()

    def close(self) -> None:
        self._connection.close()

    def _create_schema(self) -> None:
        self._connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS evidence_records (
                source TEXT NOT NULL,
                source_record_id TEXT NOT NULL,
                evidence_type TEXT NOT NULL,
                subject_type TEXT NOT NULL,
                subject_id TEXT NOT NULL,
                relation TEXT NOT NULL,
                object_type TEXT NOT NULL,
                object_id TEXT NOT NULL,
                target_symbol TEXT NOT NULL DEFAULT '',
                tier INTEGER NOT NULL,
                source_version TEXT NOT NULL DEFAULT '',
                source_group TEXT NOT NULL DEFAULT '',
                score REAL,
                effect_size REAL,
                p_value REAL,
                sample_size INTEGER,
                species TEXT NOT NULL DEFAULT '',
                tissue TEXT NOT NULL DEFAULT '',
                cell_type TEXT NOT NULL DEFAULT '',
                assay TEXT NOT NULL DEFAULT '',
                direction TEXT NOT NULL DEFAULT 'unknown',
                url TEXT NOT NULL DEFAULT '',
                license TEXT NOT NULL DEFAULT '',
                retrieved_at TEXT NOT NULL,
                payload_json TEXT NOT NULL DEFAULT '{}',
                fingerprint TEXT NOT NULL,
                recorded_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (
                    source,
                    source_record_id,
                    evidence_type,
                    subject_id,
                    relation,
                    object_id
                )
            );

            CREATE INDEX IF NOT EXISTS idx_evidence_target
                ON evidence_records(target_symbol);
            CREATE INDEX IF NOT EXISTS idx_evidence_source
                ON evidence_records(source);
            CREATE INDEX IF NOT EXISTS idx_evidence_tier
                ON evidence_records(tier);
            CREATE INDEX IF NOT EXISTS idx_evidence_relation
                ON evidence_records(relation);

            CREATE TABLE IF NOT EXISTS source_runs (
                run_id INTEGER PRIMARY KEY AUTOINCREMENT,
                source TEXT NOT NULL,
                source_version TEXT NOT NULL DEFAULT '',
                status TEXT NOT NULL,
                started_at TEXT NOT NULL,
                completed_at TEXT NOT NULL,
                record_count INTEGER NOT NULL DEFAULT 0,
                query_hash TEXT NOT NULL DEFAULT '',
                error TEXT NOT NULL DEFAULT '',
                metadata_json TEXT NOT NULL DEFAULT '{}'
            );
            """
        )
        self._connection.commit()

    def upsert_evidence(self, records: Iterable[EvidenceRecord]) -> int:
        rows = [record.to_row() for record in records]
        if not rows:
            return 0
        placeholders = ", ".join("?" for _ in EVIDENCE_COLUMNS)
        columns = ", ".join(EVIDENCE_COLUMNS)
        sql = (
            f"INSERT OR REPLACE INTO evidence_records ({columns}) "
            f"VALUES ({placeholders})"
        )
        self._connection.executemany(
            sql,
            [[row.get(column) for column in EVIDENCE_COLUMNS] for row in rows],
        )
        self._connection.commit()
        return len(rows)

    def record_source_run(
        self,
        *,
        source: str,
        source_version: str,
        status: str,
        started_at: str,
        completed_at: str,
        record_count: int,
        query_hash: str = "",
        error: str = "",
        metadata: Mapping | None = None,
    ) -> None:
        self._connection.execute(
            """
            INSERT INTO source_runs (
                source, source_version, status, started_at, completed_at,
                record_count, query_hash, error, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                source,
                source_version,
                status,
                started_at,
                completed_at,
                int(record_count),
                query_hash,
                error,
                json.dumps(metadata or {}, ensure_ascii=False, default=str),
            ),
        )
        self._connection.commit()

    def records(
        self,
        *,
        target_symbols: Iterable[str] | None = None,
        sources: Iterable[str] | None = None,
        relations: Iterable[str] | None = None,
    ) -> pd.DataFrame:
        clauses: list[str] = []
        parameters: list[str] = []
        if target_symbols is not None:
            values = sorted({str(value).upper() for value in target_symbols if value})
            if not values:
                return pd.DataFrame(columns=EVIDENCE_COLUMNS)
            clauses.append(f"target_symbol IN ({', '.join('?' for _ in values)})")
            parameters.extend(values)
        if sources is not None:
            values = sorted({str(value) for value in sources if value})
            if not values:
                return pd.DataFrame(columns=EVIDENCE_COLUMNS)
            clauses.append(f"source IN ({', '.join('?' for _ in values)})")
            parameters.extend(values)
        if relations is not None:
            values = sorted({str(value) for value in relations if value})
            if not values:
                return pd.DataFrame(columns=EVIDENCE_COLUMNS)
            clauses.append(f"relation IN ({', '.join('?' for _ in values)})")
            parameters.extend(values)
        where = f" WHERE {' AND '.join(clauses)}" if clauses else ""
        frame = pd.read_sql_query(
            f"SELECT * FROM evidence_records{where}",
            self._connection,
            params=parameters,
        )
        if not frame.empty:
            frame["tier"] = frame["tier"].map(EvidenceTier)
        return frame

    def source_runs(self) -> pd.DataFrame:
        return pd.read_sql_query(
            "SELECT * FROM source_runs ORDER BY run_id DESC",
            self._connection,
        )

    def delete_source_records(self, sources: Iterable[str]) -> int:
        """Remove records from sources that failed in the current run."""
        values = sorted({str(value) for value in sources if value})
        if not values:
            return 0
        placeholders = ", ".join("?" for _ in values)
        cursor = self._connection.execute(
            f"DELETE FROM evidence_records WHERE source IN ({placeholders})",
            values,
        )
        self._connection.commit()
        return int(cursor.rowcount or 0)

    def summary(self) -> dict:
        counts = pd.read_sql_query(
            """
            SELECT source, tier, COUNT(*) AS n
            FROM evidence_records
            GROUP BY source, tier
            ORDER BY source, tier
            """,
            self._connection,
        )
        totals = pd.read_sql_query(
            """
            SELECT COUNT(*) AS n_records,
                   COUNT(DISTINCT target_symbol) AS n_targets,
                   COUNT(DISTINCT source) AS n_sources
            FROM evidence_records
            """,
            self._connection,
        ).iloc[0]
        return {
            "database": str(self.path),
            "n_records": int(totals["n_records"]),
            "n_targets": int(totals["n_targets"]),
            "n_sources": int(totals["n_sources"]),
            "by_source_tier": counts.to_dict(orient="records"),
        }

    def export_csv(self, output_dir: str | Path) -> dict[str, str]:
        output = Path(output_dir).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        evidence_path = output / "evidence_records.csv"
        source_path = output / "source_runs.csv"
        coverage_path = output / "evidence_coverage.csv"
        self.records().to_csv(evidence_path, index=False)
        self.source_runs().to_csv(source_path, index=False)
        coverage = pd.read_sql_query(
            """
            SELECT source, source_group, tier, evidence_type, relation,
                   COUNT(*) AS n_records,
                   COUNT(DISTINCT target_symbol) AS n_targets
            FROM evidence_records
            GROUP BY source, source_group, tier, evidence_type, relation
            ORDER BY source, tier, evidence_type
            """,
            self._connection,
        )
        coverage.to_csv(coverage_path, index=False)
        return {
            "records": str(evidence_path),
            "source_runs": str(source_path),
            "coverage": str(coverage_path),
        }
