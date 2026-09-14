#!/usr/bin/env python3
"""Tests for the multi-source evidence store and prioritisation layer."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from evidence.connectors import (
    ClinVarConnector,
    DepMapLocalConnector,
    ConnectorError,
    LocalTableConnector,
    OpenTargetsConnector,
    connector_from_config,
)
from evidence.context import EvidenceContext
from evidence.hub import EvidenceHub
from evidence.models import EvidenceRecord, EvidenceTier
from evidence.scoring import benchmark_ranking, score_targets
from evidence.store import SQLiteEvidenceStore


def _record(
    *,
    source: str,
    record_id: str,
    target: str,
    relation: str = "targets",
    evidence_type: str = "direct_bioactivity",
    tier: EvidenceTier = EvidenceTier.EXPERIMENTAL,
    score: float | None = 0.8,
    subject_type: str = "compound",
    object_type: str = "target",
    source_group: str | None = None,
) -> EvidenceRecord:
    return EvidenceRecord(
        source=source,
        source_record_id=record_id,
        evidence_type=evidence_type,
        subject_type=subject_type,
        subject_id="CHEMBL1",
        relation=relation,
        object_type=object_type,
        object_id=target,
        target_symbol=target,
        tier=tier,
        score=score,
        source_group=source_group or source.lower(),
        species="Homo sapiens",
    )


class TestEvidenceModelAndStore(unittest.TestCase):
    def test_record_validation_and_roundtrip(self):
        record = _record(source="ChEMBL", record_id="A1", target="EGFR")
        record.validate()
        self.assertEqual(record.key[0], "ChEMBL")
        self.assertEqual(record.to_row()["target_symbol"], "EGFR")

    def test_store_upsert_and_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with SQLiteEvidenceStore(root / "evidence.sqlite") as store:
                inserted = store.upsert_evidence(
                    [
                        _record(
                            source="ChEMBL",
                            record_id="A1",
                            target="EGFR",
                        ),
                        _record(
                            source="OpenTargets",
                            record_id="D1",
                            target="EGFR",
                            relation="associated_with",
                            evidence_type="curated_association",
                            tier=EvidenceTier.CURATED,
                            score=0.7,
                            subject_type="disease",
                            source_group="opentargets",
                        ),
                    ]
                )
                self.assertEqual(inserted, 2)
                frame = store.records(target_symbols=["EGFR"])
                self.assertEqual(len(frame), 2)
                summary = store.summary()
                self.assertEqual(summary["n_records"], 2)
                self.assertEqual(summary["n_targets"], 1)
                removed = store.delete_source_records(["ChEMBL"])
                self.assertEqual(removed, 1)
                self.assertEqual(store.summary()["n_records"], 1)
                paths = store.export_csv(root / "export")
                self.assertTrue(Path(paths["records"]).exists())


class TestEvidenceScoring(unittest.TestCase):
    def test_missing_evidence_is_not_zero(self):
        records = [
            _record(source="ChEMBL", record_id="A1", target="EGFR", score=0.9),
            _record(
                source="OpenTargets",
                record_id="D1",
                target="EGFR",
                relation="associated_with",
                evidence_type="genetic_association",
                tier=EvidenceTier.GENETIC,
                score=0.8,
                subject_type="disease",
            ),
            _record(
                source="SwissTargetPrediction",
                record_id="P1",
                target="GPAT3",
                tier=EvidenceTier.PREDICTED,
                score=0.7,
            ),
        ]
        evidence = pd.DataFrame([record.to_row() for record in records])
        priority, matrix, ablation = score_targets(evidence)
        ranked = priority.set_index("target_symbol")
        self.assertGreater(
            ranked.loc["EGFR", "priority_score"],
            ranked.loc["GPAT3", "priority_score"],
        )
        self.assertTrue(pd.isna(ranked.loc["EGFR", "predicted_target"]))
        self.assertTrue(pd.isna(ranked.loc["GPAT3", "direct_experimental"]))
        self.assertIn("direct_experimental", ranked.loc["GPAT3", "missing_categories"])
        self.assertIn("EGFR", set(matrix["target_symbol"]))
        self.assertTrue((ablation["ablation_runs"] > 1).all())

    def test_benchmark_recall_and_auroc(self):
        priority = pd.DataFrame(
            {
                "target_symbol": ["EGFR", "GPAT3", "PPP2R2A", "APP"],
                "priority_score": [0.9, 0.8, 0.7, 0.6],
            }
        )
        result = benchmark_ranking(
            priority,
            positive_targets=["GPAT3", "PPP2R2A"],
            negative_targets=["APP"],
            top_n=2,
        )
        self.assertEqual(result["recall_at_n"], 0.5)
        self.assertIsNotNone(result["auroc"])


class TestConnectors(unittest.TestCase):
    def test_connector_registry_rejects_unknown_name(self):
        with self.assertRaises(ConnectorError):
            connector_from_config("not-a-source", {})

    def test_open_targets_parser(self):
        search_payload = {
            "data": {
                "search": {
                    "hits": [
                        {
                            "id": "MONDO_0013209",
                            "name": "metabolic dysfunction-associated steatotic liver disease",
                            "entity": "disease",
                        }
                    ]
                }
            }
        }
        target_payload = {
            "data": {
                "disease": {
                    "associatedTargets": {
                        "rows": [
                            {
                                "score": 0.82,
                                "target": {
                                    "id": "ENSG00000146648",
                                    "approvedSymbol": "EGFR",
                                },
                            }
                        ]
                    }
                }
            }
        }
        with mock.patch(
            "evidence.connectors._post_json",
            side_effect=[search_payload, target_payload],
        ):
            records = OpenTargetsConnector().collect(
                EvidenceContext(
                    disease={"name": "NAFLD"},
                    max_records_per_source=10,
                )
            )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].target_symbol, "EGFR")
        self.assertEqual(records[0].tier, EvidenceTier.CURATED)

    def test_clinvar_parser(self):
        payload = {"esearchresult": {"count": "25"}}
        with mock.patch(
            "evidence.connectors._json_get",
            return_value=payload,
        ):
            records = ClinVarConnector().collect(
                EvidenceContext(
                    disease={"name": "liver cancer"},
                    target_symbols=["TP53"],
                    max_records_per_source=10,
                )
            )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].target_symbol, "TP53")
        self.assertEqual(records[0].tier, EvidenceTier.GENETIC)
        self.assertGreater(records[0].score or 0, 0)

    def test_local_table_connector(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "local.csv"
            pd.DataFrame(
                {
                    "gene": ["EGFR", "GPAT3"],
                    "score": [0.9, 0.7],
                }
            ).to_csv(path, index=False)
            connector = LocalTableConnector(
                name="LocalBioactivity",
                path=str(path),
                source_version="snapshot-1",
                evidence_type="direct_bioactivity",
                relation="targets",
                tier=EvidenceTier.EXPERIMENTAL,
                subject_type="compound",
                subject_value="CHEMBL1",
                target_column="gene",
                score_column="score",
            )
            records = connector.collect(
                EvidenceContext(max_records_per_source=10)
            )
            self.assertEqual([row.target_symbol for row in records], ["EGFR", "GPAT3"])
            self.assertEqual(records[0].source_group, "localbioactivity")

    def test_depmap_matrix_connector(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "CRISPRGeneEffect.csv"
            pd.DataFrame(
                {
                    "DepMap_ID": ["ACH-1", "ACH-2", "ACH-3"],
                    "EGFR (1956)": [-1.2, -0.8, -1.0],
                    "GPAT3 (84803)": [-0.1, -0.2, -0.15],
                }
            ).to_csv(path, index=False)
            records = DepMapLocalConnector(path=str(path)).collect(
                EvidenceContext(
                    target_symbols=["EGFR", "GPAT3"],
                    max_records_per_source=10,
                )
            )
            self.assertEqual(
                {row.target_symbol for row in records},
                {"EGFR", "GPAT3"},
            )
            egfr = next(row for row in records if row.target_symbol == "EGFR")
            self.assertGreaterEqual(egfr.score or 0, 0.5)


class TestEvidenceHub(unittest.TestCase):
    def test_hub_collect_score_and_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_path = root / "source.csv"
            pd.DataFrame(
                {
                    "gene": ["EGFR"],
                    "score": [0.9],
                }
            ).to_csv(source_path, index=False)
            config = {
                "strict": True,
                "sources": {
                    "Local": {
                        "enabled": True,
                        "name": "LocalBioactivity",
                        "path": str(source_path),
                        "evidence_type": "direct_bioactivity",
                        "relation": "targets",
                        "tier": "EXPERIMENTAL",
                        "subject_type": "compound",
                        "subject_value": "CHEMBL1",
                        "target_column": "gene",
                        "score_column": "score",
                    }
                },
                "scoring": {},
            }
            database = root / "evidence.sqlite"
            with SQLiteEvidenceStore(database) as store:
                hub = EvidenceHub.from_config(database, config, store=store)
                context = EvidenceContext(
                    compound={"chembl_id": "CHEMBL1"},
                    target_symbols=["EGFR"],
                    allow_network=False,
                )
                status = hub.collect(context)
                self.assertEqual(status.iloc[0]["record_count"], 1)
                scored = hub.score(target_symbols=["EGFR"])
                paths = hub.export(
                    root / "output",
                    context=context,
                    source_status=status,
                    scored=scored,
                )
            self.assertTrue(Path(paths["priority"]).exists())
            manifest = json.loads(
                Path(paths["run_manifest"]).read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["store_summary"]["n_records"], 1)

    def test_evidence_cli_offline_local_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source_path = root / "source.csv"
            pd.DataFrame(
                {"gene": ["EGFR"], "score": [0.9]}
            ).to_csv(source_path, index=False)
            config = {
                "strict": True,
                "sources": {
                    "Local": {
                        "enabled": True,
                        "name": "LocalBioactivity",
                        "path": str(source_path),
                        "evidence_type": "direct_bioactivity",
                        "relation": "targets",
                        "tier": "EXPERIMENTAL",
                        "subject_type": "compound",
                        "subject_value": "CHEMBL1",
                        "target_column": "gene",
                        "score_column": "score",
                    }
                },
                "scoring": {},
            }
            config_path = root / "sources.json"
            config_path.write_text(
                json.dumps(config),
                encoding="utf-8",
            )
            output = root / "output"
            result = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve().parents[1] / "scripts" / "run_evidence_hub.py"),
                    "--config",
                    str(config_path),
                    "--output",
                    str(output),
                    "--targets",
                    "EGFR",
                    "--offline",
                ],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((output / "target_priority.csv").exists())
            self.assertTrue((output / "evidence.sqlite").exists())


if __name__ == "__main__":
    unittest.main()
