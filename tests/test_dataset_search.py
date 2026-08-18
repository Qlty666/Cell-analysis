#!/usr/bin/env python3
"""Tests for the automatic GEO dataset search script."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT) not in sys.path:
    sys.path.insert(0, str(APP_ROOT))
SCRIPTS_DIR = APP_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import search_datasets  # noqa: E402


class TestDatasetSearch(unittest.TestCase):
    def test_to_row_extracts_gse_accession(self):
        raw = {
            "Accession": "GDS1234 / GSE125449",
            "Title": "Liver cancer single cells",
            "Summary": "scRNA-seq",
            "taxon": "Homo sapiens",
            "n_samples": "8",
            "GPL": "GPL24676",
            "PDAT": "2024/01/01",
            "Type": "Expression profiling by high throughput sequencing",
        }
        row = search_datasets.to_row(
            raw,
            disease="liver cancer",
            research_direction="single cell",
        )
        self.assertEqual(row["accession"], "GSE125449")
        self.assertEqual(row["disease"], "liver cancer")
        self.assertEqual(row["research_direction"], "single cell")
        self.assertEqual(row["organism"], "Homo sapiens")
        self.assertEqual(row["samples"], "8")
        self.assertEqual(row["data_type"], "single-cell")
        self.assertIn("GSE125449", row["url"])

    def test_to_row_supports_lowercase_esummary_fields(self):
        raw = {
            "Accession": "GSE303421",
            "title": "HCC ChIP-Seq",
            "summary": "hepatocellular carcinoma",
            "taxon": "Homo sapiens",
            "n_samples": "12",
            "GPL": "34284",
            "PDAT": "2026/08/01",
            "gdsType": "Genome binding/occupancy profiling",
        }
        row = search_datasets.to_row(raw)
        self.assertEqual(row["title"], "HCC ChIP-Seq")
        self.assertEqual(row["type"], "Genome binding/occupancy profiling")
        self.assertEqual(row["samples"], "12")

    def test_to_row_classifies_bulk_rna_seq(self):
        row = search_datasets.to_row(
            {
                "Accession": "GSE299321",
                "Title": "TEAD4 in Gastric Cancer [RNA-Seq]",
                "Summary": "RNA-seq of gastric cancer cells",
                "taxon": "Homo sapiens",
            }
        )
        self.assertEqual(row["data_type"], "bulk")

    def test_to_row_defaults_to_other(self):
        row = search_datasets.to_row(
            {
                "Accession": "GSE999999",
                "Title": "Methylation profiling",
                "Summary": "array-based",
                "taxon": "Homo sapiens",
            }
        )
        self.assertEqual(row["data_type"], "other")

    def test_doc_summary_counts_samples(self):
        xml = (
            "<eSummaryResult><DocSum>"
            '<Item Name="Accession" Type="String">GSE1</Item>'
            '<Item Name="Samples" Type="List">'
            '<Item Name="Sample" Type="Structure"><Item Name="Accession">GSM1</Item></Item>'
            '<Item Name="Sample" Type="Structure"><Item Name="Accession">GSM2</Item></Item>'
            "</Item>"
            "</DocSum></eSummaryResult>"
        )
        import xml.etree.ElementTree as ET

        row = search_datasets._doc_summary(ET.fromstring(xml).find(".//DocSum"))
        self.assertEqual(row["n_samples"], "2")

    def test_filter_rows(self):
        rows = [
            {
                "accession": "GSE100001",
                "title": "HCC tumor",
                "summary": "human",
                "organism": "Homo sapiens",
                "samples": "4",
                "platform": "",
                "date": "",
                "type": "",
                "url": "",
            },
            {
                "accession": "GSE100002",
                "title": "Mouse liver",
                "summary": "mouse",
                "organism": "Mus musculus",
                "samples": "4",
                "platform": "",
                "date": "",
                "type": "",
                "url": "",
            },
        ]
        filtered = search_datasets.filter_rows(rows, organism="homo")
        self.assertEqual([r["accession"] for r in filtered], ["GSE100001"])
        filtered = search_datasets.filter_rows(rows, keyword="mouse")
        self.assertEqual([r["accession"] for r in filtered], ["GSE100002"])

    def test_filter_rows_extra_filters(self):
        rows = [
            {
                "accession": "GSE100001",
                "title": "HCC single-cell",
                "summary": "scRNA-seq",
                "data_type": "single-cell",
                "organism": "Homo sapiens",
                "samples": "8",
                "platform": "GPL24676",
                "date": "2024/01/01",
                "type": "Expression profiling by high throughput sequencing",
            },
            {
                "accession": "GSE100002",
                "title": "Methylation",
                "summary": "array",
                "data_type": "bulk",
                "organism": "Homo sapiens",
                "samples": "2",
                "platform": "GPL1",
                "date": "2020/05/05",
                "type": "Methylation profiling by array",
            },
            {
                "accession": "GSE100003",
                "title": "Liver single-cell",
                "summary": "10x",
                "data_type": "single-cell",
                "organism": "Homo sapiens",
                "samples": "20",
                "platform": "GPL24676",
                "date": "2025/06/01",
                "type": "Expression profiling by high throughput sequencing",
            },
        ]
        self.assertEqual(
            [
                r["accession"]
                for r in search_datasets.filter_rows(rows, data_type="single-cell")
            ],
            ["GSE100001", "GSE100003"],
        )
        self.assertEqual(
            [
                r["accession"]
                for r in search_datasets.filter_rows(
                    rows,
                    min_samples=5,
                    max_samples=10,
                )
            ],
            ["GSE100001"],
        )
        self.assertEqual(
            [
                r["accession"]
                for r in search_datasets.filter_rows(
                    rows,
                    start_date="2024-01-01",
                    end_date="2024-12-31",
                )
            ],
            ["GSE100001"],
        )
        self.assertEqual(
            [
                r["accession"]
                for r in search_datasets.filter_rows(
                    rows,
                    platform="GPL24676",
                )
            ],
            ["GSE100001", "GSE100003"],
        )
        self.assertEqual(
            [
                r["accession"]
                for r in search_datasets.filter_rows(
                    rows,
                    dataset_type="methylation",
                )
            ],
            ["GSE100002"],
        )

    def test_parse_date_handles_year_only(self):
        self.assertEqual(
            search_datasets._parse_date("2026"),
            (2026, 1, 1),
        )
        self.assertEqual(
            search_datasets._parse_date("2026/08/01"),
            (2026, 8, 1),
        )
        self.assertIsNone(search_datasets._parse_date("unknown"))

    def test_build_query(self):
        self.assertEqual(
            search_datasets.build_query(
                "liver cancer",
                "single cell RNA-seq",
            ),
            "liver cancer single cell RNA-seq",
        )
        self.assertEqual(
            search_datasets.build_query("", "", "custom term"),
            "custom term",
        )

    def test_search_datasets_mocks_ncbi_calls(self):
        summaries = {
            "1": {
                "Accession": "GSE100001",
                "Title": "HCC tumor",
                "Summary": "human liver",
                "taxon": "Homo sapiens",
                "n_samples": "4",
                "GPL": "GPL24676",
                "PDAT": "2024/01/01",
                "Type": "Expression profiling",
            },
            "2": {
                "Accession": "GSE100002",
                "Title": "Mouse liver",
                "Summary": "mouse",
                "taxon": "Mus musculus",
                "n_samples": "4",
                "GPL": "GPL1",
                "PDAT": "2024/02/01",
                "Type": "Expression profiling",
            },
        }
        with (
            patch.object(search_datasets, "esearch", return_value=["1", "2"]),
            patch.object(
                search_datasets,
                "esummary_many",
                return_value=summaries,
            ),
            patch.object(search_datasets.time, "sleep"),
        ):
            rows = search_datasets.search_datasets(
                "liver cancer",
                max_results=2,
                organism="homo",
                disease="liver cancer",
                research_direction="single cell",
            )
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["accession"], "GSE100001")
        self.assertEqual(rows[0]["disease"], "liver cancer")
        self.assertEqual(
            rows[0]["research_direction"],
            "single cell",
        )

    def test_write_outputs(self):
        rows = [
            {
                "accession": "GSE100001",
                "disease": "liver cancer",
                "research_direction": "single cell",
                "title": "HCC tumor",
                "summary": "human",
                "organism": "Homo sapiens",
                "samples": "4",
                "platform": "GPL24676",
                "date": "2024/01/01",
                "type": "Expression profiling",
                "url": "https://example.com/GSE100001",
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "results"
            csv_path, json_path = search_datasets.write_outputs(rows, out)
            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())
            text = csv_path.read_text(encoding="utf-8")
            self.assertIn("GSE100001", text)

    def test_write_outputs_includes_relevance_score(self):
        rows = [
            {
                "accession": "GSE100002",
                "relevance_score": 0.912,
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "results"
            csv_path, json_path = search_datasets.write_outputs(rows, out)
            text = csv_path.read_text(encoding="utf-8")
            self.assertIn("relevance_score", text)
            self.assertIn("0.912", text)


if __name__ == "__main__":
    unittest.main()
