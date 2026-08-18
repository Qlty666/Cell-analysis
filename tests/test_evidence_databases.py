#!/usr/bin/env python3
"""Unit tests for the expanded evidence database coverage."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

from docking import evidence as evidence_mod  # noqa: E402
from docking.config import ResolvedConfig  # noqa: E402


class TestExpandedDatabases(unittest.TestCase):
    def test_skill_scripts_include_expanded_database_skills(self):
        for name in (
            "string",
            "reactome",
            "pharmgkb",
            "alphafold",
            "opentargets",
        ):
            self.assertIn(name, evidence_mod.SKILL_SCRIPTS)
            self.assertTrue(
                evidence_mod.SKILL_SCRIPTS[name].exists(),
                msg=f"missing skill script for {name}",
            )

    def test_kegg_pathways_parses_rest_response(self):
        responses = [
            "hsa:7157\tTP53 tumor suppressor\n",
            "hsa:7157\tpath:hsa04115\nhsa:7157\tpath:hsa05200\n",
        ]
        with mock.patch.object(
            evidence_mod,
            "_http_text",
            side_effect=responses,
        ):
            result = evidence_mod._kegg_pathways("TP53", max_items=5)
        self.assertTrue(result["ok"])
        self.assertEqual(result["kegg_gene"], "hsa:7157")
        self.assertIn("path:hsa04115", result["pathways"])
        self.assertEqual(result["pathway_count"], 2)

    def test_gather_evidence_records_expanded_sections(self):
        with tempfile.TemporaryDirectory() as tmp:
            workdir = Path(tmp)
            cfg = ResolvedConfig(
                {
                    "workdir": str(workdir),
                    "evidence": {
                        "uniprot_accession": "P00533",
                        "target_name": "EGFR",
                        "ligand_smiles": "c1ccccc1",
                        "max_items": 3,
                    },
                },
                Path(tmp) / "cfg.json",
            )

            def fake_call(name, payload, timeout=90, log=None):
                return {"ok": True, "records": [{"id": name}]}

            with (
                mock.patch.object(
                    evidence_mod,
                    "call_skill",
                    side_effect=fake_call,
                ),
                mock.patch.object(
                    evidence_mod,
                    "_kegg_pathways",
                    return_value={
                        "ok": True,
                        "pathways": ["path:hsa04012"],
                        "pathway_count": 1,
                    },
                ),
            ):
                result = evidence_mod.gather_evidence(cfg, mock.Mock())

            summary = json.loads(
                (workdir / "evidence" / "evidence_summary.json").read_text(
                    encoding="utf-8"
                )
            )
            sections = summary["sections"]
            for name in (
                "string",
                "reactome",
                "pharmgkb",
                "alphafold",
                "opentargets",
                "kegg",
            ):
                self.assertIn(name, sections)
            self.assertEqual(result["known_ligands"], 0)

    def test_collect_gene_database_evidence_summarizes_sources(self):
        def fake_call(name, payload, timeout=90, log=None):
            if name == "opentargets":
                return {
                    "ok": True,
                    "summary": {
                        "search": {
                            "total": 50,
                            "hits": [
                                {
                                    "entity": "target",
                                    "object": {
                                        "id": "ENSG00000146648",
                                        "approvedSymbol": "EGFR",
                                    },
                                },
                                {"entity": "disease", "object": {}},
                            ],
                        }
                    },
                }
            if name == "string":
                return {
                    "ok": True,
                    "records": [
                        {
                            "preferredName_B": "GRB2",
                            "stringId_B": "9606.ENSP00000275493",
                        },
                        {
                            "preferredName_B": "EGF",
                            "stringId_B": "9606.ENSP00000265193",
                        },
                    ],
                }
            if name == "reactome":
                return {
                    "ok": True,
                    "records": [
                        {"stId": "R-HSA-177929", "id": "R-HSA-177929"}
                    ],
                }
            if name == "pharmgkb":
                return {
                    "ok": True,
                    "records": [{"id": "PA36679"}],
                }
            if name == "alphafold":
                return {
                    "ok": True,
                    "records": [{"entryId": "AF-P00533-F1"}],
                }
            return {"ok": True, "records": []}

        with (
            mock.patch.object(
                evidence_mod,
                "call_skill",
                side_effect=fake_call,
            ),
            mock.patch.object(
                evidence_mod,
                "_kegg_pathways",
                return_value={
                    "ok": True,
                    "pathways": ["path:hsa05200"],
                    "pathway_count": 1,
                },
            ),
        ):
            result = evidence_mod.collect_gene_database_evidence(
                "EGFR",
                max_items=5,
                uniprot="P00533",
            )

        self.assertEqual(result["string_partners"], 2)
        self.assertEqual(result["reactome_pathways"], 1)
        self.assertEqual(result["pharmgkb_annotations"], 1)
        self.assertEqual(result["alphafold_structures"], 1)
        self.assertEqual(result["opentargets_hits"], 1)
        self.assertEqual(result["kegg_pathways"], 1)
        self.assertIn("kegg", result["database_sources"])
        self.assertEqual(
            result["opentargets_target_ids"],
            "ENSG00000146648",
        )


if __name__ == "__main__":
    unittest.main()
