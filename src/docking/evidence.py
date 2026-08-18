#!/usr/bin/env python3
"""Gather target and ligand evidence from public bioinformatics databases."""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

from .config import ResolvedConfig
from .utils import write_json

SKILLS_ROOT = (
    Path(os.environ.get("CODEX_HOME", str(Path.home() / ".codex"))) / "skills"
)
SKILL_SCRIPTS = {
    "uniprot": SKILLS_ROOT / "uniprot-skill" / "scripts" / "rest_request.py",
    "rcsb": SKILLS_ROOT / "rcsb-pdb-skill" / "scripts" / "rest_request.py",
    "chembl": SKILLS_ROOT / "chembl-skill" / "scripts" / "rest_request.py",
    "bindingdb": SKILLS_ROOT / "bindingdb-skill" / "scripts" / "rest_request.py",
    "pubchem": SKILLS_ROOT / "pubchem-pug-skill" / "scripts" / "rest_request.py",
    "chebi": SKILLS_ROOT / "chebi-skill" / "scripts" / "rest_request.py",
    "string": SKILLS_ROOT / "string-skill" / "scripts" / "rest_request.py",
    "reactome": SKILLS_ROOT / "reactome-skill" / "scripts" / "rest_request.py",
    "pharmgkb": SKILLS_ROOT / "pharmgkb-skill" / "scripts" / "rest_request.py",
    "alphafold": SKILLS_ROOT / "alphafold-skill" / "scripts" / "rest_request.py",
    "opentargets": (
        SKILLS_ROOT / "opentargets-skill" / "scripts" / "opentargets_graphql.py"
    ),
}


def call_skill(name: str, payload: dict, timeout: int = 90, log=None) -> dict:
    script = SKILL_SCRIPTS.get(name)
    if not script or not script.exists():
        return {"ok": False, "error": {"message": f"skill script missing: {script}"}}
    try:
        proc = subprocess.run(
            [sys.executable, str(script)],
            input=json.dumps(payload),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
        )
        if proc.returncode != 0:
            return {"ok": False, "error": {"message": proc.stderr or proc.stdout}}
        return json.loads(proc.stdout or "{}")
    except Exception as exc:
        return {"ok": False, "error": {"message": str(exc)}}


def _http_text(url: str, timeout: int = 90) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", "replace")


def _kegg_pathways(gene: str, max_items: int = 10, timeout: int = 90) -> dict:
    """Resolve a gene to KEGG pathways using the public REST API."""
    if not gene:
        return {"ok": False, "error": {"message": "no gene identifier"}}
    try:
        find_url = "https://rest.kegg.jp/find/genes/" + urllib.parse.quote(gene)
        find_text = _http_text(find_url, timeout)
        kegg_gene = ""
        for line in find_text.splitlines():
            parts = line.split("\t")
            if parts and ":" in parts[0]:
                kegg_gene = parts[0]
                break
        if not kegg_gene:
            return {"ok": False, "error": {"message": "gene not found in KEGG"}}
        link_url = "https://rest.kegg.jp/link/pathway/" + kegg_gene
        link_text = _http_text(link_url, timeout)
        pathways: list[str] = []
        for line in link_text.splitlines():
            parts = line.split("\t")
            if len(parts) >= 2 and parts[1] not in pathways:
                pathways.append(parts[1])
        return {
            "ok": True,
            "kegg_gene": kegg_gene,
            "pathways": pathways[:max_items],
            "pathway_count": len(pathways),
        }
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": {"message": str(exc)}}


def _alphafold_evidence(
    uniprot: str,
    max_items: int = 10,
    timeout: int = 90,
    log=None,
) -> dict:
    if not uniprot:
        return {"ok": False, "error": {"message": "no UniProt accession"}}
    return call_skill(
        "alphafold",
        {
            "base_url": "https://alphafold.ebi.ac.uk/api",
            "path": f"prediction/{uniprot}",
            "max_items": max_items,
        },
        timeout=timeout,
        log=log,
    )


def _string_evidence(
    identifier: str,
    max_items: int = 10,
    timeout: int = 90,
    log=None,
) -> dict:
    if not identifier:
        return {"ok": False, "error": {"message": "no gene identifier"}}
    return call_skill(
        "string",
        {
            "base_url": "https://string-db.org/api/json",
            "path": "interaction_partners",
            "method": "POST",
            "form_body": {
                "identifier": identifier,
                "species": 9606,
                "caller_identity": "liver-cancer-pipeline",
                "limit": max(1, min(int(max_items), 10)),
            },
            "max_items": max_items,
        },
        timeout=timeout,
        log=log,
    )


def _reactome_evidence(
    identifier: str,
    max_items: int = 10,
    timeout: int = 90,
    log=None,
) -> dict:
    if not identifier:
        return {"ok": False, "error": {"message": "no gene identifier"}}
    return call_skill(
        "reactome",
        {
            "base_url": "https://reactome.org/ContentService",
            "path": "search/query",
            "params": {
                "query": identifier,
                "species": "Homo sapiens",
                "types": "Pathway",
            },
            "headers": {"Accept": "application/json"},
            "record_path": "results.0.entries",
            "max_items": max_items,
        },
        timeout=timeout,
        log=log,
    )


def _pharmgkb_evidence(
    gene: str,
    max_items: int = 10,
    timeout: int = 90,
    log=None,
) -> dict:
    if not gene:
        return {"ok": False, "error": {"message": "no gene name"}}
    return call_skill(
        "pharmgkb",
        {
            "base_url": "https://api.pharmgkb.org/v1/data",
            "path": "search",
            "params": {"q": gene, "type": "Gene", "limit": max_items},
            "max_items": max_items,
        },
        timeout=timeout,
        log=log,
    )


def _opentargets_evidence(
    gene: str,
    max_items: int = 10,
    timeout: int = 90,
    log=None,
) -> dict:
    if not gene:
        return {"ok": False, "error": {"message": "no gene name"}}
    return call_skill(
        "opentargets",
        {
            "query": (
                "query searchAny($q: String!) { "
                "search(queryString: $q) { total hits { entity score "
                "object { ... on Target { id approvedSymbol } } } } }"
            ),
            "variables": {"q": gene},
            "max_items": max_items,
            "max_depth": 6,
        },
        timeout=timeout,
        log=log,
    )


def _ids_from_records(res: dict, keys: list[str], limit: int = 10) -> list[str]:
    records = res.get("records") if isinstance(res, dict) else None
    if not isinstance(records, list):
        return []
    ids: list[str] = []
    for record in records:
        if not isinstance(record, dict):
            continue
        for key in keys:
            value = record.get(key)
            if value and isinstance(value, str) and value not in ids:
                ids.append(value)
                break
        if len(ids) >= limit:
            break
    return ids


def _opentargets_target_hits(res: dict, max_items: int = 10) -> tuple[int, list[str]]:
    search = (res.get("summary") or {}).get("search") or {}
    hits = search.get("hits") or []
    if not isinstance(hits, list):
        return 0, []
    ids: list[str] = []
    for hit in hits:
        if not isinstance(hit, dict) or hit.get("entity") != "target":
            continue
        obj = hit.get("object")
        if isinstance(obj, dict):
            value = obj.get("id") or obj.get("approvedSymbol") or ""
            if value and value not in ids:
                ids.append(value)
    return len(ids), ids[:max_items]


def collect_gene_database_evidence(
    gene: str,
    max_items: int = 10,
    timeout: int = 90,
    uniprot: str = "",
    ensembl: str = "",
) -> dict:
    """Collect compact per-gene counts for the full-pipeline evidence table."""
    sections = {
        "alphafold": _alphafold_evidence(uniprot, max_items, timeout),
        "string": _string_evidence(gene, max_items, timeout),
        "reactome": _reactome_evidence(gene, max_items, timeout),
        "pharmgkb": _pharmgkb_evidence(gene, max_items, timeout),
        "opentargets": _opentargets_evidence(gene, max_items, timeout),
        "kegg": _kegg_pathways(gene, max_items, timeout),
    }
    ot_hits, ot_ids = _opentargets_target_hits(
        sections["opentargets"],
        max_items,
    )
    kegg = sections["kegg"] or {}
    return {
        "string_partners": _record_count(sections["string"]),
        "string_partner_ids": ",".join(
            _ids_from_records(
                sections["string"],
                ["preferredName_B", "stringId_B"],
                max_items,
            )
        ),
        "reactome_pathways": _record_count(sections["reactome"]),
        "reactome_pathway_ids": ",".join(
            _ids_from_records(
                sections["reactome"],
                ["stId", "id"],
                max_items,
            )
        ),
        "pharmgkb_annotations": _record_count(sections["pharmgkb"]),
        "pharmgkb_ids": ",".join(
            _ids_from_records(
                sections["pharmgkb"],
                ["id", "accessionId"],
                max_items,
            )
        ),
        "alphafold_structures": _record_count(sections["alphafold"]),
        "alphafold_ids": ",".join(
            _ids_from_records(
                sections["alphafold"],
                ["entryId", "uniprotAccession", "id"],
                max_items,
            )
        ),
        "opentargets_hits": ot_hits,
        "opentargets_target_ids": ",".join(ot_ids),
        "kegg_pathways": int(kegg.get("pathway_count") or 0),
        "kegg_pathway_ids": ",".join(kegg.get("pathways") or []),
        "database_sources": ",".join(
            sorted(
                name
                for name, section in sections.items()
                if bool(
                    section.get("ok", section.get("status_code", 0) == 200)
                )
            )
        ),
    }


def gather_evidence(cfg: ResolvedConfig, log) -> dict:
    ev = cfg.data.get("evidence", {})
    acc = ev.get("uniprot_accession") or ""
    pdb_id = ev.get("pdb_id") or ""
    chembl_id = ev.get("chembl_target_id") or ""
    target_name = ev.get("target_name") or ""
    ligand_name = ev.get("ligand_name") or ""
    ligand_smiles = ev.get("ligand_smiles") or ""
    max_items = int(ev.get("max_items", 10))
    out_dir = cfg.workdir / "evidence"
    out_dir.mkdir(parents=True, exist_ok=True)

    sections: dict[str, dict] = {}
    ligand_rows: list[dict] = []

    if acc:
        res = call_skill(
            "uniprot",
            {
                "base_url": "https://rest.uniprot.org",
                "path": f"uniprotkb/{acc}",
                "params": {"format": "json"},
            },
            log=log,
        )
        sections["uniprot"] = res

    if pdb_id:
        res = call_skill(
            "rcsb",
            {
                "base_url": "https://data.rcsb.org/rest/v1",
                "path": f"core/entry/{pdb_id}",
            },
            log=log,
        )
        sections["rcsb"] = res

    chembl_res = {"ok": False, "error": {"message": "no ChEMBL target id"}}
    if not chembl_id and acc:
        for path in [
            f"target_component/uniprot/{acc}.json",
            f"target/uniprot/{acc}.json",
        ]:
            chembl_res = call_skill(
                "chembl",
                {
                    "base_url": "https://www.ebi.ac.uk/chembl/api/data",
                    "path": path,
                },
                log=log,
            )
            if chembl_res.get("ok"):
                summary = chembl_res.get("summary") or {}
                chembl_id = (
                    summary.get("target_chembl_id")
                    or summary.get("target_id")
                    or chembl_id
                )
                if chembl_id:
                    break

    if chembl_id:
        chembl_raw = out_dir / "chembl_activities_raw.json"
        res = call_skill(
            "chembl",
            {
                "base_url": "https://www.ebi.ac.uk/chembl/api/data",
                "path": "activity.json",
                "params": {
                    "target_chembl_id": chembl_id,
                    "limit": max_items,
                },
                "record_path": "activities",
                "max_items": max_items,
                "save_raw": True,
                "raw_output_path": str(chembl_raw),
            },
            log=log,
        )
        sections["chembl"] = res
        rows = _load_raw_list(chembl_raw, "activities") if chembl_raw.exists() else []
        for row in rows:
            ligand_rows.append(
                {
                    "source": "ChEMBL",
                    "ligand_id": row.get("molecule_chembl_id", ""),
                    "smiles": row.get("canonical_smiles", ""),
                    "activity": row.get("standard_value", ""),
                    "units": row.get("standard_units", ""),
                    "type": row.get("standard_type", ""),
                    "target": chembl_id,
                }
            )
    else:
        sections["chembl"] = chembl_res

    if acc:
        bdb_raw = out_dir / "bindingdb_uniprot_raw.json"
        res = call_skill(
            "bindingdb",
            {
                "base_url": "https://bindingdb.org",
                "path": "rest/getLigandsByUniprots",
                "params": {
                    "uniprot_id": acc,
                    "response": "application/json",
                },
                "max_items": max_items,
                "save_raw": True,
                "raw_output_path": str(bdb_raw),
            },
            log=log,
        )
        sections["bindingdb_uniprot"] = res
        ligand_rows.extend(_parse_bindingdb_raw(bdb_raw, acc))

    if pdb_id:
        bdb_pdb_raw = out_dir / "bindingdb_pdb_raw.json"
        res = call_skill(
            "bindingdb",
            {
                "base_url": "https://bindingdb.org",
                "path": "rest/getLigandsByPDBs",
                "params": {
                    "pdb": pdb_id,
                    "cutoff": 100,
                    "identity": 92,
                    "response": "application/json",
                },
                "max_items": max_items,
                "save_raw": True,
                "raw_output_path": str(bdb_pdb_raw),
            },
            log=log,
        )
        sections["bindingdb_pdb"] = res
        ligand_rows.extend(_parse_bindingdb_raw(bdb_pdb_raw, f"PDB:{pdb_id}"))

    if ligand_smiles:
        encoded = urllib.parse.quote(ligand_smiles, safe="")
        res = call_skill(
            "pubchem",
            {
                "base_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
                "path": (
                    "compound/smiles/"
                    f"{encoded}/property/MolecularFormula,"
                    "MolecularWeight,IUPACName,CanonicalSMILES/JSON"
                ),
                "record_path": "PropertyTable.Properties",
            },
            log=log,
        )
        sections["pubchem"] = res

    if ligand_name:
        res = call_skill(
            "chebi",
            {
                "base_url": "https://www.ebi.ac.uk/chebi/api",
                "path": "search",
                "params": {"query": ligand_name, "maximumRecords": max_items},
                "response_format": "json",
                "max_items": max_items,
            },
            log=log,
        )
        sections["chebi"] = res

    query_label = target_name or acc
    sections["alphafold"] = _alphafold_evidence(
        acc,
        max_items,
        log=log,
    )
    sections["string"] = _string_evidence(
        query_label,
        max_items,
        log=log,
    )
    sections["reactome"] = _reactome_evidence(
        query_label,
        max_items,
        log=log,
    )
    sections["pharmgkb"] = _pharmgkb_evidence(
        target_name,
        max_items,
        log=log,
    )
    sections["opentargets"] = _opentargets_evidence(
        target_name,
        max_items,
        log=log,
    )
    sections["kegg"] = _kegg_pathways(query_label, max_items)

    known_ligands = out_dir / "known_ligands.csv"
    if ligand_rows:
        with known_ligands.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "source",
                    "ligand_id",
                    "smiles",
                    "activity",
                    "units",
                    "type",
                    "target",
                ],
            )
            writer.writeheader()
            writer.writerows(ligand_rows)

    report = _render_report(ev, sections, ligand_rows, out_dir)
    (out_dir / "evidence_report.md").write_text(report, encoding="utf-8")
    write_json(out_dir / "evidence_summary.json", {"sections": sections})
    log.info(
        "evidence collected: %s sections, %s known ligands -> %s",
        len(sections),
        len(ligand_rows),
        out_dir,
    )
    return {"sections": len(sections), "known_ligands": len(ligand_rows)}


def _load_raw_list(path: Path, key: str) -> list:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return []
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return data.get(key) or []
    return []


def _load_raw_dict(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def _parse_bindingdb_raw(path: Path, target: str) -> list[dict]:
    data = _load_raw_dict(path)
    if not data:
        return []
    resp = data
    for key in [
        "getLindsByPDBsResponse",
        "getLigandsByPDBsResponse",
        "getLindsByUniprotsResponse",
        "getLigandsByUniprotsResponse",
    ]:
        if key in data and isinstance(data[key], dict):
            resp = data[key]
            break
    rows: list[dict] = []
    for item in (resp.get("affinities") or [])[:100]:
        if not isinstance(item, dict):
            continue
        rows.append(
            {
                "source": "BindingDB",
                "ligand_id": (
                    item.get("ligand_id")
                    or item.get("ligand_name")
                    or item.get("monomerid")
                    or ""
                ),
                "smiles": (
                    item.get("smiles")
                    or item.get("ligand_smiles")
                    or item.get("canonical_smiles")
                    or item.get("smile")
                    or ""
                ),
                "activity": (
                    item.get("ki")
                    or item.get("ic50")
                    or item.get("k_d")
                    or item.get("affinity")
                    or ""
                ),
                "units": item.get("units") or item.get("unit") or "nM",
                "type": (
                    item.get("measurement_type")
                    or item.get("binding_type")
                    or "Ki/IC50"
                ),
                "target": target,
            }
        )
    return rows


def _render_report(ev: dict, sections: dict, rows: list[dict], out_dir: Path) -> str:
    lines = [
        "# CADD Evidence Report",
        "",
        f"- Target name: {ev.get('target_name') or 'not set'}",
        f"- UniProt: {ev.get('uniprot_accession') or 'not set'}",
        f"- PDB: {ev.get('pdb_id') or 'not set'}",
        f"- ChEMBL target: {ev.get('chembl_target_id') or 'auto'}",
        f"- Known ligands: {len(rows)}",
        "",
        "## Sections",
        "",
    ]
    if sections:
        lines += [
            "| database | status | records |",
            "|---|---|---|",
        ]
        for name, res in sections.items():
            ok = bool(res.get("ok", res.get("status_code", 0) == 200))
            lines.append(
                f"| `{name}` | {'ok' if ok else 'failed'} | "
                f"{_record_count(res)} |"
            )
    else:
        lines.append("- no database sections requested")
    if rows:
        lines += [
            "",
            "## Known ligands",
            "",
            "| source | id | activity | units | smiles |",
            "|---|---|---|---|---|",
        ]
        for row in rows[:20]:
            lines.append(
                f"| {row['source']} | {row['ligand_id']} | "
                f"{row['activity']} | {row['units']} | {row['smiles']} |"
            )
    lines += [
        "",
        "## Suggested use",
        "",
        "Use known active ligands as positive controls during docking "
        "calibration and as training labels for ML/DL rescoring.",
        "",
    ]
    return "\n".join(lines)


def _record_count(res: dict) -> int:
    """Best-effort record count from a compact skill response."""
    if not isinstance(res, dict):
        return 0
    records = res.get("records")
    if isinstance(records, list):
        return len(records)
    if isinstance(res.get("pathway_count"), int):
        return res["pathway_count"]
    summary = res.get("summary") or {}
    if isinstance(summary, dict):
        for key in ("total", "count", "hits", "pathways"):
            value = summary.get(key)
            if isinstance(value, int):
                return value
            if isinstance(value, str) and value.isdigit():
                return int(value)
            if isinstance(value, list):
                return len(value)
    return 0
