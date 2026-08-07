#!/usr/bin/env python3
"""Gather target and ligand evidence from public bioinformatics databases."""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import urllib.parse
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


def gather_evidence(cfg: ResolvedConfig, log) -> dict:
    ev = cfg.data.get("evidence", {})
    acc = ev.get("uniprot_accession") or ""
    pdb_id = ev.get("pdb_id") or ""
    chembl_id = ev.get("chembl_target_id") or ""
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
                    or ""
                ),
                "smiles": (
                    item.get("smiles")
                    or item.get("ligand_smiles")
                    or item.get("canonical_smiles")
                    or ""
                ),
                "activity": (
                    item.get("ki")
                    or item.get("ic50")
                    or item.get("k_d")
                    or item.get("k_d")
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
    for name, res in sections.items():
        ok = bool(res.get("ok", res.get("status_code", 0) == 200))
        lines.append(f"- `{name}`: {'ok' if ok else 'failed'}")
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
