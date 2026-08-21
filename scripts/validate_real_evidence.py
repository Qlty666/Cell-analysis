#!/usr/bin/env python3
"""Validate evidence collection against ten real PDB structures."""

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

SKILLS = Path.home() / ".codex" / "skills"
SKILL_NAMES = {
    "uniprot": "uniprot-skill",
    "rcsb": "rcsb-pdb-skill",
    "chembl": "chembl-skill",
    "bindingdb": "bindingdb-skill",
    "pubchem": "pubchem-pug-skill",
    "chebi": "chebi-skill",
    "string": "string-skill",
    "reactome": "reactome-skill",
    "pharmgkb": "pharmgkb-skill",
    "alphafold": "alphafold-skill",
}
ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "dock" / "validation_real"
TARGETS = [
    ("1M17", "EGFR"),
    ("1XKK", "BRAF"),
    ("3CEJ", "MET"),
    ("4ASD", "ALK"),
    ("1KV2", "MAPK14"),
    ("1YET", "HSP90AA1"),
    ("2HYY", "KDR"),
    ("2BTR", "CDK2"),
    ("3OG7", "BRAF"),
    ("3PP0", "MAP2K1"),
]


def call_skill(skill: str, payload: dict, timeout: int = 180) -> dict:
    script = SKILLS / SKILL_NAMES[skill] / "scripts" / "rest_request.py"
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
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": {"message": f"timeout after {timeout}s"}}
    if not proc.stdout.strip():
        return {"ok": False, "error": {"message": proc.stderr or "empty output"}}
    try:
        return json.loads(proc.stdout or "{}")
    except Exception:
        return {"ok": False, "error": {"message": proc.stderr or proc.stdout}}


def fetch_bindingdb(pdb_id: str) -> dict:
    payload = {
        "base_url": "https://bindingdb.org",
        "path": "rest/getLigandsByPDBs",
        "params": {
            "pdb": pdb_id,
            "cutoff": 100,
            "identity": 92,
            "response": "application/json",
        },
        "record_path": "getLindsByPDBsResponse.affinities",
        "max_items": 10,
        "max_depth": 3,
    }
    for attempt in range(3):
        result = call_skill("bindingdb", payload)
        if result.get("ok"):
            return result
        if attempt < 2:
            time.sleep(3 * (attempt + 1))
    return result


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    summary: dict[str, dict] = {}
    print(f"{'PDB':<8}{'target':<10}{'rcsb':<6}{'ligands'}")
    for pdb, target in TARGETS:
        rcsb = call_skill(
            "rcsb",
            {
                "base_url": "https://data.rcsb.org/rest/v1",
                "path": f"core/entry/{pdb}",
            },
        )
        bdb = fetch_bindingdb(pdb)
        records = bdb.get("records") or []
        for record in records:
            if not isinstance(record, dict):
                continue
            rows.append(
                {
                    "pdb": pdb,
                    "target": target,
                    "source": "BindingDB",
                    "ligand_id": record.get("monomerid", ""),
                    "smiles": record.get("smile", ""),
                    "query": record.get("query", ""),
                }
            )
        summary[pdb] = {
            "target": target,
            "rcsb_ok": bool(rcsb.get("ok")),
            "bindingdb_ok": bool(bdb.get("ok")),
            "ligands": len(records),
        }
        print(f"{pdb:<8}{target:<10}{str(rcsb.get('ok')).lower():<6}{len(records)}")

    csv_path = OUT_DIR / "known_ligands.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pdb", "target", "source", "ligand_id", "smiles", "query"],
        )
        writer.writeheader()
        writer.writerows(rows)
    (OUT_DIR / "summary.json").write_text(
        json.dumps(
            {"targets": summary, "total_ligands": len(rows)},
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print("TOTAL_LIGANDS", len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
