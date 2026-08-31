#!/usr/bin/env python3
"""Random real-data validation for evidence collection and detect-box."""

import argparse
import csv
import json
import random
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from docking.box import detect_box_data  # noqa: E402

SKILLS = Path.home() / ".codex" / "skills"
SKILL_NAMES = {
    "uniprot": "uniprot-skill",
    "rcsb": "rcsb-pdb-skill",
    "chembl": "chembl-skill",
    "bindingdb": "bindingdb-skill",
    "pubchem": "pubchem-pug-skill",
    "chebi": "chebi-skill",
}
OUT_DIR = ROOT / "dock" / "validation_real_random"
POOL = [
    "1M17",
    "1XKK",
    "3CEJ",
    "4ASD",
    "1KV2",
    "1YET",
    "2HYY",
    "2BTR",
    "3OG7",
    "3PP0",
    "1AQ1",
    "1H1S",
    "2ITO",
    "3C4F",
    "4EKX",
    "2XHE",
    "4AG8",
    "1Y6B",
    "3LQ8",
    "5FDP",
]


def call_skill(skill: str, payload: dict, timeout: int = 180) -> dict:
    script = SKILLS / SKILL_NAMES[skill] / "scripts" / "rest_request.py"
    if not script.exists():
        return {
            "ok": False,
            "error": {"message": f"missing skill script: {script}"},
        }
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
    except FileNotFoundError as exc:
        return {"ok": False, "error": {"message": str(exc)}}
    except OSError as exc:
        return {"ok": False, "error": {"message": str(exc)}}
    except subprocess.TimeoutExpired:
        return {"ok": False, "error": {"message": f"timeout after {timeout}s"}}
    if not proc.stdout.strip():
        return {"ok": False, "error": {"message": proc.stderr or "empty output"}}
    try:
        return json.loads(proc.stdout)
    except Exception:
        return {"ok": False, "error": {"message": proc.stderr or "parse error"}}


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
        "max_items": 5,
        "max_depth": 3,
    }
    for attempt in range(3):
        result = call_skill("bindingdb", payload)
        if result.get("ok"):
            return result
        if attempt < 2:
            time.sleep(3 * (attempt + 1))
    return result


def download_pdb(pdb_id: str, dest: Path) -> None:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    req = urllib.request.Request(url, headers={"User-Agent": "Codex"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        dest.write_bytes(resp.read())


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Random real-data validation for evidence and detect-box."
    )
    parser.add_argument("--min-ok-targets", type=int, default=1)
    parser.add_argument("--min-box-ok", type=int, default=1)
    parser.add_argument("--min-ligands", type=int, default=1)
    args = parser.parse_args()

    random.seed(20260807)
    sample = random.sample(POOL, 10)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdb_dir = OUT_DIR / "pdb"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    summary: dict[str, dict] = {}

    print(f"{'PDB':<8}{'RCSB':<7}{'BindingDB':<10}{'box_mode':<18}{'center'}")
    for pdb_id in sample:
        rcsb = call_skill(
            "rcsb",
            {
                "base_url": "https://data.rcsb.org/rest/v1",
                "path": f"core/entry/{pdb_id}",
            },
        )
        bdb = fetch_bindingdb(pdb_id)
        ligand_count = len(bdb.get("records") or [])
        box: dict = {"center": None, "size": None, "mode": "failed", "error": ""}
        try:
            dest = pdb_dir / f"{pdb_id}.pdb"
            download_pdb(pdb_id, dest)
            center, size, mode = detect_box_data(dest)
            box = {"center": center, "size": size, "mode": mode, "error": ""}
        except Exception as exc:
            box["error"] = str(exc)
        rows.append(
            {
                "pdb": pdb_id,
                "rcsb_ok": bool(rcsb.get("ok")),
                "bindingdb_ligands": ligand_count,
                "box_mode": box["mode"],
                "center_x": box["center"][0] if box["center"] else "",
                "center_y": box["center"][1] if box["center"] else "",
                "center_z": box["center"][2] if box["center"] else "",
                "size": box["size"],
                "box_error": box["error"],
            }
        )
        summary[pdb_id] = {
            "rcsb_ok": bool(rcsb.get("ok")),
            "bindingdb_ligands": ligand_count,
            "box": box,
        }
        print(
            f"{pdb_id:<8}{str(rcsb.get('ok')).lower():<7}{ligand_count:<10}"
            f"{box['mode']:<18}{box['center']}"
        )

    csv_path = OUT_DIR / "results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "pdb",
                "rcsb_ok",
                "bindingdb_ligands",
                "box_mode",
                "center_x",
                "center_y",
                "center_z",
                "size",
                "box_error",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    ok_targets = sum(1 for item in summary.values() if item["rcsb_ok"])
    box_ok = sum(
        1
        for item in summary.values()
        if item["box"].get("mode") != "failed"
    )
    total_ligands = sum(r["bindingdb_ligands"] for r in rows)
    passed = (
        ok_targets >= args.min_ok_targets
        and box_ok >= args.min_box_ok
        and total_ligands >= args.min_ligands
    )
    (OUT_DIR / "summary.json").write_text(
        json.dumps(
            {
                "sample": sample,
                "summary": summary,
                "total_ligands": total_ligands,
                "ok_targets": ok_targets,
                "box_ok": box_ok,
                "passed": passed,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    print("OUTPUT", OUT_DIR)
    print("TOTAL_LIGANDS", total_ligands)
    print(
        "RESULT",
        "PASS" if passed else "FAIL",
        f"(ok_targets={ok_targets}/{args.min_ok_targets}, "
        f"box_ok={box_ok}/{args.min_box_ok}, "
        f"ligands={total_ligands}/{args.min_ligands})",
    )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
