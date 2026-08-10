#!/usr/bin/env python3
"""End-to-end automation from single-cell analysis to virtual screening.

The module wires the existing pieces together:

1. run the GEO single-cell pipeline and export a sample-level pseudobulk matrix;
2. rank significant DEGs into a compact key-gene table;
3. enrich genes with UniProt/PDB/ChEMBL evidence (network optional, cached);
4. build the virtual-knockout inputs and run multidimensional target scoring;
5. for genes with a PDB structure, collect known ligands and run the full
   AutoDock Vina pipeline in an isolated per-target workdir;
6. export the wet-lab validation plan and an integrated HTML report.

Every stage writes a marker file so a rerun resumes where it stopped.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

APP_ROOT = Path(__file__).resolve().parent.parent.parent
SRC = APP_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from docking import box, evidence as evidence_mod, pipeline as docking_pipeline  # noqa: E402
from docking.config import load_config, save_config  # noqa: E402
from docking.knockout import run_knockout  # noqa: E402
from docking.provenance import write_run_manifest  # noqa: E402
from docking.utils import DockingError, ToolNotFoundError, safe_name, write_json  # noqa: E402
from docking.validation import export_validation  # noqa: E402

from . import orchestrator  # noqa: E402

log = logging.getLogger("full_pipeline")

STAGES = [
    ("01", "single_cell", "GEO single-cell analysis (download, QC, annotation, DEG)"),
    ("02", "key_targets", "extract and rank key genes/proteins from DEGs"),
    ("03", "evidence", "enrich genes with UniProt/PDB/ChEMBL evidence"),
    ("04", "knockout_inputs", "build pseudobulk expression and knockout inputs"),
    ("05", "knockout", "virtual knockout and multidimensional target scoring"),
    ("06", "docking", "per-target virtual screening with AutoDock Vina"),
    ("07", "report", "integrated HTML report and provenance manifest"),
]

DEFAULT_GENE_BLACKLIST = [
    r"^RPL",
    r"^RPS",
    r"^MRPL",
    r"^MRPS",
    r"^MT-",
    r"^MTRNR",
    r"^SNORD",
    r"^SCGB",
    r"^IGH",
    r"^IGK",
    r"^IGL",
    r"^TRA",
    r"^TRB",
    r"^TRG",
    r"^HLA-D",
    r"^LINC",
    r"^RP[0-9]",
    r"^AC[0-9]",
    r"^AL[0-9]",
]

EVIDENCE_COLUMNS = [
    "gene",
    "entrez",
    "uniprot",
    "ensembl",
    "chembl_target_id",
    "known_ligands",
    "chembl_bioactivities",
    "pdb_structures",
    "pdb_ids",
    "off_target_paralogs",
    "safety_concern",
]


class IntegrationError(RuntimeError):
    """Raised when the integrated pipeline cannot continue."""


class PauseRequested(Exception):
    """Raised when the underlying single-cell pipeline pauses."""


def _integration_dir(workdir: Path) -> Path:
    return workdir / "outputs" / "integration"


def _stage_dir(workdir: Path) -> Path:
    return _integration_dir(workdir) / ".stages"


def _marker(workdir: Path, code: str, name: str) -> Path:
    return _stage_dir(workdir) / f"{code}_{name}.done"


def _read_json(path: Path, default=None) -> dict:
    if not path.exists():
        return default or {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else default or {}
    except Exception:
        return default or {}


def _single_cell_outputs_ready(root: Path) -> bool:
    """True when the single-cell stage produced the files later stages need."""
    return (
        (root / "results" / "pipeline_complete.json").exists()
        and (
            (root / "results" / "data" / "deg_significant.csv").exists()
            or (root / "results" / "data" / "deg_all.csv").exists()
        )
    )


def _resolve_path(value, base: Path) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def extract_key_genes(
    single_cell_root: Path,
    out_dir: Path,
    top_n: int = 50,
    keep_all: bool = False,
    blacklist_patterns: list[str] | None = None,
) -> pd.DataFrame:
    """Rank significant DEGs into a compact key-gene table."""
    data_dir = single_cell_root / "results" / "data"
    deg_path = data_dir / "deg_significant.csv"
    if not deg_path.exists():
        deg_path = data_dir / "deg_all.csv"
    if not deg_path.exists():
        raise IntegrationError(f"DEG table not found under {data_dir}")

    frame = pd.read_csv(deg_path)
    if frame.empty:
        raise IntegrationError(f"DEG table is empty: {deg_path}")

    rename = {
        "avg_log2FC": "avg_log2fc",
        "log2FoldChange": "avg_log2fc",
        "pvalue": "p_val",
        "padj": "p_val_adj",
    }
    frame = frame.rename(columns=rename)
    # Some DEG exports contain both avg_log2FC and log2FoldChange. After the
    # rename both become avg_log2fc, which turns column access into a DataFrame.
    frame = frame.loc[:, ~frame.columns.duplicated()]
    if "gene" not in frame.columns:
        raise IntegrationError(f"DEG table has no gene column: {deg_path}")

    if "significant" in frame.columns:
        flag = (
            frame["significant"]
            .astype(str)
            .str.strip()
            .str.upper()
            .isin(["TRUE", "1", "YES"])
        )
        frame = frame[flag]
    if "direction" in frame.columns and not keep_all:
        frame = frame[frame["direction"].astype(str).str.strip().isin(["Up", "Down"])]

    if "avg_log2fc" not in frame.columns:
        raise IntegrationError(f"DEG table has no log2FC column: {deg_path}")
    if "p_val_adj" not in frame.columns:
        frame["p_val_adj"] = np.nan

    frame = frame.copy()
    frame["gene"] = frame["gene"].astype(str)
    frame["avg_log2fc"] = pd.to_numeric(frame["avg_log2fc"], errors="coerce")
    frame["p_val_adj"] = pd.to_numeric(frame["p_val_adj"], errors="coerce")
    frame["abs_log2fc"] = frame["avg_log2fc"].abs()
    frame = frame.sort_values(
        ["p_val_adj", "abs_log2fc"],
        ascending=[True, False],
        na_position="last",
    ).reset_index(drop=True)
    frame["deg_rank"] = np.arange(1, len(frame) + 1)

    total_before = len(frame)
    if not keep_all:
        patterns = blacklist_patterns or DEFAULT_GENE_BLACKLIST
        if patterns:
            expr = re.compile("|".join(patterns), re.IGNORECASE)
            frame = frame[~frame["gene"].str.match(expr)]
    frame = frame.reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)

    ml_path = data_dir / "ml_feature_importance.csv"
    if ml_path.exists():
        try:
            ml = pd.read_csv(ml_path, index_col=0)
            ml.index = ml.index.astype(str)
            ml_values = ml.iloc[:, 0].astype(float)
            frame["ml_importance"] = frame["gene"].map(ml_values).fillna(0.0)
        except Exception:
            frame["ml_importance"] = 0.0
    else:
        frame["ml_importance"] = 0.0

    out_cols = [
        "rank",
        "gene",
        "direction",
        "avg_log2fc",
        "p_val_adj",
        "pct.1",
        "pct.2",
        "deg_rank",
        "ml_importance",
    ]
    for col in out_cols:
        if col not in frame.columns:
            frame[col] = ""
    frame = frame.head(top_n)[out_cols].reset_index(drop=True)
    frame["source"] = "DEG"

    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "key_genes.csv"
    frame.to_csv(csv_path, index=False)
    summary = {
        "deg_table": str(deg_path),
        "deg_total": int(total_before),
        "after_blacklist": int(len(frame)),
        "top_n": int(top_n),
        "keep_all": bool(keep_all),
        "output_csv": str(csv_path),
    }
    write_json(out_dir / "key_genes_summary.json", summary)
    log.info(
        "key targets: %s genes kept from %s DEGs -> %s",
        len(frame),
        total_before,
        csv_path,
    )
    return frame


def find_rscript() -> str:
    return orchestrator.find_rscript()


def export_pseudobulk(single_cell_root: Path, out_dir: Path) -> dict:
    """Aggregate single-cell counts by sample with the bundled R helper."""
    rscript = find_rscript()
    out_dir.mkdir(parents=True, exist_ok=True)
    script = APP_ROOT / "src" / "pipeline" / "export_pseudobulk.R"
    if not script.exists():
        raise IntegrationError(f"pseudobulk export script missing: {script}")
    log.info("exporting pseudobulk matrix (Rscript: %s)", rscript)
    proc = subprocess.run(
        [rscript, str(script), str(single_cell_root), str(out_dir)],
        cwd=str(APP_ROOT),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=1800,
    )
    if proc.returncode != 0:
        raise IntegrationError(
            "pseudobulk export failed:\n" + (proc.stderr or proc.stdout)[-3000:]
        )
    required = ["pseudobulk_expression.csv", "pseudobulk_metadata.csv"]
    missing = [name for name in required if not (out_dir / name).exists()]
    if missing:
        raise IntegrationError(f"pseudobulk export missing files: {missing}")
    return {
        "expression_csv": str(out_dir / "pseudobulk_expression.csv"),
        "metadata_csv": str(out_dir / "pseudobulk_metadata.csv"),
    }


def _http_json(url: str, payload: dict | None = None, timeout: int = 90) -> dict:
    headers = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}
    data = None
    if payload is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(payload).encode("utf-8")
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8", errors="replace"))
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"HTTP request failed: {url}: {last_error}")


def _mygene_info(gene: str, timeout: int) -> dict:
    payload = {
        "q": gene,
        "scopes": "symbol",
        "fields": "entrezgene,uniprot,ensembl.gene",
        "species": "human",
        "size": 5,
    }
    try:
        hits = _http_json(
            "https://mygene.info/v3/query",
            payload,
            timeout,
        )
    except Exception:
        return {}
    for hit in hits or []:
        if str(hit.get("query", "")).upper() != gene:
            continue
        uniprot = ""
        if isinstance(hit.get("uniprot"), dict):
            uniprot = (
                hit["uniprot"].get("Swiss-Prot")
                or hit["uniprot"].get("SWISSPROT")
                or ""
            )
        elif isinstance(hit.get("uniprot"), str):
            uniprot = hit["uniprot"]
        ensembl = ""
        if isinstance(hit.get("ensembl"), dict):
            ensembl = hit["ensembl"].get("gene") or ""
        elif isinstance(hit.get("ensembl"), str):
            ensembl = hit["ensembl"]
        return {
            "entrez": str(hit.get("entrezgene") or ""),
            "uniprot": uniprot,
            "ensembl": ensembl,
        }
    return {}


def _chembl_evidence(uniprot: str, timeout: int) -> tuple[str, int]:
    if not uniprot:
        return "", 0
    chembl_id = ""
    try:
        res = _http_json(
            "https://www.ebi.ac.uk/chembl/api/data/target.json"
            f"?target_components__accession={uniprot}&limit=50",
            timeout=timeout,
        )
        targets = res.get("targets") or []
        for target in targets:
            if str(target.get("organism", "")).lower() == "homo sapiens":
                chembl_id = target.get("target_chembl_id") or ""
                break
        if not chembl_id and targets:
            chembl_id = targets[0].get("target_chembl_id") or ""
    except Exception:
        chembl_id = ""
    if not chembl_id:
        return chembl_id, 0
    try:
        res = _http_json(
            "https://www.ebi.ac.uk/chembl/api/data/activity.json"
            f"?target_chembl_id={chembl_id}&limit=1",
            timeout=timeout,
        )
        meta = res.get("page_meta") or {}
        return chembl_id, int(meta.get("total_count") or 0)
    except Exception:
        return chembl_id, 0


def _rcsb_evidence(uniprot: str, timeout: int) -> tuple[int, list[str]]:
    if not uniprot:
        return 0, []
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_accession"
                ),
                "operator": "exact_match",
                "value": uniprot,
            },
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": 5}},
    }
    try:
        res = _http_json(
            "https://search.rcsb.org/rcsbsearch/v2/query",
            query,
            timeout,
        )
    except Exception:
        return 0, []
    total = int(res.get("total_count") or 0)
    ids = [
        item.get("identifier", "")
        for item in res.get("result_set") or []
        if item.get("identifier")
    ]
    return total, ids[:5]


def _empty_evidence(gene: str) -> dict:
    return {
        "gene": gene,
        "entrez": "",
        "uniprot": "",
        "ensembl": "",
        "chembl_target_id": "",
        "known_ligands": 0,
        "chembl_bioactivities": 0,
        "pdb_structures": 0,
        "pdb_ids": "",
        "off_target_paralogs": 0,
        "safety_concern": 0,
    }


def _evidence_for_gene(gene: str, timeout: int = 90) -> dict:
    info = _mygene_info(gene, timeout)
    uniprot = info.get("uniprot") or ""
    chembl_id, bioactivities = _chembl_evidence(uniprot, timeout)
    pdb_count, pdb_ids = _rcsb_evidence(uniprot, timeout)
    row = _empty_evidence(gene)
    row.update(
        {
            "entrez": info.get("entrez") or "",
            "uniprot": uniprot,
            "ensembl": info.get("ensembl") or "",
            "chembl_target_id": chembl_id,
            "known_ligands": bioactivities,
            "chembl_bioactivities": bioactivities,
            "pdb_structures": pdb_count,
            "pdb_ids": ",".join(pdb_ids),
        }
    )
    log.info(
        "evidence %s: ligands=%s chembl=%s pdb=%s %s",
        gene,
        bioactivities,
        chembl_id,
        pdb_count,
        pdb_ids,
    )
    return row


def ensure_gene_evidence(
    genes: list[str],
    workdir: Path,
    fetch: bool = True,
    max_workers: int = 6,
    timeout: int = 90,
) -> pd.DataFrame:
    """Return per-gene evidence, reusing the local cache when possible."""
    out_path = _integration_dir(workdir) / "gene_evidence.csv"
    cache: dict[str, dict] = {}
    if out_path.exists():
        try:
            old = pd.read_csv(out_path)
            for _, row in old.iterrows():
                cache[str(row["gene"])] = row.to_dict()
        except Exception:
            cache = {}

    missing = [gene for gene in genes if gene not in cache]
    if missing and fetch:
        log.info("fetching evidence for %s genes", len(missing))
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [
                pool.submit(_evidence_for_gene, gene, timeout)
                for gene in missing
            ]
            for future in futures:
                try:
                    row = future.result()
                    cache[str(row["gene"])] = row
                except Exception as exc:  # noqa: BLE001
                    log.warning("evidence fetch failed: %s", exc)

    rows = []
    for gene in genes:
        row = cache.get(gene) or _empty_evidence(gene)
        rows.append({key: row.get(key, _empty_evidence(gene)[key]) for key in EVIDENCE_COLUMNS})
    frame = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, index=False)
    log.info(
        "gene evidence ready: %s genes, %s with PDB structures",
        len(frame),
        int((pd.to_numeric(frame["pdb_structures"], errors="coerce").fillna(0) > 0).sum()),
    )
    return frame


def build_knockout_inputs(
    single_cell_root: Path,
    workdir: Path,
    skip_pseudobulk: bool = False,
) -> dict:
    """Write expression/metadata/prognosis/druggability inputs for knockout."""
    ko_dir = workdir / "data" / "knockout"
    ko_dir.mkdir(parents=True, exist_ok=True)
    expression_dst = ko_dir / "expression.csv"
    metadata_dst = ko_dir / "metadata.csv"

    if not _knockout_inputs_ready(expression_dst, metadata_dst):
        pseudo_dir = ko_dir / "_pseudobulk"
        if (
            not (pseudo_dir / "pseudobulk_expression.csv").exists()
            or not (pseudo_dir / "pseudobulk_metadata.csv").exists()
            or not _knockout_inputs_ready(
                pseudo_dir / "pseudobulk_expression.csv",
                pseudo_dir / "pseudobulk_metadata.csv",
            )
        ):
            if skip_pseudobulk:
                raise IntegrationError(
                    "pseudobulk files are missing and --skip-pseudobulk was set; "
                    "run the single-cell pipeline first"
                )
            if pseudo_dir.exists():
                shutil.rmtree(pseudo_dir)
            export_pseudobulk(single_cell_root, pseudo_dir)
        if expression_dst.exists():
            expression_dst.unlink()
        if metadata_dst.exists():
            metadata_dst.unlink()
        shutil.copyfile(
            pseudo_dir / "pseudobulk_expression.csv",
            expression_dst,
        )
        shutil.copyfile(
            pseudo_dir / "pseudobulk_metadata.csv",
            metadata_dst,
        )

    expression = pd.read_csv(expression_dst)
    genes = expression.iloc[:, 0].astype(str).tolist()
    metadata = pd.read_csv(metadata_dst)
    if "sample" not in metadata.columns or "condition" not in metadata.columns:
        raise IntegrationError(
            "pseudobulk metadata must contain 'sample' and 'condition' columns"
        )

    evidence_path = _integration_dir(workdir) / "gene_evidence.csv"
    evidence = (
        pd.read_csv(evidence_path)
        if evidence_path.exists()
        else pd.DataFrame([_empty_evidence(gene) for gene in genes])
    )
    evidence = evidence.drop_duplicates("gene", keep="first")

    druggability = evidence[
        [
            "gene",
            "known_ligands",
            "chembl_bioactivities",
            "pdb_structures",
            "off_target_paralogs",
            "safety_concern",
        ]
    ].copy()
    druggability.to_csv(ko_dir / "druggability.csv", index=False)
    off_target = evidence[
        ["gene", "off_target_paralogs", "safety_concern"]
    ].copy()
    off_target.to_csv(ko_dir / "off_target.csv", index=False)

    prognosis = pd.DataFrame({"gene": genes, "hr": 1.0, "p": 1.0})
    prognosis.to_csv(ko_dir / "prognosis.csv", index=False)

    summary = {
        "expression_csv": str(expression_dst),
        "metadata_csv": str(metadata_dst),
        "prognosis_csv": str(ko_dir / "prognosis.csv"),
        "druggability_csv": str(ko_dir / "druggability.csv"),
        "off_target_csv": str(ko_dir / "off_target.csv"),
        "genes": len(genes),
        "samples": int(expression.shape[1] - 1),
        "groups": metadata["condition"].drop_duplicates().tolist(),
    }
    write_json(ko_dir / "inputs_summary.json", summary)
    log.info(
        "knockout inputs ready: %s genes x %s samples -> %s",
        summary["genes"],
        summary["samples"],
        ko_dir,
    )
    return summary


def _knockout_inputs_ready(expression_path: Path, metadata_path: Path) -> bool:
    """Return False when cached knockout inputs are malformed or stale."""
    try:
        if not expression_path.exists() or not metadata_path.exists():
            return False
        expr = pd.read_csv(expression_path)
        meta = pd.read_csv(metadata_path)
        if expr.empty or meta.empty or len(expr.columns) < 2:
            return False
        numeric = expr.drop(columns=[expr.columns[0]]).apply(
            pd.to_numeric, errors="coerce"
        )
        if numeric.dropna(how="all").empty:
            return False
        if not {"sample", "condition"}.issubset(meta.columns):
            return False
        if meta["sample"].isna().any() or meta["sample"].duplicated().any():
            return False
        if meta["condition"].nunique() < 2:
            return False
        return True
    except Exception:
        return False


def run_knockout_stage(
    workdir: Path,
    docking_config: Path,
    inputs: dict,
    case_label: str | None = None,
    normal_label: str | None = None,
    ko_top_n: int | None = None,
    depmap_csv: str | None = None,
) -> dict:
    overrides = {
        "workdir": str(workdir),
        "expression_csv": inputs["expression_csv"],
        "metadata_csv": inputs["metadata_csv"],
        "prognosis_csv": inputs["prognosis_csv"],
        "druggability_csv": inputs["druggability_csv"],
        "off_target_csv": inputs["off_target_csv"],
    }
    if case_label:
        overrides["case_label"] = case_label
    if normal_label:
        overrides["normal_label"] = normal_label
    if ko_top_n:
        overrides["ko_top_n"] = int(ko_top_n)
    if depmap_csv:
        overrides["depmap_csv"] = depmap_csv

    metadata = pd.read_csv(inputs["metadata_csv"])
    if "cell_type" in metadata.columns:
        overrides["cell_type_column"] = "cell_type"

    cfg = load_config(docking_config, overrides)
    ko_summary = run_knockout(cfg, log)
    val_summary = export_validation(cfg, log)
    result = {"knockout": ko_summary, "validation": val_summary}
    write_json(
        _integration_dir(workdir) / "knockout_summary.json",
        result,
    )
    log.info(
        "knockout + validation complete: %s genes scored, %s candidates",
        ko_summary.get("genes_scored", 0),
        val_summary.get("candidates", 0),
    )
    return result


def _split_pdb_ids(value) -> list[str]:
    text = str(value or "")
    return [
        part.strip().upper()
        for part in re.split(r"[,;|\s]+", text)
        if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", part.strip())
    ]


def _download_pdb(pdb_id: str, target_dir: Path, timeout: int = 90) -> Path | None:
    out = target_dir / "data" / "receptors" / f"{pdb_id}.pdb"
    if out.exists() and out.stat().st_size > 0:
        return out
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8", errors="replace")
        if "\nATOM" not in text and not text.startswith("ATOM"):
            log.warning("PDB %s has no ATOM records", pdb_id)
            return None
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
        return out
    except Exception as exc:  # noqa: BLE001
        log.warning("PDB download failed for %s: %s", pdb_id, exc)
        return None


def _prepare_ligand_library(
    target_dir: Path,
    known_ligands_csv: Path,
    fallback_library: str | None,
    pdb_path: Path | None = None,
) -> Path | None:
    if known_ligands_csv.exists():
        try:
            df = pd.read_csv(known_ligands_csv)
            smi_col = next(
                (c for c in df.columns if c.lower() in {"smiles", "canonical_smiles"}),
                None,
            )
            if smi_col is not None:
                df = df[df[smi_col].notna()]
                df["smiles"] = df[smi_col].astype(str).str.strip()
                df = df[df["smiles"] != ""].drop_duplicates("smiles")
                df = df.head(50)
                if not df.empty:
                    id_col = "ligand_id" if "ligand_id" in df.columns else None
                    df["ID"] = [
                        safe_name(str(row.get(id_col, "")), f"ligand_{i + 1}")
                        if id_col
                        else f"ligand_{i + 1}"
                        for i, (_, row) in enumerate(df.iterrows())
                    ]
                    lib = target_dir / "data" / "ligands" / "library.csv"
                    lib.parent.mkdir(parents=True, exist_ok=True)
                    df[["ID", "smiles"]].rename(columns={"smiles": "SMILES"}).to_csv(
                        lib,
                        index=False,
                    )
                    log.info("using %s known ligands for %s", len(df), target_dir.name)
                    return lib
        except Exception as exc:  # noqa: BLE001
            log.warning("known-ligand CSV unusable: %s", exc)

    candidates = []
    if fallback_library:
        candidates.append(Path(fallback_library))
    for name in ["library.smi", "library.sdf", "library.csv"]:
        candidates.append(target_dir.parent.parent / "data" / "ligands" / name)
    for cand in candidates:
        if cand.exists():
            lib = target_dir / "data" / "ligands" / cand.name
            lib.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(cand, lib)
            log.info("using fallback ligand library: %s", cand)
            return lib
    if pdb_path is not None:
        ligands = _extract_cocrystal_ligands(pdb_path)
        if ligands:
            lib = target_dir / "data" / "ligands" / "cocrystal_library.csv"
            lib.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(ligands).rename(columns={"smiles": "SMILES"}).to_csv(
                lib,
                index=False,
            )
            log.info("using %s cocrystal ligands from %s", len(ligands), pdb_path.name)
            return lib
    return None


def _extract_cocrystal_ligands(pdb_path: Path, max_ligands: int = 5) -> list[dict]:
    """Extract non-water HETATM residues as SMILES when DB ligands are missing."""
    try:
        from rdkit import Chem
    except ImportError:
        return []
    try:
        lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    small_ions = {
        "HOH", "WAT", "DOD", "CL", "NA", "K", "MG", "CA", "ZN",
        "SO4", "PO4", "GOL", "EDO", "ACT", "FMT", "DMS", "PEG",
        "IOD", "BR", "CO", "CU", "FE", "MN", "NI",
    }
    het_atoms: list[str] = []
    conect: list[str] = []
    for line in lines:
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname in small_ions:
                continue
            het_atoms.append(line)
        elif line.startswith("CONECT"):
            conect.append(line)
    if len(het_atoms) < 3:
        return []

    groups: dict[tuple[str, str, str], list[str]] = {}
    for line in het_atoms:
        key = (line[21], line[22:26].strip(), line[17:20].strip())
        groups.setdefault(key, []).append(line)

    ligands: list[dict] = []
    for (chain, resseq, resname), atom_lines in groups.items():
        atom_ids = set()
        for line in atom_lines:
            try:
                atom_ids.add(int(line[6:11]))
            except ValueError:
                continue
        block_lines = list(atom_lines)
        for line in conect:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                ids = [int(part) for part in parts[1:] if part.isdigit()]
            except ValueError:
                continue
            if any(atom_id in atom_ids for atom_id in ids):
                block_lines.append(line)
        block = (
            "REMARK generated cocrystal ligand\n"
            + "\n".join(block_lines)
            + "\nEND\n"
        )
        try:
            mol = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=True)
            if mol is None:
                continue
            if mol.GetNumHeavyAtoms() < 6:
                continue
            smiles = Chem.MolToSmiles(mol)
            if not smiles:
                continue
            ligands.append(
                {
                    "id": f"{resname}_{chain}{resseq}",
                    "smiles": smiles,
                }
            )
            if len(ligands) >= max_ligands:
                break
        except Exception:
            continue
    return ligands


def run_target_docking(
    gene: str,
    workdir: Path,
    docking_config: Path,
    evidence: pd.DataFrame,
    ligand_library: str | None,
    force: bool = False,
) -> dict:
    row = evidence[evidence["gene"].astype(str) == gene]
    base = {
        "gene": gene,
        "status": "skipped",
        "pdb_id": "",
        "uniprot": "",
        "ligand_count": 0,
        "hits": 0,
        "best_affinity": "",
        "output_dir": "",
        "error": "",
    }
    if row.empty:
        base["error"] = "no evidence row"
        return base
    info = row.iloc[0].to_dict()
    uniprot_value = info.get("uniprot")
    uniprot = "" if uniprot_value is None or pd.isna(uniprot_value) else str(uniprot_value)
    pdb_ids = _split_pdb_ids(
        "" if info.get("pdb_ids") is None or pd.isna(info.get("pdb_ids"))
        else info.get("pdb_ids")
    )
    if not pdb_ids:
        base["error"] = "no PDB structure"
        return base

    target_dir = workdir / "work" / safe_name(gene, gene)
    target_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = None
    for pdb_id in pdb_ids:
        pdb_path = _download_pdb(pdb_id, target_dir)
        if pdb_path is not None:
            break
    if pdb_path is None:
        base["error"] = "PDB download failed for all candidates"
        return base

    cfg = load_config(
        docking_config,
        {
            "workdir": str(target_dir),
            "target_name": gene,
            "uniprot": uniprot,
            "pdb": pdb_path.stem,
        },
    )
    chembl_id = info.get("chembl_target_id")
    if chembl_id and not pd.isna(chembl_id):
        cfg.data.setdefault("evidence", {})["chembl_target_id"] = str(chembl_id)
    try:
        evidence_mod.gather_evidence(cfg, log)
    except Exception as exc:  # noqa: BLE001
        log.warning("evidence collection failed for %s: %s", gene, exc)

    known = target_dir / "evidence" / "known_ligands.csv"
    library = _prepare_ligand_library(
        target_dir,
        known,
        ligand_library,
        pdb_path=pdb_path,
    )
    if library is None:
        base["error"] = "no ligand library available"
        return base

    center, size, mode = box.detect_box_data(pdb_path)
    cfg.data["receptor"]["input"] = str(pdb_path)
    cfg.data["receptor"]["output"] = str(
        target_dir / "data" / "receptors" / f"{pdb_path.stem}.pdbqt"
    )
    cfg.data["receptor"]["center"] = center
    cfg.data["receptor"]["size"] = size
    cfg.data["receptor"]["detect_input"] = None
    cfg.data["ligand"]["input"] = str(library)
    save_config(cfg, target_dir / "config" / "docking_config.json")
    log.info(
        "docking target %s: PDB %s, box mode %s, center %s size %s",
        gene,
        pdb_path.stem,
        mode,
        center,
        size,
    )

    try:
        docking_pipeline.run_pipeline(cfg, force=force)
    except (DockingError, ToolNotFoundError) as exc:
        base["status"] = "failed"
        base["error"] = str(exc)
        base["pdb_id"] = pdb_path.stem
        base["uniprot"] = uniprot
        return base

    report_dir = cfg.output_dir / "reports"
    summary = _read_json(report_dir / "summary.json")
    ranked = report_dir / "ranked_results.csv"
    hits = 0
    best = ""
    if ranked.exists():
        try:
            ranked_df = pd.read_csv(ranked)
            hits = int((ranked_df["affinity"] <= float(cfg.get("analysis", "cutoff", -7.0))).sum())
            if "affinity" in ranked_df.columns:
                best = str(ranked_df["affinity"].min())
        except Exception:
            hits = int(summary.get("hits", 0))
            best = str(summary.get("best_affinity", ""))
    else:
        hits = int(summary.get("hits", 0))
        best = str(summary.get("best_affinity", ""))

    result = {
        "gene": gene,
        "status": "ok",
        "pdb_id": pdb_path.stem,
        "uniprot": uniprot,
        "ligand_count": int(summary.get("total_docked", 0)),
        "hits": hits,
        "best_affinity": best,
        "output_dir": str(cfg.output_dir),
        "error": "",
    }
    write_json(
        target_dir / "outputs" / "integration" / "target_summary.json",
        result,
    )
    log.info(
        "docking target %s complete: %s ligands, %s hits, best %s",
        gene,
        result["ligand_count"],
        hits,
        best,
    )
    return result


def run_docking_stage(
    workdir: Path,
    docking_config: Path,
    key_genes_csv: Path,
    evidence_csv: Path,
    max_targets: int,
    ligand_library: str | None,
    force: bool = False,
) -> dict:
    out_dir = _integration_dir(workdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = (
        pd.read_csv(key_genes_csv)["gene"]
        .astype(str)
        .head(max_targets)
        .tolist()
    )
    if not genes:
        summary = {"status": "skipped", "reason": "no key genes"}
        write_json(out_dir / "docking_summary.json", summary)
        return summary
    evidence = pd.read_csv(evidence_csv)
    rows = []
    for gene in genes:
        try:
            rows.append(
                run_target_docking(
                    gene,
                    workdir,
                    docking_config,
                    evidence,
                    ligand_library,
                    force=force,
                )
            )
        except Exception as exc:  # noqa: BLE001
            log.error("docking target %s crashed: %s", gene, exc)
            rows.append(
                {
                    "gene": gene,
                    "status": "failed",
                    "pdb_id": "",
                    "uniprot": "",
                    "ligand_count": 0,
                    "hits": 0,
                    "best_affinity": "",
                    "output_dir": "",
                    "error": str(exc),
                }
            )
    frame = pd.DataFrame(rows)
    frame.to_csv(out_dir / "docking_targets.csv", index=False)
    summary = {
        "status": "completed",
        "targets_requested": len(genes),
        "ok": int((frame["status"] == "ok").sum()),
        "failed": int((frame["status"] == "failed").sum()),
        "skipped": int((frame["status"] == "skipped").sum()),
        "total_hits": int(pd.to_numeric(frame["hits"], errors="coerce").fillna(0).sum()),
        "best_affinity": (
            str(
                frame.loc[
                    pd.to_numeric(frame["hits"], errors="coerce").fillna(0) > 0,
                    "best_affinity",
                ].min()
            )
            if len(frame)
            else ""
        ),
        "output_csv": str(out_dir / "docking_targets.csv"),
    }
    write_json(out_dir / "docking_summary.json", summary)
    log.info(
        "docking stage complete: %s ok / %s failed / %s skipped",
        summary["ok"],
        summary["failed"],
        summary["skipped"],
    )
    return summary


def _esc(value) -> str:
    import html

    return html.escape(str(value if value is not None else ""))


def _render_table(frame: pd.DataFrame, columns: list[str]) -> str:
    if frame is None or frame.empty:
        return '<p class="muted">No data.</p>'
    head = "".join(f"<th>{_esc(c)}</th>" for c in columns)
    body = ""
    for _, row in frame.head(20).iterrows():
        cells = "".join(f"<td>{_esc(row.get(c, ''))}</td>" for c in columns)
        body += f"<tr>{cells}</tr>"
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def generate_integrated_report(
    workdir: Path,
    single_cell_root: Path,
    docking_config: Path,
    ctx: dict,
) -> Path:
    out_dir = _integration_dir(workdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sc_summary = _read_json(single_cell_root / "results" / "summary.json")
    key_genes = pd.read_csv(out_dir / "key_genes.csv") if (out_dir / "key_genes.csv").exists() else pd.DataFrame()
    ko_summary = _read_json(out_dir / "knockout_summary.json")
    ko_top = pd.DataFrame()
    ko_ranked = (
        workdir / "outputs" / "run_001" / "knockout" / "ranked_knockout.csv"
    )
    if not ko_ranked.exists():
        ko_ranked = (
            workdir / "outputs" / "run_001" / "knockout" / "ranked_knockout.csv"
        )
    if ko_ranked.exists():
        ko_top = pd.read_csv(ko_ranked)
    docking_summary = _read_json(out_dir / "docking_summary.json")
    docking = pd.read_csv(out_dir / "docking_targets.csv") if (out_dir / "docking_targets.csv").exists() else pd.DataFrame()
    evidence = pd.read_csv(out_dir / "gene_evidence.csv") if (out_dir / "gene_evidence.csv").exists() else pd.DataFrame()

    sc_html = _render_table(
        pd.DataFrame(
            [
                {
                    "accession": sc_summary.get("dataset", ""),
                    "cells": sc_summary.get("n_cells_after_doublet_removal", ""),
                    "genes": sc_summary.get("n_genes", ""),
                    "deg_up": sc_summary.get("deg_up", ""),
                    "deg_down": sc_summary.get("deg_down", ""),
                }
            ]
        ),
        ["accession", "cells", "genes", "deg_up", "deg_down"],
    )
    ko_cols = [
        c
        for c in [
            "rank",
            "gene",
            "target_class",
            "target_score",
            "knockout_score",
            "druggability_score",
        ]
        if c in ko_top.columns
    ]
    dock_cols = [
        c
        for c in [
            "gene",
            "status",
            "pdb_id",
            "ligand_count",
            "hits",
            "best_affinity",
        ]
        if c in docking.columns
    ]
    ev_cols = [
        c
        for c in [
            "gene",
            "uniprot",
            "known_ligands",
            "pdb_structures",
            "pdb_ids",
        ]
        if c in evidence.columns
    ]

    def rel(path):
        try:
            return os.path.relpath(path, out_dir)
        except ValueError:
            return str(path)

    html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Integrated Discovery Pipeline Report</title>
<style>
body {{ font-family: "Segoe UI", Arial, sans-serif; margin: 24px; color: #1f2933; background: #f5f7fa; }}
h1 {{ font-size: 24px; }}
h2 {{ font-size: 18px; margin-top: 22px; }}
.card {{ background: #fff; border: 1px solid #e4e7eb; border-radius: 8px; padding: 16px; margin-bottom: 14px; }}
table {{ width: 100%; border-collapse: collapse; font-size: 12px; }}
th, td {{ border: 1px solid #e5e7eb; padding: 6px 7px; text-align: left; }}
th {{ background: #eef2f7; }}
.muted {{ color: #6b7280; }}
a {{ color: #1d4ed8; }}
</style>
</head>
<body>
<h1>Integrated Discovery Pipeline Report</h1>
<div class="card">
  <p><b>Single-cell output:</b> {rel(single_cell_root)}</p>
  <p><b>Integration output:</b> {rel(out_dir)}</p>
  <p><b>Docking summary:</b> {_esc(docking_summary)}</p>
</div>
<div class="card">
  <h2>Single-cell summary</h2>
  {sc_html}
</div>
<div class="card">
  <h2>Key genes (top 20)</h2>
  {_render_table(key_genes, ["rank", "gene", "direction", "avg_log2fc", "p_val_adj"])}
</div>
<div class="card">
  <h2>Virtual knockout targets (top 20)</h2>
  {_render_table(ko_top, ko_cols)}
</div>
<div class="card">
  <h2>Docking per target</h2>
  {_render_table(docking, dock_cols)}
</div>
<div class="card">
  <h2>Gene evidence</h2>
  {_render_table(evidence, ev_cols)}
</div>
<div class="card">
  <h2>Outputs</h2>
  <ul>
    <li><a href="{rel(out_dir / 'key_genes.csv')}">key_genes.csv</a></li>
    <li><a href="{rel(ko_ranked) if ko_ranked.exists() else '#'}">ranked_knockout.csv</a></li>
    <li><a href="{rel(out_dir / 'docking_targets.csv') if (out_dir / 'docking_targets.csv').exists() else '#'}">docking_targets.csv</a></li>
  </ul>
</div>
</body>
</html>
"""
    report_path = out_dir / "integration_report.html"
    report_path.write_text(html_text, encoding="utf-8")

    summary = {
        "single_cell": sc_summary,
        "key_genes": len(key_genes),
        "knockout": {
            "genes_scored": (ko_summary.get("knockout") or {}).get("genes_scored", 0),
            "validation_candidates": (ko_summary.get("validation") or {}).get("candidates", 0),
        },
        "docking": docking_summary,
        "evidence_genes": len(evidence),
        "report_html": str(report_path),
        "finished_at": datetime.now().isoformat(timespec="seconds"),
    }
    write_json(out_dir / "integration_summary.json", summary)
    cfg = load_config(docking_config, {"workdir": str(workdir)})
    write_run_manifest(
        out_dir,
        cfg,
        "full-pipeline",
        {
            "key_genes_csv": out_dir / "key_genes.csv",
            "gene_evidence_csv": out_dir / "gene_evidence.csv",
            "integration_summary_json": out_dir / "integration_summary.json",
        },
        summary,
    )
    log.info("integrated report generated: %s", report_path)
    return report_path


def _stage_single_cell(args, workdir: Path, ctx: dict) -> None:
    root = ctx["single_cell_root"]
    if args.skip_scrna:
        if not (root / "results" / "pipeline_complete.json").exists():
            raise IntegrationError(
                f"single-cell outputs not found under {root}; remove --skip-scrna"
            )
        log.info("using existing single-cell outputs: %s", root)
        return
    code = orchestrator.run_pipeline(
        args.force,
        args.skip_download,
        args.skip_deps,
        args.accession,
        str(root),
        args.species,
    )
    if code == 98:
        raise PauseRequested("single-cell pipeline paused; run again to resume")
    if code != 0:
        raise IntegrationError(f"single-cell pipeline exited with code {code}")


def _stage_key_targets(args, workdir: Path, ctx: dict) -> None:
    frame = extract_key_genes(
        ctx["single_cell_root"],
        _integration_dir(workdir),
        top_n=args.top_genes,
        keep_all=args.keep_all_genes,
    )
    ctx["key_genes_path"] = _integration_dir(workdir) / "key_genes.csv"
    ctx["key_genes"] = frame


def _stage_evidence(args, workdir: Path, ctx: dict) -> None:
    genes = pd.read_csv(ctx["key_genes_path"])["gene"].astype(str).tolist()
    frame = ensure_gene_evidence(
        genes,
        workdir,
        fetch=not args.skip_evidence_fetch,
        max_workers=args.evidence_workers,
        timeout=args.evidence_timeout,
    )
    ctx["evidence"] = frame
    ctx["evidence_path"] = _integration_dir(workdir) / "gene_evidence.csv"


def _stage_knockout_inputs(args, workdir: Path, ctx: dict) -> None:
    ctx["knockout_inputs"] = build_knockout_inputs(
        ctx["single_cell_root"],
        workdir,
        skip_pseudobulk=args.skip_pseudobulk,
    )


def _stage_knockout(args, workdir: Path, ctx: dict) -> None:
    if args.skip_knockout:
        summary = {"status": "skipped", "reason": "knockout disabled by arguments"}
        write_json(_integration_dir(workdir) / "knockout_summary.json", summary)
        ctx["knockout"] = summary
        return
    ctx["knockout"] = run_knockout_stage(
        workdir,
        ctx["docking_config"],
        ctx["knockout_inputs"],
        case_label=args.case_label,
        normal_label=args.normal_label,
        ko_top_n=args.ko_top_n,
        depmap_csv=args.depmap_csv,
    )


def _stage_docking(args, workdir: Path, ctx: dict) -> None:
    if args.skip_docking or args.docking_targets <= 0:
        summary = {
            "status": "skipped",
            "reason": "docking disabled by arguments",
            "ok": 0,
            "failed": 0,
            "skipped": 0,
        }
        write_json(_integration_dir(workdir) / "docking_summary.json", summary)
        ctx["docking"] = summary
        return
    if "key_genes_path" not in ctx:
        key_genes_path = _integration_dir(workdir) / "key_genes.csv"
        if not key_genes_path.exists():
            raise IntegrationError(
                "key_genes.csv missing; run stage 02 before docking"
            )
        ctx["key_genes_path"] = key_genes_path
    if "evidence_path" not in ctx:
        evidence_path = _integration_dir(workdir) / "gene_evidence.csv"
        if not evidence_path.exists():
            raise IntegrationError(
                "gene_evidence.csv missing; run stage 03 before docking"
            )
        ctx["evidence_path"] = evidence_path
    ctx["docking"] = run_docking_stage(
        workdir,
        ctx["docking_config"],
        ctx["key_genes_path"],
        ctx["evidence_path"],
        max_targets=args.docking_targets,
        ligand_library=args.ligand_library,
        force=args.force,
    )


def _stage_report(args, workdir: Path, ctx: dict) -> None:
    ctx["report"] = generate_integrated_report(
        workdir,
        ctx["single_cell_root"],
        ctx["docking_config"],
        ctx,
    )


def run_full_pipeline(args) -> int:
    """Run the integrated pipeline with per-stage resume markers."""
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    cfg_path = Path(args.docking_config).resolve()
    ctx = {
        "single_cell_root": Path(args.output).resolve(),
        "workdir": workdir,
        "docking_config": cfg_path,
    }

    for code, name, _description in STAGES:
        marker = _marker(workdir, code, name)
        if args.start_stage and code < args.start_stage:
            marker.parent.mkdir(parents=True, exist_ok=True)
            marker.write_text("skipped by start-stage", encoding="utf-8")
            log.info("stage %s %s skipped by --start-stage", code, name)
            continue
        if not args.force and marker.exists():
            if code == "01" and not _single_cell_outputs_ready(
                ctx["single_cell_root"]
            ):
                log.warning(
                    "stage 01 marker exists but single-cell outputs are missing "
                    "under %s; rerunning",
                    ctx["single_cell_root"],
                )
                marker.unlink(missing_ok=True)
            elif code == "04" and not _knockout_inputs_ready(
                workdir / "data" / "knockout" / "expression.csv",
                workdir / "data" / "knockout" / "metadata.csv",
            ):
                log.warning(
                    "stage 04 marker exists but knockout inputs are missing or "
                    "malformed; rebuilding pseudobulk"
                )
                marker.unlink(missing_ok=True)
            else:
                log.info("skip stage %s %s (already done)", code, name)
                continue
        log.info("=== stage %s %s ===", code, name)
        fn = {
            "01": _stage_single_cell,
            "02": _stage_key_targets,
            "03": _stage_evidence,
            "04": _stage_knockout_inputs,
            "05": _stage_knockout,
            "06": _stage_docking,
            "07": _stage_report,
        }[code]
        try:
            fn(args, workdir, ctx)
        except PauseRequested as exc:
            log.info(str(exc))
            return 98
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(
            datetime.now().isoformat(timespec="seconds"),
            encoding="utf-8",
        )
        log.info("stage %s %s complete", code, name)
    log.info("full pipeline complete: %s", _integration_dir(workdir))
    return 0


def load_full_config(path: Path) -> dict:
    defaults = {
        "accession": "GSE125449",
        "single_cell_output": "../liver_cancer",
        "workdir": "dock",
        "species": "auto",
        "top_genes": 50,
        "docking_targets": 3,
        "keep_all_genes": False,
        "case_label": None,
        "normal_label": None,
        "ligand_library": None,
        "ko_top_n": None,
        "depmap_csv": None,
        "evidence": {"fetch": True, "max_workers": 6, "timeout": 90},
        "gene_blacklist": DEFAULT_GENE_BLACKLIST,
    }
    raw = {}
    if path.exists():
        raw = json.loads(path.read_text(encoding="utf-8"))
    config = dict(defaults)
    config.update(raw or {})
    config["evidence"] = dict(defaults["evidence"])
    config["evidence"].update((raw.get("evidence") or {}))
    config["gene_blacklist"] = (
        raw.get("gene_blacklist") or defaults["gene_blacklist"]
    )
    return config


def _apply_defaults(args, config: dict) -> None:
    if args.accession is None:
        args.accession = config.get("accession", "GSE125449")
    if args.output is None:
        args.output = config.get("single_cell_output", "../liver_cancer")
    if args.workdir is None:
        args.workdir = config.get("workdir", "dock")
    if args.species is None:
        args.species = config.get("species", "auto")
    if args.top_genes is None:
        args.top_genes = int(config.get("top_genes", 50))
    if args.docking_targets is None:
        args.docking_targets = int(config.get("docking_targets", 3))
    if args.ligand_library is None:
        args.ligand_library = config.get("ligand_library")
    if args.case_label is None:
        args.case_label = config.get("case_label")
    if args.normal_label is None:
        args.normal_label = config.get("normal_label")
    if args.ko_top_n is None:
        args.ko_top_n = config.get("ko_top_n")
    if args.depmap_csv is None:
        args.depmap_csv = config.get("depmap_csv")
    if args.keep_all_genes is None:
        args.keep_all_genes = bool(config.get("keep_all_genes", False))
    if args.skip_evidence_fetch is None:
        args.skip_evidence_fetch = not bool(
            config.get("evidence", {}).get("fetch", True)
        )
    if args.evidence_workers is None:
        args.evidence_workers = int(config.get("evidence", {}).get("max_workers", 6))
    if args.evidence_timeout is None:
        args.evidence_timeout = int(config.get("evidence", {}).get("timeout", 90))
    if args.workdir:
        args.workdir = str(_resolve_path(args.workdir, APP_ROOT))
    if args.output:
        args.output = str(_resolve_path(args.output, APP_ROOT))
    if args.ligand_library:
        args.ligand_library = str(_resolve_path(args.ligand_library, Path.cwd()))
    if args.depmap_csv:
        args.depmap_csv = str(_resolve_path(args.depmap_csv, Path.cwd()))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_full_pipeline",
        description="Automated scRNA-seq -> key targets -> virtual screening -> knockout",
    )
    parser.add_argument(
        "--config",
        default=str(APP_ROOT / "config" / "full_pipeline_config.json"),
        help="full pipeline config JSON",
    )
    parser.add_argument("--accession", help="GEO accession (default from config)")
    parser.add_argument("--output", help="single-cell output root")
    parser.add_argument("--workdir", help="docking/integration work root")
    parser.add_argument(
        "--docking-config",
        default=str(APP_ROOT / "config" / "docking_config.json"),
        help="base docking config JSON",
    )
    parser.add_argument("--species", choices=["hs", "mm", "auto"])
    parser.add_argument("--top-genes", type=int, help="number of key genes to keep")
    parser.add_argument("--docking-targets", type=int, help="max genes to dock")
    parser.add_argument("--ligand-library", help="ligand library file (.smi/.sdf/.csv)")
    parser.add_argument("--case-label", help="case group label for knockout")
    parser.add_argument("--normal-label", help="normal group label for knockout")
    parser.add_argument("--ko-top-n", type=int, help="top N knockout report genes")
    parser.add_argument("--depmap-csv", help="DepMap gene effect CSV")
    parser.add_argument("--evidence-workers", type=int, default=None)
    parser.add_argument("--evidence-timeout", type=int, default=None)
    parser.add_argument("--skip-scrna", action="store_true")
    parser.add_argument("--skip-download", action="store_true")
    parser.add_argument("--skip-deps", action="store_true")
    parser.add_argument("--skip-evidence-fetch", action="store_true", default=None)
    parser.add_argument("--skip-pseudobulk", action="store_true")
    parser.add_argument("--skip-knockout", action="store_true")
    parser.add_argument("--skip-docking", action="store_true")
    parser.add_argument("--keep-all-genes", action="store_true", default=None)
    parser.add_argument("--force", action="store_true", help="rerun stages from scratch")
    parser.add_argument("--start-stage", default=None, help="stage code to start from")
    parser.add_argument("--list-stages", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.list_stages:
        for code, name, description in STAGES:
            print(f"{code}  {name:<20} {description}")
        return 0

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )
    config = load_full_config(Path(args.config))
    _apply_defaults(args, config)

    try:
        return run_full_pipeline(args)
    except (IntegrationError, DockingError, ToolNotFoundError) as exc:
        log.error("ERROR: %s", exc)
        return 1
    except KeyboardInterrupt:
        log.info("interrupted")
        return 130


if __name__ == "__main__":
    sys.exit(main())
