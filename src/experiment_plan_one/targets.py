"""Compound characterization and compound/disease target collection.

Network-dependent sources are queried opportunistically. Every source keeps
an explicit status record, and a failed source is never replaced with
fabricated target associations.
"""

from __future__ import annotations

import html
import http.cookiejar
import json
import logging
import math
import re
import time
from concurrent.futures import ThreadPoolExecutor
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup
from PIL import Image

from .common import (
    LOG,
    ensure_dir,
    get_json,
    read_json,
    save_figure,
    split_gene_symbol,
    write_json,
)

PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
    "MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,"
    "IUPACName,InChIKey/JSON"
)
PUBCHEM_3D_SDF_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
    "?record_type=3d"
)
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
SWISS_TARGET_BASE = "https://www.swisstargetprediction.ch"
OPEN_TARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"


def compound_properties(cid: int | str) -> dict[str, Any]:
    """Fetch authoritative 2D identifiers from PubChem."""
    payload = get_json(PUBCHEM_PROPERTY_URL.format(cid=cid))
    rows = ((payload.get("PropertyTable") or {}).get("Properties") or [])
    if not rows:
        raise RuntimeError(f"PubChem returned no compound properties for CID {cid}")
    row = rows[0]
    return {
        "pubchem_cid": str(row.get("CID") or cid),
        "molecular_formula": row.get("MolecularFormula") or "",
        "molecular_weight": _number(row.get("MolecularWeight")),
        "canonical_smiles": row.get("SMILES") or row.get("CanonicalSMILES") or "",
        "isomeric_smiles": row.get("IsomericSMILES") or "",
        "iupac_name": row.get("IUPACName") or "",
        "inchi_key": row.get("InChIKey") or "",
    }


def _number(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def rdkit_descriptors(smiles: str) -> dict[str, Any]:
    """Calculate transparent physicochemical descriptors with RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"RDKit is unavailable: {exc}") from exc
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    return {
        "rdkit_logp": float(Crippen.MolLogP(molecule)),
        "h_bond_donors": int(Lipinski.NumHDonors(molecule)),
        "h_bond_acceptors": int(Lipinski.NumHAcceptors(molecule)),
        "tpsa_angstrom2": float(rdMolDescriptors.CalcTPSA(molecule)),
        "rotatable_bonds": int(Lipinski.NumRotatableBonds(molecule)),
        "heavy_atoms": int(molecule.GetNumHeavyAtoms()),
        "rings": int(rdMolDescriptors.CalcNumRings(molecule)),
        "fraction_csp3": float(rdMolDescriptors.CalcFractionCSP3(molecule)),
        "qed": float(QED.qed(molecule)),
        "mol_wt_rdkit": float(Descriptors.MolWt(molecule)),
    }


def write_compound_figures(
    smiles: str,
    properties: dict[str, Any],
    descriptors: dict[str, Any],
    out_dir: Path,
) -> dict[str, Path]:
    """Render 2D structure, 3D conformer and a compact property panel."""
    from matplotlib import pyplot as plt
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw

    ensure_dir(out_dir)
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")

    path_2d = out_dir / "fig1b_compound_2d.png"
    drawer = Draw.MolDraw2DCairo(2800, 1900)
    drawer.drawOptions().useBWAtomPalette()
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
    drawer.FinishDrawing()
    path_2d.write_bytes(drawer.GetDrawingText())
    with Image.open(path_2d) as image:
        image.save(path_2d, format="PNG", dpi=(600, 600))
    svg_drawer = Draw.MolDraw2DSVG(2800, 1900)
    svg_drawer.drawOptions().useBWAtomPalette()
    Draw.rdMolDraw2D.PrepareAndDrawMolecule(svg_drawer, molecule)
    svg_drawer.FinishDrawing()
    path_2d.with_suffix(".svg").write_text(
        svg_drawer.GetDrawingText(),
        encoding="utf-8",
    )

    path_3d = out_dir / "fig1c_compound_3d.png"
    molecule_3d = Chem.AddHs(molecule)
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260913
    status = AllChem.EmbedMolecule(molecule_3d, params)
    if status != 0:
        AllChem.EmbedMolecule(molecule_3d, useRandomCoords=True, randomSeed=20260913)
    try:
        AllChem.MMFFOptimizeMolecule(molecule_3d, maxIters=1000)
    except Exception:
        pass
    conformer = molecule_3d.GetConformer()
    coords = np.asarray(
        [
            [conformer.GetAtomPosition(index).x, conformer.GetAtomPosition(index).y, conformer.GetAtomPosition(index).z]
            for index in range(molecule_3d.GetNumAtoms())
        ]
    )
    colors = {
        "C": "#59636e",
        "N": "#2f6bb3",
        "O": "#d1495b",
        "H": "#d8dde3",
        "S": "#d9a21b",
        "F": "#59a14f",
        "Cl": "#59a14f",
        "Br": "#8b5a2b",
        "I": "#7f5aa2",
    }
    atom_colors = [
        colors.get(atom.GetSymbol(), "#777777") for atom in molecule_3d.GetAtoms()
    ]
    sizes = [90 if atom.GetSymbol() != "H" else 30 for atom in molecule_3d.GetAtoms()]
    fig = plt.figure(figsize=(5.8, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    for bond in molecule_3d.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ax.plot(
            coords[[begin, end], 0],
            coords[[begin, end], 1],
            coords[[begin, end], 2],
            color="#a9b0b8",
            linewidth=1.2,
        )
    ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=atom_colors, s=sizes)
    ax.set_axis_off()
    ax.view_init(elev=18, azim=32)
    ax.set_title("6PPD-Q 3D conformer")
    save_figure(fig, path_3d)

    path_props = out_dir / "fig1c_physicochemical_properties.png"
    values = [
        ("Molecular weight", properties.get("molecular_weight"), "g/mol"),
        ("LogP (RDKit)", descriptors.get("rdkit_logp"), ""),
        ("H-bond donors", descriptors.get("h_bond_donors"), ""),
        ("H-bond acceptors", descriptors.get("h_bond_acceptors"), ""),
        ("TPSA", descriptors.get("tpsa_angstrom2"), "A^2"),
        ("Rotatable bonds", descriptors.get("rotatable_bonds"), ""),
        ("Heavy atoms", descriptors.get("heavy_atoms"), ""),
        ("QED", descriptors.get("qed"), ""),
    ]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.axis("off")
    y = 0.94
    for label, value, unit in values:
        shown = "NA" if value is None else f"{float(value):.3g}"
        ax.text(0.02, y, label, ha="left", va="top", fontsize=10)
        ax.text(
            0.98,
            y,
            f"{shown} {unit}".strip(),
            ha="right",
            va="top",
            fontsize=10,
            fontweight="bold",
        )
        y -= 0.12
    ax.text(
        0.02,
        0.02,
        "Descriptors are calculated locally with RDKit; SwissADME values are not "
        "substituted.",
        ha="left",
        va="bottom",
        fontsize=7.5,
        color="#5d6670",
    )
    save_figure(fig, path_props)
    path_combined = out_dir / "fig1c_compound_3d_properties.png"
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.2, 4.4),
        gridspec_kw={"width_ratios": [1.12, 0.88]},
    )
    with Image.open(path_3d) as image:
        axes[0].imshow(image)
    axes[0].axis("off")
    axes[0].set_title("6PPD-Q 3D conformer", fontsize=10, fontweight="bold")
    axes[1].axis("off")
    y = 0.96
    for label, value, unit in values:
        shown = "NA" if value is None else f"{float(value):.3g}"
        axes[1].text(0.0, y, label, ha="left", va="top", fontsize=8.5)
        axes[1].text(
            1.0,
            y,
            f"{shown} {unit}".strip(),
            ha="right",
            va="top",
            fontsize=8.5,
            fontweight="bold",
        )
        y -= 0.115
    axes[1].text(
        0.0,
        0.02,
        "Local RDKit descriptors",
        ha="left",
        va="bottom",
        fontsize=7,
        color="#5d6670",
    )
    fig.tight_layout()
    save_figure(fig, path_combined)
    return {
        "structure_2d": path_2d,
        "structure_3d": path_3d,
        "properties": path_props,
        "structure_3d_properties": path_combined,
    }


def _fetch_swiss_target(
    smiles: str,
    organism: str = "Homo_sapiens",
    *,
    timeout: int = 300,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Submit one molecule and parse the SwissTargetPrediction result table."""
    cookiejar = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(cookiejar))
    headers = {
        "User-Agent": "Mozilla/5.0 (compatible; liver-cancer-pipeline/1.6)",
        "Accept": "text/html,application/xhtml+xml",
    }
    opener.open(f"{SWISS_TARGET_BASE}/index.php", timeout=60).read()
    body = urllib.parse.urlencode(
        {"organism": organism, "smiles": smiles, "ioi": "2"}
    ).encode()
    request = urllib.request.Request(
        f"{SWISS_TARGET_BASE}/predict.php",
        data=body,
        headers={**headers, "Referer": f"{SWISS_TARGET_BASE}/index.php"},
    )
    response = opener.open(request, timeout=timeout)
    submitted_url = response.geturl()
    submitted = response.read().decode("utf-8", "replace")
    match = re.search(r'location\.replace\("([^"]+)', submitted)
    if not match:
        reason = _html_message(submitted) or "SwissTargetPrediction did not return a job URL"
        return pd.DataFrame(), {"status": "failed", "reason": reason}
    result_url = urllib.parse.urljoin(SWISS_TARGET_BASE, match.group(1))
    last_html = ""
    for _ in range(60):
        with opener.open(result_url, timeout=timeout) as result_response:
            last_html = result_response.read().decode("utf-8", "replace")
        if "Target prediction" in last_html and "probability" in last_html.lower():
            break
        if "Data available" in last_html or "resultTable" in last_html:
            break
        time.sleep(2)
    frame = _parse_swiss_target_table(last_html)
    status = {
        "status": "completed" if not frame.empty else "empty",
        "job_url": result_url,
    }
    return frame, status


def _parse_swiss_target_table(page: str) -> pd.DataFrame:
    soup = BeautifulSoup(page, "html.parser")
    records: list[dict[str, Any]] = []
    for row in soup.find_all("tr"):
        cells = row.find_all("td")
        if len(cells) < 6:
            continue
        gene = ""
        for anchor in cells[1].find_all("a", href=True):
            parsed = urllib.parse.urlparse(anchor["href"])
            query = urllib.parse.parse_qs(parsed.query)
            candidate = split_gene_symbol((query.get("gene") or [""])[0])
            if candidate and candidate != "N/A":
                gene = candidate
                break
        if not gene:
            gene = split_gene_symbol(cells[1].get_text(" ", strip=True))
        target_name = _clean_target_name(cells[0].get_text(" ", strip=True))
        uniprot = cells[2].get_text(" ", strip=True)
        probability_text = cells[5].get_text(" ", strip=True)
        probability_match = re.search(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", probability_text)
        probability = float(probability_match.group(1)) if probability_match else np.nan
        if not gene:
            continue
        records.append(
            {
                "target_name": target_name,
                "gene": gene,
                "uniprot": uniprot,
                "probability": probability,
            }
        )
    output = pd.DataFrame(records)
    if output.empty:
        return output
    output = output[output["gene"] != ""].drop_duplicates("gene", keep="first")
    return output.reset_index(drop=True)


def _clean_target_name(value: str) -> str:
    text = re.sub(r"\s+", " ", str(value or "")).strip()
    return text


def _extract_gene_symbol(value: str) -> str:
    text = _clean_target_name(value)
    candidates = re.findall(r"\b[A-Z][A-Z0-9-]{1,14}\b", text)
    stop = {
        "HUMAN",
        "HOMO",
        "SAPIENS",
        "PROTEIN",
        "RECEPTOR",
        "KINASE",
        "CHAIN",
        "ISOFORM",
        "UNIPROT",
    }
    for candidate in candidates:
        if candidate not in stop and not candidate.startswith("P0"):
            return candidate
    return ""


def _html_message(page: str) -> str:
    soup = BeautifulSoup(page, "html.parser")
    for script in soup.find_all("script"):
        text = script.get_text(" ", strip=True)
        match = re.search(r'alert\("([^"]+)', text)
        if match:
            return html.unescape(match.group(1))
    text = soup.get_text(" ", strip=True)
    return text[-500:] if text else ""


def _fetch_stitch(cid: str) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Query the legacy STITCH API for a PubChem CID."""
    identifier = f"CIDs{int(cid):09d}"
    url = (
        "https://stitch-db.org/api/tsv/interactors?"
        + urllib.parse.urlencode(
            {
                "identifiers": identifier,
                "species": "9606",
                "limit": "100",
            }
        )
    )
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            text = response.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", "replace")
        return pd.DataFrame(), {
            "status": "no_results" if exc.code == 400 and "no results" in detail else "failed",
            "http_status": exc.code,
            "reason": detail[-500:],
            "identifier": identifier,
        }
    except Exception as exc:  # noqa: BLE001
        return pd.DataFrame(), {"status": "failed", "reason": str(exc), "identifier": identifier}
    frame = _read_tsv_text(text)
    if frame.empty:
        return frame, {"status": "empty", "identifier": identifier}
    gene_column = next(
        (
            column
            for column in frame.columns
            if str(column).lower() in {"preferredname_b", "preferred_name_b", "gene", "target"}
        ),
        None,
    )
    if gene_column is None:
        return pd.DataFrame(), {
            "status": "failed",
            "reason": "STITCH response has no gene column",
            "columns": list(frame.columns),
        }
    score_column = next(
        (
            column
            for column in frame.columns
            if "score" in str(column).lower() or "confidence" in str(column).lower()
        ),
        None,
    )
    output = pd.DataFrame(
        {
            "gene": frame[gene_column].map(split_gene_symbol),
            "score": pd.to_numeric(frame[score_column], errors="coerce") if score_column else np.nan,
            "source": "STITCH",
        }
    )
    output = output[output["gene"] != ""].drop_duplicates("gene")
    return output.reset_index(drop=True), {"status": "completed", "identifier": identifier}


def _read_tsv_text(text: str) -> pd.DataFrame:
    from io import StringIO

    if not text.strip():
        return pd.DataFrame()
    try:
        return pd.read_csv(StringIO(text), sep="\t")
    except Exception:
        return pd.DataFrame()


def _chembl_similarity_targets(
    smiles: str,
    *,
    threshold: int = 40,
    max_molecules: int = 10,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Use ChEMBL similarity hits as an auditable ligand-based target source."""
    encoded = urllib.parse.quote(smiles, safe="")
    url = f"{CHEMBL_BASE}/similarity/{encoded}/{threshold}.json?limit={max_molecules}"
    try:
        payload = get_json(url, timeout=120)
    except Exception as exc:  # noqa: BLE001
        return pd.DataFrame(), {"status": "failed", "reason": str(exc)}
    molecules = payload.get("molecules") or []
    records: list[dict[str, Any]] = []
    molecule_ids: list[str] = []
    similarity: dict[str, float] = {}
    for molecule in molecules:
        molecule_id = str(molecule.get("molecule_chembl_id") or "")
        if not molecule_id:
            continue
        molecule_ids.append(molecule_id)
        similarity[molecule_id] = float(molecule.get("similarity") or 0.0)
        if len(molecule_ids) >= max_molecules:
            break
    if not molecule_ids:
        return pd.DataFrame(), {"status": "empty", "reason": "no similar ChEMBL molecules"}

    for molecule_id in molecule_ids:
        activity_url = (
            f"{CHEMBL_BASE}/activity.json?molecule_chembl_id={molecule_id}"
            "&limit=200&pchembl_value__isnull=false"
        )
        try:
            activity_payload = get_json(activity_url, timeout=120)
        except Exception as exc:  # noqa: BLE001
            LOG.debug("ChEMBL activities failed for %s: %s", molecule_id, exc)
            continue
        for activity in activity_payload.get("activities") or []:
            target_id = str(activity.get("target_chembl_id") or "")
            pchembl = _number(activity.get("pchembl_value"))
            standard_value = _number(activity.get("standard_value"))
            if not target_id:
                continue
            records.append(
                {
                    "molecule_chembl_id": molecule_id,
                    "tanimoto": similarity.get(molecule_id, np.nan),
                    "target_chembl_id": target_id,
                    "target_pref_name": activity.get("target_pref_name") or "",
                    "standard_type": activity.get("standard_type") or "",
                    "standard_value": standard_value,
                    "standard_units": activity.get("standard_units") or "",
                    "pchembl_value": pchembl,
                }
            )
    if not records:
        return pd.DataFrame(), {
            "status": "empty",
            "reason": "similar ChEMBL molecules had no target bioactivities",
            "n_similar_molecules": len(molecule_ids),
        }
    activities = pd.DataFrame(records)
    # Resolve the target gene symbol once per ChEMBL target.
    gene_map: dict[str, str] = {}
    for target_id in activities["target_chembl_id"].dropna().unique():
        try:
            target = get_json(f"{CHEMBL_BASE}/target/{target_id}.json", timeout=60)
        except Exception:
            continue
        if str(target.get("organism") or "") != "Homo sapiens":
            continue
        gene = ""
        for component in target.get("target_components") or []:
            synonyms = component.get("target_component_synonyms") or []
            for synonym in synonyms:
                if str(synonym.get("syn_type") or "").upper().startswith(
                    "GENE_SYMBOL"
                ):
                    gene = split_gene_symbol(synonym.get("component_synonym", ""))
                    if gene:
                        break
            if gene:
                break
            for accession in component.get("target_component_xrefs") or []:
                if str(accession.get("xref_src_db", "")).lower() in {
                    "hgnc",
                    "gene symbol",
                }:
                    gene = split_gene_symbol(accession.get("xref_id", ""))
                    if gene:
                        break
            if gene:
                break
        gene_map[target_id] = gene
    activities["gene"] = activities["target_chembl_id"].map(gene_map)
    activities = activities[activities["gene"] != ""]
    if activities.empty:
        return pd.DataFrame(), {
            "status": "empty",
            "reason": "ChEMBL targets could not be mapped to gene symbols",
        }
    normalized = _normalize_human_symbols(activities["gene"].unique().tolist())
    activities["gene"] = activities["gene"].map(normalized)
    activities = activities[activities["gene"].notna() & (activities["gene"] != "")]
    grouped = (
        activities.groupby("gene", as_index=False)
        .agg(
            source=("gene", lambda _: "ChEMBL_similarity"),
            chembl_molecules=("molecule_chembl_id", lambda values: ";".join(sorted(set(values)))),
            targets=("target_chembl_id", lambda values: ";".join(sorted(set(values)))),
            max_tanimoto=("tanimoto", "max"),
            max_pchembl=("pchembl_value", "max"),
            activity_records=("gene", "size"),
        )
        .sort_values(["max_tanimoto", "max_pchembl"], ascending=False)
    )
    if grouped.empty:
        return pd.DataFrame(), {
            "status": "empty",
            "reason": "similarity hits had no human target-gene annotation",
            "n_similar_molecules": len(molecule_ids),
        }
    return grouped, {
        "status": "completed",
        "n_similar_molecules": len(molecule_ids),
        "n_activities": len(activities),
    }


def _normalize_human_symbols(symbols: list[str], workers: int = 8) -> dict[str, str | None]:
    """Resolve legacy aliases to current human HGNC symbols with MyGene.info."""
    symbols = [
        str(symbol).strip()
        for symbol in symbols
        if isinstance(symbol, str) and str(symbol).strip()
    ]

    def resolve(symbol: str) -> tuple[str, str | None]:
        params = urllib.parse.urlencode(
            {
                "q": symbol,
                "scopes": "symbol,alias",
                "species": "human",
                "fields": "symbol",
                "size": 1,
            }
        )
        try:
            payload = get_json(f"https://mygene.info/v3/query?{params}", timeout=30, retries=2)
        except Exception:
            return symbol, None
        if isinstance(payload, dict):
            hits = [payload]
        elif isinstance(payload, list):
            hits = payload
        else:
            hits = []
        for hit in hits:
            if not isinstance(hit, dict):
                continue
            value = split_gene_symbol(hit.get("symbol") or hit.get("query") or "")
            if value:
                return symbol, value
        return symbol, None

    unique = sorted(set(symbols))
    with ThreadPoolExecutor(max_workers=max(1, min(workers, len(unique) or 1))) as executor:
        return dict(executor.map(resolve, unique))


def collect_compound_targets(
    compound: dict[str, Any],
    out_dir: Path,
    *,
    cache_dir: Path | None = None,
    max_swiss_targets: int = 500,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Collect and combine accessible compound-target sources."""
    ensure_dir(out_dir)
    sources_dir = ensure_dir(out_dir / "sources")
    if cache_dir is not None:
        ensure_dir(cache_dir)
    statuses: dict[str, Any] = {}
    frames: dict[str, pd.DataFrame] = {}

    swiss_cache = (cache_dir or out_dir) / "swisstargetprediction_targets.csv"
    swiss_status_cache = (cache_dir or out_dir) / "swisstargetprediction_status.json"
    if swiss_cache.exists():
        swiss_frame = pd.read_csv(swiss_cache)
        statuses["SwissTargetPrediction"] = {
            "status": "cached",
            "source": str(swiss_cache),
        }
    elif swiss_status_cache.exists():
        swiss_frame = pd.DataFrame()
        statuses["SwissTargetPrediction"] = read_json(
            swiss_status_cache,
            {"status": "cached_empty"},
        )
    else:
        swiss_frame, swiss_status = _fetch_swiss_target(
            str(compound["canonical_smiles"])
        )
        statuses["SwissTargetPrediction"] = swiss_status
        write_json(swiss_status_cache, swiss_status)
        if not swiss_frame.empty:
            swiss_frame = swiss_frame.sort_values(
                "probability", ascending=False, na_position="last"
            ).head(max_swiss_targets)
            swiss_frame["source"] = "SwissTargetPrediction"
            swiss_frame.to_csv(swiss_cache, index=False)
    if not swiss_frame.empty:
        frames["SwissTargetPrediction"] = swiss_frame[
            ["gene", "source", "probability"]
        ].copy()

    stitch_cache = (cache_dir or out_dir) / "stitch_targets.csv"
    stitch_status_cache = (cache_dir or out_dir) / "stitch_status.json"
    if stitch_cache.exists():
        stitch_frame = pd.read_csv(stitch_cache)
        statuses["STITCH"] = {"status": "cached", "source": str(stitch_cache)}
    elif stitch_status_cache.exists():
        stitch_frame = pd.DataFrame()
        statuses["STITCH"] = read_json(
            stitch_status_cache,
            {"status": "cached_empty"},
        )
    else:
        stitch_frame, stitch_status = _fetch_stitch(str(compound["pubchem_cid"]))
        statuses["STITCH"] = stitch_status
        write_json(stitch_status_cache, stitch_status)
        if not stitch_frame.empty:
            stitch_frame.to_csv(stitch_cache, index=False)
    if not stitch_frame.empty:
        frames["STITCH"] = stitch_frame[["gene", "source", "score"]].rename(
            columns={"score": "probability"}
        )

    chembl_cache = (cache_dir or out_dir) / "chembl_similarity_targets.csv"
    chembl_status_cache = (cache_dir or out_dir) / "chembl_similarity_status.json"
    if chembl_cache.exists():
        chembl_frame = pd.read_csv(chembl_cache)
        statuses["ChEMBL_similarity"] = {
            "status": "cached",
            "source": str(chembl_cache),
        }
    elif chembl_status_cache.exists():
        chembl_frame = pd.DataFrame()
        statuses["ChEMBL_similarity"] = read_json(
            chembl_status_cache,
            {"status": "cached_empty"},
        )
    else:
        chembl_frame, chembl_status = _chembl_similarity_targets(
            str(compound["canonical_smiles"])
        )
        statuses["ChEMBL_similarity"] = chembl_status
        write_json(chembl_status_cache, chembl_status)
        if not chembl_frame.empty:
            chembl_frame.to_csv(chembl_cache, index=False)
    if not chembl_frame.empty:
        chembl_subset = chembl_frame.copy()
        chembl_subset["probability"] = pd.to_numeric(
            chembl_subset.get("max_pchembl", np.nan), errors="coerce"
        )
        frames["ChEMBL_similarity"] = chembl_subset[
            ["gene", "source", "probability"]
        ].copy()

    for name in ("SwissTargetPrediction", "STITCH", "ChEMBL_similarity"):
        frame = frames.get(name)
        if frame is None:
            pd.DataFrame(columns=["gene", "source", "probability"]).to_csv(
                sources_dir / f"{name}.csv",
                index=False,
            )
            continue
        frame = frame.copy()
        frame["gene"] = frame["gene"].map(lambda value: split_gene_symbol(value))
        frame = frame[frame["gene"] != ""].drop_duplicates("gene")
        frame.to_csv(sources_dir / f"{name}.csv", index=False)
        frames[name] = frame

    all_rows: list[dict[str, Any]] = []
    all_genes = sorted(set().union(*(set(frame["gene"]) for frame in frames.values())))
    for gene in all_genes:
        gene_sources: list[str] = []
        scores: list[float] = []
        for name, frame in frames.items():
            rows = frame[frame["gene"] == gene]
            if rows.empty:
                continue
            gene_sources.append(name)
            value = pd.to_numeric(rows["probability"], errors="coerce").max()
            if pd.notna(value):
                scores.append(float(value))
        all_rows.append(
            {
                "gene": gene,
                "source_count": len(gene_sources),
                "sources": ";".join(gene_sources),
                "best_score": max(scores) if scores else np.nan,
            }
        )
    combined = pd.DataFrame(
        all_rows,
        columns=["gene", "source_count", "sources", "best_score"],
    ).sort_values(["source_count", "best_score", "gene"], ascending=[False, False, True])
    combined.to_csv(out_dir / "compound_targets.csv", index=False)
    write_json(
        out_dir / "compound_target_sources.json",
        {
            "status": "completed" if not combined.empty else "empty",
            "sources": statuses,
            "targets": int(len(combined)),
        },
    )
    return combined.reset_index(drop=True), statuses


def open_targets_disease_targets(
    disease_name: str,
    *,
    page_size: int = 500,
    min_score: float = 0.0,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Fetch disease-associated targets from Open Targets GraphQL."""
    query = """
    query searchDisease($q: String!) {
      search(queryString: $q, entityNames: ["disease"]) {
        hits {
          id
          name
          entity
        }
      }
    }
    """
    payload = _graphql(query, {"q": disease_name})
    hits = (((payload.get("data") or {}).get("search") or {}).get("hits") or [])
    disease = next((hit for hit in hits if hit.get("entity") == "disease"), None)
    if not disease:
        return pd.DataFrame(), {
            "status": "empty",
            "reason": f"disease not found: {disease_name}",
        }
    query_targets = """
    query diseaseTargets($id: String!, $size: Int!) {
      disease(efoId: $id) {
        id
        name
        associatedTargets(page: {index: 0, size: $size}) {
          count
          rows {
            score
            target {
              id
              approvedSymbol
            }
          }
        }
      }
    }
    """
    payload = _graphql(
        query_targets,
        {"id": disease["id"], "size": int(page_size)},
    )
    disease_data = (payload.get("data") or {}).get("disease") or {}
    rows = ((disease_data.get("associatedTargets") or {}).get("rows") or [])
    records = [
        {
            "gene": split_gene_symbol((row.get("target") or {}).get("approvedSymbol")),
            "ensembl_id": (row.get("target") or {}).get("id") or "",
            "score": float(row.get("score") or 0.0),
            "source": "OpenTargets",
        }
        for row in rows
    ]
    frame = pd.DataFrame(records)
    if frame.empty:
        return frame, {
            "status": "empty",
            "disease_id": disease.get("id"),
            "disease_name": disease.get("name"),
        }
    frame = frame[frame["gene"] != ""]
    frame = frame[pd.to_numeric(frame["score"], errors="coerce") >= min_score]
    frame = (
        frame.sort_values("score", ascending=False)
        .drop_duplicates("gene")
        .reset_index(drop=True)
    )
    return frame, {
        "status": "completed",
        "disease_id": disease.get("id"),
        "disease_name": disease.get("name"),
        "count": int(len(frame)),
    }


def gwas_catalog_disease_targets(
    disease_terms: list[str],
    *,
    max_records: int = 1000,
    timeout: int = 120,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Collect disease genes from GWAS Catalog as an open, auditable source."""
    try:
        from evidence.connectors import GWASCatalogConnector
        from evidence.context import EvidenceContext
    except Exception as exc:  # noqa: BLE001
        return pd.DataFrame(), {
            "status": "failed",
            "reason": f"evidence connector unavailable: {exc}",
        }
    rows: list[dict[str, Any]] = []
    statuses: dict[str, Any] = {}
    for term in disease_terms:
        try:
            records = GWASCatalogConnector().collect(
                EvidenceContext(
                    disease={"name": term},
                    max_records_per_source=max_records,
                    timeout_seconds=timeout,
                    allow_network=True,
                )
            )
            statuses[term] = {"status": "completed", "count": len(records)}
        except Exception as exc:  # noqa: BLE001
            statuses[term] = {"status": "failed", "reason": str(exc)}
            continue
        for record in records:
            rows.append(
                {
                    "gene": split_gene_symbol(record.target_symbol),
                    "score": float(record.score or 0.0),
                    "source": "GWAS_Catalog",
                }
            )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame, {
            "status": "empty",
            "terms": statuses,
        }
    frame = (
        frame[frame["gene"] != ""]
        .groupby("gene", as_index=False)
        .agg(
            score=("score", "max"),
            source=("source", lambda _: "GWAS_Catalog"),
        )
        .sort_values(["score", "gene"], ascending=[False, True])
    )
    return frame, {
        "status": "completed",
        "count": int(len(frame)),
        "terms": statuses,
    }


def clinvar_disease_targets(
    target_symbols: list[str],
    *,
    max_records: int = 500,
    timeout: int = 120,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Collect ClinVar evidence for candidate disease genes."""
    try:
        from evidence.connectors import ClinVarConnector
        from evidence.context import EvidenceContext
    except Exception as exc:  # noqa: BLE001
        return pd.DataFrame(), {
            "status": "failed",
            "reason": f"evidence connector unavailable: {exc}",
        }
    symbols = [
        split_gene_symbol(value)
        for value in target_symbols
        if split_gene_symbol(value)
    ]
    if not symbols:
        return pd.DataFrame(), {
            "status": "empty",
            "reason": "no target symbols supplied",
        }
    try:
        records = ClinVarConnector().collect(
            EvidenceContext(
                disease={"name": "NAFLD"},
                target_symbols=symbols,
                max_records_per_source=max_records,
                timeout_seconds=timeout,
                allow_network=True,
            )
        )
    except Exception as exc:  # noqa: BLE001
        return pd.DataFrame(), {"status": "failed", "reason": str(exc)}
    rows = [
        {
            "gene": split_gene_symbol(record.target_symbol),
            "score": float(record.score or 0.0),
            "source": "ClinVar",
        }
        for record in records
    ]
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame, {"status": "empty"}
    frame = (
        frame[frame["gene"] != ""]
        .groupby("gene", as_index=False)
        .agg(
            score=("score", "max"),
            source=("source", lambda _: "ClinVar"),
        )
        .sort_values(["score", "gene"], ascending=[False, True])
    )
    return frame, {"status": "completed", "count": int(len(frame))}


def _graphql(query: str, variables: dict[str, Any]) -> dict[str, Any]:
    body = json.dumps({"query": query, "variables": variables}).encode("utf-8")
    request = urllib.request.Request(
        OPEN_TARGETS_GRAPHQL,
        data=body,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "Mozilla/5.0",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8", "replace"))


def load_disease_source_file(
    name: str,
    path: str | Path | None,
    *,
    gene_column: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load a user-supplied disease target source without guessing values."""
    if path in (None, ""):
        return pd.DataFrame(), {
            "status": "unavailable",
            "reason": (
                f"{name} bulk download requires licensed/credentialed access; "
                "provide a local CSV/TSV to enable this source"
            ),
        }
    source = Path(path).expanduser()
    if not source.is_absolute():
        source = (Path.cwd() / source).resolve()
    if not source.exists():
        return pd.DataFrame(), {"status": "failed", "reason": f"file not found: {source}"}
    sep = "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(source, sep=sep, dtype=str)
    gene_column = gene_column or _guess_gene_column(frame)
    if not gene_column or gene_column not in frame.columns:
        return pd.DataFrame(), {
            "status": "failed",
            "reason": f"no gene column in {source}: {list(frame.columns)}",
        }
    output = pd.DataFrame(
        {
            "gene": frame[gene_column].map(split_gene_symbol),
            "score": 1.0,
            "source": name,
        }
    )
    output = output[output["gene"] != ""].drop_duplicates("gene")
    return output.reset_index(drop=True), {
        "status": "completed",
        "source": str(source),
        "count": int(len(output)),
    }


def load_compound_target_file(
    name: str,
    path: str | Path | None,
    *,
    gene_column: str | None = None,
    score_column: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Load a user-supplied compound-target prediction table."""
    if path in (None, ""):
        return pd.DataFrame(), {
            "status": "not_configured",
            "reason": "no local compound-target prediction table supplied",
        }
    source = Path(path).expanduser()
    if not source.is_absolute():
        source = (Path.cwd() / source).resolve()
    if not source.exists():
        return pd.DataFrame(), {
            "status": "failed",
            "reason": f"file not found: {source}",
        }
    sep = "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","
    frame = pd.read_csv(source, sep=sep, dtype=str)
    gene_column = gene_column or _guess_gene_column(frame)
    if not gene_column or gene_column not in frame.columns:
        return pd.DataFrame(), {
            "status": "failed",
            "reason": f"no gene column in {source}: {list(frame.columns)}",
        }
    if score_column and score_column in frame.columns:
        score = pd.to_numeric(frame[score_column], errors="coerce")
    else:
        score = pd.Series(1.0, index=frame.index, dtype=float)
    output = pd.DataFrame(
        {
            "gene": frame[gene_column].map(split_gene_symbol),
            "score": score.fillna(0.0),
            "source": name,
        }
    )
    output = (
        output[output["gene"] != ""]
        .groupby("gene", as_index=False)
        .agg(
            score=("score", "max"),
            source=("source", lambda _: name),
        )
    )
    return output.reset_index(drop=True), {
        "status": "completed",
        "source": str(source),
        "count": int(len(output)),
    }


def _guess_gene_column(frame: pd.DataFrame) -> str | None:
    aliases = {
        "gene",
        "gene_symbol",
        "symbol",
        "hgnc",
        "hugo",
        "approved_symbol",
    }
    for column in frame.columns:
        if str(column).strip().lower().replace(" ", "_") in aliases:
            return str(column)
    return None


def combine_disease_sources(
    sources: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    genes = sorted(set().union(*(set(frame["gene"]) for frame in sources.values())))
    rows = []
    for gene in genes:
        present = sorted(name for name, frame in sources.items() if gene in set(frame["gene"]))
        scores = []
        for name in present:
            values = pd.to_numeric(sources[name].loc[sources[name]["gene"] == gene, "score"], errors="coerce")
            if not values.empty and pd.notna(values.max()):
                scores.append(float(values.max()))
        rows.append(
            {
                "gene": gene,
                "source_count": len(present),
                "sources": ";".join(present),
                "best_score": max(scores) if scores else np.nan,
            }
        )
    return pd.DataFrame(
        rows,
        columns=["gene", "source_count", "sources", "best_score"],
    ).sort_values(["source_count", "best_score", "gene"], ascending=[False, False, True])


def make_venn_figure(
    sets: dict[str, set[str]],
    output: Path,
    *,
    title: str,
) -> Path:
    """Draw a publication-readable 2- or 3-set Venn-like diagram."""
    from matplotlib.patches import Circle

    labels = list(sets)[:3]
    values = [sets[label] for label in labels]
    if len(values) < 2:
        raise ValueError("Venn diagram requires at least two gene sets")
    fig, ax = plt.subplots(figsize=(6.8, 5.5))
    ax.set_aspect("equal")
    ax.axis("off")
    colors = ["#2f6bb3", "#e07a3f", "#4d9b6a"]
    if len(values) == 2:
        centers = [(-0.55, 0.0), (0.55, 0.0)]
        radius = 1.1
        regions = {
            ("A",): len(values[0] - values[1]),
            ("B",): len(values[1] - values[0]),
            ("AB",): len(values[0] & values[1]),
        }
        positions = {
            ("A",): (-1.15, 0),
            ("B",): (1.15, 0),
            ("AB",): (0, 0),
        }
    else:
        centers = [
            (0.0, 0.95),
            (-0.95, -0.55),
            (0.95, -0.55),
        ]
        radius = 1.28
        regions = _venn_counts(values)
        positions = {
            "A": (-1.55, 1.2),
            "B": (-1.7, -1.25),
            "C": (1.7, -1.25),
            "AB": (-0.55, 0.5),
            "AC": (0.55, 0.5),
            "BC": (0.0, -1.1),
            "ABC": (0.0, 0.1),
        }
    for center, color in zip(centers, colors):
        ax.add_patch(
            Circle(center, radius, facecolor=color, edgecolor=color, alpha=0.15, linewidth=2)
        )
    for key, count in regions.items():
        x, y = positions.get(key, (0, 0))
        ax.text(x, y, str(count), ha="center", va="center", fontsize=12, fontweight="bold")
    for (x, y), label in zip(centers, labels):
        ax.text(x, y + radius + 0.12, label, ha="center", va="bottom", fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlim(-2.3, 2.3)
    ax.set_ylim(-2.1, 2.4)
    output.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, output)
    return output


def make_workflow_figure(output: Path) -> Path:
    """Draw the computational workflow from the supplied experiment plan."""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    boxes = [
        (0.04, 0.66, 0.20, 0.17, "6PPD-Q\nstructure + target prediction", "#dbe8f0"),
        (0.29, 0.66, 0.20, 0.17, "NAFLD\npublic target sets", "#e8e1d0"),
        (0.54, 0.66, 0.17, 0.17, "Intersection\ncandidate targets", "#dce8dc"),
        (0.76, 0.66, 0.20, 0.17, "STRING PPI\nhub ranking", "#e6ddec"),
        (0.10, 0.33, 0.22, 0.17, "Bulk microarray\ntraining and validation", "#dbe8f0"),
        (0.39, 0.33, 0.22, 0.17, "11-model CV\n+ SHAP", "#e8e1d0"),
        (0.68, 0.33, 0.22, 0.17, "Mouse scRNA\n+ human snRNA", "#dce8dc"),
        (0.24, 0.05, 0.22, 0.15, "Molecular docking\nVina", "#e6ddec"),
        (0.54, 0.05, 0.22, 0.15, "100 ns MD\n+ MM-PBSA prep", "#eadfdc"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    centers: dict[int, tuple[float, float]] = {}
    for index, (x, y, width, height, label, color) in enumerate(boxes):
        patch = FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            facecolor=color,
            edgecolor="#62717d",
            linewidth=1.1,
        )
        ax.add_patch(patch)
        ax.text(
            x + width / 2,
            y + height / 2,
            label,
            ha="center",
            va="center",
            fontsize=6.4,
        )
        centers[index] = (x + width / 2, y + height / 2)
    arrows = [
        (0, 2),
        (1, 2),
        (2, 3),
        (2, 4),
        (3, 5),
        (4, 5),
        (5, 6),
        (5, 7),
        (7, 8),
    ]
    for source, target in arrows:
        start = centers[source]
        end = centers[target]
        if abs(start[0] - end[0]) < 1e-6 and abs(start[1] - end[1]) < 1e-6:
            continue
        arrow = FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.1,
            color="#7b8794",
            connectionstyle="arc3,rad=0.08",
            shrinkA=35,
            shrinkB=35,
        )
        ax.add_patch(arrow)
    ax.set_title("Experiment plan one computational workflow", fontsize=13, fontweight="bold")
    output.parent.mkdir(parents=True, exist_ok=True)
    save_figure(fig, output)
    return output


def _venn_counts(values: list[set[str]]) -> dict[str, int]:
    a, b, c = values
    return {
        "A": len(a - b - c),
        "B": len(b - a - c),
        "C": len(c - a - b),
        "AB": len((a & b) - c),
        "AC": len((a & c) - b),
        "BC": len((b & c) - a),
        "ABC": len(a & b & c),
    }
