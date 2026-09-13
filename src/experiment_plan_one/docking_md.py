"""Docking and GROMACS/MM-PBSA preparation for experiment plan one."""

from __future__ import annotations

import json
import logging
import math
import os
import shutil
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .common import LOG, download_file, ensure_dir, get_json, save_figure, write_json

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_FILE = "https://files.rcsb.org/download/{pdb_id}.pdb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files/{entry_id}-model_v4.pdb"


def uniprot_for_gene(gene: str, *, species: int = 9606) -> dict[str, Any]:
    """Resolve one gene symbol to a reviewed UniProt entry."""
    query = f"gene_exact:{gene} AND organism_id:{species} AND reviewed:true"
    params = urllib.parse.urlencode(
        {
            "query": query,
            "format": "json",
            "size": 5,
            "fields": "accession,id,protein_name,gene_names,sequence",
        }
    )
    payload = get_json(f"{UNIPROT_SEARCH}?{params}", timeout=90)
    results = payload.get("results") or []
    if not results:
        raise RuntimeError(f"no reviewed UniProt entry for {gene}")
    entry = results[0]
    sequences = entry.get("sequence") or {}
    return {
        "gene": gene.upper(),
        "uniprot_accession": entry.get("primaryAccession") or "",
        "uniprot_entry": entry.get("uniProtkbId") or "",
        "protein_name": _protein_name(entry),
        "sequence_length": int(sequences.get("length") or 0),
        "sequence": sequences.get("value") or "",
    }


def _protein_name(entry: dict[str, Any]) -> str:
    description = entry.get("proteinDescription") or {}
    recommended = description.get("recommendedName") or {}
    full_name = recommended.get("fullName") or {}
    return full_name.get("value") or entry.get("uniProtkbId") or ""


def best_pdb_for_target(uniprot: str, *, min_resolution: float = 3.8) -> dict[str, Any] | None:
    """Return the best high-resolution experimental PDB entry from RCSB."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": 30},
            "results_content_type": ["experimental"],
            "sort": [
                {"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}
            ],
        },
    }
    request = urllib.request.Request(
        RCSB_SEARCH,
        data=json.dumps(query).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = json.loads(response.read().decode("utf-8", "replace"))
    except Exception as exc:  # noqa: BLE001
        LOG.warning("RCSB search failed for %s: %s", uniprot, exc)
        return None
    identifiers = [
        row.get("identifier")
        for row in payload.get("result_set") or []
        if row.get("identifier")
    ]
    if not identifiers:
        return None
    pdb_id = str(identifiers[0]).upper()
    entry = get_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", timeout=60)
    info = entry.get("rcsb_entry_info") or {}
    methods = (entry.get("exptl") or [{}])[0].get("method") or ""
    resolution_values = info.get("resolution_combined") or []
    resolution = min(resolution_values) if resolution_values else None
    if resolution is not None and resolution > min_resolution:
        return None
    return {
        "pdb_id": pdb_id,
        "method": methods,
        "resolution_angstrom": resolution,
        "title": (entry.get("struct") or {}).get("title") or "",
    }


def alphafold_structure(uniprot: str) -> dict[str, Any] | None:
    try:
        payload = get_json(ALPHAFOLD_API.format(accession=uniprot), timeout=90)
    except Exception as exc:  # noqa: BLE001
        LOG.warning("AlphaFold lookup failed for %s: %s", uniprot, exc)
        return None
    if not isinstance(payload, list) or not payload:
        return None
    entry = payload[0]
    return {
        "entry_id": entry.get("entryId") or "",
        "pdb_url": entry.get("pdbUrl") or ALPHAFOLD_FILE.format(entry_id=entry.get("entryId") or ""),
        "mean_plddt": entry.get("globalMetricValue"),
    }


def write_ligand_sdf(
    cid: str,
    smiles: str,
    output_path: Path,
) -> Path:
    """Write a 3D SDF for 6PPD-Q, preferring the PubChem conformer."""
    ensure_dir(output_path.parent)
    if output_path.exists():
        return output_path
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    try:
        download_file(url, output_path, timeout=120)
        if output_path.stat().st_size > 100:
            return output_path
    except Exception:
        pass
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260913
    if AllChem.EmbedMolecule(molecule, params) != 0:
        raise RuntimeError("RDKit failed to generate the 6PPD-Q 3D conformer")
    AllChem.MMFFOptimizeMolecule(molecule, maxIters=2000)
    writer = Chem.SDWriter(str(output_path))
    writer.write(molecule)
    writer.close()
    return output_path


def _prepare_receptor_with_existing_tools(
    receptor_pdb: Path,
    output_path: Path,
    cfg,
    log,
) -> None:
    from docking.receptor import prepare_receptor

    prepare_receptor(cfg, log)
    if not output_path.exists():
        raise RuntimeError(f"receptor preparation did not create {output_path}")


def run_docking_for_targets(
    targets: list[dict[str, str]],
    compound: dict[str, Any],
    output_dir: Path,
    *,
    exhaustiveness: int = 16,
    cpu: int = 4,
    timeout: int = 3600,
) -> dict[str, Any]:
    """Prepare and dock each target, keeping each target in its own workdir."""
    from docking.analysis import analyze_results
    from docking.config import ResolvedConfig
    from docking.docking import run_docking
    from docking.ligands import prepare_ligands
    from docking.receptor import prepare_receptor
    from docking.utils import setup_logging

    ensure_dir(output_dir)
    ligand = write_ligand_sdf(
        str(compound["pubchem_cid"]),
        str(compound["canonical_smiles"]),
        output_dir / "ligands" / "6PPD-Q.sdf",
    )
    rows: list[dict[str, Any]] = []
    target_records: list[dict[str, Any]] = []
    for target in targets:
        gene = str(target["gene"]).upper()
        target_dir = ensure_dir(output_dir / "targets" / gene)
        data: dict[str, Any] = {
            "workdir": str(target_dir),
            "output_dir": "outputs/run_001",
            "receptor": {
                "input": "",
                "output": "data/receptors/receptor.pdbqt",
                "center": [0.0, 0.0, 0.0],
                "size": [24.0, 24.0, 24.0],
                "flexible": [],
            },
            "ligand": {
                "input": str(ligand),
                "output_dir": "data/ligands/prepared",
                "smiles_column": "SMILES",
                "id_column": "ID",
                "ph": 7.4,
                "remove_salts": True,
                "neutralize": True,
                "max_heavy_atoms": 80,
                "max_rotatable_bonds": 20,
                "max_ligands": None,
                "conformers": 1,
                "seed": 42,
                "engine": "auto",
            },
            "docking": {
                "engine": "vina",
                "executable": _find_vina(),
                "scoring": "vina",
                "exhaustiveness": int(exhaustiveness),
                "num_modes": 9,
                "energy_range": 3.0,
                "cpu": int(cpu),
                "max_workers": 1,
                "seed": 42,
                "seeds": [],
                "replicates": 3,
                "timeout_seconds": int(timeout),
                "resume": True,
            },
            "analysis": {
                "cutoff": -7.0,
                "moderate_cutoff": -5.0,
                "top_n": 3,
                "figures": True,
                "diversity": False,
                "tanimoto_cutoff": 0.7,
            },
            "md_simulation": {
                "mode": "prepare",
                "top_n": 1,
                "prod_steps": 50_000_000,
                "dt_ps": 0.002,
                "temperature": 310.0,
                "pressure": 1.0,
                "protein_forcefield": "amber99sb-ildn",
                "water": "tip3p",
                "box_padding_nm": 1.2,
                "ion_concentration": 0.15,
                "timeout_seconds": 172800,
                "cpu": max(1, cpu),
                "gpu": True,
                "figures": True,
                "equil_steps": 250_000,
                "mmpbsa_command": None,
            },
            "knockout": {},
            "insilico_knockout": {},
        }

        structure = _resolve_structure(gene, target_dir, target)
        if structure is None:
            target_records.append(
                {
                    **target,
                    "status": "skipped",
                    "reason": "no experimental or AlphaFold structure",
                }
            )
            continue
        receptor_pdb = structure["receptor_pdb"]
        center, size, cofactor = _docking_box(Path(structure["raw_pdb"]))
        data["receptor"]["input"] = str(receptor_pdb)
        data["receptor"]["center"] = center
        data["receptor"]["size"] = size
        receptor_copy = target_dir / "data" / "receptors" / "receptor.pdb"
        ensure_dir(receptor_copy.parent)
        shutil.copy2(receptor_pdb, receptor_copy)
        config_path = target_dir / "docking_config.json"
        write_json(config_path, data)
        cfg = ResolvedConfig(data, config_path)
        log = setup_logging(str(target_dir / "docking.log"))
        try:
            prepare_receptor(cfg, log)
            prepare_ligands(cfg, log)
            run_docking(cfg, log)
            summary = analyze_results(cfg, log)
        except Exception as exc:  # noqa: BLE001
            LOG.warning("docking failed for %s: %s", gene, exc)
            target_records.append(
                {
                    **target,
                    "status": "failed",
                    "reason": str(exc),
                    "structure": structure,
                }
            )
            continue
        results_path = cfg.results_path()
        energies = pd.read_csv(results_path) if results_path.exists() else pd.DataFrame()
        energy_column = next(
            (
                column
                for column in ("affinity_kcal_mol", "affinity", "binding_affinity")
                if column in energies.columns
            ),
            None,
        )
        best_energy = (
            float(pd.to_numeric(energies[energy_column], errors="coerce").min())
            if not energies.empty and energy_column
            else np.nan
        )
        rows.append(
            {
                "gene": gene,
                "uniprot_accession": target.get("uniprot_accession", ""),
                "structure_source": structure["source"],
                "pdb_id": structure.get("pdb_id", ""),
                "best_affinity_kcal_mol": best_energy,
                "status": "completed",
                "cofactor_box": cofactor,
                "workdir": str(target_dir),
                "results_csv": str(results_path),
            }
        )
        target_records.append(
            {
                **target,
                **rows[-1],
                "config": str(config_path),
                "structure": structure,
            }
        )
    score_frame = pd.DataFrame(rows)
    score_frame.to_csv(output_dir / "docking_scores.csv", index=False)
    if not score_frame.empty:
        _plot_docking_heatmap(score_frame, output_dir / "fig5c_docking_affinity_heatmap.png")
        best_target = score_frame.sort_values("best_affinity_kcal_mol").iloc[0]
        best_dir = Path(str(best_target["workdir"]))
        _write_docking_pose_figure(
            best_dir,
            output_dir / "fig5a_docking_pose_3d.png",
            str(best_target["gene"]),
        )
        _write_interaction_figure(best_dir, output_dir / "fig5b_interaction_schematic.png", str(best_target["gene"]))
    write_json(
        output_dir / "docking_summary.json",
        {
            "targets_requested": len(targets),
            "targets_completed": len(rows),
            "targets": target_records,
            "exhaustiveness": exhaustiveness,
            "seed_replicates": 3,
        },
    )
    return {
        "scores": output_dir / "docking_scores.csv",
        "records": target_records,
        "completed": score_frame,
    }


def _find_vina() -> str:
    found = shutil.which("vina") or shutil.which("vina.exe")
    if found:
        return found
    local = Path(__file__).resolve().parents[2] / "dock" / "tools" / "vina.exe"
    return str(local) if local.exists() else "vina"


def _resolve_structure(gene: str, target_dir: Path, target: dict[str, str]) -> dict[str, Any] | None:
    accession = str(target.get("uniprot_accession") or "").strip()
    protein_name = str(target.get("protein_name") or "")
    if not accession:
        try:
            info = uniprot_for_gene(gene)
            accession = info["uniprot_accession"]
            protein_name = info["protein_name"]
            target.update(info)
        except Exception as exc:  # noqa: BLE001
            LOG.warning("UniProt mapping failed for %s: %s", gene, exc)
            return None
    pdb = best_pdb_for_target(accession)
    if pdb:
        raw_pdb = target_dir / "structure" / f"{pdb['pdb_id']}.pdb"
        download_file(RCSB_FILE.format(pdb_id=pdb["pdb_id"]), raw_pdb, timeout=180)
        protein_pdb = _strip_pdb(raw_pdb, target_dir / "structure" / f"{pdb['pdb_id']}_protein.pdb")
        return {
            "source": "RCSB PDB",
            **pdb,
            "raw_pdb": raw_pdb,
            "receptor_pdb": protein_pdb,
            "protein_name": protein_name,
        }
    alphafold = alphafold_structure(accession)
    if not alphafold:
        return None
    raw_pdb = target_dir / "structure" / f"{gene}_AlphaFold.pdb"
    download_file(alphafold["pdb_url"], raw_pdb, timeout=180)
    protein_pdb = _strip_pdb(raw_pdb, target_dir / "structure" / f"{gene}_AlphaFold_protein.pdb")
    return {
        "source": "AlphaFold",
        **alphafold,
        "raw_pdb": raw_pdb,
        "receptor_pdb": protein_pdb,
        "protein_name": protein_name,
    }


def _strip_pdb(input_path: Path, output_path: Path) -> Path:
    ensure_dir(output_path.parent)
    with input_path.open("r", encoding="utf-8", errors="replace") as source, output_path.open(
        "w", encoding="utf-8"
    ) as destination:
        for line in source:
            record = line[0:6].strip()
            if record in {"ATOM", "TER", "END", "MODEL", "ENDMDL"}:
                destination.write(line)
    return output_path


def _docking_box(receptor_pdb: Path) -> tuple[list[float], list[float], str]:
    atoms, hetero = _parse_pdb_coordinates(receptor_pdb)
    if atoms.size == 0:
        raise ValueError(f"no atoms in receptor PDB: {receptor_pdb}")
    if hetero.size:
        center = np.median(hetero, axis=0)
        size = [22.0, 22.0, 22.0]
        cofactor = "cocrystallized ligand/cofactor"
    else:
        center = (atoms.min(axis=0) + atoms.max(axis=0)) / 2.0
        size = [24.0, 24.0, 24.0]
        cofactor = "blind docking"
    return [round(float(value), 3) for value in center], size, cofactor


def _parse_pdb_coordinates(path: Path) -> tuple[np.ndarray, np.ndarray]:
    protein: list[list[float]] = []
    hetero: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                coordinate = [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]
            except ValueError:
                continue
            if line.startswith("HETATM"):
                residue = line[17:20].strip().upper()
                if residue not in {
                    "HOH",
                    "WAT",
                    "NA",
                    "K",
                    "CL",
                    "MG",
                    "CA",
                    "ZN",
                    "MN",
                    "FE",
                }:
                    hetero.append(coordinate)
            else:
                protein.append(coordinate)
    return np.asarray(protein), np.asarray(hetero)


def _plot_docking_heatmap(scores: pd.DataFrame, output: Path) -> None:
    frame = scores.dropna(subset=["best_affinity_kcal_mol"]).sort_values("best_affinity_kcal_mol")
    if frame.empty:
        return
    fig, ax = plt.subplots(figsize=(6.2, max(3.5, len(frame) * 0.45)))
    image = ax.imshow(
        frame[["best_affinity_kcal_mol"]].to_numpy(),
        cmap="YlGnBu_r",
        aspect="auto",
    )
    ax.set_xticks([0])
    ax.set_xticklabels(["Vina affinity"])
    ax.set_yticks(range(len(frame)))
    ax.set_yticklabels(frame["gene"])
    for row, value in enumerate(frame["best_affinity_kcal_mol"]):
        ax.text(0, row, f"{value:.2f}", ha="center", va="center", color="white", fontsize=8)
    ax.set_title("6PPD-Q docking affinity", fontweight="bold")
    fig.colorbar(image, ax=ax, label="kcal/mol", shrink=0.7)
    save_figure(fig, output)


def _write_interaction_figure(target_dir: Path, output: Path, gene: str) -> None:
    receptor = target_dir / "data" / "receptors" / "receptor.pdb"
    pose_table = target_dir / "outputs" / "run_001" / "docked" / "results.csv"
    if not receptor.exists() or not pose_table.exists():
        return
    results = pd.read_csv(pose_table)
    energy_column = next(
        (
            column
            for column in ("affinity_kcal_mol", "affinity", "binding_affinity")
            if column in results.columns
        ),
        None,
    )
    if energy_column:
        results["_energy"] = pd.to_numeric(results[energy_column], errors="coerce")
        results = results.sort_values("_energy")
    pose_path = Path(str(results.iloc[0].get("pose_file") or ""))
    if not pose_path.is_absolute():
        pose_path = target_dir / pose_path
    if not pose_path.exists():
        return
    ligand_atoms = _parse_pdbqt_coordinates(pose_path)
    contacts = _receptor_contacts(receptor, ligand_atoms)
    contacts.to_csv(output.with_suffix(".csv"), index=False)
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.axis("off")
    ax.text(
        0.5,
        0.96,
        f"{gene}-6PPD-Q predicted contacts",
        ha="center",
        va="top",
        fontsize=12,
        fontweight="bold",
        transform=ax.transAxes,
    )
    ax.text(
        0.5,
        0.5,
        "Ligand",
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        bbox={"boxstyle": "round,pad=0.7", "facecolor": "#d7e4ee", "edgecolor": "#385a73"},
        transform=ax.transAxes,
    )
    if contacts.empty:
        ax.text(0.5, 0.42, "No contact within 4.5 A", ha="center", transform=ax.transAxes)
    else:
        top = contacts.head(14)
        for index, row in top.iterrows():
            angle = 2 * math.pi * index / max(len(top), 1)
            x = 0.5 + 0.38 * math.cos(angle)
            y = 0.5 + 0.38 * math.sin(angle)
            ax.annotate(
                f"{row['residue']}{row['residue_number']}",
                xy=(0.5, 0.5),
                xytext=(x, y),
                ha="center",
                va="center",
                fontsize=8,
                arrowprops={"arrowstyle": "-", "color": "#7b8794", "lw": 0.8},
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#a7b1bb"},
                xycoords=ax.transAxes,
                textcoords=ax.transAxes,
            )
    ax.text(
        0.5,
        0.02,
        "Contacts are geometric approximations; Discovery Studio interaction typing was not used.",
        ha="center",
        va="bottom",
        fontsize=7.5,
        color="#5d6670",
        transform=ax.transAxes,
    )
    save_figure(fig, output)


def _write_docking_pose_figure(target_dir: Path, output: Path, gene: str) -> None:
    receptor = target_dir / "data" / "receptors" / "receptor.pdb"
    pose_table = target_dir / "outputs" / "run_001" / "docked" / "results.csv"
    if not receptor.exists() or not pose_table.exists():
        return
    results = pd.read_csv(pose_table)
    energy_column = next(
        (
            column
            for column in ("affinity_kcal_mol", "affinity", "binding_affinity")
            if column in results.columns
        ),
        None,
    )
    if energy_column:
        results["_energy"] = pd.to_numeric(results[energy_column], errors="coerce")
        results = results.sort_values("_energy")
    pose_path = Path(str(results.iloc[0].get("pose_file") or ""))
    if not pose_path.is_absolute():
        pose_path = target_dir / pose_path
    ligand = _parse_pdbqt_coordinates(pose_path)
    if ligand.size == 0:
        return
    protein_atoms: list[np.ndarray] = []
    with receptor.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            try:
                coordinate = [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]
            except ValueError:
                continue
            if np.linalg.norm(np.asarray(coordinate) - ligand.mean(axis=0)) <= 18.0:
                protein_atoms.append(np.asarray(coordinate))
    protein = np.asarray(protein_atoms)
    fig = plt.figure(figsize=(6.4, 5.8))
    ax = fig.add_subplot(111, projection="3d")
    if protein.size:
        ax.scatter(
            protein[:, 0],
            protein[:, 1],
            protein[:, 2],
            s=10,
            c="#7b8794",
            alpha=0.55,
            depthshade=False,
            label="Protein C-alpha",
        )
    ax.scatter(
        ligand[:, 0],
        ligand[:, 1],
        ligand[:, 2],
        s=36,
        c="#cf4f45",
        depthshade=False,
        label="6PPD-Q",
    )
    ax.set_xlabel("x (A)")
    ax.set_ylabel("y (A)")
    ax.set_zlabel("z (A)")
    ax.set_title(f"{gene}-6PPD-Q best docking pose", fontweight="bold")
    ax.legend(loc="upper right")
    ax.view_init(elev=18, azim=42)
    save_figure(fig, output)


def _parse_pdbqt_coordinates(path: Path) -> np.ndarray:
    coordinates: list[list[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                coordinates.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
            except ValueError:
                continue
    return np.asarray(coordinates)


def _receptor_contacts(receptor: Path, ligand: np.ndarray, cutoff: float = 4.5) -> pd.DataFrame:
    protein_atoms: list[tuple[str, int, str, list[float]]] = []
    with receptor.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            try:
                coordinate = [
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ]
            except ValueError:
                continue
            residue = line[17:20].strip()
            chain = line[21:22].strip()
            try:
                residue_number = int(line[22:26])
            except ValueError:
                continue
            protein_atoms.append((residue, residue_number, chain, coordinate))
    if not protein_atoms or ligand.size == 0:
        return pd.DataFrame(columns=["residue", "residue_number", "chain", "min_distance"])
    tree = cKDTree(ligand)
    records: dict[tuple[str, int, str], float] = {}
    for residue, residue_number, chain, coordinate in protein_atoms:
        distance = float(tree.query(coordinate)[0])
        if distance <= cutoff:
            key = (residue, residue_number, chain)
            records[key] = min(records.get(key, float("inf")), distance)
    return (
        pd.DataFrame(
            [
                {
                    "residue": residue,
                    "residue_number": residue_number,
                    "chain": chain,
                    "min_distance": distance,
                }
                for (residue, residue_number, chain), distance in records.items()
            ]
        )
        .sort_values("min_distance")
        .reset_index(drop=True)
        if records
        else pd.DataFrame(columns=["residue", "residue_number", "chain", "min_distance"])
    )


def prepare_or_run_md(
    docking_dir: Path,
    *,
    run: bool,
    gpu: bool,
    cpu: int = 4,
    timeout: int = 172800,
) -> dict[str, Any]:
    """Run the project's GROMACS wrapper for the best docked target."""
    from docking.config import ResolvedConfig
    from docking.md_simulation import run_md_simulation
    from docking.utils import setup_logging

    scores_path = docking_dir / "docking_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(f"docking scores not found: {scores_path}")
    scores = pd.read_csv(scores_path)
    if scores.empty:
        return {"status": "skipped", "reason": "no completed docking target"}
    best = scores.sort_values("best_affinity_kcal_mol").iloc[0]
    target_dir = Path(str(best["workdir"]))
    config_path = target_dir / "docking_config.json"
    data = json.loads(config_path.read_text(encoding="utf-8"))
    data["md_simulation"].update(
        {
            "mode": "auto" if run else "prepare",
            "top_n": 1,
            "prod_steps": 50_000_000,
            "dt_ps": 0.002,
            "temperature": 310.0,
            "pressure": 1.0,
            "timeout_seconds": int(timeout),
            "cpu": int(cpu),
            "gpu": bool(gpu),
            "mmpbsa_command": data["md_simulation"].get("mmpbsa_command"),
        }
    )
    write_json(config_path, data)
    cfg = ResolvedConfig(data, config_path)
    log = setup_logging(str(target_dir / "md.log"))
    result = run_md_simulation(cfg, log)
    write_json(docking_dir / "md_execution.json", result)
    return {
        "status": result.get("status") or ("completed" if run else "prepared"),
        "target": str(best["gene"]),
        "mode": "auto" if run else "prepare",
        "simulation_ns": 100.0,
        "temperature_k": 310.0,
        "pressure_bar": 1.0,
        "gpu": bool(gpu),
        "result": result,
    }
