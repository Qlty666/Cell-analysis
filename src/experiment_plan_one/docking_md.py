"""Docking and GROMACS/MM-PBSA preparation for experiment plan one."""

from __future__ import annotations

import json
import math
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .common import (
    LOG,
    download_file,
    ensure_dir,
    get_json,
    save_figure,
    sha256_file,
    write_json,
)

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
    exhaustiveness: int = 100,
    num_modes: int = 20,
    replicates: int = 3,
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
    run_manifest: list[dict[str, Any]] = []
    ligand_id = str(compound.get("pubchem_cid") or "ligand")
    for target in targets:
        gene = str(target["gene"]).upper()
        target_dir = ensure_dir(output_dir / "targets" / gene)
        docking_run_id = f"{gene}|{ligand_id}|dock|replicate_set"
        md_run_id = f"{gene}|{ligand_id}|md|run_001"
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
                "num_modes": int(num_modes),
                "energy_range": 3.0,
                "cpu": int(cpu),
                "max_workers": 1,
                "seed": 42,
                "seeds": [],
                "replicates": int(replicates),
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
                "target_id": gene,
                "ligand_id": ligand_id,
                "run_id": md_run_id,
                "mode": "prepare",
                "top_n": 1,
                "prod_steps": 50_000_000,
                "dt_ps": 0.002,
                "temperature": 310.0,
                "pressure": 1.0,
                "protein_forcefield": "amber14sb",
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
        run_manifest.append(
            {
                "target_id": gene,
                "ligand_id": ligand_id,
                "docking_run_id": docking_run_id,
                "md_run_id": md_run_id,
                "workdir": str(target_dir),
                "config": str(config_path),
                "config_sha256": sha256_file(config_path),
                "receptor_sha256": sha256_file(receptor_copy),
                "ligand_sha256": sha256_file(ligand),
                "random_seed": 42,
                "docking_replicates": int(replicates),
                "structure": structure,
            }
        )
        cfg = ResolvedConfig(data, config_path)
        log = setup_logging(str(target_dir / "docking.log"))
        try:
            prepare_receptor(cfg, log)
            prepare_ligands(cfg, log)
            run_docking(cfg, log)
            analyze_results(cfg, log)
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
    write_json(
        output_dir / "docking_run_manifest.json",
        {
            "schema_version": 1,
            "targets": run_manifest,
            "selection_rule": (
                "prepare MD only for the first pre-specified target in the "
                "priority order; Vina scores are not compared across proteins"
            ),
        },
    )
    if not score_frame.empty:
        _plot_docking_heatmap(score_frame, output_dir / "fig5c_docking_affinity_heatmap.png")
        # Vina scores are not comparable across different protein pockets.
        best_target = score_frame.iloc[0]
        best_dir = Path(str(best_target["workdir"]))
        _write_docking_pose_figure(
            best_dir,
            output_dir / "fig5a_docking_pose_3d.png",
            str(best_target["gene"]),
        )
        interaction_status = _write_interaction_figure(
            best_dir,
            output_dir / "fig5b_interaction_schematic.png",
            str(best_target["gene"]),
        )
        write_json(
            output_dir / "fig5b_status.json",
            {
                **interaction_status,
                "limitations": (
                    "PLIP interactions describe the selected docking pose; "
                    "they are not experimental binding evidence."
                    if interaction_status.get("validated")
                    else (
                        "Aromatic plane, protonation and charge evidence were "
                        "not validated with PLIP; candidate labels must not "
                        "be reported as confirmed interactions."
                    )
                ),
            },
        )
    write_json(
        output_dir / "docking_summary.json",
        {
            "targets_requested": len(targets),
            "targets_completed": len(rows),
            "targets": target_records,
            "exhaustiveness": exhaustiveness,
            "seed_replicates": int(replicates),
            "run_manifest": str(output_dir / "docking_run_manifest.json"),
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
    frame = scores.dropna(subset=["best_affinity_kcal_mol"]).copy()
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
    values = frame["best_affinity_kcal_mol"].to_numpy(dtype=float)
    value_min = float(np.nanmin(values))
    value_max = float(np.nanmax(values))
    for row, value in enumerate(values):
        normalized = (
            (float(value) - value_min) / (value_max - value_min)
            if value_max > value_min
            else 0.5
        )
        rgba = plt.get_cmap("YlGnBu_r")(normalized)
        luminance = (
            0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
        )
        ax.text(
            0,
            row,
            f"{value:.2f}",
            ha="center",
            va="center",
            color="#111111" if luminance > 0.55 else "white",
            fontsize=7,
            fontweight="bold",
        )
    ax.set_title("6PPD-Q docking affinity by pre-specified target order", fontweight="bold")
    ax.text(
        0,
        -0.24,
        "Vina scores are pocket-specific and are not ranked across proteins.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        color="#4a5560",
    )
    fig.colorbar(image, ax=ax, label="kcal/mol", shrink=0.7)
    save_figure(fig, output)


def _find_plip() -> str | None:
    discovered = (
        shutil.which("plip")
        or shutil.which("plip.exe")
    )
    if discovered:
        return discovered
    executable = "plip.exe" if sys.platform.startswith("win") else "plip"
    candidates = [
        Path(sys.executable).parent / executable,
        Path(sys.executable).parent / "Scripts" / executable,
        Path(sys.executable).parent / "bin" / executable,
    ]
    local = next((path for path in candidates if path.exists()), None)
    return str(local) if local else None


def _write_ligand_pdb(
    receptor: Path,
    ligand_atoms: list[dict[str, Any]],
    output: Path,
) -> None:
    lines: list[str] = []
    with receptor.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM"):
                lines.append(line.rstrip("\n"))
    for serial, atom in enumerate(ligand_atoms, start=1):
        coordinate = np.asarray(atom["coordinate"], dtype=float)
        element = str(atom.get("element") or "C").upper()[:1]
        atom_name = re.sub(r"[^A-Za-z0-9]", "", str(atom.get("atom") or element))
        if not atom_name:
            atom_name = f"{element}{serial}"
        lines.append(
            "HETATM"
            f"{serial:5d} {atom_name[:4]:<4} LIG Z{1:4d}    "
            f"{coordinate[0]:8.3f}{coordinate[1]:8.3f}{coordinate[2]:8.3f}"
            f"  1.00  0.00          {element:>2}"
        )
    lines.append("END")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _plip_interactions(
    receptor: Path,
    ligand_atoms: list[dict[str, Any]],
    output_dir: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    executable = _find_plip()
    if not executable:
        return pd.DataFrame(), {
            "status": "blocked",
            "reason": "PLIP executable is not installed in this environment",
        }
    work_dir = ensure_dir(output_dir / "plip")
    complex_path = work_dir / "complex.pdb"
    _write_ligand_pdb(receptor, ligand_atoms, complex_path)
    report_dir = ensure_dir(work_dir / "report")
    process = subprocess.run(
        [
            executable,
            "-f",
            str(complex_path),
            "-o",
            str(report_dir),
            "-x",
            "-q",
            "--nofix",
            "--maxthreads",
            "1",
            "--name",
            "plip_report",
        ],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=600,
    )
    xml_path = report_dir / "plip_report.xml"
    if process.returncode != 0 or not xml_path.exists():
        detail = (process.stderr or process.stdout or "").strip()
        return pd.DataFrame(), {
            "status": "failed",
            "reason": detail[-1200:] or "PLIP did not produce an XML report",
            "complex": str(complex_path),
        }
    from plip.exchange.xml import PlipXML

    report = PlipXML(str(xml_path))
    type_map = {
        "hydrogen_bond": (
            "hbonds",
            lambda item: getattr(item, "dist", np.nan),
            lambda item: getattr(item, "don_angle", np.nan),
        ),
        "hydrophobic_contact": (
            "hydrophobics",
            lambda item: getattr(item, "dist", np.nan),
            lambda item: np.nan,
        ),
        "water_bridge": (
            "wbridges",
            lambda item: getattr(item, "dist", np.nan),
            lambda item: getattr(item, "don_angle", np.nan),
        ),
        "salt_bridge": (
            "sbridges",
            lambda item: getattr(item, "dist", np.nan),
            lambda item: np.nan,
        ),
        "pi_stacking": (
            "pi_stacks",
            lambda item: getattr(item, "dist", np.nan),
            lambda item: getattr(item, "angle", np.nan),
        ),
        "pi_cation": (
            "pi_cations",
            lambda item: getattr(item, "dist", np.nan),
            lambda item: np.nan,
        ),
        "halogen_bond": (
            "halogens",
            lambda item: getattr(item, "dist", np.nan),
            lambda item: getattr(item, "don_angle", np.nan),
        ),
        "metal_complex": (
            "metal_complexes",
            lambda item: getattr(item, "dist", np.nan),
            lambda item: np.nan,
        ),
    }
    records: list[dict[str, Any]] = []
    sites = list(report.bsites.values())
    for interaction_type, (attribute, distance, angle) in type_map.items():
        items = [
            item
            for site in sites
            for item in getattr(site, attribute, [])
        ]
        for item in items:
            records.append(
                {
                    "interaction_type": interaction_type,
                    "ligand_atom": str(
                        getattr(item, "ligcarbonidx", "")
                        or getattr(item, "donoridx", "")
                        or getattr(item, "lig_idx_list", "")
                    ),
                    "residue": str(getattr(item, "restype", "")),
                    "residue_number": int(getattr(item, "resnr", 0) or 0),
                    "chain": str(getattr(item, "reschain", "")),
                    "distance_angstrom": distance(item),
                    "angle_degrees": angle(item),
                    "evidence_class": "PLIP_3.0.1",
                    "evidence_note": "validated by PLIP chemistry rules",
                }
            )
    frame = pd.DataFrame(records)
    return frame, {
        "status": "completed" if not frame.empty else "valid_negative",
        "method": "PLIP 3.0.1",
        "complex": str(complex_path),
        "xml_report": str(xml_path),
        "n_interactions": int(len(frame)),
    }


def _write_interaction_figure(
    target_dir: Path,
    output: Path,
    gene: str,
) -> dict[str, Any]:
    receptor = target_dir / "data" / "receptors" / "receptor.pdb"
    pose_table = target_dir / "outputs" / "run_001" / "docked" / "results.csv"
    if not receptor.exists() or not pose_table.exists():
        return {
            "status": "not_run",
            "reason": "receptor or docking result table is missing",
        }
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
        return {
            "status": "not_run",
            "reason": "top docking pose file is missing",
        }
    ligand_atoms = _parse_pdbqt_atoms(pose_path)
    contacts, plip_status = _plip_interactions(
        receptor,
        ligand_atoms,
        target_dir / "outputs" / "interaction_validation",
    )
    validated = not contacts.empty
    plip_ran = str(plip_status.get("status") or "") in {
        "completed",
        "valid_negative",
    }
    if not validated and not plip_ran:
        contacts = _typed_interactions(receptor, ligand_atoms)
    contacts.to_csv(output.with_suffix(".csv"), index=False)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
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
        ax.text(
            0.5,
            0.42,
            (
                "No PLIP-validated interaction"
                if plip_ran
                else "No contact within 4.5 A"
            ),
            ha="center",
            transform=ax.transAxes,
        )
    else:
        top = (
            contacts.sort_values("distance_angstrom")
            .groupby("interaction_type", group_keys=False)
            .head(5)
            .sort_values(
                [
                    "interaction_type",
                    "distance_angstrom",
                ]
            )
            .head(18)
        )
        type_colors = {
            "hydrogen_bond": "#2f6bb3",
            "hydrophobic_contact": "#d29b32",
            "water_bridge": "#2b9eb3",
            "salt_bridge": "#c0392b",
            "pi_stacking": "#6c5fa7",
            "pi_cation": "#8e6bb3",
            "halogen_bond": "#7b6d5d",
            "metal_complex": "#6d7f3c",
            "hydrogen_bond_geometric_candidate": "#2f6bb3",
            "polar_contact_candidate": "#6b8eb5",
            "hydrophobic_contact": "#d29b32",
            "salt_bridge_candidate": "#c0392b",
            "aromatic_contact_candidate": "#6c5fa7",
            "vdw_contact": "#7b8794",
        }
        for position, (_, row) in enumerate(top.iterrows()):
            angle = 2 * math.pi * position / max(len(top), 1)
            x = 0.5 + 0.38 * math.cos(angle)
            y = 0.5 + 0.38 * math.sin(angle)
            interaction = str(row["interaction_type"])
            ax.annotate(
                f"{row['residue']}{row['residue_number']}\n{interaction}",
                xy=(0.5, 0.5),
                xytext=(x, y),
                ha="center",
                va="center",
                fontsize=7.3,
                arrowprops={
                    "arrowstyle": "-",
                    "color": type_colors.get(interaction, "#7b8794"),
                    "lw": 1.35,
                },
                bbox={
                    "boxstyle": "round,pad=0.3",
                    "facecolor": "white",
                    "edgecolor": type_colors.get(interaction, "#a7b1bb"),
                },
                zorder=6,
                xycoords=ax.transAxes,
                textcoords=ax.transAxes,
            )
        handles = [
            plt.Line2D(
                [0],
                [0],
                color=color,
                lw=2,
                label=interaction.replace("_", " "),
            )
            for interaction, color in type_colors.items()
            if interaction in set(contacts["interaction_type"])
        ]
        if handles:
            ax.legend(
                handles=handles,
                loc="lower left",
                frameon=False,
                fontsize=7.2,
            )
    ax.text(
        0.5,
        0.02,
        (
            "Interactions validated by PLIP 3.0.1."
            if validated
            else "PLIP 3.0.1 completed without a validated interaction."
            if plip_ran
            else "Interaction labels are distance/geometry candidates; PLIP "
            "did not provide a validated interaction table."
        ),
        ha="center",
        va="bottom",
        fontsize=7.8,
        color="#5d6670",
        transform=ax.transAxes,
    )
    save_figure(fig, output)
    return {
        **plip_status,
        "figure": str(output),
        "machine_readable_table": str(output.with_suffix(".csv")),
        "validated": validated,
    }


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
    ligand_atoms = _parse_pdbqt_atoms(pose_path)
    if not ligand_atoms:
        return
    ligand = np.asarray([atom["coordinate"] for atom in ligand_atoms])
    protein_atoms: list[dict[str, Any]] = []
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
            if np.linalg.norm(np.asarray(coordinate) - ligand.mean(axis=0)) <= 12.0:
                protein_atoms.append(
                    {
                        "coordinate": np.asarray(coordinate),
                        "atom": line[12:16].strip(),
                        "residue": line[17:20].strip(),
                        "residue_number": int(line[22:26]),
                        "chain": line[21:22].strip(),
                        "element": (
                            line[76:78].strip().upper()
                            or re.sub(r"[^A-Za-z]", "", line[12:16])[:1].upper()
                        ),
                    }
                )
    protein = np.asarray([atom["coordinate"] for atom in protein_atoms])
    fig = plt.figure(figsize=(6.6, 5.6))
    ax = fig.add_subplot(111, projection="3d")
    if protein.size:
        ax.scatter(
            protein[:, 0],
            protein[:, 1],
            protein[:, 2],
            s=10,
            c="#7b8794",
            alpha=0.30,
            depthshade=False,
            label="Protein binding-site atoms",
        )
    ca_trace: dict[str, list[np.ndarray]] = {}
    for atom in protein_atoms:
        if atom["atom"] == "CA":
            ca_trace.setdefault(str(atom["chain"]), []).append(atom["coordinate"])
    for trace in ca_trace.values():
        if len(trace) >= 2:
            points = np.asarray(trace)
            ax.plot(
                points[:, 0],
                points[:, 1],
                points[:, 2],
                color="#5e6b78",
                alpha=0.55,
                linewidth=1.0,
            )
    element_colors = {
        "C": "#d3544f",
        "N": "#4169a8",
        "O": "#c0392b",
        "S": "#d6a62f",
        "H": "#d9dde1",
        "F": "#4e9b6e",
        "CL": "#4e9b6e",
    }
    for atom in ligand_atoms:
        coordinate = atom["coordinate"]
        element = str(atom["element"]).upper()
        ax.scatter(
            [coordinate[0]],
            [coordinate[1]],
            [coordinate[2]],
            s=42,
            c=element_colors.get(element, "#cf4f45"),
            depthshade=False,
            edgecolors="white",
            linewidths=0.4,
        )
    for left in range(len(ligand)):
        for right in range(left + 1, len(ligand)):
            distance = float(np.linalg.norm(ligand[left] - ligand[right]))
            if distance <= 1.9:
                ax.plot(
                    ligand[[left, right], 0],
                    ligand[[left, right], 1],
                    ligand[[left, right], 2],
                    color="#7c3f3a",
                    linewidth=1.4,
                )
    pocket = _receptor_contacts(receptor, ligand, cutoff=4.5)
    if not pocket.empty:
        labels = []
        for row in pocket.head(3).itertuples(index=False):
            labels.append(f"{row.residue}{row.residue_number}")
            match = next(
                (
                    atom
                    for atom in protein_atoms
                    if atom["residue"] == row.residue
                    and int(atom["residue_number"]) == int(row.residue_number)
                ),
                None,
            )
            if match is not None:
                coordinate = match["coordinate"]
                ax.text(
                    coordinate[0],
                    coordinate[1],
                    coordinate[2],
                    f"{row.residue}{row.residue_number}",
        fontsize=7,
                    color="#39424c",
                )
        fig.text(
            0.02,
            0.02,
            "Pocket residues: " + ", ".join(labels),
            fontsize=7.5,
            color="#5d6670",
        )
    ligand_handle = plt.Line2D(
        [0],
        [0],
        marker="o",
        color="none",
        markerfacecolor="#cf4f45",
        markersize=6,
        label="6PPD-Q",
    )
    protein_handle = plt.Line2D(
        [0],
        [0],
        color="#5e6b78",
        linewidth=1.4,
        label="Protein backbone trace",
    )
    ax.set_xlabel("x (Angstrom)", labelpad=2)
    ax.set_ylabel("y (Angstrom)", labelpad=2)
    ax.set_zlabel("z (Angstrom)", labelpad=2)
    ax.set_title(f"{gene}-6PPD-Q best docking pose", fontweight="bold", pad=8)
    ax.set_box_aspect((1.0, 1.0, 0.82))
    ax.grid(False)
    if ligand.size:
        center = ligand.mean(axis=0)
        radius = max(4.5, float(np.ptp(ligand, axis=0).max()) * 0.9)
        if protein.size:
            distances = np.linalg.norm(protein - center, axis=1)
            radius = max(radius, float(np.percentile(distances, 92)) + 1.5)
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)
        ax.set_zlim(center[2] - radius, center[2] + radius)
    ax.legend(
        handles=[protein_handle, ligand_handle],
        loc="upper right",
        frameon=False,
        fontsize=7,
    )
    ax.view_init(elev=18, azim=42)
    save_figure(fig, output)


def _parse_pdbqt_coordinates(path: Path) -> np.ndarray:
    atoms = _parse_pdbqt_atoms(path)
    return np.asarray([atom["coordinate"] for atom in atoms])


def _parse_pdbqt_atoms(path: Path) -> list[dict[str, Any]]:
    atoms: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                coordinate = np.asarray(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ]
                )
            except ValueError:
                continue
            atom_type = re.sub(r"[^A-Za-z]", "", line[76:79]).upper()
            atoms.append(
                {
                    "coordinate": coordinate,
                    "atom": line[12:16].strip(),
                    "element": atom_type[:1] or line[12:16].strip()[:1].upper(),
                    "atom_type": atom_type,
                }
            )
    return atoms


def _typed_interactions(
    receptor: Path,
    ligand_atoms: list[dict[str, Any]],
    cutoff: float = 4.5,
) -> pd.DataFrame:
    if not ligand_atoms:
        return pd.DataFrame(
            columns=[
                "interaction_type",
                "ligand_atom",
                "residue",
                "residue_number",
                "chain",
                "distance_angstrom",
            ]
        )
    protein_atoms: list[dict[str, Any]] = []
    with receptor.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            try:
                coordinate = np.asarray(
                    [
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ]
                )
                residue_number = int(line[22:26])
            except ValueError:
                continue
            protein_atoms.append(
                {
                    "coordinate": coordinate,
                    "atom": line[12:16].strip().upper(),
                    "residue": line[17:20].strip().upper(),
                    "residue_number": residue_number,
                    "chain": line[21:22].strip(),
                    "element": (
                        line[76:78].strip().upper()
                        or re.sub(r"[^A-Za-z]", "", line[12:16])[:1].upper()
                    ),
                }
            )
    if not protein_atoms:
        return pd.DataFrame()
    ligand_coordinates = np.asarray(
        [atom["coordinate"] for atom in ligand_atoms]
    )
    tree = cKDTree(ligand_coordinates)
    records: dict[
        tuple[str, str, str, int, str],
        dict[str, Any],
    ] = {}
    aromatic_residues = {"PHE", "TYR", "TRP", "HIS"}
    negative_atoms = {
        ("ASP", "OD1"),
        ("ASP", "OD2"),
        ("GLU", "OE1"),
        ("GLU", "OE2"),
    }
    positive_atoms = {
        ("LYS", "NZ"),
        ("ARG", "NE"),
        ("ARG", "NH1"),
        ("ARG", "NH2"),
        ("HIS", "ND1"),
        ("HIS", "NE2"),
    }
    def donor_angle(
        donor: dict[str, Any],
        acceptor_coordinate: np.ndarray,
    ) -> float | None:
        donor_coordinate = np.asarray(donor["coordinate"], dtype=float)
        hydrogens = [
            atom
            for atom in ligand_atoms
            if str(atom.get("element")).upper() == "H"
            and float(np.linalg.norm(atom["coordinate"] - donor_coordinate)) <= 1.25
        ]
        if not hydrogens:
            return None
        hydrogen = min(
            hydrogens,
            key=lambda atom: float(
                np.linalg.norm(atom["coordinate"] - donor_coordinate)
            ),
        )
        vector_a = donor_coordinate - hydrogen["coordinate"]
        vector_b = acceptor_coordinate - hydrogen["coordinate"]
        denominator = float(np.linalg.norm(vector_a) * np.linalg.norm(vector_b))
        if denominator <= 0:
            return None
        cosine = float(np.dot(vector_a, vector_b) / denominator)
        return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))

    for protein in protein_atoms:
        distance, ligand_index = tree.query(protein["coordinate"])
        distance = float(distance)
        if distance > cutoff:
            continue
        ligand = ligand_atoms[int(ligand_index)]
        ligand_element = str(ligand["element"]).upper()
        protein_element = str(protein["element"]).upper()
        interactions: list[tuple[str, str, float | None, str]] = []
        if (
            protein_element in {"N", "O", "S"}
            and ligand_element in {"N", "O", "S"}
            and distance <= 3.5
        ):
            ligand_donor = ligand
            protein_donor = protein
            angle = donor_angle(
                ligand_donor,
                np.asarray(protein["coordinate"], dtype=float),
            )
            protein_hydrogen_candidates = [
                atom
                for atom in protein_atoms
                if str(atom.get("element")).upper() == "H"
                and float(
                    np.linalg.norm(
                        atom["coordinate"]
                        - np.asarray(protein_donor["coordinate"], dtype=float)
                    )
                )
                <= 1.25
            ]
            if angle is None and protein_hydrogen_candidates:
                hydrogen = min(
                    protein_hydrogen_candidates,
                    key=lambda atom: float(
                        np.linalg.norm(
                            atom["coordinate"]
                            - np.asarray(protein_donor["coordinate"], dtype=float)
                        )
                    ),
                )
                vector_a = (
                    np.asarray(protein_donor["coordinate"], dtype=float)
                    - hydrogen["coordinate"]
                )
                vector_b = (
                    np.asarray(ligand["coordinate"], dtype=float)
                    - hydrogen["coordinate"]
                )
                denominator = float(
                    np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
                )
                if denominator > 0:
                    angle = float(
                        np.degrees(
                            np.arccos(
                                np.clip(
                                    float(
                                        np.dot(vector_a, vector_b) / denominator
                                    ),
                                    -1.0,
                                    1.0,
                                )
                            )
                        )
                    )
            if angle is not None and 120.0 <= angle <= 180.0:
                interactions.append(
                    (
                        "hydrogen_bond_geometric_candidate",
                        "geometry_and_distance",
                        angle,
                        "donor-H-acceptor angle passes the geometric candidate rule",
                    )
                )
            else:
                interactions.append(
                    (
                        "polar_contact_candidate",
                        "distance_only",
                        angle,
                        "polar atoms are close; explicit hydrogen-bond geometry was not verifiable",
                    )
                )
        if (
            (protein["residue"], protein["atom"]) in negative_atoms
            and ligand_element in {"N", "O"}
            and distance <= 4.0
        ) or (
            (protein["residue"], protein["atom"]) in positive_atoms
            and ligand_element in {"N", "O"}
            and distance <= 4.0
        ):
            interactions.append(
                (
                    "salt_bridge_candidate",
                    "charged_group_distance",
                    None,
                    "charged functional groups are close; not experimentally confirmed",
                )
            )
        if (
            protein["residue"] in aromatic_residues
            and ligand_element in {"C", "N", "O"}
            and distance <= 4.5
        ):
            interactions.append(
                (
                    "aromatic_contact_candidate",
                    "ring_identity_and_distance",
                    None,
                    "aromatic residue and ligand atom are close; plane geometry is not claimed",
                )
            )
        if (
            protein_element == "C"
            and ligand_element == "C"
            and distance <= 4.5
        ):
            interactions.append(
                (
                    "hydrophobic_contact",
                    "carbon_distance",
                    None,
                    "carbon-carbon close contact",
                )
            )
        if not interactions:
            interactions.append(
                (
                    "vdw_contact",
                    "distance_only",
                    None,
                    "no stronger typed contact candidate was detected",
                )
            )
        for interaction, evidence, angle, note in interactions:
            key = (
                interaction,
                str(ligand["atom"]),
                str(protein["residue"]),
                int(protein["residue_number"]),
                str(protein["chain"]),
            )
            current = records.get(key)
            if current is None or distance < float(current["distance_angstrom"]):
                records[key] = {
                    "distance_angstrom": distance,
                    "evidence_class": evidence,
                    "angle_degrees": (
                        round(float(angle), 3)
                        if angle is not None and np.isfinite(angle)
                        else np.nan
                    ),
                    "evidence_note": note,
                }
    rows = [
        {
            "interaction_type": key[0],
            "ligand_atom": key[1],
            "residue": key[2],
            "residue_number": key[3],
            "chain": key[4],
            **value,
        }
        for key, value in records.items()
    ]
    return (
        pd.DataFrame(rows)
        .sort_values(
            ["interaction_type", "distance_angstrom", "residue_number"]
        )
        .reset_index(drop=True)
        if rows
        else pd.DataFrame()
    )


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
    protein_forcefield: str = "amber14sb",
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
    best = scores.iloc[0]
    gene = str(best["gene"]).upper()
    manifest_path = docking_dir / "docking_run_manifest.json"
    manifest = (
        json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest_path.exists()
        else {}
    )
    manifest_entry = next(
        (
            entry
            for entry in (manifest.get("targets") or [])
            if str(entry.get("target_id") or "").upper() == gene
        ),
        {},
    )
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
            "protein_forcefield": str(protein_forcefield),
            "timeout_seconds": int(timeout),
            "cpu": int(cpu),
            "gpu": bool(gpu),
            "mmpbsa_command": data["md_simulation"].get("mmpbsa_command"),
            "target_id": gene,
            "ligand_id": str(
                manifest_entry.get("ligand_id")
                or data["md_simulation"].get("ligand_id")
                or ""
            ),
            "run_id": str(
                manifest_entry.get("md_run_id")
                or data["md_simulation"].get("run_id")
                or f"{gene}|md|run_001"
            ),
        }
    )
    write_json(config_path, data)
    cfg = ResolvedConfig(data, config_path)
    log = setup_logging(str(target_dir / "md.log"))
    result = run_md_simulation(cfg, log)
    write_json(docking_dir / "md_execution.json", result)
    run_dir = _select_manifest_run_dir(
        cfg.md_dir(),
        ligand_id=str(
            manifest_entry.get("ligand_id")
            or data["md_simulation"].get("ligand_id")
            or ""
        ),
    )
    production_ns = float(result.get("prod_ns") or 0.0)
    if str(result.get("mode")) == "prepare":
        run_status = "prepared"
    elif int(result.get("completed") or 0) > 0:
        run_status = "completed"
    else:
        run_status = "failed"
    if run_status == "completed" and production_ns < 100.0 - 1e-9:
        run_status = "incomplete_duration"
    write_json(
        docking_dir / "md_run_manifest.json",
        {
            "schema_version": 1,
            "target_id": gene,
            "ligand_id": str(manifest_entry.get("ligand_id") or ""),
            "run_id": str(
                manifest_entry.get("md_run_id")
                or data["md_simulation"].get("run_id")
                or f"{gene}|md|run_001"
            ),
            "run_dir": str(run_dir) if run_dir else "",
            "status": run_status,
            "production_ns": production_ns,
            "configured_production_ns": 100.0,
            "random_seed": int(data["md_simulation"].get("gen_seed", 42)),
            "receptor_sha256": sha256_file(
                Path(data["receptor"]["input"])
            ),
            "ligand_sha256": sha256_file(
                Path(data["ligand"]["input"])
            ),
            "selection_rule": (
                "first pre-specified target in docking_scores.csv; "
                "cross-protein Vina scores are not ranked"
            ),
            "result": result,
        },
    )
    return {
        "status": result.get("status") or ("completed" if run else "prepared"),
        "target": str(best["gene"]),
        "mode": "auto" if run else "prepare",
        "simulation_ns": 100.0,
        "production_ns": production_ns,
        "run_dir": str(run_dir) if run_dir else "",
        "temperature_k": 310.0,
        "pressure_bar": 1.0,
        "gpu": bool(gpu),
        "result": result,
    }


def _select_manifest_run_dir(md_dir: Path, *, ligand_id: str) -> Path | None:
    if not md_dir.exists():
        return None
    candidates = sorted(path for path in md_dir.iterdir() if path.is_dir())
    exact = [path for path in candidates if path.name == ligand_id]
    if len(exact) == 1:
        return exact[0]
    return candidates[0] if len(candidates) == 1 else None
