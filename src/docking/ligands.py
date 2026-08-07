#!/usr/bin/env python3
"""Ligand library preparation: standardization, 3D generation and PDBQT export."""

from __future__ import annotations

import csv
import sys
import tempfile
from pathlib import Path

from .config import ResolvedConfig
from .utils import (
    DockingError,
    ToolNotFoundError,
    find_script,
    find_tool,
    run_command,
    safe_name,
    tail,
    tool_command,
    write_json,
)

MANIFEST_FIELDS = [
    "id",
    "smiles",
    "heavy_atoms",
    "rotatable_bonds",
    "pdbqt",
    "status",
    "error",
]


def prepare_ligands(cfg: ResolvedConfig, log):
    rdkit = _require_rdkit()
    Chem = rdkit.Chem
    input_path = cfg.ligand_input()
    if not input_path.exists():
        raise DockingError(f"ligand library not found: {input_path}")

    out_dir = cfg.ligand_output_dir()
    out_dir.mkdir(parents=True, exist_ok=True)
    max_ligands = cfg.get("ligand", "max_ligands") or None
    max_heavy = int(cfg.get("ligand", "max_heavy_atoms", 60))
    max_rot = int(cfg.get("ligand", "max_rotatable_bonds", 15))
    seen_ids: set[str] = set()
    rows: list[dict] = []
    ok_count = 0

    prepared_sdf = out_dir / "prepared_library.sdf"
    writer = Chem.SDWriter(str(prepared_sdf))
    try:
        for item in iter_input_molecules(input_path, cfg):
            if max_ligands and ok_count >= max_ligands:
                break
            idx = item["index"]
            base_id = safe_name(item.get("id") or f"mol_{idx}", f"mol_{idx}")
            lig_id = base_id
            suffix = 2
            while lig_id in seen_ids:
                lig_id = f"{base_id}_{suffix}"
                suffix += 1
            seen_ids.add(lig_id)

            mol = standardize_mol(item["mol"], cfg)
            if mol is None:
                rows.append(_row(lig_id, item["smiles"], "", "failed", "standardization failed"))
                continue

            heavy = mol.GetNumHeavyAtoms()
            from rdkit.Chem import rdMolDescriptors

            rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
            if heavy > max_heavy:
                rows.append(
                    _row(lig_id, item["smiles"], "", "filtered",
                         f"{heavy} heavy atoms > {max_heavy}")
                )
                continue
            if rot > max_rot:
                rows.append(
                    _row(lig_id, item["smiles"], "", "filtered",
                         f"{rot} rotatable bonds > {max_rot}")
                )
                continue

            embedded, conf_id, energy = embed_mol(mol, cfg, log)
            if embedded is None:
                rows.append(_row(lig_id, item["smiles"], "", "failed", "3D embedding failed"))
                continue

            embedded.SetProp("_Name", lig_id)
            embedded.SetProp("LigandID", lig_id)
            embedded.SetProp("SourceSMILES", item["smiles"])
            if energy is not None:
                embedded.SetProp("conformer_energy", f"{energy:.3f}")

            pdbqt_path = out_dir / f"{lig_id}.pdbqt"
            status, error = write_pdbqt(embedded, conf_id, pdbqt_path, cfg, log)
            writer.write(embedded, confId=conf_id)
            if status == "ok":
                ok_count += 1
                log.info("prepared %s (%s)", lig_id, item["smiles"])
            else:
                log.warning("PDBQT failed for %s: %s", lig_id, error)
            rows.append(
                _row(lig_id, item["smiles"], str(pdbqt_path), status, error,
                     heavy, rot)
            )
    finally:
        writer.close()

    manifest_path = out_dir / "manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as fh:
        writer_csv = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS)
        writer_csv.writeheader()
        for row in rows:
            writer_csv.writerow(row)

    summary = {
        "input": str(input_path),
        "prepared_sdf": str(prepared_sdf),
        "total": len(rows),
        "ok": ok_count,
        "failed": sum(1 for r in rows if r["status"] in ("failed", "filtered")),
        "manifest": str(manifest_path),
    }
    write_json(out_dir / "summary.json", summary)
    log.info(
        "ligand preparation complete: %s ok, %s total",
        ok_count,
        len(rows),
    )
    return summary


def _row(
    lig_id: str,
    smiles: str,
    pdbqt: str,
    status: str,
    error: str,
    heavy: int = 0,
    rot: int = 0,
) -> dict:
    return {
        "id": lig_id,
        "smiles": smiles,
        "heavy_atoms": heavy,
        "rotatable_bonds": rot,
        "pdbqt": pdbqt,
        "status": status,
        "error": error,
    }


def _require_rdkit():
    try:
        import rdkit  # type: ignore
        from rdkit import Chem, rdBase  # noqa: F401
    except ImportError as exc:
        raise DockingError(
            "RDKit is required for ligand preparation; run "
            "launchers\\install_dock_dependencies.bat first"
        ) from exc
    return rdkit


def iter_input_molecules(path: Path, cfg: ResolvedConfig):
    import pandas as pd
    from rdkit import Chem

    suffix = path.suffix.lower()
    if suffix in (".sdf", ".mol"):
        supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
        for idx, mol in enumerate(supplier, 1):
            if mol is None:
                continue
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
            yield {
                "index": idx,
                "id": name or f"mol_{idx}",
                "mol": mol,
                "smiles": Chem.MolToSmiles(mol),
            }
    elif suffix == ".smi":
        with path.open("r", encoding="utf-8") as fh:
            for idx, line in enumerate(fh, 1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                smi = parts[0]
                name = parts[1] if len(parts) > 1 else f"mol_{idx}"
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                yield {
                    "index": idx,
                    "id": name,
                    "mol": mol,
                    "smiles": smi,
                }
    elif suffix == ".csv":
        df = pd.read_csv(path)
        id_col = cfg.get("ligand", "id_column", "ID")
        smi_col = cfg.get("ligand", "smiles_column", "SMILES")
        if smi_col not in df.columns:
            raise DockingError(
                f"column '{smi_col}' not found in {path.name}; "
                "set ligand.smiles_column in the config"
            )
        for idx, (_, row) in enumerate(df.iterrows(), 1):
            smi = str(row[smi_col]).strip()
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            name = str(row.get(id_col, "")) if id_col in df.columns else ""
            yield {
                "index": idx,
                "id": name or f"mol_{idx}",
                "mol": mol,
                "smiles": smi,
            }
    else:
        raise DockingError(
            f"unsupported ligand format '{suffix}'; use .sdf, .smi or .csv"
        )


def standardize_mol(mol, cfg: ResolvedConfig):
    from rdkit import Chem
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError:
        from rdkit.Chem import rdMolStandardize

    remove_salts = cfg.get("ligand", "remove_salts", True)
    neutralize = cfg.get("ligand", "neutralize", True)
    try:
        mol = Chem.RemoveHs(mol)
        if remove_salts:
            chooser = rdMolStandardize.LargestFragmentChooser(
                preferOrganic=True
            )
            mol = chooser.choose(mol)
        if neutralize:
            mol = rdMolStandardize.Uncharger().uncharge(mol)
        return mol
    except Exception:
        return None


def embed_mol(mol, cfg: ResolvedConfig, log):
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdForceFieldHelpers

    seed = int(cfg.get("ligand", "seed", 42))
    n_conf = max(1, min(int(cfg.get("ligand", "conformers", 1)), 20))
    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.maxIterations = 1000
    conf_ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_conf, params=params)
    if not conf_ids:
        conf_ids = AllChem.EmbedMultipleConfs(
            mol_h, numConfs=1, params=AllChem.ETKDGv3()
        )
    if not conf_ids:
        return None, None, None

    energies: dict[int, float] = {}
    mmff_props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_h)
    for conf_id in conf_ids:
        mmff_ok = AllChem.MMFFOptimizeMolecule(mol_h, confId=conf_id, maxIters=1000)
        if mmff_props is not None:
            ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(
                mol_h, mmff_props, confId=conf_id
            )
            if ff is not None:
                energies[conf_id] = ff.CalcEnergy()
                continue
        if mmff_ok != 0:
            AllChem.UFFOptimizeMolecule(mol_h, confId=conf_id, maxIters=1000)
            uff = AllChem.UFFGetMoleculeForceField(mol_h, confId=conf_id)
            if uff is not None:
                energies[conf_id] = uff.CalcEnergy()

    best_id = min(conf_ids, key=lambda cid: energies.get(cid, float("inf")))
    return mol_h, best_id, energies.get(best_id)


def write_pdbqt(mol, conf_id: int, out_path: Path, cfg: ResolvedConfig, log):
    from rdkit import Chem

    engine = str(cfg.get("ligand", "engine", "auto")).lower()
    if engine == "none":
        return "skipped", "PDBQT export disabled (ligand.engine=none)"

    meeko_script = find_script("mk_prepare_ligand.py")
    obabel = find_tool("obabel")

    def _write_sdf(tmp_dir: Path) -> Path:
        tmp_sdf = tmp_dir / "ligand.sdf"
        writer = Chem.SDWriter(str(tmp_sdf))
        writer.write(mol, confId=conf_id)
        writer.close()
        return tmp_sdf

    if engine in ("auto", "meeko") and meeko_script:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                tmp_sdf = _write_sdf(Path(tmp))
                cmd = tool_command(meeko_script) + [
                    "-i",
                    str(tmp_sdf),
                    "-o",
                    str(out_path),
                ]
                result = run_command(cmd, timeout=120)
                if result.returncode != 0:
                    raise DockingError(
                        "mk_prepare_ligand.py failed: "
                        f"{tail(result.stderr or result.stdout)}"
                    )
                if out_path.exists() and out_path.stat().st_size > 0:
                    return "ok", ""
                raise DockingError("mk_prepare_ligand.py produced an empty output")
        except DockingError as exc:
            if engine == "meeko":
                return "failed", str(exc)
            log.debug("Meeko ligand prep failed, falling back to Open Babel: %s", exc)

    if engine in ("auto", "openbabel") and obabel:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                tmp_sdf = _write_sdf(Path(tmp))
                cmd = [
                    obabel,
                    str(tmp_sdf),
                    "-O",
                    str(out_path),
                    "--partialcharge",
                    "gasteiger",
                    "-p",
                    str(cfg.get("ligand", "ph", 7.4)),
                ]
                result = run_command(cmd, timeout=120)
                if result.returncode != 0:
                    return "failed", tail(result.stderr or result.stdout)
                if out_path.exists() and out_path.stat().st_size > 0:
                    return "ok", ""
                return "failed", "obabel produced an empty output"
        except DockingError as exc:
            if engine == "openbabel":
                return "failed", str(exc)

    return (
        "failed",
        "no PDBQT preparation tool available; install Meeko or Open Babel",
    )
