"""Optional PoseBusters validation for SDF docking poses."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

from docking.utils import write_json


def _convert_pdbqt_to_sdf(
    pdbqt_path: Path,
    sdf_dir: Path,
    gene: str,
) -> tuple[Path | None, str]:
    try:
        from meeko import PDBQTMolecule
        from meeko.rdkit_mol_create import RDKitMolCreate
    except ImportError:
        return None, "Meeko is not installed"
    try:
        molecule = PDBQTMolecule.from_file(
            str(pdbqt_path),
            skip_typing=False,
        )
        sdf_text, failures = RDKitMolCreate.write_sd_string(molecule)
    except Exception as exc:  # noqa: BLE001
        return None, f"PDBQT conversion failed: {exc}"
    if not sdf_text:
        return None, "PDBQT conversion produced no poses"
    gene_dir = sdf_dir / gene
    gene_dir.mkdir(parents=True, exist_ok=True)
    output = gene_dir / f"{pdbqt_path.stem}.sdf"
    output.write_text(sdf_text, encoding="utf-8")
    note = f"{len(failures)} pose(s) failed conversion" if failures else ""
    return output, note


def _find_receptor(pose_path: Path) -> Path | None:
    for ancestor in pose_path.parents:
        receptor_dir = ancestor / "data" / "receptors"
        if receptor_dir.exists():
            receptors = sorted(receptor_dir.glob("*.pdb"))
            if receptors:
                return receptors[0]
    return None


def run_pose_qc(workdir: Path, out_dir: Path) -> dict:
    """Validate SDF poses when PoseBusters and receptor structures are available."""
    out_dir.mkdir(parents=True, exist_ok=True)
    work_root = workdir / "work"
    sdf_files = sorted(work_root.glob("*/**/*.sdf"))
    pdbqt_files = sorted(work_root.glob("*/**/*.pdbqt"))
    results_path = out_dir / "pose_qc_results.csv"
    empty_columns = [
        "gene",
        "pose_file",
        "pose_index",
        "pb_valid",
        "source",
        "pose_format",
        "source_pose_file",
        "receptor_file",
    ]

    def gene_for(path: Path) -> str:
        try:
            relative = path.relative_to(work_root)
            return relative.parts[0] if relative.parts else "unknown"
        except ValueError:
            return path.parents[1].name if len(path.parents) > 1 else "unknown"

    # Keep the original pose and its receptor together.  A converted PDBQT
    # file lives under outputs/integration, so looking for data/receptors
    # from the converted SDF itself would silently lose the receptor.
    pose_items: list[dict[str, object]] = []
    conversion_notes: list[str] = []
    if sdf_files:
        for sdf in sdf_files:
            pose_items.append(
                {
                    "gene": gene_for(sdf),
                    "pose_file": sdf,
                    "source_pose_file": sdf,
                    "receptor": _find_receptor(sdf),
                    "pose_format": "sdf",
                }
            )
    else:
        for pdbqt in pdbqt_files:
            gene = gene_for(pdbqt)
            receptor = _find_receptor(pdbqt)
            sdf, note = _convert_pdbqt_to_sdf(
                pdbqt,
                out_dir / "poses_sdf",
                gene,
            )
            if sdf is not None:
                pose_items.append(
                    {
                        "gene": gene,
                        "pose_file": sdf,
                        "source_pose_file": pdbqt,
                        "receptor": receptor,
                        "pose_format": "pdbqt_converted",
                    }
                )
            if note:
                conversion_notes.append(f"{pdbqt.name}: {note}")

    if not pose_items:
        if not sdf_files and not pdbqt_files:
            summary = {
                "status": "skipped",
                "reason": "no SDF or PDBQT pose files found",
                "sdf_files": 0,
                "poses": 0,
                "pb_valid": 0,
                "pb_valid_rate": None,
                "gate_passed": False,
            }
        else:
            summary = {
                "status": "unavailable",
                "reason": (
                    "PDBQT conversion unavailable: "
                    + ("; ".join(conversion_notes) or "Meeko not available")
                ),
                "sdf_files": 0,
                "poses": 0,
                "pb_valid": 0,
                "pb_valid_rate": None,
                "gate_passed": False,
            }
        pd.DataFrame(columns=empty_columns).to_csv(results_path, index=False)
        write_json(out_dir / "pose_qc_summary.json", summary)
        return summary

    if importlib.util.find_spec("posebusters") is None:
        summary = {
            "status": "unavailable",
            "reason": "PoseBusters is not installed",
            "sdf_files": len(pose_items),
            "pose_format": (
                str(pose_items[0]["pose_format"])
                if len({item["pose_format"] for item in pose_items}) == 1
                else "mixed"
            ),
            "poses": 0,
            "pb_valid": 0,
            "pb_valid_rate": None,
            "gate_passed": False,
        }
        pd.DataFrame(columns=empty_columns).to_csv(results_path, index=False)
        write_json(out_dir / "pose_qc_summary.json", summary)
        return summary

    from posebusters import PoseBusters

    bust = PoseBusters(config="dock")
    rows: list[dict] = []
    for item in pose_items:
        gene = str(item["gene"])
        sdf = Path(item["pose_file"])
        receptor = item["receptor"]
        pose_format = str(item["pose_format"])
        source_pose_file = Path(item["source_pose_file"])
        if receptor is None:
            rows.append(
                {
                    "gene": gene,
                    "pose_file": str(sdf),
                    "pose_index": -1,
                    "pb_valid": False,
                    "source": "no_receptor",
                    "pose_format": pose_format,
                    "source_pose_file": str(source_pose_file),
                    "receptor_file": "",
                }
            )
            continue
        try:
            frame = bust.bust(
                mol_pred=str(sdf),
                mol_cond=str(Path(receptor)),
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "gene": gene,
                    "pose_file": str(sdf),
                    "pose_index": -1,
                    "pb_valid": False,
                    "source": f"error:{exc}",
                    "pose_format": pose_format,
                    "source_pose_file": str(source_pose_file),
                    "receptor_file": str(receptor),
                }
            )
            continue
        bool_columns = [
            column
            for column in frame.select_dtypes(include="bool").columns
            if not column.lower().startswith("rmsd")
        ]
        for index, row in frame.iterrows():
            valid = bool(row[bool_columns].all()) if bool_columns else False
            rows.append(
                {
                    "gene": gene,
                    "pose_file": str(sdf),
                    "pose_index": int(index),
                    "pb_valid": valid,
                    "source": "posebusters",
                    "pose_format": pose_format,
                    "source_pose_file": str(source_pose_file),
                    "receptor_file": str(receptor),
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(results_path, index=False)
    valid_count = int(result["pb_valid"].sum()) if "pb_valid" in result else 0
    total = int(len(result))
    rate = valid_count / total if total else 0.0
    summary = {
        "status": "completed",
        "reason": "",
        "sdf_files": len(pose_items),
        "pose_format": (
            str(pose_items[0]["pose_format"])
            if len({item["pose_format"] for item in pose_items}) == 1
            else "mixed"
        ),
        "poses": total,
        "pb_valid": valid_count,
        "pb_valid_rate": rate,
        "gate_passed": bool(total > 0 and rate >= 0.50),
    }
    write_json(out_dir / "pose_qc_summary.json", summary)
    return summary
