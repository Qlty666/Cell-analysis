"""Final Figure-5 trajectory panels parsed from GROMACS/MM-PBSA outputs."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .common import LOG, ensure_dir, save_figure, write_json


def _parse_xvg(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "@")):
            continue
        try:
            rows.append([float(value) for value in line.split()])
        except ValueError:
            continue
    return np.asarray(rows, dtype=float) if rows else None


def _latest_run_dir(docking_dir: Path) -> Path | None:
    roots = [
        path
        for path in docking_dir.rglob("06_md")
        if path.is_dir()
    ]
    if not roots:
        return None
    root = max(roots, key=lambda path: path.stat().st_mtime)
    runs = [path for path in root.iterdir() if path.is_dir()]
    if not runs:
        return None
    return max(runs, key=lambda path: path.stat().st_mtime)


def _metrics(docking_dir: Path) -> pd.DataFrame:
    candidates = sorted(docking_dir.rglob("md_simulation_results.csv"))
    if not candidates:
        return pd.DataFrame()
    return pd.read_csv(candidates[-1])


def _empty_panel(output: Path, title: str, reason: str) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.axis("off")
    ax.text(
        0.5,
        0.62,
        title,
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.40,
        "Not available",
        ha="center",
        va="center",
        fontsize=11,
        color="#b24c3c",
        fontweight="bold",
    )
    ax.text(
        0.5,
        0.22,
        reason,
        ha="center",
        va="center",
        fontsize=8.5,
        color="#5d6670",
    )
    save_figure(fig, output)


def _time_series_panel(
    data: np.ndarray | None,
    output: Path,
    *,
    title: str,
    ylabel: str,
    color: str,
    scale: float = 1.0,
) -> bool:
    import matplotlib.pyplot as plt

    if data is None or not len(data):
        _empty_panel(
            output,
            title,
            "The corresponding GROMACS XVG file was not found.",
        )
        return False
    values = data[:, 1] * scale
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(data[:, 0], values, color=color, linewidth=1.35)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight="bold")
    ax.grid(alpha=0.22)
    mean = float(np.nanmean(values))
    std = float(np.nanstd(values, ddof=1)) if len(values) > 1 else 0.0
    ax.text(
        0.02,
        0.96,
        f"mean={mean:.3f}; SD={std:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    save_figure(fig, output)
    return True


def _rmsf_panel(
    data: np.ndarray | None,
    output: Path,
    binding_site: Path,
) -> bool:
    import matplotlib.pyplot as plt

    if data is None or not len(data):
        _empty_panel(
            output,
            "关键残基 RMSF",
            "The residue RMSF XVG file was not found.",
        )
        return False
    residues = data[:, 0]
    values = data[:, 1] * 10.0
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.plot(residues, values, color="#7d3c98", linewidth=1.2)
    ax.set_xlabel("Residue index")
    ax.set_ylabel("RMSF (Angstrom)")
    ax.set_title("Protein residue RMSF", fontweight="bold")
    ax.grid(alpha=0.22)
    if binding_site.exists():
        numbers = [
            int(value)
            for value in re.findall(r"\d+", binding_site.read_text(encoding="utf-8"))
        ]
        if numbers:
            mask = np.isin(residues.astype(int), numbers)
            if mask.any():
                ax.scatter(
                    residues[mask],
                    values[mask],
                    color="#c0392b",
                    s=18,
                    label="Binding-site residues",
                    zorder=3,
                )
                ax.legend(frameon=False, fontsize=8)
    save_figure(fig, output)
    return True


def _parse_mmpbsa_decomposition(run_dir: Path) -> pd.DataFrame:
    candidates = list(run_dir.rglob("*DECOMP*MMPBSA*.dat")) + list(
        run_dir.rglob("FINAL_DECOMP_MMPBSA.dat")
    )
    if not candidates:
        return pd.DataFrame()
    path = candidates[0]
    records: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        numbers = []
        for token in reversed(tokens):
            try:
                numbers.append(float(token))
            except ValueError:
                break
        if len(numbers) < 2:
            continue
        label = " ".join(tokens[: max(1, len(tokens) - len(numbers))])
        if not re.search(r"\d", label):
            continue
        records.append(
            {
                "residue": label,
                "delta_kj_mol": numbers[-1],
            }
        )
    return pd.DataFrame(records)


def _mmpbsa_panel(
    docking_dir: Path,
    run_dir: Path | None,
    metrics: pd.DataFrame,
    output: Path,
) -> dict[str, Any]:
    import matplotlib.pyplot as plt

    decomposition = (
        _parse_mmpbsa_decomposition(run_dir)
        if run_dir is not None
        else pd.DataFrame()
    )
    if not decomposition.empty:
        frame = (
            decomposition.sort_values("delta_kj_mol")
            .head(20)
            .sort_values("delta_kj_mol")
        )
        fig, ax = plt.subplots(figsize=(7.2, 5.0))
        ax.barh(
            frame["residue"],
            frame["delta_kj_mol"],
            color="#456b8c",
        )
        ax.set_xlabel("MM-PBSA contribution (kJ/mol)")
        ax.set_title("MM-PBSA residue energy decomposition", fontweight="bold")
        save_figure(fig, output)
        return {
            "status": "completed",
            "mode": "residue_decomposition",
            "n_residues": int(len(decomposition)),
        }

    total = None
    if not metrics.empty:
        for column in (
            "mmpbsa_delta_total_kj_mol",
            "mmpbsa_delta_g",
            "mmpbsa_total_kj_mol",
        ):
            if column in metrics.columns:
                values = pd.to_numeric(metrics[column], errors="coerce").dropna()
                if not values.empty:
                    total = float(values.iloc[0])
                    break
    if total is not None and np.isfinite(total):
        fig, ax = plt.subplots(figsize=(7.2, 4.4))
        ax.bar(["MM-PBSA total"], [total], color="#456b8c")
        ax.axhline(0, color="#333333", linewidth=0.8)
        ax.set_ylabel("Delta G (kJ/mol)")
        ax.set_title("MM-PBSA binding free energy", fontweight="bold")
        ax.text(
            0,
            total,
            f"{total:.2f}",
            ha="center",
            va="bottom" if total >= 0 else "top",
            fontsize=9,
        )
        save_figure(fig, output)
        return {
            "status": "completed",
            "mode": "total_energy_only",
            "delta_total_kj_mol": total,
            "note": "Residue decomposition was not present in the external output.",
        }
    _empty_panel(
        output,
        "MM-PBSA 结合自由能分解",
        "No parsed MM-PBSA decomposition or total-energy output was found. "
        "Configure md.mmpbsa_command and run the 100 ns trajectory.",
    )
    return {
        "status": "unavailable",
        "mode": "not_available",
        "reason": "gmx_MMPBSA output not found",
    }


def generate_plan_md_figures(
    docking_dir: Path,
    output_dir: Path,
) -> dict[str, Any]:
    """Generate final Figure-5d-h panels from real trajectory files."""
    output_dir = ensure_dir(output_dir)
    run_dir = _latest_run_dir(docking_dir)
    metrics = _metrics(docking_dir)
    if run_dir is None:
        for panel, title in (
            ("d_rmsd", "蛋白主链 RMSD"),
            ("e_ligand_rmsd", "配体 RMSD"),
            ("f_rmsf", "关键残基 RMSF"),
            ("g_rg", "回旋半径 Rg"),
        ):
            _empty_panel(
                output_dir / f"fig5{panel}.png",
                title,
                "No completed GROMACS run directory was found. "
                "Set md.run=true to execute the 100 ns workflow.",
            )
        mmpbsa = _mmpbsa_panel(
            docking_dir,
            None,
            metrics,
            output_dir / "fig5h_mmpbsa.png",
        )
        result = {
            "status": "not_run",
            "run_dir": "",
            "panels": {
                "d_rmsd": False,
                "e_ligand_rmsd": False,
                "f_rmsf": False,
                "g_rg": False,
                "h_mmpbsa": mmpbsa["status"] == "completed",
            },
            "mmpbsa": mmpbsa,
        }
        write_json(output_dir / "md_plan_figures.json", result)
        return result

    panels = {
        "d_rmsd": _time_series_panel(
            _parse_xvg(run_dir / "rmsd_protein.xvg"),
            output_dir / "fig5d_rmsd.png",
            title="Protein backbone RMSD",
            ylabel="RMSD (nm)",
            color="#1665c0",
        ),
        "e_ligand_rmsd": _time_series_panel(
            _parse_xvg(run_dir / "rmsd_ligand.xvg"),
            output_dir / "fig5e_ligand_rmsd.png",
            title="Ligand RMSD",
            ylabel="RMSD (nm)",
            color="#c0392b",
        ),
        "f_rmsf": _rmsf_panel(
            _parse_xvg(run_dir / "rmsf_protein_residue.xvg"),
            output_dir / "fig5f_rmsf.png",
            run_dir / "binding_site_residues.txt",
        ),
        "g_rg": _time_series_panel(
            _parse_xvg(run_dir / "gyrate_protein.xvg"),
            output_dir / "fig5g_rg.png",
            title="Protein radius of gyration",
            ylabel="Rg (Angstrom)",
            color="#2e7d32",
            scale=10.0,
        ),
    }
    mmpbsa = _mmpbsa_panel(
        docking_dir,
        run_dir,
        metrics,
        output_dir / "fig5h_mmpbsa.png",
    )
    panels["h_mmpbsa"] = mmpbsa["status"] == "completed"
    result = {
        "status": "completed" if any(panels.values()) else "unavailable",
        "run_dir": str(run_dir),
        "panels": panels,
        "mmpbsa": mmpbsa,
    }
    write_json(output_dir / "md_plan_figures.json", result)
    LOG.info("Figure 5 MD panels generated: %s", result["panels"])
    return result
