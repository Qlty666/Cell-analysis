"""Configuration helpers for the standalone molecular docking board."""

from __future__ import annotations

import json
from pathlib import Path

from docking.config import ResolvedConfig, load_config as load_docking_config

APP_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_CONFIG = APP_ROOT / "config" / "molecular_docking_config.json"


def load_config(
    config_path: str | Path | None = None,
    overrides: dict | None = None,
) -> ResolvedConfig:
    """Load a molecular docking config with docking defaults applied."""
    return load_docking_config(
        Path(config_path).resolve() if config_path else DEFAULT_CONFIG,
        overrides,
    )


def save_config(cfg: ResolvedConfig, path: Path) -> None:
    """Save only the fields the standalone docking board manages."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "name": cfg.data.get("name", "molecular_docking"),
        "workdir": str(cfg.workdir),
        "output_dir": _rel(cfg.output_dir, cfg.workdir),
        "receptor": {
            "input": _rel(cfg.receptor_input(), cfg.workdir),
            "output": _rel(cfg.receptor_output(), cfg.workdir),
            "detect_input": cfg.get("receptor", "detect_input"),
            "center": cfg.receptor_center(),
            "size": cfg.receptor_size(),
            "flexible": [_rel(item, cfg.workdir) for item in cfg.receptor_flexible()],
        },
        "ligand": {
            "input": _rel(cfg.ligand_input(), cfg.workdir),
            "output_dir": _rel(cfg.ligand_output_dir(), cfg.workdir),
            "smiles_column": cfg.get("ligand", "smiles_column", "SMILES"),
            "id_column": cfg.get("ligand", "id_column", "ID"),
            "ph": cfg.get("ligand", "ph", 7.4),
            "remove_salts": cfg.get("ligand", "remove_salts", True),
            "neutralize": cfg.get("ligand", "neutralize", True),
            "max_heavy_atoms": cfg.get("ligand", "max_heavy_atoms", 60),
            "max_rotatable_bonds": cfg.get("ligand", "max_rotatable_bonds", 15),
            "max_ligands": cfg.get("ligand", "max_ligands"),
            "conformers": cfg.get("ligand", "conformers", 1),
            "seed": cfg.get("ligand", "seed", 42),
            "engine": cfg.get("ligand", "engine", "auto"),
        },
        "docking": {
            "engine": cfg.get("docking", "engine", "vina"),
            "executable": cfg.get("docking", "executable", "vina"),
            "scoring": cfg.get("docking", "scoring", "vina"),
            "exhaustiveness": cfg.get("docking", "exhaustiveness", 8),
            "num_modes": cfg.get("docking", "num_modes", 9),
            "energy_range": cfg.get("docking", "energy_range", 3.0),
            "cpu": cfg.get("docking", "cpu", 4),
            "max_workers": cfg.get("docking", "max_workers", 4),
            "seed": cfg.get("docking", "seed", 42),
            "timeout_seconds": cfg.get("docking", "timeout_seconds", 600),
            "resume": cfg.get("docking", "resume", True),
        },
        "analysis": {
            "cutoff": cfg.get("analysis", "cutoff", -7.0),
            "top_n": cfg.get("analysis", "top_n", 100),
            "figures": cfg.get("analysis", "figures", True),
            "diversity": cfg.get("analysis", "diversity", True),
            "tanimoto_cutoff": cfg.get("analysis", "tanimoto_cutoff", 0.7),
        },
        "redock": {
            "enabled": cfg.get("redock", "enabled", True),
            "top_n": cfg.get("redock", "top_n", 20),
            "exhaustiveness": cfg.get("redock", "exhaustiveness", 32),
            "num_modes": cfg.get("redock", "num_modes", 9),
            "energy_range": cfg.get("redock", "energy_range", 3.0),
            "max_workers": cfg.get("redock", "max_workers", 4),
            "timeout_seconds": cfg.get("redock", "timeout_seconds", 600),
            "resume": cfg.get("redock", "resume", True),
        },
        "report": {
            "top_n": cfg.get("report", "top_n", 20),
        },
    }
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )


def _rel(path: Path, base: Path) -> str:
    path = Path(path).resolve()
    base = Path(base).resolve()
    try:
        value = path.relative_to(base).as_posix()
    except ValueError:
        value = str(path)
    return value
