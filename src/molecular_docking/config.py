"""Configuration helpers for the standalone molecular docking board."""

from __future__ import annotations

from pathlib import Path

from docking.config import ResolvedConfig
from docking.config import load_config as load_docking_config
from docking.config import save_config as _save_config

APP_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_CONFIG = APP_ROOT / "config" / "molecular_docking_config.json"

# Top-level config keys the standalone board manages; everything else keeps
# its docking-pipeline defaults and is intentionally not written out.
MANAGED_SECTIONS = [
    "name",
    "workdir",
    "output_dir",
    "receptor",
    "ligand",
    "docking",
    "analysis",
    "redock",
    "report",
]


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
    _save_config(
        cfg,
        Path(path),
        sections=MANAGED_SECTIONS,
        posix_paths=True,
        strict_relative=True,
    )
