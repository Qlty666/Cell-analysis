#!/usr/bin/env python3
"""Load, validate and serialize docking configuration."""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path

try:
    import yaml
except ImportError:  # YAML configs are optional; JSON is the default format.
    yaml = None

DEFAULTS = {
    "name": "liver_cancer_vs",
    "workdir": ".",
    "output_dir": "outputs/run_001",
    "receptor": {
        "input": "data/receptors/receptor.pdb",
        "output": "data/receptors/receptor.pdbqt",
        "detect_input": None,
        "center": [0.0, 0.0, 0.0],
        "size": [22.5, 22.5, 22.5],
        "flexible": [],
    },
    "ligand": {
        "input": "data/ligands/library.sdf",
        "output_dir": "data/ligands/prepared",
        "smiles_column": "SMILES",
        "id_column": "ID",
        "ph": 7.4,
        "remove_salts": True,
        "neutralize": True,
        "max_heavy_atoms": 60,
        "max_rotatable_bonds": 15,
        "max_ligands": None,
        "conformers": 1,
        "seed": 42,
        "engine": "auto",
    },
    "docking": {
        "engine": "vina",
        "executable": "vina",
        "scoring": "vina",
        "exhaustiveness": 8,
        "num_modes": 9,
        "energy_range": 3.0,
        "cpu": 4,
        "max_workers": 4,
        "seed": 42,
        "timeout_seconds": 600,
        "resume": True,
    },
    "analysis": {
        "cutoff": -7.0,
        "top_n": 100,
        "figures": True,
        "diversity": True,
        "tanimoto_cutoff": 0.7,
    },
    "evidence": {
        "uniprot_accession": None,
        "pdb_id": None,
        "chembl_target_id": None,
        "target_name": None,
        "ligand_name": None,
        "ligand_smiles": None,
        "max_items": 10,
    },
    "ml": {
        "model": "rf",
        "task": "auto",
        "label_column": "active",
        "training_csv": "data/ml/training.csv",
        "test_size": 0.2,
        "random_state": 42,
        "epochs": 80,
        "hidden_size": 128,
    },
    "md": {
        "export_dir": "outputs/md",
        "top_n": 10,
    },
    "redock": {
        "enabled": True,
        "top_n": 20,
        "exhaustiveness": 32,
        "num_modes": 9,
        "energy_range": 3.0,
        "max_workers": 4,
        "timeout_seconds": 600,
        "resume": True,
    },
    "report": {
        "top_n": 20,
    },
    "knockout": {
        "enabled": True,
        "expression_csv": "data/knockout/expression.csv",
        "metadata_csv": None,
        "ppi_network_csv": None,
        "depmap_csv": None,
        "prognosis_csv": None,
        "druggability_csv": None,
        "off_target_csv": None,
        "group_column": "condition",
        "cell_type_column": None,
        "case_label": None,
        "normal_label": None,
        "top_n": 50,
        "figures": True,
        "max_genes": 2000,
        "max_samples": 5000,
        "corr_cutoff": 0.7,
        "liver_lineage": "liver",
        "off_target_penalty": 0.05,
        "disease_up_genes": None,
        "disease_down_genes": None,
        "pathway_genes": None,
        "weights": {
            "expression": 0.35,
            "proliferation": 0.25,
            "network": 0.20,
            "depmap": 0.20,
        },
        "target_weights": {
            "base": 0.35,
            "reversal": 0.20,
            "pathway": 0.15,
            "specificity": 0.10,
            "prognosis": 0.10,
            "druggability": 0.10,
            "ppi_hub": 0.10,
        },
    },
    "network_toxicology": {
        "compound_name": None,
        "disease_name": None,
        "compound_targets_csv": None,
        "target_sources": None,
        "disease_genes_csv": None,
        "disease_gene_column": None,
        "ppi_network_csv": None,
        "output_dir": "outputs/run_001/network_toxicology",
        "venn": True,
    },
    "faers": {
        "input_csv": None,
        "drug_column": "drug",
        "event_column": "event",
        "count_column": None,
        "min_count": 3,
        "output_dir": "outputs/run_001/faers",
    },
    "validation": {
        "top_n": 10,
        "output_dir": "outputs/run_001/validation",
    },
}


def deep_merge(base: dict, extra: dict) -> dict:
    out = copy.deepcopy(base)
    for key, value in (extra or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


class ResolvedConfig:
    def __init__(self, data: dict, config_path: Path):
        self.data = data
        self.config_path = Path(config_path).resolve()
        self.root = self.config_path.parent
        env_workdir = os.environ.get("DOCK_WORKDIR")
        if env_workdir:
            data["workdir"] = env_workdir
        self.workdir = self._resolve(data.get("workdir") or ".", self.root)
        self.output_dir = self._resolve(
            data.get("output_dir") or "outputs/run_001", self.workdir
        )

    def _resolve(self, value, base: Path | None = None) -> Path:
        path = Path(str(value)).expanduser()
        if not path.is_absolute():
            path = (base or self.workdir) / path
        return path.resolve()

    def get(self, section: str, key: str, default=None):
        return self.data.get(section, {}).get(key, default)

    def receptor_input(self) -> Path:
        return self._resolve(self.get("receptor", "input"))

    def receptor_output(self) -> Path:
        return self._resolve(self.get("receptor", "output"))

    def receptor_center(self) -> list[float]:
        return [float(v) for v in self.get("receptor", "center", [0.0, 0.0, 0.0])]

    def receptor_size(self) -> list[float]:
        return [float(v) for v in self.get("receptor", "size", [22.5, 22.5, 22.5])]

    def receptor_flexible(self) -> list[Path]:
        return [
            self._resolve(str(item), self.workdir)
            for item in self.get("receptor", "flexible", [])
        ]

    def ligand_input(self) -> Path:
        return self._resolve(self.get("ligand", "input"))

    def ligand_output_dir(self) -> Path:
        return self._resolve(self.get("ligand", "output_dir"))

    def docked_dir(self) -> Path:
        return self.output_dir / "docked"

    def reports_dir(self) -> Path:
        return self.output_dir / "results"

    def analysis_dir(self) -> Path:
        return self.reports_dir() / "01_analysis"

    def redock_dir(self) -> Path:
        return self.reports_dir() / "02_redock"

    def ml_dir(self) -> Path:
        return self.reports_dir() / "03_ml"

    def knockout_dir(self) -> Path:
        return self.reports_dir() / "04_knockout"

    def validation_dir(self) -> Path:
        return self.reports_dir() / "05_validation"

    def logs_dir(self) -> Path:
        return self.output_dir / "logs"

    def stage_dir(self) -> Path:
        return self.output_dir / ".stages"

    def results_path(self) -> Path:
        return self.docked_dir() / "results.csv"

    def manifest_path(self) -> Path:
        return self.ligand_output_dir() / "manifest.csv"

    def ml_training_csv(self) -> Path:
        return self._resolve(self.get("ml", "training_csv", "data/ml/training.csv"))

    def md_export_dir(self) -> Path:
        return self._resolve(self.get("md", "export_dir", "outputs/md"))

    def validate(self) -> None:
        for name, values in [
            ("center", self.receptor_center()),
            ("size", self.receptor_size()),
        ]:
            if not isinstance(values, (list, tuple)) or len(values) != 3:
                raise ValueError(f"receptor.{name} must be a list of 3 numbers")
        if any(v <= 0 for v in self.receptor_size()):
            raise ValueError("receptor.size values must be positive")
        for section, key in [
            ("docking", "exhaustiveness"),
            ("docking", "num_modes"),
            ("docking", "cpu"),
            ("docking", "max_workers"),
        ]:
            if int(self.get(section, key, 0)) < 1:
                raise ValueError(f"{section}.{key} must be >= 1")
        if float(self.get("docking", "energy_range", 3.0)) <= 0:
            raise ValueError("docking.energy_range must be positive")


def load_config(
    config_path: str | Path | None = None,
    overrides: dict | None = None,
) -> ResolvedConfig:
    if config_path is None:
        config_path = (
            Path(__file__).resolve().parent.parent.parent
            / "config"
            / "docking_config.json"
        )
    config_path = Path(config_path).resolve()
    raw = {}
    if config_path.exists():
        text = config_path.read_text(encoding="utf-8")
        if config_path.suffix.lower() in (".yaml", ".yml"):
            if yaml is None:
                raise RuntimeError("PyYAML is required to load YAML docking configs")
            raw = yaml.safe_load(text) or {}
        else:
            raw = json.loads(text) or {}
    data = deep_merge(DEFAULTS, raw)
    cfg = ResolvedConfig(data, config_path)
    if overrides:
        cfg = apply_overrides(cfg, overrides)
    cfg.validate()
    return cfg


def apply_overrides(cfg: ResolvedConfig, overrides: dict) -> ResolvedConfig:
    data = cfg.data
    for key, dotted in [
        ("workdir", "workdir"),
        ("outdir", "output_dir"),
        ("receptor", "receptor/input"),
        ("ligand", "ligand/input"),
    ]:
        if overrides.get(key) is not None:
            set_dotted(data, dotted, overrides[key])

    section_map = {
        "center": ("receptor", "center"),
        "size": ("receptor", "size"),
        "flexible": ("receptor", "flexible"),
        "exhaustiveness": ("docking", "exhaustiveness"),
        "num_modes": ("docking", "num_modes"),
        "energy_range": ("docking", "energy_range"),
        "cpu": ("docking", "cpu"),
        "max_workers": ("docking", "max_workers"),
        "seed": ("docking", "seed"),
        "scoring": ("docking", "scoring"),
        "executable": ("docking", "executable"),
        "max_ligands": ("ligand", "max_ligands"),
        "ligand_engine": ("ligand", "engine"),
        "conformers": ("ligand", "conformers"),
        "cutoff": ("analysis", "cutoff"),
        "top_n": ("analysis", "top_n"),
        "figures": ("analysis", "figures"),
        "diversity": ("analysis", "diversity"),
        "model": ("ml", "model"),
        "label_column": ("ml", "label_column"),
        "training_csv": ("ml", "training_csv"),
        "epochs": ("ml", "epochs"),
        "hidden_size": ("ml", "hidden_size"),
        "uniprot": ("evidence", "uniprot_accession"),
        "pdb": ("evidence", "pdb_id"),
        "chembl_target": ("evidence", "chembl_target_id"),
        "target_name": ("evidence", "target_name"),
        "ligand_name": ("evidence", "ligand_name"),
        "ligand_smiles": ("evidence", "ligand_smiles"),
        "max_items": ("evidence", "max_items"),
        "expression_csv": ("knockout", "expression_csv"),
        "metadata_csv": ("knockout", "metadata_csv"),
        "depmap_csv": ("knockout", "depmap_csv"),
        "prognosis_csv": ("knockout", "prognosis_csv"),
        "druggability_csv": ("knockout", "druggability_csv"),
        "off_target_csv": ("knockout", "off_target_csv"),
        "cell_type_column": ("knockout", "cell_type_column"),
        "group_column": ("knockout", "group_column"),
        "case_label": ("knockout", "case_label"),
        "normal_label": ("knockout", "normal_label"),
        "ko_top_n": ("knockout", "top_n"),
        "validation_top_n": ("validation", "top_n"),
        "compound_name": ("network_toxicology", "compound_name"),
        "disease_name": ("network_toxicology", "disease_name"),
        "compound_targets_csv": ("network_toxicology", "compound_targets_csv"),
        "disease_genes_csv": ("network_toxicology", "disease_genes_csv"),
        "disease_gene_column": ("network_toxicology", "disease_gene_column"),
        "network_output_dir": ("network_toxicology", "output_dir"),
        "faers_input": ("faers", "input_csv"),
        "faers_drug_column": ("faers", "drug_column"),
        "faers_event_column": ("faers", "event_column"),
        "faers_count_column": ("faers", "count_column"),
        "faers_min_count": ("faers", "min_count"),
    }
    for key, (section, field) in section_map.items():
        if overrides.get(key) is not None:
            data.setdefault(section, {})[field] = overrides[key]
    if overrides.get("ppi_network_csv") is not None:
        data.setdefault("knockout", {})["ppi_network_csv"] = overrides[
            "ppi_network_csv"
        ]
        data.setdefault("network_toxicology", {})["ppi_network_csv"] = overrides[
            "ppi_network_csv"
        ]
    return ResolvedConfig(data, cfg.config_path)


def set_dotted(data: dict, dotted: str, value) -> None:
    parts = dotted.split("/")
    target = data
    for part in parts[:-1]:
        target = target.setdefault(part, {})
    target[parts[-1]] = value


def save_config(cfg: ResolvedConfig, path: Path) -> None:
    data = {
        "name": cfg.data.get("name", "virtual_screening"),
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
            "max_ligands": cfg.get("ligand", "max_ligands", None),
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
        "evidence": {
            "uniprot_accession": cfg.get("evidence", "uniprot_accession"),
            "pdb_id": cfg.get("evidence", "pdb_id"),
            "chembl_target_id": cfg.get("evidence", "chembl_target_id"),
            "target_name": cfg.get("evidence", "target_name"),
            "ligand_name": cfg.get("evidence", "ligand_name"),
            "ligand_smiles": cfg.get("evidence", "ligand_smiles"),
            "max_items": cfg.get("evidence", "max_items", 10),
        },
        "ml": {
            "model": cfg.get("ml", "model", "rf"),
            "task": cfg.get("ml", "task", "auto"),
            "label_column": cfg.get("ml", "label_column", "active"),
            "training_csv": cfg.get("ml", "training_csv", "data/ml/training.csv"),
            "test_size": cfg.get("ml", "test_size", 0.2),
            "random_state": cfg.get("ml", "random_state", 42),
            "epochs": cfg.get("ml", "epochs", 80),
            "hidden_size": cfg.get("ml", "hidden_size", 128),
        },
        "md": {
            "export_dir": cfg.get("md", "export_dir", "outputs/md"),
            "top_n": cfg.get("md", "top_n", 10),
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
        "knockout": {
            "enabled": cfg.get("knockout", "enabled", True),
            "expression_csv": cfg.get(
                "knockout", "expression_csv", "data/knockout/expression.csv"
            ),
            "metadata_csv": cfg.get("knockout", "metadata_csv"),
            "ppi_network_csv": cfg.get("knockout", "ppi_network_csv"),
            "depmap_csv": cfg.get("knockout", "depmap_csv"),
            "prognosis_csv": cfg.get("knockout", "prognosis_csv"),
            "druggability_csv": cfg.get("knockout", "druggability_csv"),
            "off_target_csv": cfg.get("knockout", "off_target_csv"),
            "group_column": cfg.get("knockout", "group_column", "condition"),
            "cell_type_column": cfg.get("knockout", "cell_type_column"),
            "case_label": cfg.get("knockout", "case_label"),
            "normal_label": cfg.get("knockout", "normal_label"),
            "top_n": cfg.get("knockout", "top_n", 50),
            "figures": cfg.get("knockout", "figures", True),
            "max_genes": cfg.get("knockout", "max_genes", 2000),
            "max_samples": cfg.get("knockout", "max_samples", 5000),
            "corr_cutoff": cfg.get("knockout", "corr_cutoff", 0.7),
            "liver_lineage": cfg.get("knockout", "liver_lineage", "liver"),
            "off_target_penalty": cfg.get("knockout", "off_target_penalty", 0.05),
            "disease_up_genes": cfg.get("knockout", "disease_up_genes"),
            "disease_down_genes": cfg.get("knockout", "disease_down_genes"),
            "pathway_genes": cfg.get("knockout", "pathway_genes"),
            "weights": cfg.get(
                "knockout",
                "weights",
                {
                    "expression": 0.35,
                    "proliferation": 0.25,
                    "network": 0.20,
                    "depmap": 0.20,
                },
            ),
            "target_weights": cfg.get(
                "knockout",
                "target_weights",
                {
                    "base": 0.35,
                    "reversal": 0.20,
                    "pathway": 0.15,
                    "specificity": 0.10,
                    "prognosis": 0.10,
                    "druggability": 0.10,
                    "ppi_hub": 0.10,
                },
            ),
        },
        "network_toxicology": {
            "compound_name": cfg.get("network_toxicology", "compound_name"),
            "disease_name": cfg.get("network_toxicology", "disease_name"),
            "compound_targets_csv": cfg.get(
                "network_toxicology", "compound_targets_csv"
            ),
            "target_sources": cfg.get("network_toxicology", "target_sources"),
            "disease_genes_csv": cfg.get(
                "network_toxicology", "disease_genes_csv"
            ),
            "disease_gene_column": cfg.get(
                "network_toxicology", "disease_gene_column"
            ),
            "ppi_network_csv": cfg.get("network_toxicology", "ppi_network_csv"),
            "output_dir": cfg.get(
                "network_toxicology",
                "output_dir",
                "outputs/run_001/network_toxicology",
            ),
            "venn": cfg.get("network_toxicology", "venn", True),
        },
        "faers": {
            "input_csv": cfg.get("faers", "input_csv"),
            "drug_column": cfg.get("faers", "drug_column", "drug"),
            "event_column": cfg.get("faers", "event_column", "event"),
            "count_column": cfg.get("faers", "count_column"),
            "min_count": cfg.get("faers", "min_count", 3),
            "output_dir": cfg.get(
                "faers", "output_dir", "outputs/run_001/faers"
            ),
        },
        "validation": {
            "top_n": cfg.get("validation", "top_n", 10),
            "output_dir": cfg.get(
                "validation", "output_dir", "outputs/run_001/validation"
            ),
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _rel(path: Path, base: Path) -> str:
    try:
        return os.path.relpath(path, base)
    except ValueError:
        return str(path)
