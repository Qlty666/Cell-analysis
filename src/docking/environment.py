#!/usr/bin/env python3
"""Environment checks for the docking pipeline."""

from __future__ import annotations

import importlib.metadata as metadata
import importlib.util
import sys

from .utils import find_script, find_tool

PACKAGE_INFO = {
    "numpy": ("numpy", "https://numpy.org/"),
    "pandas": ("pandas", "https://pandas.pydata.org/"),
    "yaml": ("PyYAML", "https://pyyaml.org/"),
    "matplotlib": ("matplotlib", "https://matplotlib.org/"),
    "openpyxl": ("openpyxl", "https://openpyxl.readthedocs.io/"),
    "rdkit": ("rdkit", "https://www.rdkit.org/"),
    "meeko": ("meeko", "https://github.com/forlilab/Meeko"),
    "gemmi": ("gemmi", "https://project-gemmi.github.io/"),
    "openbabel": (
        "openbabel-wheel",
        "https://github.com/openbabel/openbabel-wheel",
    ),
    "sklearn": ("scikit-learn", "https://scikit-learn.org/"),
    "joblib": ("joblib", "https://joblib.readthedocs.io/"),
    "torch": ("torch", "https://pytorch.org/"),
    "versioneer": (
        "versioneer",
        "https://github.com/python-versioneer/python-versioneer",
    ),
}

TOOL_INFO = [
    (
        "vina",
        "AutoDock Vina",
        "https://github.com/ccsb-scripps/AutoDock-Vina/releases",
        "conda install -c conda-forge autodock-vina",
    ),
    (
        "obabel",
        "Open Babel",
        "https://openbabel.org/",
        "python -m pip install openbabel-wheel",
    ),
    (
        "mk_prepare_ligand.py",
        "Meeko ligand prep",
        "https://github.com/forlilab/Meeko",
        "python -m pip install meeko",
    ),
    (
        "mk_prepare_receptor.py",
        "Meeko receptor prep",
        "https://github.com/forlilab/Meeko",
        "python -m pip install meeko",
    ),
    (
        "prepare_receptor4.py",
        "AutoDockTools",
        "https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3",
        "python -m pip install ./dock/tools/AutoDockTools_py3",
    ),
]


def check_environment() -> list[dict]:
    checks = [
        {
            "name": "Python",
            "kind": "software",
            "version": sys.version.split()[0],
            "ok": sys.version_info >= (3, 10),
            "url": "https://www.python.org/",
            "install": "",
            "detail": sys.version.split()[0],
        }
    ]
    for module, (package, url) in PACKAGE_INFO.items():
        found = importlib.util.find_spec(module) is not None
        version = ""
        if found:
            try:
                version = metadata.version(package)
            except Exception:
                version = "installed"
        checks.append(
            {
                "name": package,
                "kind": "package",
                "version": version or "missing",
                "ok": found,
                "url": url,
                "install": f"python -m pip install {package}",
                "detail": version or "missing",
            }
        )
    for tool, label, url, install in TOOL_INFO:
        path = find_script(tool) if tool.endswith(".py") else find_tool(tool)
        checks.append(
            {
                "name": label,
                "kind": "software",
                "version": path or "not found",
                "ok": path is not None,
                "url": url,
                "install": install,
                "detail": path or "not found",
            }
        )
    return checks
