#!/usr/bin/env python3
"""Check whether the current computer can run the single-cell pipeline."""

import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REQUIRED_PYTHON = (3, 10)
REQUIRED_R = (4, 5, 0)

R_PACKAGES = [
    "Seurat",
    "Matrix",
    "data.table",
    "dplyr",
    "ggplot2",
    "patchwork",
    "jsonlite",
    "ggrepel",
    "pheatmap",
    "scDblFinder",
    "SingleCellExperiment",
    "BiocParallel",
    "clusterProfiler",
    "org.Hs.eg.db",
    "org.Mm.eg.db",
    "enrichplot",
    "DESeq2",
]


def find_rscript() -> str | None:
    if shutil.which("Rscript"):
        return shutil.which("Rscript")
    if shutil.which("Rscript.exe"):
        return shutil.which("Rscript.exe")
    for base in [
        Path(r"C:\Program Files\R"),
        Path(r"C:\Program Files\Microsoft\R Open"),
        Path(r"C:\Users") / Path.home().name / r"AppData\Local\Programs\R",
    ]:
        if not base.exists():
            continue
        candidates = sorted(base.glob("R-*/bin/Rscript.exe"), reverse=True)
        if candidates:
            return str(candidates[0])
    return None


def parse_version(text: str) -> tuple:
    parts = []
    for part in text.strip().split("."):
        digits = "".join(ch for ch in part if ch.isdigit())
        if digits:
            parts.append(int(digits))
        else:
            break
    return tuple(parts)


def main() -> int:
    ok = True
    checks = []

    def report(name: str, passed: bool, detail: str) -> None:
        nonlocal ok
        checks.append((name, passed, detail))
        if not passed:
            ok = False

    py = sys.version_info
    report(
        "Python",
        py >= REQUIRED_PYTHON,
        f"{py.major}.{py.minor}.{py.micro} (required >= "
        f"{REQUIRED_PYTHON[0]}.{REQUIRED_PYTHON[1]})",
    )

    rscript = find_rscript()
    report("Rscript", rscript is not None, rscript or "Rscript not found in PATH")

    r_version = ""
    if rscript:
        try:
            result = subprocess.run(
                [rscript, "-e", "cat(as.character(getRversion()))"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=60,
            )
            r_version = result.stdout.strip()
            report(
                "R version",
                parse_version(r_version) >= REQUIRED_R,
                f"{r_version} (required >= "
                f"{'.'.join(str(x) for x in REQUIRED_R)})",
            )
        except Exception as exc:
            report("R version", False, f"could not query R version: {exc}")

    curl = shutil.which("curl") or shutil.which("curl.exe")
    report("curl", curl is not None, curl or "curl not found in PATH")

    if rscript:
        pkg_script = (
            "pkgs <- c("
            + ",".join(f'"{p}"' for p in R_PACKAGES)
            + "); ip <- installed.packages(); "
            "for(p in pkgs) cat(p, '=', "
            "ifelse(p %in% rownames(ip), ip[p,'Version'], 'MISSING'), '\\n')"
        )
        try:
            result = subprocess.run(
                [rscript, "-e", pkg_script],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=180,
            )
            missing = []
            for line in result.stdout.splitlines():
                if "=" in line:
                    name, version = line.split("=", 1)
                    name = name.strip()
                    version = version.strip()
                    if version == "MISSING":
                        missing.append(name)
            report(
                "R packages",
                not missing,
                "all installed" if not missing else "missing: " + ", ".join(missing),
            )
        except Exception as exc:
            report("R packages", False, f"could not query packages: {exc}")
    else:
        report("R packages", False, "Rscript unavailable")

    try:
        req = urllib.request.Request(
            "https://ftp.ncbi.nlm.nih.gov/",
            method="HEAD",
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(req, timeout=15) as resp:
            reachable = resp.status < 400
        report("NCBI GEO access", reachable, "reachable" if reachable else "unreachable")
    except Exception as exc:
        report("NCBI GEO access", False, str(exc))

    print()
    print("Environment check result")
    print("=" * 60)
    for name, passed, detail in checks:
        status = "OK " if passed else "FAIL"
        print(f"[{status}] {name}: {detail}")
    print("=" * 60)
    print("Result:", "PASS" if ok else "FAIL")
    if not ok:
        print("Run install_pipeline_dependencies.py to install missing R packages.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
