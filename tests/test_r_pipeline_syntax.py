#!/usr/bin/env python3
"""Smoke tests for the R analysis pipeline scripts.

These tests are deliberately cheap: they parse every R script and check that
the pipeline driver can be sourced far enough to define its helpers. They do
not run any analysis stage and do not need input data.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent
if str(APP_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(APP_ROOT / "src"))

RSCRIPT = shutil.which("Rscript")
R_ANALYSIS_DIR = APP_ROOT / "src" / "analysis"
R_MODULES_DIR = R_ANALYSIS_DIR / "R"
R_PIPELINE = R_ANALYSIS_DIR / "analysis_pipeline.R"

# Helpers the driver and its modules must expose once sourced.
REQUIRED_DRIVER_FUNCTIONS = (
    "read_h5ad_matrix",
    "run_stage",
    "stage_allowed",
    "normalize_ensembl_ids",
)


def _r_string(value: str) -> str:
    """Return an R string literal for *value* (paths are made POSIX first)."""
    escaped = value.replace("\\", "/").replace("'", "\\'")
    return f"'{escaped}'"


def _r_scripts() -> list[Path]:
    scripts: list[Path] = []
    for folder in (R_ANALYSIS_DIR, APP_ROOT / "src" / "pipeline"):
        if folder.is_dir():
            scripts.extend(sorted(folder.glob("*.R")))
    return scripts


def _module_scripts() -> list[Path]:
    if not R_MODULES_DIR.is_dir():
        return []
    return sorted(R_MODULES_DIR.glob("*.R"))


def _run_r(
    code: str,
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    timeout: float = 300,
) -> subprocess.CompletedProcess:
    merged = os.environ.copy()
    if env:
        merged.update(env)
    return subprocess.run(
        [RSCRIPT, "-e", code],
        cwd=str(cwd or APP_ROOT),
        env=merged,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )


def _run_r_script(
    code: str,
    *,
    cwd: Path,
    env: dict[str, str] | None = None,
    trailing_args: tuple[str, ...] = (),
    timeout: float = 300,
) -> subprocess.CompletedProcess:
    """Run *code* from a script file so Rscript forwards trailing args.

    ``Rscript -e`` silently drops everything after the expression, so the
    ``--start-stage=`` argument the driver reads from ``commandArgs()`` is only
    visible in file mode.
    """
    script = Path(cwd) / "_r_smoke_driver.R"
    script.write_text(code, encoding="utf-8")
    cmd = [RSCRIPT, str(script)]
    if trailing_args:
        cmd += ["--args", *trailing_args]
    merged = os.environ.copy()
    if env:
        merged.update(env)
    return subprocess.run(
        cmd,
        cwd=str(cwd),
        env=merged,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )


@unittest.skipUnless(RSCRIPT, "Rscript not installed")
class TestRPipelineSyntax(unittest.TestCase):
    def test_every_r_script_parses(self):
        scripts = _r_scripts()
        self.assertTrue(scripts, "no R scripts found under src/analysis or src/pipeline")
        for script in scripts:
            with self.subTest(script=script.name):
                proc = _run_r(
                    f"invisible(parse({_r_string(str(script))}))",
                    timeout=120,
                )
                self.assertEqual(
                    proc.returncode,
                    0,
                    msg=(
                        f"Rscript failed to parse {script}:\n"
                        f"{proc.stderr or proc.stdout}"
                    ),
                )

    def test_driver_sources_and_defines_helpers(self):
        checks = ", ".join(
            f'"{name}" = exists("{name}", mode = "function")'
            for name in REQUIRED_DRIVER_FUNCTIONS
        )
        code = (
            f"source({_r_string(str(R_PIPELINE))}, local = globalenv())\n"
            f"flags <- c({checks})\n"
            'cat("DRIVER_FUNCS", paste(names(flags), flags, sep = "="), "\\n")\n'
            "if (!all(flags)) quit(status = 3)\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            env = {}
            if R_MODULES_DIR.is_dir():
                env["LIVER_R_MODULES_DIR"] = str(R_MODULES_DIR)
            # --start-stage=99 disables every stage block, so sourcing only
            # exercises parsing, module loading and helper definition.
            proc = _run_r_script(
                code,
                cwd=Path(tmp),
                env=env,
                trailing_args=("--start-stage=99",),
                timeout=300,
            )
            self.assertEqual(
                proc.returncode,
                0,
                msg=(
                    "sourcing the pipeline driver failed:\n"
                    f"{proc.stderr or proc.stdout}"
                ),
            )
            self.assertIn("DRIVER_FUNCS", proc.stdout)
            for name in REQUIRED_DRIVER_FUNCTIONS:
                self.assertIn(f"{name}=TRUE", proc.stdout)

    def test_modules_source_standalone(self):
        modules = _module_scripts()
        if not modules:
            self.skipTest("no src/analysis/R module directory yet")
        for module in modules:
            with self.subTest(module=module.name):
                # Run from a file: this build's ``Rscript -e`` evaluates only
                # the first statement of a multi-statement string.
                code = (
                    f"source({_r_string(str(module))}, local = globalenv())\n"
                    "n <- length(lsf.str(envir = globalenv()))\n"
                    'cat("MODULE_FUNCS", n, "\\n")\n'
                    "if (n < 1) quit(status = 4)\n"
                )
                with tempfile.TemporaryDirectory() as tmp:
                    proc = _run_r_script(code, cwd=Path(tmp), timeout=300)
                self.assertEqual(
                    proc.returncode,
                    0,
                    msg=(
                        f"sourcing module {module} failed:\n"
                        f"{proc.stderr or proc.stdout}"
                    ),
                )
                self.assertIn("MODULE_FUNCS", proc.stdout)
                count = int(proc.stdout.split("MODULE_FUNCS", 1)[1].split()[0])
                self.assertGreaterEqual(
                    count, 1, msg=f"{module} defines no functions"
                )

    def test_snapshot_layout_resolves_modules(self):
        """The orchestrator snapshot (<logs>/pipeline_analysis.R + <logs>/R)."""
        modules = _module_scripts()
        if not modules:
            self.skipTest("no src/analysis/R module directory yet")
        with tempfile.TemporaryDirectory() as tmp:
            logs = Path(tmp) / "logs"
            logs.mkdir()
            shutil.copy2(R_PIPELINE, logs / "pipeline_analysis.R")
            shutil.copytree(R_MODULES_DIR, logs / "R")
            code = (
                f"source({_r_string(str(logs / 'pipeline_analysis.R'))}, "
                "local = globalenv()); "
                'cat("MODULES_DIR", liver_modules_dir, "\\n"); '
                "if (!exists('read_h5ad_matrix', mode = 'function')) "
                "quit(status = 5)"
            )
            # The env override is disabled on purpose: only the
            # <script_dir>/R branch can resolve the copied modules.
            proc = _run_r_script(
                code,
                cwd=logs,
                env={
                    "LIVER_R_MODULES_DIR": "",
                    "LIVER_ROOT": str(Path(tmp) / "out"),
                },
                trailing_args=("--start-stage=99",),
                timeout=300,
            )
            self.assertEqual(
                proc.returncode,
                0,
                msg=(
                    "sourcing the snapshot driver failed:\n"
                    f"{proc.stderr or proc.stdout}"
                ),
            )
            self.assertIn("MODULES_DIR", proc.stdout)
            expected = str((logs / "R").resolve()).replace("\\", "/").lower()
            actual = proc.stdout.replace("\\", "/").lower()
            self.assertIn(expected, actual)


if __name__ == "__main__":
    unittest.main()
