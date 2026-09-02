---
name: liver-expression-analysis
description: Run local GEO/ArrayExpress/BioStudies expression analysis (single-cell, bulk RNA-seq, or microarray) with the Liver Cancer Bioinformatics project scripts. Use only when working in that project or when LIVER_PIPELINE_ROOT points to it.
---

# Liver Expression Analysis

Run the project's expression pipeline rather than reimplementing it.

## Locate the project

- Use the current working directory or its nearest ancestor that contains both `AGENTS.md` and `config/full_pipeline_config.json`.
- If that search fails and `LIVER_PIPELINE_ROOT` is set, use that directory.
- If neither is available, ask the user for the project path before running commands.

## Quick Start

```bash
python scripts/liverbio.py expression GSE125449 --output ../liver_cancer --species auto
```

The equivalent raw command is `python scripts/run_pipeline.py`. Prefer the `liverbio` dispatcher so the feature entry point stays consistent.

## Operating rules

- Use `--help` on the underlying script to discover every available option before inventing parameters.
- Pass `--output` as a project-independent absolute or repository-relative path.
- Keep `--species auto` unless the user states human (`hs`) or mouse (`mm`).
- Rerunning the same command without `--force` resumes from the last completed checkpoint. Use `--force` only when the user explicitly asks for a full rerun.
- Do not change `src/analysis/*.R`, `src/pipeline/orchestrator.py`, or data-download modules unless the user asks for code changes.
- If R or Python dependencies are missing, run `python scripts/liverbio.py doctor pipeline` and follow the launcher installers.

## Useful options

The pipeline supports environment-controlled features documented in `README.md`: QC thresholds, ML model selection, GSEA limits, figure skipping, and dataset mode detection. Read the README or `python scripts/run_pipeline.py --help` when a user request depends on those choices.

## References

- `README.md`: installation, CLI usage, outputs, and environment variables.
- `config/full_pipeline_config.json`: integrated-run defaults.
