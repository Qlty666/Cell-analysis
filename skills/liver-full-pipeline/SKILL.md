---
name: liver-full-pipeline
description: Run the integrated Liver Cancer Bioinformatics workflow from expression analysis through key targets, evidence, knockout, docking, cell feedback, and report generation. Use only inside that project or when LIVER_PIPELINE_ROOT points to it.
---

# Liver Full Pipeline

Run the project's end-to-end integrated pipeline instead of composing ad hoc commands.

## Locate the project

- Use the current working directory or its nearest ancestor containing `AGENTS.md` and `config/full_pipeline_config.json`.
- If that search fails and `LIVER_PIPELINE_ROOT` is set, use that directory.
- If neither is available, ask the user for the project path before running commands.

## Quick Start

```bash
python scripts/liverbio.py full \
  --accession GSE125449 \
  --output ../liver_cancer \
  --workdir ../liver_cancer_full
```

The `--output` directory is the single-cell output root used by stage 01 and downstream stages; `--workdir` is the directory that receives the integrated `outputs/` tree. They are separate and should not be confused.

## Operating rules

- Run `python scripts/liverbio.py full --list-stages` or `--dry-run` before a long first run when the user wants to confirm scope.
- Do not pass `--force` by default; the stage markers preserve completed work and rerun only stale stages.
- Respect the stage order 01-08 and use `--start-stage` only when the user needs to resume from a known stage.
- Report missing R, Python, Vina, or evidence-skill dependencies instead of continuing with fabricated output.
- Keep user-provided labels (`--case-label`, `--normal-label`) verbatim; do not infer tumor/normal when the metadata is ambiguous.

## Useful options

The full pipeline accepts `--skip-scrna`, `--skip-download`, `--skip-evidence-fetch`, `--skip-knockout`, `--skip-docking`, and `--skip-cell-feedback` to run partial workflows. Additional options and config defaults are documented by `python scripts/run_full_pipeline.py --help` and `README.md`.

## References

- `README.md`: integrated pipeline usage and stage outputs.
- `config/full_pipeline_config.json`: default parameters, QC gate, and feedback settings.
