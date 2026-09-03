---
name: liver-virtual-screening
description: Run local AutoDock Vina virtual screening, evidence collection, ML rescoring, knockout, network toxicology, FAERS signal detection, and validation-export commands in the Liver Cancer Bioinformatics project. Use only when working in that project or when LIVER_PIPELINE_ROOT points to it.
---

# Liver Virtual Screening

Run the project's docking workflows through its tested CLI rather than recreating the docking logic.

## Locate the project

- Use the current working directory or its nearest ancestor containing `AGENTS.md`, `config/docking_config.json`, and `scripts/run_docking.py`.
- If that search fails and `LIVER_PIPELINE_ROOT` is set, use that directory.
- If neither is available, ask the user for the project path before running commands.

## Quick Start

```bash
python scripts/liverbio.py docking pipeline --config config/docking_config.json
```

`liverbio docking` forwards every argument to `scripts/run_docking.py`, so docking subcommands such as `evidence`, `prepare-receptor`, `prepare-ligands`, `dock`, `analyze`, `redock`, `report`, `virtual-knockout`, `network`, and `faers` work through the same dispatcher.

## Operating rules

- Run `python scripts/liverbio.py docking init` before a fresh docking workspace when the workdir skeleton does not exist.
- Use the `pipeline` subcommand with checkpoint/resume support by default; use individual subcommands only when the user asks to run one stage.
- Pass existing inputs (receptor, ligand library, CSVs, config) as real paths and do not fabricate missing targets or docking results.
- Check `check-env` or `check-cadd` before reporting an environment problem; evidence collection depends on installed Codex database skills under `%USERPROFILE%\.codex\skills`.
- Do not edit docking source modules unless the user explicitly asks for code changes.

## Useful options

The CLI uses subcommands, so run `python scripts/run_docking.py <subcommand> --help` to inspect a stage. `README.md` documents each subcommand and required inputs.

## References

- `README.md`: command reference and workflow description.
- `VIRTUAL_SCREENING_REQUIREMENTS.md`: external software and Codex skills required by CADD.
- `config/docking_config.json`: current docking defaults.
