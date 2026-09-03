---
name: liver-dataset-search
description: Search and download NCBI GEO, EBI BioStudies/ArrayExpress, or Expression Atlas datasets with the Liver Cancer Bioinformatics project dataset-search CLI. Use only inside that project or when LIVER_PIPELINE_ROOT points to it.
---

# Liver Dataset Search

Use the project's dataset search script; do not hand-build GEO or EBI queries from memory.

## Locate the project

- Use the current working directory or its nearest ancestor containing `AGENTS.md` and `scripts/search_datasets.py`.
- If that search fails and `LIVER_PIPELINE_ROOT` is set, use that directory.
- If neither is available, ask the user for the project path before running commands.

## Quick Start

```bash
python scripts/liverbio.py datasets \
  --disease "liver cancer" \
  --research-direction "single cell RNA-seq" \
  --max-results 20
```

`liverbio datasets` forwards arguments to `scripts/search_datasets.py`.

## Operating rules

- Require at least one of `--query`, `--disease`, or `--research-direction`; the CLI rejects otherwise.
- Save results to the default `data_cache/dataset_search/` output unless the user asks for another directory.
- Treat `run_supported` as a candidate marker, not a guarantee: file availability is verified only during download.
- Use `--download <GSE,...>` or `--download-top N` only when the user explicitly wants data downloaded.
- Do not invent accession counts, sample counts, relevance scores, or download success.

## Useful options

The CLI supports organism, data type, sample count, date, platform, and dataset-type filters plus optional ML reranking. Read `python scripts/search_datasets.py --help` before adding options the user did not mention.

## References

- `README.md`: dataset-search workflow and download details.
- `docs/software_guide.md`: end-user suite guide.
