# Multi-source target evidence hub

The evidence hub collects target associations into a local SQLite database,
preserves source-level provenance, distinguishes missing evidence from a
negative result, and produces coverage-aware target rankings.

## Why it exists

A union of target lists is not sufficient for publication-grade prioritisation.
Two databases may contain the same underlying evidence, predicted targets may be
ranked beside direct biochemical measurements, and a target absent from a
database is not evidence that the association is false.

The hub therefore stores one row per source assertion and evaluates:

- direct experimental compound-target evidence;
- curated compound-target evidence;
- predicted compound-target evidence;
- human genetic disease association;
- curated disease association;
- liver expression or cell-context evidence;
- functional dependency evidence;
- pathway evidence;
- experimental or predicted structure evidence.

## Canonical record

Every record contains:

```text
source, source_version, source_record_id
subject_type, subject_id, relation, object_type, object_id
target_symbol, evidence_type, tier
score, effect_size, p_value, sample_size
species, tissue, cell_type, assay, direction
retrieved_at, license, payload_json, fingerprint
```

The current open connectors are:

| Connector | Evidence family | Default tier |
|---|---|---|
| Open Targets | disease-target association | curated |
| ChEMBL | direct compound bioactivity | experimental |
| BindingDB | direct binding measurement | experimental |
| PubChem BioAssay | bioassay activity | experimental |
| GWAS Catalog | human genetic association | genetic |
| ClinVar | clinical variant association | genetic |
| GTEx | liver expression context | context |
| Human Protein Atlas | liver expression context | context |
| DepMap local snapshot | cancer-cell dependency | context |

CTD, Tox21, ToxCast/CompTox, LINCS, DisGeNET and licensed resources such as
GeneCards, OMIM, TTD and DrugBank can be imported through the generic local
table connector. Licensed data are never downloaded or redistributed by the
repository.

`config/evidence_local_sources.example.json` also contains disabled templates
for ClinVar, gnomAD, cBioPortal, OncoKB, CIViC and ClinicalTrials.gov exports.
When a user has the corresponding authorized or public snapshot, the template
can be copied into a local evidence config and enabled without changing code.

## Scoring

Scores are calculated from available categories only. Missing categories remain
`NaN`, are listed in `missing_categories`, and reduce `coverage_ratio`.

```text
priority_score =
    evidence_weight * weighted_mean(available category scores)
    + coverage_weight * coverage_ratio
```

Within a category, records are first collapsed to one maximum score per
`source_group`, preventing one source with many rows from dominating another
independent source. `source_ablation.csv` reports rank changes after removing
each source group.

The default category weights are in `config/evidence_sources.json`. They are
analysis assumptions, not universal biological constants. They must be
pre-specified or evaluated through sensitivity analysis before publication.

## CLI

```bash
python scripts/run_evidence_hub.py \
  --disease "metabolic dysfunction-associated steatotic liver disease" \
  --targets PNPLA3,TM6SF2,GPAT3,PPP2R2A \
  --output ../evidence_run
```

Compound-target evidence:

```bash
python scripts/run_evidence_hub.py \
  --pubchem-cid 154926030 \
  --smiles "CC(C)CC(C)NC1=CC(=O)C(=CC1=O)NC2=CC=CC=C2" \
  --disease NAFLD \
  --output ../evidence_run
```

The unified entry point is:

```bash
liverbio evidence --config config/evidence_sources.json \
  --disease NAFLD --targets GPAT3,PPP2R2A --output ../evidence_run
```

## Outputs

| File | Purpose |
|---|---|
| `evidence.sqlite` | reusable local evidence database |
| `evidence_records.csv` | one row per source assertion |
| `source_runs.csv` | connector status, timing, version and query hash |
| `evidence_coverage.csv` | coverage by source, tier and relation |
| `target_priority.csv` | target ranking with coverage and sensitivity metrics |
| `target_evidence_matrix.csv` | target-by-category score matrix |
| `source_ablation.csv` | rank stability after source removal |
| `evidence_run_manifest.json` | command context, software and store summary |
| `target_priority.md` | readable top-target table |

## Experiment-plan-one integration

`config/experiment_plan_one.json` contains an optional `evidence` stage between
`disease` and `ppi`. It writes results under `02b_evidence/`.

```json
{
  "evidence": {
    "enabled": true,
    "config_file": "config/evidence_sources.json",
    "target_scope": "intersection",
    "allow_network": true
  }
}
```

Set `allow_network` to `false` for an offline run using cached/local sources.
The stage records each connector failure and continues unless `strict=true`.
Set `LIVER_CONTACT_EMAIL` before large automated runs so remote services can
contact the operator about API traffic.

## Full-pipeline integration

The full automated pipeline now uses the evidence hub in stage 03. Stage 02
writes `candidate_universe.csv` before the legacy `key_genes.csv` Top-N view.
Stage 03 collects evidence for `evidence.max_targets` candidates, writes the
hub to `outputs/integration/evidence_hub/`, and produces
`outputs/integration/target_priority.csv`. Stage 05 merges the heuristic
perturbation score into `integrated_target_priority.csv`, and stage 06 uses
that integrated ranking before falling back to `target_priority.csv` or
`key_genes.csv`.

Set `evidence.hub_config`, `evidence.disease_name`, `evidence.max_targets`,
and `evidence.allow_network` in `config/full_pipeline_config.json`. Missing
evidence remains `not_found` or `not_queried` in the target priority table and
is never converted to a negative result. The integrated report also writes a
`publication_readiness` block that distinguishes exploratory results from
paper-supporting or publication-grade evidence gates.

The full pipeline now keeps a DEG-derived `candidate_universe.csv` and writes
an evidence-expanded `candidate_universe_evidence_expanded.csv`. Disease and
genetic targets discovered by Open Targets or GWAS Catalog are unioned into the
expanded universe, controlled by `evidence.candidate_expansion_max_targets`.
Legacy target-level ChEMBL, PDB, AlphaFold, Reactome and KEGG evidence is
converted into the canonical evidence store before scoring. This makes the
chemical, structural and pathway axes available even when no particular
compound was supplied by the user.

Positive and negative benchmark targets can be supplied through
`evidence.benchmark_positive_targets`, `evidence.benchmark_negative_targets`
and `evidence.benchmark_top_n`, or through the corresponding CLI options.
`evidence.strict=true` now propagates failures to the full-pipeline stage
instead of being silently downgraded to a warning. A changed local evidence
file also changes the evidence-store fingerprint, so stale records are not
silently reused. If a configured source fails in the current run while strict
mode is off, its previous records are removed before scoring so that a stale
success cannot masquerade as current evidence.

## Scientific limits

- Open Targets aggregates many underlying sources. Do not add its underlying
  databases again as independent evidence without source-level de-duplication.
- Predicted and experimental evidence remain separate categories.
- Target coverage depends on correct gene, protein, compound and disease
  identifiers. Identifier harmonisation is recorded but does not replace manual
  review.
- Ranking stability is not biological validation. Publication-level conclusions
  still require endpoint-matched external cohorts and orthogonal experiments.
