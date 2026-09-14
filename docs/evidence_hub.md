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
| GTEx | liver expression context | context |
| Human Protein Atlas | liver expression context | context |
| DepMap local snapshot | cancer-cell dependency | context |

CTD, Tox21, ToxCast/CompTox, LINCS, DisGeNET and licensed resources such as
GeneCards, OMIM, TTD and DrugBank can be imported through the generic local
table connector. Licensed data are never downloaded or redistributed by the
repository.

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

## Scientific limits

- Open Targets aggregates many underlying sources. Do not add its underlying
  databases again as independent evidence without source-level de-duplication.
- Predicted and experimental evidence remain separate categories.
- Target coverage depends on correct gene, protein, compound and disease
  identifiers. Identifier harmonisation is recorded but does not replace manual
  review.
- Ranking stability is not biological validation. Publication-level conclusions
  still require endpoint-matched external cohorts and orthogonal experiments.
