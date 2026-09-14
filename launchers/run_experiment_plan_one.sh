#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python "$ROOT/scripts/run_experiment_plan_one.py" \
  --config "$ROOT/config/experiment_plan_one.json" \
  --output-root "$ROOT/../experiment_plan_one_results" "$@"
