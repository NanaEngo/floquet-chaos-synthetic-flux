#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"
export PYTHONPATH="${ROOT}/scripts:${PYTHONPATH:-}"
mkdir -p "$ROOT/results"
# Production defaults: 9 phases x 3 replicates, long QR window, 24 CPU workers.
set +e
"$PYTHON_BIN" "$ROOT/scripts/flux_grid.py" --output "$ROOT/results/flux_grid.json" --points "${FLUX_POINTS:-9}" --replicates "${FLUX_REPLICATES:-3}" --n-steps "${GRID_LYAPUNOV_STEPS:-2000000}" --transient-steps "${GRID_TRANSIENT_STEPS:-200000}" --workers "${GRID_WORKERS:-24}" --backend python --max-rate-difference 5e-4 --max-divergence-residual 5e-3 --max-block-std 5e-3
grid_rc=$?
set -e
if [ "$grid_rc" -ne 0 ]; then
  echo "Flux grid did not pass its strict spectrum gate; retaining FAIL records and continuing to the independent covariance/Fisher stage."
fi
"$PYTHON_BIN" "$ROOT/scripts/covariance_fisher.py" --output "$ROOT/results/covariance_fisher.json"
"$PYTHON_BIN" "$ROOT/scripts/generate_extended_figures.py" --results "$ROOT/results/covariance_fisher.json" --output-dir "$ROOT/figures"
"$PYTHON_BIN" "$ROOT/scripts/hash_sources.py" --root "$ROOT" --output "$ROOT/results/source_hashes.json"
"$PYTHON_BIN" "$ROOT/scripts/finalize_run.py" --results "$ROOT/results" --workspace "$ROOT"
echo "Extended validation completed."
