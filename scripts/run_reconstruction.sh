#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS="${ROOT}/results"
FIGURES="${ROOT}/figures"
PYTHON_BIN="${PYTHON_BIN:-python}"

mkdir -p "$RESULTS" "$FIGURES"
export PYTHONPATH="${ROOT}/scripts:${PYTHONPATH:-}"

"$PYTHON_BIN" "${ROOT}/scripts/model_gate.py" --output-dir "$RESULTS"
"$PYTHON_BIN" "${ROOT}/scripts/floquet_lyapunov_gate.py" \
  --output-dir "$RESULTS" \
  --n-steps "${LYAPUNOV_STEPS:-2000000}" \
  --transient-steps "${LYAPUNOV_TRANSIENT:-200000}" \
  --dt "${LYAPUNOV_DT:-0.01}" \
  --qr-interval "${LYAPUNOV_QR_INTERVAL:-10}"
"$PYTHON_BIN" "${ROOT}/scripts/generate_figures.py" \
  --results "${RESULTS}/floquet_lyapunov.json" \
  --output-dir "$FIGURES"
"$PYTHON_BIN" "${ROOT}/scripts/finalize_run.py" \
  --results "$RESULTS" \
  --workspace "$ROOT"

echo "Reconstruction pipeline completed with PASS status."
