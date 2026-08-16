#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FL_QOM_ROOT="${FL_QOM_ROOT:-/home/nanaengo/FL_QOM}"
ENV_FILE="${ROOT}/environment.yml"

command -v mamba >/dev/null 2>&1 || {
  echo "ERROR: mamba is required; install it before running this script." >&2
  exit 2
}
[[ -d "$FL_QOM_ROOT" ]] || {
  echo "ERROR: FL_QOM root does not exist: $FL_QOM_ROOT" >&2
  exit 2
}

if ! mamba env list | awk '{print $1}' | grep -qx qom; then
  echo "Creating qom from ${ENV_FILE}"
  mamba env create --file "$ENV_FILE"
else
  echo "qom already exists; updating it from ${ENV_FILE}."
  mamba env update --name qom --file "$ENV_FILE"
fi

mamba run -n qom python "${ROOT}/scripts/validate_environment.py" \
  --fl-qom-root "$FL_QOM_ROOT" \
  --output "${ROOT}/results/environment_manifest.json"

# PyQt5 is an optional GUI dependency and is intentionally excluded from the
# headless HPC validation; all numerical/physics tests remain included.
mamba run -n qom python -m pytest "$FL_QOM_ROOT/tests" -q \
  --disable-warnings --maxfail=1 --ignore="$FL_QOM_ROOT/tests/test_gui.py"

echo "Environment verification completed."
