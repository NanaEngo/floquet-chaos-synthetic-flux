#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT/manuscript"
TEX=Chaos_Floquet_SyntheticFlux.tex
pdflatex -interaction=nonstopmode -halt-on-error "$TEX" > build_1.log 2>&1
bibtex Chaos_Floquet_SyntheticFlux > bibtex.log 2>&1
pdflatex -interaction=nonstopmode -halt-on-error "$TEX" > build_2.log 2>&1
pdflatex -interaction=nonstopmode -halt-on-error "$TEX" > build_3.log 2>&1
if grep -nE '^!|Emergency stop|Fatal error|undefined|multiply defined' build_3.log; then
  echo "ERROR: manuscript build has fatal or unresolved-reference diagnostics." >&2
  exit 1
fi
printf '%s\n' "Built $ROOT/manuscript/Chaos_Floquet_SyntheticFlux.pdf"
