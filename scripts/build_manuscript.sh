#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT/submission_SciPost"
MAIN=SciPost_Chaos_Floquet_SyntheticFlux.tex
SI=si_Chaos_Floquet_SyntheticFlux.tex
# Build the SI first (its aux feeds the main's external xr references), then the
# main, then repeat so the bidirectional cross-references resolve.
pdflatex -interaction=nonstopmode -halt-on-error "$SI" > build_si_1.log 2>&1
bibtex si_Chaos_Floquet_SyntheticFlux > build_si_bib.log 2>&1 || true
pdflatex -interaction=nonstopmode -halt-on-error "$SI" > build_si_2.log 2>&1
pdflatex -interaction=nonstopmode -halt-on-error "$MAIN" > build_main_1.log 2>&1
bibtex SciPost_Chaos_Floquet_SyntheticFlux > build_main_bib.log 2>&1 || true
pdflatex -interaction=nonstopmode -halt-on-error "$MAIN" > build_main_2.log 2>&1
pdflatex -interaction=nonstopmode -halt-on-error "$MAIN" > build_main_3.log 2>&1
pdflatex -interaction=nonstopmode -halt-on-error "$SI" > build_si_3.log 2>&1
if grep -nE '^!|Emergency stop|Fatal error|undefined|multiply defined' build_main_3.log build_si_3.log; then
  echo "ERROR: manuscript build has fatal or unresolved-reference diagnostics." >&2
  exit 1
fi
printf '%s\n' "Built $ROOT/submission_SciPost/SciPost_Chaos_Floquet_SyntheticFlux.pdf and si_Chaos_Floquet_SyntheticFlux.pdf"
