# Chaos manuscript draft

`Chaos_Floquet_SyntheticFlux.tex` is the first scientific-article draft for the
isolated reconstruction workspace. It is deliberately written as an article—not as
a job report—and is centred on the question of whether gauge-invariant synthetic
flux produces a reproducible, noise-resolved change in dissipative dynamics.

## Current status

`INTERNAL DRAFT — NOT SUBMISSION READY`

The Results section contains the validated stable-orbit pilot, the three-phase
covariance/Fisher pilot, the long-window production grid, and a bounded theoretical
bad-cavity sensitivity screen. The grid covers nine phases and three independent
initial conditions per phase; all 27 strict Floquet--Lyapunov gates pass, with no
positive exponent in the weakly coupled pilot grid. The full-spectrum
classification was then extended to strong coupling (g = 0.3): a
drive- and flux-controlled transition to hyperchaos (up to four positive
exponents, Kaplan--Yorke dimension up to 4.8) is resolved and shown in
`figures/hyperchaos_transition_map.{pdf,png}` (`fig:hyperchaos-map`). The
transition is drive- and coupling-gated; the synthetic-flux phase modulates the
hyperchaos order only in the strong-coupling, strong-drive regime. The reduced screen reports a
16-candidate, three-point local fixed-point Pareto pattern under an explicitly
assumed literature-informed closure. It is not a physical calibration or optimum. A provisional finite-time
basin-capture validation (three non-dominated candidates, five independent initial
conditions per parameter replica, vectorized RK4 independently compared with
DOP853) has since passed for all three candidates
(`fast_basin_capture_validation_v2.json`, status
`PROVISIONAL_REDUCED_MODEL_BASIN_CAPTURE_VALIDATION`); this is a finite-time,
finite-ensemble check, not a global basin-of-attraction proof. The draft remains internal because
no matched-reference sensing gain has been established and declarations and final
venue checks remain to be completed. The eight cited references have now been checked against the publisher records
(Crossref) and stored with complete journal metadata in
`Chaos_Floquet_SyntheticFlux.bib`.

## AIP journal format (converted 2026-08-16)

Both documents now use the official AIP submission class
`\documentclass[aip,cha,reprint]{revtex4-2}` (the class underlying the AIP/Overleaf
submission template for *Chaos*), with `aipnum4-2.bst` reference style, native
`\keywords{}`, and `alt={...}` accessibility text on every figure. The
affiliation and corresponding-author email are placeholders to be completed by
the authors before submission.

## Build

Both documents cross-reference each other through `xr`, so build the SI first,
then the main text, then repeat (and rerun bibtex on the main text):

```bash
cd manuscript
pdflatex -interaction=nonstopmode si_Chaos_Floquet_SyntheticFlux.tex
pdflatex -interaction=nonstopmode Chaos_Floquet_SyntheticFlux.tex
bibtex Chaos_Floquet_SyntheticFlux
pdflatex -interaction=nonstopmode Chaos_Floquet_SyntheticFlux.tex
pdflatex -interaction=nonstopmode Chaos_Floquet_SyntheticFlux.tex
pdflatex -interaction=nonstopmode si_Chaos_Floquet_SyntheticFlux.tex
pdflatex -interaction=nonstopmode Chaos_Floquet_SyntheticFlux.tex
```

Known benign REVTeX warnings: "float is stuck" / "Deferred float stuck during
\clearpage" notices for the two page-6 figures (they are placed correctly) and
microtype's footnote-patch notice; the SI's `[h]` tables are converted to `[ht]`.
The bibliography is formatted by `aipnum4-2.bst` (AIP style) from
publisher-verified metadata; DOIs are retained for machine resolution. The draft
contains no numerical result that is not yet bound to an immutable result
manifest. The provisional finite-time basin-capture validation
(`PROVISIONAL_REDUCED_MODEL_BASIN_CAPTURE_VALIDATION`) is reported in the manuscript
as a finite-time, finite-ensemble check only and does not support a global
basin-of-attraction claim unless its joint criteria are met.
