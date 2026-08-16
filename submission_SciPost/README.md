# SciPost Physics — submission package (Track D)

This directory is a **format conversion** of the canonical AIP/Chaos manuscript
(`../manuscript/Chaos_Floquet_SyntheticFlux.tex`) to the official **SciPost
Physics** template, following the editorial decision to target SciPost Physics
(IF 5.29 in 2024; Genuine Open Access — **no APC, no subscription fees**, CC BY 4.0).

## Files

| File | Role |
|---|---|
| `SciPost_Chaos_Floquet_SyntheticFlux.tex` | Main text (SciPost template) |
| `si_Chaos_Floquet_SyntheticFlux.tex` | Supplemental Material (SciPost template) |
| `SciPost.cls` | Official SciPost class (v2019-08, option `Phys`) |
| `SciPost_bibstyle.bst` | Official SciPost bibliography style |
| `Chaos_Floquet_SyntheticFlux.bib` | Shared bibliography (36 entries) |

Figures are resolved from `../figures/` (unchanged paths).

## Build (exact order — bidirectional `xr` cross-references)

```bash
cd submission_SciPost
pdflatex -interaction=nonstopmode si_Chaos_Floquet_SyntheticFlux.tex   # SI aux
bibtex   si_Chaos_Floquet_SyntheticFlux                                  # SI bibliography
pdflatex -interaction=nonstopmode si_Chaos_Floquet_SyntheticFlux.tex
pdflatex -interaction=nonstopmode SciPost_Chaos_Floquet_SyntheticFlux.tex
bibtex   SciPost_Chaos_Floquet_SyntheticFlux                             # main bibliography
pdflatex -interaction=nonstopmode SciPost_Chaos_Floquet_SyntheticFlux.tex
pdflatex -interaction=nonstopmode SciPost_Chaos_Floquet_SyntheticFlux.tex
pdflatex -interaction=nonstopmode si_Chaos_Floquet_SyntheticFlux.tex    # SI reads MAIN aux
```

Verified state: main **19 pp, 0 errors, 0 undefined refs, 36 citations**; SI **6 pp, 0 errors, 0 undefined, 9 citations** (the SI carries its own bibliography: anchors Mayor2025/Mathew2020 plus the methodological references).

## What changed vs the AIP version

1. **Documentclass** `revtex4-2 (aip,cha,reprint)` → `SciPost (submission, Phys)`.
2. **Title/author/affiliation** → SciPost centred blocks (no `\maketitle`/`\keywords`).
3. **Abstract** rewritten for the SciPost "context → problem → methods → results →
   conclusions → outlook" structure (156 words, boldface, headline-first).
4. **Table of contents** added (paper > 6 pages, per template guideline).
5. **Introduction** strengthened with an explicit statement of the *novel synergetic
   link* across synthetic-gauge physics, nonlinear dynamics, and quantum-limited
   measurement (this is the expectation the paper is argued to satisfy).
6. **Supplemental Material bibliography** added: the two anchors (Mayor2025,
   Mathew2020) are now proper `\cite` references (previously inline text with
   DOIs), and the methodological sources (QR/Benettin, Cram\'er--Rao/measurement,
   master-equation noise model) are cited — 9 SI citations total.
7. Two latent LaTeX bugs fixed in **both** the converted and canonical files:
   an unbalanced `$--$` in the ensemble-size bootstrap range, and unescaped
   underscores in two file paths in the SI.

## Acceptance-criteria mapping (SciPost Physics)

SciPost Physics requires **at least one** *Expectation* and **all** *General
acceptance criteria*.

| Expectation | Status |
|---|---|
| Novel and synergetic link between different research areas | ✅ **Primary argument** — synthetic gauge / nonlinear dynamics (hyperchaos) / optomechanics / quantum-limited sensing, jointly certified on one orbit |
| Open a new pathway with multi-pronged follow-up | ✅ the Floquet–Lyapunov consistency protocol is a reusable diagnostic; calibrated-device UQ is the named follow-up |

| General criterion | Status |
|---|---|
| Clear, jargon-free, unambiguous | ✅ (reframed abstract + intro) |
| Reproducible derivations in appendices | ✅ (SI carries the derivations) |
| Representative, complete citations | ✅ 36 refs (grown from 8) |
| Conclusion with objective reach/limitations + outlook | ✅ (existing Conclusions are explicit) |
| Detailed abstract + introduction | ✅ (rewritten) |
| Reproducibility resources (code/data in repository) | ⚠️ **OPEN** — scripts/manifests exist but a public release (Zenodo/GitHub) with DOI is still required before/at submission |

## Remaining author-controlled items

- Author list, affiliations, and corresponding-author email: **filled** (Stella
  Rolande Mbokop Tchounda, Carolle Tchodimou, Philippe Djorwe, Sifeu Takougang
  Kingni, Serge Guy Nana Engo — corresponding author, from the DQS_260518 archive).
- CRediT author contributions: **filled** (SRMT/CT simulations and analysis;
  PD/STK theory and feedback; SGNE conception, supervision, and writing).
- Competing interests: **filled** (none).
- Funding information: **filled** (no external funding).
- Data/code repository release with a DOI (the one substantive blocker above).
- The SI is currently a **separate document**; SciPost also accepts it merged as
  appendices — a merge is optional and can be done at production time.
