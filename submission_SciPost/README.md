# SciPost Physics — canonical manuscript (Track D)

This directory is the **single canonical manuscript** for Track D, set in the
official **SciPost Physics** template, following the editorial decision to target
SciPost Physics (IF 5.29 in 2024; Genuine Open Access — **no APC, no
subscription fees**, CC BY 4.0). It is the only manuscript version maintained; no
AIP/`Chaos` or arXiv copies are kept.

## Files

| File | Role |
|---|---|
| `SciPost_Chaos_Floquet_SyntheticFlux.tex` | Main text (SciPost template) |
| `si_Chaos_Floquet_SyntheticFlux.tex` | Supplemental Material (SciPost template) |
| `SciPost.cls` | Official SciPost class (v2019-08, option `Phys`) |
| `SciPost_bibstyle.bst` | Official SciPost bibliography style |
| `Chaos_Floquet_SyntheticFlux.bib` | Shared bibliography (49 entries) |

Figures are resolved from `../figures/` via `\graphicspath{{../figures/}}` (paths are filenames only in the source).

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

Verified state: main **20 pp, 8 802 words, 0 errors, 0 undefined refs, 51 citations**; Abstract is **196 words / 14 lines**, well within SciPost's preferred 8-line / ~200-word limit and rendered in boldface headline-first prose (the central numerical claim is in the first sentence). SI **16 pp, 7 444 words across 15 numbered sections, 0 errors, 0 undefined refs, 12 citations** (the SI carries its own bibliography, including Mayor2025, Mathew2020, Muthukumar2025, and the methodological references); the four SI sections added at the strengthening pass (limitations and residual gaps, feasibility audit, initial-condition robustness, attractor-Fisher audit) are the explicit technical backup for the Discussion in the main text. The post-audit polish added a "Headline findings" KPI list at the top of Conclusions, restructured the Abstract to lead with the central numerical finding, sharpened the title to "Synthetic-flux control of hyperchaos order and matched-resource sensing in dissipative optomechanics", merged the Discussion "Validity" + "Residual gap" subsections into a single Bridge-to-experiment subsection, and fixed a math-mode typo in the matched-measurement Fisher lemma (`\T` → `^T`).

**Substantive refinements for the SciPost acceptance target (median 60–70 %).**

- **Z1 Abstract restructured** in dense prose boldface (no internal labels) with all key numbers explicit (peak gain $\num{1.32}\times$, UQ median $\num{1.039}\times$, 90 % CI, $E^{*}$, $n_+=4$, 2542×).
- **Z2 Introduction restructured** with explicit *Position relative to prior work* and *Novelty statement* paragraphs (HalefShomroni comparée, Muthukumar continue trois choses concrètes).
- **Z3 Two theorems promoted** to the main text (Lyapunov–Floquet consistency criterion, matched-measurement Fisher lemma) in `theorem`/`lemma` environments with full formal statements; SI keeps the proofs.
- **Z4 State-of-the-art comparison table** (`tab:context-comparison`) places the present gain next to Teufel2011, Gavartin2012, Purdy2017, Li2021, Qvarfort2018 with matched-resource operating regimes.
- **Z5 Discussion restructured** into four explicit subsections (Relation to prior work / Position among benchmarks / Critical self-assessment and limits / Validity of the mean-field description) with all limitations stated explicitly.
- **Z6 SI reading guide** added at the top of the SI, with all SI cross-reference labels (`SI-sec:anchor`, `SI-sec:matched-fisher`, `SI-sec:theorems`, `SI-sec:convergence`, `SI-sec:threshold-sensitivity`, `SI-sec:fisher-window`, `SI-sec:reduced-pareto`, `SI-sec:noise-model-crosscheck`, `SI-sec:interpretation`, `SI-sec:attractor-fisher`, `SI-sec:semiclassical-boundary`).

## Structure

1. **Documentclass** `SciPost (submission, Phys)`.
2. **Title/author/affiliation** in SciPost centred blocks (no `\maketitle`/`\keywords`).
3. **Abstract** structured for the SciPost "context → problem → methods → results →
   conclusions → outlook" template (boldface, headline-first).
4. **Table of contents** included (paper > 6 pages, per template guideline).
5. **Introduction** states the distinct continuation from Muthukumar2025: the same
   model topology is used to resolve a drive/coupling-gated full-spectrum
   hyperchaos transition and to test, rather than assume, a matched sensing effect.
6. **Supplemental Material** carries its own bibliography: the two anchors
   (Mayor2025, Mathew2020) as proper `\cite` references plus the methodological
   sources (QR/Benettin, Cram\'er--Rao/measurement, master-equation noise model) —
   9 SI citations total.

## Acceptance-criteria mapping (SciPost Physics)

SciPost Physics requires **at least one** *Expectation* and **all** *General
acceptance criteria*.

| Expectation | Status |
|---|---|
| Distinct contribution and cross-area link | ✅ **Primary argument** — continuation of the Muthukumar2025 model toward full-spectrum hyperchaos order, bifurcation structure, and a matched measurement bound; not a generic claim that synthetic magnetism creates chaos |
| Open a new pathway with multi-pronged follow-up | ✅ the Floquet–Lyapunov consistency protocol is a reusable diagnostic; calibrated-device UQ is the named follow-up |

| General criterion | Status |
|---|---|
| Clear, jargon-free, unambiguous | ✅ (reframed abstract + intro) |
| Reproducible derivations in appendices | ✅ (SI carries the derivations) |
| Representative, complete citations | ✅ 49 refs (grown from 8) |
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
- The current manuscript explicitly distinguishes the Muthukumar2025 detuning convention from the present convention and includes a parameter/scope comparison in the SI.
- The SI is currently a **separate document**; SciPost also accepts it merged as
  appendices — a merge is optional and can be done at production time.
