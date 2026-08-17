# Public code and data deposit

**Paper:** *Floquet–Lyapunov consistency and matched sensing bounds for
flux-controlled hyperchaos in dissipative optomechanics*

**Authors:** Stella Rolande Mbokop Tchounda, Carolle Tchodimou,
Philippe Djorwe, Sifeu Takougang Kingni, Serge Guy Nana Engo
(corresponding author).

**Target venue:** SciPost Physics (diamond open access).

This repository is the public, reproducible code and data package for the
manuscript. Every numerical claim in the paper traces to a deterministic script
in `scripts/` and a machine-readable artifact in `results/`; every figure in
`figures/` is generated from those artifacts and carries a sidecar
`*.manifest.json` provenance record.

## Contents

| Directory | Contents |
|---|---|
| `scripts/` | Deterministic generators, validators, and figure scripts (Python). |
| `results/` | Machine-readable JSON outputs, run manifests, `checksums.sha256`, `source_hashes.json`. |
| `figures/` | Vector (PDF) and raster (PNG) figures plus `*.manifest.json` provenance. |
| `submission_SciPost/` | SciPost Physics submission package — the canonical manuscript (main + SI, class, bib style) and compiled PDFs. |
| `environment.yml` | Conda environment specification (`fl_qom`). |

## Environment

The environment is specified in `environment.yml`:

```bash
mamba env create -f environment.yml
mamba activate fl_qom
```

The Floquet infrastructure package **FL_QOM** is a separate repository and is
required by the scripts. Install it from its own repository before running the
generators:

```bash
pip install -e /path/to/FL_QOM
```

## Reproducing the headline results

All commands are run from the repository root with the `fl_qom` environment active
(`mamba run -n fl_qom python ...`). The key generators are:

| Claim | Script → artifact |
|---|---|
| Model/Jacobian/trace/flux gate | `scripts/model_gate.py` → `results/model_gate.json` |
| Floquet–Lyapunov consistency (reference orbit) | `scripts/floquet_lyapunov_gate.py` → `results/floquet_lyapunov.json` |
| Long-window flux grid (27 runs) | `scripts/flux_grid.py` → `results/flux_grid.json` |
| Full-spectrum weak-coupling transition map | `scripts/full_spectrum_transition_map.py` → `results/full_spectrum_transition_map.json` |
| Strong-coupling hyperchaos transition map | `scripts/full_spectrum_hyperchaos_map.py` → `results/full_spectrum_hyperchaos_map_g03.json` |
| Hyperchaos-order seed robustness | `scripts/hyperchaos_seed_robustness.py` → `results/hyperchaos_seed_robustness.json` |
| Finite-time convergence probe | `scripts/hyperchaos_convergence_probe.py` → `results/hyperchaos_convergence_probe.json` |
| Drive-axis bifurcation continuation (Neimark–Sacker) | `scripts/drive_axis_bifurcation_scan.py` → `results/drive_axis_bifurcation_scan.json` |
| Torus observation | `scripts/torus_observation.py` → `results/torus_observation.json` |
| Attractor structure (unstable seed vs attractor) | `scripts/attractor_structure.py` → `results/attractor_structure.json` |
| Correlation dimension (Grassberger–Procaccia) | `scripts/correlation_dimension.py` → `results/correlation_dimension.json` |
| Covariance/Fisher pilot | `scripts/covariance_fisher.py` → `results/covariance_fisher.json` |
| Matched measurement-Fisher reference (force sensing) | `scripts/matched_fisher_reference.py` → `results/matched_fisher_reference.json` |
| Matched-Fisher UQ (Monte Carlo) | `scripts/matched_fisher_uq.py` → `results/matched_fisher_uq.json` |
| Noise-dependent observability | `scripts/noise_observability.py` → `results/noise_observability.json` |
| Amplitude-quadrature noise cross-check (master equation) | `scripts/langevin_master_crosscheck.py` → `results/langevin_master_crosscheck.json` |
| Semiclassical validity at the reference point | `scripts/semiclassical_validity.py` → `results/semiclassical_validity.json` |
| Strong-coupling reachability (2542× bound) | `scripts/strong_coupling_reachability.py` → `results/strong_coupling_reachability.json` |
| Reduced-model (adiabatic) Pareto screen | `scripts/pareto_robust_optimizer.py` → `results/physical_parameter_reduced_bad_cavity_pareto_v4.json` |

## Figures

Each figure is generated from its `results/` artifact by the matching
`scripts/plot_*.py` script. The figure list and the mapping to the manuscript
figures are in `figures/README.md` and the individual `figures/*.manifest.json`
records.

## Integrity and provenance

- `results/checksums.sha256` — SHA-256 of every result file.
- `results/source_hashes.json` — hashes of the generator scripts.
- `results/run_manifest.json` and the per-analysis `*_run_manifest*.json` —
  command line, environment, and design metadata.

## Honest scope of the claims

This package reproduces a model-level, numerically validated study. It is **not**
a demonstration of a quantum advantage, an SQL-beating sensitivity, Wigner
negativity, topological protection, or an exceptional-point sensing gain. The
matched force-sensing gain is order-unity (median ≈ 1.02–1.04; peak ≤ 1.32×),
and the strong-coupling hyperchaos regime lies ≈2542× beyond the anchored
optomechanical coupling, so it is reported as a model-level demonstration rather
than a device prediction.

## License

Code: MIT License. Data and figures: CC BY 4.0 (consistent with SciPost's
open-access policy). See `LICENSE` and `LICENSE-DATA`.

## Citation

A Zenodo DOI is minted upon publication. Until then, cite the manuscript and
this repository URL.
