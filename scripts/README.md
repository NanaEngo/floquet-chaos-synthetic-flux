# Reconstruction scripts

The scripts in this directory are the only approved entry points for the new
Chaos reconstruction. They use the six-dimensional two-mechanical-mode model in
`reconstruction_core.py`, with normalized units declared in the internal
mathematical derivations.

## Order of execution

```bash
# after activating fl_qom
bash scripts/run_reconstruction.sh
```

The pipeline executes:

1. `model_gate.py` — analytic Jacobian versus central finite differences, divergence, and gauge-invariant flux;
2. `floquet_lyapunov_gate.py` — drive-locked periodic orbit, monodromy, Floquet rates, QR/Benettin spectrum, and divergence balance;
3. `generate_figures.py` — figures only if the JSON result has `status=PASS`;
4. `finalize_run.py` — provenance files and checksums.

The extended stages are run with `bash scripts/run_extended_validation.sh`; production defaults are 9 phases, 3 replicates, 2,000,000 integration steps, 200,000 transient steps, and 24 CPU workers:

- `flux_grid.py` — prespecified long-window flux grid with three independent initial-condition replicates, deterministic multiprocessing, and strict fail-closed thresholds;
- `covariance_fisher.py` — periodic covariance and classical Fisher information for the declared cavity-quadrature record;
- `generate_flux_grid_figure.py` — provenance-bound grid diagnostics figure, generated only from a `PASS` JSON;
- `lyapunov_numba.py` — optional equivalence-tested Numba backend; the validated production result uses the reference Python/NumPy backend because its small QR blocks are faster here.

Supporting scripts:

- `validate_environment.py` records Python/dependency versions, FL_QOM revision, and the pyproject hash.
- `validate_physical_manifest.py` fails closed until SI calibration, admissible ranges, objective, constraints, and uncertainties are supplied; legacy candidate values are not accepted automatically.
- `pareto_robust_optimizer.py` implements the robust Pareto design. It requires the accepted manifest by default; `--allow-provisional-literature-scenario` is an explicit exploratory override and labels output `PROVISIONAL`. Use `--initial-condition-replicates >= 3` for the robustness gate; a failed orbit under any replicate is retained as a failure, not silently dropped.
- `convert_drive_calibration.py` converts an explicitly defined coherent input amplitude to input power and normalized additive drive. It requires all SI inputs and the `coherent_input_amplitude_sqrt_photons_per_second` convention; incomplete or quadrature-phase input returns `NOT_COMPUTED`.
- `generate_pareto_figure.py` creates a diagnostic figure only from the provisional machine-generated JSON and never feeds the manuscript automatically.
- `bootstrap_qom.sh` creates/verifies the `fl_qom` environment, installs FL_QOM editable, runs the upstream tests, and records the environment manifest.
- `reconstruction_core.py` contains the equations and numerical primitives; it does not write results itself.

## Failure policy

A missing dependency, failed model gate, non-converged periodic orbit, mismatch
between Floquet and Lyapunov diagnostics, or failed divergence balance stops the
pipeline. No figure generator may fall back to mock, representative, or manually
entered data.

All output JSON files belong in `results/` and all generated figures belong in
`figures/`. They must be accompanied by provenance and should be committed only
after scientific review of the actual values. The provisional screening command and
limits are recorded in `results/physical_parameter_run_manifest.json`.
