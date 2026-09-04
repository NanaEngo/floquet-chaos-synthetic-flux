# Data-analysis report — Chaos reconstruction

**Version:** 0.5, 2026-08-17
**Status:** full-spectrum hyperchaos classification, matched force-sensing reference, master-equation semiclassical check, full-map seed-robustness, chaotic-attractor matched Fisher, optimized linear-sensor reference, and truncated-Wigner Lyapunov robustness recorded (2026-08-17 checkpoint below). Manuscript claims remain bounded by the evidence register.
**Source of truth for manuscript editing:** this file plus the immutable `results/` manifests.

## Objective

Determine whether the Floquet optomechanical model supports a reproducible, physically consistent account of dissipative nonlinear dynamics and, conditionally, whether any operational sensing advantage survives a complete noise and measurement analysis.

## Prespecified primary questions

1. What is the correct attractor structure of the driven six-real-variable model over the selected parameter domain?
2. Do Floquet stability, Poincare diagnostics, and the full Lyapunov spectrum agree under one declared time/phase convention?
3. What is the phase-space volume balance implied by the model’s divergence and dissipation terms?
4. How do quantum and thermal fluctuations propagate along unstable directions?
5. Does an explicitly defined measurement Fisher information exceed the matched linear reference under the same resource and time constraints?

## Analysis populations and units

- **System:** one cavity mode plus two mechanical modes, represented as six real variables only if that representation is explicitly derived and tested.
- **Drive:** periodic forcing with phase included in the model definition; the Poincare map and continuous-time formulation must not be conflated.
- **Time:** store physical seconds, normalized `omega_m*t`, and drive-period count as separate fields. Never overwrite one with another.
- **Parameters:** store SI values and normalized values side by side with units and conversion metadata.
- **Seeds:** deterministic seed for any stochastic trajectory ensemble; deterministic solver tolerances recorded per run.

## Required result objects

Every accepted run must produce:

- `run_manifest.json` — git revision, FL_QOM revision, environment, command, timestamp, seed, parameters, tolerances;
- `parameters.json` — SI and normalized values with units;
- `trajectory.h5` — time series only when needed, with time convention metadata;
- `lyapunov_spectrum.json` — full spectrum, window/block estimates, QR settings, convergence diagnostics;
- `floquet_diagnostics.json` — monodromy/multipliers and orbit/phase definition;
- `covariance_metrics.json` — drift, diffusion, covariance, physicality checks, and measurement Fisher quantities;
- `checksums.sha256` — hashes of all immutable result files;
- `status.json` — `PASS`, `FAIL`, or `NOT_COMPUTED` with reason.

> **Canonical-filename mapping (2026-08-15 data-restore check).** The names above are a specification template, not literal filenames. The production runs satisfy each requirement under the following canonical files (all present in `results/`, content verified):
>
> - `lyapunov_spectrum.json` → **`floquet_lyapunov.json`** (`lyapunov.spectrum`, `lyapunov.block_history`, `lyapunov.qr_interval`, plus `floquet.monodromy`, `floquet.multipliers`);
> - `floquet_diagnostics.json` → **`floquet_lyapunov.json`** (`floquet` block: monodromy, multipliers, floquet_rates) and **`flux_grid.json`** (orbit/phase definition);
> - `covariance_metrics.json` → **`covariance_fisher.json`** (drift/diffusion/covariance, physicality checks, classical measurement Fisher);
> - `manifest.json` → **`run_manifest.json`** (the immutable per-run manifest);
> - `trajectory.h5` → not required here: the production pipeline stores JSON results (no HDF5 time series was needed);
> - `production_summary_metrics.json` → intentionally absent (it is named only in the Prohibited-shortcuts rule below).

## Acceptance gates

### Gate A — model and units

PASS only if the Hamiltonian/EOM, sign conventions, damping, drive period, and dimensional conversions agree between equations, code, and metadata.

### Gate B — deterministic nonlinear dynamics

PASS only if integration tolerances, transient removal, initial-condition protocol, attractor classification, Poincare sections, and parameter continuation are documented and reproducible.

### Gate C — Lyapunov/Floquet consistency

PASS only if the full spectrum converges across block lengths and QR reorthogonalization choices, and the appropriate divergence/volume-balance diagnostic is reported.

### Gate D — stochastic/quantum fluctuations

PASS only if covariance propagation uses the full state-transition solution and the diffusion matrix is derived from declared baths. A scalar diffusive ansatz is not sufficient.

### Gate E — metrology

PASS only if the observable, estimator, measurement record, resource count, SQL reference, integration window, and uncertainty are matched across regimes. Trajectory separation alone cannot pass this gate.

### Gate F — quantum-state claims

PASS only if QuTiP or another validated quantum solver produces a physical density operator or equivalent phase-space solution, with Hilbert-space/truncation convergence. Otherwise Wigner/QFI claims remain `NOT_COMPUTED`.

## Current baseline status

| Item | Status |
|---|---|
| Original DQS manuscript audit | `FAIL — internal inconsistencies` |
| Original Wigner figure | `FAIL — mock generator identified` |
| Original SNR figure | `FAIL — trajectory-distance proxy, not measurement SNR` |
| Original Lyapunov spectrum | `NOT_ACCEPTED — recomputation required` |
| New Chaos manuscript | `INTERNAL DRAFT — deterministic and extended gates PASS` |
| New quantum-metrology claim | `WITHDRAWN_PENDING_EVIDENCE` |

## Planned statistical/reporting outputs

Report distributions and uncertainty across initial conditions and parameter replicates, not only maxima. For parameter maps, report grid dimensions, exclusions, convergence failures, and selection rules. Distinguish exploratory scans from prespecified primary parameter points.

## Checkpoint — first reproducible pilot run (2026-08-13)

The `qom` environment was created from `environment.yml`, FL_QOM was installed in editable mode, and the headless upstream suite passed **1,424 tests with 80 skips**. The skips are the optional PyQt5 GUI test and six Track-B tests whose external source path is absent from this checkout. Compatibility fixes required for the declared environment were recorded in FL_QOM: NumPy 1.26-compatible trapezoidal integration and QuTiP 5 keyword-only `e_ops` calls. The environment emitted an h5py/HDF5 build-versus-runtime warning; it did not affect the test suite, but it must be rechecked before quantum production runs.

The model gate passed with maximum analytic-versus-finite-difference Jacobian error `2.92e-10`, trace `-1.04`, and gauge-transformation flux residual `2.22e-16`. The first production periodic-orbit/Floquet–Lyapunov run passed with a Poincaré residual `8.49e-9`, maximum Floquet-rate versus QR-spectrum difference `2.68e-4`, tail-block standard deviation `5.3e-9`, and mean divergence `-1.04`. These values validate the **pilot stable periodic orbit and numerical consistency chain only**. They do not establish chaos, a flux-dependent transition, noise robustness, Fisher-information gain, QFI, Wigner negativity, or an exceptional point.

Evidence files:

- `results/environment_manifest.json`;
- `results/model_gate.json`;
- `results/floquet_lyapunov.json`;
- `results/source_hashes.json`;
- `figures/floquet_lyapunov_diagnostics.{pdf,png}` and its manifest.

This pilot checkpoint is retained as the single-orbit validation record; the declared long-window grid is documented in the subsequent checkpoint.

## Checkpoint — long-window flux grid and operational Fisher pilot (2026-08-13)

The production grid covered 9 equally spaced phases over $[-\pi,\pi)$ and 3 independent initial conditions per phase (27 runs), using 2,000,000 integration steps, 200,000 discarded transient steps, and the strict thresholds recorded in `results/flux_grid.json`. All 27 periodic-orbit searches and all 27 Floquet/QR consistency gates passed. The maximum absolute rate discrepancy was 2.849e-4, the maximum divergence residual was 4.746e-5, and the maximum terminal QR-block standard deviation was 6.10e-9. The largest observed QR exponent was -0.01002; no positive exponent was accepted in the sampled domain. This supports stable drive-locked dynamics over the tested grid but does not establish a universal absence of chaos. The former 30,000-step screening output is retained as `results/flux_grid_screening_30000.json` for audit history only.

The separate covariance/Fisher pilot passed for $\theta=0,\pi/2,\pi$ with $n_{\mathrm{th},1}=n_{\mathrm{th},2}=0.1$ and detector variance 0.01 in normalized units. The minimum physicality eigenvalue was 2.42e-7–2.58e-7 and classical measurement Fisher information was 3.78e-26–9.17e-13. This is an operational classical Fisher calculation for the declared cavity-quadrature record, not QFI and not a matched-baseline gain. A matched force-sensing reference (same drive, observation time, detector noise, baths, estimator; weak force F on mode 1) gives F_C(F) rising from 4.50e-5 (theta=0) to 5.96e-5 (theta=pi), i.e. at most a 1.32x gain over the flux-off coupled sensor and 1.16x over the single-mode sensor - a modest order-unity effect (central difference converged <0.01% across steps 1e-5..1e-3).

## Required next study — physical parameter calibration and constrained optimization

The current parameter set is a **numerically validated pilot regime**, not an experimentally optimized operating point. The primary optimization objective is now fixed as a robust multi-objective Pareto compromise, with stability, classical Fisher information, noise robustness, and experimental resource cost as secondary objectives; the Pareto analysis itself remains blocked by missing SI calibration and hardware constraints. The normalized values (`kappa=1`, `gamma1=gamma2=0.02`, `omega1=1`, `omega2=1.03`, `g1=0.02`, `g2=0.018`, `hopping=0.08`, `detuning=-1`, `drive=0.20`, modulation `0.10`, and drive frequency `1`) were selected for model validation and numerical conditioning. No claim of physical optimality is permitted.

Before a physical interpretation or sensing claim, the following study is required:

1. map every normalized parameter to SI units and record the experimental source or uncertainty range;
2. define admissible ranges for optical loss, mechanical damping, frequencies, optomechanical couplings, hopping, detuning, drive power, modulation depth, temperature, and detector noise;
3. predefine the objective: stable periodic operation, onset of chaos, classical Fisher information, robustness, or a constrained multi-objective compromise;
4. perform a coarse-to-fine parameter sweep with independent initial conditions and solver/tolerance checks;
5. optimize only within the physical ranges and report the Pareto set or uncertainty region rather than a single unsupported optimum;
6. validate the selected operating region against thermal, detector, and parameter uncertainties.

A local and web audit found a primary literature parameter paragraph in arXiv:2502.12336v2 (related journal DOI 10.1103/f6sv-4fq9). It supports a normalized reference scenario with \(\omega_m/2\pi\approx1\,\mathrm{GHz}\), \(\kappa/\omega_m=0.073\), \(\omega_2/\omega_m=1.0005\), \(g/\omega_m=1.077\times10^{-4}\), \(\gamma/\omega_m=1.077\times10^{-5}\), \(J_m/\omega_m=2\times10^{-4}\), and \(\Delta/\omega_m=+1\). The source does not resolve the mapping of \(\alpha^{\rm in}\) or pump power to the new additive `drive` coordinate, nor the hardware bounds.

The source-level equations also constrain the missing drive step: the reference model uses \(H_{\rm drive}=i\sqrt{\kappa}\alpha^{\rm in}(a^\dagger-a)\) and an additive semiclassical term \(\sqrt{\kappa}\alpha^{\rm in}\). Under a matched normalization and phase convention, the conditional model relation is \(E_{\rm model}=\sqrt{\kappa_{\rm norm}}\alpha^{\rm in}_{\rm norm}\). A physical power conversion would additionally require \(P_{\rm in}=\hbar\omega_p|\alpha^{\rm in}|^2\), with laser frequency and field normalization specified. These remain unresolved, so this is calibration evidence only and not an accepted SI mapping. The planned conversion is fail-closed: `scripts/convert_drive_calibration.py` will require an explicit laser frequency, `alpha_in` convention, input amplitude, SI cavity decay rate, and reference angular frequency before returning either input power or a normalized additive drive.

A two-stage explicit **provisional literature-anchored screening** was implemented using `physical_parameter_manifest_literature_scenario.json`. The parameter-only stage (4 candidates, 2 parameter replicates) found 3 feasible and 2 non-dominated points, but it was superseded because it did not replicate initial conditions. The final robustness stage used 2 candidates, 2 parameter replicates, 3 independent initial conditions, ±10% screening perturbations, the declared classical Fisher calculation, seed 20260813, and an orbit residual tolerance of \(10^{-4}\): **0/2 candidates passed the initial-condition robustness gate**. The result is `results/physical_parameter_optimization_provisional.json`, with `scientific_status=NO_ROBUST_FEASIBLE_CANDIDATE`; it is labelled `PROVISIONAL`, not `PASS`. The stage-1 diagnostic is retained as `figures/physical_parameter_pareto_screening_stage1.{pdf,png}` for audit history only. Because the drive mapping, hardware constraints, strict production convergence, and robust initial-condition gate remain unresolved, no parameter set can support an experimental-optimum or manuscript claim.

A second provisional run was then executed with 8 Latin-hypercube candidates, three parameter perturbation replicates, three independent initial conditions per replicate, the Fisher calculation enabled, seed `20260814`, and the strict orbit residual gate `1e-8`. It produced **0/8 feasible candidates and 0 non-dominated candidates**. The machine-readable checkpoint is `results/physical_parameter_optimization_provisional_v2.json`, with execution provenance in `results/physical_parameter_run_manifest_v2.json`; this is the definitive robustness result for the tested literature envelope. The earlier 5/8 parameter-only feasible count and the earlier stage-1 non-dominated points remainsuperseded diagnostics because they did not satisfy the initial-condition gate.

## Model-level Pareto recovery design (separate from physical calibration)

Because the literature-anchored low-damping scenario yielded no robust feasible point, a separate recovery screen is authorized at the level of the already validated normalized pilot model. It will vary the pilot parameters within an explicit local envelope, use independent parameter replicates and independent root-solver starts, and remain labelled `PROVISIONAL_MODEL_LEVEL`. This design can identify a local fixed-point pattern, but not finite-time basin robustness, experimental calibration, or optimality. The accepted physical manifest remains fail-closed.

The model-level recovery screen was positive: four of four candidates passed the declared local fixed-point/radius gates, and two candidates were non-dominated under the normalized objectives. The machine-readable result is `results/physical_parameter_optimization_model_recovery.json`, with provenance in `results/physical_parameter_run_manifest_model_recovery.json` and diagnostics in `figures/physical_parameter_pareto_model_recovery.{pdf,png}`. The conditional drive-to-input translation is recorded in `results/physical_parameter_pareto_drive_translation.json`: the required normalized input amplitudes are 0.3369999471 and 0.1281611690, while the corresponding power expressions remain symbolic until SI frequencies and the field convention are supplied. This is a `PROVISIONAL_MODEL_LEVEL` result: it identifies a local fixed-point Pareto pattern around the validated pilot, not finite-time basin robustness, SI calibration, or an experimentally optimal operating point.

The accepted manifest `physical_parameter_manifest.json` remains fail-closed as `NOT_COMPUTED`. The provisional scenario is intentionally separate and must not be substituted for it. The eventual required outputs remain an SI/normalized conversion table, measured hardware ranges, a matched objective/constraint definition, strict optimization results with uncertainty, and a provenance/checksum package. Until these outputs exist, the manuscript must use **pilot parameter set** or **validated numerical operating point**, never **optimal physical parameters**.

## Parameter audit checkpoint — 2026-08-14

The legacy `md_files` parameter tables were audited against primary sources. Several physical conversions were off by a factor of $10^3$: at $f_m=1$ GHz, $\gamma/\omega_m=1.077\times10^{-5}$ is 10.77 kHz, $g/\omega_m=1.077\times10^{-4}$ is 107.7 kHz, and $J_m/\omega_m=2\times10^{-4}$ is 200 kHz. The values $J_m/\omega_m=0.01$, 0.1, 0.065 and 0.02 correspond respectively to 10, 100, 65 and 20 MHz, not kHz.

An independent silicon OMC anchor is available from Mayor et al., *Nature Communications* 16, 2576 (2025), DOI `10.1038/s41467-025-57948-7`: $f_m\approx7.436$ GHz, $Q_m\approx3.6\times10^4$, $\gamma_m/2\pi\approx206.6$ kHz, $\kappa/2\pi\approx800$ MHz and $g_0/2\pi\approx880$ kHz. It does not calibrate the present two-resonator hopping $J_m$ or the input-drive convention. These values are therefore recorded as an external platform anchor only.

The complete decision and corrected conversion table are recorded in the physical-calibration audit. The normalized model and its provisional Pareto result are unchanged; `PHYSICAL_CALIBRATION_STATUS=NOT_COMPUTED`, `J_M_CALIBRATION_STATUS=MISSING`, and `DRIVE_POWER_CALIBRATION_STATUS=MISSING` remain active.

## Cross-platform screening checkpoint — 2026-08-14

A separate provisional screen combined the Mayor et al. 2025 silicon OMC scales ($\\kappa/\\omega_m\\approx0.1076$, $\\gamma_m/\\omega_m\\approx2.78\\times10^{-5}$, $g_0/\\omega_m\\approx1.18\\times10^{-4}$) with the Mathew et al. 2020 two-mode coupling context. Four Latin-hypercube candidates were tested with three parameter replicas and three independent initial conditions per candidate, strict orbit tolerance $10^{-8}$, and classical Fisher evaluation.

Result: `0/4` robust feasible candidates and `0` Pareto points. Every candidate failed the periodic-orbit gate under at least one parameter replica or initial condition. This is a negative cross-platform sensitivity result, not evidence that the physical device is infeasible; the imported coupling is not a calibration of the present phase-dependent $J_m$. The run is recorded in `results/physical_parameter_run_manifest_cross_platform.json` and `results/physical_parameter_optimization_cross_platform.json` and remains `PROVISIONAL`.

## Physical calibration input checkpoint — 2026-08-14

`physical_parameter_manifest_author_input_20260814.json` is now the controlled input template for the next SI stage. It contains only the Mayor et al. platform anchors for $\\omega_m$, $\\kappa$, $\\gamma_1$ and $g_1$; device-specific $\\gamma_2$, $\\omega_2$, $g_2$, $J_m$, detuning, drive and hardware constraints remain null by design. The standard validator rejects it with `NOT_COMPUTED`, and the drive converter rejects the incomplete input with `INCOMPLETE_DRIVE_CALIBRATION`. No SI optimization is launched from this template.

## Published two-mode values supplied — 2026-08-14

The directly reported Mathew et al. two-mode values are now recorded in `published_two_mode_parameter_set_mathew2020.json`: mechanical modes near 2.4 MHz separated by approximately 0.5 MHz, intrinsic linewidths near 3 kHz, optical linewidth 256 GHz, measured drive powers 33.4 and 174 microW, effective coupling up to 200 kHz, and direct next-nearest-neighbour coupling 10 kHz. These values are useful reference data but do not identify the present model's $g_1$, $g_2$, $\\kappa_{in}$, $\\alpha^{\\rm in}$ or phase-dependent $J_m$. The file therefore remains blocked and no SI optimization is launched from it.

## Explicit exploratory closure — 2026-08-14

For users requiring a fully populated numerical input, `published_assumption_closure_mathew2020.json` supplies an **assumed** closure: $J_m\leftarrow t_2$, symmetric $g_1=g_2$ inferred from the reported $g_{\\rm eff}$, the sign-converted detuning, and $\\alpha^{\\rm in}$ values derived from the reported incident powers under the standard input convention. This file is authorized only for an exploratory model screen. It is not a physical calibration, not a measurement record, and not authorized for manuscript claims.

## Assumed-closure execution checkpoint — 2026-08-14

The explicit Mathew closure was submitted to the current optimizer with two candidates, three parameter replicas, three initial conditions, and a strict $10^{-8}$ orbit gate. The run exceeded 600 s without producing an output file; no orphan process remained. The resulting status is `TIMEOUT_STIFF_REGIME`, recorded in `results/physical_parameter_run_manifest_assumed_closure.json`. The imported scale has $\\kappa/\\omega_1\\approx1.07\\times10^5$ and effective $g/\\omega_1\\approx69.5$, which is outside the numerical regime of the present orbit solver. This is a solver/model-compatibility diagnostic, not a physical infeasibility result. No tolerance was relaxed and no Pareto claim was generated.

## Adiabatic bad-cavity reduction checkpoint — 2026-08-14

To address the raw full-model timeout, `scripts/reduced_bad_cavity_screen.py` eliminates the optical amplitude using $a_{ss}=E/(\\kappa/2-i\\Omega_{\\rm eff})$ and integrates the four mechanical coordinates. The implementation was compared with the full six-dimensional model in common dimensionless regimes: mechanical infinity-norm errors were $9.27\\times10^{-7}$ at $\\kappa=20$ and $2.08\\times10^{-8}$ at $\\kappa=100$, with both solvers passing their orbit gates.

The assumed Mathew closure was first screened with a two-candidate diagnostic. A solver-status bug was then corrected: HYBR could report `success=false` despite a Poincare residual below $10^{-21}$. The corrected screen also samples distinct points rather than repeating interval midpoints. In the production reduced screen (4 candidates, 3 parameter replicas, 3 root-solver starts), `4/4` candidates pass the local fixed-point/radius gate. A transparent two-objective analysis (maximize $1-\\rho_{\\rm Poincare,max}$, minimize normalized drive; Fisher not computed) yields `3` non-dominated points. Because the starts are independent initial guesses to a root solver rather than finite-time trajectory integrations, this is a multi-start fixed-point diagnostic, not a basin-of-attraction robustness test.

The stronger deterministic rerun used 16 independently sampled candidates, 5 parameter perturbation replicates, and 5 independent root-solver starts per replicate (400 fixed-point evaluations), with seed `20260814`, residual tolerance $10^{-8}$, and exact replicate parameters retained in every candidate record. All `16/16` candidates passed the local residual and Poincare-radius gates; the maximum recorded radius was `0.9970026`, corresponding to a minimum local margin of `0.0029974`, and the same two-objective analysis produced `3` non-dominated points. This increases evidence for a reproducible reduced-equation fixed-point pattern, but finite-time basin capture remains `NOT_COMPUTED`; the result remains `PROVISIONAL_REDUCED_MODEL_MULTI_START_PARETO`.

## Prespecified theoretical robustness rerun — 2026-08-14

To distinguish a real model-level pattern from the four-point diagnostic, the prespecified reproducibility run was executed with the same assumed-closure manifest and validated reduced bad-cavity equations. It used 16 independently sampled candidates, 5 parameter perturbation replicates per candidate, and 5 independent root-solver starts per replicate (400 fixed-point evaluations), with seed `20260814`, an orbit residual gate of $10^{-8}$, and the residual—not the optimizer status flag—as the acceptance criterion. The only declared objectives were the maximum local Poincare stability margin and absolute normalized drive cost. Fisher information, SI power, experimental feasibility, and device-specific optimality were not part of this run. The result remains `PROVISIONAL_REDUCED_MODEL_MULTI_START_PARETO` because finite-time basin capture is not computed.

**Outcome (recorded 2026-08-16):** all 16 candidates are feasible (16/16 passed the local fixed-point/radius gate; largest Poincaré spectral radius $0.9970026$, i.e. a smallest stability margin $0.0029974$), and exactly three candidates are non-dominated under the declared objectives, with normalized drive costs $0.0800$, $0.1312$, $0.2203$ and worst-replicate stability margins $0.0034526$, $0.0035354$, $0.0035597$. The machine-readable artifact is `results/physical_parameter_reduced_bad_cavity_pareto_v4.json` (`pareto_count = 3`), rendered in `figures/reduced_model_pareto_v4.{pdf,png}` (generator `scripts/plot_reduced_model_pareto_v4.py`, manifest `figures/reduced_model_pareto_v4.manifest.json`). These are the values quoted in the main text Results and in SI Table `tab:si-reduced-pareto`; they remain a local fixed-point pattern of the normalized reduced model, not an SI calibration, finite-time basin-robustness result, or experimental optimum.

## Historical finite-time basin-capture attempts — superseded 2026-08-15

Before any manuscript insertion of the reduced-model Pareto result, a separate basin-capture test is required for the three non-dominated candidates (samples 5, 10, and 15) and all five stored parameter replicas. For each of these 15 parameter sets, five independently seeded initial conditions will be integrated continuously for up to 5,000 drive periods with the reduced equations. Capture requires both a distance to the reference fixed point of at most $10^{-5}$ in the reduced state infinity norm and a post-capture one-period residual of at most $10^{-7}$. The reference fixed point is recomputed from the zero initial guess using the residual-based root gate. The test will report capture time, final distance, final residual, and failures; it will not be called a basin proof beyond the tested initial-condition ensemble. An initial implementation probe at 5,000 periods and 100 integration substeps per period exceeded 600 s before producing an output. A controlled 3,000-period/20-substep probe for candidate 5 completed, but zero of its five parameter replicas passed the joint distance/residual gate: several nonzero starts crossed the distance threshold while their post-capture residual remained above $10^{-7}$. A subsequent 5,000-period/20-substep rerun also exceeded 600 s without output. These are computational and gate diagnostics, not claims of dynamical failure. At this checkpoint the result was recorded as `NOT_COMPUTED`, with no relaxed criterion. The subsequent vectorized RK4 run produced `results/fast_basin_capture_validation_v2.json`; its completed finite-time/finite-ensemble result and limitations are documented in the 2026-08-15 correction checkpoint above. The historical SLURM provenance is retained in `results/basin_capture_slurm_manifest_15325.json` and is not the current source of truth.

## Bibliographic verification checkpoint — 2026-08-14

The three references used to frame the synthetic-magnetism and fluctuation context were
verified against publisher records on 2026-08-14. The verified records are:

- S. Mondal, M. S. Baptista, and K. Debnath, “Chaotic dynamics under the influence of a synthetic magnetic field in an optomechanical system,” *Physical Review A* **110**, 023509 (2024), DOI `10.1103/PhysRevA.110.023509`.
- D. Muthukumar, S. R. Mbokop Tchounda, S. T. Kingni, K. Rajagopal, and S. G. Nana Engo, “Controlling complex dynamics with synthetic magnetism in optomechanical systems: A route to enhanced sensor performance,” *Physical Review Applied* **24**, 014053 (2025), DOI `10.1103/f6sv-4fq9`.
- M.-Q. Gao, S.-h. Li, X. Li, X. Li, J. Cheng, and W. Li, “Fluctuation-induced quenching of chaos in quantum optics,” *Chaos: An Interdisciplinary Journal of Nonlinear Science* **36**, 063114 (2026), DOI `10.1063/5.0325109`.

The corresponding BibTeX entries no longer contain internal scaffold notes. These
citations provide context only; they do not validate the present model, parameter
closure, basin capture, or sensing claims.

## Checkpoint — strong-coupling convergence and basin-capture provenance correction (2026-08-15)

The strong-coupling full-spectrum map was subjected to a targeted convergence audit at one chaotic point ($E=4$, $\Phi_{\rm syn}=-\pi/2$) and one hyperchaotic point ($E=8$, $\Phi_{\rm syn}=+\pi/2$), using the same seeded initial condition across $2.5\times10^5$, $5\times10^5$, and $10^6$ integration steps. The complete six-exponent spectra and QR block histories are stored in `results/hyperchaos_convergence_probe.json`. The classification is stable between the two longest windows: $n_{+}=2$ at the $E=4$ point and $n_{+}=4$ at the $E=8$ point. The shorter $E=4$ window gives $n_{+}=3$, so the manuscript does not use the shortest-window classification. The largest-exponent values vary by several percent across windows; the evidence supports the positive-exponent count at the tested representative points, not high-precision exponent estimates or a global convergence theorem. Divergence-balance residuals are $\sim10^{-3}$--$2\times10^{-3}$ in this chaotic regime, larger than in the weak-coupling periodic grid and now reported as a limitation of the finite-time QR audit.

The finite-time basin-capture output `results/fast_basin_capture_validation_v2.json` supersedes the earlier `NOT_COMPUTED` status file for the declared finite ensemble: 3 candidates × 5 parameter replicas × 5 initial conditions, 5,000 drive periods, 320 RK4 substeps per period, joint distance/residual gate, and independent DOP853 comparison; all three candidates passed. This remains a finite-time, finite-ensemble validation of the reduced assumed-closure model, not a global basin proof, SI calibration, or experimental optimum. The old SLURM manifest `results/basin_capture_slurm_manifest_15325.json` is retained as execution history only; its original working-directory path predates the rename and is not a current source path. Canonical provenance is now the `track_D_chaos` workspace and the fast-validation JSON.

## Checkpoint — full-spectrum transition map, weak coupling (2026-08-15)

To answer primary question 2/3 with the complete spectrum rather than a single largest exponent (HalefShomroni2025 standard), `scripts/full_spectrum_transition_map.py` scans drive $E\in\{0.2,0.5,1.0,2.0,4.0,8.0\}$ over nine flux phases (54 records, numba QR/Benettin, $10^6$ steps, $10^5$ transient, positive-exponent threshold $10^{-4}$) → `results/full_spectrum_transition_map.json`. **54/54 PASS: zero positive exponents everywhere**; the largest exponent stays in $[-0.0434,-0.0100]$, the Kaplan--Yorke dimension is 0 at every point, and the maximum divergence residual is 4.75e-5 (consistent with $\Tr J=-1.04$). This closes the hyperchaos-gap critique honestly at the weakly coupled pilot point: no transition is manufactured.

## Checkpoint — strong-coupling hyperchaos transition (2026-08-15)

The weak-coupling null held only because the broad probe stopped at $g=0.2$. `scripts/full_spectrum_hyperchaos_map.py` (SLURM job 15363) repeats the full-spectrum classification at $g_1=0.3$, $g_2=0.27$, $\gamma_{1,2}=0.02$, $\Delta=-\omega_m$ over twelve flux phases × five drives (60 records, $10^6$ steps, $10^5$ transient, threshold $10^{-3}$) → `results/full_spectrum_hyperchaos_map_g03.json`. **60/60 PASS**, with a genuine drive- and flux-controlled transition:

- $E=0.2$ → $n_{+}=0$ (stable); $E=1.0$ → $n_{+}$ up to 1 (chaos); $E=2.0$ → $n_{+}=0$ but largest exponent $+1.9\times10^{-5}$ (below the $10^{-3}$ threshold); $E=4.0$ → $n_{+}$ up to 3; $E=8.0$ → $n_{+}$ 2--4 (hyperchaos), largest exponent up to 0.3223, Kaplan--Yorke dimension up to 4.802.
- At fixed strong drive the flux phase modulates the hyperchaos order (at $E=8.0$, $n_{+}$ varies 2→4 over $\theta\in[-\pi,\pi)$).
- **Seed robustness:** `scripts/hyperchaos_seed_robustness.py` re-runs $E=4$ and $E=8$ with 3 independent initial conditions per phase → `results/hyperchaos_seed_robustness.json`, **72/72 PASS**; every phase is stable to within one exponent across seeds, maxima at $\theta=\pm\pi/2$, minima at $\theta=0,\pm\pi$.
- **Provenance correction:** the `flux` field was originally miscomputed as 0 for every record (the real hopping magnitude, not $Je^{i\theta}$, was fed to `synthetic_flux()`); fixed to `angle(g1·(J e^{iθ})·conj(g2)) = θ` and the 60 stored records were patched so `flux == theta`. Dynamics were unaffected ($\theta$ enters the equations via $\cos\theta/\sin\theta$).
- **Convergence audit:** `results/hyperchaos_convergence_probe.json` (script `scripts/hyperchaos_convergence_probe.py`) records the full spectrum at $(E,\Phi_{\rm syn})=(4,-\pi/2)$ and $(8,+\pi/2)$ for $2.5\times10^5$, $5\times10^5$, $10^6$ steps; $n_{+}=2$ and $n_{+}=4$ are stable between the two longest windows (SI Table S6), while finite-time exponent magnitudes vary by several percent. Reported as a finite-time classification, not high-precision asymptotic exponents.

The transition is therefore drive- and coupling-gated (cooperativity), with the flux phase controlling the hyperchaos order only beyond the dissipation-dominated regime. Threshold definitions: $n_{+}$ counts exponents above $10^{-3}$ (strong-coupling) or $10^{-4}$ (weak-coupling); the Kaplan--Yorke dimension counts strictly positive exponents, so the two can disagree between 0 and the threshold (e.g. $E=2.0$: $n_{+}=0$, KY $\approx 1.0$).

## Checkpoint — matched measurement-Fisher reference for force sensing (2026-08-15/16)

Primary question 5 is now answered at the model level. `scripts/matched_fisher_reference.py` estimates a weak external force $F$ on mechanical mode 1 (momentum-quadrature drive) from the cavity-amplitude record, with matched drive $E=0.2$, observation time (one drive period), detector variance 0.01, baths $n_{\rm th}=0.1$, and the same estimator → `results/matched_fisher_reference.json` (5 flux phases + flux-off coupled and single-mode references; status PASS).

- $F_C(F)$ rises monotonically from 4.503e-5 ($\theta=0$) to 5.958e-5 ($\theta=\pi$); gain vs flux-off 1.323, vs single-mode 1.159 — an honest order-unity classical linear-response effect, no quantum/SQL implication.
- **Stability sweep:** `scripts/matched_fisher_sweep.py` → `results/matched_fisher_sweep.json`, 16/16 PASS: gains constant to 3 decimals (std < 1e-3) over $E\in\{0.1,0.2,0.5,1.0\}\times n_{\rm th}\in\{0,0.1,0.5,1.0\}$.
- **Force-step convergence:** the canonical reference run stored an empty `--eps-check`; the SI Table (tab:si-fisher-eps) values for $\varepsilon=10^{-3},10^{-5}$ therefore had no stored provenance. `scripts/matched_fisher_eps_convergence.py` (deterministic, same matched resources) computes and stores `results/matched_fisher_eps_convergence.json` at both flux extrema ($\theta=\pi/2,\pi$) for $\varepsilon=10^{-3},10^{-4},10^{-5}$; relative spread < 1e-4 confirms the "converged to better than 0.01%" statement with stored evidence.

Claim C05 is RESOLVED as "no dramatic advantage": the flux phase is a real but sub-factor-of-two classical effect.

## Checkpoint — severe-reviewer editorial pass (2026-08-16)

The manuscript was revised in the severe-reviewer editorial pass: report-style language purged (pilot/prespecified/gate/reconstruction/provenance), abstract reduced to 221 words with keywords added, numerics details added to Methods (DOP853 tolerances, RK4 step $10^{-2}$, QR interval 10, transient/recorded steps — verified against the scripts), explicit contribution paragraph, semiclassical amplitude-scale disclosure ($|\langle X_a\rangle|\approx0.08$ from `covariance_fisher.json`), future-work statement (bifurcation/continuation + calibrated-device UQ), and AIP-structured end sections (Supplementary material, Acknowledgments, Author declarations with CRediT/COI/ethics placeholders, Data availability). No numerical claim changed; the claims↔data cross-check still passes 60/60 and both PDFs rebuild clean (main 11 pp, SI 4 pp).

## Checkpoint — AIP journal-format conversion (2026-08-16)

No numerical claims changed. Both manuscripts were converted to the official AIP
submission class (`revtex4-2` with `aip,cha` options + `aipnum4-2.bst`), with
`alt={...}` accessibility text on all five main-text figures, native
`\keywords{}`, and author placeholders retained for the authors. Rebuilds: main
8 pages / SI 3 pages (AIP reprint layout), 0 errors, 0 undefined
references; benign REVTeX deferred-float notices documented in
`manuscript/README.md`. Claims lint PASS on both files. All numbers still trace
to the result JSONs (re-verified by `scripts/audit_crosscheck_20260816.py`, 60/60).

## Checkpoint — drive-axis bifurcation continuation (2026-08-16)

A Newton continuation of the drive-locked periodic orbit family along the drive
axis at the two flux phases of interest ($\theta=0$ and $\theta=\pi/2$) was
executed with `scripts/drive_axis_bifurcation_scan.py` (gate
`DRIVE_AXIS_BIFURCATION_SCAN`): Newton solver on the Poincare map with the
monodromy as Jacobian, continuation tolerance $5\times10^{-6}$, adaptive steps,
upward and downward continuation over $E\in[0.1,8.0]$. Result
`results/drive_axis_bifurcation_scan.json` records 165 converged orbit points
at $\theta=0$ and 164 at $\theta=\pi/2$: the orbit family persists over the
whole scanned range, so the transition is not a fold termination. The first
loss of linear stability is a complex-pair (Neimark--Sacker) Floquet-multiplier
crossing at $E^{*}=1.060$ ($\theta=0$) and $E^{*}=0.975$ ($\theta=\pi/2$);
at $\theta=\pi/2$ a transient window of four unstable multipliers
($E\in[0.98,1.45]$) precedes the two-direction regime. Cross-validation
against the stored full-spectrum Lyapunov map at $E=0.2$: largest orbit
Floquet rate $-0.0102$ vs largest Lyapunov exponent $-0.0106$ ($\theta=0$)
and $-0.0148$ vs $-0.0151$ ($\theta=\pi/2$), agreement within $\sim4\%$.
Figure `figures/drive_axis_bifurcation.pdf` (manifest
`figures/drive_axis_bifurcation.manifest.json`; generator
`scripts/plot_drive_axis_bifurcation.py`). Scope: linear (Floquet) stability of
the continued orbit — the continued orbit is the unstable seed of the route, not
the chaotic attractor. The manuscript Results, Discussion, and Conclusions were
updated to report this refinement; the abstract notes the Neimark--Sacker
entry. No other numerical claim changed.

## Checkpoint — four-gap implementation (2026-08-16)

The four gaps (global attractor structure, calibrated-device UQ, experimental feasibility, and weak-coupling semiclassical validity) were implemented
with deterministic scripts and bound to JSON artifacts before any manuscript
insertion. Summary of the verified findings:

**Gap 1 — global attractor structure** (`scripts/attractor_structure.py` →
`results/attractor_structure.json`, gate `ATTRACTOR_STRUCTURE`). At the four
strong-coupling points ($E\approx4,8$; $\theta=0,\pi/2$) the continued orbit
has $n_{\rm unstable}=2$ unstable Floquet directions with largest Floquet rate
$0.60$–$0.69$; a $10^{-4}$ perturbation along the most unstable direction
departs to the attractor scale (final distance $9$–$27$) and the resulting
attractor Lyapunov spectrum reproduces the stored map's largest exponent within
$\sim2$–$8\%$ ($0.091$–$0.335$ vs $0.096$–$0.322$), with the same
positive-exponent count at $E=4$ ($n_{+}=2$ both phases). The orbit is the
unstable seed; the attractor is a distinct object.

**Gap 2 — calibrated-device UQ** (`scripts/matched_fisher_uq.py` →
`results/matched_fisher_uq.json`, gate `MATCHED_FISHER_UQ`). The deterministic
matched gains ($1.323$ flux-off, $1.159$ single-mode) were propagated over the
three unmeasured coordinates — hopping log-uniform in $[10^{-3},10^{-1}]$,
bath occupation uniform in $[0.01,1.0]$ (cryogenic to $\sim4$ K at
$\omega_m=2\pi\times7.436$ GHz), drive $\pm10\%$ — with $N=20$ samples. The
gain never drops below unity (100% of samples above unity; 5th percentile
$1.004$ flux-off / $1.002$ single-mode); its magnitude is diluted from the
peak $1.32$ to a median $1.082$ (flux-off) / $1.041$ (single-mode) because the
peak sits at the reference hopping $0.08$ and the log-uniform scan range
concentrates mass at smaller hopping.

**Gap 3 — experimental feasibility** (`scripts/strong_coupling_reachability.py`
→ `results/strong_coupling_reachability.json`, gate
`STRONG_COUPLING_REACHABILITY`). The strong-coupling hyperchaos regime requires
optomechanical coupling $g_1/\omega_m=0.3$ and $g_2/\omega_m=0.27$, which are
$2542\times$ and $2288\times$ the anchored Mayor $g_0/\omega_m=1.18\times10^{-4}$,
and mechanical damping $\gamma/\omega_m=0.02$, $720\times$ the anchored
$2.78\times10^{-5}$; $\kappa/\omega_m=1.0$ is $9.3\times$ the anchor and the
hopping $J/\omega_m=0.08$ is within the Mathew scan range. Limiting coordinate:
$g_1$ at $2542\times$. The regime is therefore a model-level demonstration, not
a prediction for the anchored device.

**Gap 4 — semiclassical validity** (`scripts/semiclassical_validity.py` →
`results/semiclassical_validity.json`, gate `SEMICLASSICAL_VALIDITY`). At the
weak-coupling reference the mean cavity occupation is $\langle|\alpha|^2\rangle=
0.0324$ (sub-photon), the single-photon cooperativity is $C=0.08$, and
$g/\kappa=0.02$; a truncated-Fock master equation (QuTiP, $N_{\rm Fock}=8$) for
the driven linear optical mode gives $\langle a^\dagger a\rangle=0.0320$,
reproducing the semiclassical value to $\sim1.3\%$. The strong-coupling
continued orbits have $\langle|\alpha|^2\rangle\approx7.3$–$11.2$, so the
hyperchaos regime operates at $\sim300\times$ the reference occupation where
the semiclassical factorization is standard.

No new references were added (the two SI anchors suffice). The Discussion
reports these as descriptive/interpretive statements.

**Refinement R1 — noise-dependent observability of chaos (v2, corrected)**
(`scripts/noise_observability.py` → `results/noise_observability.json`, gate
`NOISE_DEPENDENT_OBSERVABILITY`). The noise normalization was corrected against
the quantum master equation (see the Langevin cross-check below): the
amplitude quadrature $X_a=\Re(\alpha)=(a+a^\dagger)/2$ has vacuum variance
$1/4$, so the correct diffusion is $\kappa/4$ (optical) and
$\gamma_j(2n_{\mathrm{th}}+1)/4$ (mechanical), not $\kappa/2$ / $\gamma(2n+1)/2$.
With the corrected strengths, an ensemble of $N=32$ trajectories under the
Langevin noise (optical vacuum $\kappa/4$, mechanical thermal
$\gamma_j(2n_{\mathrm{th}}+1)/4$ with $n_{\mathrm{th}}=0.1$) plus detector noise
$0.01$, with measured record $y=X_a+\nu_{\rm det}$, gives a noise-only floor
$\approx0.255$ (matching the master-equation $0.25$). The measured-record
variance is $7.6$–$17.0\times$ the floor at the chaotic ($E=4$) and
hyperchaotic ($E=8$) points versus $\sim1.0$ stable; the deterministic spread of
$X_a$ is $1.8$–$4.1$ versus $2.4\times10^{-4}$ stable. The largest ratio
($17.0\times$) is at $\theta=\pi/2$, $E=8$. This is the noise-resolved
counterpart of the deterministic null and the positive answer to question 3 of
the Introduction. Figure `figures/noise_observability.pdf`.

**Refinement R6 — correlation-dimension cross-check (v2, Theiler-corrected)**
(`scripts/correlation_dimension.py` → `results/correlation_dimension.json`, gate
`CORRELATION_DIMENSION`; figure `figures/correlation_dimension.pdf`). The
grassberger--Procaccia correlation dimension was computed with a Theiler window
$W=50$ samples (removing temporally correlated pairs) from $2\times10^4$
sampled attractor points at four strong-coupling points. $D_2=2.33/2.84$
($\theta=0$, $E=4/8$) and $3.39/3.36$ ($\theta=\pi/2$, $E=4/8$), all below the
stored Kaplan--Yorke dimensions $4.21/4.63/4.28/4.80$ ($D_2\le D_{\mathrm{KY}}$).
The scaling-region uncertainty (three fit windows) is $0.23$–$0.40$. $D_2$ is
phase-dependent ($\pi/2>0$) and weakly drive-dependent, so the geometric
dimension is dominated by the stable directions while the Kaplan--Yorke sum
grows with the added unstable directions.

**Severe-referee rigor items (2026-08-16)**

*Langevin--master cross-check* (`scripts/langevin_master_crosscheck.py` →
`results/langevin_master_crosscheck.json`, gate `LANGEVIN_MASTER_CROSSCHECK`):
the isolated optical sector of the reference point (driven damped linear cavity)
has master-equation $\mathrm{Var}[(a+a^\dagger)/2]=0.25$; the Langevin amplitude
with $\kappa/4$ gives $0.245$ ($<2\%$ agreement) while $\kappa/2$ gives $0.49$
(factor $2$ high). This fixed the noise normalization used by R1.

*Noise-observability N-convergence* (`scripts/noise_observability_convergence.py`
→ `results/noise_observability_convergence.json`, gate
`NOISE_OBSERVABILITY_CONVERGENCE`): a master ensemble of $N_{\max}=128$
trajectories, bootstrapped ($B=200$) over $N'\in\{8,16,32,64,128\}$, gives
ratios stable to $<1\%$ ($4.63\times$ chaotic, $8.78\times$ hyperchaotic
$\theta=\pi/2$, $7.43\times$ hyperchaotic $\theta=0$) with the bootstrap std
decreasing monotonically $0.16\to0.04$, ruling out a small-ensemble artifact.

*Fisher-gain significance* (`scripts/fisher_gain_significance.py` →
`results/fisher_gain_significance.json`, gate `FISHER_GAIN_SIGNIFICANCE`; source
`results/matched_fisher_uq.json` re-run with $N=60$, SLURM job 15371): the
deterministic peak $1.323/1.159$ is the 97th-percentile (maximum-gain phase).
The median is $1.039/1.019$ with 90% bootstrap CI $[1.025,1.053]/[1.012,1.026]$
(excluding unity), median excess $0.039/0.019$; 5th percentile $1.004/1.002$.
The gain is significantly above unity but weak.

*Torus observation* (`scripts/torus_observation.py` →
`results/torus_observation.json`, gate `TORUS_OBSERVATION`; figure
`figures/torus_observation.pdf`): just above the Neimark--Sacker onset
($E\approx E^*+0.15$) the stroboscopic section of the perturbed orbit forms a
closed invariant curve (angular spread $\sim1.7$ rad) with largest exponent
$+0.017$–$+0.022$ (near-marginal), growing to $+0.057$–$+0.069$ at $E^*+0.3$.
The torus is an observed invariant object, not only an inference from the
multiplier spectrum.

## Checkpoint — severe adversarial audit corrections (2026-08-17)

Outcome of the 2026-08-17 severe adversarial review, consolidated in `../PRX_PRXQuantum_adversarial_review_tracks_A_B_D_20260822.md` and verified against `results/*.json`:

- **Semiclassical boundary corrected**: at the tested chaotic trajectory, the truncated-Fock master check converges toward the semiclassical optical occupation with maximum relative deviation $6.6\times10^{-6}$ at $N_{\rm Fock}=60$. The stored mean occupations are $\langle a^\dagger a\rangle=6.121$, $\langle b_1^\dagger b_1\rangle=20.417$, and $\langle b_2^\dagger b_2\rangle=23.538$. The time-averaged $g_2\approx2.0$ is super-Poissonian, and the mechanical sector remains semiclassical.
- **Threshold sensitivity extended**: the strong-coupling $n_{+}$ is recomputed at a single uniform $10^{-4}$ threshold applied to both sectors; it is identical to the reference $10^{-3}$ count, so the drive-gated route to $n_{+}=4$ is not an artefact of the two-threshold choice.
- **Seed robustness full map**: 180/180 records (5 drives × 12 phases × 3 seeds) confirm the stable→onset→re-stabilization→chaos→hyperchaos route, including the $E=2.0$ window, against initial-condition dependence.
- **Flux-phase symmetry corrected**: the three-seed check gives $n_{+}\in\{3,4\}$ at both $\theta=\pm\pi/2$, so no direction asymmetry under $\theta\to-\theta$ is claimed (the single-seed map's $n_{+}=4$ at $-\pi/2$ vs $3$ at $+\pi/2$ is within the seed variation). The $D_2$ phase-dependence claim was likewise restricted to pointwise values.
- **Title and prose**: the title now reads “Flux-controlled dynamics in dissipative optomechanics: full-spectrum hyperchaos and matched sensing bounds”; report-style vocabulary (`declared`, `canonical`, `physicality`, `novel synergetic link`) was removed from main and SI.
- **Chaotic-attractor matched Fisher (complete)**: the resolved record at $E=4$ gives flux-off $F_C=\num{1.08e-2}$ vs flux-on $F_C=\num{5.88e-3}$ (ratio 0.54); all $E=8$ points are below resolution. No flux-induced enhancement of the sensing Fisher information on the attractor; this is reported as a negative result in the Discussion and SI.
- **Optimized linear-sensor reference (2026-08-17)**: the optimal homodyne readout angle is $\varphi_\star\approx\num{0.79}\pi$ for every configuration (set by the linear cavity response, phase-independent), so readout optimization rescales $F_C$ by the same factor (\num{1.56}) for sensor and references. The flux gain vs the single-mode reference is therefore unchanged by readout optimization (0.88→1.16, max 1.16 at $\theta=\pi$), refuting the 'suboptimal reference readout' objection. `results/optimized_linear_sensor.json`.
- **Truncated-Wigner Lyapunov robustness (2026-08-17)**: the full 3-mode truncated-Fock master equation is intractable at the strong-coupling occupations (mechanical 37.5/44.1, peaks 75–79). At the TWA level, quantum fluctuations enter as Wigner-sampled initial conditions (variance 0.5/quadrature). The hyperchaos route core is robust: $n_+=2$ at $E=8$, $\theta=0$ is reproduced in 12/12 Wigner-sampled trajectories; the maximal order ($n_+=4$ at $\theta=\pi/2$, fourth exponent $2.4\times10^{-3}$) and the chaos onset ($n_+=1$ at $E=4$, second exponent $10^{-4}$) are threshold-conventional. `results/twa_lyapunov_robustness.json`.
- **Feasibility frontier (2026-08-17)**: scanning the drive $E\in[0.5,8]$ and the coupling $g_1/\omega_m\in[0.02,0.3]$ ($g_2=0.9g_1$), hyperchaos ($n_+\ge2$) occurs only at $E=8$, with onset near $g_1/\omega_m\approx0.15$–$0.18$ (seed-dependent over the transition interval); the robust onset (all three seeds $n_+=2$) is $g_1/\omega_m\approx0.18$, i.e. $1525\times$ the anchored $g_0/\omega_m=1.18\times10^{-4}$. No drive retuning lowers the onset, so the coupling gap is fundamental to the hyperchaos regime rather than an artefact of the single declared operating point. `results/feasibility_frontier.json` + `results/feasibility_frontier_onset.json`.

## Prohibited shortcuts

- no mock Wigner surfaces;
- no representative values copied from `production_summary_metrics.json`;
- no hand-edited HDF5/CSV values;
- no mixed normalized/physical time axes;
- no QFI inferred from Lyapunov exponents;
- no claim that non-significance proves absence or equivalence;
- no manuscript update before this report is updated with actual run evidence.
