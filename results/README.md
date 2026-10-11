# Results directory

Only machine-generated outputs from the reconstruction may be stored here.
The current package is the full production set behind the Chaos, Solitons &
Fractals resubmission (89 JSON artifacts; `checksums.sha256` is the
authoritative inventory). The headline result is the drive- and
coupling-gated full-spectrum hyperchaos route with a null matched-sensing
outcome, not a periodic-orbit pilot.

## Gate and provenance core

- `environment_manifest.json` — fl_qom environment, FL_QOM revision, worktree state, and dependency versions;
- `model_gate.json` — equations/Jacobian/trace/gauge-flux gate;
- `floquet_lyapunov.json` — periodic orbit, monodromy, Floquet rates, QR spectrum, and convergence history;
- `parameters.json` — normalized pilot parameters;
- `status.json` — aggregate gate status;
- `run_manifest.json` — provenance links and execution metadata;
- `source_hashes.json` — hashes of the reconstruction scripts;
- `checksums.sha256` — hashes of the result files;
- `flux_grid.json` — canonical 27-run long-window production grid (9 phases × 3 independent replicates), all strict QR/Floquet gates passing;
- `covariance_fisher.json` — three-phase periodic covariance and classical measurement-Fisher pilot;
- `physical_parameter_optimization_provisional.json` — superseded provisional robustness screen (0/2 candidates passed the three-initial-condition gate); no physical optimum;
- `physical_parameter_optimization_provisional_v2.json` — definitive eight-candidate literature-envelope robustness screen (0/8 feasible, 0 non-dominated); no physical optimum;
- `physical_parameter_optimization_model_recovery.json` — positive normalized-pilot recovery screen (4/4 feasible, 2 non-dominated); `PROVISIONAL_MODEL_LEVEL`, not SI-calibrated;
- `physical_parameter_run_manifest_model_recovery.json` — provenance for the positive normalized-pilot recovery screen;
- `physical_parameter_pareto_drive_translation.json` — conditional normalized input-amplitude translation; symbolic power only, not SI calibration;
- `physical_parameter_stability_screening_provisional.json` — superseded parameter-only stage retained for audit;
- `../data/physical_parameter_manifest_literature_scenario.json` — source and assumptions for the provisional screen (`data/`, outside `results/`, not covered by `checksums.sha256`);
- `physical_parameter_run_manifest.json` — command, environment, design, and limits for the superseded screen;
- `physical_parameter_run_manifest_v2.json` — exact provenance and interpretation of the definitive eight-candidate robustness screen.

The full-spectrum maps, seed/window convergence audits, basin and
attractor-geometry analyses, semiclassical master-equation cross-checks,
matched-Fisher sensing results, and reduced-model comparisons are listed
individually in `checksums.sha256`.

All listed gates have `extended_status=PASS`. The production grid and
covariance/Fisher calculation pass their declared gates; the manuscript
reports the resulting hyperchaos classification and the null flux-induced
sensing gain, and still claims no QFI, Wigner negativity, or exceptional
point. The former 30,000-step screening output is retained as
`flux_grid_screening_30000.json` for audit history and is not used by the
manuscript.

The full QR block history is retained in `floquet_lyapunov.json` for auditability;
its large size is intentional. No hand-edited values, mock data, representative
summaries, screenshots, or figures belong here. Failed and incomplete runs may be
retained for audit but must be marked explicitly and must never feed a manuscript
generator. The provisional Pareto output is diagnostic evidence only; the final
initial-condition gate found no feasible candidate, so it does not upgrade claim C09
or authorize the words “optimal” or “experimentally feasible.”
