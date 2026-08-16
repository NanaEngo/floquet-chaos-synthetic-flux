#!/usr/bin/env python3
"""Grouped DOP853 finite-time basin-capture validation."""
from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from reconstruction_core import ModelParameters
from reduced_bad_cavity_screen import find_periodic_orbit_root, integrate_period
from fast_basin_capture import rhs_batch


def parameter_from_dict(values: dict) -> ModelParameters:
    allowed = {"kappa", "gamma1", "gamma2", "omega1", "omega2", "g1", "g2",
               "hopping", "detuning", "drive", "drive_modulation", "drive_frequency", "theta"}
    return ModelParameters(**{key: float(value) for key, value in values.items() if key in allowed})


def integrate_grouped(initial_states: np.ndarray, p: ModelParameters, max_periods: int,
                      max_step_fraction: int) -> np.ndarray:
    n = initial_states.shape[0]
    def fun(t: float, flat: np.ndarray) -> np.ndarray:
        states = flat.reshape(n, 4)
        phi = p.drive_frequency * t
        return rhs_batch(states, phi, p).ravel()
    t_eval = np.arange(1, max_periods + 1, dtype=float) * p.period
    sol = solve_ivp(
        fun,
        (0.0, max_periods * p.period),
        np.asarray(initial_states, dtype=float).ravel(),
        method="DOP853",
        t_eval=t_eval,
        rtol=1e-10,
        atol=1e-12,
        max_step=p.period / float(max_step_fraction),
    )
    if not sol.success or sol.y.shape[1] != max_periods:
        raise RuntimeError(sol.message)
    return np.asarray(sol.y.T, dtype=float).reshape(max_periods, n, 4)


def validate_parameter_set(p: ModelParameters, rng: np.random.Generator, *, initial_conditions: int,
                           max_periods: int, max_step_fraction: int, capture_tol: float,
                           residual_tol: float) -> dict:
    reference = find_periodic_orbit_root(np.zeros(4), p, max_iter=80, residual_tol=1e-8)
    if reference["status"] != "PASS":
        return {"status": "FAIL", "reason": "reference root failed", "reference": reference}
    reference_state = np.asarray(reference["y0"], dtype=float)
    states0 = np.zeros((initial_conditions, 4), dtype=float)
    if initial_conditions > 1:
        states0[1:] = rng.normal(0.0, 0.1, size=(initial_conditions - 1, 4))
    trajectory = integrate_grouped(states0, p, max_periods, max_step_fraction)
    distances = np.max(np.abs(trajectory - reference_state[None, None, :]), axis=2)
    records = []
    for i in range(initial_conditions):
        candidate_indices = np.flatnonzero(distances[:, i] <= capture_tol)
        passed = None
        first_distance = None if len(candidate_indices) == 0 else int(candidate_indices[0] + 1)
        minimum_residual = None
        for idx in candidate_indices:
            state = trajectory[idx, i]
            post = integrate_period(state, p)
            residual = float(np.linalg.norm(post.y[:, -1] - state, ord=np.inf))
            minimum_residual = residual if minimum_residual is None else min(minimum_residual, residual)
            if residual <= residual_tol:
                passed = {"initial_condition": i, "status": "PASS", "capture_period": int(idx + 1),
                          "first_distance_threshold_period": first_distance,
                          "distance_at_capture": float(distances[idx, i]),
                          "post_capture_residual": residual}
                break
        if passed is None:
            records.append({"initial_condition": i, "status": "FAIL",
                            "reason": "joint distance/residual gate not reached",
                            "first_distance_threshold_period": first_distance,
                            "final_distance": float(distances[-1, i]),
                            "minimum_post_capture_residual": minimum_residual})
        else:
            records.append(passed)
    return {"status": "PASS" if all(r["status"] == "PASS" for r in records) else "FAIL",
            "reference_residual": float(reference["residual"]),
            "initial_condition_records": records}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, nargs="+", default=[5, 10, 15])
    parser.add_argument("--initial-condition-replicates", type=int, default=5)
    parser.add_argument("--max-periods", type=int, default=5000)
    parser.add_argument("--max-step-fraction", type=int, default=20)
    parser.add_argument("--capture-tol", type=float, default=1e-5)
    parser.add_argument("--residual-tol", type=float, default=1e-7)
    parser.add_argument("--seed", type=int, default=20260815)
    args = parser.parse_args()
    source = json.loads(args.screen.read_text(encoding="utf-8"))
    candidates = {int(item["sample"]): item for item in source["candidate_records"]}
    missing = sorted(set(args.samples) - candidates.keys())
    if missing:
        raise ValueError(f"missing samples: {missing}")
    rng = np.random.default_rng(args.seed)
    results = []
    for sample in sorted(set(args.samples)):
        item = candidates[sample]
        reps = []
        for rep in range(args.initial_condition_replicates):
            p = parameter_from_dict(item["replicate_parameters"][str(rep)])
            reps.append({"replicate": rep, "parameters": p.to_dict(),
                         "capture": validate_parameter_set(
                             p, rng, initial_conditions=args.initial_condition_replicates,
                             max_periods=args.max_periods, max_step_fraction=args.max_step_fraction,
                             capture_tol=args.capture_tol, residual_tol=args.residual_tol)})
        results.append({"sample": sample, "replicate_results": reps,
                       "candidate_pass": all(r["capture"]["status"] == "PASS" for r in reps)})
    out = {
        "status": "PASS" if all(r["candidate_pass"] for r in results) else "FAIL",
        "scientific_status": "PROVISIONAL_REDUCED_MODEL_BASIN_CAPTURE_VALIDATION",
        "kind": "grouped_dopri853_finite_time_basin_capture_validation",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "source_screen": str(args.screen.resolve()),
        "design": {"samples": sorted(set(args.samples)), "candidate_count": len(results),
                   "parameter_replicates_per_candidate": args.initial_condition_replicates,
                   "initial_conditions_per_replicate": args.initial_condition_replicates,
                   "max_periods": args.max_periods, "max_step_fraction": args.max_step_fraction,
                   "capture_tolerance": args.capture_tol,
                   "post_capture_residual_tolerance": args.residual_tol, "seed": args.seed,
                   "integrator": "grouped DOP853; five initial conditions integrated in one system per parameter set"},
        "candidate_results": results,
        "candidate_pass_count": sum(r["candidate_pass"] for r in results),
        "interpretation": "Finite-time capture validation for the tested ensemble only; it does not prove global basin structure, SI calibration, or experimental optimality.",
        "limitations": ["adiabatic optical elimination", "assumed J_m and g_i proxies", "finite ensemble", "Fisher information not computed", "normalized drive"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: out[k] for k in ("status", "scientific_status", "candidate_pass_count")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
