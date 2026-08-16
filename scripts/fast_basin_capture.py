#!/usr/bin/env python3
"""Vectorized finite-time basin-capture validation for the reduced equations."""
from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters
from reduced_bad_cavity_screen import find_periodic_orbit_root, integrate_period


def rhs_batch(states: np.ndarray, phi: float, p: ModelParameters) -> np.ndarray:
    b1r, b1i, b2r, b2i = states.T
    omega_eff = p.detuning + 2.0 * p.g1 * b1r + 2.0 * p.g2 * b2r
    e = p.drive * (1.0 + p.drive_modulation * np.cos(phi))
    n = e * e / np.maximum((0.5 * p.kappa) ** 2 + omega_eff * omega_eff, np.finfo(float).tiny)
    c, s = np.cos(p.theta), np.sin(p.theta)
    return np.column_stack((
        -0.5 * p.gamma1 * b1r + p.omega1 * b1i + p.hopping * (s * b2r + c * b2i),
        p.g1 * n - p.omega1 * b1r - 0.5 * p.gamma1 * b1i - p.hopping * c * b2r + p.hopping * s * b2i,
        -0.5 * p.gamma2 * b2r + p.omega2 * b2i - p.hopping * s * b1r + p.hopping * c * b1i,
        p.g2 * n - p.omega2 * b2r - 0.5 * p.gamma2 * b2i - p.hopping * c * b1r - p.hopping * s * b1i,
    ))


def rk4_period_batch(states: np.ndarray, p: ModelParameters, substeps: int) -> np.ndarray:
    h = p.period / float(substeps)
    y = np.asarray(states, dtype=float).copy()
    t = 0.0
    for _ in range(substeps):
        phi = p.drive_frequency * t
        k1 = rhs_batch(y, phi, p)
        k2 = rhs_batch(y + 0.5 * h * k1, phi + 0.5 * p.drive_frequency * h, p)
        k3 = rhs_batch(y + 0.5 * h * k2, phi + 0.5 * p.drive_frequency * h, p)
        k4 = rhs_batch(y + h * k3, phi + p.drive_frequency * h, p)
        y += h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        t += h
    return y


def parameter_from_dict(values: dict) -> ModelParameters:
    allowed = {"kappa", "gamma1", "gamma2", "omega1", "omega2", "g1", "g2",
               "hopping", "detuning", "drive", "drive_modulation", "drive_frequency", "theta"}
    return ModelParameters(**{key: float(value) for key, value in values.items() if key in allowed})


def validate_candidate(p: ModelParameters, rng: np.random.Generator, *, initial_conditions: int,
                       max_periods: int, substeps: int, capture_tol: float,
                       residual_tol: float, residual_check_interval: int) -> dict:
    reference = find_periodic_orbit_root(np.zeros(4), p, max_iter=80, residual_tol=1e-8)
    if reference["status"] != "PASS":
        return {"status": "FAIL", "reason": "reference root failed", "reference": reference}
    reference_state = np.asarray(reference["y0"], dtype=float)
    states = np.zeros((initial_conditions, 4), dtype=float)
    if initial_conditions > 1:
        states[1:] = rng.normal(0.0, 0.1, size=(initial_conditions - 1, 4))
    first_distance_period = np.full(initial_conditions, -1, dtype=int)
    capture_period = np.full(initial_conditions, -1, dtype=int)
    post_residual = np.full(initial_conditions, np.nan)
    min_distance = np.full(initial_conditions, np.inf)
    for period in range(1, max_periods + 1):
        states = rk4_period_batch(states, p, substeps)
        distances = np.max(np.abs(states - reference_state[None, :]), axis=1)
        min_distance = np.minimum(min_distance, distances)
        for i in range(initial_conditions):
            if first_distance_period[i] < 0 and distances[i] <= capture_tol:
                first_distance_period[i] = period
            if capture_period[i] >= 0 or distances[i] > capture_tol:
                continue
            if period % residual_check_interval != 0 and period != first_distance_period[i]:
                continue
            check = integrate_period(states[i], p)
            residual = float(np.linalg.norm(check.y[:, -1] - states[i], ord=np.inf))
            post_residual[i] = residual
            if residual <= residual_tol:
                capture_period[i] = period
        if np.all(capture_period >= 0):
            break
    records = []
    for i in range(initial_conditions):
        passed = capture_period[i] >= 0
        records.append({
            "initial_condition": i,
            "status": "PASS" if passed else "FAIL",
            "first_distance_threshold_period": None if first_distance_period[i] < 0 else int(first_distance_period[i]),
            "capture_period": None if capture_period[i] < 0 else int(capture_period[i]),
            "minimum_distance": float(min_distance[i]),
            "post_capture_residual": None if not np.isfinite(post_residual[i]) else float(post_residual[i]),
        })
    return {
        "status": "PASS" if all(r["status"] == "PASS" for r in records) else "FAIL",
        "reference_residual": float(reference["residual"]),
        "initial_condition_records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, nargs="+", default=[5, 10, 15])
    parser.add_argument("--initial-condition-replicates", type=int, default=5)
    parser.add_argument("--max-periods", type=int, default=5000)
    parser.add_argument("--substeps", type=int, default=40)
    parser.add_argument("--capture-tol", type=float, default=1e-5)
    parser.add_argument("--residual-tol", type=float, default=1e-7)
    parser.add_argument("--residual-check-interval", type=int, default=25)
    parser.add_argument("--seed", type=int, default=20260815)
    args = parser.parse_args()
    source = json.loads(args.screen.read_text(encoding="utf-8"))
    candidates = {int(item["sample"]): item for item in source["candidate_records"]}
    missing = sorted(set(args.samples) - candidates.keys())
    if missing:
        raise ValueError(f"missing samples: {missing}")
    rng = np.random.default_rng(args.seed)
    candidate_results = []
    for sample in sorted(set(args.samples)):
        item = candidates[sample]
        reps = []
        for rep in range(args.initial_condition_replicates):
            values = item["replicate_parameters"].get(str(rep))
            if values is None:
                raise ValueError(f"missing replicate parameters for sample={sample}, rep={rep}")
            p = parameter_from_dict(values)
            reps.append({"replicate": rep, "parameters": p.to_dict(), "capture": validate_candidate(
                p, rng, initial_conditions=args.initial_condition_replicates,
                max_periods=args.max_periods, substeps=args.substeps,
                capture_tol=args.capture_tol, residual_tol=args.residual_tol,
                residual_check_interval=args.residual_check_interval)})
        candidate_results.append({"sample": sample, "replicate_results": reps,
                                 "candidate_pass": all(r["capture"]["status"] == "PASS" for r in reps)})
    out = {
        "status": "PASS" if all(c["candidate_pass"] for c in candidate_results) else "FAIL",
        "scientific_status": "PROVISIONAL_REDUCED_MODEL_BASIN_CAPTURE_VALIDATION",
        "kind": "vectorized_finite_time_basin_capture_validation",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "source_screen": str(args.screen.resolve()),
        "design": {
            "samples": sorted(set(args.samples)),
            "candidate_count": len(candidate_results),
            "parameter_replicates_per_candidate": args.initial_condition_replicates,
            "initial_conditions_per_replicate": args.initial_condition_replicates,
            "max_periods": args.max_periods,
            "rk4_substeps_per_period": args.substeps,
            "capture_tolerance": args.capture_tol,
            "post_capture_residual_tolerance": args.residual_tol,
            "residual_check_interval": args.residual_check_interval,
            "seed": args.seed,
            "integrator": "vectorized classical RK4; independently compared with DOP853",
        },
        "candidate_results": candidate_results,
        "candidate_pass_count": sum(c["candidate_pass"] for c in candidate_results),
        "interpretation": "Finite-time capture validation for the tested ensemble only; it does not prove global basin structure, SI calibration, or experimental optimality.",
        "limitations": ["adiabatic optical elimination", "assumed J_m and g_i proxies", "finite ensemble", "Fisher information not computed", "normalized drive"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: out[k] for k in ("status", "scientific_status", "candidate_pass_count")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
