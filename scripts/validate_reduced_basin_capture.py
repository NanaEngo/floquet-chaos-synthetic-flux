#!/usr/bin/env python3
"""Finite-time basin-capture test for selected reduced-model Pareto candidates."""
from __future__ import annotations

import argparse
import json
import platform
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from reconstruction_core import ModelParameters
from reduced_bad_cavity_screen import (
    find_periodic_orbit_root,
    integrate_period,
)


def integrate_periods(y0: np.ndarray, p: ModelParameters, n_periods: int,
                      max_step_fraction: int) -> np.ndarray:
    """Integrate one continuous trajectory and return one state per period."""
    if n_periods <= 0 or max_step_fraction <= 0:
        raise ValueError("n_periods and max_step_fraction must be positive")

    def fun(t: float, y: np.ndarray) -> np.ndarray:
        from reduced_bad_cavity_screen import rhs_reduced
        return rhs_reduced(y, p.drive_frequency * t, p)

    final_time = n_periods * p.period
    t_eval = np.arange(1, n_periods + 1, dtype=float) * p.period
    sol = solve_ivp(
        fun,
        (0.0, final_time),
        np.asarray(y0, dtype=float),
        method="DOP853",
        t_eval=t_eval,
        rtol=1e-10,
        atol=1e-12,
        max_step=p.period / float(max_step_fraction),
    )
    if not sol.success or sol.y.shape[1] != n_periods:
        raise RuntimeError(sol.message)
    return np.asarray(sol.y.T, dtype=float)


def parameter_from_dict(values: dict) -> ModelParameters:
    allowed = {
        "kappa", "gamma1", "gamma2", "omega1", "omega2", "g1", "g2",
        "hopping", "detuning", "drive", "drive_modulation", "drive_frequency",
        "theta",
    }
    return ModelParameters(**{key: float(value) for key, value in values.items() if key in allowed})


def test_parameter_set(p: ModelParameters, rng: np.random.Generator, *,
                      initial_conditions: int, max_periods: int,
                      capture_tol: float, residual_tol: float,
                      max_step_fraction: int) -> dict:
    reference = find_periodic_orbit_root(np.zeros(4), p, max_iter=80, residual_tol=1e-8)
    if reference["status"] != "PASS":
        return {"status": "FAIL", "reason": "reference root did not pass", "reference": reference}
    reference_state = np.asarray(reference["y0"], dtype=float)
    records = []
    for initial_condition in range(initial_conditions):
        y0 = np.zeros(4) if initial_condition == 0 else rng.normal(0.0, 0.1, 4)
        trajectory = integrate_periods(y0, p, max_periods, max_step_fraction)
        distances = np.max(np.abs(trajectory - reference_state[None, :]), axis=1)
        passed_record = None
        first_distance_period = None
        minimum_post_residual = None
        for capture_index in np.flatnonzero(distances <= capture_tol):
            capture_index = int(capture_index)
            if first_distance_period is None:
                first_distance_period = capture_index + 1
            capture_state = trajectory[capture_index]
            post = integrate_period(capture_state, p)
            post_residual = float(np.linalg.norm(post.y[:, -1] - capture_state, ord=np.inf))
            minimum_post_residual = post_residual if minimum_post_residual is None else min(minimum_post_residual, post_residual)
            if post_residual <= residual_tol:
                passed_record = {
                    "initial_condition": initial_condition,
                    "status": "PASS",
                    "capture_period": capture_index + 1,
                    "first_distance_threshold_period": first_distance_period,
                    "final_distance_at_capture": float(distances[capture_index]),
                    "post_capture_residual": post_residual,
                }
                break
        if passed_record is not None:
            records.append(passed_record)
        else:
            records.append({
                "initial_condition": initial_condition,
                "status": "FAIL",
                "reason": "joint distance/residual capture threshold not reached",
                "first_distance_threshold_period": first_distance_period,
                "final_distance": float(distances[-1]),
                "minimum_post_capture_residual": minimum_post_residual,
            })
    return {
        "status": "PASS" if all(record["status"] == "PASS" for record in records) else "FAIL",
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
    parser.add_argument("--capture-tol", type=float, default=1e-5)
    parser.add_argument("--residual-tol", type=float, default=1e-7)
    parser.add_argument("--max-step-fraction", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260815)
    args = parser.parse_args()
    source = json.loads(args.screen.read_text(encoding="utf-8"))
    wanted = set(args.samples)
    candidates = {int(item["sample"]): item for item in source["candidate_records"]}
    missing = sorted(wanted - candidates.keys())
    if missing:
        raise ValueError(f"Pareto candidate records missing from screen: {missing}")

    rng = np.random.default_rng(args.seed)
    candidate_results = []
    for sample in sorted(wanted):
        candidate = candidates[sample]
        replicate_results = []
        for rep in range(args.initial_condition_replicates):
            values = candidate["replicate_parameters"].get(str(rep))
            if values is None:
                raise ValueError(f"missing exact parameters for sample={sample}, replicate={rep}")
            p = parameter_from_dict(values)
            replicate_results.append({
                "replicate": rep,
                "parameters": p.to_dict(),
                "capture": test_parameter_set(
                    p, rng,
                    initial_conditions=args.initial_condition_replicates,
                    max_periods=args.max_periods,
                    capture_tol=args.capture_tol,
                    residual_tol=args.residual_tol,
                    max_step_fraction=args.max_step_fraction,
                ),
            })
        candidate_results.append({
            "sample": sample,
            "replicate_results": replicate_results,
            "candidate_pass": all(r["capture"]["status"] == "PASS" for r in replicate_results),
        })

    out = {
        "status": "PASS" if all(item["candidate_pass"] for item in candidate_results) else "FAIL",
        "scientific_status": "PROVISIONAL_REDUCED_MODEL_BASIN_CAPTURE_VALIDATION",
        "kind": "finite_time_basin_capture_validation",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "source_screen": str(args.screen.resolve()),
        "design": {
            "samples": sorted(wanted),
            "candidate_count": len(candidate_results),
            "parameter_replicates_per_candidate": args.initial_condition_replicates,
            "initial_conditions_per_replicate": args.initial_condition_replicates,
            "max_periods": args.max_periods,
            "capture_tolerance": args.capture_tol,
            "post_capture_residual_tolerance": args.residual_tol,
            "max_step_fraction": args.max_step_fraction,
            "seed": args.seed,
            "reference_orbit": "root from zero initial guess, residual gate 1e-8",
        },
        "candidate_results": candidate_results,
        "candidate_pass_count": sum(item["candidate_pass"] for item in candidate_results),
        "interpretation": "Finite-time capture validation for the tested 3-candidate, 5-replica, 5-initial-condition ensemble only; it does not prove global basin structure, SI calibration, or experimental optimality.",
        "limitations": [
            "adiabatic optical elimination",
            "assumed J_m and g_i proxies",
            "finite ensemble rather than global basin proof",
            "Fisher information not computed",
            "drive remains a normalized model coordinate",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: out[key] for key in ("status", "scientific_status", "candidate_pass_count")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
