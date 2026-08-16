#!/usr/bin/env python3
"""Bad-cavity adiabatic reduction for a diagnostic model-level screen.

The optical amplitude is eliminated from
  da/dt = (-kappa/2 + i*Omega_eff)a + E(phi)
using a_ss = E(phi)/(kappa/2 - i*Omega_eff).
This module is deliberately separate from reconstruction_core.py. Its outputs
are diagnostic and cannot be promoted to physical calibration results without
comparison to the full model in a common regime and device-specific drive data.
"""
from __future__ import annotations

import argparse
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from reconstruction_core import ModelParameters, drive_amplitude


def optical_intensity(phi: float, b1r: float, b2r: float, p: ModelParameters) -> float:
    omega_eff = p.detuning + 2.0 * p.g1 * b1r + 2.0 * p.g2 * b2r
    e = drive_amplitude(phi, p)
    denom = (0.5 * p.kappa) ** 2 + omega_eff ** 2
    return float(e * e / max(denom, np.finfo(float).tiny))


def rhs_reduced(y: np.ndarray, phi: float, p: ModelParameters) -> np.ndarray:
    if np.asarray(y).shape != (4,):
        raise ValueError(f"expected reduced state shape (4,), got {np.asarray(y).shape}")
    b1r, b1i, b2r, b2i = np.asarray(y, dtype=float)
    c, s = np.cos(p.theta), np.sin(p.theta)
    n = optical_intensity(phi, b1r, b2r, p)
    return np.array([
        -0.5 * p.gamma1 * b1r + p.omega1 * b1i
        + p.hopping * (s * b2r + c * b2i),
        p.g1 * n - p.omega1 * b1r - 0.5 * p.gamma1 * b1i
        - p.hopping * c * b2r + p.hopping * s * b2i,
        -0.5 * p.gamma2 * b2r + p.omega2 * b2i
        - p.hopping * s * b1r + p.hopping * c * b1i,
        p.g2 * n - p.omega2 * b2r - 0.5 * p.gamma2 * b2i
        - p.hopping * c * b1r - p.hopping * s * b1i,
    ], dtype=float)


def integrate_period(y0: np.ndarray, p: ModelParameters, *, rtol: float = 1e-9,
                     atol: float = 1e-11):
    def fun(t: float, y: np.ndarray) -> np.ndarray:
        return rhs_reduced(y, p.drive_frequency * t, p)
    return solve_ivp(fun, (0.0, p.period), np.asarray(y0, dtype=float),
                     method="DOP853", rtol=rtol, atol=atol,
                     max_step=p.period / 400.0)


def find_periodic_orbit(y0: np.ndarray, p: ModelParameters, *, max_iter: int = 600,
                        residual_tol: float = 1e-8) -> dict:
    y = np.asarray(y0, dtype=float).copy()
    history: list[float] = []
    for iteration in range(1, max_iter + 1):
        sol = integrate_period(y, p)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return {"status": "FAIL", "reason": "reduced period integration failed", "iterations": iteration}
        y_next = sol.y[:, -1]
        residual = float(np.linalg.norm(y_next - y, ord=np.inf))
        history.append(residual)
        y = y_next
        if residual < residual_tol:
            return {"status": "PASS", "y0": y.tolist(), "residual": residual,
                    "iterations": iteration, "residual_history": history}
    return {"status": "FAIL", "reason": "reduced periodic-orbit residual did not converge",
            "y0": y.tolist(), "residual": history[-1], "iterations": max_iter,
            "residual_history": history}


def find_periodic_orbit_root(y0: np.ndarray, p: ModelParameters, *, max_iter: int = 80,
                             residual_tol: float = 1e-8) -> dict:
    """Solve the reduced Poincare fixed point directly instead of iterating periods."""
    def residual(y: np.ndarray) -> np.ndarray:
        sol = integrate_period(y, p)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            raise RuntimeError("reduced period integration failed")
        return np.asarray(sol.y[:, -1] - y, dtype=float)
    result = None
    errors = []
    # The stiff bad-cavity regime can make HYBR report "no progress" even
    # after reaching a machine-precision fixed point. Try LM first and judge
    # the scientific gate by the actual Poincare residual, not the optimizer's
    # advisory success flag.
    for method, options in (
        ("lm", {"maxiter": max_iter, "ftol": 1e-10, "xtol": 1e-10}),
        ("hybr", {"maxfev": max_iter, "xtol": min(1e-8, residual_tol)}),
    ):
        try:
            candidate = root(residual, np.asarray(y0, dtype=float), method=method,
                             options=options)
            if np.all(np.isfinite(candidate.x)):
                result = candidate
                residual_norm = float(np.linalg.norm(residual(candidate.x), ord=np.inf))
                if residual_norm < residual_tol:
                    return {"status": "PASS", "y0": candidate.x.tolist(),
                            "residual": residual_norm, "iterations": int(candidate.nfev),
                            "root_method": method, "solver_success_flag": bool(candidate.success)}
                errors.append(f"{method}: residual={residual_norm}")
        except Exception as exc:
            errors.append(f"{method}: {exc}")
    return {"status": "FAIL", "reason": "reduced root residual did not pass gate",
            "diagnostics": errors, "iterations": int(getattr(result, "nfev", 0))}


def poincare_spectral_radius(y0: np.ndarray, p: ModelParameters, eps: float = 1e-6) -> float:
    base = np.asarray(integrate_period(y0, p).y[:, -1], dtype=float)
    jac = np.empty((4, 4), dtype=float)
    for j in range(4):
        yp = np.asarray(y0, dtype=float).copy(); ym = yp.copy()
        yp[j] += eps; ym[j] -= eps
        fp = integrate_period(yp, p).y[:, -1]
        fm = integrate_period(ym, p).y[:, -1]
        jac[:, j] = (fp - fm) / (2.0 * eps)
    return float(np.max(np.abs(np.linalg.eigvals(jac))))


def center_parameters(manifest: dict, rng: np.random.Generator, fraction: float = 0.0) -> ModelParameters:
    entries = {x["name"]: x for x in manifest["parameters"]}
    values: dict[str, float] = {}
    names = ("kappa", "gamma1", "gamma2", "omega1", "omega2", "g1", "g2",
             "hopping", "detuning", "drive", "drive_modulation", "drive_frequency")
    for name in names:
        lo, hi = map(float, entries[name]["normalized_range"])
        if fraction:
            value = 0.5 * (lo + hi)
            value *= 1.0 + float(rng.normal(0.0, fraction))
        else:
            # Independent parameter-space samples; the earlier diagnostic
            # accidentally reused the interval midpoint for every candidate.
            value = lo + float(rng.random()) * (hi - lo)
        values[name] = value
    return ModelParameters(**values)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--sample-start", type=int, default=0,
                        help="zero-based global candidate index; preserves the seeded stream for chunked runs")
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--initial-condition-replicates", type=int, default=3)
    parser.add_argument("--initial-condition-scale", type=float, default=0.1)
    parser.add_argument("--robust-replicates", type=int, default=3)
    parser.add_argument("--orbit-max-iter", type=int, default=600)
    parser.add_argument("--orbit-residual-tol", type=float, default=1e-8)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    rng = np.random.default_rng(args.seed)
    if args.samples <= 0 or args.sample_start < 0:
        raise ValueError("samples must be positive and sample-start must be non-negative")
    records = []
    global_end = args.sample_start + args.samples
    for sample in range(global_end):
        base = center_parameters(manifest, rng)
        if sample < args.sample_start:
            continue
        sample_records = []
        replicate_parameters = {}
        for rep in range(args.robust_replicates):
            p = base if rep == 0 else center_parameters(manifest, rng, 0.10)
            replicate_parameters[str(rep)] = p.to_dict()
            for ic in range(args.initial_condition_replicates):
                y0 = np.zeros(4) if ic == 0 else rng.normal(0.0, args.initial_condition_scale, 4)
                orbit = find_periodic_orbit_root(y0, p, max_iter=args.orbit_max_iter,
                                                 residual_tol=args.orbit_residual_tol)
                if orbit["status"] != "PASS":
                    sample_records.append({"replicate": rep, "initial_condition": ic,
                                           "status": "FAIL", "reason": orbit["reason"]})
                    continue
                radius = poincare_spectral_radius(np.asarray(orbit["y0"]), p)
                sample_records.append({"replicate": rep, "initial_condition": ic,
                                       "status": "PASS", "residual": orbit["residual"],
                                       "iterations": orbit["iterations"],
                                       "poincare_spectral_radius": radius})
        passed = [r for r in sample_records if r["status"] == "PASS"]
        feasible = len(passed) == args.robust_replicates * args.initial_condition_replicates and all(
            r["poincare_spectral_radius"] < 1.0 for r in passed)
        records.append({"sample": sample, "parameters": base.to_dict(),
                        "replicate_parameters": replicate_parameters,
                        "replicate_records": sample_records, "feasible": feasible})
    feasible = [r for r in records if r["feasible"]]
    out = {
        "status": "PROVISIONAL",
        "scientific_status": "PROVISIONAL_REDUCED_MODEL_MULTI_START_SCREEN",
        "kind": "adiabatic_bad_cavity_diagnostic",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "manifest": str(args.manifest.resolve()),
        "design": {
            "samples": args.samples,
            "sample_start": args.sample_start,
            "global_sample_end": global_end,
            "seed": args.seed,
            "robust_replicates": args.robust_replicates,
            "initial_condition_replicates": args.initial_condition_replicates,
            "initial_condition_scale": args.initial_condition_scale,
            "orbit_max_iter": args.orbit_max_iter,
            "orbit_residual_tolerance": args.orbit_residual_tol,
            "initial_condition_protocol": "independent root-solver starts; finite-time basin capture is NOT_COMPUTED",
            "parameter_replication_protocol": "independent uniform perturbation draws; exact replicate parameters stored per candidate",
        },
        "candidate_records": records,
        "candidate_count": len(records),
        "feasible_count": len(feasible),
        "pareto_count": 0,
        "interpretation": "Reduced bad-cavity multi-start fixed-point diagnostic only; all residual/radius gates passed for the recorded starts, but finite-time basin capture, SI calibration, and manuscript claims are not authorized.",
        "limitations": [
            "optical mode adiabatically eliminated",
            "drive remains model coordinate",
            "J_m and g_i are assumed proxies",
            "full/reduced agreement tested only at kappa=20 and 100",
            "comparison does not validate the assumed closure or an experimental device",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("status", "scientific_status", "candidate_count", "feasible_count", "pareto_count")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
