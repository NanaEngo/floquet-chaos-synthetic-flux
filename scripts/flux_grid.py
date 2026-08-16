#!/usr/bin/env python3
"""Explore flux-dependent dynamics without silently promoting exploratory points."""
from __future__ import annotations

import argparse
import json
import platform
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, find_periodic_orbit, lyapunov_qr, monodromy
try:
    from lyapunov_numba import lyapunov_qr_numba
except ImportError:  # documented fallback for non-production technical tests
    lyapunov_qr_numba = None


def evaluate_record(task: tuple[int, float, int, int, int, int, float, int, float, float, float, str]) -> dict:
    """Evaluate one independent phase/initial-condition replicate."""
    theta_index, theta, replicate, seed, n_steps, transient_steps, dt, qr_interval, max_diff, max_div, max_block, backend = task
    p = ModelParameters(theta=float(theta))
    rng = np.random.default_rng(seed)
    x0 = np.zeros(6) if replicate == 0 else rng.normal(0.0, 0.02, 6)
    orbit = find_periodic_orbit(x0, p, max_iter=600)
    record = {
        "theta_index": theta_index,
        "theta": float(theta),
        "flux": float(theta),
        "replicate": replicate,
        "seed": seed,
        "orbit_status": orbit.get("status"),
        "orbit_residual": orbit.get("residual"),
    }
    if orbit.get("status") != "PASS":
        record.update({"status": "FAIL", "reason": orbit.get("reason", "orbit failed")})
        return record
    try:
        floquet = monodromy(np.asarray(orbit["x0"]), p)
        if backend == "numba":
            if lyapunov_qr_numba is None:
                raise RuntimeError("Numba backend requested but unavailable")
            lyap = lyapunov_qr_numba(
                np.asarray(orbit["x0"]), p, n_steps=n_steps, dt=dt,
                transient_steps=transient_steps, qr_interval=qr_interval,
            )
        else:
            lyap = lyapunov_qr(
                np.asarray(orbit["x0"]), p, n_steps=n_steps, dt=dt,
                transient_steps=transient_steps, qr_interval=qr_interval,
            )
        floquet_rates = np.sort(np.asarray(floquet["floquet_rates"], dtype=float))
        spectrum = np.sort(np.asarray(lyap["spectrum"], dtype=float))
        rate_difference = float(np.max(np.abs(floquet_rates - spectrum)))
        divergence_residual = float(abs(np.sum(spectrum) - lyap["mean_divergence"]))
        blocks = np.asarray(lyap["block_history"], dtype=float)
        tail = blocks[-10:] if len(blocks) >= 10 else blocks
        block_std = float(np.max(np.std(tail, axis=0))) if len(tail) > 1 else float("inf")
        record.update(
            {
                "status": "PASS"
                if rate_difference < max_diff
                and divergence_residual < max_div
                and block_std < max_block
                else "FAIL",
                "floquet_rates": floquet_rates.tolist(),
                "lyapunov_spectrum": spectrum.tolist(),
                "max_rate_difference": rate_difference,
                "divergence_residual": divergence_residual,
                "tail_block_std": block_std,
                "is_chaotic": bool(spectrum[-1] > 1e-6),
                "backend": backend,
                "gate_thresholds": {
                    "max_rate_difference": max_diff,
                    "max_divergence_residual": max_div,
                    "max_block_std": max_block,
                },
            }
        )
    except Exception as exc:  # keep one failed record from hiding other records
        record.update({"status": "FAIL", "reason": f"{type(exc).__name__}: {exc}"})
    return record


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--points", type=int, default=9)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--n-steps", type=int, default=2_000_000)
    parser.add_argument("--transient-steps", type=int, default=200_000)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--qr-interval", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260813)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--backend", choices=("numba", "python"), default="numba")
    parser.add_argument("--max-rate-difference", type=float, default=5e-4)
    parser.add_argument("--max-divergence-residual", type=float, default=5e-3)
    parser.add_argument("--max-block-std", type=float, default=5e-3)
    args = parser.parse_args()
    if args.points < 3 or args.replicates < 2:
        raise SystemExit("points>=3 and replicates>=2 are required")
    if args.workers < 1:
        raise SystemExit("workers must be >= 1")

    thetas = np.linspace(-np.pi, np.pi, args.points, endpoint=False)
    tasks = []
    for theta_index, theta in enumerate(thetas):
        for replicate in range(args.replicates):
            record_seed = args.seed + theta_index * args.replicates + replicate
            tasks.append(
                (
                    theta_index,
                    float(theta),
                    replicate,
                    record_seed,
                    args.n_steps,
                    args.transient_steps,
                    args.dt,
                    args.qr_interval,
                    args.max_rate_difference,
                    args.max_divergence_residual,
                    args.max_block_std,
                    args.backend,
                )
            )

    records = []
    if args.workers == 1:
        records = [evaluate_record(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(evaluate_record, task): i for i, task in enumerate(tasks)}
            completed = {}
            for future in as_completed(futures):
                completed[futures[future]] = future.result()
            records = [completed[i] for i in range(len(tasks))]

    for record in records:
        print(json.dumps(record, sort_keys=True), flush=True)
    output = {
        "status": "PASS" if records and all(r.get("status") == "PASS" for r in records) else "FAIL",
        "kind": "exploratory_flux_grid",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "points": args.points,
        "replicates": args.replicates,
        "seed": args.seed,
        "workers": args.workers,
        "backend": args.backend,
        "n_steps": args.n_steps,
        "transient_steps": args.transient_steps,
        "dt": args.dt,
        "qr_interval": args.qr_interval,
        "gate_thresholds": {
            "max_rate_difference": args.max_rate_difference,
            "max_divergence_residual": args.max_divergence_residual,
            "max_block_std": args.max_block_std,
        },
        "records": records,
        "interpretation": "Grid diagnostics only; no chaos or sensing claim is accepted from this file without independent review.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"], "records": len(records), "output": str(args.output)}, indent=2))
    return 0 if output["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
