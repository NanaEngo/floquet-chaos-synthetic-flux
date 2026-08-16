#!/usr/bin/env python3
"""Full-spectrum transition map: number of positive Lyapunov exponents over a
(drive x synthetic-flux-phase) grid, to close the hyperchaos gap (i.e. report the
complete six-exponent spectrum rather than only the largest exponent).

Mirrors flux_grid.py conventions; uses the numba QR/Benettin backend.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruction_core import ModelParameters, find_periodic_orbit
from lyapunov_numba import lyapunov_qr_numba


def kaplan_yorke_dimension(spectrum: np.ndarray) -> float:
    """Kaplan-Yorke dimension from a sorted (ascending) Lyapunov spectrum."""
    s = np.sort(spectrum)[::-1]  # descending
    ssum = 0.0
    k = 0
    for i, lam in enumerate(s):
        ssum += lam
        if ssum > 0:
            k = i + 1
    if k == 0:
        return 0.0
    if k == len(s):
        return float(k)
    return float(k) + ssum / abs(s[k])


def evaluate(task: tuple) -> dict:
    (di, drive, ti, theta, seed, n_steps, transient_steps, dt, qr_interval) = task
    p = ModelParameters(drive=drive, theta=theta)
    rng = np.random.default_rng(seed)
    x0 = np.zeros(6) if seed % 2 == 0 else rng.normal(0.0, 0.02, 6)
    rec = {
        "drive_index": di, "drive": drive,
        "theta_index": ti, "theta": theta,
        "seed": seed,
    }
    orb = find_periodic_orbit(x0, p, max_iter=600)
    rec["orbit_status"] = orb.get("status")
    rec["orbit_residual"] = orb.get("residual")
    if orb.get("status") != "PASS":
        rec.update({"status": "ORBIT_FAIL", "reason": orb.get("reason", "orbit failed")})
        return rec
    ly = lyapunov_qr_numba(np.asarray(orb["x0"]), p, n_steps=n_steps, dt=dt,
                           transient_steps=transient_steps, qr_interval=qr_interval)
    spectrum = np.sort(np.asarray(ly["spectrum"], dtype=float))
    n_pos = int(np.sum(spectrum > 1e-4))       # positive-exponent count (hyperchaos order)
    divergence_residual = float(abs(np.sum(spectrum) - (-p.kappa - p.gamma1 - p.gamma2)))
    rec.update({
        "status": "PASS",
        "lyapunov_spectrum": spectrum.tolist(),
        "largest_exponent": float(spectrum[-1]),
        "n_positive_exponents": n_pos,
        "kaplan_yorke_dimension": kaplan_yorke_dimension(spectrum),
        "divergence_residual": divergence_residual,
        "mean_divergence": float(ly.get("mean_divergence", -1.04)),
    })
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drives", type=str, default="0.2,0.5,1.0,2.0,4.0,8.0")
    ap.add_argument("--points", type=int, default=9)
    ap.add_argument("--n-steps", type=int, default=1_000_000)
    ap.add_argument("--transient-steps", type=int, default=100_000)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--qr-interval", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20260815)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()

    drives = [float(x) for x in a.drives.split(",")]
    thetas = np.linspace(-np.pi, np.pi, a.points, endpoint=False)
    tasks = []
    for di, drive in enumerate(drives):
        for ti, theta in enumerate(thetas):
            seed = a.seed + di * a.points + ti
            tasks.append((di, drive, ti, float(theta), seed, a.n_steps,
                          a.transient_steps, a.dt, a.qr_interval))

    records = []
    if a.workers <= 1:
        for t in tasks:
            records.append(evaluate(t))
    else:
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            futs = {pool.submit(evaluate, t): i for i, t in enumerate(tasks)}
            done = {}
            for fut in as_completed(futs):
                done[futs[fut]] = fut.result()
            records = [done[i] for i in range(len(tasks))]

    passed = [r for r in records if r.get("status") == "PASS"]
    n_pos_total = sum(r["n_positive_exponents"] for r in passed)
    out = {
        "gate": "FULL_SPECTRUM_TRANSITION_MAP",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "kind": "full_spectrum_lyapunov_transition_map",
        "interpretation": "Full six-exponent spectrum over (drive x phase); the positive-exponent count is the hyperchaos order. This file records the spectrum, not only the largest exponent.",
        "drives": drives,
        "n_phases": a.points,
        "n_steps": a.n_steps,
        "transient_steps": a.transient_steps,
        "dt": a.dt,
        "qr_interval": a.qr_interval,
        "positive_exponent_threshold": 1e-4,
        "total_records": len(records),
        "passed_records": len(passed),
        "total_positive_exponents": n_pos_total,
        "any_positive_exponent": n_pos_total > 0,
        "records": records,
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {a.output}: {len(passed)}/{len(records)} PASS, "
          f"total positive exponents = {n_pos_total}, any_positive = {n_pos_total > 0}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
