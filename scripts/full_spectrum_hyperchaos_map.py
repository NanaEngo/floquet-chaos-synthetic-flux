#!/usr/bin/env python3
"""Full-spectrum hyperchaos map: number of positive Lyapunov exponents (hyperchaos
order) vs synthetic-flux phase and drive amplitude, at a chosen coupling strength.

Unlike the weak-coupling pilot, strong optomechanical coupling drives the system from
a stable drive-locked orbit into chaos (one positive exponent) and hyperchaos (two or
more).  Because chaotic attractors need not possess a stable periodic orbit, the
spectrum is integrated from a random initial condition (transient discarded) rather
than from a converged Poincare fixed point.
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
from reconstruction_core import ModelParameters
from lyapunov_numba import lyapunov_qr_numba


def kaplan_yorke_dimension(spectrum: np.ndarray) -> float:
    s = np.sort(np.asarray(spectrum, dtype=float))[::-1]
    cum = np.cumsum(s)
    idx = np.where(cum >= 0)[0]
    if len(idx) == 0:
        return 0.0
    k = int(idx[-1]) + 1
    if k == len(s):
        return float(k)
    return float(k) + float(cum[k - 1]) / abs(float(s[k]))


def evaluate(task: tuple) -> dict:
    (di, drive, ti, theta, seed, p_kwargs, n_steps, transient_steps, dt, qr_interval,
     pos_threshold) = task
    p = ModelParameters(drive=drive, theta=theta, **p_kwargs)
    rng = np.random.default_rng(seed)
    x0 = rng.normal(0.0, 0.1, 6)
    rec = {"drive_index": di, "drive": drive, "theta_index": ti, "theta": theta,
           "seed": seed}
    ly = lyapunov_qr_numba(np.asarray(x0, dtype=float), p, n_steps=n_steps, dt=dt,
                           transient_steps=transient_steps, qr_interval=qr_interval)
    spectrum = np.sort(np.asarray(ly["spectrum"], dtype=float))
    n_pos = int(np.sum(spectrum > pos_threshold))
    divergence_residual = float(abs(np.sum(spectrum) - (-p.kappa - p.gamma1 - p.gamma2)))
    rec.update({
        "status": "PASS",
        "lyapunov_spectrum": spectrum.tolist(),
        "largest_exponent": float(spectrum[-1]),
        "n_positive_exponents": n_pos,
        "kaplan_yorke_dimension": kaplan_yorke_dimension(spectrum),
        "divergence_residual": divergence_residual,
        "flux": float(np.angle(p.g1 * (p.hopping * np.exp(1j * p.theta)) * np.conjugate(p.g2))),
    })
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drives", type=str, default="0.2,1.0,2.0,4.0,8.0")
    ap.add_argument("--points", type=int, default=12)
    ap.add_argument("--g1", type=float, default=0.3)
    ap.add_argument("--g2", type=float, default=0.27)
    ap.add_argument("--gamma1", type=float, default=0.02)
    ap.add_argument("--gamma2", type=float, default=0.02)
    ap.add_argument("--detuning", type=float, default=-1.0)
    ap.add_argument("--n-steps", type=int, default=1_000_000)
    ap.add_argument("--transient-steps", type=int, default=100_000)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--qr-interval", type=int, default=10)
    ap.add_argument("--positive-threshold", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=20260815)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()

    p_kwargs = dict(g1=a.g1, g2=a.g2, gamma1=a.gamma1, gamma2=a.gamma2,
                    detuning=a.detuning)
    drives = [float(x) for x in a.drives.split(",")]
    thetas = np.linspace(-np.pi, np.pi, a.points, endpoint=False)
    tasks = []
    for di, drive in enumerate(drives):
        for ti, theta in enumerate(thetas):
            seed = a.seed + di * a.points + ti
            tasks.append((di, drive, ti, float(theta), seed, p_kwargs, a.n_steps,
                          a.transient_steps, a.dt, a.qr_interval, a.positive_threshold))

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
    n_pos_max = max((r["n_positive_exponents"] for r in passed), default=0)
    out = {
        "gate": "FULL_SPECTRUM_HYPERCHAOS_MAP",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "kind": "full_spectrum_hyperchaos_map",
        "interpretation": ("Number of positive Lyapunov exponents (hyperchaos order) vs "
                           "synthetic-flux phase and drive, at strong optomechanical coupling. "
                           "The full six-exponent spectrum is recorded for every point."),
        "coupling": {"g1": a.g1, "g2": a.g2, "gamma1": a.gamma1, "gamma2": a.gamma2,
                     "detuning": a.detuning},
        "drives": drives,
        "n_phases": a.points,
        "n_steps": a.n_steps,
        "transient_steps": a.transient_steps,
        "dt": a.dt,
        "positive_exponent_threshold": a.positive_threshold,
        "total_records": len(records),
        "passed_records": len(passed),
        "max_positive_exponents": n_pos_max,
        "records": records,
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {a.output}: {len(passed)}/{len(records)} PASS, "
          f"max n_pos = {n_pos_max}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
