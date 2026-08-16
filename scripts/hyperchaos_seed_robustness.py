#!/usr/bin/env python3
"""Seed-robustness check for the strong-coupling hyperchaos transition map.

The canonical map (full_spectrum_hyperchaos_map_g03.json) uses one distinct seed
per (drive, phase) point, so the flux-phase dependence of the hyperchaos order
(n_pos) at strong drive rests on single realizations. This script re-runs the
full six-exponent spectrum at the hyperchaotic drive values (E = 4.0, 8.0) with
n_seeds independent initial conditions per phase, and reports per-phase n_pos
statistics (min/mean/max) to show the phase modulation is robust to initial
conditions.
"""
from __future__ import annotations
import argparse, json, platform, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

from full_spectrum_hyperchaos_map import evaluate


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drives", type=str, default="4.0,8.0")
    ap.add_argument("--points", type=int, default=12)
    ap.add_argument("--n-seeds", type=int, default=3)
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
    ap.add_argument("--seed-base", type=int, default=20260815)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()

    p_kwargs = dict(g1=a.g1, g2=a.g2, gamma1=a.gamma1, gamma2=a.gamma2,
                    detuning=a.detuning)
    drives = [float(x) for x in a.drives.split(",")]
    thetas = np.linspace(-np.pi, np.pi, a.points, endpoint=False)

    tasks = []
    for di, drive in enumerate(drives):
        for ti, theta in enumerate(thetas):
            for si in range(a.n_seeds):
                seed = a.seed_base + di * a.points * a.n_seeds + ti * a.n_seeds + si
                tasks.append((di, drive, ti, float(theta), seed, p_kwargs, a.n_steps,
                              a.transient_steps, a.dt, a.qr_interval,
                              a.positive_threshold))

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

    # per-phase statistics
    summary = {}
    for di, drive in enumerate(drives):
        for ti, theta in enumerate(thetas):
            rs = [r for r in records
                  if r.get("drive_index") == di and r.get("theta_index") == ti]
            npos = [r["n_positive_exponents"] for r in rs if r.get("status") == "PASS"]
            lam = [r["largest_exponent"] for r in rs if r.get("status") == "PASS"]
            summary[f"{drive:.1f}"] = summary.get(f"{drive:.1f}", [])
            summary[f"{drive:.1f}"].append({
                "theta": float(theta),
                "n_pos_min": min(npos) if npos else None,
                "n_pos_mean": float(np.mean(npos)) if npos else None,
                "n_pos_max": max(npos) if npos else None,
                "n_pos_values": npos,
                "largest_exponent_min": min(lam) if lam else None,
                "largest_exponent_max": max(lam) if lam else None,
                "seeds_passed": len(npos),
            })

    out = {
        "gate": "HYPERCHAOS_SEED_ROBUSTNESS",
        "kind": "full_spectrum_hyperchaos_seed_robustness",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "coupling": {"g1": a.g1, "g2": a.g2, "gamma1": a.gamma1, "gamma2": a.gamma2,
                     "detuning": a.detuning},
        "drives": drives, "n_phases": a.points, "n_seeds_per_point": a.n_seeds,
        "positive_exponent_threshold": a.positive_threshold,
        "total_records": len(records), "passed_records": len(passed),
        "summary_by_drive": summary,
        "records": records,
        "interpretation": ("Per-phase hyperchaos-order (n_pos) statistics over "
                           "independent initial conditions at strong coupling. The "
                           "flux-phase modulation of n_pos is robust if the min/max "
                           "spread at each phase is small and the phase-to-phase "
                           "variation exceeds the seed spread."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"wrote {a.output}: {len(passed)}/{len(records)} PASS")
    for drive in drives:
        print(f"  E={drive:.1f} n_pos per phase (min~max over seeds):")
        for row in summary[f"{drive:.1f}"]:
            print(f"    theta={row['theta']/np.pi:+.2f}pi: "
                  f"n_pos {row['n_pos_values']} (min={row['n_pos_min']}, mean={row['n_pos_mean']:.2f}, max={row['n_pos_max']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
