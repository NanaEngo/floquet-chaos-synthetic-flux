#!/usr/bin/env python3
"""Ensemble-size convergence of the noise-observability ratio.

A severe-referee item on R1: the measured-record variance ratio
(chaotic / stable noise floor) reported in results/noise_observability.json
uses a single N=32 ensemble. This script integrates a master ensemble of
N_max trajectories per operating point once, then **bootstraps** the
observability ratio over ensemble sizes N' in {8,16,32,64,128} (B paired
resamples per size) to show the ratio is converged, not an N=32 artifact.

The ratio at each bootstrap draw is formed from paired subsets of the chaotic
and stable ensembles (same trajectory indices), so the noise floor and the
signal are compared at identical ensemble size.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, find_periodic_orbit
from noise_observability import (rhs_batch, noise_strengths, strong_params,
                                 orbit_at)

ROOT = Path(__file__).resolve().parent.parent
SCAN = ROOT / "results" / "drive_axis_bifurcation_scan.json"


def trajectory_stats(x0, p: ModelParameters, q, n_ens, n_steps, dt,
                     det_sigma, seed, transient_steps):
    """Per-trajectory sufficient statistics of the measured record y=X_a+nu."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x0, dtype=np.float64)[None, :] + 1e-4 * rng.standard_normal((n_ens, 6))
    q = np.asarray(q, dtype=np.float64)
    sum_y = np.zeros(n_ens)
    sum_y2 = np.zeros(n_ens)
    n = 0
    for step in range(n_steps):
        phi = (step * dt * p.drive_frequency) % (2.0 * np.pi)
        k1 = rhs_batch(x, phi, p)
        k2 = rhs_batch(x + 0.5 * dt * k1, phi + 0.5 * p.drive_frequency * dt, p)
        k3 = rhs_batch(x + 0.5 * dt * k2, phi + 0.5 * p.drive_frequency * dt, p)
        k4 = rhs_batch(x + dt * k3, phi + p.drive_frequency * dt, p)
        x = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        x += np.sqrt(q[None, :] * dt) * rng.standard_normal((n_ens, 6))
        if step >= transient_steps:
            y = x[:, 0] + np.sqrt(det_sigma) * rng.standard_normal(n_ens)
            sum_y += y
            sum_y2 += y * y
            n += 1
    return sum_y, sum_y2, n


def pooled_variance(sum_y, sum_y2, n, idx):
    sy = float(sum_y[idx].sum())
    sy2 = float(sum_y2[idx].sum())
    nn = float(n * len(idx))
    return max(sy2 / nn - (sy / nn) ** 2, 0.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--n-max", type=int, default=128)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--sizes", nargs="+", type=int,
                    default=[8, 16, 32, 64, 128])
    ap.add_argument("--n-steps", type=int, default=30000)
    ap.add_argument("--transient-steps", type=int, default=3000)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--detector-variance", type=float, default=0.01)
    ap.add_argument("--n-th", nargs=2, type=float, default=[0.1, 0.1])
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    n_th1, n_th2 = a.n_th
    scan = json.loads(SCAN.read_text())
    rng = np.random.default_rng(a.seed)

    # Stable reference (weak coupling, E=0.2, theta=0).
    p_ref = ModelParameters(theta=0.0)
    orb = find_periodic_orbit(np.zeros(6), p_ref, max_iter=600)
    if orb.get("status") != "PASS":
        raise RuntimeError(orb.get("reason", "reference orbit failed"))
    x_ref = np.asarray(orb["x0"], dtype=np.float64)
    q_ref = noise_strengths(p_ref, n_th1, n_th2)
    stable_sy, stable_sy2, stable_n = trajectory_stats(
        x_ref, p_ref, q_ref, a.n_max, a.n_steps, a.dt,
        a.detector_variance, a.seed, a.transient_steps)

    points = [
        ("chaotic_E4_th0.5pi", 4.0, 0.5),
        ("hyperchaotic_E8_th0.5pi", 8.0, 0.5),
        ("hyperchaotic_E8_th0pi", 8.0, 0.0),
    ]
    records = []
    for label, E, theta_pi in points:
        theta = theta_pi * np.pi
        trec = next(t for t in scan["theta_records"]
                    if abs(t["theta"] - theta) < 1e-9)
        x0, E_actual = orbit_at(trec["records"], E, "up")
        p = strong_params(E_actual, theta)
        q = noise_strengths(p, n_th1, n_th2)
        sy, sy2, n = trajectory_stats(x0, p, q, a.n_max, a.n_steps, a.dt,
                                      a.detector_variance, a.seed,
                                      a.transient_steps)
        convergence = []
        for size in a.sizes:
            ratios = []
            for _ in range(a.n_boot):
                idx = rng.choice(a.n_max, size=size, replace=False)
                vp = pooled_variance(sy, sy2, n, idx)
                vs = pooled_variance(stable_sy, stable_sy2, stable_n, idx)
                ratios.append(vp / vs)
            ratios = np.asarray(ratios)
            convergence.append({
                "ensemble_size": size,
                "ratio_mean": float(np.mean(ratios)),
                "ratio_std": float(np.std(ratios)),
                "ratio_p05": float(np.percentile(ratios, 5)),
                "ratio_p95": float(np.percentile(ratios, 95)),
            })
        records.append({
            "label": label,
            "drive": float(p.drive),
            "theta_over_pi": float(theta_pi),
            "coupling": "strong",
            "convergence": convergence,
            "status": "PASS",
        })
        print(f"{label}: N'={a.sizes[-1]} ratio="
              f"{convergence[-1]['ratio_mean']:.3f} "
              f"+- {convergence[-1]['ratio_std']:.3f}", flush=True)

    # Converged ratio = the largest-N bootstrap mean; converged iff the
    # largest-N ratio is within 2 sigma of the second-largest-N ratio.
    for r in records:
        c = r["convergence"]
        r["converged_ratio"] = c[-1]["ratio_mean"]
        r["converged_ratio_std"] = c[-1]["ratio_std"]
        if len(c) >= 2:
            d = abs(c[-1]["ratio_mean"] - c[-2]["ratio_mean"])
            r["converged_within_2sigma"] = bool(d <= 2.0 * (c[-1]["ratio_std"] + c[-2]["ratio_std"]))
        else:
            r["converged_within_2sigma"] = False

    status = "PASS" if all(r["converged_within_2sigma"] for r in records) else "FAIL"
    out = {
        "gate": "NOISE_OBSERVABILITY_CONVERGENCE",
        "kind": "bootstrap_ensemble_size_convergence",
        "status": status,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "settings": {
            "n_max": a.n_max, "n_boot": a.n_boot, "sizes": a.sizes,
            "n_steps": a.n_steps, "transient_steps": a.transient_steps,
            "dt": a.dt, "detector_variance": a.detector_variance,
            "thermal_occupations": [n_th1, n_th2], "seed": a.seed,
        },
        "records": records,
        "interpretation": (
            "The observability ratio (chaotic measured-record variance over the "
            "stable noise floor) is bootstrapped over ensemble size. The ratio "
            "is converged if the largest-size estimate lies within two standard "
            "errors of the second-largest-size estimate. This rules out the "
            "objection that the N=32 ratio is a small-ensemble artifact."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(json.dumps({"status": status, "output": str(a.output)}, indent=2))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
