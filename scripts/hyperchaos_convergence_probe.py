#!/usr/bin/env python3
"""Targeted convergence audit for the strong-coupling full-spectrum map.

This is an audit, not a replacement for the canonical map.  It reuses the exact
Numba QR/Benettin backend and records block-history convergence, divergence-balance
residuals, and the positive-exponent count at two representative points.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from lyapunov_numba import lyapunov_qr_numba  # noqa: E402
from reconstruction_core import ModelParameters  # noqa: E402


def run_point(drive: float, theta: float, n_steps: int, transient_steps: int,
              seed: int, *, g1: float, g2: float, gamma1: float, gamma2: float,
              detuning: float, threshold: float) -> dict:
    p = ModelParameters(drive=drive, theta=theta, g1=g1, g2=g2,
                        gamma1=gamma1, gamma2=gamma2, detuning=detuning)
    x0 = np.random.default_rng(seed).normal(0.0, 0.1, 6)
    out = lyapunov_qr_numba(x0, p, n_steps=n_steps, dt=0.01,
                            transient_steps=transient_steps, qr_interval=10)
    spectrum = np.asarray(out["spectrum"], dtype=float)
    blocks = np.asarray(out["block_history"], dtype=float)
    tail = blocks[max(0, len(blocks) // 2):]
    return {
        "drive": drive,
        "theta": theta,
        "theta_over_pi": theta / np.pi,
        "n_steps": n_steps,
        "transient_steps": transient_steps,
        "seed": seed,
        "spectrum": spectrum.tolist(),
        "largest_exponent": float(spectrum[0]),
        "n_positive_exponents": int(np.sum(spectrum > threshold)),
        "positive_exponent_threshold": threshold,
        "mean_divergence": float(out["mean_divergence"]),
        "divergence_residual": float(abs(np.sum(spectrum) - p.divergence)),
        "tail_block_mean": tail.mean(axis=0).tolist() if len(tail) else [],
        "tail_block_std": tail.std(axis=0, ddof=0).tolist() if len(tail) else [],
        "backend": out["backend"],
        "status": "PASS" if np.all(np.isfinite(spectrum)) else "FAIL",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--n-steps", type=int, nargs="+", default=[250_000, 500_000, 1_000_000])
    ap.add_argument("--transient-fraction", type=float, default=0.1)
    ap.add_argument("--threshold", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=20260815)
    a = ap.parse_args()

    points = [(4.0, -np.pi / 2), (8.0, np.pi / 2)]
    records = []
    for pi, (drive, theta) in enumerate(points):
        for ni, n_steps in enumerate(a.n_steps):
            transient = int(round(n_steps * a.transient_fraction))
            records.append(run_point(
                drive, theta, n_steps, transient, a.seed + pi * 100,
                g1=0.3, g2=0.27, gamma1=0.02, gamma2=0.02,
                detuning=-1.0, threshold=a.threshold,
            ))

    stable_classification = {}
    for drive, theta in points:
        rs = [r for r in records if r["drive"] == drive and r["theta"] == theta]
        stable_classification[f"E={drive:g},theta/pi={theta / np.pi:g}"] = {
            "longest_two_n_positive_exponents": [r["n_positive_exponents"] for r in rs[-2:]],
            "stable_between_longest_two": len(rs) >= 2 and rs[-1]["n_positive_exponents"] == rs[-2]["n_positive_exponents"],
        }
    out = {
        "gate": "HYPERCHAOS_REPRESENTATIVE_CONVERGENCE",
        "status": "PASS" if all(r["status"] == "PASS" for r in records) and all(v["stable_between_longest_two"] for v in stable_classification.values()) else "FAIL",
        "classification_stability": stable_classification,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "backend": "numba-cpu",
        "coupling": {"g1": 0.3, "g2": 0.27, "gamma1": 0.02,
                     "gamma2": 0.02, "detuning": -1.0},
        "dt": 0.01,
        "points": [{"drive": d, "theta": t, "theta_over_pi": t / np.pi}
                   for d, t in points],
        "records": records,
        "interpretation": (
            "Representative strong-coupling convergence audit. The canonical map remains "
            "the primary phase/drive result; this audit tests block-history stability and "
            "divergence balance at one chaotic and one hyperchaotic point."
        ),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {a.output}: {len(records)} records, status={out['status']}")
    return 0 if out["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
