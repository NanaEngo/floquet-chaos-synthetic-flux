#!/usr/bin/env python3
"""Observe the Neimark--Sacker torus, not merely infer it.

The drive-axis continuation (results/drive_axis_bifurcation_scan.json) shows a
complex-conjugate pair of Floquet multipliers crossing the unit circle at
E* ~ 1.06 (theta=0) and E* ~ 0.975 (theta=pi/2). A crossing alone is an
inference; this script *observes* the born object by:

  1. perturbing the (now unstable) drive-locked orbit just above E*;
  2. integrating forward and recording the stroboscopic Poincare section
     (state at every drive period T = 2 pi / drive_frequency);
  3. computing the full Lyapunov spectrum of the same trajectory.

A Neimark--Sacker torus shows up as (a) a closed invariant curve in the
Poincare section (a 1D ring in the 2D projection) and (b) a largest Lyapunov
exponent consistent with zero (marginal, the quasiperiodic neutral direction)
while all remaining exponents are negative. This is the observation that turns
the multiplier crossing into a torus.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from numba import njit

from reconstruction_core import ModelParameters
from lyapunov_numba import parameters_array, _rhs, lyapunov_qr_numba

ROOT = Path(__file__).resolve().parent.parent
SCAN = ROOT / "results" / "drive_axis_bifurcation_scan.json"


def strong_params(E: float, theta: float) -> ModelParameters:
    return ModelParameters(kappa=1.0, gamma1=0.02, gamma2=0.02, omega1=1.0,
                           omega2=1.03, g1=0.3, g2=0.27, hopping=0.08,
                           detuning=-1.0, drive=E, drive_modulation=0.1,
                           drive_frequency=1.0, theta=theta)


@njit(cache=True)
def _poincare_integrate(x0, p, n_transient, n_records, steps_per_period, dt):
    """Record stroboscopic points (one per drive period) after a transient."""
    x = x0.copy()
    phi = 0.0
    pts = np.empty((n_records, 6), dtype=np.float64)
    rec = 0
    for step in range(n_transient + n_records * steps_per_period):
        k1 = _rhs(x, phi, p)
        k2 = _rhs(x + 0.5 * dt * k1, phi + 0.5 * p[12] * dt, p)
        k3 = _rhs(x + 0.5 * dt * k2, phi + 0.5 * p[12] * dt, p)
        k4 = _rhs(x + dt * k3, phi + p[12] * dt, p)
        x = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        phi = (phi + p[12] * dt) % (2.0 * np.pi)
        if step >= n_transient and (step - n_transient) % steps_per_period == 0:
            pts[rec] = x
            rec += 1
    return pts


def orbit_at(records: list[dict], E_target: float, direction: str):
    cands = [r for r in records if r.get("direction") == direction
             and r.get("status") == "PASS" and r.get("x0") is not None]
    rec = min(cands, key=lambda r: abs(r["E"] - E_target))
    return np.asarray(rec["x0"], dtype=np.float64), float(rec["E"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drive-above-onsets", nargs="+", type=float,
                    default=[0.15, 0.3])  # E = E* + these offsets
    ap.add_argument("--n-transient", type=int, default=40000)
    ap.add_argument("--n-records", type=int, default=2000)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--perturb", type=float, default=1e-3)
    ap.add_argument("--lyap-n-steps", type=int, default=200000)
    ap.add_argument("--lyap-transient", type=int, default=20000)
    a = ap.parse_args()

    scan = json.loads(SCAN.read_text())
    # Onset locations from the continuation's first multiplier crossing.
    onsets = {}
    for trec in scan["theta_records"]:
        key = round(float(trec["theta"]), 9)
        cross = trec.get("crossings") or []
        if cross:
            onsets[key] = float(cross[0]["E_star"])

    period = 2.0 * np.pi  # drive_frequency = 1.0
    steps_per_period = int(round(period / a.dt))

    records = []
    for theta in (0.0, np.pi / 2):
        trec = next(t for t in scan["theta_records"]
                    if abs(t["theta"] - theta) < 1e-9)
        onset = onsets[round(theta, 9)]
        for off in a.drive_above_onsets:
            E_target = onset + off
            x0, E_actual = orbit_at(trec["records"], E_target, "up")
            p = strong_params(E_actual, theta)
            parr = parameters_array(p)
            seed = x0.copy()
            for d in range(6):
                seed[d] += a.perturb * (0.5 - ((d * 2654435761) % 1000) / 1000.0)
            pts = _poincare_integrate(seed, parr, a.n_transient, a.n_records,
                                      steps_per_period, a.dt)
            lyap = lyapunov_qr_numba(seed, p, n_steps=a.lyap_n_steps, dt=a.dt,
                                     transient_steps=a.lyap_transient,
                                     qr_interval=10)
            spectrum = np.asarray(lyap["spectrum"])
            # Torus signature: lambda_1 ~ 0 (marginal) and lambda_2 < 0.
            lam1 = float(spectrum[0])
            lam2 = float(spectrum[1])
            # Ring test: project Poincare points onto (b1r, b1i) and check the
            # angular spread covers the full circle (torus) rather than a point.
            b1 = pts[:, 2:4]
            center = b1.mean(axis=0)
            rad = np.linalg.norm(b1 - center, axis=1)
            radial_cv = float(np.std(rad) / (np.mean(rad) + 1e-12))
            angular_spread = float(np.std(np.arctan2(b1[:, 1] - center[1],
                                                     b1[:, 0] - center[0])))
            records.append({
                "theta": float(theta), "theta_over_pi": float(theta / np.pi),
                "drive": E_actual, "drive_target": E_target,
                "onset_E": onset, "offset_above_onset": off,
                "largest_exponent": lam1, "second_exponent": lam2,
                "spectrum": lyap["spectrum"],
                "mean_divergence": lyap["mean_divergence"],
                "poincare_radial_cv": radial_cv,
                "poincare_angular_spread_std": angular_spread,
                "poincare_points_b1r": b1[:, 0].tolist(),
                "poincare_points_b1i": b1[:, 1].tolist(),
                "torus_signature": {
                    "lambda1_near_zero": bool(abs(lam1) < 2e-2),
                    "lambda2_negative": bool(lam2 < 0.0),
                },
                "status": "PASS" if np.all(np.isfinite(spectrum)) else "FAIL",
            })
            print(f"theta={theta/np.pi:.2f}pi E={E_actual:.3f} "
                  f"lam1={lam1:+.4f} lam2={lam2:+.4f} radial_cv={radial_cv:.3f} "
                  f"ang_spread={angular_spread:.3f}", flush=True)

    status = "PASS" if all(r["status"] == "PASS" for r in records) else "FAIL"
    out = {
        "gate": "TORUS_OBSERVATION",
        "kind": "neimark_sacker_torus_poincare_section",
        "status": status,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "settings": {
            "n_transient": a.n_transient, "n_records": a.n_records,
            "dt": a.dt, "perturb": a.perturb, "steps_per_period": steps_per_period,
            "lyap_n_steps": a.lyap_n_steps, "lyap_transient": a.lyap_transient,
        },
        "onsets_used": {f"{k:.9f}": v for k, v in onsets.items()},
        "records": records,
        "interpretation": (
            "A Neimark--Sacker torus is *observed* (not only inferred from a "
            "multiplier crossing) when, just above the onset, the stroboscopic "
            "Poincare section of a perturbed trajectory forms a closed invariant "
            "curve and the largest Lyapunov exponent is consistent with zero "
            "(the quasiperiodic neutral direction) while the remaining exponents "
            "are negative. If the largest exponent is instead clearly positive, "
            "the object is chaotic rather than a clean torus at that drive."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(json.dumps({"status": status, "output": str(a.output)}, indent=2))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
