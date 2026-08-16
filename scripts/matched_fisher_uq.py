#!/usr/bin/env python3
"""Uncertainty quantification of the matched force-sensing Fisher gains.

Gap 2 of GAP_IMPLEMENTATION_20260816.md. The matched measurement-Fisher gains
(1.323 vs flux-off, 1.159 vs single-mode) are deterministic: one hopping, one
bath occupation, one drive. This script propagates uncertainty in the three
*unmeasured* coordinates to the gain:

  * inter-resonator hopping J_m  -- log-uniform over the declared scan range
    [1e-3, 1e-1] (Mathew et al. scale), because no two-resonator GHz value exists;
  * mechanical bath occupation n_th -- uniform over [0.01, 1.0], spanning
    cryogenic to ~4 K at the anchored omega_m = 2pi*7.436 GHz;
  * drive amplitude E -- uniform over +/-10% around the nominal 0.2.

The measured coordinates (kappa, gamma, g0) and the detector noise are held at
their reference values; the gain is evaluated at the maximum-gain phase
theta = pi. The drive-locked orbit is found with a damped-Newton solver on the
Poincare map (monodromy as Jacobian) for speed. Classical measurement Fisher
information, not QFI.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, integrate_one_period, monodromy
from covariance_fisher import periodic_covariance

ROOT = Path(__file__).resolve().parent.parent


def newton_orbit(p: ModelParameters, *, max_iter: int = 40,
                 tol: float = 1e-8) -> np.ndarray:
    """Damped Newton on F(x) = P(x) - x with J = M - I (monodromy Jacobian)."""
    x = np.zeros(6)
    for _ in range(max_iter):
        sol = integrate_one_period(x, p)
        x1 = np.asarray(sol.y[:, -1], dtype=float)
        F = x1 - x
        if float(np.linalg.norm(F, ord=np.inf)) < tol:
            return x1
        M = np.asarray(monodromy(x, p)["monodromy"], dtype=float)
        J = M - np.eye(6)
        try:
            delta = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            delta = np.linalg.pinv(J) @ (-F)
        x = x + delta
    raise RuntimeError("Newton orbit solve did not converge")


def fisher_force_fast(theta: float, hopping: float, eps: float, *,
                      n_th1: float, n_th2: float,
                      detector_variance: float, drive: float) -> dict:
    """Matched classical measurement Fisher information with a Newton orbit solve."""
    def obs(F: float):
        p = ModelParameters(theta=float(theta), hopping=float(hopping),
                            force=float(F), drive=float(drive))
        x0 = newton_orbit(p)
        return periodic_covariance(x0, p, n_th1=n_th1, n_th2=n_th2,
                                   detector_variance=detector_variance)

    plus = obs(+eps)
    minus = obs(-eps)
    center = obs(0.0)
    dm = (plus["mean_signal"] - minus["mean_signal"]) / (2.0 * eps)
    dv = (plus["mean_variance"] - minus["mean_variance"]) / (2.0 * eps)
    var = center["mean_variance"]
    fc = float(dm * dm / var + 0.5 * dv * dv / (var * var))
    phys = min(center["min_quantum_physicality_eigenvalue"],
               plus["min_quantum_physicality_eigenvalue"],
               minus["min_quantum_physicality_eigenvalue"])
    return {
        "classical_fisher_information": fc,
        "min_physicality_eigenvalue": float(phys),
        "status": "PASS" if (phys >= -1e-8 and np.isfinite(fc)) else "FAIL",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--n-samples", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--eps", type=float, default=1e-4)
    ap.add_argument("--theta", type=float, default=np.pi)
    ap.add_argument("--hopping-range", nargs=2, type=float, default=[1e-3, 1e-1])
    ap.add_argument("--n-th-range", nargs=2, type=float, default=[0.01, 1.0])
    ap.add_argument("--drive-range", nargs=2, type=float, default=[0.18, 0.22])
    ap.add_argument("--detector-variance", type=float, default=0.01)
    ap.add_argument("--omega-m-ghz", type=float, default=7.436)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    h_lo, h_hi = a.hopping_range
    n_lo, n_hi = a.n_th_range
    d_lo, d_hi = a.drive_range

    base_hopping = 0.08
    base_sensor = fisher_force_fast(a.theta, base_hopping, a.eps,
                                    n_th1=0.1, n_th2=0.1,
                                    detector_variance=a.detector_variance, drive=0.2)
    base_refA = fisher_force_fast(0.0, base_hopping, a.eps,
                                  n_th1=0.1, n_th2=0.1,
                                  detector_variance=a.detector_variance, drive=0.2)
    base_refB = fisher_force_fast(0.0, 0.0, a.eps,
                                  n_th1=0.1, n_th2=0.1,
                                  detector_variance=a.detector_variance, drive=0.2)
    baseline = {
        "gain_vs_flux_off": base_sensor["classical_fisher_information"] / base_refA["classical_fisher_information"],
        "gain_vs_single_mode": base_sensor["classical_fisher_information"] / base_refB["classical_fisher_information"],
    }

    records = []
    for i in range(a.n_samples):
        hopping = float(np.exp(rng.uniform(np.log(h_lo), np.log(h_hi))))
        n_th = float(rng.uniform(n_lo, n_hi))
        drive = float(rng.uniform(d_lo, d_hi))
        try:
            sensor = fisher_force_fast(a.theta, hopping, a.eps,
                                       n_th1=n_th, n_th2=n_th,
                                       detector_variance=a.detector_variance, drive=drive)
            refA = fisher_force_fast(0.0, hopping, a.eps,
                                     n_th1=n_th, n_th2=n_th,
                                     detector_variance=a.detector_variance, drive=drive)
            refB = fisher_force_fast(0.0, 0.0, a.eps,
                                     n_th1=n_th, n_th2=n_th,
                                     detector_variance=a.detector_variance, drive=drive)
            gA = sensor["classical_fisher_information"] / refA["classical_fisher_information"] if refA["classical_fisher_information"] != 0 else float("nan")
            gB = sensor["classical_fisher_information"] / refB["classical_fisher_information"] if refB["classical_fisher_information"] != 0 else float("nan")
            records.append({
                "sample": i, "hopping": hopping, "n_th": n_th, "drive": drive,
                "gain_vs_flux_off": float(gA), "gain_vs_single_mode": float(gB),
                "status": "PASS" if (np.isfinite(gA) and np.isfinite(gB)) else "FAIL",
            })
            print(f"sample {i}/{a.n_samples} done (hop={hopping:.4f}, n_th={n_th:.3f}, drive={drive:.3f})", flush=True)
        except Exception as exc:
            records.append({"sample": i, "hopping": hopping, "n_th": n_th,
                            "drive": drive, "status": "FAIL",
                            "reason": f"{type(exc).__name__}: {exc}"})
            print(f"sample {i}/{a.n_samples} FAIL ({type(exc).__name__})", flush=True)

    ok = [r for r in records if r.get("status") == "PASS"]
    gA = np.array([r["gain_vs_flux_off"] for r in ok])
    gB = np.array([r["gain_vs_single_mode"] for r in ok])

    def stats(x):
        if x.size == 0:
            return {"n": 0}
        return {
            "n": int(x.size),
            "median": float(np.median(x)),
            "p05": float(np.percentile(x, 5)),
            "p95": float(np.percentile(x, 95)),
            "fraction_above_unity": float(np.mean(x > 1.0)),
        }

    out = {
        "gate": "MATCHED_FISHER_UQ",
        "kind": "classical_measurement_fisher_uncertainty_quantification",
        "status": "PASS" if ok else "FAIL",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "seed": a.seed,
        "sensing_phase": a.theta,
        "sampling": {
            "hopping_loguniform": [h_lo, h_hi],
            "bath_occupation_uniform": [n_lo, n_hi],
            "drive_uniform": [d_lo, d_hi],
            "bath_temperature_note": (
                f"n_th in [{n_lo},{n_hi}] spans cryogenic to ~4 K at "
                f"omega_m = 2pi*{a.omega_m_ghz} GHz"),
        },
        "baseline_deterministic": baseline,
        "gain_vs_flux_off": stats(gA),
        "gain_vs_single_mode": stats(gB),
        "records": records,
        "interpretation": (
            "The order-unity matched force-sensing gain is propagated over the "
            "unmeasured coordinates (hopping, bath occupation, drive). The gain "
            "is robust if the Monte Carlo median remains near the deterministic "
            "value and the 5th percentile stays above or near unity; it is not "
            "robust if the spread pushes the gain below unity for a substantial "
            "fraction of samples. Classical Fisher information, not QFI."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": out["status"], "n_ok": len(ok),
                      "gain_vs_flux_off": stats(gA),
                      "gain_vs_single_mode": stats(gB),
                      "output": str(a.output)}, indent=2))
    return 0 if out["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
