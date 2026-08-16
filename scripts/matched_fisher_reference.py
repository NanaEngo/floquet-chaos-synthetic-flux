#!/usr/bin/env python3
"""Matched measurement-Fisher reference for optomechanical force sensing.

Computes the classical measurement Fisher information F_C(F) for estimating a
weak external force F on mechanical mode 1, read out via the cavity amplitude
quadrature X_a = Re(alpha) (the same record as the covariance/Fisher pilot).

Configurations (all matched in input power/drive, observation time = one drive
period, bandwidth/detector noise, thermal baths, and estimator):
  * synthetic-flux sensor: theta in {0, pi/4, pi/2, 3pi/4, pi}, hopping = 0.08
  * reference A (flux off): theta = 0,  hopping = 0.08  (same coupled system)
  * reference B (linear single-mode): hopping = 0  (J = 0, uncoupled modes)

This is a classical measurement Fisher information, not QFI; the only thing that
varies between a flux configuration and its reference is the synthetic-flux phase
(or the inter-resonator coupling for the single-mode reference).
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

from reconstruction_core import ModelParameters, find_periodic_orbit
from covariance_fisher import periodic_covariance


def fisher_force(theta: float, hopping: float, eps: float, *,
                 n_th1: float = 0.1, n_th2: float = 0.1,
                 detector_variance: float = 0.01, drive: float = 0.2) -> dict:
    """Central-difference F_C for a weak force F on mode 1, at fixed (theta, J)."""
    def obs(F: float):
        p = ModelParameters(theta=float(theta), hopping=float(hopping),
                            force=float(F), drive=float(drive))
        orbit = find_periodic_orbit(np.zeros(6), p, max_iter=600)
        if orbit.get("status") != "PASS":
            raise RuntimeError(orbit.get("reason", "periodic orbit failed"))
        return periodic_covariance(np.asarray(orbit["x0"]), p,
                                   n_th1=n_th1, n_th2=n_th2,
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
        "theta": float(theta), "hopping": float(hopping), "force_eps": eps,
        "d_mean_dF": float(dm), "d_variance_dF": float(dv),
        "mean_variance": float(var),
        "classical_fisher_information": fc,
        "min_physicality_eigenvalue": float(phys),
        "status": "PASS" if (phys >= -1e-8 and np.isfinite(fc)) else "FAIL",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--thetas", nargs="+", type=float,
                    default=[0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    ap.add_argument("--eps", type=float, default=1e-4)
    ap.add_argument("--eps-check", nargs="*", type=float,
                    default=[1e-3, 1e-4, 1e-5])
    ap.add_argument("--hopping", type=float, default=0.08)
    ap.add_argument("--n-th1", type=float, default=0.1)
    ap.add_argument("--n-th2", type=float, default=0.1)
    ap.add_argument("--detector-variance", type=float, default=0.01)
    ap.add_argument("--drive", type=float, default=0.2)
    a = ap.parse_args()

    noise = dict(n_th1=a.n_th1, n_th2=a.n_th2, detector_variance=a.detector_variance)

    # Flux scan
    records = []
    for th in a.thetas:
        try:
            rec = fisher_force(th, a.hopping, a.eps, drive=a.drive, **noise)
            records.append(rec)
        except Exception as exc:
            records.append({"theta": float(th), "hopping": a.hopping,
                            "status": "FAIL",
                            "reason": f"{type(exc).__name__}: {exc}"})

    # Reference A: flux off (theta = 0, same hopping)
    refA = fisher_force(0.0, a.hopping, a.eps, drive=a.drive, **noise)

    # Reference B: single-mode linear reference (hopping = 0)
    refB = fisher_force(0.0, 0.0, a.eps, drive=a.drive, **noise)

    # eps convergence check on the theta = pi/2 record (primary reported value)
    eps_check = []
    for e in a.eps_check:
        try:
            r = fisher_force(np.pi / 2, a.hopping, e, drive=a.drive, **noise)
            eps_check.append({"eps": e, "fisher": r["classical_fisher_information"],
                              "status": r["status"]})
        except Exception as exc:
            eps_check.append({"eps": e, "status": "FAIL",
                              "reason": f"{type(exc).__name__}: {exc}"})

    # Gains
    gains = []
    for rec in records:
        if rec.get("status") != "PASS":
            gains.append({"theta": rec.get("theta"), "status": "FAIL"})
            continue
        gA = rec["classical_fisher_information"] / refA["classical_fisher_information"] if refA["classical_fisher_information"] != 0 else float("nan")
        gB = rec["classical_fisher_information"] / refB["classical_fisher_information"] if refB["classical_fisher_information"] != 0 else float("nan")
        gains.append({"theta": rec["theta"], "gain_vs_flux_off": float(gA),
                      "gain_vs_single_mode": float(gB)})

    out = {
        "gate": "MATCHED_MEASUREMENT_FISHER_REFERENCE",
        "kind": "classical_measurement_fisher_force_sensing",
        "status": "PASS" if all(r.get("status") == "PASS" for r in records) else "FAIL",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "sensing_parameter": "external force F on mechanical mode 1 (momentum-quadrature drive d(beta1_i)/dt += F)",
        "measurement": {
            "observable": "cavity amplitude quadrature X_a = Re(alpha)",
            "record": "y(t) = X_a(t) + nu_det(t)",
            "detector_noise_variance": a.detector_variance,
            "thermal_occupations": [a.n_th1, a.n_th2],
            "observation_time": "one drive period (matched across all configurations)",
        },
        "matched_resources": {
            "drive": a.drive, "observation_time": "1 drive period",
            "bandwidth": "detector noise variance (matched)",
            "baths": [a.n_th1, a.n_th2], "estimator": "eq:fisher (mean + variance terms)",
        },
        "force_eps": a.eps,
        "eps_convergence_check_theta_pi_over_2": eps_check,
        "references": {
            "flux_off_theta_0_hopping_0.08": refA,
            "single_mode_hopping_0": refB,
        },
        "records": records,
        "gains": gains,
        "interpretation": "Matched classical measurement Fisher information for force sensing. A gain > 1 relative to both references would indicate a flux-mediated sensing advantage; a gain ~ 1 or < 1 indicates no advantage. Not QFI and not a claim of quantum advantage.",
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": out["status"], "output": str(a.output),
                      "records": len(records)}, indent=2))
    return 0 if out["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
