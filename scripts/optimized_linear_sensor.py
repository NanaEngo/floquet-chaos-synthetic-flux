#!/usr/bin/env python3
"""Optimized linear-sensor Fisher reference for the matched force-sensing bound.

The matched Fisher pilot (matched_fisher_reference.py) compares the synthetic-flux
sensor against two references that all share the *same* readout, the cavity
amplitude quadrature X_a = Re(alpha).  A severe reviewer can object that a fixed
readout is not the best possible *linear* sensor: the force signal may live in a
rotated homodyne quadrature X_phi = Re(alpha) cos phi + Im(alpha) sin phi, so a
single-mode (J = 0) reference read out at an optimal angle phi* is a stronger
baseline than the matched X_a readout.

This script computes the classical measurement Fisher information F_C(F) for the
external force F on mechanical mode 1, for an arbitrary cavity quadrature angle
phi, reusing the exact Gaussian-likelihood estimator of the pilot

    F_C = (d mu / dF)^2 / V  +  (1/2) (dV / dF)^2 / V^2,

with mu = time-averaged rotated signal and V = time-averaged rotated cavity
variance plus the detector noise (the phi = 0 case reduces exactly to the X_a
readout of the pilot).  It scans phi over [0, pi) for three configurations:

  * single-mode linear reference (J = 0)          -- the optimized linear sensor;
  * flux-off reference (theta = 0, J = 0.08);
  * flux-on sensor  (theta = pi/2, J = 0.08).

The result is the optimized readout angle and F_C for each configuration, plus
the flux gain against the *optimized* single-mode reference (the honest bound
on the sensing advantage).  This is classical measurement Fisher information,
not QFI, and no quantum-advantage claim is made.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from reconstruction_core import (ModelParameters, find_periodic_orbit,
                                 integrate_one_period, jacobian_phase, monodromy)


def cavity_statistics(p: ModelParameters, *, n_th1: float = 0.1,
                      n_th2: float = 0.1) -> dict:
    """Return the cavity mean mu(t) = [ar(t), ai(t)] and covariance V_cav(t)
    over one drive period, around the drive-locked orbit of ``p``.

    Returns
    -------
    dict with ``times`` (shape (N,)), ``mu`` (shape (2, N)) and ``V_cav``
    (shape (N, 2, 2)).
    """
    orbit = find_periodic_orbit(np.zeros(6), p, max_iter=600)
    if orbit.get("status") != "PASS":
        raise RuntimeError(orbit.get("reason", "periodic orbit failed"))
    x0 = np.asarray(orbit["x0"], dtype=float)

    dense = integrate_one_period(x0, p, dense=True)
    if not dense.success or dense.sol is None:
        raise RuntimeError("orbit interpolation failed")

    n = 6
    q = np.array([p.kappa / 2, p.kappa / 2,
                  p.gamma1 * (2 * n_th1 + 1) / 2, p.gamma1 * (2 * n_th1 + 1) / 2,
                  p.gamma2 * (2 * n_th2 + 1) / 2, p.gamma2 * (2 * n_th2 + 1) / 2],
                 dtype=float)
    Q = np.diag(q)

    def fun(t, v):
        x = dense.sol(t)
        A = jacobian_phase(x, p.drive_frequency * t, p)
        V = v.reshape(n, n, order="F")
        return (A @ V + V @ A.T + Q).ravel(order="F")

    zero = solve_ivp(fun, (0.0, p.period), np.zeros(n * n), method="DOP853",
                     rtol=1e-9, atol=1e-11, max_step=p.period / 200.0)
    if not zero.success:
        raise RuntimeError(zero.message)
    qT = zero.y[:, -1]
    M = np.asarray(monodromy(x0, p)["monodromy"], dtype=float)
    K = np.eye(n * n) - np.kron(M, M)
    v0 = np.linalg.solve(K, qT)
    V0 = v0.reshape(n, n, order="F")
    V0 = 0.5 * (V0 + V0.T)

    full = solve_ivp(fun, (0.0, p.period), V0.ravel(order="F"), method="DOP853",
                     rtol=1e-9, atol=1e-11, dense_output=True,
                     max_step=p.period / 200.0)
    if not full.success or full.sol is None:
        raise RuntimeError(full.message)

    times = np.linspace(0.0, p.period, 1001)
    mu = np.asarray([dense.sol(t) for t in times]).T          # (6, N)
    Vs = np.asarray([full.sol(t).reshape(n, n, order="F") for t in times])  # (N,6,6)
    return {
        "times": times,
        "mu": mu[0:2, :],              # (2, N): [ar(t), ai(t)]
        "V_cav": Vs[:, 0:2, 0:2],      # (N, 2, 2): cavity covariance sub-block
    }


def fisher_at_angle(stats_plus: dict, stats_minus: dict, stats_center: dict,
                    detector_variance: float, angle: float, eps: float) -> dict:
    """Fisher information for the rotated cavity quadrature at ``angle``.

    Uses central differences in the force F with step ``eps``; the phi = 0 case
    reproduces the X_a = Re(alpha) readout of the matched pilot exactly.
    """
    w = np.array([np.cos(angle), np.sin(angle)])

    def mean_signal(s):
        return float(np.mean(np.einsum("i,it->t", w, s["mu"])))

    def mean_variance(s):
        var = np.einsum("i,nij,j->n", w, s["V_cav"], w) + detector_variance
        return float(np.mean(var))

    mu_c = mean_signal(stats_center)
    d_mu = (mean_signal(stats_plus) - mean_signal(stats_minus)) / (2.0 * eps)
    var_c = mean_variance(stats_center)
    d_var = (mean_variance(stats_plus) - mean_variance(stats_minus)) / (2.0 * eps)
    fc = d_mu * d_mu / var_c + 0.5 * d_var * d_var / (var_c * var_c)
    return {
        "angle": float(angle),
        "angle_over_pi": float(angle / np.pi),
        "mean_signal": float(mu_c),
        "mean_variance": float(var_c),
        "d_mean_dF": float(d_mu),
        "d_variance_dF": float(d_var),
        "classical_fisher_information": float(fc),
    }


def scan_configuration(theta: float, hopping: float, *, drive: float, eps: float,
                       n_th1: float, n_th2: float, detector_variance: float,
                       n_angles: int = 360) -> dict:
    """Scan the readout angle for one (theta, hopping) configuration."""
    def stats(F):
        p = ModelParameters(theta=theta, hopping=hopping, force=F, drive=drive)
        return cavity_statistics(p, n_th1=n_th1, n_th2=n_th2)

    plus = stats(+eps)
    minus = stats(-eps)
    center = stats(0.0)

    angles = np.linspace(0.0, np.pi, n_angles, endpoint=False)
    recs = [fisher_at_angle(plus, minus, center, detector_variance, a, eps)
            for a in angles]
    best = max(recs, key=lambda r: r["classical_fisher_information"])
    xa = fisher_at_angle(plus, minus, center, detector_variance, 0.0, eps)
    return {
        "theta": float(theta), "theta_over_pi": float(theta / np.pi),
        "hopping": float(hopping),
        "Xa_readout_fisher": xa["classical_fisher_information"],
        "optimized_readout": best,
        "scan": recs,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--eps", type=float, default=1e-4)
    ap.add_argument("--n-th1", type=float, default=0.1)
    ap.add_argument("--n-th2", type=float, default=0.1)
    ap.add_argument("--detector-variance", type=float, default=0.01)
    ap.add_argument("--drive", type=float, default=0.2)
    ap.add_argument("--hopping", type=float, default=0.08)
    ap.add_argument("--n-angles", type=int, default=360)
    ap.add_argument("--thetas", nargs="+", type=float,
                    default=[0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    a = ap.parse_args()

    noise = dict(n_th1=a.n_th1, n_th2=a.n_th2, detector_variance=a.detector_variance)

    # The optimized single-mode (J=0) linear reference: the strongest linear baseline.
    single = scan_configuration(0.0, 0.0, drive=a.drive, eps=a.eps, n_angles=a.n_angles,
                                **noise)
    single_opt = single["optimized_readout"]["classical_fisher_information"]

    # Every flux phase, read out at the matched X_a angle and at its optimum.
    flux_phases = []
    for th in a.thetas:
        cfg = scan_configuration(float(th), a.hopping, drive=a.drive, eps=a.eps,
                                 n_angles=a.n_angles, **noise)
        flux_phases.append({
            "theta": float(th),
            "theta_over_pi": float(th / np.pi),
            "Xa_readout_fisher": cfg["Xa_readout_fisher"],
            "optimized_angle_over_pi": cfg["optimized_readout"]["angle_over_pi"],
            "optimized_fisher": cfg["optimized_readout"]["classical_fisher_information"],
            "gain_matched_Xa_vs_optimized_linear":
                float(cfg["Xa_readout_fisher"] / single_opt),
            "gain_optimized_vs_optimized_linear":
                float(cfg["optimized_readout"]["classical_fisher_information"] / single_opt),
        })

    max_matched = max(flux_phases, key=lambda r: r["gain_matched_Xa_vs_optimized_linear"])
    max_optimized = max(flux_phases, key=lambda r: r["gain_optimized_vs_optimized_linear"])

    out = {
        "gate": "OPTIMIZED_LINEAR_SENSOR_REFERENCE",
        "kind": "classical_measurement_fisher_optimized_readout",
        "status": "PASS",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "sensing_parameter": "external force F on mechanical mode 1",
        "estimator": "Gaussian-likelihood Fisher: (d mu/dF)^2/V + 0.5 (dV/dF)^2/V^2",
        "readout": "rotated cavity quadrature X_phi = Re(alpha) cos phi + Im(alpha) sin phi",
        "parameters": {"drive": a.drive, "eps": a.eps, "n_angles": a.n_angles,
                       "n_th1": a.n_th1, "n_th2": a.n_th2,
                       "detector_variance": a.detector_variance,
                       "hopping": a.hopping},
        "single_mode_linear_reference": {
            "Xa_readout_fisher": single["Xa_readout_fisher"],
            "optimized_angle_over_pi": single["optimized_readout"]["angle_over_pi"],
            "optimized_fisher": single["optimized_readout"]["classical_fisher_information"],
        },
        "flux_phases": flux_phases,
        "max_gain_matched_Xa_vs_optimized_linear": max_matched,
        "max_gain_optimized_vs_optimized_linear": max_optimized,
        "interpretation": (
            "The single-mode (J=0) linear reference read out at its optimal "
            "homodyne angle is the strongest non-chaotic linear baseline. A flux "
            "gain close to or below unity against this optimized reference "
            "removes the 'suboptimal reference readout' objection to the "
            "order-unity Fisher gain. Classical measurement Fisher information "
            "only; not QFI and no quantum-advantage claim."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({
        "status": "PASS",
        "single_Xa": single["Xa_readout_fisher"],
        "single_opt_fisher": single_opt,
        "single_opt_angle_pi": single["optimized_readout"]["angle_over_pi"],
        "flux_phases": [
            {"theta_over_pi": r["theta_over_pi"],
             "gain_matched": r["gain_matched_Xa_vs_optimized_linear"],
             "gain_optimized": r["gain_optimized_vs_optimized_linear"]}
            for r in flux_phases],
        "output": str(a.output),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
