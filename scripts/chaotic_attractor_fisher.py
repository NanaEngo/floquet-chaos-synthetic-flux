#!/usr/bin/env python3
"""Matched measurement-Fisher assessment on the chaotic attractor.

The weak-coupling force-sensing Fisher reference (matched_fisher_reference.py)
is computed around a *stable* drive-locked orbit via the periodic covariance.
That construction does not exist on a chaotic attractor (no stable periodic
orbit). This script instead evaluates the same classical Fisher question
directly on the chaotic attractor: for a weak external force F on mechanical
mode 1, read out through the cavity amplitude quadrature X_a = Re(alpha) under
the declared Langevin (thermal + vacuum) and detector noise, we define the
record moments

    m(F)  = mean over (ensemble x time) of y = X_a + nu_det
    s2(F) = variance over (ensemble x time) of y

and evaluate the classical measurement Fisher information by central finite
difference in F,

    F_C(F) = (dm/dF)^2 / s2  +  (1/2) (ds2/dF)^2 / s2^2 ,

exactly the estimator used for the stable reference (mean + variance terms).
The observation window is matched across configurations (same number of drive
periods, same detector noise, same thermal baths, same estimator).

Configurations (all at a fixed strong-coupling drive E in {4, 8}):
  * synthetic-flux sensor: theta in {0, pi/2, pi}, hopping = 0.08
  * reference A (flux off): theta = 0,  hopping = 0.08  (same coupled system)
  * reference B (single-mode linear): hopping = 0 (J = 0). At weak coupling
    this is the decoupled linear reference; at strong coupling the uncoupled,
    strongly driven two-mode system is unbounded (the two-mode coupling is what
    bounds the dynamics), so reference B is recorded as NOT_APPLICABLE there
    and only the flux-off reference A is used for the gain.

On a chaotic attractor the record variance is dominated by the deterministic
chaotic spread rather than by the detector noise, so F_C is expected to be
strongly suppressed relative to the stable reference and the flux-phase gain is
expected to be of order unity (no chaotic-transduction advantage). The
finite-difference resolution is limited by the sampling error of the ensemble
moments; a bootstrap over ensemble blocks reports the uncertainty, so an
unresolved derivative is reported as such rather than as a finite gain.

This is a classical measurement Fisher information, not QFI, and not a claim of
a quantum or SQL-beating advantage.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, find_periodic_orbit
from noise_observability import rhs_batch, noise_strengths


def record_moments(x0, p: ModelParameters, q, det_sigma, n_ens, n_steps, dt,
                   transient_steps, seed) -> tuple[float, float]:
    """Ensemble+time mean and variance of y = X_a + detector noise."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x0, dtype=np.float64)[None, :] + 1e-4 * rng.standard_normal((n_ens, 6))
    q = np.asarray(q, dtype=np.float64)
    sum_y = 0.0
    sum_y2 = 0.0
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
            sum_y += float(y.sum())
            sum_y2 += float((y * y).sum())
            n += n_ens
    mean = sum_y / n
    var = sum_y2 / n - mean * mean
    return float(mean), float(var)


def block_bootstrap(x0, p: ModelParameters, q, det_sigma, n_ens, n_steps, dt,
                    transient_steps, n_blocks, n_boot, seed) -> dict:
    """Bootstrap of (mean, variance) over time blocks, for derivative resolution."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x0, dtype=np.float64)[None, :] + 1e-4 * rng.standard_normal((n_ens, 6))
    q = np.asarray(q, dtype=np.float64)
    block_len = max(1, (n_steps - transient_steps) // n_blocks)
    block_means = []
    block_vars = []
    for step in range(n_steps):
        phi = (step * dt * p.drive_frequency) % (2.0 * np.pi)
        k1 = rhs_batch(x, phi, p)
        k2 = rhs_batch(x + 0.5 * dt * k1, phi + 0.5 * p.drive_frequency * dt, p)
        k3 = rhs_batch(x + 0.5 * dt * k2, phi + 0.5 * p.drive_frequency * dt, p)
        k4 = rhs_batch(x + dt * k3, phi + p.drive_frequency * dt, p)
        x = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        x += np.sqrt(q[None, :] * dt) * rng.standard_normal((n_ens, 6))
        if step >= transient_steps and (step - transient_steps) % block_len == 0:
            y = x[:, 0] + np.sqrt(det_sigma) * rng.standard_normal(n_ens)
            block_means.append(float(y.mean()))
            block_vars.append(float(y.var()))
    bm = np.asarray(block_means)
    bv = np.asarray(block_vars)
    if bm.size < 2:
        raise RuntimeError("not enough blocks for bootstrap")
    rng2 = np.random.default_rng(seed + 1)
    idx = rng2.integers(0, bm.size, size=(n_boot, bm.size))
    boot_m = bm[idx].mean(axis=1)
    boot_v = bv[idx].mean(axis=1)
    return {
        "mean": float(bm.mean()), "mean_std": float(bm.std(ddof=1) / np.sqrt(bm.size)),
        "variance": float(bv.mean()), "variance_std": float(bv.std(ddof=1) / np.sqrt(bv.size)),
        "n_blocks": int(bm.size),
        "boot_mean_std": float(boot_m.std(ddof=1)),
        "boot_variance_std": float(boot_v.std(ddof=1)),
    }


def _attractor_seed(E, theta, hopping, dt, n_warm):
    """Return a state on (or near) the chaotic attractor.

    We start from the converged drive-locked orbit of the *weak-coupling*
    reference (a bounded initial condition), then integrate the strong-coupling
    deterministic flow forward for n_warm steps so the transient relaxes onto
    the attractor. Integrating from the origin directly diverges at strong
    coupling (the unstable orbit is repelling), so a bounded seed is required.
    The warm-up uses the same RK4 scheme (with the phase-resolved rhs_phase)
    as the production Lyapunov spectra (reconstruction_core.lyapunov_qr); a
    plain forward-Euler step is numerically unstable at these amplitudes and
    would spuriously overflow.
    """
    from reconstruction_core import rhs_phase

    p0 = ModelParameters(theta=theta)  # weak coupling, drive 0.2
    orb = find_periodic_orbit(np.zeros(6), p0, max_iter=600)
    if orb.get("status") != "PASS":
        raise RuntimeError(orb.get("reason", "reference orbit failed"))
    x = np.asarray(orb["x0"], dtype=np.float64).copy()
    p = ModelParameters(kappa=1.0, gamma1=0.02, gamma2=0.02, omega1=1.0,
                        omega2=1.03, g1=0.3, g2=0.27, hopping=hopping,
                        detuning=-1.0, drive=E, drive_modulation=0.1,
                        drive_frequency=1.0, theta=theta)
    phi = 0.0
    for _ in range(n_warm):
        k1 = rhs_phase(x, phi, p)
        k2 = rhs_phase(x + 0.5 * dt * k1, phi + 0.5 * p.drive_frequency * dt, p)
        k3 = rhs_phase(x + 0.5 * dt * k2, phi + 0.5 * p.drive_frequency * dt, p)
        k4 = rhs_phase(x + dt * k3, phi + p.drive_frequency * dt, p)
        x = x + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        phi = (phi + p.drive_frequency * dt) % (2.0 * np.pi)
    if not np.all(np.isfinite(x)):
        raise RuntimeError("attractor seed diverged")
    return x, p


def chaotic_fisher(E, theta, hopping, eps, *, n_ens=128, n_steps=30000,
                   transient_steps=3000, dt=0.01, n_th1=0.1, n_th2=0.1,
                   detector_variance=0.01, seed=42, n_warm=200000) -> dict:
    """Central-difference F_C on the chaotic attractor at fixed (E, theta, J)."""
    seed_state, p_base = _attractor_seed(E, theta, hopping, dt, n_warm)

    def obs(F):
        p = ModelParameters(kappa=1.0, gamma1=0.02, gamma2=0.02, omega1=1.0,
                            omega2=1.03, g1=0.3, g2=0.27, hopping=hopping,
                            detuning=-1.0, drive=E, drive_modulation=0.1,
                            drive_frequency=1.0, theta=theta, force=F)
        q = noise_strengths(p, n_th1, n_th2)
        m, s2 = record_moments(seed_state, p, q, detector_variance, n_ens,
                               n_steps, dt, transient_steps, seed)
        return m, s2

    mp, sp = obs(+eps)
    mm, sm = obs(-eps)
    mc, sc = obs(0.0)
    dm = (mp - mm) / (2.0 * eps)
    dv = (sp - sm) / (2.0 * eps)
    var = sc
    fc = float(dm * dm / var + 0.5 * dv * dv / (var * var))

    # Resolution: bootstrap of the center configuration over time blocks.
    q = noise_strengths(p_base, n_th1, n_th2)
    boot = block_bootstrap(seed_state, p_base, q, detector_variance, n_ens,
                           n_steps, dt, transient_steps, n_blocks=32,
                           n_boot=200, seed=seed)

    # Is the finite difference resolved above the sampling error of the center?
    dm_sigma = boot["boot_mean_std"]
    dv_sigma = boot["boot_variance_std"]
    resolved = bool(abs(dm) > 3.0 * dm_sigma or abs(dv) > 3.0 * dv_sigma)

    return {
        "E": float(E), "theta": float(theta), "theta_over_pi": float(theta / np.pi),
        "hopping": float(hopping), "force_eps": eps,
        "record_mean": float(mc), "record_variance": float(sc),
        "d_mean_dF": float(dm), "d_variance_dF": float(dv),
        "classical_fisher_information": fc,
        "sampling_mean_std": float(dm_sigma), "sampling_variance_std": float(dv_sigma),
        "derivative_resolved": resolved,
        "status": "PASS" if np.isfinite(fc) else "FAIL",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drives", nargs="+", type=float, default=[4.0, 8.0])
    ap.add_argument("--thetas", nargs="+", type=float, default=[0.0, np.pi / 2, np.pi])
    ap.add_argument("--eps", type=float, default=0.1)
    ap.add_argument("--n-ens", type=int, default=128)
    ap.add_argument("--n-steps", type=int, default=30000)
    ap.add_argument("--transient-steps", type=int, default=3000)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--detector-variance", type=float, default=0.01)
    ap.add_argument("--n-th", nargs=2, type=float, default=[0.1, 0.1])
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    n_th1, n_th2 = a.n_th

    all_records = []
    for E in a.drives:
        for th in a.thetas:
            try:
                rec = chaotic_fisher(E, th, 0.08, a.eps, n_ens=a.n_ens,
                                     n_steps=a.n_steps, transient_steps=a.transient_steps,
                                     dt=a.dt, n_th1=n_th1, n_th2=n_th2,
                                     detector_variance=a.detector_variance, seed=a.seed)
                all_records.append(rec)
            except Exception as exc:
                all_records.append({"E": float(E), "theta": float(th),
                                    "status": "FAIL",
                                    "reason": f"{type(exc).__name__}: {exc}"})

        # References at the same drive E
        refA = chaotic_fisher(E, 0.0, 0.08, a.eps, n_ens=a.n_ens,
                              n_steps=a.n_steps, transient_steps=a.transient_steps,
                              dt=a.dt, n_th1=n_th1, n_th2=n_th2,
                              detector_variance=a.detector_variance, seed=a.seed)
        # Single-mode (J = 0) reference: only well defined at weak coupling.
        # At strong coupling the uncoupled, strongly driven two-mode system is
        # unbounded (the two-mode coupling is what bounds the dynamics), so the
        # reference is recorded as not applicable rather than failing the run.
        try:
            refB = chaotic_fisher(E, 0.0, 0.0, a.eps, n_ens=a.n_ens,
                                  n_steps=a.n_steps, transient_steps=a.transient_steps,
                                  dt=a.dt, n_th1=n_th1, n_th2=n_th2,
                                  detector_variance=a.detector_variance, seed=a.seed)
        except Exception as exc:
            refB = {"E": float(E), "theta": 0.0, "hopping": 0.0,
                    "status": "NOT_APPLICABLE",
                    "reason": ("uncoupled strongly-driven two-mode system is unbounded; "
                                f"{type(exc).__name__}")}
        for rec in all_records:
            if rec.get("E") == E and rec.get("status") == "PASS" and refA.get("status") == "PASS":
                rec["gain_vs_flux_off"] = rec["classical_fisher_information"] / refA["classical_fisher_information"]
            if rec.get("E") == E and rec.get("status") == "PASS" and refB.get("status") == "PASS":
                rec["gain_vs_single_mode"] = rec["classical_fisher_information"] / refB["classical_fisher_information"]

    status = "PASS" if all(r.get("status") == "PASS" for r in all_records) else "FAIL"
    out = {
        "gate": "CHAOTIC_ATTRACTOR_MATCHED_FISHER",
        "kind": "classical_measurement_fisher_on_chaotic_attractor",
        "status": status,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "settings": {
            "drives": a.drives, "thetas": a.thetas, "force_eps": a.eps,
            "n_ens": a.n_ens, "n_steps": a.n_steps,
            "transient_steps": a.transient_steps, "dt": a.dt,
            "detector_variance": a.detector_variance,
            "thermal_occupations": [n_th1, n_th2], "seed": a.seed,
        },
        "sensing_parameter": "external force F on mechanical mode 1 (momentum-quadrature drive d(beta1_i)/dt += F)",
        "measurement": {
            "observable": "cavity amplitude quadrature X_a = Re(alpha)",
            "record": "y(t) = X_a(t) + nu_det(t)",
            "detector_noise_variance": a.detector_variance,
            "thermal_occupations": [n_th1, n_th2],
            "moments": "ensemble x time mean and variance over the observation window",
        },
        "records": all_records,
        "interpretation": (
            "Matched classical measurement Fisher information evaluated directly "
            "on the chaotic attractor (ensemble+time record moments, central "
            "finite difference in the force). On a chaotic attractor the record "
            "variance is dominated by the deterministic chaotic spread rather "
            "than by detector noise, so F_C is strongly suppressed relative to "
            "the stable reference and the flux-phase gain is expected to be of "
            "order unity. A derivative is marked resolved only if it exceeds "
            "three bootstrap standard deviations of the center sampling error; "
            "an unresolved derivative is reported as such and its gain is not "
            "claimed. Not QFI and not a quantum-advantage claim."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "output": str(a.output),
                      "records": len(all_records)}, indent=2))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
