#!/usr/bin/env python3
"""Noise-dependent observability of the chaos signature in the measured record.

Refinement R1 of NOVELTY_REFINEMENT_ANALYSIS_20260816.md, implementing
P4_P5_MATHEMATICAL_DERIVATIONS.md §7.3: a map of whether the deterministic
chaotic signature remains observable after thermal, vacuum, and detector noise
are propagated.

For each operating point (stable reference, chaotic, hyperchaotic) we integrate
an ensemble of N trajectories under the declared Langevin noise
(optical vacuum kappa/4 per amplitude quadrature, mechanical thermal
gamma_j (2 n_th + 1)/4 per amplitude quadrature, both with n_th = 0.1) with the
measurement record y = X_a + nu_det (detector variance 0.01), and report the
measured-record variance as the declared classifier. The deterministic spread
(no noise) is reported alongside, so the chaos signature is observable whenever
the deterministic spread exceeds the noise-only floor set by the stable
reference.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, find_periodic_orbit

ROOT = Path(__file__).resolve().parent.parent
SCAN = ROOT / "results" / "drive_axis_bifurcation_scan.json"


def rhs_batch(x, phi, p: ModelParameters) -> np.ndarray:
    """Vectorized right-hand side for an (N, 6) state batch."""
    ar, ai, b1r, b1i, b2r, b2i = x[:, 0], x[:, 1], x[:, 2], x[:, 3], x[:, 4], x[:, 5]
    c, s = np.cos(p.theta), np.sin(p.theta)
    omega_eff = p.detuning + 2.0 * p.g1 * b1r + 2.0 * p.g2 * b2r
    e = p.drive * (1.0 + p.drive_modulation * np.cos(phi))
    out = np.empty_like(x)
    out[:, 0] = -0.5 * p.kappa * ar - omega_eff * ai + e
    out[:, 1] = omega_eff * ar - 0.5 * p.kappa * ai
    out[:, 2] = -0.5 * p.gamma1 * b1r + p.omega1 * b1i + p.hopping * (s * b2r + c * b2i)
    out[:, 3] = p.g1 * (ar * ar + ai * ai) - p.omega1 * b1r - 0.5 * p.gamma1 * b1i - p.hopping * c * b2r + p.hopping * s * b2i + p.force
    out[:, 4] = -0.5 * p.gamma2 * b2r + p.omega2 * b2i - p.hopping * s * b1r + p.hopping * c * b1i
    out[:, 5] = p.g2 * (ar * ar + ai * ai) - p.omega2 * b2r - 0.5 * p.gamma2 * b2i - p.hopping * c * b1r - p.hopping * s * b1i
    return out


def ensemble_variance(x0, p: ModelParameters, q, n_ens, n_steps, dt,
                      det_sigma, seed, transient_steps):
    """Ensemble+time variance of the measured record y = X_a + detector noise."""
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
    return float(var)


def noise_strengths(p: ModelParameters, n_th1: float, n_th2: float) -> np.ndarray:
    """Diffusion diagonal Q for the Langevin equation (per amplitude quadrature).

    The state is x = (Re alpha, Im alpha, Re beta_1, Im beta_1, Re beta_2,
    Im beta_2), i.e. the *amplitude* quadrature (a+a^dag)/2 of each mode. That
    quadrature has vacuum variance 1/4 and thermal variance (2 n_th + 1)/4, so
    the correct diffusion coefficients are kappa/4 (optical) and
    gamma_j (2 n_th_j + 1)/4 (mechanical) -- NOT kappa/2, which is the
    diffusion of the *normalized* quadrature (a+a^dag)/sqrt(2) (variance 1/2).
    This normalization is fixed by results/langevin_master_crosscheck.json:
    kappa/4 reproduces the master-equation variance of (a+a^dag)/2 to <2%,
    while kappa/2 is a factor 2 too large.
    """
    return np.array([
        p.kappa / 4.0, p.kappa / 4.0,
        p.gamma1 * (2.0 * n_th1 + 1.0) / 4.0,
        p.gamma1 * (2.0 * n_th1 + 1.0) / 4.0,
        p.gamma2 * (2.0 * n_th2 + 1.0) / 4.0,
        p.gamma2 * (2.0 * n_th2 + 1.0) / 4.0,
    ], dtype=np.float64)


def strong_params(E: float, theta: float) -> ModelParameters:
    return ModelParameters(kappa=1.0, gamma1=0.02, gamma2=0.02, omega1=1.0,
                           omega2=1.03, g1=0.3, g2=0.27, hopping=0.08,
                           detuning=-1.0, drive=E, drive_modulation=0.1,
                           drive_frequency=1.0, theta=theta)


def orbit_at(records: list[dict], E_target: float, direction: str):
    cands = [r for r in records if r.get("direction") == direction
             and r.get("status") == "PASS" and r.get("x0") is not None]
    rec = min(cands, key=lambda r: abs(r["E"] - E_target))
    return np.asarray(rec["x0"], dtype=np.float64), float(rec["E"])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--n-ens", type=int, default=32)
    ap.add_argument("--n-steps", type=int, default=50000)
    ap.add_argument("--transient-steps", type=int, default=5000)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--detector-variance", type=float, default=0.01)
    ap.add_argument("--n-th", nargs=2, type=float, default=[0.1, 0.1])
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()

    n_th1, n_th2 = a.n_th
    scan = json.loads(SCAN.read_text())

    points = []
    for th in (0.0, np.pi / 2):
        p = ModelParameters(theta=th)  # weak coupling, drive=0.2
        orb = find_periodic_orbit(np.zeros(6), p, max_iter=600)
        if orb.get("status") != "PASS":
            raise RuntimeError(orb.get("reason", "reference orbit failed"))
        points.append((f"stable_E0.2_th{th/np.pi:.2f}pi", p,
                       np.asarray(orb["x0"], dtype=np.float64)))

    for E in (4.0, 8.0):
        for th in (0.0, np.pi / 2):
            trec = next(t for t in scan["theta_records"]
                        if abs(t["theta"] - th) < 1e-9)
            x0, E_actual = orbit_at(trec["records"], E, "up")
            p = strong_params(E_actual, th)
            label = "hyperchaotic" if E == 8.0 else "chaotic"
            points.append((f"{label}_E{E_actual:.2f}_th{th/np.pi:.2f}pi", p, x0))

    records = []
    for label, p, x0 in points:
        q = noise_strengths(p, n_th1, n_th2)
        var_det = ensemble_variance(x0, p, np.zeros(6), a.n_ens, a.n_steps,
                                    a.dt, 0.0, a.seed, a.transient_steps)
        var_noisy = ensemble_variance(x0, p, q, a.n_ens, a.n_steps, a.dt,
                                      a.detector_variance, a.seed, a.transient_steps)
        records.append({
            "label": label,
            "drive": float(p.drive),
            "theta": float(p.theta),
            "theta_over_pi": float(p.theta / np.pi),
            "coupling": "weak" if p.g1 < 0.05 else "strong",
            "deterministic_spread_var": var_det,
            "measured_record_var": var_noisy,
            "status": "PASS" if np.isfinite(var_noisy) else "FAIL",
        })

    stable_noisy = [r for r in records if r["label"].startswith("stable")]
    noise_floor = float(np.mean([r["measured_record_var"] for r in stable_noisy]))
    for r in records:
        r["observability_ratio_vs_noise_floor"] = r["measured_record_var"] / noise_floor

    status = "PASS" if all(r["status"] == "PASS" for r in records) else "FAIL"
    out = {
        "gate": "NOISE_DEPENDENT_OBSERVABILITY",
        "kind": "noise_resolved_chaos_observability",
        "status": status,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "settings": {
            "n_ens": a.n_ens, "n_steps": a.n_steps,
            "transient_steps": a.transient_steps, "dt": a.dt,
            "detector_variance": a.detector_variance,
            "thermal_occupations": [n_th1, n_th2], "seed": a.seed,
        },
        "classifier": "measured-record variance of y = X_a + detector noise",
        "noise_floor_measured_record_var": noise_floor,
        "records": records,
        "interpretation": (
            "The deterministic spread of the cavity amplitude quadrature is "
            "compared with the measured-record variance under the declared "
            "Langevin (thermal + vacuum) and detector noise. A chaos signature "
            "is observable whenever the deterministic spread exceeds the "
            "noise-only floor set by the stable reference (observability ratio "
            "well above unity); at the stable reference the ratio is unity by "
            "construction. Deterministic spread alone is not a metrological "
            "gain and not QFI."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "noise_floor": noise_floor,
                      "output": str(a.output)}, indent=2))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
