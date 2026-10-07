#!/usr/bin/env python3
"""Local susceptibility--stretching covariation on the chaotic attractor.

Refinement of the null matched-Fisher result on the chaotic attractor
(SciPost referee demand iii): *why* does the synthetic-flux phase reshape the
hyperchaos order without changing the matched measurement Fisher information?

For each declared point (drive E, flux phase theta) the script integrates one
long trajectory with the fixed-step RK4/QR backend of ``lyapunov_numba`` and,
alongside the state, the force-response vector

    dw/dt = J(x(t)) w + e_F ,   e_F = (0,0,0,1,0,0),   w = 0 at each window start,

which is the exact linear response of the state to a constant weak force F
added to the beta_1 momentum equation (``reconstruction_core.rhs_phase``, the
``+ p.force`` term).  Because w is reset at every window boundary, the
window-averaged observable susceptibility

    S_k = d<X_a>/dF |_window k = mean_t w_0(t) over window k

is a *local* (window-restricted) linear-response coefficient, free of the
finite-difference nonlinearity of a two-trajectory estimator.

Per window k the script also records
  * lambda_k  : the six local finite-time Lyapunov exponents (QR block
                accumulators over the window) -- the local stretching;
  * sigma_k   : the standard deviation of the record X_a = Re(alpha) inside
                the window -- the local record spread ("noise");
  * S_k       : the local force-to-observable susceptibility ("signal").

The matched-Fisher lemma's shared-stretching premise requires signal and noise
to respond to the same local factor lambda_k.  The script reports the Spearman
covariations of lambda_k with |S_k| and sigma_k, the log-log response slopes of
ln|S_k| and ln sigma_k against lambda_k, the distribution of the local Fisher
proxy f_k = S_k^2 / sigma_k^2, and two-sample KS tests of that proxy across
flux phases at fixed drive.  The numbers are descriptive; no pass/fail claim
is attached to their signs.

All results are model-level; no quantum or SQL claim is made.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from numba import njit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lyapunov_numba import _jac, _rhs, parameters_array  # noqa: E402
from reconstruction_core import ModelParameters  # noqa: E402

E_FORCE = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float64)


@njit(cache=True)
def _rhs_f(x, phi, p, force):
    """Model right-hand side with the weak force F on the beta_1 momentum."""
    out = _rhs(x, phi, p)
    out[3] += force
    return out


@njit(cache=True)
def _run(x0, p, force, n_steps, transient_steps, qr_interval, window_steps, dt):
    """Integrate state + reset-per-window force response + QR stretching.

    Returns per-window records: (6 local FTLEs, mean X_a, std X_a, mean w_0).
    """
    x = x0.copy()
    w = np.zeros(6, dtype=np.float64)
    Q = np.eye(6)
    accum = np.zeros(6, dtype=np.float64)
    phi = 0.0
    two_pi = 2.0 * np.pi

    n_windows = (n_steps - transient_steps) // window_steps
    win_lam = np.empty((n_windows, 6), dtype=np.float64)
    win_mean_ar = np.empty(n_windows, dtype=np.float64)
    win_std_ar = np.empty(n_windows, dtype=np.float64)
    win_mean_w = np.empty(n_windows, dtype=np.float64)

    sum_ar = 0.0
    sum_ar2 = 0.0
    sum_w = 0.0
    n_in_win = 0
    win_count = 0

    for step in range(n_steps):
        # --- RK4 on the augmented state (x, w) -----------------------------
        k1x = _rhs_f(x, phi, p, force)
        k1w = _jac(x, p).dot(w) + E_FORCE
        x2 = x + 0.5 * dt * k1x
        phi2 = phi + 0.5 * p[12] * dt
        k2x = _rhs_f(x2, phi2, p, force)
        k2w = _jac(x2, p).dot(w + 0.5 * dt * k1w) + E_FORCE
        x3 = x + 0.5 * dt * k2x
        k3x = _rhs_f(x3, phi2, p, force)
        k3w = _jac(x3, p).dot(w + 0.5 * dt * k2w) + E_FORCE
        x4 = x + dt * k3x
        phi4 = phi + p[12] * dt
        k4x = _rhs_f(x4, phi4, p, force)
        k4w = _jac(x4, p).dot(w + dt * k3w) + E_FORCE
        x_new = x + dt * (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0
        w_new = w + dt * (k1w + 2.0 * k2w + 2.0 * k3w + k4w) / 6.0

        # --- QR tangent flow (same formulas as lyapunov_numba._integrate) ---
        amat = _jac(x, p)
        amat_new = _jac(x_new, p)
        qtrial = Q + 0.5 * dt * (amat @ Q + amat_new @ (Q + dt * amat @ Q))

        x = x_new
        w = w_new
        phi = (phi + p[12] * dt) % two_pi

        if step < transient_steps:
            Q = np.eye(6)
            accum[:] = 0.0
            w[:] = 0.0
            continue

        sum_ar += x[0]
        sum_ar2 += x[0] * x[0]
        sum_w += w[0]
        n_in_win += 1

        if (step - transient_steps + 1) % qr_interval == 0:
            Q, R = np.linalg.qr(qtrial)
            for j in range(6):
                if R[j, j] < 0.0:
                    Q[:, j] *= -1.0
                val = abs(R[j, j])
                if val < np.finfo(np.float64).tiny:
                    val = np.finfo(np.float64).tiny
                accum[j] += np.log(val)
        else:
            Q = qtrial

        if (step - transient_steps + 1) % window_steps == 0:
            win_time = window_steps * dt
            mean_ar = sum_ar / n_in_win
            var_ar = max(sum_ar2 / n_in_win - mean_ar * mean_ar, 0.0)
            win_lam[win_count, :] = accum / win_time
            win_mean_ar[win_count] = mean_ar
            win_std_ar[win_count] = np.sqrt(var_ar)
            win_mean_w[win_count] = sum_w / n_in_win
            win_count += 1
            # reset for the next window: new linear-response initial condition
            sum_ar = 0.0
            sum_ar2 = 0.0
            sum_w = 0.0
            n_in_win = 0
            accum[:] = 0.0
            w[:] = 0.0

    return (win_lam[:win_count], win_mean_ar[:win_count],
            win_std_ar[:win_count], win_mean_w[:win_count])


def analyse_point(drive: float, theta_over_pi: float, seed: int, args) -> dict:
    """Run one (drive, phase) point and summarise the window records."""
    from scipy.stats import spearmanr

    p = ModelParameters(drive=drive, theta=theta_over_pi * np.pi,
                        g1=0.3, g2=0.27, gamma1=0.02, gamma2=0.02,
                        detuning=-1.0)
    x0 = np.random.default_rng(seed).normal(0.0, 0.1, 6)
    lam, mean_ar, std_ar, mean_w = _run(
        x0, parameters_array(p), args.force, args.n_steps, args.transient_steps,
        args.qr_interval, args.window_steps, args.dt)

    sig = np.abs(mean_w)
    noise = std_ar
    ok = (np.all(np.isfinite(lam)) and np.all(np.isfinite(sig))
          and np.all(np.isfinite(noise)) and np.all(sig > 0.0)
          and np.all(noise > 0.0))

    def _rho(a, b):
        r, pv = spearmanr(a, b)
        return {"spearman_rho": float(r), "p_value": float(pv)}

    with np.errstate(divide="ignore"):
        f_k = (mean_w ** 2) / (noise ** 2)  # local Fisher proxy
    lam1 = lam[:, 0]
    f_cv = float(np.std(f_k) / np.mean(f_k))
    ok = ok and np.isfinite(f_cv)
    slope_logS = float(np.polyfit(lam1, np.log(sig), 1)[0])
    slope_logs = float(np.polyfit(lam1, np.log(noise), 1)[0])

    summary = {
        "n_windows": int(len(lam)),
        "window_time_units": float(args.window_steps * args.dt),
        "rho_lambda_absS": _rho(lam1, sig),
        "rho_lambda_sigma": _rho(lam1, noise),
        "rho_absS_sigma": _rho(sig, noise),
        "slope_log_absS_vs_lambda": slope_logS,
        "slope_log_sigma_vs_lambda": slope_logs,
        "fisher_proxy": {
            "median": float(np.median(f_k)),
            "q25": float(np.percentile(f_k, 25)),
            "q75": float(np.percentile(f_k, 75)),
            "cv": f_cv,
        },
        "fraction_negative_response": float(np.mean(mean_w < 0.0)),
        "leading_exponent_window_median": float(np.median(lam1)),
        "leading_exponent_window_mean": float(np.mean(lam1)),
        "status": "PASS" if ok else "FAIL",
    }
    return {
        "drive": drive,
        "theta_over_pi": theta_over_pi,
        "seed": seed,
        "summary": summary,
        "windows": {
            "lambda_1": lam1.tolist(),
            "lambda_all": lam.tolist(),
            "mean_ar": mean_ar.tolist(),
            "sigma_ar": std_ar.tolist(),
            "susceptibility_S": mean_w.tolist(),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--points", type=str, default="4:0,4:0.5,8:0,8:0.5",
                    help="Comma-separated drive:theta_over_pi pairs")
    ap.add_argument("--n-steps", type=int, default=1_000_000)
    ap.add_argument("--transient-steps", type=int, default=100_000)
    ap.add_argument("--window-steps", type=int, default=2000)
    ap.add_argument("--qr-interval", type=int, default=10)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--force", type=float, default=1e-6,
                    help="Not used for the response itself (w is the exact "
                         "unit-force linear response); kept as an explicit "
                         "record of the linearisation scale.")
    ap.add_argument("--seed", type=int, default=20260815)
    a = ap.parse_args()

    points = []
    for pair in a.points.split(","):
        d_str, t_str = pair.strip().split(":")
        points.append((float(d_str), float(t_str)))

    results = []
    for pi, (drive, theta) in enumerate(points):
        results.append(analyse_point(drive, theta, a.seed + pi * 7, a))
        s = results[-1]["summary"]
        print(f"E={drive:g} theta/pi={theta:g}: windows={s['n_windows']} "
              f"rho(lam,|S|)={s['rho_lambda_absS']['spearman_rho']:+.3f} "
              f"rho(lam,sigma)={s['rho_lambda_sigma']['spearman_rho']:+.3f} "
              f"rho(|S|,sigma)={s['rho_absS_sigma']['spearman_rho']:+.3f} "
              f"f_k cv={s['fisher_proxy']['cv']:.3f} status={s['status']}")

    # Cross-phase comparison of the local Fisher-proxy distribution: does the
    # flux phase change the per-window signal-to-noise budget?
    from scipy.stats import ks_2samp
    logf = {}
    for r in results:
        w = r["windows"]
        f = (np.abs(np.array(w["susceptibility_S"])) ** 2
             / np.array(w["sigma_ar"]) ** 2)
        logf[(r["drive"], r["theta_over_pi"])] = np.log(f)
    drives = sorted({r["drive"] for r in results})
    thetas = sorted({r["theta_over_pi"] for r in results})
    cross_phase = []
    for d in drives:
        for i in range(len(thetas)):
            for j in range(i + 1, len(thetas)):
                ks = ks_2samp(logf[(d, thetas[i])], logf[(d, thetas[j])])
                cross_phase.append({
                    "drive": d,
                    "theta_a": thetas[i], "theta_b": thetas[j],
                    "ks_statistic": float(ks.statistic),
                    "p_value": float(ks.pvalue),
                    "median_log_f_a": float(np.median(logf[(d, thetas[i])])),
                    "median_log_f_b": float(np.median(logf[(d, thetas[j])])),
                })
                print(f"KS log f  E={d:g} theta {thetas[i]:g} vs {thetas[j]:g}: "
                      f"stat={ks.statistic:.3f} p={ks.pvalue:.3g}")

    status = ("PASS" if all(r["summary"]["status"] == "PASS" for r in results)
              else "FAIL")
    out = {
        "gate": "LOCAL_FISHER_STRETCHING",
        "kind": "local_susceptibility_stretching_covariation",
        "status": status,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "backend": "numba-cpu",
        "coupling": {"g1": 0.3, "g2": 0.27, "gamma1": 0.02,
                     "gamma2": 0.02, "detuning": -1.0},
        "method": {
            "response": ("exact linear response: w' = J(x) w + e_F with "
                         "e_F = (0,0,0,1,0,0), w reset to zero at every "
                         "window boundary; S_k = window mean of w_0"),
            "noise": "window standard deviation of the record X_a = Re(alpha)",
            "stretching": "QR block accumulators of the tangent flow per window",
            "fisher_proxy": "f_k = S_k^2 / sigma_k^2 (local, window-restricted)",
            "dt": a.dt,
            "n_steps": a.n_steps,
            "transient_steps": a.transient_steps,
            "window_steps": a.window_steps,
            "qr_interval": a.qr_interval,
        },
        "results": results,
        "cross_phase_tests": cross_phase,
        "interpretation": (
            "Window-resolved covariation of the local force-to-observable "
            "susceptibility, the local record spread, and the local finite-time "
            "Lyapunov stretching along the chaotic attractor, plus log-log "
            "response slopes of both quantities against the stretching and "
            "two-sample tests of the local Fisher-proxy distribution across "
            "flux phases at fixed drive. The matched-Fisher lemma's "
            "shared-stretching premise requires signal and noise to be "
            "functions of the same local stretching; the reported correlations "
            "and slopes quantify how strongly each responds. This is a "
            "model-level classical linear-response statement on the tested "
            "trajectory ensemble; it is not a quantum, SQL, or device claim.",
        ),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(f"wrote {a.output}: status={status}")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
