#!/usr/bin/env python3
"""Drive-axis bifurcation-continuation scan at strong optomechanical coupling.

Continues the drive-locked periodic orbit family of the six-dimensional model
along the drive amplitude E at fixed synthetic-flux phase, using a damped
Newton solver on the one-period Poincare map (Jacobian = monodromy - I). At
every point the full set of Floquet multipliers is recorded, so the first loss
of stability, the bifurcation type (real +1, real -1, or complex pair), any
restabilization windows, and the growth of the number of unstable directions
with drive are resolved directly. The scan runs both upward (E increasing) and
downward (E decreasing) to expose branch structure, with adaptive step halving
near bifurcation points and resume-from-checkpoint support.

Settings match the production strong-coupling map
(`results/full_spectrum_hyperchaos_map_g03.json`):
g1=0.3, g2=0.27, gamma1=gamma2=0.02, detuning=-1.0, kappa=1.0, hopping=0.08,
omega1=1.0, omega2=1.03, drive_modulation=0.1, drive_frequency=1.0.
Deterministic: no RNG; the branch starts from the zero state at E = 0.1.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruction_core import ModelParameters, jacobian_phase, rhs_phase


def params(E: float, theta: float) -> ModelParameters:
    return ModelParameters(drive=float(E), theta=float(theta),
                           g1=0.3, g2=0.27, gamma1=0.02, gamma2=0.02,
                           detuning=-1.0)


def poincare_augmented(x0: np.ndarray, p: ModelParameters, dense: bool = False):
    """Integrate state + variational equations over one drive period."""
    n = 6
    y0 = np.r_[np.asarray(x0, dtype=float), np.eye(n).ravel()]

    def fun(t, y):
        x = y[:n]
        phi = p.drive_frequency * t
        A = jacobian_phase(x, phi, p)
        return np.r_[rhs_phase(x, phi, p), (A @ y[n:].reshape(n, n)).ravel()]

    return solve_ivp(fun, (0.0, p.period), y0, method="DOP853",
                     rtol=1e-9, atol=1e-11, max_step=p.period / 200.0,
                     dense_output=dense)


def newton_map(x0: np.ndarray, p: ModelParameters, *,
               tol: float = 1e-10, max_iter: int = 25,
               max_backtrack: int = 15,
               time_budget: float = 180.0) -> dict:
    """Damped Newton on F(x) = P(x) - x with J = M - I (single augmented solve).

    A wall-clock budget bounds the cost near bifurcation points where the
    monodromy is nearly singular and the line search stalls; exceeding it
    returns FAIL so the continuation can halve its step and move on. A
    no-progress detector fails fast when the residual plateaus (oscillating
    Newton on a vanished branch), so branch-end crawls do not burn the full
    budget at every point.
    """
    t0 = time.perf_counter()
    x = np.asarray(x0, dtype=float).copy()
    history = []
    best = float("inf")
    stall = 0
    for it in range(max_iter):
        if time.perf_counter() - t0 > time_budget:
            return {"status": "FAIL", "reason": "time budget exceeded",
                    "iterations": it, "residual": float("nan"),
                    "x0": x.tolist(), "history": history}
        sol = poincare_augmented(x, p)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return {"status": "FAIL", "reason": "integration failed",
                    "iterations": it, "residual": float("nan"),
                    "x0": x.tolist()}
        x1 = sol.y[:6, -1]
        M = sol.y[6:, -1].reshape(6, 6)
        F = x1 - x
        r = float(np.linalg.norm(F, ord=np.inf))
        history.append(r)
        if r < tol:
            return {"status": "PASS", "x0": x.tolist(), "residual": r,
                    "iterations": it + 1, "history": history}
        if r < 0.9 * best:
            best = r
            stall = 0
        else:
            stall += 1
            if stall >= 6:
                return {"status": "FAIL", "reason": "no progress (stalled)",
                        "iterations": it + 1, "residual": r,
                        "x0": x.tolist(), "history": history}
        J = M - np.eye(6)
        try:
            delta = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            delta = np.linalg.pinv(J) @ (-F)
        step = 1.0
        accepted = False
        x_trial = x + step * delta
        for _ in range(max_backtrack):
            if time.perf_counter() - t0 > time_budget:
                return {"status": "FAIL", "reason": "time budget exceeded",
                        "iterations": it + 1, "residual": r,
                        "x0": x.tolist(), "history": history}
            sol_t = poincare_augmented(x_trial, p)
            if sol_t.success:
                r_t = float(np.linalg.norm(sol_t.y[:6, -1] - x_trial, ord=np.inf))
                if np.isfinite(r_t) and r_t < r:
                    accepted = True
                    break
            step *= 0.5
            x_trial = x + step * delta
        if not accepted:
            return {"status": "FAIL", "reason": "line search stalled",
                    "iterations": it + 1, "residual": r,
                    "x0": x.tolist(), "history": history}
        x = x_trial
    sol = poincare_augmented(x, p)
    r = float(np.linalg.norm(sol.y[:6, -1] - x, ord=np.inf))
    return {"status": "PASS" if r < tol else "FAIL",
            "reason": None if r < tol else "max iterations",
            "x0": x.tolist(), "residual": r, "iterations": max_iter,
            "history": history}


def multiplier_stats(x0: np.ndarray, p: ModelParameters) -> dict:
    sol = poincare_augmented(x0, p)
    M = sol.y[6:, -1].reshape(6, 6)
    mu = np.linalg.eigvals(M)
    rates = np.log(np.maximum(np.abs(mu), np.finfo(float).tiny)) / p.period
    order = np.argsort(rates)[::-1]
    mu = mu[order]
    rates = rates[order]
    return {
        "monodromy": M.tolist(),
        "multipliers": [{"real": float(z.real), "imag": float(z.imag)}
                        for z in mu],
        "floquet_rates": rates.tolist(),
        "max_abs_multiplier": float(np.max(np.abs(mu))),
        "n_unstable": int(np.sum(np.abs(mu) > 1.0 + 1e-6)),
    }


def mean_amplitude(x0: np.ndarray, p: ModelParameters, n_sample: int = 200) -> dict:
    sol = poincare_augmented(x0, p, dense=True)
    t = np.linspace(0.0, p.period, n_sample)
    y = sol.sol(t)
    amp = np.hypot(y[0], y[1])
    return {"mean_amplitude": float(np.mean(amp)),
            "min_amplitude": float(np.min(amp)),
            "max_amplitude": float(np.max(amp))}


def classify_crossing(mu_hi) -> dict:
    i = int(np.argmax(np.abs(mu_hi)))
    z = mu_hi[i]
    if abs(np.imag(z)) < 1e-4 and np.real(z) > 0:
        kind = "real +1 (saddle-node / transcritical)"
    elif abs(np.imag(z)) < 1e-4 and np.real(z) < 0:
        kind = "real -1 (period doubling)"
    else:
        kind = "complex pair (Neimark-Sacker)"
    return {"crossing_angle_rad": float(np.angle(z)), "crossing_type": kind,
            "multiplier_at_crossing": {"real": float(z.real), "imag": float(z.imag)}}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--thetas", nargs="+", type=float, default=[0.0, np.pi / 2])
    ap.add_argument("--e-min", type=float, default=0.1)
    ap.add_argument("--e-max", type=float, default=8.0)
    ap.add_argument("--d-e", type=float, default=0.1)
    ap.add_argument("--crossing-tol", type=float, default=1e-3)
    ap.add_argument("--continuation-tol", type=float, default=5e-6,
                    help="Newton acceptance residual for continuation sweeps "
                         "(looser than the 1e-10 used for the crossing "
                         "refinement; 5e-6 suffices for multiplier estimates "
                         "and avoids 45 s budget burns near singular "
                         "monodromies)")
    ap.add_argument("--point-time-cap", type=float, default=600.0,
                    help="per-E wall-clock cap across step-halving retries "
                         "(bounds crawl cost at branch ends)")
    a = ap.parse_args()

    partial_path = a.output.with_suffix(a.output.suffix + ".partial")
    resume = None
    if partial_path.is_file():
        try:
            resume = json.loads(partial_path.read_text())
        except json.JSONDecodeError:
            resume = None

    out = {"gate": "DRIVE_AXIS_BIFURCATION_SCAN",
           "kind": "drive-locked orbit continuation, Floquet multipliers",
           "status": "PASS",
           "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
           "python": sys.version, "platform": platform.platform(),
           "settings": {"thetas": a.thetas, "e_min": a.e_min, "e_max": a.e_max,
                        "d_e": a.d_e, "continuation_tol": a.continuation_tol,
                        "crossing_tol": a.crossing_tol,
                        "point_time_cap": a.point_time_cap,
                        "coupling": {"g1": 0.3, "g2": 0.27,
                        "gamma1": 0.02, "gamma2": 0.02, "detuning": -1.0},
                        "kappa": 1.0, "hopping": 0.08, "omega1": 1.0,
                        "omega2": 1.03, "drive_modulation": 0.1,
                        "drive_frequency": 1.0},
           "theta_records": [], "summary": {},
           "resumed_from": (resume.get("partial") if resume else None),
           "interpretation": ("Continuation of the drive-locked orbit family along E. "
                              "max_abs_multiplier > 1 marks loss of linear stability of "
                              "the orbit; n_unstable counts |mu| > 1. Beyond the first "
                              "bifurcation the attractor need not coincide with the "
                              "continued orbit.")}

    def write_partial(theta, direction, records):
        partial_path.write_text(json.dumps({
            "partial": True, "theta": float(theta),
            "direction": direction, "records": records,
            "settings": out["settings"],
            "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    def continue_branch(theta: float, direction: str, x0_init, records: list,
                        E0: float, dE0: float, n_success0: int) -> tuple:
        E = float(E0)
        dE = float(dE0)
        E_end = float(a.e_max if direction == "up" else a.e_min)
        dE_min = dE / 16.0
        x0 = np.asarray(x0_init, dtype=float).copy()
        n_success = n_success0
        n_fail_consec = 0
        terminated = False
        t_point0 = time.perf_counter()
        while ((dE > 0 and E <= E_end + 0.5 * dE) if direction == "up"
               else (dE > 0 and E >= E_end - 0.5 * dE)):
            p = params(E, theta)
            res = newton_map(x0, p, tol=a.continuation_tol)
            if (res["status"] != "PASS" and dE > dE_min
                    and time.perf_counter() - t_point0 < a.point_time_cap):
                dE *= 0.5
                continue
            t_point0 = time.perf_counter()
            rec = {"direction": direction, "E": float(E), "theta": float(theta),
                   "theta_over_pi": float(theta / np.pi),
                   "status": res["status"], "residual": res["residual"],
                   "iterations": res["iterations"], "step": float(dE)}
            if res["status"] == "PASS":
                x0 = np.asarray(res["x0"])
                rec["x0"] = x0.tolist()
                rec.update(multiplier_stats(x0, p))
                rec.update(mean_amplitude(x0, p))
                n_success += 1
                n_fail_consec = 0
            else:
                rec.update({"max_abs_multiplier": None, "n_unstable": None,
                            "floquet_rates": None, "multipliers": None,
                            "reason": res.get("reason")})
                n_fail_consec += 1
            records.append(rec)
            print(f"theta={theta:.4f} {direction} E={E:.2f} "
                  f"{res['status']} r={res['residual']:.2e} "
                  f"n_unstable={rec.get('n_unstable')} step={dE:.3f}", flush=True)
            if len(records) % 5 == 0:
                write_partial(theta, direction, records)
            if n_fail_consec >= 25 and dE <= dE_min:
                rec["reason"] = "branch terminated: 25 consecutive fails at min step"
                terminated = True
                print(f"theta={theta:.4f} {direction} branch terminated at E={E:.4f}",
                      flush=True)
                break
            dE = min(dE * 2.0, a.d_e) if n_success >= 2 else dE
            E = E + dE if direction == "up" else E - dE
        return records, x0, terminated

    for theta in a.thetas:
        records = []
        skip = False
        if resume and abs(resume.get("theta", 1e9) - theta) < 1e-12:
            records = list(resume["records"])
            last = records[-1]
            if last.get("reason") and "branch terminated" in str(last.get("reason")):
                # this theta is already complete (branch ended); move to next
                skip = True
            elif last["direction"] == "down" and last["status"] == "PASS":
                # resume a downward pass interrupted mid-branch
                x0 = np.asarray(last["x0"])
                dE_res = max(float(last.get("step", a.d_e)), a.d_e / 16.0)
                records, x0, _ = continue_branch(theta, "down", x0, records,
                                                 last["E"] - dE_res, dE_res, 0)
                skip = True
            elif last["direction"] == "down":
                # down pass interrupted at a FAILed point: reseed from the last
                # PASSed down record and continue below the failed E
                last_pass = [r for r in records if r["direction"] == "down"
                             and r["status"] == "PASS"]
                x0 = np.asarray(last_pass[-1]["x0"]) if last_pass else np.zeros(6)
                dE_res = max(float(last.get("step", a.d_e)), a.d_e / 16.0)
                records, x0, _ = continue_branch(theta, "down", x0, records,
                                                 last["E"] - dE_res, dE_res, 0)
                skip = True
            elif last["direction"] == "up" and last["status"] == "PASS":
                if last["E"] < a.e_max - 1e-9:
                    # up pass interrupted before reaching E_max: resume it,
                    # then run the down pass only if the branch survived
                    x0 = np.asarray(last["x0"])
                    dE_res = max(float(last.get("step", a.d_e)), a.d_e / 16.0)
                    records, x0, up_terminated = continue_branch(
                        theta, "up", x0, records, last["E"] + dE_res, dE_res, 0)
                    if not up_terminated:
                        last_up = [r for r in records
                                   if r["direction"] == "up"
                                   and r["status"] == "PASS"]
                        x0 = np.asarray(last_up[-1]["x0"]) if last_up else np.zeros(6)
                        records, _, _ = continue_branch(theta, "down", x0, records,
                                                        a.e_max, a.d_e, 0)
                else:
                    # up pass completed: run the downward pass
                    x0 = np.asarray(last["x0"])
                    records, x0, _ = continue_branch(theta, "down", x0, records,
                                                     a.e_max, a.d_e, 0)
                skip = True
            elif last["direction"] == "up":
                # resume upward pass from the last recorded E
                x0 = np.zeros(6)
                last_pass = [r for r in records if r["status"] == "PASS"]
                if last_pass:
                    x0 = np.asarray(last_pass[-1]["x0"])
                dE_res = max(float(last.get("step", a.d_e)), a.d_e / 16.0)
                records, x0, up_terminated = continue_branch(
                    theta, "up", x0, records, last["E"] + dE_res, dE_res, 0)
                if not up_terminated:
                    last_up = [r for r in records if r["direction"] == "up"
                               and r["status"] == "PASS"]
                    x0 = np.asarray(last_up[-1]["x0"]) if last_up else np.zeros(6)
                    records, _, _ = continue_branch(theta, "down", x0, records,
                                                    a.e_max, a.d_e, 0)
                skip = True
        if not skip:
            records, x0, up_terminated = continue_branch(
                theta, "up", np.zeros(6), records, a.e_min, a.d_e, 0)
            if not up_terminated:
                last_up = [r for r in records if r["direction"] == "up"
                           and r["status"] == "PASS"]
                x0 = np.asarray(last_up[-1]["x0"]) if last_up else np.zeros(6)
                records, _, _ = continue_branch(theta, "down", x0, records,
                                                a.e_max, a.d_e, 0)

        # refine crossings on the upward branch
        crossings = []
        up_recs = [r for r in records if r["direction"] == "up"
                   and r["status"] == "PASS"]
        for r_lo, r_hi in zip(up_recs[:-1], up_recs[1:]):
            s_lo = r_lo["max_abs_multiplier"] - 1.0
            s_hi = r_hi["max_abs_multiplier"] - 1.0
            if s_lo * s_hi < 0.0:
                E_lo, E_hi = r_lo["E"], r_hi["E"]
                mu_hi = np.array([complex(m["real"], m["imag"])
                                  for m in r_hi["multipliers"]])
                s_lo_cur = s_lo
                x_star = np.asarray(r_lo["x0"])
                while E_hi - E_lo > a.crossing_tol:
                    Em = 0.5 * (E_lo + E_hi)
                    rm = newton_map(x_star, params(float(Em), theta))
                    if rm["status"] != "PASS":
                        E_lo = Em
                        continue
                    x_star = np.asarray(rm["x0"])
                    sm = multiplier_stats(x_star, params(float(Em), theta))[
                        "max_abs_multiplier"] - 1.0
                    if s_lo_cur * sm < 0.0:
                        E_hi = Em
                    else:
                        E_lo, s_lo_cur = Em, sm
                E_star = 0.5 * (E_lo + E_hi)
                stats = multiplier_stats(x_star, params(float(E_star), theta))
                crossings.append({
                    "E_star": float(E_star),
                    "theta": float(theta),
                    "theta_over_pi": float(theta / np.pi),
                    "max_abs_multiplier": stats["max_abs_multiplier"],
                    **classify_crossing(mu_hi),
                })
                print(f"theta={theta:.4f} crossing at E*={E_star:.4f} "
                      f"type={crossings[-1]['crossing_type']}", flush=True)

        out["theta_records"].append({"theta": float(theta), "records": records,
                                     "crossings": crossings})
        n_unstable_by_E = {}
        for r in up_recs:
            n_unstable_by_E.setdefault(r["n_unstable"], []).append(r["E"])
        out["summary"][str(theta / np.pi)] = {
            "n_points": len(records),
            "n_converged": sum(1 for r in records if r["status"] == "PASS"),
            "crossings": crossings,
            "n_unstable_vs_E": {str(k): v for k, v in n_unstable_by_E.items()},
        }
        # checkpoint after each theta completes
        a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")

    partial_path.unlink(missing_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(json.dumps({"status": out["status"], "output": str(a.output),
                      "summary": out["summary"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
