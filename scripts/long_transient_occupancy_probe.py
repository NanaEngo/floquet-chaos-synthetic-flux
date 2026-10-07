#!/usr/bin/env python3
"""Long-transient occupation probe at a strong-coupling point.

The 2.1-period analysis-window occupations stored in
``results/strong_coupling_master_check_E*.json`` are window-precise means of a
chaotic record.  The manuscript quotes them only as *scales* (order 5e2 for the
mechanical modes at the hyperchaotic point), and the SI cites a
"longer-transient probe" for that statement.  This script is that probe: it
integrates the deterministic mean-field equations (``reconstruction_core.
rhs_phase``) from the origin with a transient ``n_attr`` drive periods (default
400, versus 40 in the production master check), then averages the mode
occupations |beta_j|^2 over a separate ``n_win``-period window (default 200).

The probe is fully deterministic (no RNG): same DOP853 tolerances and
max-step as ``strong_coupling_master_check.chaotic_trajectory``.  Window-to-
window variability is reported as the min/max of block means over ``n_blocks``
contiguous blocks, so the manuscript can state boundedness rather than a
single number.

Units: occupations are dimensionless (mean photon/phonon numbers); times are
in units of the inverse drive frequency (drive_frequency = 1, period 2*pi).

Args (CLI):
    --output        path of the JSON artifact to write.
    --drive         drive amplitude E of the checked point (default 8.0).
    --theta-over-pi synthetic-flux phase in units of pi (default 0.5).
    --n-attr        transient length in drive periods (default 400).
    --n-win         averaging window length in drive periods (default 200).
    --n-t           number of samples in the averaging window (default 20000).
    --n-blocks      number of contiguous blocks for variability (default 10).
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

from reconstruction_core import ModelParameters, rhs_phase


def make_point(drive: float, theta_over_pi: float) -> ModelParameters:
    """Build model parameters at (drive, theta) with the production values."""
    return ModelParameters(kappa=1.0, gamma1=0.02, gamma2=0.02, omega1=1.0,
                           omega2=1.03, g1=0.3, g2=0.27, hopping=0.08,
                           detuning=-1.0, drive=drive, drive_modulation=0.1,
                           drive_frequency=1.0,
                           theta=theta_over_pi * np.pi)


def long_transient_occupations(p: ModelParameters, n_attr: int, n_win: int,
                               n_t: int, n_blocks: int) -> dict:
    """Integrate a long transient, then average occupations over n_win periods.

    Args:
        p: model parameters of the checked point (dimensionless).
        n_attr: transient length in drive periods (discarded).
        n_win: averaging window length in drive periods.
        n_t: number of uniformly spaced samples stored in the window.
        n_blocks: number of contiguous blocks used for block-mean variability.

    Returns:
        dict with per-mode window means, block-mean min/max, and the final
        state of the transient (for restartability checks).
    """
    T = p.period
    t_attr = np.linspace(0.0, n_attr * T, 4000)
    sol = solve_ivp(lambda t, x: rhs_phase(x, p.drive_frequency * t, p),
                    (0.0, n_attr * T), np.zeros(6), method="DOP853",
                    rtol=1e-9, atol=1e-11, max_step=T / 200.0, t_eval=t_attr)
    if not sol.success:
        raise RuntimeError(sol.message)
    x0 = sol.y[:, -1]

    t2 = np.linspace(0.0, n_win * T, n_t)
    sol2 = solve_ivp(lambda t, x: rhs_phase(x, p.drive_frequency * t, p),
                     (0.0, n_win * T), x0, method="DOP853",
                     rtol=1e-9, atol=1e-11, max_step=T / 200.0, t_eval=t2)
    if not sol2.success:
        raise RuntimeError(sol2.message)
    y = sol2.y

    occ = {
        "optical": y[0] ** 2 + y[1] ** 2,
        "mechanical_1": y[2] ** 2 + y[3] ** 2,
        "mechanical_2": y[4] ** 2 + y[5] ** 2,
    }
    out: dict = {"window_means": {}, "block_means": {}}
    block_edges = np.linspace(0, n_t, n_blocks + 1, dtype=int)
    for name, vals in occ.items():
        block_means = [float(vals[block_edges[i]:block_edges[i + 1]].mean())
                       for i in range(n_blocks)]
        out["window_means"][name] = float(vals.mean())
        out["block_means"][name] = {
            "min": float(min(block_means)),
            "max": float(max(block_means)),
            "values": block_means,
        }
    out["transient_final_state"] = [float(v) for v in x0]
    out["all_finite"] = bool(all(np.all(np.isfinite(v)) for v in occ.values()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drive", type=float, default=8.0)
    ap.add_argument("--theta-over-pi", type=float, default=0.5)
    ap.add_argument("--n-attr", type=int, default=400)
    ap.add_argument("--n-win", type=int, default=200)
    ap.add_argument("--n-t", type=int, default=20000)
    ap.add_argument("--n-blocks", type=int, default=10)
    a = ap.parse_args()

    p = make_point(a.drive, a.theta_over_pi)
    res = long_transient_occupations(p, a.n_attr, a.n_win, a.n_t, a.n_blocks)

    # Self-check: occupations must be finite, positive and inside the
    # documented bands (a diverging or collapsed trajectory fails loudly
    # instead of silently replacing the manuscript's scale statement).
    bands = {"optical": (1.0, 100.0),
             "mechanical_1": (50.0, 5000.0),
             "mechanical_2": (50.0, 5000.0)}
    for name, mean in res["window_means"].items():
        lo, hi = bands[name]
        assert np.isfinite(mean) and mean > 0.0, f"{name} not positive finite"
        assert lo < mean < hi, f"{name}={mean} outside band ({lo}, {hi})"
    assert res["all_finite"], "trajectory not finite"

    out = {
        "gate": "LONG_TRANSIENT_OCCUPANCY_PROBE",
        "kind": "long_transient_mean_field_occupations",
        "status": "PASS",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "point": {"drive": p.drive, "theta": p.theta},
        "settings": {"n_attr_periods": a.n_attr, "n_win_periods": a.n_win,
                     "n_samples": a.n_t, "n_blocks": a.n_blocks,
                     "rtol": 1e-9, "atol": 1e-11, "max_step_periods": 200},
        "occupations": res,
        "interpretation": (
            "Deterministic mean-field occupations averaged over a 200-period "
            "window after a 400-period transient at the hyperchaotic headline "
            "point. The window means bound the mechanical occupations as "
            "scales of order 5e2; block means report window-to-window "
            "variability of the chaotic record. Semiclassical mean-field "
            "quantities only, not quantum occupations."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(json.dumps({"status": out["status"], "output": str(a.output),
                      "window_means": res["window_means"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
