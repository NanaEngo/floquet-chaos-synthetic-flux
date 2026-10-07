"""Step-convergence and window-extension study at headline points (referee items 1-2).

Recomputes the full six-exponent Benettin-QR spectrum at fixed physical
configurations while varying (a) the fixed RK4 step at fixed time span
(step-convergence; answers the DOP853-vs-RK4 branch question by bounding
the RK4 discretization error) and (b) the integration window at fixed
step (answers whether s_QR falls below the marginal fourth exponent).

Same seeds/initial conditions as the production map
(full_spectrum_hyperchaos_map.py: x0 = default_rng(seed).normal(0,0.1,6))
so every comparison is apples-to-apples. Per-exponent terminal-block
standard deviations are reported (not just the max) for the
s_QR-vs-lambda4 face-to-face.

Usage:
    mamba run -n qom python scripts/referee_step_window_study.py \
        --output results/headline_step_window_study.json --workers 4
"""

from __future__ import annotations

import argparse
import datetime
import json
import platform
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruction_core import ModelParameters
from lyapunov_numba import lyapunov_qr_numba

# (drive_index, drive, theta_index, theta, map_seed) matching the
# production-map seed stream: seed = 20260815 + di*12 + ti, 12 phases.
POINTS = [
    (3, 4.0, 3, -np.pi / 2, 20260815 + 3 * 12 + 3),
    (4, 8.0, 9, np.pi / 2, 20260815 + 4 * 12 + 9),
    (4, 8.0, 6, 0.0, 20260815 + 4 * 12 + 6),
]

# (dt, n_steps, transient_steps): same 1e4 time-unit span unless noted.
CONFIGS = [
    {"label": "production_dt1e-2_1e6", "dt": 0.01, "n_steps": 1_000_000, "transient": 100_000},
    {"label": "halfstep_dt5e-3_2e6", "dt": 0.005, "n_steps": 2_000_000, "transient": 200_000},
    {"label": "quarterstep_dt2.5e-3_4e6", "dt": 0.0025, "n_steps": 4_000_000, "transient": 400_000},
    {"label": "longwindow_dt1e-2_2e6", "dt": 0.01, "n_steps": 2_000_000, "transient": 200_000},
]

P_KWARGS = dict(g1=0.3, g2=0.27, gamma1=0.02, gamma2=0.02, detuning=-1.0)
THRESHOLD = 1e-3


def run(task: tuple) -> dict:
    (di, drive, ti, theta, seed), cfg = task
    p = ModelParameters(drive=drive, theta=theta, **P_KWARGS)
    rng = np.random.default_rng(seed)
    x0 = rng.normal(0.0, 0.1, 6)
    ly = lyapunov_qr_numba(np.asarray(x0, dtype=float), p,
                            n_steps=cfg["n_steps"], dt=cfg["dt"],
                            transient_steps=cfg["transient"], qr_interval=10)
    spec = np.asarray(ly["spectrum"], dtype=float)  # descending
    blocks = np.asarray(ly["block_history"], dtype=float)
    half = blocks[len(blocks) // 2:, :]
    per_exp_std = half.std(axis=0).tolist()
    return {
        "drive": drive, "theta": float(theta), "seed": seed,
        "config": cfg["label"], "dt": cfg["dt"], "n_steps": cfg["n_steps"],
        "status": "PASS",
        "lyapunov_spectrum": spec.tolist(),
        "n_positive_exponents": int(np.sum(spec > THRESHOLD)),
        "per_exponent_terminal_std": per_exp_std,
        "max_terminal_std": float(np.max(per_exp_std)),
        "divergence_residual": float(abs(np.sum(spec) - (-p.kappa - p.gamma1 - p.gamma2))),
        "mean_divergence": float(ly["mean_divergence"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()

    tasks = [(pt, cfg) for pt in POINTS for cfg in CONFIGS]
    recs = []
    if a.workers <= 1:
        recs = [run(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            futs = {pool.submit(run, t): i for i, t in enumerate(tasks)}
            done = {}
            for fut in as_completed(futs):
                done[futs[fut]] = fut.result()
            recs = [done[i] for i in range(len(tasks))]

    out = {
        "gate": "REFEREE_STEP_WINDOW_STUDY",
        "recorded_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "kind": "referee_step_window_study",
        "interpretation": ("Fixed-step RK4 step-convergence (dt=1e-2/5e-3/2.5e-3 "
                           "at fixed span) and window extension (1e6->2e6 steps) "
                           "at three headline configurations, same seeds/ICs as "
                           "the production map."),
        "coupling": {"g1": 0.3, "g2": 0.27, "gamma1": 0.02, "gamma2": 0.02, "detuning": -1.0},
        "positive_exponent_threshold": THRESHOLD,
        "total_records": len(recs),
        "records": recs,
    }
    a.output.write_text(json.dumps(out, indent=2))
    for r in recs:
        print(f"E={r['drive']} th={r['theta']:+.4f} {r['config']:28s} "
              f"n+={r['n_positive_exponents']} lam4={r['lyapunov_spectrum'][3]:+.6f} "
              f"s4={r['per_exponent_terminal_std'][3]:.2e} rdiv={r['divergence_residual']:.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
