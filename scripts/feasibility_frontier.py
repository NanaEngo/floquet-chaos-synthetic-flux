#!/usr/bin/env python3
"""Experimental-feasibility frontier: minimum optomechanical coupling for
hyperchaos.

The strong-coupling hyperchaos regime is placed against the Mayor/Mathew anchor
by the single-point reachability check (scripts/strong_coupling_reachability.py),
which reports g1 ~ 2542x the anchored single-photon coupling g0/omega_m = 1.18e-4.
A single point does not show whether a different drive (or a retuned dissipation)
could reach hyperchaos at a substantially smaller coupling.  This script maps the
hyperchaos onset in the (drive E, coupling g1) plane: for each drive the full
six-exponent Lyapunov spectrum is computed on a logarithmic grid of g1 (with
g2 = 0.9 g1, the ratio of the declared hyperchaos point), and the smallest g1 at
which n_+ >= 1 (chaos) and n_+ >= 2 (hyperchaos) is recorded.

The result is a *feasibility frontier*, not a device design: it answers whether
the coupling gap to the anchor is fundamental to the hyperchaos regime or an
artefact of the single declared operating point.  No pump frequency or
input-field normalization is used, so no input power is reported.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reconstruction_core import ModelParameters
from lyapunov_numba import lyapunov_qr_numba

# Mayor et al., Nat. Commun. 16, 2576 (2025): g0/omega_m = 880 kHz / 7.436 GHz.
G0_OVER_OMEGA_M = 1.18e-4


def evaluate(task: tuple) -> dict:
    (di, drive, gi, g1, p_kwargs, n_steps, transient_steps, dt, qr_interval,
     pos_threshold, seed) = task
    p = ModelParameters(drive=drive, g1=g1, g2=0.9 * g1, **p_kwargs)
    rng = np.random.default_rng(seed)
    x0 = rng.normal(0.0, 0.1, 6)
    ly = lyapunov_qr_numba(np.asarray(x0, dtype=float), p, n_steps=n_steps, dt=dt,
                           transient_steps=transient_steps, qr_interval=qr_interval)
    spectrum = np.sort(np.asarray(ly["spectrum"], dtype=float))
    n_pos = int(np.sum(spectrum > pos_threshold))
    return {
        "status": "PASS",
        "drive_index": di, "drive": drive,
        "g1_index": gi, "g1": g1, "g2": 0.9 * g1,
        "lyapunov_spectrum": spectrum.tolist(),
        "largest_exponent": float(spectrum[-1]),
        "n_positive_exponents": n_pos,
        "divergence_residual": float(abs(np.sum(spectrum)
                                         + p.kappa + p.gamma1 + p.gamma2)),
    }


def onset_coupling(records: list, order: int) -> dict:
    """Smallest g1 at which n_+ >= order, or None if never reached."""
    passing = [r for r in records if r.get("status") == "PASS"
               and r["n_positive_exponents"] >= order]
    if not passing:
        return {"reached": False, "g1": None}
    best = min(passing, key=lambda r: r["g1"])
    return {"reached": True, "g1": float(best["g1"]),
            "n_positive": best["n_positive_exponents"],
            "drive": float(best["drive"])}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drives", type=str, default="0.5,1.0,2.0,4.0,8.0")
    ap.add_argument("--g1-list", type=str,
                    default="0.02,0.0294,0.0433,0.0637,0.0937,0.138,0.203,0.3")
    ap.add_argument("--gamma1", type=float, default=0.02)
    ap.add_argument("--gamma2", type=float, default=0.02)
    ap.add_argument("--kappa", type=float, default=1.0)
    ap.add_argument("--detuning", type=float, default=-1.0)
    ap.add_argument("--hopping", type=float, default=0.08)
    ap.add_argument("--theta", type=float, default=0.0)
    ap.add_argument("--n-steps", type=int, default=1_000_000)
    ap.add_argument("--transient-steps", type=int, default=100_000)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--qr-interval", type=int, default=10)
    ap.add_argument("--positive-threshold", type=float, default=1e-3)
    ap.add_argument("--seed-base", type=int, default=20260817)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()

    p_kwargs = dict(gamma1=a.gamma1, gamma2=a.gamma2, kappa=a.kappa,
                    detuning=a.detuning, hopping=a.hopping, theta=a.theta)
    drives = [float(x) for x in a.drives.split(",")]
    g1s = [float(x) for x in a.g1_list.split(",")]

    tasks = []
    for di, drive in enumerate(drives):
        for gi, g1 in enumerate(g1s):
            seed = a.seed_base + di * len(g1s) + gi
            tasks.append((di, drive, gi, g1, p_kwargs, a.n_steps,
                          a.transient_steps, a.dt, a.qr_interval,
                          a.positive_threshold, seed))

    records = []
    with ProcessPoolExecutor(max_workers=max(1, a.workers)) as pool:
        futs = {pool.submit(evaluate, t): i for i, t in enumerate(tasks)}
        out = {}
        for fut in as_completed(futs):
            out[futs[fut]] = fut.result()
        records = [out[i] for i in range(len(tasks))]

    # Per-drive onset.
    per_drive = []
    for di, drive in enumerate(drives):
        rs = [r for r in records if r.get("drive_index") == di]
        chaos = onset_coupling(rs, 1)
        hyper = onset_coupling(rs, 2)
        per_drive.append({
            "drive": drive,
            "chaos_onset_n_ge_1": chaos,
            "hyperchaos_onset_n_ge_2": hyper,
        })

    # Global minimum coupling for hyperchaos across the whole envelope.
    all_hyper = [r for r in records if r.get("status") == "PASS"
                 and r["n_positive_exponents"] >= 2]
    global_min = min(all_hyper, key=lambda r: r["g1"]) if all_hyper else None

    global_hyper = onset_coupling([r for r in records if r.get("status") == "PASS"], 2)
    ratio = (global_hyper["g1"] / G0_OVER_OMEGA_M) if global_hyper["reached"] else None

    out = {
        "gate": "FEASIBILITY_FRONTIER",
        "kind": "minimum_coupling_for_hyperchaos",
        "status": "PASS",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "anchor_g0_over_omega_m": G0_OVER_OMEGA_M,
        "anchor_sources": ["Mayor et al., Nat. Commun. 16, 2576 (2025)"],
        "parameters": {**p_kwargs, "g2": "0.9 * g1"},
        "drives": drives,
        "g1_grid": g1s,
        "n_steps": a.n_steps,
        "transient_steps": a.transient_steps,
        "positive_exponent_threshold": a.positive_threshold,
        "per_drive_onset": per_drive,
        "global_minimum_hyperchaos_coupling": (
            {"g1": global_min["g1"], "drive": global_min["drive"],
             "n_positive": global_min["n_positive_exponents"]}
            if global_min else None),
        "global_hyperchaos_onset": global_hyper,
        "ratio_to_anchor": ratio,
        "records": records,
        "interpretation": (
            "Minimum optomechanical coupling at which n_+ >= 2 (hyperchaos) is "
            "reached, as a function of the drive, with all other coordinates at "
            "the declared hyperchaos values. If the ratio to the anchor "
            "g0/omega_m = 1.18e-4 remains of order 1e3 across the whole envelope, "
            "the coupling gap is fundamental to the hyperchaos regime rather "
            "than an artefact of the single declared operating point. No input "
            "power is reported (drive is a model coordinate)."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"wrote {a.output}")
    for pd in per_drive:
        ch = pd["chaos_onset_n_ge_1"]
        hy = pd["hyperchaos_onset_n_ge_2"]
        print(f"  E={pd['drive']:.1f}: chaos onset g1="
              f"{ch['g1'] if ch['reached'] else 'none'} | hyperchaos onset g1="
              f"{hy['g1'] if hy['reached'] else 'none'}")
    if global_hyper["reached"]:
        print(f"  global hyperchaos onset g1={global_hyper['g1']:.3f} "
              f"= {ratio:.0f}x the anchor g0/omega_m={G0_OVER_OMEGA_M:.2e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
