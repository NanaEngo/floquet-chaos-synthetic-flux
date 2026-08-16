#!/usr/bin/env python3
"""Force-step (eps) convergence of the matched measurement-Fisher reference.

Stores the central-difference force-step dependence of F_C(F) at the two flux
extrema (theta = pi/2 and theta = pi) for eps in {1e-3, 1e-4, 1e-5}, under the
same matched resources as `matched_fisher_reference.py` (drive 0.2, one drive
period, detector variance 0.01, baths n_th = 0.1).

This artifact exists because the canonical stored run of
`matched_fisher_reference.py` was executed with an empty `--eps-check`, so the
eps-convergence rows of SI Table (tab:si-fisher-eps) previously had no stored
provenance. The values are deterministic (no RNG), so they reproduce exactly
the table entries computed on 2026-08-15.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from matched_fisher_reference import fisher_force


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--eps-list", nargs="+", type=float, default=[1e-3, 1e-4, 1e-5])
    ap.add_argument("--thetas", nargs="+", type=float,
                    default=[np.pi / 2, np.pi])
    ap.add_argument("--hopping", type=float, default=0.08)
    ap.add_argument("--n-th1", type=float, default=0.1)
    ap.add_argument("--n-th2", type=float, default=0.1)
    ap.add_argument("--detector-variance", type=float, default=0.01)
    ap.add_argument("--drive", type=float, default=0.2)
    a = ap.parse_args()

    noise = dict(n_th1=a.n_th1, n_th2=a.n_th2, detector_variance=a.detector_variance)
    records = []
    ok = True
    for th in a.thetas:
        for eps in a.eps_list:
            try:
                r = fisher_force(float(th), a.hopping, eps, drive=a.drive, **noise)
                records.append({"theta": float(th), "theta_over_pi": float(th / np.pi),
                                "force_eps": eps,
                                "classical_fisher_information": r["classical_fisher_information"],
                                "d_mean_dF": r["d_mean_dF"], "d_variance_dF": r["d_variance_dF"],
                                "min_physicality_eigenvalue": r["min_physicality_eigenvalue"],
                                "status": r["status"]})
                ok = ok and r["status"] == "PASS"
            except Exception as exc:
                records.append({"theta": float(th), "theta_over_pi": float(th / np.pi),
                                "force_eps": eps, "status": "FAIL",
                                "reason": f"{type(exc).__name__}: {exc}"})
                ok = False

    # Relative spread per theta (max-min)/mean over the three eps values.
    spread = {}
    for th in a.thetas:
        sub = [r["classical_fisher_information"] for r in records
               if abs(r["theta"] - th) < 1e-12 and r["status"] == "PASS"]
        if len(sub) == 3:
            mean = float(np.mean(sub))
            spread[float(th / np.pi)] = {
                "max_minus_min": float(max(sub) - min(sub)),
                "relative_spread": float((max(sub) - min(sub)) / mean),
                "converged_below_1e-4_relative": float((max(sub) - min(sub)) / mean) < 1e-4,
            }

    out = {
        "gate": "MATCHED_FISHER_EPS_CONVERGENCE",
        "kind": "classical_measurement_fisher_force_sensing_eps_convergence",
        "status": "PASS" if ok else "FAIL",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "force_eps_list": a.eps_list,
        "thetas": a.thetas,
        "matched_resources": {
            "drive": a.drive, "observation_time": "1 drive period",
            "detector_noise_variance": a.detector_variance,
            "baths": [a.n_th1, a.n_th2],
            "estimator": "eq:fisher (mean + variance terms)",
        },
        "records": records,
        "relative_spread": spread,
        "interpretation": ("Central-difference force-step convergence of the matched "
                           "measurement Fisher information at the two flux extrema. "
                           "Not QFI and not a claim of quantum advantage."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": out["status"], "output": str(a.output),
                      "records": len(records), "spread": spread}, indent=2))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
