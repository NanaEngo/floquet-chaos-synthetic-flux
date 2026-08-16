#!/usr/bin/env python3
"""Common-regime consistency test for the bad-cavity reduction."""
from __future__ import annotations

import json
import platform
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, find_periodic_orbit
from reduced_bad_cavity_screen import find_periodic_orbit_root


def run_case(kappa: float) -> dict:
    p = replace(ModelParameters(), kappa=kappa, drive=0.2, drive_modulation=0.1,
                drive_frequency=1.0)
    full = find_periodic_orbit(np.zeros(6), p, max_iter=300, residual_tol=1e-8)
    reduced = find_periodic_orbit_root(np.zeros(4), p, max_iter=60, residual_tol=1e-8)
    record = {"kappa": kappa, "full_status": full["status"],
              "reduced_status": reduced["status"]}
    if full["status"] == "PASS" and reduced["status"] == "PASS":
        full_mechanical = np.asarray(full["x0"])[2:]
        reduced_mechanical = np.asarray(reduced["y0"])
        record["mechanical_inf_error"] = float(np.linalg.norm(full_mechanical - reduced_mechanical, ord=np.inf))
        record["full_residual"] = full["residual"]
        record["reduced_residual"] = reduced["residual"]
    else:
        record["reason"] = "one of the two periodic-orbit solves did not pass"
    return record


def main() -> int:
    cases = [run_case(20.0), run_case(100.0)]
    out = {
        "status": "PASS" if all(c["full_status"] == "PASS" and c["reduced_status"] == "PASS" for c in cases) else "INCOMPLETE",
        "kind": "full_reduced_common_regime_comparison",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "cases": cases,
        "interpretation": "The comparison checks implementation consistency only; it does not validate the experimental parameter mapping or the assumed Mathew closure."
    }
    path = Path(__file__).resolve().parents[1] / "results" / "reduced_full_common_regime_comparison.json"
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
