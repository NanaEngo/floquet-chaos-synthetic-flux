#!/usr/bin/env python3
"""Convergence check for the vectorized reduced-equation RK4 map."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters
from reduced_bad_cavity_screen import integrate_period
from fast_basin_capture import parameter_from_dict, rk4_period_batch


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--screen", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, nargs="+", default=[5, 10, 15])
    parser.add_argument("--replicates", type=int, default=5)
    parser.add_argument("--states", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260816)
    parser.add_argument("--substeps", type=int, nargs="+", default=[20, 40, 80])
    args = parser.parse_args()
    source = json.loads(args.screen.read_text(encoding="utf-8"))
    candidates = {int(item["sample"]): item for item in source["candidate_records"]}
    rng = np.random.default_rng(args.seed)
    records = []
    for sample in sorted(set(args.samples)):
        item = candidates[sample]
        for rep in range(args.replicates):
            p = parameter_from_dict(item["replicate_parameters"][str(rep)])
            states = rng.normal(0.0, 0.1, size=(args.states, 4))
            reference = np.vstack([integrate_period(state, p).y[:, -1] for state in states])
            errors = {}
            for substeps in args.substeps:
                approx = rk4_period_batch(states, p, substeps)
                errors[str(substeps)] = float(np.max(np.abs(approx - reference)))
            records.append({"sample": sample, "replicate": rep, "errors_inf_norm": errors})
    summary = {str(s): max(r["errors_inf_norm"][str(s)] for r in records) for s in args.substeps}
    out = {
        "status": "PASS" if summary[str(max(args.substeps))] < 1e-7 else "INCOMPLETE",
        "scientific_status": "REDUCED_RK4_DOP853_CONVERGENCE",
        "kind": "one_period_integrator_comparison",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_screen": str(args.screen.resolve()),
        "design": {"samples": sorted(set(args.samples)), "replicates": args.replicates,
                   "states_per_case": args.states, "seed": args.seed, "substeps": args.substeps,
                   "reference": "DOP853 rtol=1e-9 atol=1e-11 max_step=T/400"},
        "maximum_error_by_substeps": summary,
        "records": records,
        "interpretation": "This validates the one-period map approximation used by the vectorized basin screen; it does not validate SI calibration or global basin structure.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": out["status"], "maximum_error_by_substeps": summary}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
