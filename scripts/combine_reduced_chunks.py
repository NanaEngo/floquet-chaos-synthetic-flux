#!/usr/bin/env python3
"""Combine chunked outputs from the prespecified reduced-model screen."""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("inputs", nargs="+", type=Path)
    args = parser.parse_args()

    chunks = [json.loads(path.read_text(encoding="utf-8")) for path in args.inputs]
    if not chunks:
        raise ValueError("at least one chunk is required")
    records = []
    for chunk in chunks:
        records.extend(chunk["candidate_records"])
    indices = [int(record["sample"]) for record in records]
    if len(indices) != len(set(indices)):
        raise ValueError(f"duplicate candidate indices: {indices}")
    records.sort(key=lambda record: record["sample"])
    out = {
        "status": "PROVISIONAL",
        "scientific_status": "PROVISIONAL_REDUCED_MODEL_MULTI_START_SCREEN",
        "kind": "adiabatic_bad_cavity_diagnostic_chunked",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_chunks": [str(path.resolve()) for path in args.inputs],
        "manifest": chunks[0]["manifest"],
        "design": {
            "samples": len(records),
            "global_candidate_indices": indices,
            "seed": chunks[0]["design"]["seed"],
            "robust_replicates": chunks[0]["design"]["robust_replicates"],
            "initial_condition_replicates": chunks[0]["design"]["initial_condition_replicates"],
            "initial_condition_scale": chunks[0]["design"]["initial_condition_scale"],
            "orbit_max_iter": chunks[0]["design"]["orbit_max_iter"],
            "orbit_residual_tolerance": chunks[0]["design"]["orbit_residual_tolerance"],
            "initial_condition_protocol": "independent root-solver starts; finite-time basin capture is NOT_COMPUTED",
            "parameter_replication_protocol": "independent uniform perturbation draws; exact replicate parameters stored per candidate",
            "sampling": "independent uniform stream, executed as deterministic chunks",
        },
        "candidate_records": records,
        "candidate_count": len(records),
        "feasible_count": sum(bool(record["feasible"]) for record in records),
        "pareto_count": 0,
        "interpretation": "Reduced bad-cavity multi-start fixed-point diagnostic only; all residual/radius gates passed for the recorded starts, but finite-time basin capture, SI calibration, and manuscript claims are not authorized.",
        "limitations": [
            "optical mode adiabatically eliminated",
            "drive remains model coordinate",
            "J_m and g_i are assumed proxies",
            "full/reduced agreement tested only at kappa=20 and 100",
            "comparison does not validate the assumed closure or an experimental device",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: out[key] for key in ("candidate_count", "feasible_count", "pareto_count")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
