#!/usr/bin/env python3
"""Pareto analysis for the reduced-model diagnostic.

Only two measured objectives are available in the reduced screen:
maximize robust stability margin and minimize normalized drive cost. Fisher
information is intentionally not imputed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def dominates(a: dict, b: dict) -> bool:
    return (
        a["robust_stability_margin"] >= b["robust_stability_margin"]
        and a["drive_cost"] <= b["drive_cost"]
        and (
            a["robust_stability_margin"] > b["robust_stability_margin"]
            or a["drive_cost"] < b["drive_cost"]
        )
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.input.read_text(encoding="utf-8"))
    records = []
    for item in source["candidate_records"]:
        radii = [r["poincare_spectral_radius"] for r in item["replicate_records"] if r["status"] == "PASS"]
        if not item["feasible"] or not radii:
            continue
        record = {
            "sample": item["sample"],
            "parameters": item["parameters"],
            "robust_stability_margin": 1.0 - max(radii),
            "minimum_stability_margin": 1.0 - max(radii),
            "drive_cost": abs(item["parameters"]["drive"]),
            "maximum_poincare_spectral_radius": max(radii),
            "minimum_poincare_spectral_radius": min(radii),
        }
        records.append(record)
    pareto = [r for r in records if not any(dominates(other, r) for other in records)]
    out = {
        "status": "PROVISIONAL",
        "scientific_status": "PROVISIONAL_REDUCED_MODEL_MULTI_START_PARETO",
        "kind": "two_objective_reduced_model_pareto",
        "source": str(args.input.resolve()),
        "objectives": {
            "maximize": "robust_stability_margin",
            "minimize": "drive_cost",
            "fisher_information": "NOT_COMPUTED"
        },
        "candidate_count": len(records),
        "pareto_count": len(pareto),
        "candidates": records,
        "pareto_candidates": pareto,
        "interpretation": "Provisional reduced-model multi-start front only. The gate concerns fixed-point residuals and local Poincare radii from independent root starts; finite-time basin capture, SI calibration, experimental optimality, and Fisher information are not established."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("scientific_status", "candidate_count", "pareto_count")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
