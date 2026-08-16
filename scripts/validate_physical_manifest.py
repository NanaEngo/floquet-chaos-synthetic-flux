#!/usr/bin/env python3
"""Validate physical calibration completeness without inventing missing values."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    data = json.loads(args.manifest.read_text(encoding="utf-8"))
    errors = []
    if data.get("status") != "PASS":
        errors.append(f"manifest status is {data.get('status')!r}, expected 'PASS'")
    normalization = data.get("normalization", {})
    if normalization.get("reference_frequency_si_hz") is None:
        errors.append("reference_frequency_si_hz is missing")
    if normalization.get("reference_frequency_source") is None:
        errors.append("reference_frequency_source is missing")
    objective = data.get("objective", {})
    if objective.get("selected_primary_objective") is None:
        errors.append("selected_primary_objective is missing")
    if not objective.get("constraints"):
        errors.append("objective constraints are missing")
    for parameter in data.get("parameters", []):
        name = parameter.get("name", "<unnamed>")
        for field in ("si_value", "si_units", "admissible_si_range", "experimental_source", "uncertainty"):
            if parameter.get(field) is None:
                errors.append(f"{name}.{field} is missing")
    if errors:
        print("NOT_COMPUTED: physical calibration manifest is incomplete")
        for error in errors:
            print(f"- {error}")
        return 1
    print("PASS: physical calibration manifest is complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
