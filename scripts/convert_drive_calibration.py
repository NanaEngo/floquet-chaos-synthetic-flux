#!/usr/bin/env python3
"""Fail-closed conversion of coherent input amplitude to model drive.

The input convention is deliberately explicit:
  alpha_in_magnitude_sqrt_per_s = |alpha_in| in sqrt(photons / second)
  omega_p_rad_s                = laser angular frequency
  kappa_rad_s                  = cavity decay rate
  omega_ref_rad_s              = normalization angular frequency
  alpha_in_phase_rad           = phase relative to the model's real drive axis

Under the standard input-output convention:
    P_in = hbar * omega_p * |alpha_in|**2

After nondimensionalizing time by omega_ref:
    alpha_in_norm = alpha_in / sqrt(omega_ref)
    kappa_norm = kappa / omega_ref
    E_model = sqrt(kappa_norm) * alpha_in_norm

The current reconstruction accepts a real additive drive only.  Therefore a
non-zero quadrature phase is rejected rather than silently projected away.
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

HBAR = 1.054571817e-34  # J s, exact SI-defined value
REQUIRED = (
    "alpha_in_magnitude_sqrt_per_s",
    "alpha_in_phase_rad",
    "omega_p_rad_s",
    "kappa_rad_s",
    "omega_ref_rad_s",
)


def _finite_positive(name: str, value: object) -> float | None:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return None
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        return None
    return value


def convert(payload: dict) -> dict:
    errors: list[str] = []
    if payload.get("alpha_in_definition") != "coherent_input_amplitude_sqrt_photons_per_second":
        errors.append(
            "alpha_in_definition must be 'coherent_input_amplitude_sqrt_photons_per_second'"
        )

    values: dict[str, float] = {}
    for name in REQUIRED:
        raw = payload.get(name)
        if name == "alpha_in_phase_rad":
            if not isinstance(raw, (int, float)) or isinstance(raw, bool):
                errors.append(f"{name} is missing or not numeric")
            else:
                value = float(raw)
                if not math.isfinite(value):
                    errors.append(f"{name} is not finite")
                else:
                    values[name] = value
        else:
            value = _finite_positive(name, raw)
            if value is None:
                errors.append(f"{name} is missing, non-positive, or not finite")
            else:
                values[name] = value

    if errors:
        return {
            "status": "NOT_COMPUTED",
            "scientific_status": "INCOMPLETE_DRIVE_CALIBRATION",
            "errors": errors,
            "interpretation": "No power or model drive was calculated; missing conventions or values must be supplied explicitly.",
        }

    phase = values["alpha_in_phase_rad"]
    # The reconstruction's drive coordinate is real.  Do not discard a quadrature
    # component silently; the author must either supply an in-phase amplitude or
    # extend the model before using a complex input field.
    phase_mod = math.atan2(math.sin(phase), math.cos(phase))
    if not math.isclose(math.sin(phase_mod), 0.0, abs_tol=1.0e-12):
        return {
            "status": "NOT_COMPUTED",
            "scientific_status": "COMPLEX_DRIVE_PHASE_UNSUPPORTED_BY_REAL_MODEL",
            "errors": [
                "alpha_in_phase_rad is not 0 or pi; the current model accepts only a real in-phase drive",
            ],
            "phase_rad": phase,
        }

    alpha = values["alpha_in_magnitude_sqrt_per_s"]
    omega_p = values["omega_p_rad_s"]
    kappa = values["kappa_rad_s"]
    omega_ref = values["omega_ref_rad_s"]
    kappa_norm = kappa / omega_ref
    alpha_norm = alpha / math.sqrt(omega_ref)
    drive_norm = math.sqrt(kappa_norm) * alpha_norm
    if math.cos(phase_mod) < 0.0:
        drive_norm = -drive_norm
    power_w = HBAR * omega_p * alpha * alpha

    return {
        "status": "PASS",
        "scientific_status": "CONDITIONAL_STANDARD_INPUT_OUTPUT_CONVERSION",
        "convention": {
            "alpha_in_definition": payload["alpha_in_definition"],
            "power_relation": "P_in = hbar * omega_p * |alpha_in|^2",
            "model_relation": "E_model = sqrt(kappa_norm) * alpha_in_norm",
            "real_drive_phase_requirement": "alpha_in_phase_rad must be 0 or pi",
        },
        "inputs": {
            **values,
            "hbar_j_s": HBAR,
        },
        "derived": {
            "kappa_normalized": kappa_norm,
            "alpha_in_normalized": alpha_norm,
            "input_power_w": power_w,
            "model_drive_normalized": drive_norm,
        },
        "provenance": {
            "source": "arXiv:2502.12336v2 Sections 2.2-2.3.3; standard input-output power convention",
            "source_url": "https://arxiv.org/html/2502.12336",
            "calibration_scope": "Conversion only; not evidence that the setup is experimentally calibrated or optimal.",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        payload = json.loads(args.input.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        result = {
            "status": "NOT_COMPUTED",
            "scientific_status": "INVALID_CALIBRATION_INPUT",
            "errors": [str(exc)],
        }
        rc = 2
    else:
        result = convert(payload)
        result["recorded_at_utc"] = datetime.now(timezone.utc).isoformat()
        result["python"] = sys.version
        result["platform"] = platform.platform()
        rc = 0 if result["status"] == "PASS" else 2
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
