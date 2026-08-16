#!/usr/bin/env python3
"""Gate N1/N2: model, Jacobian, units and gauge-invariant flux checks."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, finite_difference_jacobian, jacobian_phase, rhs_phase, synthetic_flux


def revision(path: Path) -> str | None:
    try:
        return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260813)
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    p = ModelParameters()
    states = [rng.normal(0.0, 0.2, 6) for _ in range(5)]
    phases = rng.uniform(0.0, 2.0*np.pi, 5)
    errors = []
    for x, phi in zip(states, phases):
        analytic = jacobian_phase(x, phi, p)
        numeric = finite_difference_jacobian(x, phi, p)
        errors.append(float(np.max(np.abs(analytic - numeric))))
    # Gauge transformation: g1 -> g1 exp(i chi1), J -> J exp(i(chi2-chi1)),
    # g2 -> g2 exp(i chi2), so g1*J*conj(g2) is invariant.
    g1, hop, g2 = 0.02*np.exp(0.3j), 0.08*np.exp(0.7j), 0.018*np.exp(-0.2j)
    base_flux = synthetic_flux(g1, hop, g2)
    chi1, chi2 = 1.1, -0.4
    transformed_flux = synthetic_flux(g1*np.exp(1j*chi1), hop*np.exp(1j*(chi2-chi1)), g2*np.exp(1j*chi2))
    trace = float(np.trace(jacobian_phase(states[0], phases[0], p)))
    expected_trace = p.divergence
    report = {
        "status": "PASS" if max(errors) < 1e-6 and abs(trace-expected_trace) < 1e-12 and abs(np.angle(np.exp(1j*(base_flux-transformed_flux)))) < 1e-12 else "FAIL",
        "gate": "N1_N2_MODEL_JACOBIAN_FLUX",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "workspace_revision": revision(Path(__file__).resolve().parents[1]),
        "seed": args.seed,
        "parameters": p.to_dict(),
        "jacobian_max_abs_errors": errors,
        "jacobian_max_abs_error": max(errors),
        "trace": trace,
        "expected_trace": expected_trace,
        "flux_base": base_flux,
        "flux_gauge_transformed": transformed_flux,
        "flux_phase_residual": float(abs(np.angle(np.exp(1j*(base_flux-transformed_flux))))),
    }
    path = out / "model_gate.json"
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "PASS" else 1

if __name__ == "__main__":
    raise SystemExit(main())
