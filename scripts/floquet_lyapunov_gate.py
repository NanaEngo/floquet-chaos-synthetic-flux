#!/usr/bin/env python3
"""Gate N3: periodic orbit, monodromy and Lyapunov consistency."""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, find_periodic_orbit, lyapunov_qr, monodromy


def git_revision(path: Path) -> str | None:
    try:
        return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    # The weakly damped periodic orbit needs a long production window for
    # individual QR exponents to resolve nearly degenerate Floquet rates.
    parser.add_argument("--n-steps", type=int, default=1_000_000)
    parser.add_argument("--transient-steps", type=int, default=100_000)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--qr-interval", type=int, default=10)
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    p = ModelParameters()
    x0 = np.zeros(6, dtype=float)
    # The weakly damped mechanical modes require more than 100 Poincare
    # iterations to reach the declared 1e-8 fixed-point residual.
    orbit = find_periodic_orbit(x0, p, max_iter=600)
    report = {
        "gate": "N3_FLOQUET_LYAPUNOV",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "workspace_revision": git_revision(Path(__file__).resolve().parents[1]),
        "parameters": p.to_dict(),
        "periodic_orbit": orbit,
        "status": "FAIL",
    }
    if orbit.get("status") != "PASS":
        report["reason"] = orbit.get("reason", "periodic orbit not converged")
        (out / "floquet_lyapunov.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2, sort_keys=True))
        return 1
    x_periodic = np.asarray(orbit["x0"], dtype=float)
    floquet = monodromy(x_periodic, p)
    lyap = lyapunov_qr(x_periodic, p, n_steps=args.n_steps, dt=args.dt,
                       transient_steps=args.transient_steps, qr_interval=args.qr_interval)
    rates = np.sort(np.real(np.asarray(floquet["floquet_rates"])))
    spectrum = np.sort(np.asarray(lyap["spectrum"]))
    comparison = np.abs(rates - spectrum)
    divergence_residual = abs(float(np.sum(spectrum)) - float(lyap["mean_divergence"]))
    block = np.asarray(lyap["block_history"], dtype=float)
    block_tail = block[max(0, len(block)-10):]
    block_uncertainty = np.max(np.std(block_tail, axis=0)) if len(block_tail) > 1 else float("inf")
    report.update({
        "floquet": floquet,
        "lyapunov": lyap,
        "comparison": {
            "floquet_rates_sorted": rates.tolist(),
            "lyapunov_spectrum_sorted": spectrum.tolist(),
            "max_rate_spectrum_abs_difference": float(np.max(comparison)),
            "lyapunov_sum": float(np.sum(spectrum)),
            "mean_divergence": float(lyap["mean_divergence"]),
            "divergence_residual": divergence_residual,
            "tail_block_max_std": float(block_uncertainty),
        },
    })
    # This pilot is a valid gate only if the periodic orbit and independent
    # diagnostics agree. Thresholds are explicit and may be tightened later.
    passed = (orbit["residual"] < 1e-8 and np.max(comparison) < 5e-4
              and divergence_residual < 5e-3 and block_uncertainty < 5e-3)
    report["status"] = "PASS" if passed else "FAIL"
    if not passed:
        report["reason"] = "Floquet/Lyapunov/divergence convergence gate failed"
    (out / "floquet_lyapunov.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if passed else 1

if __name__ == "__main__":
    raise SystemExit(main())
