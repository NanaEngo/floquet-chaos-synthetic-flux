#!/usr/bin/env python3
"""Drive-amplitude and bath-temperature sweep of the matched force-sensing gain.

Confirms whether the ~1.32x (vs flux-off) / ~1.16x (vs single-mode) enhancement
of the classical measurement Fisher information found at (drive E=0.2, n_th=0.1)
is stable across operating points. For each (drive, n_th) grid point we compute
F_C(force) at theta=pi (flux on), theta=0 (flux off, J=0.08), and J=0 (single
mode), all under matched detector noise, observation time, and estimator.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

from matched_fisher_reference import fisher_force


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drives", type=str, default="0.1,0.2,0.5,1.0")
    ap.add_argument("--n-ths", type=str, default="0.0,0.1,0.5,1.0")
    ap.add_argument("--eps", type=float, default=1e-4)
    ap.add_argument("--detector-variance", type=float, default=0.01)
    a = ap.parse_args()

    drives = [float(x) for x in a.drives.split(",")]
    n_ths = [float(x) for x in a.n_ths.split(",")]

    records = []
    for E in drives:
        for nth in n_ths:
            noise = dict(n_th1=nth, n_th2=nth,
                         detector_variance=a.detector_variance)
            row = {"drive": E, "n_th": nth}
            try:
                flux = fisher_force(np.pi, 0.08, a.eps, drive=E, **noise)
                ref_off = fisher_force(0.0, 0.08, a.eps, drive=E, **noise)
                ref_single = fisher_force(0.0, 0.0, a.eps, drive=E, **noise)
                row.update({
                    "status": "PASS" if all(r["status"] == "PASS" for r in (flux, ref_off, ref_single)) else "FAIL",
                    "F_C_flux_theta_pi": flux["classical_fisher_information"],
                    "F_C_flux_off_theta_0": ref_off["classical_fisher_information"],
                    "F_C_single_mode_J0": ref_single["classical_fisher_information"],
                    "gain_vs_flux_off": flux["classical_fisher_information"] / ref_off["classical_fisher_information"],
                    "gain_vs_single_mode": flux["classical_fisher_information"] / ref_single["classical_fisher_information"],
                    "min_physicality_eigenvalue": min(flux["min_physicality_eigenvalue"], ref_off["min_physicality_eigenvalue"], ref_single["min_physicality_eigenvalue"]),
                })
            except Exception as exc:
                row.update({"status": "FAIL", "reason": f"{type(exc).__name__}: {exc}"})
            records.append(row)
            print(f"  E={E:4} n_th={nth:4}: "
                  f"gain_vs_flux_off={row.get('gain_vs_flux_off', float('nan')):.4f} "
                  f"gain_vs_single={row.get('gain_vs_single_mode', float('nan')):.4f} "
                  f"[{row['status']}]", flush=True)

    passed = [r for r in records if r.get("status") == "PASS"]
    gains_off = [r["gain_vs_flux_off"] for r in passed]
    gains_single = [r["gain_vs_single_mode"] for r in passed]
    out = {
        "gate": "MATCHED_FISHER_DRIVE_TEMPERATURE_SWEEP",
        "kind": "classical_measurement_fisher_sweep",
        "status": "PASS" if passed and all(r.get("status") == "PASS" for r in records) else "PARTIAL" if passed else "FAIL",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "sensing_parameter": "external force F on mechanical mode 1 (momentum-quadrature drive)",
        "eps": a.eps, "detector_variance": a.detector_variance,
        "drives": drives, "n_ths": n_ths,
        "gain_vs_flux_off_stats": {
            "min": min(gains_off), "max": max(gains_off),
            "mean": float(np.mean(gains_off)), "std": float(np.std(gains_off)),
        },
        "gain_vs_single_mode_stats": {
            "min": min(gains_single), "max": max(gains_single),
            "mean": float(np.mean(gains_single)), "std": float(np.std(gains_single)),
        },
        "records": records,
        "interpretation": "Gain stability across drive amplitude and bath temperature. The ~1.3x (flux-off) enhancement is stable if gain_vs_flux_off varies only weakly across the grid.",
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nwrote {a.output}: {len(passed)}/{len(records)} PASS")
    if gains_off:
        print(f"gain_vs_flux_off: min={min(gains_off):.3f} max={max(gains_off):.3f} "
              f"mean={np.mean(gains_off):.3f} std={np.std(gains_off):.3f}")
    if gains_single:
        print(f"gain_vs_single_mode: min={min(gains_single):.3f} max={max(gains_single):.3f} "
              f"mean={np.mean(gains_single):.3f} std={np.std(gains_single):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
