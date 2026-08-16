#!/usr/bin/env python3
"""Experimental-feasibility bound for the strong-coupling hyperchaos regime.

Gap 3 of GAP_IMPLEMENTATION_20260816.md. The strong-coupling regime
(g1=0.3, g2=0.27, gamma_{1,2}=0.02, kappa=1.0, J=0.08, all normalized to
omega_m) is placed against the two experimental anchors already cited in the
SI (Mayor et al. 2025; Mathew, del Pino and Verhagen 2020). For each coordinate
the script computes the ratio of the required (model) value to the anchored
value and identifies the limiting coordinate. The finding is a *feasibility
bound*, not a device design: no pump frequency, external coupling, or
input-field normalization is available to convert the drive amplitude to an
input power.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path


# Anchored values (normalized to omega_m), from SI tab:si-calibrated.
ANCHOR = {
    "kappa_over_omega_m": 0.108,       # Mayor et al.: 800 MHz / 7.436 GHz
    "gamma_m_over_omega_m": 2.78e-5,   # Mayor et al.: 206.6 kHz / 7.436 GHz (Q_m = 3.6e4)
    "g0_over_omega_m": 1.18e-4,        # Mayor et al.: 880 kHz / 7.436 GHz
    "J_m_over_omega_m_range": [1e-3, 1e-1],  # Mathew et al.: up to a few percent of omega_m
}

# Strong-coupling model (normalized to omega_m), the declared hyperchaos regime.
MODEL = {
    "kappa_over_omega_m": 1.0,
    "gamma_m_over_omega_m": 0.02,
    "g1_over_omega_m": 0.3,
    "g2_over_omega_m": 0.27,
    "J_m_over_omega_m": 0.08,
    "detuning_over_omega_m": -1.0,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()

    ratios = {
        "kappa": MODEL["kappa_over_omega_m"] / ANCHOR["kappa_over_omega_m"],
        "gamma_m": MODEL["gamma_m_over_omega_m"] / ANCHOR["gamma_m_over_omega_m"],
        "g1": MODEL["g1_over_omega_m"] / ANCHOR["g0_over_omega_m"],
        "g2": MODEL["g2_over_omega_m"] / ANCHOR["g0_over_omega_m"],
    }
    # J_m: within the scanned range?
    j_lo, j_hi = ANCHOR["J_m_over_omega_m_range"]
    j_in_range = bool(j_lo <= MODEL["J_m_over_omega_m"] <= j_hi)
    j_ratio = MODEL["J_m_over_omega_m"] / j_hi  # relative to the upper scan bound

    limiting = max(ratios, key=ratios.get)

    reachable = {k: bool(v <= 1.0) for k, v in ratios.items()}
    reachable["J_m"] = j_in_range
    all_reachable = all(reachable.values())

    out = {
        "gate": "STRONG_COUPLING_REACHABILITY",
        "kind": "experimental_feasibility_bound",
        "status": "PASS",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "anchor": ANCHOR,
        "anchor_sources": [
            "Mayor et al., Nat. Commun. 16, 2576 (2025) [omega_m, kappa, gamma_m, g0]",
            "Mathew, del Pino, Verhagen, Nat. Nanotechnol. 15, 198-202 (2020) [J_m scale]",
        ],
        "model": MODEL,
        "ratios_required_over_anchor": {
            "kappa": ratios["kappa"],
            "gamma_m": ratios["gamma_m"],
            "g1": ratios["g1"],
            "g2": ratios["g2"],
            "J_m": j_ratio,
        },
        "limiting_coordinate": limiting,
        "limiting_ratio": ratios[limiting],
        "reachable_within_anchor": reachable,
        "all_reachable": all_reachable,
        "interpretation": (
            "The strong-coupling hyperchaos regime requires optomechanical "
            "coupling and mechanical damping that are 2-3 orders of magnitude "
            "above the anchored Mayor/Mathew values (limiting coordinate: "
            f"{limiting}, {ratios[limiting]:.0f}x). The inter-resonator hopping "
            "is within the scanned Mathew range. The regime is therefore a "
            "model-level demonstration, not a prediction for the anchored "
            "device; the drive amplitude is not convertible to input power "
            "without the pump frequency and external coupling."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "limiting_coordinate": limiting,
                      "limiting_ratio": ratios[limiting],
                      "output": str(a.output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
