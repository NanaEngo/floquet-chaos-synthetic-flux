#!/usr/bin/env python3
"""Figure: matched measurement-Fisher reference for force sensing.

Plots F_C(force) versus synthetic-flux phase with the two matched references
(flux-off and single-mode) overlaid, so the operational sensing question
(advantage or no advantage) is visible at a glance.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path,
                    default=ROOT / "results" / "matched_fisher_reference.json")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "figures" / "matched_fisher_reference")
    a = ap.parse_args()

    d = json.loads(a.input.read_text(encoding="utf-8"))
    records = [r for r in d["records"] if r.get("status") == "PASS"]
    refA = d["references"]["flux_off_theta_0_hopping_0.08"]
    refB = d["references"]["single_mode_hopping_0"]
    gains = {g["theta"]: g for g in d["gains"]}

    thetas = np.array([r["theta"] for r in records]) / np.pi
    fc = np.array([r["classical_fisher_information"] for r in records])
    gA = np.array([gains[r["theta"]]["gain_vs_flux_off"] for r in records])
    gB = np.array([gains[r["theta"]]["gain_vs_single_mode"] for r in records])

    fig, (ax, axg) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    ax.plot(thetas, fc, "-o", ms=5, lw=1.6, color="#1f77b4", label="synthetic-flux sensor")
    ax.axhline(refA["classical_fisher_information"], ls="--", lw=1.2,
               color="#2ca02c", label="reference A: flux off ($\\theta=0$)")
    ax.axhline(refB["classical_fisher_information"], ls=":", lw=1.4,
               color="#d62728", label="reference B: single mode ($J=0$)")
    ax.set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    ax.set_ylabel(r"measurement Fisher $F_C(F)$")
    ax.text(0.02, 0.98, "(a)", transform=ax.transAxes,
            va="top", ha="left", fontsize=11)
    ax.legend(fontsize=8, loc="best")

    axg.plot(thetas, gA, "-o", ms=5, lw=1.6, color="#2ca02c", label="vs flux off")
    axg.plot(thetas, gB, "-s", ms=5, lw=1.6, color="#d62728", label="vs single mode")
    axg.axhline(1.0, ls="--", lw=1.0, color="0.5")
    axg.set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    axg.set_ylabel("gain $F_C/F_C^{\\rm ref}$")
    axg.text(0.02, 0.98, "(b)", transform=axg.transAxes,
             va="top", ha="left", fontsize=11)
    axg.legend(fontsize=8, loc="best")
    axg.set_ylim(0.8, 1.3)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{a.output}.{ext}", bbox_inches="tight", dpi=200)
    plt.close(fig)

    manifest = {
        "gate": "MATCHED_FISHER_REFERENCE_FIGURE",
        "input": str(a.input),
        "outputs": [f"{a.output}.pdf", f"{a.output}.png"],
        "note": "(a) F_C(force) vs synthetic-flux phase with flux-off and single-mode references overlaid; (b) gain relative to each reference. A gain <= 1 indicates no flux-mediated sensing advantage.",
    }
    (ROOT / "figures" / "matched_fisher_reference.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
