#!/usr/bin/env python3
"""Figure: drive-axis bifurcation structure at strong optomechanical coupling.

Panels:
(a) maximum Floquet multiplier |mu|_max of the drive-locked orbit versus drive
    amplitude E for theta = 0 and theta = pi/2 (upward and downward
    continuation); the dashed line marks the unit circle and markers locate the
    refined first bifurcation points E*.
(b) largest Floquet rate of the orbit versus E, with the largest Lyapunov
    exponent of the stored strong-coupling map at the five canonical drives
    overlaid (per-phase records for theta = 0 and theta = pi/2) as a
    cross-validation of the two independent diagnostics.
(c) number of unstable multipliers (|mu| > 1) versus E.
(d) mean cavity amplitude <|alpha|> over the orbit period versus E.

Provenance-bound: reads only the two immutable result JSONs listed in the
output manifest.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

try:
    import scienceplots  # noqa: F401
    plt.style.use(["science", "no-latex"])
except Exception:
    plt.style.use("seaborn-v0_8-whitegrid")

ROOT = Path(__file__).resolve().parent.parent
SCAN = ROOT / "results" / "drive_axis_bifurcation_scan.json"
LYAP = ROOT / "results" / "full_spectrum_hyperchaos_map_g03.json"


def main() -> int:
    scan = json.loads(SCAN.read_text())
    lyap = json.loads(LYAP.read_text())

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    colors = {0.0: "#1f77b4", np.pi / 2: "#d62728"}
    labels = {0.0: r"$\theta=0$", np.pi / 2: r"$\theta=\pi/2$"}

    for theta_rec in scan["theta_records"]:
        th = theta_rec["theta"]
        col = colors[th]
        up = [r for r in theta_rec["records"] if r["direction"] == "up"]
        dn = [r for r in theta_rec["records"] if r["direction"] == "down"]
        for recs, ls in ((up, "-"), (dn, "--")):
            Es = [r["E"] for r in recs if r["status"] == "PASS"]
            mus = [r["max_abs_multiplier"] for r in recs if r["status"] == "PASS"]
            axes[0, 0].plot(Es, mus, ls, color=col, lw=1.1,
                            label=labels[th] if ls == "-" else None)
            nus = [r["n_unstable"] for r in recs if r["status"] == "PASS"]
            axes[1, 0].plot(Es, nus, ls, color=col, lw=1.1)
            amps = [r["mean_amplitude"] for r in recs if r["status"] == "PASS"]
            axes[1, 1].plot(Es, amps, ls, color=col, lw=1.1)
        for cr in theta_rec["crossings"]:
            axes[0, 0].axvline(cr["E_star"], color=col, lw=0.8, ls=":", alpha=0.7)
            axes[0, 0].annotate(f"$E^*={cr['E_star']:.2f}$", xy=(cr["E_star"], 1.0),
                                xytext=(cr["E_star"] + 0.15, 1.06),
                                fontsize=8, color=col)

    axes[0, 0].axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.6)
    axes[0, 0].set_ylabel(r"max $|\mu|$")
    axes[0, 0].set_xlabel(r"drive $E$")
    axes[0, 0].legend(fontsize=8, frameon=True)
    axes[0, 0].set_ylim(0.9, None)

    # (b) largest Floquet rate + Lyapunov cross-validation at canonical drives
    for th in (0.0, np.pi / 2):
        col = colors[th]
        up = [r for r in scan["theta_records"] if r["theta"] == th][0]
        up = [r for r in up["records"] if r["direction"] == "up" and r["status"] == "PASS"]
        axes[0, 1].plot([r["E"] for r in up], [r["floquet_rates"][0] for r in up],
                        color=col, lw=1.1, label=labels[th] + " (orbit)")
        ly_pts = [r for r in lyap["records"]
                  if abs(r["theta"] - th) < 1e-9]
        axes[0, 1].plot([r["drive"] for r in ly_pts], [r["largest_exponent"] for r in ly_pts],
                        "o", color=col, ms=4, alpha=0.85,
                        label=labels[th] + r" ($\lambda_{\max}$, map)")
    axes[0, 1].axhline(0.0, color="k", lw=0.8, ls="--", alpha=0.6)
    axes[0, 1].set_ylabel("largest rate")
    axes[0, 1].set_xlabel(r"drive $E$")
    axes[0, 1].legend(fontsize=7.5, frameon=True)

    axes[1, 0].set_ylabel(r"$n_{\rm unstable}$ ($|\mu|>1$)")
    axes[1, 0].set_xlabel(r"drive $E$")
    axes[1, 0].set_yticks([0, 1, 2, 3, 4])

    axes[1, 1].set_ylabel(r"$\langle|\alpha|\rangle$")
    axes[1, 1].set_xlabel(r"drive $E$")

    for ax in axes.ravel():
        ax.set_xlim(0.0, 8.0)

    fig.tight_layout()
    out_pdf = ROOT / "figures" / "drive_axis_bifurcation.pdf"
    out_png = ROOT / "figures" / "drive_axis_bifurcation.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    manifest = {
        "gate": "DRIVE_AXIS_BIFURCATION_FIGURE",
        "input": [str(SCAN), str(LYAP)],
        "outputs": [str(out_pdf), str(out_png)],
        "note": ("(a) max Floquet multiplier vs drive for theta=0 and theta=pi/2 "
                 "(up/down continuation; E* markers); (b) largest Floquet rate with "
                 "the stored full-spectrum Lyapunov largest exponent overlaid at the "
                 "canonical drives; (c) number of unstable multipliers; (d) mean "
                 "cavity amplitude over the orbit period."),
    }
    (ROOT / "figures" / "drive_axis_bifurcation.manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "outputs": [str(out_pdf), str(out_png)]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
