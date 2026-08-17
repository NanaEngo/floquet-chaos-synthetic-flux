#!/usr/bin/env python3
"""Full-spectrum transition map figure (weak coupling): (a) broken-axis six-exponent
spectrum vs flux phase separating the cavity and mechanical scales, (b) largest
exponent vs drive with the phase-to-phase spread shaded, (c) largest-exponent
heatmap over the drive-phase plane (everywhere negative).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


def load(p: Path):
    with open(p) as f:
        return json.load(f)


def main() -> int:
    # high-resolution phase data (drive=0.2) from the canonical flux grid
    fg = load(ROOT / "results" / "flux_grid.json")
    fg_recs = [r for r in fg["records"] if r.get("status") == "PASS"]
    thetas = sorted({round(r["theta"], 6) for r in fg_recs})
    spectra_by_phase = {}
    for t in thetas:
        recs = [r for r in fg_recs if round(r["theta"], 6) == t]
        mat = np.array([r["lyapunov_spectrum"] for r in recs])
        spectra_by_phase[t] = mat.mean(axis=0)

    # extended (drive x phase) full-spectrum map
    tm = load(ROOT / "results" / "full_spectrum_transition_map.json")
    tm_recs = [r for r in tm["records"] if r.get("status") == "PASS"]
    drives = sorted({r["drive"] for r in tm_recs})
    tm_thetas = sorted({round(r["theta"], 6) for r in tm_recs})

    x = np.array(thetas) / np.pi
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, 6))

    fig = plt.figure(figsize=(13.5, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1.0, 1.15], wspace=0.35)

    # ---- (a) broken-axis full spectrum vs phase ----
    gs_a = gs[0, 0].subgridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.12)
    ax_m = fig.add_subplot(gs_a[0])              # mechanical exponents (zoom)
    ax_c = fig.add_subplot(gs_a[1], sharex=ax_m)  # cavity exponents

    # four mechanical exponents (lambda3..lambda6)
    for k in range(2, 6):
        vals = [spectra_by_phase[t][k] for t in thetas]
        ax_m.plot(x, vals, "-o", ms=3.5, lw=1.5, color=colors[k],
                  label=rf"$\lambda_{{{k+1}}}$")
    ax_m.axhline(-0.01, color="0.5", lw=0.7, ls=":")
    ax_m.set_ylim(-0.0106, -0.0099)
    ax_m.set_ylabel(r"Lyapunov exponent $\lambda_k$")
    ax_m.legend(fontsize=7, ncol=4, frameon=False, loc="upper center",
                handletextpad=0.3, columnspacing=0.8)
    ax_m.text(0.02, 0.98, "(a)", transform=ax_m.transAxes,
              va="top", ha="left", fontsize=11)

    # two cavity exponents (lambda1, lambda2)
    for k in range(2):
        vals = [spectra_by_phase[t][k] for t in thetas]
        ax_c.plot(x, vals, "-o", ms=3.5, lw=1.5, color=colors[k],
                  label=rf"$\lambda_{{{k+1}}}$")
    ax_c.axhline(-0.5, color="0.5", lw=0.7, ls=":")
    ax_c.set_ylim(-0.505, -0.495)
    ax_c.set_xlabel(r"synthetic flux phase $\Phi_{\rm syn}/\pi$")
    ax_c.legend(fontsize=7, ncol=2, frameon=False, loc="upper center",
                handletextpad=0.3, columnspacing=0.8)

    # family labels at the reference lines (blended axes-fraction x / data y)
    blend_m = matplotlib.transforms.blended_transform_factory(
        ax_m.transAxes, ax_m.transData)
    blend_c = matplotlib.transforms.blended_transform_factory(
        ax_c.transAxes, ax_c.transData)
    ax_m.text(0.985, -0.00995, r"mechanical $(-\gamma/2)$", transform=blend_m,
              va="bottom", ha="right", color="0.3", fontsize=8)
    ax_c.text(0.985, -0.5005, r"cavity $(-\kappa/2)$", transform=blend_c,
              va="bottom", ha="right", color="0.3", fontsize=8)

    # broken-axis: hide the facing spines and draw the break marks
    ax_m.spines["bottom"].set_visible(False)
    ax_c.spines["top"].set_visible(False)
    ax_m.xaxis.tick_top()
    ax_m.tick_params(labeltop=False)
    ax_c.xaxis.tick_bottom()
    d = 0.02
    kw = dict(color="k", lw=1.0, clip_on=False)
    ax_m.plot([-d, d], [-d, d], transform=ax_m.transAxes, **kw)
    ax_c.plot([-d, d], [1 - d, 1 + d], transform=ax_c.transAxes, **kw)

    # ---- (b) largest exponent vs drive ----
    axb = fig.add_subplot(gs[0, 1])
    drive_largest = {}
    drive_lowest = {}
    for d in drives:
        vals = [r["largest_exponent"] for r in tm_recs if r["drive"] == d]
        drive_largest[d] = float(np.max(vals))
        drive_lowest[d] = float(np.min(vals))
    xs = [d for d in drives]
    ymax = [drive_largest[d] for d in drives]
    ymin = [drive_lowest[d] for d in drives]
    axb.fill_between(xs, ymin, ymax, color="crimson", alpha=0.18,
                     lw=0, label="phase spread")
    axb.plot(xs, ymax, "-o", color="crimson", lw=1.8,
             label="largest over phases")
    axb.axhline(0.0, color="black", lw=0.8, ls="--")
    axb.set_xlabel(r"drive amplitude $E$")
    axb.set_ylabel(r"largest exponent $\lambda_{\max}$")
    axb.set_xscale("log")
    axb.set_xticks(drives)
    axb.set_xticklabels([f"{d:g}" for d in drives])
    axb.xaxis.set_minor_locator(plt.NullLocator())
    axb.set_ylim(-0.052, 0.006)
    axb.text(0.02, 0.98, "(b)", transform=axb.transAxes,
             va="top", ha="left", fontsize=11)
    axb.text(0.5, 0.90, "no positive exponent", transform=axb.transAxes,
             ha="center", color="crimson", fontsize=8)
    axb.legend(fontsize=7, frameon=False, loc="lower left")

    # ---- (c) largest-exponent heatmap over the drive-phase plane ----
    axc = fig.add_subplot(gs[0, 2])
    Z = np.full((len(drives), len(tm_thetas)), np.nan)
    for di, d in enumerate(drives):
        for ti, t in enumerate(tm_thetas):
            rs = [r for r in tm_recs if r["drive"] == d and round(r["theta"], 6) == t]
            if rs:
                Z[di, ti] = rs[0]["largest_exponent"]
    im = axc.imshow(Z, aspect="auto", origin="lower",
                    extent=[tm_thetas[0] / np.pi, tm_thetas[-1] / np.pi,
                            drives[0], drives[-1]],
                    cmap="Blues", vmin=-0.05, vmax=0, interpolation="nearest")
    axc.set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    axc.set_ylabel("drive amplitude $E$")
    axc.set_yscale("log")
    axc.set_yticks(drives)
    axc.set_yticklabels([f"{d:g}" for d in drives])
    axc.yaxis.set_minor_locator(plt.NullLocator())
    axc.text(0.02, 0.98, "(c)", transform=axc.transAxes,
             va="top", ha="left", fontsize=11)
    cb = fig.colorbar(im, ax=axc, ticks=[-0.04, -0.03, -0.02, -0.01, 0])
    cb.set_label(r"largest exponent $\lambda_{\max}$")
    axc.text(0.5, 0.5, "no positive\nexponent", transform=axc.transAxes,
             ha="center", va="center", color="0.15", fontsize=8)

    fig.savefig(FIG / "full_spectrum_transition_map.pdf", bbox_inches="tight")
    fig.savefig(FIG / "full_spectrum_transition_map.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    # manifest
    manifest = {
        "gate": "FULL_SPECTRUM_TRANSITION_MAP_FIGURE",
        "inputs": ["results/flux_grid.json", "results/full_spectrum_transition_map.json"],
        "outputs": ["figures/full_spectrum_transition_map.pdf",
                    "figures/full_spectrum_transition_map.png"],
        "n_flux_phases": len(thetas),
        "n_drives": len(drives),
        "any_positive_exponent": bool(tm.get("any_positive_exponent", False)),
        "note": "Figure (a) shows all six exponents on a broken vertical scale: the two "
                "cavity-dominated exponents near -kappa/2 and the four mechanical exponents "
                "near -gamma/2; (b) the largest (most marginal) exponent vs drive with the "
                "phase-to-phase spread shaded; (c) the largest exponent over the drive-phase "
                "plane, everywhere negative (no positive exponent across the tested domain).",
    }
    (FIG / "full_spectrum_transition_map.manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print("wrote", FIG / "full_spectrum_transition_map.pdf")
    print("any_positive_exponent =", manifest["any_positive_exponent"],
          "| n_flux_phases =", manifest["n_flux_phases"],
          "| n_drives =", manifest["n_drives"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
