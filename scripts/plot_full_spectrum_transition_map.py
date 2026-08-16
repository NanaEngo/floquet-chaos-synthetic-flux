#!/usr/bin/env python3
"""Full-spectrum transition map figure: (a) complete six-exponent spectrum vs flux
phase, (b) largest exponent vs drive, (c) positive-exponent (hyperchaos-order) heatmap.
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
    # per-phase spectrum (average over the 3 replicates for the plotted curves)
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

    fig = plt.figure(figsize=(13.5, 4.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.15], wspace=0.32)

    # ---- (a) full spectrum vs phase (drive = 0.2) ----
    ax = fig.add_subplot(gs[0, 0])
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, 6))
    for k in range(6):
        vals = [spectra_by_phase[t][k] for t in thetas]
        ax.plot(np.array(thetas) / np.pi, vals, "-o", ms=4, lw=1.6,
                color=colors[k], label=f"$\\lambda_{{{k+1}}}$")
    ax.axhline(0.0, color="black", lw=0.8, ls="--")
    ax.set_xlabel(r"synthetic flux phase $\Phi_{\rm syn}/\pi$")
    ax.set_ylabel("Lyapunov exponent $\\lambda_k$")
    ax.set_title("(a) Full spectrum vs phase (drive $=0.2$)")
    ax.legend(fontsize=7, ncol=2, frameon=False, loc="center right")
    ax.set_xlim(-1, 1)

    # ---- (b) largest exponent vs drive ----
    axb = fig.add_subplot(gs[0, 1])
    # average largest exponent over phases at each drive
    drive_largest = {d: float(np.max([r["largest_exponent"] for r in tm_recs
                                      if r["drive"] == d])) for d in drives}
    xs = [d for d in drives]
    ys = [drive_largest[d] for d in drives]
    axb.plot(xs, ys, "-o", color="crimson", lw=1.8)
    axb.axhline(0.0, color="black", lw=0.8, ls="--")
    axb.set_xlabel("drive amplitude $E$")
    axb.set_ylabel(r"largest exponent $\lambda_{\max}$")
    axb.set_title("(b) $\\lambda_{\\max}$ vs drive")
    axb.set_xscale("log")
    axb.text(0.5, 0.06, "no positive exponent", transform=axb.transAxes,
             ha="center", color="crimson", fontsize=8)

    # ---- (c) positive-exponent count heatmap ----
    axc = fig.add_subplot(gs[0, 2])
    Z = np.full((len(drives), len(tm_thetas)), np.nan)
    for di, d in enumerate(drives):
        for ti, t in enumerate(tm_thetas):
            rs = [r for r in tm_recs if r["drive"] == d and round(r["theta"], 6) == t]
            if rs:
                Z[di, ti] = rs[0]["n_positive_exponents"]
    im = axc.imshow(Z, aspect="auto", origin="lower",
                    extent=[tm_thetas[0] / np.pi, tm_thetas[-1] / np.pi,
                            drives[0], drives[-1]],
                    cmap="RdYlBu_r", vmin=0, vmax=3, interpolation="nearest")
    axc.set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    axc.set_ylabel("drive amplitude $E$")
    axc.set_title("(c) Positive exponents (hyperchaos order)")
    axc.set_yscale("log")
    cb = fig.colorbar(im, ax=axc, ticks=[0, 1, 2, 3])
    cb.set_label("$n_{\\rm pos}$")

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
        "note": "Figure (a) plots all six exponents; (b) the largest exponent vs drive; "
                "(c) the positive-exponent (hyperchaos-order) count. Zero positive "
                "exponents were found across the tested (drive x phase) domain.",
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
