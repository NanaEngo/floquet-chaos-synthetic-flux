#!/usr/bin/env python3
"""Hyperchaos transition map figure: full spectrum vs phase, positive-exponent
(hyperchaos-order) count vs phase, and its heatmap over the drive-phase plane.
"""
from __future__ import annotations

import json
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
    d = load(ROOT / "results" / "full_spectrum_hyperchaos_map_g03.json")
    recs = [r for r in d["records"] if r.get("status") == "PASS"]
    drives = d["drives"]
    thetas = sorted({round(r["theta"], 6) for r in recs})
    th_pi = np.array(thetas) / np.pi

    def rec_at(drive, theta):
        for r in recs:
            if r["drive"] == drive and round(r["theta"], 6) == round(theta, 6):
                return r
        return None

    # ---- (a) full spectrum vs phase at E=4.0 ----
    E4 = 4.0
    spec_mat = np.array([rec_at(E4, t)["lyapunov_spectrum"] for t in thetas])

    fig = plt.figure(figsize=(14.0, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.2], wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, 6))
    for k in range(6):
        vals = spec_mat[:, k]
        ax.plot(th_pi, vals, "-o", ms=4, lw=1.6, color=colors[k],
                label=f"$\\lambda_{{{k+1}}}$")
    ax.axhline(0.0, color="black", lw=0.8, ls="--")
    ax.fill_between(th_pi, 0, np.max(spec_mat, axis=1), color="crimson",
                    alpha=0.08)
    ax.set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    ax.set_ylabel("Lyapunov exponent $\\lambda_k$")
    ax.set_title("(a) Full spectrum vs phase ($g=0.3$, $E=4$)")
    ax.legend(fontsize=7, ncol=2, frameon=False, loc="upper center")
    ax.set_xlim(-1, 1)

    # ---- (b) n_pos vs phase ----
    axb = fig.add_subplot(gs[0, 1])
    for dr, marker, label in [(4.0, "o-", "$E=4$"), (8.0, "s--", "$E=8$"),
                              (0.2, "^-", "$E=0.2$")]:
        ys = [rec_at(dr, t)["n_positive_exponents"] for t in thetas]
        axb.plot(th_pi, ys, marker, lw=1.6, ms=5, label=label)
    axb.set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    axb.set_ylabel("positive exponents $n_{\\rm pos}$")
    axb.set_title("(b) Hyperchaos order vs phase")
    axb.legend(fontsize=8, frameon=False)
    axb.set_ylim(-0.3, 4.6)
    axb.set_xlim(-1, 1)
    axb.set_yticks(range(0, 5))

    # ---- (c) heatmap n_pos vs (drive, phase) ----
    axc = fig.add_subplot(gs[0, 2])
    Z = np.full((len(drives), len(thetas)), np.nan)
    for di, dr in enumerate(drives):
        for ti, t in enumerate(thetas):
            r = rec_at(dr, t)
            if r:
                Z[di, ti] = r["n_positive_exponents"]
    im = axc.imshow(Z, aspect="auto", origin="lower",
                    extent=[th_pi[0], th_pi[-1], drives[0], drives[-1]],
                    cmap="inferno", vmin=0, vmax=4, interpolation="nearest")
    axc.set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    axc.set_ylabel("drive amplitude $E$")
    axc.set_title("(c) Hyperchaos order $n_{\\rm pos}(E,\\Phi)$")
    axc.set_yscale("log")
    cb = fig.colorbar(im, ax=axc, ticks=[0, 1, 2, 3, 4])
    cb.set_label("$n_{\\rm pos}$")

    fig.savefig(FIG / "hyperchaos_transition_map.pdf", bbox_inches="tight")
    fig.savefig(FIG / "hyperchaos_transition_map.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    # manifest
    man = {
        "gate": "HYPERCHAOS_TRANSITION_MAP_FIGURE",
        "input": "results/full_spectrum_hyperchaos_map_g03.json",
        "outputs": ["figures/hyperchaos_transition_map.pdf",
                    "figures/hyperchaos_transition_map.png"],
        "coupling": d.get("coupling"),
        "max_positive_exponents": d.get("max_positive_exponents"),
        "note": ("(a) full six-exponent spectrum vs flux phase at g=0.3, E=4; "
                 "(b) hyperchaos order vs phase at E=0.2/4/8; (c) heatmap of the "
                 "hyperchaos order over the drive-phase plane. The number of positive "
                 "exponents is controlled by the synthetic-flux phase."),
    }
    (FIG / "hyperchaos_transition_map.manifest.json").write_text(
        json.dumps(man, indent=2) + "\n", encoding="utf-8")
    print("wrote", FIG / "hyperchaos_transition_map.pdf")
    print("max n_pos =", d.get("max_positive_exponents"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
