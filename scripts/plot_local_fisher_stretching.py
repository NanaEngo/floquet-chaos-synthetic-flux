#!/usr/bin/env python3
"""Local susceptibility--stretching covariation figure.

Three panels from ``results/local_fisher_stretching.json``: (a) the local
force-to-observable susceptibility |S_k| against the local finite-time
Lyapunov exponent at the most hyperchaotic point; (b) the local record spread
sigma_k against the same local exponent (nearly flat response); (c) grouped
box plots of the local Fisher proxy log10(S_k^2/sigma_k^2) at all six
(drive, flux-phase) points.
"""
from __future__ import annotations
import argparse, json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "local_fisher_stretching.json"
FIG = ROOT / "figures"


def theta_label(t: float) -> str:
    return {0.0: r"$\theta=0$", 0.5: r"$\theta=\pi/2$",
            1.0: r"$\theta=\pi$"}.get(t, rf"$\theta={t:.2g}\pi$")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=FIG / "local_fisher_stretching.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        import scienceplots  # noqa: F401
        plt.style.use(["science", "no-latex"])
    except Exception:
        plt.style.use("seaborn-v0_8-whitegrid")

    data = json.loads(RES.read_text())
    rec = {(r["drive"], r["theta_over_pi"]): r for r in data["results"]}
    rep = rec[(8.0, 0.5)]

    fig, (axa, axb, axc) = plt.subplots(
        1, 3, figsize=(10.2, 3.3),
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.25]})

    lam = np.array(rep["windows"]["lambda_1"])
    sig = np.abs(np.array(rep["windows"]["susceptibility_S"]))
    sgd = np.array(rep["windows"]["sigma_ar"])

    axa.scatter(lam, sig, s=6, alpha=0.55, color="#1f77b4", linewidths=0)
    sl = np.polyfit(lam, np.log10(sig), 1)
    xs = np.linspace(lam.min(), lam.max(), 50)
    axa.plot(xs, 10 ** (sl[0] * xs + sl[1]), color="#d62728", lw=1.2)
    axa.set_xlabel(r"local $\lambda_k$  (1/time)")
    axa.set_ylabel(r"$|S_k| = |\partial_F \langle X_a\rangle_k|$")
    axa.set_yscale("log")
    axa.set_title("(a) signal vs. local stretching", fontsize=9)
    axa.annotate(rf"$\mathrm{{d}}\ln|S|/\mathrm{{d}}\lambda={sl[0]*np.log(10):.1f}$",
                 xy=(0.03, 0.94), xycoords="axes fraction", fontsize=8)

    axb.scatter(lam, sgd, s=6, alpha=0.55, color="#2ca02c", linewidths=0)
    sl2 = np.polyfit(lam, np.log(sgd), 1)[0]
    axb.set_xlabel(r"local $\lambda_k$  (1/time)")
    axb.set_ylabel(r"record spread $\sigma_k$")
    axb.set_title("(b) noise vs. local stretching", fontsize=9)
    axb.annotate(rf"$\mathrm{{d}}\ln\sigma/\mathrm{{d}}\lambda={sl2:.2f}$",
                 xy=(0.03, 0.94), xycoords="axes fraction", fontsize=8)

    drives = (4.0, 8.0)
    thetas = (0.0, 0.5, 1.0)
    groups, labels = [], []
    for d in drives:
        for t in thetas:
            w = rec[(d, t)]["windows"]
            f = (np.abs(np.array(w["susceptibility_S"])) ** 2
                 / np.array(w["sigma_ar"]) ** 2)
            groups.append(np.log10(f))
            labels.append(theta_label(t))
    bp = axc.boxplot(groups, tick_labels=labels, showfliers=False, widths=0.6,
                     patch_artist=True, medianprops=dict(color="black"))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor("#4c72b0" if i < 3 else "#dd8452")
        patch.set_alpha(0.75)
    axc.set_xticklabels(labels, rotation=0)
    axc.set_xlabel("flux phase", labelpad=2)
    axc.set_ylabel(r"$\log_{10}$ local Fisher proxy $S_k^2/\sigma_k^2$")
    axc.set_title("(c) local Fisher budget by phase", fontsize=9)
    axc.axvline(3.5, color="0.6", lw=0.8, ls="--")
    axc.text(2.0, 0.97, r"$E=4$", transform=axc.get_xaxis_transform(),
             ha="center", va="top", fontsize=8)
    axc.text(5.0, 0.97, r"$E=8$", transform=axc.get_xaxis_transform(),
             ha="center", va="top", fontsize=8)

    fig.tight_layout()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.output)
    png = a.output.with_suffix(".png")
    fig.savefig(png, dpi=300)
    manifest = {
        "figure": "local_fisher_stretching",
        "source_json": str(RES.relative_to(ROOT)),
        "generator": str(Path(__file__).relative_to(ROOT)),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": [a.output.name, png.name],
        "caption_alt": (
            "Three panels of the window-resolved local analysis at the "
            "strong-coupling points. (a) The local force-to-observable "
            "susceptibility against the local finite-time Lyapunov exponent at "
            "E=8, theta=pi/2; the signal responds strongly to local "
            "stretching. (b) The local record spread against the same local "
            "exponent; the noise responds weakly because the within-window "
            "record saturates at the attractor scale. (c) Grouped box plots of "
            "the local Fisher proxy logarithm at all six drive and flux-phase "
            "points, comparing the flux-phase dependence of the per-window "
            "signal-to-noise budget."),
    }
    (FIG / "local_fisher_stretching.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"wrote {a.output} and {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
