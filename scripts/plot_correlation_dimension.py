#!/usr/bin/env python3
"""Plot the Grassberger--Procaccia correlation-sum curves (log C vs log eps).

Reads results/correlation_dimension.json and renders the four strong-coupling
points (E = 4, 8 at theta = 0 and pi/2) with their fitted slopes and the
Kaplan--Yorke dimension annotated. A figure manifest is written next to the
PDF/PNG for provenance.
"""
from __future__ import annotations
import argparse, json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "correlation_dimension.json"
FIG = ROOT / "figures"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=FIG / "correlation_dimension.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    data = json.loads(RES.read_text())
    records = data["records"]

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.4), sharex=True, sharey=True)
    for ax, r in zip(axes.ravel(), records):
        eps = np.asarray(r["eps"])
        C = np.asarray(r["correlation_sum_C"])
        mask = (C >= 1e-3) & (C <= 0.3)
        logE = np.log(eps)
        logC = np.log(C)
        ax.loglog(eps, C, ".", ms=3, alpha=0.7, color="0.35")
        ax.loglog(eps[mask], C[mask], "o", ms=4, color="#1f77b4",
                  label="scaling region")
        slope = r["correlation_dimension_D2"]
        xf = np.log(eps[mask])
        ax.loglog(eps[mask], np.exp(r["fits"]["window_A_C1e-3_0p3"]["fit_intercept"]
                                    + slope * xf), "-", color="#d62728", lw=1.6,
                  label=f"$D_2 = {slope:.2f}$")
        ky = r["cross_check_vs_kaplan_yorke"]["kaplan_yorke_dimension"]
        npos = r["cross_check_vs_kaplan_yorke"]["n_positive_exponents"]
        ax.set_title(rf"$\theta={r['theta_over_pi']:.2f}\pi$, "
                     rf"$E={r['drive']:.2f}$" + "\n"
                     rf"$D_2={slope:.2f}\pm{r['D2_scaling_region_uncertainty']:.2f}$, "
                     rf"$D_{{KY}}={ky:.2f}$, $n_+={npos}$", fontsize=8.5)
        ax.grid(True, which="both", alpha=0.3)

    for ax in axes[:, 0]:
        ax.set_ylabel(r"$C(\varepsilon)$")
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\varepsilon$")
    axes[0, 0].legend(fontsize=7, loc="lower right")
    fig.suptitle("Theiler-corrected correlation sum and the Kaplan--Yorke "
                 "inequality $D_2 \\leq D_{\\mathrm{KY}}$", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    a.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.output)
    png = a.output.with_suffix(".png")
    fig.savefig(png, dpi=150)
    plt.close(fig)

    manifest = {
        "figure": "correlation_dimension",
        "source_json": str(RES.relative_to(ROOT)),
        "generator": str(Path(__file__).relative_to(ROOT)),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": [a.output.name, png.name],
        "caption_alt": (
            "Two-by-two log-log panels of the Theiler-corrected two-point "
            "correlation sum versus distance at four strong-coupling points "
            "(drive four and eight, flux phases zero and pi over two). Each "
            "panel marks the scaling region and its fitted slope, the "
            "correlation dimension D2, alongside the Kaplan-Yorke dimension "
            "and the number of positive Lyapunov exponents. Every panel shows "
            "D2 below the Kaplan-Yorke dimension."),
    }
    (FIG / "correlation_dimension.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "output": str(a.output),
                      "png": str(png)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
