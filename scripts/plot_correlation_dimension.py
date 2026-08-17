#!/usr/bin/env python3
"""Correlation-dimension figure.

Two panels: (a) the Theiler-corrected two-point correlation sum C(eps) versus
distance at the most hyperchaotic strong-coupling point, with the scaling region
and its fitted slope (the correlation dimension D2); (b) the correlation
dimension D2 (with scaling-region uncertainty) against the Kaplan--Yorke
dimension at the four strong-coupling points, showing D2 <= D_KY everywhere.
"""
from __future__ import annotations
import argparse, json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "correlation_dimension.json"
FIG = ROOT / "figures"


def theta_label(theta_over_pi: float) -> str:
    if abs(theta_over_pi) < 0.01:
        return r"$\theta=0$"
    if abs(theta_over_pi - 0.5) < 0.01:
        return r"$\theta=\pi/2$"
    return rf"$\theta={theta_over_pi:.2f}\pi$"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=FIG / "correlation_dimension.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = json.loads(RES.read_text())
    records = data["records"]

    # representative point: the most hyperchaotic (largest positive-exponent count)
    rep = max(records,
              key=lambda r: r["cross_check_vs_kaplan_yorke"]["n_positive_exponents"])

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(9.0, 3.9),
                                   gridspec_kw={"width_ratios": [1.15, 1.0]})

    # ---- (a) correlation-sum scaling at the representative point ----
    eps = np.asarray(rep["eps"])
    C = np.asarray(rep["correlation_sum_C"])
    C_floor = 1e-10  # exclude the no-pairs underflow (C ~ 1e-300)
    plot_mask = C >= C_floor
    mask = (C >= 1e-3) & (C <= 0.3)  # scaling region (fit window A)
    slope = rep["correlation_dimension_D2"]
    ky = rep["cross_check_vs_kaplan_yorke"]["kaplan_yorke_dimension"]
    npos = rep["cross_check_vs_kaplan_yorke"]["n_positive_exponents"]
    unc = rep["D2_scaling_region_uncertainty"]

    axa.loglog(eps[plot_mask], C[plot_mask], ".", ms=3, alpha=0.7, color="0.35")
    axa.loglog(eps[mask], C[mask], "o", ms=4, color="#1f77b4",
               label="scaling region")
    xf = np.log(eps[mask])
    axa.loglog(eps[mask],
               np.exp(rep["fits"]["window_A_C1e-3_0p3"]["fit_intercept"]
                      + slope * xf),
               "-", color="#d62728", lw=1.6,
               label=rf"fit, $D_2={slope:.2f}$")
    axa.text(0.02, 0.98, "(a)", transform=axa.transAxes,
             va="top", ha="left", fontsize=11)
    axa.text(0.03, 0.82,
             rf"{theta_label(rep['theta_over_pi'])}, $E={rep['drive']:.0f}$" + "\n"
             rf"$D_2={slope:.2f}\pm{unc:.2f}$, $D_{{\mathrm{{KY}}}}={ky:.2f}$, $n_+={npos}$",
             transform=axa.transAxes, va="top", ha="left", fontsize=8.5)
    axa.set_xlabel(r"distance $\varepsilon$")
    axa.set_ylabel(r"correlation sum $C(\varepsilon)$")
    axa.legend(fontsize=7, loc="lower right", frameon=False)

    # ---- (b) D2 vs D_KY at the four points ----
    D2 = np.array([r["correlation_dimension_D2"] for r in records])
    D2e = np.array([r["D2_scaling_region_uncertainty"] for r in records])
    DKY = np.array([r["cross_check_vs_kaplan_yorke"]["kaplan_yorke_dimension"]
                    for r in records])
    x = np.arange(len(records))
    width = 0.34
    axb.bar(x - width / 2, D2, width, yerr=D2e, capsize=3,
            color="#1f77b4", label=r"$D_2$ (correlation)")
    axb.bar(x + width / 2, DKY, width,
            color="#d62728", label=r"$D_{\mathrm{KY}}$ (Kaplan--Yorke)")
    labels = [rf"{theta_label(r['theta_over_pi'])}, $E={r['drive']:.0f}$"
              for r in records]
    axb.set_xticks(x)
    axb.set_xticklabels(labels, fontsize=8)
    axb.set_ylabel("dimension")
    axb.set_ylim(0, 5.6)
    axb.text(0.02, 0.98, "(b)", transform=axb.transAxes,
             va="top", ha="left", fontsize=11)
    axb.legend(fontsize=8, loc="upper left", frameon=False)

    fig.tight_layout()

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
            "Two-panel correlation-dimension cross-check. Left: the "
            "Theiler-corrected two-point correlation sum versus distance at the "
            "most hyperchaotic point, with the scaling region and its fitted "
            "slope giving the correlation dimension D2. Right: grouped bars of "
            "the correlation dimension D2 with its scaling-region uncertainty "
            "against the Kaplan-Yorke dimension at the four strong-coupling "
            "points, showing D2 below the Kaplan-Yorke dimension at every point."),
    }
    (FIG / "correlation_dimension.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "output": str(a.output),
                      "png": str(png)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
