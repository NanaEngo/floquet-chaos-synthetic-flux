#!/usr/bin/env python3
"""Figure: noise-dependent observability of the chaos signature.

Panels:
(a) measured-record variance of y = X_a + detector noise versus drive amplitude
    for theta = 0 and theta = pi/2, with the noise-only floor marked;
(b) the observability ratio (measured variance / noise floor) versus drive,
    with the unit line marking the noise floor.
"""
from __future__ import annotations
import json, sys
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
DATA = ROOT / "results" / "noise_observability.json"


def main() -> int:
    d = json.loads(DATA.read_text())
    floor = d["noise_floor_measured_record_var"]

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
    colors = {0.0: "#1f77b4", np.pi / 2: "#d62728"}
    labels = {0.0: r"$\theta=0$", np.pi / 2: r"$\theta=\pi/2$"}

    for th in (0.0, np.pi / 2):
        recs = [r for r in d["records"] if abs(r["theta"] - th) < 1e-9]
        recs.sort(key=lambda r: r["drive"])
        drives = [r["drive"] for r in recs]
        var = [r["measured_record_var"] for r in recs]
        ratio = [r["observability_ratio_vs_noise_floor"] for r in recs]
        axes[0].plot(drives, var, "o-", color=colors[th], label=labels[th], lw=1.2, ms=4)
        axes[1].plot(drives, ratio, "o-", color=colors[th], label=labels[th], lw=1.2, ms=4)

    axes[0].axhline(floor, color="k", ls="--", lw=0.9, alpha=0.6,
                    label="noise floor")
    axes[1].axhline(1.0, color="k", ls="--", lw=0.9, alpha=0.6)
    axes[0].set_xlabel(r"drive $E$")
    axes[0].set_ylabel(r"measured-record variance")
    axes[1].set_xlabel(r"drive $E$")
    axes[1].set_ylabel(r"observability ratio")
    axes[0].legend(fontsize=7, frameon=True)
    for ax in axes:
        ax.set_xlim(0.0, 8.5)
    axes[1].set_yscale("log")

    fig.tight_layout()
    out_pdf = ROOT / "figures" / "noise_observability.pdf"
    out_png = ROOT / "figures" / "noise_observability.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=200)
    manifest = {
        "gate": "NOISE_OBSERVABILITY_FIGURE",
        "input": str(DATA),
        "outputs": [str(out_pdf), str(out_png)],
        "note": ("(a) measured-record variance of y = X_a + detector noise vs drive "
                 "for theta=0 and pi/2 with the noise-only floor marked; (b) the "
                 "observability ratio (measured variance / noise floor) vs drive."),
    }
    (ROOT / "figures" / "noise_observability.manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "outputs": [str(out_pdf), str(out_png)]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
