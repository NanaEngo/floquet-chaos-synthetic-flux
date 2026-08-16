#!/usr/bin/env python3
"""Plot the observed Neimark--Sacker torus (Poincare sections).

Reads results/torus_observation.json and renders, for each flux phase, the
stroboscopic Poincare section of the perturbed trajectory just above the onset
(the born torus as an invariant closed curve) together with the largest
Lyapunov exponent as a function of drive.
"""
from __future__ import annotations
import argparse, json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "torus_observation.json"
FIG = ROOT / "figures"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, default=FIG / "torus_observation.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    data = json.loads(RES.read_text())
    records = data["records"]

    # For each phase, the point closest to onset is the cleanest torus section.
    by_phase: dict[str, list[dict]] = {}
    for r in records:
        by_phase.setdefault(str(r["theta_over_pi"]), []).append(r)

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 7.0))
    for col, phase in enumerate(sorted(by_phase, key=float)):
        rs = sorted(by_phase[phase], key=lambda r: r["drive"])
        # Poincare section at the closest-to-onset drive.
        r0 = rs[0]
        ax = axes[0, col]
        b1r = np.asarray(r0["poincare_points_b1r"])
        b1i = np.asarray(r0["poincare_points_b1i"])
        ax.scatter(b1r, b1i, s=2, alpha=0.5, color="#1f77b4")
        ax.set_title(rf"$\theta={r0['theta_over_pi']:.2f}\pi$, "
                     rf"$E={r0['drive']:.3f}$ ($\approx E^*+{r0['offset_above_onset']:.2f}$)"
                     + "\n" + rf"$\lambda_1={r0['largest_exponent']:+.3f}$, "
                     rf"$\lambda_2={r0['second_exponent']:+.3f}$", fontsize=8.5)
        ax.set_xlabel(r"$\mathrm{Re}\,\beta_1$")
        ax.set_ylabel(r"$\mathrm{Im}\,\beta_1$")
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)

        # Largest exponent vs drive.
        ax2 = axes[1, col]
        Es = [r["drive"] for r in rs]
        lams = [r["largest_exponent"] for r in rs]
        ax2.plot(Es, lams, "o-", color="#d62728")
        ax2.axhline(0.0, color="0.5", lw=0.8, ls="--")
        ax2.axvline(r0["onset_E"], color="0.5", lw=0.8, ls=":",
                    label=f"onset $E^*={r0['onset_E']:.3f}$")
        ax2.set_xlabel("drive $E$")
        ax2.set_ylabel(r"largest Lyapunov exponent $\lambda_1$")
        ax2.set_title("Near-marginal just above onset, growing with drive",
                      fontsize=8.5)
        ax2.legend(fontsize=7)

    axes[0, 0].set_ylabel(r"$\mathrm{Im}\,\beta_1$")
    axes[0, 1].set_ylabel(r"$\mathrm{Im}\,\beta_1$")
    fig.suptitle("Neimark--Sacker torus: stroboscopic section and largest "
                 "exponent", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    a.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.output)
    png = a.output.with_suffix(".png")
    fig.savefig(png, dpi=150)
    plt.close(fig)

    manifest = {
        "figure": "torus_observation",
        "source_json": str(RES.relative_to(ROOT)),
        "generator": str(Path(__file__).relative_to(ROOT)),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": [a.output.name, png.name],
        "caption_alt": (
            "Two-by-two figure of the observed Neimark-Sacker torus. Top row: "
            "stroboscopic Poincare sections of the perturbed drive-locked orbit "
            "just above the onset, projected onto the first mechanical mode, "
            "showing an invariant closed curve at flux phases zero and pi over "
            "two. Bottom row: the largest Lyapunov exponent versus drive, "
            "near-marginal just above the onset and growing with drive, with "
            "the onset drive marked."),
    }
    (FIG / "torus_observation.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "output": str(a.output),
                      "png": str(png)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
