#!/usr/bin/env python3
"""Plot the reduced-model (adiabatic) two-objective Pareto front (v4).

Reads ``results/physical_parameter_reduced_bad_cavity_pareto_v4.json`` and
renders the sixteen feasible candidates in (drive cost, stability margin)
space, with the three non-dominated points highlighted and joined by the
front. The declared objectives are ``minimize drive cost`` and ``maximize
robust stability margin``; Fisher information was not computed in this screen.
A figure manifest is written next to the PDF/PNG for provenance.
"""
from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results" / "physical_parameter_reduced_bad_cavity_pareto_v4.json"
FIG = ROOT / "figures"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=FIG / "reduced_model_pareto_v4.pdf")
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("seaborn-v0_8-whitegrid")

    data = json.loads(RES.read_text(encoding="utf-8"))
    if data.get("scientific_status") != "PROVISIONAL_REDUCED_MODEL_MULTI_START_PARETO":
        raise SystemExit(f"refusing to plot status={data.get('scientific_status')!r}")

    cands = data["candidates"]
    pareto = sorted(data["pareto_candidates"], key=lambda c: c["drive_cost"])

    cost_all = np.array([c["drive_cost"] for c in cands], float)
    margin_all = np.array([c["minimum_stability_margin"] for c in cands], float)
    cost_p = np.array([c["drive_cost"] for c in pareto], float)
    margin_p = np.array([c["minimum_stability_margin"] for c in pareto], float)

    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.scatter(cost_all, margin_all, s=60, color="0.55", alpha=0.6,
               edgecolor="black", linewidth=0.4, label="Feasible (16)")
    ax.plot(cost_p, margin_p, "--", color="crimson", lw=1.2,
            label="Non-dominated front")
    ax.scatter(cost_p, margin_p, marker="*", s=230, facecolor="none",
               edgecolor="crimson", linewidth=1.5, label="Non-dominated (3)")
    for c, m in zip(cost_p, margin_p):
        ax.annotate(f"{c:.4f}", (c, m), textcoords="offset points",
                    xytext=(7, 7), fontsize=8.5, color="crimson")

    ax.set_xlabel("Normalized drive cost (drive amplitude)")
    ax.set_ylabel(r"Worst-replicate stability margin $1-\max|\mu|$")
    ax.text(0.03, 0.03, "PROVISIONAL REDUCED MODEL; not SI-calibrated",
            transform=ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.5"})
    ax.legend(loc="lower right", frameon=True, fontsize=8)
    fig.tight_layout()

    a.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.output, bbox_inches="tight")
    png = a.output.with_suffix(".png")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "figure": "reduced_model_pareto_v4",
        "source_json": str(RES.relative_to(ROOT)),
        "generator": str(Path(__file__).relative_to(ROOT)),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "outputs": [a.output.name, png.name],
        "caption_alt": (
            "Scatter of sixteen feasible candidates of the reduced adiabatic "
            "model in normalized drive cost versus worst-replicate stability "
            "margin, with the three non-dominated candidates marked by crimson "
            "stars and joined by a dashed front. The front trades increasing "
            "drive cost for increasing stability margin. Model-level only."),
    }
    (FIG / "reduced_model_pareto_v4.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    print(json.dumps({"status": "PASS", "candidates": len(cands),
                      "pareto": len(pareto), "output": str(a.output),
                      "png": str(png)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
