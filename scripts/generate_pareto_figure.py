#!/usr/bin/env python3
"""Generate a provenance-bound diagnostic plot for the provisional Pareto screen."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    if data.get("status") not in {"PASS", "PROVISIONAL"}:
        raise SystemExit(f"refusing to plot status={data.get('status')!r}")
    records = [r for r in data.get("candidate_records", []) if r.get("status") == "PASS" and r.get("feasible")]
    if not records:
        raise SystemExit("no feasible records to plot")
    fisher = np.array([r["metrics"]["classical_fisher_information"] for r in records], float)
    margin = np.array([r["metrics"]["robustness_margin"] for r in records], float)
    cost = np.array([r["metrics"]["resource_cost"] for r in records], float)
    pareto_keys = {
        json.dumps(r.get("normalized_parameters", {}), sort_keys=True)
        for r in data.get("pareto_candidates", [])
    }
    pareto_mask = np.array([
        json.dumps(r.get("normalized_parameters", {}), sort_keys=True) in pareto_keys
        for r in records
    ], dtype=bool)
    model_level = "model_recovery" in str(data.get("manifest", ""))
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    sc = ax.scatter(margin, fisher, c=cost, cmap="viridis", s=75, alpha=0.55,
                    edgecolor="black", linewidth=0.4, label="Feasible")
    if np.any(pareto_mask):
        ax.scatter(margin[pareto_mask], fisher[pareto_mask], marker="*", s=180,
                   facecolor="none", edgecolor="crimson", linewidth=1.2,
                   label="Non-dominated")
    ax.set_yscale("log")
    ax.set_xlabel("Worst-replicate stability margin (normalized rate)")
    ax.set_ylabel("Classical measurement Fisher information")
    ax.set_title("Provisional normalized-model Pareto recovery" if model_level else "Provisional literature-anchored Pareto screening")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Effective normalized drive coordinate")
    annotation = ("PROVISIONAL_MODEL_LEVEL; not SI-calibrated"
                  if model_level else "PROVISIONAL; not an experimental optimum")
    ax.text(0.02, 0.02, annotation, transform=ax.transAxes,
            fontsize=8, bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.5"})
    ax.legend(loc="best", frameon=True, fontsize=8)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({"status": "PASS", "records": len(records), "pareto_records": int(np.sum(pareto_mask)), "output": str(args.output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
