#!/usr/bin/env python3
"""Generate the production flux-grid diagnostic figure from PASS JSON only."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.results.read_text(encoding="utf-8"))
    if data.get("status") != "PASS":
        raise SystemExit("ERROR: flux-grid figure requires a PASS result")
    records = data.get("records", [])
    if len(records) != data.get("points", 0) * data.get("replicates", 0):
        raise SystemExit("ERROR: incomplete flux-grid records")
    if any(record.get("status") != "PASS" for record in records):
        raise SystemExit("ERROR: all plotted flux-grid records must be PASS")

    theta = np.array([record["theta"] for record in records])
    max_exp = np.array([max(record["lyapunov_spectrum"]) for record in records])
    rate_diff = np.array([record["max_rate_difference"] for record in records])
    unique_theta = np.array(sorted(set(theta)))
    grouped_exp = [max_exp[np.isclose(theta, value)] for value in unique_theta]
    grouped_diff = [rate_diff[np.isclose(theta, value)] for value in unique_theta]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)
    axes[0].boxplot(grouped_exp, positions=unique_theta/np.pi, widths=0.07, manage_ticks=False)
    axes[0].axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    axes[0].set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    axes[0].set_ylabel(r"Largest QR exponent $\lambda_{\max}$")
    axes[0].text(0.02, 0.98, "(a)", transform=axes[0].transAxes,
                 va="top", ha="left", fontsize=11)
    axes[0].tick_params(axis="x", labelrotation=35)
    axes[1].boxplot(grouped_diff, positions=unique_theta/np.pi, widths=0.07, manage_ticks=False)
    axes[1].axhline(data["gate_thresholds"]["max_rate_difference"], color="#b2182b", linewidth=0.9, linestyle="--")
    axes[1].set_xlabel(r"flux phase $\Phi_{\rm syn}/\pi$")
    axes[1].set_ylabel(r"$\max|\rho_i-\lambda_i|$")
    axes[1].text(0.02, 0.98, "(b)", transform=axes[1].transAxes,
                 va="top", ha="left", fontsize=11)
    axes[1].tick_params(axis="x", labelrotation=35)
    stem = args.output_dir / "flux_grid_diagnostics"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    manifest = {
        "status": "PASS",
        "source": str(args.results.resolve()),
        "source_status": data["status"],
        "outputs": [str(stem.with_suffix(".pdf")), str(stem.with_suffix(".png"))],
        "points": data["points"],
        "replicates": data["replicates"],
        "no_mock_data": True,
        "interpretation": "Diagnostics for the sampled production grid; not evidence that chaos is absent outside the grid.",
    }
    (args.output_dir / "flux_grid_diagnostics.manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
