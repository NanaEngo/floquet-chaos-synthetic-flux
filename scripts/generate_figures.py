#!/usr/bin/env python3
"""Generate reconstruction figures from immutable, passed result JSON only."""
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
    result_path = args.results.resolve()
    if not result_path.is_file():
        raise SystemExit(f"ERROR: missing result file: {result_path}")
    data = json.loads(result_path.read_text(encoding="utf-8"))
    if data.get("status") != "PASS":
        raise SystemExit("ERROR: figure generation requires status=PASS")
    multiplier_records = data["floquet"]["multipliers"]
    multipliers = np.asarray([complex(z["real"], z["imag"]) for z in multiplier_records], dtype=complex)
    rates = np.asarray(data["comparison"]["floquet_rates_sorted"], dtype=float)
    spectrum = np.asarray(data["comparison"]["lyapunov_spectrum_sorted"], dtype=float)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    unit = np.linspace(0, 2*np.pi, 400)
    axes[0].plot(np.cos(unit), np.sin(unit), color="0.65", lw=1, label="unit circle")
    axes[0].scatter(multipliers.real, multipliers.imag, s=24, color="#1f77b4")
    axes[0].axhline(0, color="0.85", lw=0.7); axes[0].axvline(0, color="0.85", lw=0.7)
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_xlabel("Re($\\mu$)"); axes[0].set_ylabel("Im($\\mu$)")
    axes[0].text(0.02, 0.98, "(a)", transform=axes[0].transAxes,
                 va="top", ha="left", fontsize=11)
    axes[0].legend(frameon=False, fontsize=8)
    idx = np.arange(len(rates))
    width = 0.38
    axes[1].bar(idx-width/2, rates, width, label="Floquet rates")
    axes[1].bar(idx+width/2, spectrum, width, label="QR spectrum")
    axes[1].axhline(0, color="0.25", lw=0.8)
    axes[1].set_xlabel("Mode index (sorted)")
    axes[1].set_ylabel("Rate (normalised time$^{-1}$)")
    axes[1].text(0.02, 0.98, "(b)", transform=axes[1].transAxes,
                 va="top", ha="left", fontsize=11)
    axes[1].legend(frameon=False, fontsize=8)
    stem = args.output_dir / "floquet_lyapunov_diagnostics"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    manifest = {"status": "PASS", "source": str(result_path), "outputs": [str(stem.with_suffix('.pdf')), str(stem.with_suffix('.png'))], "no_mock_data": True}
    (args.output_dir / "floquet_lyapunov_diagnostics.manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
