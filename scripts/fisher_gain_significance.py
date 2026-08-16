#!/usr/bin/env python3
"""Statistical significance of the matched force-sensing Fisher gains.

A severe-referee item on the sensing result: the abstract-level "gain 1.323"
is a single deterministic point at the maximum-gain phase, while the
uncertainty-quantified distribution (results/matched_fisher_uq.json) has a
much lower median and a 5th percentile barely above unity. This script
bootstraps the gain distribution to put confidence intervals on the median,
the 5th percentile, and the fraction of the sampled parameter range above
unity, and reports the one-sided test of whether the gain is distinguishable
from 1.0 at the conservative (5th-percentile) end.

Pure post-processing of matched_fisher_uq.json; no new simulations.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
UQ = ROOT / "results" / "matched_fisher_uq.json"


def bootstrap_ci(x: np.ndarray, stat, *, n_boot: int, seed: int, alpha: float):
    rng = np.random.default_rng(seed)
    vals = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, x.size, size=x.size)
        vals[b] = stat(x[idx])
    lo = float(np.percentile(vals, 100 * alpha / 2))
    hi = float(np.percentile(vals, 100 * (1 - alpha / 2)))
    return {"point": float(stat(x)), "ci_low": lo, "ci_high": hi,
            "alpha": alpha, "bootstrap_std": float(np.std(vals))}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--n-boot", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--alpha", type=float, default=0.10)
    a = ap.parse_args()

    data = json.loads(UQ.read_text())
    baseline = data["baseline_deterministic"]

    def stats_key(name: str):
        g = np.array([r[name] for r in data["records"]
                      if r.get("status") == "PASS"])
        g = g[np.isfinite(g)]
        med = bootstrap_ci(g, lambda x: float(np.median(x)),
                           n_boot=a.n_boot, seed=a.seed, alpha=a.alpha)
        p05 = bootstrap_ci(g, lambda x: float(np.percentile(x, 5)),
                           n_boot=a.n_boot, seed=a.seed + 1, alpha=a.alpha)
        p95 = bootstrap_ci(g, lambda x: float(np.percentile(x, 95)),
                           n_boot=a.n_boot, seed=a.seed + 2, alpha=a.alpha)
        frac = bootstrap_ci(g, lambda x: float(np.mean(x > 1.0)),
                            n_boot=a.n_boot, seed=a.seed + 3, alpha=a.alpha)
        # One-sided: is the conservative (5th percentile) end above unity?
        p05_above_unity = bool(p05["ci_low"] > 1.0)
        # Excess of the median over unity with its CI.
        med_excess = {
            "point": med["point"] - 1.0,
            "ci_low": med["ci_low"] - 1.0,
            "ci_high": med["ci_high"] - 1.0,
        }
        return {
            "n": int(g.size),
            "deterministic_baseline": baseline[name],
            "baseline_percentile_rank": float(np.mean(g <= baseline[name]) * 100),
            "median": med,
            "median_excess_over_unity": med_excess,
            "p05": p05,
            "p95": p95,
            "fraction_above_unity": frac,
            "p05_ci_low_above_unity": p05_above_unity,
            "gain_distinguishable_from_unity_at_5pct": p05_above_unity,
        }

    out = {
        "gate": "FISHER_GAIN_SIGNIFICANCE",
        "kind": "bootstrap_confidence_intervals_classical_fisher_gain",
        "status": "PASS",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "settings": {"n_boot": a.n_boot, "seed": a.seed, "alpha": a.alpha,
                     "source": str(UQ.relative_to(ROOT))},
        "gain_vs_flux_off": stats_key("gain_vs_flux_off"),
        "gain_vs_single_mode": stats_key("gain_vs_single_mode"),
        "interpretation": (
            "The deterministic 1.323 (flux-off) / 1.159 (single-mode) gain is "
            "the value at the maximum-gain phase; across the sampled unmeasured "
            "coordinates the gain distribution is much lower (median ~1.08/1.04, "
            "5th percentile ~1.004/1.002). The bootstrap places a confidence "
            "interval on these statistics. The gain is 'distinguishable from "
            "unity at the conservative end' only if the 5th-percentile CI lower "
            "bound exceeds 1.0; otherwise the enhancement is real but weak and "
            "not robustly separable from unity at the conservative end. "
            "Classical Fisher information, not QFI."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
