#!/usr/bin/env python3
"""Audit cross-check: verify every quantitative claim in the Track D manuscript
against the immutable result JSONs (2026-08-16 audit).

Usage:
  mamba run -n qom python scripts/audit_crosscheck_20260816.py

Exit code 0 if every checked claim passes; 1 otherwise.
"""
from __future__ import annotations
import json
from pathlib import Path

R = Path(__file__).resolve().parent.parent / "results"


def load(name: str):
    with open(R / name) as f:
        return json.load(f)


fails = 0


def check(label: str, ok: bool, detail: str = ""):
    global fails
    if not ok:
        fails += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f"  ({detail})" if detail else ""))


# ---- model_gate.json ----
print("== model_gate.json")
mg = load("model_gate.json")
check("status", mg["status"] == "PASS")
check("jacobian 2.92e-10", abs(mg["jacobian_max_abs_error"] - 2.92e-10) < 0.01e-10, str(mg["jacobian_max_abs_error"]))
check("trace -1.04", abs(mg["trace"] - (-1.04)) < 1e-12)
check("flux residual < 3e-16", mg["flux_phase_residual"] < 3e-16, str(mg["flux_phase_residual"]))

# ---- floquet_lyapunov.json ----
print("== floquet_lyapunov.json")
fl = load("floquet_lyapunov.json")
po = fl.get("periodic_orbit", {})
check("orbit residual 8.49e-9 (<9e-9)", abs(po.get("residual", 0) - 8.485457714301273e-9) < 1e-12, str(po.get("residual")))
c = fl["comparison"]
check("max rate diff 2.68e-4", abs(c["max_rate_spectrum_abs_difference"] - 2.678e-4) < 0.01e-4, str(c["max_rate_spectrum_abs_difference"]))
check("tail block std 5.3e-9", abs(c["tail_block_max_std"] - 5.294459447847316e-9) < 0.1e-9)
check("mean divergence -1.04", abs(c["mean_divergence"] - (-1.04)) < 1e-3)
check("divergence residual 4.75e-5", abs(c["divergence_residual"] - 4.745978429410158e-5) < 0.01e-5)

# ---- flux_grid.json ----
print("== flux_grid.json")
fg = load("flux_grid.json")
recs = fg["records"]
check("27 records", len(recs) == 27, f"n={len(recs)}")
check("9 phases x 3 seeds", len({r["theta"] for r in recs}) == 9 and len({r["seed"] for r in recs}) == 27)
check("all status PASS", {r["status"] for r in recs} == {"PASS"})
check("orbit status all PASS", {r["orbit_status"] for r in recs} == {"PASS"})
check("no chaos", not any(r["is_chaotic"] for r in recs))
check("largest QR exponent -0.01002", abs(max(r["lyapunov_spectrum"][-1] for r in recs) - (-0.010022961746636227)) < 1e-6,
      str(max(r["lyapunov_spectrum"][-1] for r in recs)))
check("max Floquet-QR diff 2.849e-4", abs(max(r["max_rate_difference"] for r in recs) - 2.8493199670920344e-4) < 1e-9,
      str(max(r["max_rate_difference"] for r in recs)))
check("max block std 6.10e-9", abs(max(r["tail_block_std"] for r in recs) - 6.102306531569741e-9) < 1e-12)
check("max divergence residual 4.75e-5", abs(max(r["divergence_residual"] for r in recs) - 4.74597625501616e-5) < 1e-9)

# ---- full_spectrum_transition_map.json ----
print("== full_spectrum_transition_map.json")
ft = load("full_spectrum_transition_map.json")
recs = ft["records"]
check("54 records", len(recs) == 54, f"n={len(recs)}")
check("54/54 PASS", ft["passed_records"] == 54 == ft["total_records"])
check("largest exponent in [-0.043, -0.010]",
      min(r["largest_exponent"] for r in recs) >= -0.0435 and max(r["largest_exponent"] for r in recs) <= -0.0100,
      f"[{min(r['largest_exponent'] for r in recs)}, {max(r['largest_exponent'] for r in recs)}]")
check("zero positive exponents", all(r["n_positive_exponents"] == 0 for r in recs))
check("KY dim 0 everywhere", all(r["kaplan_yorke_dimension"] == 0.0 for r in recs))

# ---- full_spectrum_hyperchaos_map_g03.json ----
print("== full_spectrum_hyperchaos_map_g03.json")
fh = load("full_spectrum_hyperchaos_map_g03.json")
recs = fh["records"]
check("60 records, 60/60 PASS", len(recs) == 60 and fh["passed_records"] == 60 == fh["total_records"])
check("coupling g1=0.3 g2=0.27", fh["coupling"]["g1"] == 0.3 and fh["coupling"]["g2"] == 0.27)
check("12 phases x 5 drives", len({r["theta"] for r in recs}) == 12 and len({r["drive"] for r in recs}) == 5)
by_drive = {}
for r in recs:
    by_drive.setdefault(r["drive"], set()).add(r["n_positive_exponents"])
check("E=0.2: n_pos=0", by_drive[0.2] == {0}, str(by_drive[0.2]))
check("E=1.0: n_pos<=1", by_drive[1.0] <= {0, 1}, str(by_drive[1.0]))
check("E=4.0: n_pos up to 3", max(by_drive[4.0]) == 3, str(by_drive[4.0]))
check("E=8.0: n_pos 2..4", by_drive[8.0] == {2, 3, 4}, str(by_drive[8.0]))
check("max largest exponent 0.3223", abs(max(r["largest_exponent"] for r in recs) - 0.32226279689568815) < 1e-4,
      str(max(r["largest_exponent"] for r in recs)))
check("max KY dim 4.802", abs(max(r["kaplan_yorke_dimension"] for r in recs) - 4.801606545905644) < 1e-3,
      str(max(r["kaplan_yorke_dimension"] for r in recs)))
check("flux == theta (patched provenance)", all(abs((r.get("flux") or 0) - r["theta"]) < 1e-9 for r in recs))

# ---- hyperchaos_seed_robustness.json ----
print("== hyperchaos_seed_robustness.json")
hs = load("hyperchaos_seed_robustness.json")
recs = hs["records"]
check("72 records, 72/72 PASS", len(recs) == 72 and hs["passed_records"] == 72)
e8 = [r for r in recs if r["drive"] == 8.0]
per_phase = {}
for r in e8:
    per_phase.setdefault(r["theta"], set()).add(r["n_positive_exponents"])
check("E=8 n_pos 2..4", set().union(*per_phase.values()) == {2, 3, 4})
check("E=8 spread <= 1 per phase", all(len(v) <= 2 and max(v) - min(v) <= 1 for v in per_phase.values()))
check("E=8 maxima at +-pi/2, minima at 0,+-pi",
      max(per_phase[-1.5707963267948966]) >= 3 and max(per_phase[1.5707963267948966]) >= 3
      and max(per_phase[0.0]) == 2 and max(per_phase[-3.141592653589793]) == 2)

# ---- hyperchaos_convergence_probe.json ----
print("== hyperchaos_convergence_probe.json")
hc = load("hyperchaos_convergence_probe.json")
recs = hc["records"]
expected = [
    (4.0, -0.5, 250000, 3, 0.132142, 1.01e-3, 6.69e-3),
    (4.0, -0.5, 500000, 2, 0.127291, 1.03e-3, 9.89e-3),
    (4.0, -0.5, 1000000, 2, 0.123729, 1.03e-3, 1.05e-2),
    (8.0, 0.5, 250000, 4, 0.324142, 1.90e-3, 3.30e-3),
    (8.0, 0.5, 500000, 4, 0.315008, 1.92e-3, 5.88e-3),
    (8.0, 0.5, 1000000, 4, 0.317277, 1.84e-3, 3.90e-3),
]
for r, (E, tp, steps, np_, le, dr, s) in zip(recs, expected):
    s_max = max(r["tail_block_std"])
    ok = (r["drive"] == E and r["theta_over_pi"] == tp and r["n_steps"] == steps
          and r["n_positive_exponents"] == np_
          and abs(r["largest_exponent"] - le) < 1e-6
          and abs(r["divergence_residual"] - dr) < 1e-5
          and abs(s_max - s) < 1e-4)
    check(f"probe (E={E}, theta/pi={tp}, steps={steps})", ok)
check("classification stable (longest two windows)", hc["classification_stability"]["E=4,theta/pi=-0.5"]["stable_between_longest_two"]
      and hc["classification_stability"]["E=8,theta/pi=0.5"]["stable_between_longest_two"])

# ---- covariance_fisher.json ----
print("== covariance_fisher.json")
cf = load("covariance_fisher.json")
recs = cf["records"]
mines = [r["center"]["min_quantum_physicality_eigenvalue"] for r in recs]
fishers = [r["classical_fisher_information"] for r in recs]
check("min phys eig 2.42e-7..2.58e-7", min(mines) >= 2.42e-7 and max(mines) <= 2.59e-7,
      f"[{min(mines)}, {max(mines)}]")
check("classical Fisher 3.78e-26..9.17e-13", min(fishers) >= 3.77e-26 and max(fishers) <= 9.18e-13,
      f"[{min(fishers)}, {max(fishers)}]")

# ---- matched_fisher_reference.json ----
print("== matched_fisher_reference.json")
mfr = load("matched_fisher_reference.json")
recs = mfr["records"]
check("5 flux records PASS", len(recs) == 5 and all(r["status"] == "PASS" for r in recs))
fc0 = recs[0]["classical_fisher_information"]
fcp = recs[-1]["classical_fisher_information"]
check("F_C 4.50e-5 (theta=0) -> 5.96e-5 (theta=pi)",
      abs(fc0 - 4.503008099534932e-5) < 1e-9 and abs(fcp - 5.9581655778556957e-5) < 1e-9,
      f"[{fc0}, {fcp}]")
gains = {round(g["theta"], 3): (g["gain_vs_flux_off"], g["gain_vs_single_mode"]) for g in mfr["gains"]}
check("gain vs flux-off 1.323 at theta=pi", abs(gains[3.142][0] - 1.3231523120002942) < 1e-9)
check("gain vs single-mode 1.159 at theta=pi", abs(gains[3.142][1] - 1.1590095550739294) < 1e-9)
check("references present", "flux_off_theta_0_hopping_0.08" in mfr["references"] and "single_mode_hopping_0" in mfr["references"])
eps_check = mfr.get("eps_convergence_check_theta_pi_over_2", [])
print(f"  INFO: eps_convergence_check_theta_pi_over_2 = {eps_check} (empty => SI Table tab:si-fisher-eps eps=1e-3/1e-5 rows not stored)")

# ---- matched_fisher_sweep.json ----
print("== matched_fisher_sweep.json")
mfs = load("matched_fisher_sweep.json")
check("16/16 sweep PASS", len(mfs["records"]) == 16 and all(r["status"] == "PASS" for r in mfs["records"]))
check("gain_vs_flux_off 1.323 std<1e-3",
      abs(mfs["gain_vs_flux_off_stats"]["mean"] - 1.323) < 1e-3 and mfs["gain_vs_flux_off_stats"]["std"] < 1e-3)
check("gain_vs_single_mode 1.159 std<1e-3",
      abs(mfs["gain_vs_single_mode_stats"]["mean"] - 1.159) < 1e-3 and mfs["gain_vs_single_mode_stats"]["std"] < 1e-3)

# ---- reduced_full_common_regime_comparison.json ----
print("== reduced_full_common_regime_comparison.json")
rc = load("reduced_full_common_regime_comparison.json")
by_k = {c["kappa"]: c for c in rc["cases"]}
check("kappa=20 diff 9.27e-7", abs(by_k[20.0]["mechanical_inf_error"] - 9.270964667811068e-7) < 1e-12)
check("kappa=100 diff 2.08e-8", abs(by_k[100.0]["mechanical_inf_error"] - 2.0768420381518217e-8) < 1e-12)

# ---- Pareto v4 ----
print("== physical_parameter_reduced_bad_cavity_pareto_v4.json")
p4 = load("physical_parameter_reduced_bad_cavity_pareto_v4.json")
pc = sorted(p4["pareto_candidates"], key=lambda x: x["drive_cost"])
costs = [c["drive_cost"] for c in pc]
margins = [c["minimum_stability_margin"] for c in pc]
check("3 non-dominated candidates", len(pc) == 3)
check("costs 0.0800/0.1312/0.2203", all(abs(a - b) < 1e-6 for a, b in zip(costs, [0.08002407081169909, 0.13120321582088784, 0.22025850962654964])), str(costs))
check("margins 0.0034526/0.0035354/0.0035597", all(abs(a - b) < 1e-7 for a, b in zip(margins, [0.0034525598353587483, 0.0035354157489547733, 0.0035596603297002316])), str(margins))

# ---- screen v4 (radius claim) ----
print("== physical_parameter_reduced_bad_cavity_screen_robust_v4.json")
s4 = load("physical_parameter_reduced_bad_cavity_screen_robust_v4.json")
check("16/16 candidates feasible", s4["candidate_count"] == 16 and s4["feasible_count"] == 16)
radii = [rr["poincare_spectral_radius"] for cr in s4["candidate_records"] for rr in cr["replicate_records"]]
check("max radius 0.9970026 -> margin 0.0029974", abs(max(radii) - 0.9970026418196289) < 1e-9,
      str(max(radii)))
check("5 replicas x 5 starts per candidate", all(len(cr["replicate_records"]) == 25 for cr in s4["candidate_records"]))

# ---- fast_basin_capture_validation_v2.json ----
print("== fast_basin_capture_validation_v2.json")
fb = load("fast_basin_capture_validation_v2.json")
ics = [ic for cr in fb["candidate_results"] for rr in cr["replicate_results"] for ic in rr["capture"]["initial_condition_records"]]
check("75/75 joint-gate passes", len(ics) == 75 and all(ic["status"] == "PASS" for ic in ics), f"n={len(ics)}")

print()
print(f"RESULT: {'ALL CHECKS PASS' if fails == 0 else f'{fails} CHECK(S) FAILED'}")
raise SystemExit(1 if fails else 0)
