#!/usr/bin/env python3
"""Fail-closed robust Pareto optimizer for the calibrated FL_QOM model.

The optimizer refuses to run until the physical manifest contains an explicit SI
calibration, normalized ranges, uncertainties, objective constraints, and a declared
primary objective. It therefore cannot silently optimize the current pilot values or
legacy candidate values.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, find_periodic_orbit, lyapunov_qr, monodromy
from covariance_fisher import periodic_covariance

REQUIRED_PARAMETERS = (
    "kappa", "gamma1", "gamma2", "omega1", "omega2", "g1", "g2",
    "hopping", "detuning", "drive", "drive_modulation", "drive_frequency",
)
MAXIMIZE = ("stability_margin", "classical_fisher_information", "robustness_margin")
MINIMIZE = ("resource_cost",)


def calibration_errors(manifest: dict, *, allow_provisional: bool = False) -> list[str]:
    """Return blocking errors without silently promoting exploratory values.

    A normal run requires a complete SI calibration.  The explicit provisional mode
    is only for a literature-anchored numerical screening and requires a source,
    normalized ranges, and constraints; its output is labelled PROVISIONAL and cannot
    support an experimental-optimum claim.
    """
    errors: list[str] = []
    status = manifest.get("status")
    if status != "PASS":
        if not (allow_provisional and status == "PROVISIONAL"):
            errors.append(f"manifest status is {status!r}, expected 'PASS'")
    normalization = manifest.get("normalization", {})
    if normalization.get("reference_frequency_si_hz") is None:
        errors.append("normalization.reference_frequency_si_hz is missing")
    if normalization.get("reference_frequency_source") is None:
        errors.append("normalization.reference_frequency_source is missing")
    if allow_provisional and not manifest.get("source_reference", {}).get("url"):
        errors.append("provisional source_reference.url is missing")
    objective = manifest.get("objective", {})
    if objective.get("selected_primary_objective") != "multiobjective_pareto_compromise":
        errors.append("primary objective is not the selected robust Pareto objective")
    if not objective.get("constraints"):
        errors.append("objective.constraints is empty")
    entries = {item.get("name"): item for item in manifest.get("parameters", [])}
    for name in REQUIRED_PARAMETERS:
        item = entries.get(name)
        if item is None:
            errors.append(f"parameter {name} is absent")
            continue
        fields = ("normalized_range",)
        if not allow_provisional:
            fields = ("si_value", "si_units", "admissible_si_range", "experimental_source", "uncertainty", "normalized_range")
        for field in fields:
            if item.get(field) is None:
                errors.append(f"{name}.{field} is missing")
    return errors


def lhs(n: int, d: int, rng: np.random.Generator) -> np.ndarray:
    """Generate a reproducible Latin-hypercube matrix in [0, 1]."""
    result = np.empty((n, d), dtype=float)
    for j in range(d):
        result[:, j] = (rng.permutation(n) + rng.random(n)) / n
    return result


def build_parameters(values: dict[str, float], theta: float) -> ModelParameters:
    return ModelParameters(theta=theta, **{name: float(values[name]) for name in REQUIRED_PARAMETERS})


def fisher_at_theta(p: ModelParameters, theta: float, eps: float, n_th1: float, n_th2: float, detector_variance: float,
                     orbit_max_iter: int, orbit_residual_tol: float) -> float:
    def observation(phase: float) -> dict:
        phase_params = replace(p, theta=float(phase))
        orbit = find_periodic_orbit(np.zeros(6), phase_params, max_iter=orbit_max_iter,
                                     residual_tol=orbit_residual_tol)
        if orbit.get("status") != "PASS":
            raise RuntimeError(orbit.get("reason", "periodic orbit failed"))
        return periodic_covariance(
            np.asarray(orbit["x0"]), phase_params,
            n_th1=n_th1, n_th2=n_th2, detector_variance=detector_variance,
        )
    plus = observation(theta + eps)
    minus = observation(theta - eps)
    center = observation(theta)
    dm = (plus["mean_signal"] - minus["mean_signal"]) / (2 * eps)
    dv = (plus["mean_variance"] - minus["mean_variance"]) / (2 * eps)
    variance = center["mean_variance"]
    value = float(dm * dm / variance + 0.5 * dv * dv / (variance * variance))
    if not np.isfinite(value):
        raise RuntimeError("classical Fisher information is non-finite")
    return value


def evaluate_candidate(values: dict[str, float], args: argparse.Namespace, rng: np.random.Generator) -> dict:
    p = build_parameters(values, args.theta)
    perturbations = [np.zeros(len(REQUIRED_PARAMETERS))]
    if args.robust_replicates > 1:
        perturbations.extend(rng.normal(0.0, args.uncertainty_fraction, (args.robust_replicates - 1, len(REQUIRED_PARAMETERS))))
    margins: list[float] = []
    center_margins: list[float] = []
    lyap_max: list[float] = []
    fisher_values: list[float] = []
    initial_conditions = [np.zeros(6)]
    if args.initial_condition_replicates > 1:
        initial_conditions.extend(
            rng.normal(0.0, args.initial_condition_scale,
                       (args.initial_condition_replicates - 1, 6))
        )
    for perturbation_index, perturbation in enumerate(perturbations):
        perturbed = {
            name: values[name] * (1.0 + float(perturbation[i]))
            for i, name in enumerate(REQUIRED_PARAMETERS)
        }
        pp = build_parameters(perturbed, args.theta)
        for x_initial in initial_conditions:
            orbit = find_periodic_orbit(np.asarray(x_initial), pp,
                                         max_iter=args.orbit_max_iter,
                                         residual_tol=args.orbit_residual_tol)
            if orbit.get("status") != "PASS":
                return {"status": "FAIL", "reason": "periodic orbit failed under robustness replicate or initial condition"}
            floquet = monodromy(np.asarray(orbit["x0"]), pp)
            floquet_max = float(max(np.real(floquet["floquet_rates"])))
            margin = -floquet_max
            margins.append(margin)
            if perturbation_index == 0:
                center_margins.append(margin)
            if args.lyapunov_steps > 0:
                lyap = lyapunov_qr(np.asarray(orbit["x0"]), pp, n_steps=args.lyapunov_steps,
                                   transient_steps=args.lyapunov_transient, dt=args.dt, qr_interval=args.qr_interval)
                lyap_max.append(float(max(lyap["spectrum"])))
        if args.compute_fisher:
            fisher_values.append(fisher_at_theta(pp, args.theta, args.fisher_eps, args.n_th1, args.n_th2,
                                                  args.detector_variance, args.orbit_max_iter,
                                                  args.orbit_residual_tol))
    metrics = {
        "stability_margin": float(min(center_margins)),
        "robustness_margin": float(min(margins)),
        "floquet_max_rate": float(-min(margins)),
        "resource_cost": float(abs(values["drive"])),
    }
    if lyap_max:
        metrics["lyapunov_max"] = float(max(lyap_max))
    if fisher_values:
        metrics["classical_fisher_information"] = float(min(fisher_values))
    else:
        metrics["classical_fisher_information"] = None
    return {"status": "PASS", "normalized_parameters": values, "metrics": metrics, "robust_replicates": len(perturbations)}


def satisfies_constraints(metrics: dict, constraints: list[dict]) -> bool:
    for constraint in constraints:
        metric = constraint.get("metric")
        operator = constraint.get("operator")
        limit = constraint.get("value")
        value = metrics.get(metric)
        if value is None or limit is None:
            return False
        if operator == "<=" and not value <= limit:
            return False
        if operator == ">=" and not value >= limit:
            return False
        if operator == "<" and not value < limit:
            return False
        if operator == ">" and not value > limit:
            return False
        if operator == "==" and not np.isclose(value, limit):
            return False
    return True


def dominates(a: dict, b: dict) -> bool:
    av, bv = a["metrics"], b["metrics"]
    a_values = [av["stability_margin"], av["classical_fisher_information"], av["robustness_margin"], -av["resource_cost"]]
    b_values = [bv["stability_margin"], bv["classical_fisher_information"], bv["robustness_margin"], -bv["resource_cost"]]
    if any(value is None for value in a_values + b_values):
        return False
    return all(x >= y for x, y in zip(a_values, b_values)) and any(x > y for x, y in zip(a_values, b_values))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-provisional-literature-scenario", action="store_true",
                        help="Run an explicitly labelled literature-anchored screening; never an experimental calibration.")
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260813)
    parser.add_argument("--theta", type=float, default=0.0)
    parser.add_argument("--robust-replicates", type=int, default=3)
    parser.add_argument("--orbit-max-iter", type=int, default=600,
                        help="Maximum Poincare fixed-point iterations; low-damping literature regimes may need more.")
    parser.add_argument("--orbit-residual-tol", type=float, default=1e-8,
                        help="Periodic-orbit residual gate; provisional screening may use a looser value, but production must tighten it.")
    parser.add_argument("--uncertainty-fraction", type=float, default=0.02)
    parser.add_argument("--initial-condition-replicates", type=int, default=1,
                        help="Independent initial conditions per parameter replicate; production robustness should use at least 3.")
    parser.add_argument("--initial-condition-scale", type=float, default=0.1,
                        help="Dimensionless Gaussian scale for additional initial conditions.")
    parser.add_argument("--compute-fisher", action="store_true")
    parser.add_argument("--fisher-eps", type=float, default=1e-4)
    parser.add_argument("--n-th1", type=float, default=0.1)
    parser.add_argument("--n-th2", type=float, default=0.1)
    parser.add_argument("--detector-variance", type=float, default=0.01)
    parser.add_argument("--lyapunov-steps", type=int, default=0)
    parser.add_argument("--lyapunov-transient", type=int, default=0)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--qr-interval", type=int, default=10)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    errors = calibration_errors(manifest, allow_provisional=args.allow_provisional_literature_scenario)
    objective_names = set(manifest.get("objective", {}).get("secondary_objectives", []))
    if "classical_measurement_fisher_information" in objective_names and not args.compute_fisher:
        errors.append("classical Fisher objective is declared but --compute-fisher was not supplied")
    if errors:
        output = {
            "status": "NOT_COMPUTED",
            "kind": "robust_pareto_optimization",
            "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version,
            "platform": platform.platform(),
            "manifest": str(args.manifest.resolve()),
            "errors": errors,
            "interpretation": "No optimization was run because physical calibration, ranges, uncertainties, or constraints are incomplete.",
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(output, indent=2, sort_keys=True))
        return 1

    entries = {item["name"]: item for item in manifest["parameters"]}
    ranges = np.array([entries[name]["normalized_range"] for name in REQUIRED_PARAMETERS], dtype=float)
    if np.any(~np.isfinite(ranges)) or np.any(ranges[:, 1] <= ranges[:, 0]):
        raise SystemExit("invalid normalized ranges")
    rng = np.random.default_rng(args.seed)
    unit_design = lhs(args.samples, len(REQUIRED_PARAMETERS), rng)
    candidates = []
    for row in unit_design:
        values = {name: float(ranges[i, 0] + row[i] * (ranges[i, 1] - ranges[i, 0])) for i, name in enumerate(REQUIRED_PARAMETERS)}
        candidate = evaluate_candidate(values, args, rng)
        if candidate.get("status") == "PASS":
            candidate["feasible"] = satisfies_constraints(candidate["metrics"], manifest["objective"]["constraints"])
        candidates.append(candidate)
    feasible = [candidate for candidate in candidates if candidate.get("feasible")]
    pareto = [candidate for candidate in feasible if not any(dominates(other, candidate) for other in feasible)]
    output = {
        "status": ("PROVISIONAL" if manifest.get("status") == "PROVISIONAL" else ("PASS" if pareto else "FAIL")),
        "scientific_status": (
            "PROVISIONAL_PARETO_FOUND" if pareto else "NO_ROBUST_FEASIBLE_CANDIDATE"
        ) if manifest.get("status") == "PROVISIONAL" else ("PASS" if pareto else "NO_FEASIBLE_CANDIDATE"),
        "kind": "robust_pareto_optimization",
        "calibration_status": manifest.get("status"),
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "manifest": str(args.manifest.resolve()),
        "seed": args.seed,
        "samples": args.samples,
        "robust_replicates": args.robust_replicates,
        "initial_condition_replicates": args.initial_condition_replicates,
        "initial_condition_scale": args.initial_condition_scale,
        "orbit_max_iter": args.orbit_max_iter,
        "orbit_residual_tol": args.orbit_residual_tol,
        "objective": manifest["objective"],
        "candidate_records": candidates,
        "candidate_count": len(candidates),
        "feasible_count": len(feasible),
        "pareto_count": len(pareto),
        "pareto_candidates": pareto,
        "interpretation": (
            "Provisional literature-anchored numerical Pareto screening; not an experimental optimum."
            if manifest.get("status") == "PROVISIONAL"
            else "Pareto set under the declared calibrated normalized ranges and constraints."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"], "feasible": len(feasible), "pareto": len(pareto), "output": str(args.output)}, indent=2))
    return 0 if output["status"] in {"PASS", "PROVISIONAL"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
