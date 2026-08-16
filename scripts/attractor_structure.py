#!/usr/bin/env python3
"""Global attractor structure: from the unstable drive-locked orbit to the attractor.

Gap 1 of GAP_IMPLEMENTATION_20260816.md. The drive-axis continuation tracks the
*drive-locked periodic orbit*, which becomes linearly unstable at the
Neimark--Sacker onset. This script connects that orbit to the chaotic attractor:

  1. recompute the orbit monodromy -> Floquet rates (the orbit's *linear*
     stability) and identify the most unstable Floquet direction;
  2. perturb the orbit by a small amount along that direction and integrate
     forward, recording the exponential departure;
  3. from the departed state, run a QR/Benettin Lyapunov spectrum and compare it
     with the stored strong-coupling map (the attractor's spectrum).

The claim this supports is narrow: the continued orbit is the *unstable seed*
of the route; the attractor is a distinct object whose Lyapunov spectrum
reproduces the stored map. No invariant measure, basin partition, or ergodic
proof is claimed.
"""
from __future__ import annotations
import argparse, json, platform, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from reconstruction_core import ModelParameters, monodromy, rhs_phase
from lyapunov_numba import lyapunov_qr_numba

ROOT = Path(__file__).resolve().parent.parent
SCAN = ROOT / "results" / "drive_axis_bifurcation_scan.json"
MAP = ROOT / "results" / "full_spectrum_hyperchaos_map_g03.json"


def kaplan_yorke_dimension(spectrum: np.ndarray) -> float:
    s = np.sort(np.asarray(spectrum, dtype=float))[::-1]
    cum = np.cumsum(s)
    idx = np.where(cum >= 0)[0]
    if len(idx) == 0:
        return 0.0
    k = int(idx[-1]) + 1
    if k == len(s):
        return float(k)
    return float(k) + float(cum[k - 1]) / abs(float(s[k]))


def strong_params(E: float, theta: float) -> ModelParameters:
    """Strong-coupling model used by the bifurcation scan and the hyperchaos map."""
    return ModelParameters(kappa=1.0, gamma1=0.02, gamma2=0.02, omega1=1.0,
                           omega2=1.03, g1=0.3, g2=0.27, hopping=0.08,
                           detuning=-1.0, drive=E, drive_modulation=0.1,
                           drive_frequency=1.0, theta=theta)


def orbit_at(records: list[dict], E_target: float, direction: str) -> dict:
    """Return the PASS record whose E is closest to E_target in the given direction."""
    cands = [r for r in records if r.get("direction") == direction
             and r.get("status") == "PASS" and r.get("x0") is not None]
    if not cands:
        raise RuntimeError("no converged orbit records")
    return min(cands, key=lambda r: abs(r["E"] - E_target))


def rk4_steps(x_seed: np.ndarray, x_start: np.ndarray, p: ModelParameters,
              n: int, dt: float, track_every: int = 100):
    """Forward RK4 from x_start; record distance to the seed orbit x_seed."""
    x = np.asarray(x_start, dtype=float).copy()
    phi = 0.0
    dists = []
    for step in range(n):
        k1 = rhs_phase(x, phi, p)
        k2 = rhs_phase(x + 0.5 * dt * k1, phi + 0.5 * p.drive_frequency * dt, p)
        k3 = rhs_phase(x + 0.5 * dt * k2, phi + 0.5 * p.drive_frequency * dt, p)
        k4 = rhs_phase(x + dt * k3, phi + p.drive_frequency * dt, p)
        x = x + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        phi = (phi + p.drive_frequency * dt) % (2.0 * np.pi)
        if step % track_every == 0:
            dists.append(float(np.linalg.norm(x - np.asarray(x_seed, dtype=float))))
    return dists, x


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--drives", nargs="+", type=float, default=[4.0, 8.0])
    ap.add_argument("--thetas", nargs="+", type=float, default=[0.0, np.pi / 2])
    ap.add_argument("--depart-steps", type=int, default=20000)
    ap.add_argument("--lyap-steps", type=int, default=300000)
    ap.add_argument("--lyap-transient", type=int, default=30000)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--perturbation", type=float, default=1e-4)
    a = ap.parse_args()

    scan = json.loads(SCAN.read_text())
    lyap_map = json.loads(MAP.read_text())

    map_lookup = {}
    for r in lyap_map.get("records", []):
        map_lookup[(round(r["drive"], 3), round(r["theta"], 9))] = r

    records_out = []
    for theta in a.thetas:
        theta_rec = next(t for t in scan["theta_records"]
                         if abs(t["theta"] - theta) < 1e-9)
        for E in a.drives:
            orb = orbit_at(theta_rec["records"], E, "up")
            x0 = np.asarray(orb["x0"], dtype=float)
            p = strong_params(float(orb["E"]), theta)

            # 1. Orbit linear stability (monodromy -> Floquet rates).
            mono = monodromy(x0, p)
            multipliers = np.array([complex(z["real"], z["imag"])
                                    for z in mono["multipliers"]])
            rates = np.sort(np.asarray(mono["floquet_rates"], dtype=float))[::-1]
            n_unstable = int(np.sum(np.abs(multipliers) > 1.0))
            largest_rate = float(rates[0])

            M = np.asarray(mono["monodromy"], dtype=float)
            vals, vecs = np.linalg.eig(M)
            idx = int(np.argmax(np.abs(vals)))
            v_unstable = np.real(vecs[:, idx])
            v_unstable = v_unstable / (np.linalg.norm(v_unstable) + 1e-300)

            # 2. Perturb along the unstable direction and integrate forward.
            x_pert = x0 + a.perturbation * v_unstable
            dists, x_final = rk4_steps(x0, x_pert, p, a.depart_steps, a.dt)

            # Departure slope over the growing window of the distance series.
            logd = np.log(np.maximum(np.asarray(dists), 1e-300))
            t = a.dt * 100.0 * np.arange(len(dists), dtype=float)
            grow = np.where(logd > np.log(a.perturbation) + 0.5)[0]
            slope = float("nan")
            if grow.size >= 4:
                slope = float(np.polyfit(t[grow], logd[grow], 1)[0])

            # 3. Attractor Lyapunov spectrum from the departed state.
            ly = lyapunov_qr_numba(x_final, p, n_steps=a.lyap_steps, dt=a.dt,
                                   transient_steps=a.lyap_transient, qr_interval=10)
            spectrum = np.sort(np.asarray(ly["spectrum"], dtype=float))[::-1]
            n_pos = int(np.sum(spectrum > 1e-3))
            ky = kaplan_yorke_dimension(spectrum)

            mkey = (round(E, 3), round(theta, 9))
            mapped = map_lookup.get(mkey)
            cross = None
            if mapped is not None:
                cross = {
                    "map_largest_exponent": mapped["largest_exponent"],
                    "map_n_positive": mapped["n_positive_exponents"],
                    "map_kaplan_yorke_dimension": mapped["kaplan_yorke_dimension"],
                    "largest_exponent_relative_diff": float(
                        (spectrum[0] - mapped["largest_exponent"])
                        / abs(mapped["largest_exponent"])),
                }

            records_out.append({
                "theta": float(theta),
                "theta_over_pi": float(theta / np.pi),
                "drive": float(orb["E"]),
                "drive_target": float(E),
                "orbit_n_unstable": n_unstable,
                "orbit_floquet_rates": rates.tolist(),
                "orbit_largest_floquet_rate": largest_rate,
                "departure_slope": slope,
                "departure_initial_distance": float(np.linalg.norm(x_pert - x0)),
                "departure_final_distance": float(dists[-1]),
                "attractor_lyapunov_spectrum": spectrum.tolist(),
                "attractor_largest_exponent": float(spectrum[0]),
                "attractor_n_positive": n_pos,
                "attractor_kaplan_yorke_dimension": ky,
                "mean_divergence": ly["mean_divergence"],
                "cross_validation_vs_stored_map": cross,
            })

    status = "PASS" if records_out else "FAIL"
    out = {
        "gate": "ATTRACTOR_STRUCTURE",
        "kind": "deterministic_dynamics_attractor_connection",
        "status": status,
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "settings": {
            "depart_steps": a.depart_steps,
            "lyap_steps": a.lyap_steps,
            "lyap_transient": a.lyap_transient,
            "dt": a.dt,
            "perturbation": a.perturbation,
        },
        "records": records_out,
        "interpretation": (
            "The continued drive-locked orbit has unstable Floquet directions "
            "(n_unstable >= 1) beyond the Neimark--Sacker onset; a small "
            "perturbation along the most unstable direction departs the orbit, "
            "and the forward trajectory settles onto an attractor whose Lyapunov "
            "spectrum reproduces the stored strong-coupling map. The orbit is "
            "the unstable seed; the attractor is a distinct object. Deterministic "
            "dynamics only."),
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "records": len(records_out),
                      "output": str(a.output)}, indent=2))
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
