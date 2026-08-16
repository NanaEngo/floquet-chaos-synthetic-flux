"""Numba backend for the deterministic QR/Benettin integration.

The formulas intentionally mirror reconstruction_core.py.  This module contains
no model changes; it only removes Python interpreter overhead from the fixed-step
integration loop used by long validation windows.
"""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _rhs(x, phi, p):
    ar, ai, b1r, b1i, b2r, b2i = x
    c = np.cos(p[7])
    s = np.sin(p[7])
    omega_eff = p[8] + 2.0 * p[5] * b1r + 2.0 * p[6] * b2r
    e = p[10] * (1.0 + p[11] * np.cos(phi))
    out = np.empty(6, dtype=np.float64)
    out[0] = -0.5 * p[0] * ar - omega_eff * ai + e
    out[1] = omega_eff * ar - 0.5 * p[0] * ai
    out[2] = -0.5 * p[1] * b1r + p[3] * b1i + p[9] * (s * b2r + c * b2i)
    out[3] = p[5] * (ar * ar + ai * ai) - p[3] * b1r - 0.5 * p[1] * b1i - p[9] * c * b2r + p[9] * s * b2i
    out[4] = -0.5 * p[2] * b2r + p[4] * b2i + p[9] * (-s * b1r + c * b1i)
    out[5] = p[6] * (ar * ar + ai * ai) - p[4] * b2r - 0.5 * p[2] * b2i - p[9] * c * b1r - p[9] * s * b1i
    return out


@njit(cache=True)
def _jac(x, p):
    ar, ai, b1r, _b1i, b2r, _b2i = x
    c = np.cos(p[7])
    s = np.sin(p[7])
    omega_eff = p[8] + 2.0 * p[5] * b1r + 2.0 * p[6] * b2r
    out = np.empty((6, 6), dtype=np.float64)
    out[0, :] = np.array([-0.5*p[0], -omega_eff, -2*p[5]*ai, 0.0, -2*p[6]*ai, 0.0])
    out[1, :] = np.array([omega_eff, -0.5*p[0], 2*p[5]*ar, 0.0, 2*p[6]*ar, 0.0])
    out[2, :] = np.array([0.0, 0.0, -0.5*p[1], p[3], p[9]*s, p[9]*c])
    out[3, :] = np.array([2*p[5]*ar, 2*p[5]*ai, -p[3], -0.5*p[1], -p[9]*c, p[9]*s])
    out[4, :] = np.array([0.0, 0.0, -p[9]*s, p[9]*c, -0.5*p[2], p[4]])
    out[5, :] = np.array([2*p[6]*ar, 2*p[6]*ai, -p[9]*c, -p[9]*s, -p[4], -0.5*p[2]])
    return out


@njit(cache=True)
def _integrate(x0, p, n_steps, dt, transient_steps, qr_interval):
    x = x0.copy()
    phi = 0.0
    qmat = np.eye(6)
    accum = np.zeros(6)
    n_blocks = (n_steps - transient_steps) // qr_interval
    blocks = np.empty((n_blocks, 6), dtype=np.float64)
    block_count = 0
    div_accum = 0.0
    total = 0.0
    for step in range(n_steps):
        k1 = _rhs(x, phi, p)
        k2 = _rhs(x + 0.5 * dt * k1, phi + 0.5 * p[12] * dt, p)
        k3 = _rhs(x + 0.5 * dt * k2, phi + 0.5 * p[12] * dt, p)
        k4 = _rhs(x + dt * k3, phi + p[12] * dt, p)
        x_new = x + dt * (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0
        amat = _jac(x, p)
        amat_new = _jac(x_new, p)
        qtrial = qmat + 0.5 * dt * (amat @ qmat + amat_new @ (qmat + dt * amat @ qmat))
        x = x_new
        phi = (phi + p[12] * dt) % (2.0 * np.pi)
        if step < transient_steps:
            qmat = np.eye(6)
            continue
        div_accum += np.trace(amat) * dt
        if (step - transient_steps + 1) % qr_interval == 0:
            qmat, rmat = np.linalg.qr(qtrial)
            for j in range(6):
                if rmat[j, j] < 0.0:
                    qmat[:, j] *= -1.0
            for j in range(6):
                value = abs(rmat[j, j])
                if value < np.finfo(np.float64).tiny:
                    value = np.finfo(np.float64).tiny
                accum[j] += np.log(value)
            total += qr_interval * dt
            blocks[block_count, :] = accum / total
            block_count += 1
        else:
            qmat = qtrial
    return accum / total, blocks[:block_count, :], total, div_accum / ((n_steps - transient_steps) * dt)


def parameters_array(p) -> np.ndarray:
    """Serialize ModelParameters in the fixed order used by the JIT backend."""
    return np.array([
        p.kappa, p.gamma1, p.gamma2, p.omega1, p.omega2, p.g1, p.g2,
        p.theta, p.detuning, p.hopping, p.drive, p.drive_modulation,
        p.drive_frequency,
    ], dtype=np.float64)


def lyapunov_qr_numba(x0, p, *, n_steps: int, dt: float, transient_steps: int, qr_interval: int) -> dict:
    """Return the same fields as reconstruction_core.lyapunov_qr."""
    if n_steps <= transient_steps or qr_interval <= 0:
        raise ValueError("n_steps must exceed transient_steps and qr_interval > 0")
    spectrum, blocks, total, mean_divergence = _integrate(
        np.asarray(x0, dtype=np.float64), parameters_array(p), n_steps, dt, transient_steps, qr_interval
    )
    return {
        "spectrum": np.sort(spectrum)[::-1].tolist(),
        "block_history": blocks.tolist(),
        "total_time": float(total),
        "mean_divergence": float(mean_divergence),
        "qr_interval": qr_interval,
        "dt": dt,
        "transient_steps": transient_steps,
        "n_steps": n_steps,
        "backend": "numba-cpu",
    }
