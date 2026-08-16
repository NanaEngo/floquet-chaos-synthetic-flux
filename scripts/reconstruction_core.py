#!/usr/bin/env python3
"""Core equations and diagnostics for the Chaos reconstruction.

All quantities in this module are dimensionless and use the declared reference
frequency.  The state is x=(Re(alpha), Im(alpha), Re(beta_1), Im(beta_1),
Re(beta_2), Im(beta_2)).  No result is considered publication evidence until a
caller writes a manifest and passes the corresponding gate.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class ModelParameters:
    kappa: float = 1.0
    gamma1: float = 0.02
    gamma2: float = 0.02
    omega1: float = 1.0
    omega2: float = 1.03
    g1: float = 0.02
    g2: float = 0.018
    hopping: float = 0.08
    detuning: float = -1.0
    drive: float = 0.20
    drive_modulation: float = 0.10
    drive_frequency: float = 1.0
    theta: float = 0.0
    force: float = 0.0

    def to_dict(self) -> dict[str, float]:
        return asdict(self)

    @property
    def period(self) -> float:
        if self.drive_frequency <= 0:
            raise ValueError("drive_frequency must be positive")
        return 2.0 * np.pi / self.drive_frequency

    @property
    def divergence(self) -> float:
        return -self.kappa - self.gamma1 - self.gamma2


def synthetic_flux(g1: complex, hopping: complex, g2: complex) -> float:
    """Return the gauge-invariant loop phase arg(g1*J*conj(g2))."""
    product = g1 * hopping * np.conjugate(g2)
    if abs(product) == 0:
        raise ValueError("loop couplings must be non-zero")
    return float(np.angle(product))


def drive_amplitude(phi: float, p: ModelParameters) -> float:
    return p.drive * (1.0 + p.drive_modulation * np.cos(phi))


def rhs_phase(x: np.ndarray, phi: float, p: ModelParameters) -> np.ndarray:
    """Right-hand side derived from one fixed Hamiltonian/QLE convention."""
    if np.asarray(x).shape != (6,):
        raise ValueError(f"expected state shape (6,), got {np.asarray(x).shape}")
    ar, ai, b1r, b1i, b2r, b2i = np.asarray(x, dtype=float)
    c, s = np.cos(p.theta), np.sin(p.theta)
    omega_eff = p.detuning + 2.0 * p.g1 * b1r + 2.0 * p.g2 * b2r
    e = drive_amplitude(phi, p)

    return np.array([
        -0.5 * p.kappa * ar - omega_eff * ai + e,
        omega_eff * ar - 0.5 * p.kappa * ai,
        -0.5 * p.gamma1 * b1r + p.omega1 * b1i
            + p.hopping * (s * b2r + c * b2i),
        p.g1 * (ar * ar + ai * ai) - p.omega1 * b1r
            - 0.5 * p.gamma1 * b1i
            - p.hopping * c * b2r + p.hopping * s * b2i
            + p.force,
        -0.5 * p.gamma2 * b2r + p.omega2 * b2i
            - p.hopping * s * b1r + p.hopping * c * b1i,
        p.g2 * (ar * ar + ai * ai) - p.omega2 * b2r
            - 0.5 * p.gamma2 * b2i
            - p.hopping * c * b1r - p.hopping * s * b1i,
    ], dtype=float)


def jacobian_phase(x: np.ndarray, phi: float, p: ModelParameters) -> np.ndarray:
    """Analytic Jacobian of :func:`rhs_phase` with respect to the six states."""
    ar, ai, b1r, _b1i, b2r, _b2i = np.asarray(x, dtype=float)
    c, s = np.cos(p.theta), np.sin(p.theta)
    omega_eff = p.detuning + 2.0 * p.g1 * b1r + 2.0 * p.g2 * b2r
    return np.array([
        [-0.5*p.kappa, -omega_eff, -2*p.g1*ai, 0, -2*p.g2*ai, 0],
        [omega_eff, -0.5*p.kappa, 2*p.g1*ar, 0, 2*p.g2*ar, 0],
        [0, 0, -0.5*p.gamma1, p.omega1, p.hopping*s, p.hopping*c],
        [2*p.g1*ar, 2*p.g1*ai, -p.omega1, -0.5*p.gamma1, -p.hopping*c, p.hopping*s],
        [0, 0, -p.hopping*s, p.hopping*c, -0.5*p.gamma2, p.omega2],
        [2*p.g2*ar, 2*p.g2*ai, -p.hopping*c, -p.hopping*s, -p.omega2, -0.5*p.gamma2],
    ], dtype=float)


def phase_rhs_augmented(z: np.ndarray, p: ModelParameters) -> np.ndarray:
    x = np.asarray(z[:6], dtype=float)
    phi = float(z[6])
    return np.r_[rhs_phase(x, phi, p), p.drive_frequency]


def integrate_one_period(x0: np.ndarray, p: ModelParameters, *, rtol: float = 1e-9,
                         atol: float = 1e-11, dense: bool = False):
    """Integrate the forced system from phase zero through one drive period."""
    def fun(t, x):
        return rhs_phase(x, p.drive_frequency * t, p)
    return solve_ivp(fun, (0.0, p.period), np.asarray(x0, dtype=float),
                     method="DOP853", rtol=rtol, atol=atol,
                     dense_output=dense, max_step=p.period / 200.0)


def find_periodic_orbit(x0: np.ndarray, p: ModelParameters, *, max_iter: int = 100,
                        residual_tol: float = 1e-8, rtol: float = 1e-9,
                        atol: float = 1e-11) -> dict:
    """Find a drive-locked orbit by fixed-point iteration of the Poincare map."""
    x = np.asarray(x0, dtype=float).copy()
    history = []
    for iteration in range(1, max_iter + 1):
        sol = integrate_one_period(x, p, rtol=rtol, atol=atol)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])):
            return {"status": "FAIL", "reason": "period integration failed", "iterations": iteration}
        x_next = sol.y[:, -1]
        residual = float(np.linalg.norm(x_next - x, ord=np.inf))
        history.append(residual)
        x = x_next
        if residual < residual_tol:
            return {"status": "PASS", "x0": x.tolist(), "residual": residual,
                    "iterations": iteration, "residual_history": history}
    return {"status": "FAIL", "reason": "periodic-orbit residual did not converge",
            "x0": x.tolist(), "residual": history[-1], "iterations": max_iter,
            "residual_history": history}


def monodromy(x0: np.ndarray, p: ModelParameters, *, rtol: float = 1e-9,
              atol: float = 1e-11) -> dict:
    """Integrate the state and variational equation over one drive period."""
    n = 6
    y0 = np.r_[np.asarray(x0, dtype=float), np.eye(n).ravel()]
    def fun(t, y):
        x = y[:n]
        phi = p.drive_frequency * t
        A = jacobian_phase(x, phi, p)
        return np.r_[rhs_phase(x, phi, p), (A @ y[n:].reshape(n, n)).ravel()]
    sol = solve_ivp(fun, (0.0, p.period), y0, method="DOP853", rtol=rtol,
                    atol=atol, max_step=p.period / 200.0)
    if not sol.success:
        raise RuntimeError(sol.message)
    M = sol.y[n:, -1].reshape(n, n)
    multipliers = np.linalg.eigvals(M)
    rates = np.log(np.maximum(np.abs(multipliers), np.finfo(float).tiny)) / p.period
    multiplier_records = [{"real": float(z.real), "imag": float(z.imag)} for z in multipliers]
    return {"monodromy": M.tolist(), "multipliers": multiplier_records,
            "floquet_rates": rates.tolist(), "solver_nfev": int(sol.nfev)}


def lyapunov_qr(x0: np.ndarray, p: ModelParameters, *, n_steps: int = 10000,
                dt: float = 0.01, transient_steps: int = 1000,
                qr_interval: int = 10) -> dict:
    """Compute a full six-dimensional QR/Benettin spectrum with block history."""
    if n_steps <= transient_steps or qr_interval <= 0:
        raise ValueError("n_steps must exceed transient_steps and qr_interval > 0")
    x = np.asarray(x0, dtype=float).copy()
    phi = 0.0
    Q = np.eye(6)
    accum = np.zeros(6)
    blocks = []
    div_accum = 0.0
    total = 0

    def rk4_state(x_, phi_):
        k1 = rhs_phase(x_, phi_, p)
        k2 = rhs_phase(x_ + 0.5*dt*k1, phi_ + 0.5*p.drive_frequency*dt, p)
        k3 = rhs_phase(x_ + 0.5*dt*k2, phi_ + 0.5*p.drive_frequency*dt, p)
        k4 = rhs_phase(x_ + dt*k3, phi_ + p.drive_frequency*dt, p)
        return x_ + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0

    for step in range(n_steps):
        x_new = rk4_state(x, phi)
        A = jacobian_phase(x, phi, p)
        A_new = jacobian_phase(x_new, phi + p.drive_frequency*dt, p)
        Q_trial = Q + 0.5 * dt * (A @ Q + A_new @ (Q + dt*A @ Q))
        x, phi = x_new, (phi + p.drive_frequency*dt) % (2.0*np.pi)
        if step < transient_steps:
            Q = np.eye(6)
            continue
        # Accumulate divergence over the same post-transient interval used
        # for the QR spectrum, not over the discarded transient.
        div_accum += float(np.trace(A)) * dt
        if (step - transient_steps + 1) % qr_interval == 0:
            Q, R = np.linalg.qr(Q_trial)
            signs = np.sign(np.diag(R)); signs[signs == 0] = 1.0
            Q *= signs[np.newaxis, :]
            diag = np.maximum(np.abs(np.diag(R)), np.finfo(float).tiny)
            accum += np.log(diag)
            total += qr_interval * dt
            blocks.append((accum / total).tolist())
        else:
            Q = Q_trial
    if not blocks:
        raise RuntimeError("no QR blocks were accumulated")
    spectrum = np.sort(accum / total)[::-1]
    return {            "spectrum": spectrum.tolist(), "block_history": blocks,
            "total_time": total, "mean_divergence": div_accum / ((n_steps - transient_steps)*dt),

            "qr_interval": qr_interval, "dt": dt,
            "transient_steps": transient_steps, "n_steps": n_steps}


def finite_difference_jacobian(x: np.ndarray, phi: float, p: ModelParameters,
                               eps: float = 1e-7) -> np.ndarray:
    J = np.zeros((6, 6), dtype=float)
    for col in range(6):
        xp = np.asarray(x, dtype=float).copy(); xm = xp.copy()
        xp[col] += eps; xm[col] -= eps
        J[:, col] = (rhs_phase(xp, phi, p) - rhs_phase(xm, phi, p)) / (2*eps)
    return J
