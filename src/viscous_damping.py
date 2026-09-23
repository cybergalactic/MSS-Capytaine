"""Derive constant diagonal viscous damping from total-damping targets."""

from __future__ import annotations

import logging

import numpy as np
from scipy.optimize import brentq


LOG = logging.getLogger(__name__)
DOFS = ("surge", "sway", "heave", "roll", "pitch", "yaw")


def _parameter(settings: dict, name: str, *, positive: bool) -> float:
    try:
        number = float(settings[name])
    except KeyError as exc:
        raise ValueError(f"total_damping.{name} is required") from exc
    except (TypeError, ValueError) as exc:
        raise ValueError(f"total_damping.{name} must be a number") from exc
    if not np.isfinite(number) or (number <= 0 if positive else number < 0):
        condition = "positive" if positive else "nonnegative"
        raise ValueError(f"total_damping.{name} must be finite and {condition}")
    return number


def diagonal_viscous_damping(
    settings: dict | None,
    omega: np.ndarray,
    rigid_body_mass: np.ndarray,
    zero_added_mass: np.ndarray,
    zero_radiation_damping: np.ndarray,
    finite_added_mass: np.ndarray,
    finite_radiation_damping: np.ndarray,
    hydrostatic_stiffness: np.ndarray,
) -> np.ndarray:
    """Return a CG-referenced 6x6 Bv, constant over the export frequency grid.

    For surge, sway and yaw, use A(0) and assume B(0) = 0. For heave,
    roll and pitch, solve Cii = omega_n**2 * (MRBii + Aii(omega_n)) on the
    finite-frequency grid, then subtract Bii(omega_n) from the target Dii.
    A(omega_n) and B(omega_n) are linearly interpolated; the artificial
    10 rad/s infinite-frequency endpoint is never used for interpolation.
    """
    if settings is None:
        return np.zeros((6, 6), dtype=float)
    if not isinstance(settings, dict):
        raise ValueError("total_damping must be an object or null")

    expected_fields = {
        "surge_time_constant_s", "sway_time_constant_s", "heave_damping_ratio",
        "roll_damping_ratio", "pitch_damping_ratio", "yaw_time_constant_s",
    }
    unknown_fields = set(settings) - expected_fields
    if unknown_fields:
        raise ValueError(
            "Unknown total_damping entries: " + ", ".join(sorted(unknown_fields))
        )

    time_constants = {
        0: _parameter(settings, "surge_time_constant_s", positive=True),
        1: _parameter(settings, "sway_time_constant_s", positive=True),
        5: _parameter(settings, "yaw_time_constant_s", positive=True),
    }
    damping_ratios = {
        2: _parameter(settings, "heave_damping_ratio", positive=False),
        3: _parameter(settings, "roll_damping_ratio", positive=False),
        4: _parameter(settings, "pitch_damping_ratio", positive=False),
    }

    omega = np.asarray(omega, dtype=float)
    rigid_body_mass = np.asarray(rigid_body_mass, dtype=float)
    zero_added_mass = np.asarray(zero_added_mass, dtype=float)
    zero_radiation_damping = np.asarray(zero_radiation_damping, dtype=float)
    finite_added_mass = np.asarray(finite_added_mass, dtype=float)
    finite_radiation_damping = np.asarray(finite_radiation_damping, dtype=float)
    hydrostatic_stiffness = np.asarray(hydrostatic_stiffness, dtype=float)
    if omega.ndim != 1 or len(omega) == 0 or np.any(np.diff(omega) <= 0):
        raise ValueError("Finite frequencies must be a nonempty increasing vector")
    if not np.all(np.isfinite(omega)) or omega[0] <= 0:
        raise ValueError("Finite frequencies must be positive and finite")
    expected_shapes = (
        (rigid_body_mass, (6, 6), "MRB"),
        (zero_added_mass, (6, 6), "A(0)"),
        (zero_radiation_damping, (6, 6), "B(0)"),
        (finite_added_mass, (len(omega), 6, 6), "A(omega)"),
        (finite_radiation_damping, (len(omega), 6, 6), "B(omega)"),
        (hydrostatic_stiffness, (6, 6), "C"),
    )
    for matrix, shape, name in expected_shapes:
        if matrix.shape != shape or not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} must have shape {shape} and finite entries")

    grid = np.concatenate(([0.0], omega))
    viscous_diagonal = np.zeros(6, dtype=float)

    # Unrestored DOFs: their target Dii is Mii(0)/Ti, with Bii(0) assumed zero.
    for dof, time_constant in time_constants.items():
        effective_mass = rigid_body_mass[dof, dof] + zero_added_mass[dof, dof]
        if effective_mass <= 0:
            raise ValueError(f"MRB + A(0) must be positive in {DOFS[dof]}")
        viscous_diagonal[dof] = effective_mass / time_constant
        if abs(zero_radiation_damping[dof, dof]) > 1e-6 * viscous_diagonal[dof]:
            LOG.warning(
                "B(0) in %s is not negligible; its time-constant target "
                "assumes B(0) = 0", DOFS[dof],
            )

    # Restored DOFs: the natural frequency is estimated separately per DOF.
    for dof, damping_ratio in damping_ratios.items():
        stiffness = hydrostatic_stiffness[dof, dof]
        if stiffness <= 0:
            raise ValueError(f"Cii must be positive in {DOFS[dof]}")
        added_mass_curve = np.concatenate(
            ([zero_added_mass[dof, dof]], finite_added_mass[:, dof, dof])
        )
        radiation_curve = np.concatenate(
            ([zero_radiation_damping[dof, dof]], finite_radiation_damping[:, dof, dof])
        )
        mass_curve = rigid_body_mass[dof, dof] + added_mass_curve
        if np.any(mass_curve <= 0):
            raise ValueError(f"MRB + A(omega) must be positive in {DOFS[dof]}")

        def constraint(frequency: float) -> float:
            mass = rigid_body_mass[dof, dof] + np.interp(
                frequency, grid, added_mass_curve
            )
            return frequency**2 * mass - stiffness

        residuals = grid**2 * mass_curve - stiffness
        crossings = np.flatnonzero(residuals[1:] >= 0)
        if len(crossings) == 0:
            raise ValueError(
                f"Natural frequency in {DOFS[dof]} exceeds the finite "
                f"frequency grid ({omega[-1]:g} rad/s); extend omega_rad_s"
            )
        upper = int(crossings[0]) + 1
        natural_frequency = brentq(constraint, grid[upper - 1], grid[upper])
        effective_mass = rigid_body_mass[dof, dof] + np.interp(
            natural_frequency, grid, added_mass_curve
        )
        potential_damping = float(np.interp(
            natural_frequency, grid, radiation_curve
        ))
        target_damping = 2.0 * damping_ratio * natural_frequency * effective_mass
        remainder = target_damping - potential_damping
        if remainder < 0:
            LOG.warning(
                "%s: B(omega_n)=%.6g exceeds target Dii=%.6g at "
                "omega_n=%.6g rad/s; setting Bvii=0",
                DOFS[dof], potential_damping, target_damping, natural_frequency,
            )
        viscous_diagonal[dof] = max(0.0, remainder)

    return np.diag(viscous_diagonal)
