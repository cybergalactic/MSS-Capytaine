"""Derive constant diagonal damping added to potential-flow radiation damping."""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq


DOFS = ("surge", "sway", "heave", "roll", "pitch", "yaw")


def _parameter(settings: dict, name: str) -> float:
    try:
        number = float(settings[name])
    except KeyError as exc:
        raise ValueError(f"viscous_damping.{name} is required") from exc
    except (TypeError, ValueError) as exc:
        raise ValueError(f"viscous_damping.{name} must be a number") from exc
    if not np.isfinite(number) or number < 0:
        raise ValueError(f"viscous_damping.{name} must be finite and nonnegative")
    return number


def diagonal_viscous_damping(
    settings: dict | None,
    omega: np.ndarray,
    rigid_body_mass: np.ndarray,
    zero_added_mass: np.ndarray,
    finite_added_mass: np.ndarray,
    hydrostatic_stiffness: np.ndarray,
) -> np.ndarray:
    """Return a CG-referenced 6x6 Bv, constant over the export frequency grid.

    Surge, sway and yaw use Bvii = (MRBii + Aii(0)) / Ti. A zero time
    constant disables that entry. Heave, roll and pitch use
    Bvii = 2*zeta_v*omega_n*(MRBii + Aii(omega_n)), where omega_n solves
    Cii = omega_n**2 * (MRBii + Aii(omega_n)). Zero additional ratio
    disables that entry. Only finite-frequency added mass is interpolated.
    """
    if settings is None:
        return np.zeros((6, 6), dtype=float)
    if not isinstance(settings, dict):
        raise ValueError("viscous_damping must be an object or null")

    expected_fields = {
        "surge_viscous_time_constant_s", "sway_viscous_time_constant_s",
        "heave_additional_damping_ratio", "roll_additional_damping_ratio",
        "pitch_additional_damping_ratio", "yaw_viscous_time_constant_s",
    }
    unknown_fields = set(settings) - expected_fields
    if unknown_fields:
        raise ValueError(
            "Unknown viscous_damping entries: " + ", ".join(sorted(unknown_fields))
        )

    time_constants = {
        0: _parameter(settings, "surge_viscous_time_constant_s"),
        1: _parameter(settings, "sway_viscous_time_constant_s"),
        5: _parameter(settings, "yaw_viscous_time_constant_s"),
    }
    damping_ratios = {
        2: _parameter(settings, "heave_additional_damping_ratio"),
        3: _parameter(settings, "roll_additional_damping_ratio"),
        4: _parameter(settings, "pitch_additional_damping_ratio"),
    }

    omega = np.asarray(omega, dtype=float)
    rigid_body_mass = np.asarray(rigid_body_mass, dtype=float)
    zero_added_mass = np.asarray(zero_added_mass, dtype=float)
    finite_added_mass = np.asarray(finite_added_mass, dtype=float)
    hydrostatic_stiffness = np.asarray(hydrostatic_stiffness, dtype=float)
    if omega.ndim != 1 or len(omega) == 0 or np.any(np.diff(omega) <= 0):
        raise ValueError("Finite frequencies must be a nonempty increasing vector")
    if not np.all(np.isfinite(omega)) or omega[0] <= 0:
        raise ValueError("Finite frequencies must be positive and finite")
    expected_shapes = (
        (rigid_body_mass, (6, 6), "MRB"),
        (zero_added_mass, (6, 6), "A(0)"),
        (finite_added_mass, (len(omega), 6, 6), "A(omega)"),
        (hydrostatic_stiffness, (6, 6), "C"),
    )
    for matrix, shape, name in expected_shapes:
        if matrix.shape != shape or not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} must have shape {shape} and finite entries")

    grid = np.concatenate(([0.0], omega))
    viscous_diagonal = np.zeros(6, dtype=float)

    # Unrestored DOFs: time constants describe added damping only.
    for dof, time_constant in time_constants.items():
        if time_constant == 0:
            continue
        effective_mass = rigid_body_mass[dof, dof] + zero_added_mass[dof, dof]
        if effective_mass <= 0:
            raise ValueError(f"MRB + A(0) must be positive in {DOFS[dof]}")
        viscous_diagonal[dof] = effective_mass / time_constant

    # Restored DOFs: the natural frequency is estimated separately per DOF.
    for dof, damping_ratio in damping_ratios.items():
        if damping_ratio == 0:
            continue
        stiffness = hydrostatic_stiffness[dof, dof]
        if stiffness <= 0:
            raise ValueError(f"Cii must be positive in {DOFS[dof]}")
        added_mass_curve = np.concatenate(
            ([zero_added_mass[dof, dof]], finite_added_mass[:, dof, dof])
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
        viscous_diagonal[dof] = (
            2.0 * damping_ratio * natural_frequency * effective_mass
        )

    return np.diag(viscous_diagonal)
