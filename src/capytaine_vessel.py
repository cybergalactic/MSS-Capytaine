"""Export a CG-referenced Capytaine solve as an MSS-style vessel.mat."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import savemat
import xarray as xr


DOFS = ("Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw")

# Reflection in the x-z plane:
#
#   surge  +     roll  -
#   sway   -     pitch +
#   heave  +     yaw   -
#
# The same parity applies to forces/moments and rigid-body motions.
SYMMETRY_SIGNS = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0])


def _mirror_headings(
    values: np.ndarray,
    headings: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Mirror 0...180 deg data to a full 0...360 deg directional set.

    Parameters
    ----------
    values :
        Complex array with shape

            frequency x heading x DOF

        for headings from 0 through pi.

    headings :
        Heading angles in radians, including both 0 and pi.

    Returns
    -------
    values_full :
        Data for the full directional set.

    headings_full :
        Corresponding headings in radians.

    Notes
    -----
    For a port-starboard symmetric hull,

        beta -> 2*pi - beta

    and the sway, roll and yaw components change sign.
    """

    values = np.asarray(values)
    headings = np.asarray(headings, dtype=float)

    if values.ndim != 3 or values.shape[2] != 6:
        raise ValueError(
            "Directional data must have shape frequency x heading x 6"
        )

    if values.shape[1] != len(headings):
        raise ValueError(
            "Heading dimension of directional data does not match headings"
        )

    if len(headings) < 2:
        raise ValueError(
            "At least two headings are required for symmetry mirroring"
        )

    if not np.isclose(headings[0], 0.0):
        raise ValueError(
            "First solved heading must be 0 degrees"
        )

    if not np.isclose(headings[-1], np.pi):
        raise ValueError(
            "Last solved heading must be 180 degrees"
        )

    # Sway, roll, and yaw are odd under port-starboard reflection. Their
    # force and motion RAOs must vanish for head and following seas; small
    # BEM residuals can otherwise be strongly amplified at low frequency.
    values = values.copy()
    odd_dofs = SYMMETRY_SIGNS < 0
    values[:, 0, odd_dofs] = 0.0
    values[:, -1, odd_dofs] = 0.0

    # Example:
    #
    # solved:
    #   0, 10, ..., 170, 180
    #
    # mirrored:
    #   190, 200, ..., 340, 350
    #
    # Hence reverse the interior headings 170,...,10.

    mirrored_headings = (
        2.0 * np.pi - headings[-2:0:-1]
    )

    mirrored_values = (
        values[:, -2:0:-1, :]
        * SYMMETRY_SIGNS[None, None, :]
    )

    headings_full = np.concatenate(
        (headings, mirrored_headings)
    )

    values_full = np.concatenate(
        (values, mirrored_values),
        axis=1,
    )

    return values_full, headings_full


def _rao_cells(
    values: np.ndarray,
) -> dict[str, np.ndarray]:
    """Store six frequency x heading x speed arrays using MSS cell layout."""

    # Capytaine uses exp(-i*w*t), while MSS uses exp(+i*w*t).
    values = np.conj(values)

    amp = np.empty((1, 6), dtype=object)
    phase = np.empty((1, 6), dtype=object)

    for dof in range(6):

        amp[0, dof] = (
            np.abs(values[:, :, dof])[:, :, None]
        )

        phase[0, dof] = (
            np.mod(
                np.angle(values[:, :, dof]),
                2.0 * np.pi,
            )[:, :, None]
        )

    return {
        "amp": amp,
        "phase": phase,
    }


def write_vessel(
    dataset: xr.Dataset,
    hydrostatics: dict,
    vertices: np.ndarray,
    faces: np.ndarray,
    name: str,
    rho: float,
    gravity: float,
    path: Path,
    zero_added_mass: np.ndarray,
    zero_radiation_damping: np.ndarray,
    infinite_added_mass: np.ndarray,
    infinite_radiation_damping: np.ndarray,
    viscous_damping_matrix: np.ndarray,
) -> None:
    """Write zero-speed Capytaine results as an MSS vessel structure.

    Hydrodynamic quantities use MSS FSD axes and are referenced to CG.

    The hydrodynamic frequency vector is

        vessel.freqs = [0, finite Capytaine frequencies, 10]

    The data stored at 10 rad/s are Capytaine's omega=infinity radiation
    solution. This endpoint is included in A, B, Bv, and C, but not in RAOs.

    RAOs contain only the positive finite frequencies.

    The Capytaine solution is assumed to contain headings from 0 to 180 deg
    for a port-starboard symmetric monohull. These are mirrored to the full
    0 to 350 deg directional set in the MATLAB vessel structure.

    Bv is a configured constant diagonal matrix, repeated at every
    coefficient frequency. Second-order drift forces are not included.
    """

    # ------------------------------------------------------------------
    # Coordinate-system checks
    # ------------------------------------------------------------------

    if (
        dataset.attrs.get("matrix_reference") != "CG"
        or "MSS FSD"
        not in dataset.attrs.get("coordinate_system", "")
    ):
        raise ValueError(
            "MSS vessel export requires FSD data referenced to CG"
        )

    # ------------------------------------------------------------------
    # Finite frequencies and solved headings
    # ------------------------------------------------------------------

    omega = np.asarray(
        dataset.omega.values,
        dtype=float,
    )

    headings = np.asarray(
        dataset.wave_direction.values,
        dtype=float,
    )

    if (
        omega.ndim != 1
        or len(omega) == 0
        or np.any(omega <= 0)
        or np.any(omega >= 10.0)
        or not np.all(np.isfinite(omega))
    ):
        raise ValueError(
            "RAO frequencies must be positive, finite, and below 10 rad/s"
        )

    # Hydrodynamic matrices additionally contain zero and infinity. Infinity
    # is represented by the conventional practical endpoint at 10 rad/s.
    frequencies = np.concatenate(
        ([0.0], omega, [10.0])
    )

    # ------------------------------------------------------------------
    # Rigid-body and hydrostatic quantities
    # ------------------------------------------------------------------

    mass_matrix = np.asarray(
        dataset.inertia_matrix.values,
        dtype=float,
    )

    mass = float(
        hydrostatics["body_mass_kg"]
    )

    if not np.allclose(
        np.diag(mass_matrix)[:3],
        mass,
        rtol=1e-5,
    ):
        raise ValueError(
            "The first three mass-matrix diagonal entries "
            "must equal body mass"
        )

    cg = np.asarray(
        hydrostatics["center_of_mass_m"],
        dtype=float,
    )

    cb = np.asarray(
        hydrostatics["center_of_buoyancy_m"],
        dtype=float,
    )

    # vertices are still expressed in the original Capytaine mesh axes.
    length = float(
        np.ptp(vertices[:, 0])
    )

    breadth = float(
        2.0 * np.max(np.abs(vertices[:, 1]))
    )

    draft = float(
        -np.min(vertices[:, 2])
    )

    if min(length, breadth, draft) <= 0:
        raise ValueError(
            "Offset hull must have positive length, breadth, and draft"
        )

    volume = float(
        hydrostatics["volume_m3"]
    )

    stiffness = np.asarray(
        dataset.hydrostatic_stiffness.values,
        dtype=float,
    )

    main = {
        "name": name,
        "g": gravity,
        "rho": rho,

        "Lpp": length,
        "Lwl": length,
        "B": breadth,
        "T": draft,

        "m": mass,
        "nabla": volume,

        "C_B": (
            volume
            / (length * breadth * draft)
        ),

        # MSS FSD coordinates:
        # x forward, y starboard, z down,
        # waterline origin.
        "CG": cg,
        "CB": cb,

        "k44": np.sqrt(
            mass_matrix[3, 3] / mass
        ),

        "k55": np.sqrt(
            mass_matrix[4, 4] / mass
        ),

        "k66": np.sqrt(
            mass_matrix[5, 5] / mass
        ),

        "GM_T": (
            stiffness[3, 3]
            / (mass * gravity)
        ),

        "GM_L": (
            stiffness[4, 4]
            / (mass * gravity)
        ),
    }

    # ------------------------------------------------------------------
    # Added mass
    #
    # MATLAB:
    #
    #   A(:,:,1,1)     = A(0)
    #   A(:,:,2:end-1) = finite-frequency A(w)
    #   A(:,:,end,1)   = A(infinity), labeled as 10 rad/s
    # ------------------------------------------------------------------

    zero_added_mass = np.asarray(
        zero_added_mass,
        dtype=float,
    )

    if zero_added_mass.shape != (6, 6):
        raise ValueError(
            "zero_added_mass must be a 6 x 6 matrix"
        )

    infinite_added_mass = np.asarray(
        infinite_added_mass,
        dtype=float,
    )

    if infinite_added_mass.shape != (6, 6):
        raise ValueError(
            "infinite_added_mass must be a 6 x 6 matrix"
        )

    added_mass_finite = np.asarray(
        dataset.added_mass.values,
        dtype=float,
    ).transpose(1, 2, 0)

    added_mass = np.concatenate(
        (
            zero_added_mass[:, :, None],
            added_mass_finite,
            infinite_added_mass[:, :, None],
        ),
        axis=2,
    )[:, :, :, None]

    # ------------------------------------------------------------------
    # Potential damping
    #
    # MATLAB:
    #
    #   B(:,:,1,1)     = B(0)
    #   B(:,:,2:end-1) = finite-frequency B(w)
    #   B(:,:,end,1)   = B(infinity), labeled as 10 rad/s
    # ------------------------------------------------------------------

    zero_radiation_damping = np.asarray(
        zero_radiation_damping,
        dtype=float,
    )

    if zero_radiation_damping.shape != (6, 6):
        raise ValueError(
            "zero_radiation_damping must be a 6 x 6 matrix"
        )

    infinite_radiation_damping = np.asarray(
        infinite_radiation_damping,
        dtype=float,
    )

    if infinite_radiation_damping.shape != (6, 6):
        raise ValueError(
            "infinite_radiation_damping must be a 6 x 6 matrix"
        )

    damping_finite = np.asarray(
        dataset.radiation_damping.values,
        dtype=float,
    ).transpose(1, 2, 0)

    damping = np.concatenate(
        (
            zero_radiation_damping[:, :, None],
            damping_finite,
            infinite_radiation_damping[:, :, None],
        ),
        axis=2,
    )[:, :, :, None]

    # ------------------------------------------------------------------
    # Hydrostatic restoring
    #
    # C is frequency independent, but MSS stores it over the same
    # frequency dimension as A and B.
    # ------------------------------------------------------------------

    restoring = np.repeat(
        stiffness[:, :, None, None],
        len(frequencies),
        axis=2,
    )

    # ------------------------------------------------------------------
    # Viscous damping
    #
    # The same CG-referenced diagonal matrix applies at all frequencies.
    # ------------------------------------------------------------------

    viscous_damping_matrix = np.asarray(viscous_damping_matrix, dtype=float)
    if (
        viscous_damping_matrix.shape != (6, 6)
        or not np.all(np.isfinite(viscous_damping_matrix))
        or not np.allclose(
            viscous_damping_matrix,
            np.diag(np.diag(viscous_damping_matrix)),
        )
        or np.any(np.diag(viscous_damping_matrix) < 0)
    ):
        raise ValueError("Bv must be a finite, nonnegative diagonal 6x6 matrix")
    viscous_damping = np.repeat(
        viscous_damping_matrix[:, :, None], len(frequencies), axis=2
    )

    # ------------------------------------------------------------------
    # Force RAOs
    #
    # Only finite positive frequencies are included.
    # ------------------------------------------------------------------

    force = np.asarray(
        dataset.excitation_force.transpose(
            "omega",
            "wave_direction",
            "influenced_dof",
        ).values
    )

    force, headings_full = _mirror_headings(
        force,
        headings,
    )

    force_rao = _rao_cells(
        force
    )

    force_rao["w"] = omega

    # ------------------------------------------------------------------
    # Motion RAOs
    #
    # Only finite positive frequencies are included.
    # ------------------------------------------------------------------

    motion = np.asarray(
        dataset.motion_rao.transpose(
            "omega",
            "wave_direction",
            "radiating_dof",
        ).values
    )

    motion, headings_motion = _mirror_headings(
        motion,
        headings,
    )

    if not np.allclose(
        headings_motion,
        headings_full,
    ):
        raise RuntimeError(
            "Force and motion RAO headings are inconsistent"
        )

    motion_rao = _rao_cells(
        motion
    )

    motion_rao["w"] = omega

    # ------------------------------------------------------------------
    # MSS vessel structure
    # ------------------------------------------------------------------

    vessel = {
        "main": main,

        "MRB": mass_matrix,

        "A": added_mass,
        "B": damping,
        "Bv": viscous_damping,
        "C": restoring,

        # A, B, Bv and C use this frequency vector.
        "freqs": frequencies,

        # Full 360-degree directional set.
        "headings": headings_full,

        "velocities": np.array([0.0]),

        # RAO.w contains finite frequencies only.
        "forceRAO": force_rao,
        "motionRAO": motion_rao,

        "hydrodynamic_reference": "CG",
        "hydrodynamic_axes": "MSS FSD",
        "hydrodynamic_source": (
            "Capytaine monohull offset points"
        ),

        "viscous_damping_model": (
            "constant diagonal Bv calibrated at zero/natural frequencies"
            if np.any(np.diag(viscous_damping_matrix) > 0)
            else "none; Bv is zero"
        ),

        "second_order_drift": (
            "not computed"
        ),
    }

    # ------------------------------------------------------------------
    # MATLAB output
    # ------------------------------------------------------------------

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    savemat(
        path,
        {"vessel": vessel},
        do_compression=True,
    )

    # ------------------------------------------------------------------
    # Panel geometry for MATLAB/Python inspection
    # ------------------------------------------------------------------

    np.savez(
        path.parent
        / "generated_hull_panels.npz",
        vertices=vertices,
        faces=faces,
    )

    savemat(
        path.parent
        / "generated_hull_panels.mat",
        {
            "vertices": vertices,
            "faces": faces + 1,
        },
        do_compression=True,
    )
