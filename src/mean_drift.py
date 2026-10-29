"""Zero-speed, diagonal mean-drift loads from Capytaine Kochin functions.

All calculations here use the solver frame. Coefficients are real signed
loads per squared incident wave amplitude, not per squared wave height.
Only surge, sway and yaw are computed; cross-frequency QTFs are not exported.
"""

import numpy as np
import xarray as xr

DRIFT_DOFS = ("Surge", "Sway", "Yaw")


def kochin_angles(step_degrees=0.5):
    """Include the periodic seam and neighbors needed by the yaw derivative."""
    if isinstance(step_degrees, bool):
        raise ValueError("mean_drift_theta_step_deg must be a positive number")
    try:
        step = float(step_degrees)
    except (TypeError, ValueError) as error:
        raise ValueError("mean_drift_theta_step_deg must be a positive number") from error
    if not np.isfinite(step) or not 0 < step <= 10:
        raise ValueError("mean_drift_theta_step_deg must be in (0, 10]")
    intervals = round(360 / step)
    if not np.isclose(intervals * step, 360, rtol=0, atol=1e-10):
        raise ValueError("mean_drift_theta_step_deg must divide 360 degrees exactly")
    return np.arange(-1, intervals + 2) * (2*np.pi / intervals)


def angular_difference(dataset, motion, center_of_mass, fine):
    """Compare against half as many angular intervals, without another solve.

    Return a relative norm difference per physical component. This checks
    angular integration only; it does not check hull-panel convergence.
    """
    intervals = dataset.sizes['theta'] - 3
    if intervals % 2:
        return None
    core = dataset.isel(theta=slice(1, -1, 2))
    step = float(core.theta[1] - core.theta[0])
    coarse = xr.concat([
        core.isel(theta=[-2]).assign_coords(theta=[-step]), core,
        core.isel(theta=[1]).assign_coords(theta=[2*np.pi+step])], dim='theta',
        data_vars='minimal', coords='minimal', compat='equals')
    values = compute_mean_drift(coarse, motion, center_of_mass).values
    difference = np.linalg.norm(values-fine.values, axis=(0, 1))
    scale = np.linalg.norm(fine.values, axis=(0, 1))
    return difference / np.maximum(scale, np.finfo(float).eps)


def compute_mean_drift(dataset, motion, center_of_mass):
    """Return (omega, wave_direction, drift_dof) coefficients about CG.

    Capytaine's far-field moment is about the horizontal coordinate origin.
    Translate that moment to CG before applying any axis transformation.
    Evaluate one incident heading at a time: only the same-heading diagonal
    coefficients are needed, avoiding a full directional QTF allocation.
    The supplied motions must be expressed in the solver frame, with the
    same incident-wave phase origin as the Kochin functions.
    """
    from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force

    if float(dataset.forward_speed) != 0 or not np.isinf(float(dataset.water_depth)):
        raise ValueError("Mean-drift export requires zero speed and infinite water depth")
    cg = np.asarray(center_of_mass, dtype=float)
    if cg.shape != (3,) or not np.all(np.isfinite(cg)):
        raise ValueError("Mean-drift reference must contain three finite coordinates")
    motion = motion.transpose("omega", "wave_direction", "radiating_dof")
    if not np.array_equal(motion.omega, dataset.omega) or not np.array_equal(
            motion.wave_direction, dataset.wave_direction):
        raise ValueError("Mean-drift motion and Kochin frequency/heading grids must match")
    if not np.all(np.isfinite(motion)):
        raise ValueError("Mean-drift motion RAOs must be finite")
    values = np.empty((dataset.sizes['omega'], dataset.sizes['wave_direction'], 3))
    keys = ('drift_force_surge', 'drift_force_sway', 'drift_force_yaw')
    for heading in range(dataset.sizes['wave_direction']):
        selection = dict(wave_direction=slice(heading, heading + 1))
        result = far_field_mean_drift_force(motion.isel(**selection), dataset.isel(**selection))
        for component, key in enumerate(keys):
            raw = result[key].transpose('omega', 'wave_direction_k', 'wave_direction_l').values[:, 0, 0]
            if not np.all(np.isfinite(raw)) or np.any(
                    np.abs(raw.imag) > 1e-9*np.maximum(1, np.abs(raw.real))):
                raise ValueError("Diagonal mean-drift coefficients must be finite and real")
            values[:, heading, component] = raw.real
    values[:, :, 2] -= cg[0]*values[:, :, 1] - cg[1]*values[:, :, 0]
    return xr.DataArray(values, dims=('omega', 'wave_direction', 'drift_dof'),
                        coords=dict(omega=dataset.omega, wave_direction=dataset.wave_direction,
                                    drift_dof=list(DRIFT_DOFS)), name='mean_drift_force')
