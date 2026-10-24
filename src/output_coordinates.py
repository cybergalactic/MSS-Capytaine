"""Internal transforms between solver axes and the public MSS FSD convention."""

import numpy as np
import xarray as xr


DOF_SIGNS = np.array([-1.0, 1.0, -1.0, -1.0, 1.0, -1.0])
POSITION_SIGNS = np.array([-1.0, 1.0, -1.0])


def transform_positions(values: np.ndarray) -> np.ndarray:
    """Rotate positions in either direction (the transform is its own inverse)."""
    return np.asarray(values, dtype=float) * POSITION_SIGNS


def transform_matrix(values: np.ndarray) -> np.ndarray:
    """Rotate a generalized 6x6 matrix in either direction at the same origin."""
    return DOF_SIGNS[:, None] * np.asarray(values, dtype=float) * DOF_SIGNS[None, :]


def to_mss_fsd(dataset: xr.Dataset) -> xr.Dataset:
    """Rotate generalized results by 180 degrees about y, preserving CG.

    The offset mesh frame has x aft, y starboard, z up. MSS FSD uses x forward,
    y starboard, z down. Wave headings are propagation directions in each frame.
    """
    converted = dataset.copy(deep=True)
    sign_force = xr.DataArray(DOF_SIGNS, dims=["influenced_dof"],
                              coords={"influenced_dof": dataset.influenced_dof})
    sign_motion = xr.DataArray(DOF_SIGNS, dims=["radiating_dof"],
                               coords={"radiating_dof": dataset.radiating_dof})
    for name, value in dataset.data_vars.items():
        signs = 1
        if "influenced_dof" in value.dims:
            signs = signs * sign_force
        if "radiating_dof" in value.dims:
            signs = signs * sign_motion
        if name == "center_of_buoyancy":
            signs = xr.DataArray(POSITION_SIGNS, dims=["space_coordinate"],
                                 coords={"space_coordinate": dataset.space_coordinate})
        converted[name] = value * signs

    headings = np.mod(np.pi - np.asarray(dataset.wave_direction), 2 * np.pi)
    headings[np.isclose(headings, 2 * np.pi, atol=1e-12)] = 0.0
    converted = converted.assign_coords(wave_direction=headings).sortby("wave_direction")
    converted.attrs.update({
        "coordinate_system": "MSS FSD: x forward, y starboard, z down",
        "coordinate_origin": dataset.attrs.get(
            "coordinate_origin",
            "midships, centerline, design waterline",
        ),
        "matrix_reference": "CG", "force_reference": "CG",
        "motion_reference": "CG",
        "wave_direction_definition": "propagation angle from FSD +x toward +y",
        "source_frame": "internal solver frame",
    })
    return converted
