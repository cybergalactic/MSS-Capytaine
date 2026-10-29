"""Zero-speed, six-DOF Capytaine solver with MATLAB vessel export."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np

from .output_coordinates import transform_matrix, transform_positions
from .damping import validate_damping as _viscous_damping_parameters


DOFS = ("Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw")
MSS_WAVE_DIRECTIONS_DEG = np.arange(0.0, 181.0, 10.0)
LOG = logging.getLogger(__name__)
RECIPROCITY_WARNING_THRESHOLD = 0.01
RECIPROCITY_GLOBAL_WARNING_THRESHOLD = 1e-3


def _vector3(value: object, name: str) -> tuple[float, float, float]:
    array = np.asarray(value, dtype=float)
    if array.shape != (3,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain three finite coordinates")
    return tuple(float(x) for x in array)


def _positive(value: object, name: str) -> float:
    number = float(value)
    if not np.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _check_transverse_stability(config: dict, gm_t: float) -> None:
    """Reject unstable roll or disagreement with an independent GM reference."""
    z_g = float(config["center_of_mass_m"][2])
    axes = (f"Configured zG={z_g:g} m. JSON uses MSS FSD (z positive down); "
            "surface-vessel centers are relative to the waterline, not keel heights. ")
    if not np.isfinite(gm_t) or gm_t <= 0:
        raise ValueError(
            f"Calculated GM_T={gm_t:g} m is not positive. " + axes +
            "Check the CG sign, origin, hull geometry and loading condition before solving."
        )
    if "GM_T_m" not in config:
        if "GM_T_tolerance_m" in config:
            raise ValueError("GM_T_tolerance_m requires GM_T_m")
        return
    reference = _positive(config["GM_T_m"], "GM_T_m")
    tolerance = _positive(config.get("GM_T_tolerance_m", max(0.05, 0.02*reference)),
                          "GM_T_tolerance_m")
    if abs(gm_t-reference) > tolerance:
        implied_z = z_g + reference-gm_t
        raise ValueError(
            f"GM_T consistency check failed: calculated {gm_t:.4g} m, "
            f"reference GM_T_m={reference:.4g} m, tolerance {tolerance:.4g} m. "
            + axes + f"For this hull and reference GM, zG should be approximately "
            f"{implied_z:.4g} m in FSD. Check the sign and origin; do not alter GM "
            "merely to make the check pass."
        )


def _mass_matrix(
    config: dict,
    base: Path,
    mass: float,
) -> np.ndarray:
    """Read a CSV mass matrix at CG or build it from gyration radii.

    Optional R46_m is a signed cross-radius: I46=m*R46_m*abs(R46_m).
    Omitting it retains the original diagonal inertia.
    """
    matrix_file = config.get("inertia_matrix_csv")
    if matrix_file:
        if "R46_m" in config:
            raise ValueError("R46_m is used with gyration radii, not inertia_matrix_csv")
        matrix = np.loadtxt((base / matrix_file).resolve(), delimiter=",")
    else:
        required = ("R44_m", "R55_m", "R66_m")
        if not all(name in config for name in required):
            raise ValueError(
                "Provide inertia_matrix_csv or R44_m, R55_m, and R66_m"
            )
        radii = np.array(
            [_positive(config[name], name) for name in required]
        )
        matrix = np.diag(
            np.concatenate((np.full(3, mass), mass * radii**2))
        )
        cross_radius = config.get("R46_m", 0.0)
        if (isinstance(cross_radius, bool) or not isinstance(cross_radius, (int, float))
                or not np.isfinite(cross_radius)):
            raise ValueError("R46_m must be a finite signed radius in metres")
        matrix[3, 5] = matrix[5, 3] = mass * cross_radius * abs(cross_radius)

    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (6, 6) or not np.all(np.isfinite(matrix)):
        raise ValueError(
            "The mass matrix must contain 6 rows of 6 finite numbers"
        )
    if not np.allclose(matrix, matrix.T, rtol=1e-8, atol=1e-8):
        raise ValueError("The mass matrix must be symmetric")
    if np.min(np.linalg.eigvalsh(matrix)) <= 0.0:
        raise ValueError("The mass matrix must be positive definite")
    if not np.allclose(np.diag(matrix)[:3], mass, rtol=1e-5):
        raise ValueError(
            "The first three mass-matrix diagonal entries must equal mass_kg"
        )
    # JSON inertia is specified in MSS FSD at CG, including cross terms.
    return transform_matrix(matrix)


def _symmetrize_hydrodynamic_matrix(
    matrix: np.ndarray,
    name: str,
    omega: np.ndarray,
) -> np.ndarray:
    """Enforce zero-speed reciprocity and warn about excessive asymmetry."""
    matrix = np.asarray(matrix, dtype=float)
    omega = np.asarray(omega, dtype=float)
    if matrix.ndim != 3 or matrix.shape[1:] != (6, 6):
        raise ValueError(f"{name} must have shape (n_omega, 6, 6)")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must contain only finite values")
    if omega.shape != (matrix.shape[0],):
        raise ValueError(f"{name} frequencies do not match its matrix data")

    transpose = matrix.transpose(0, 2, 1)
    matrix_norm = np.linalg.norm(matrix, axis=(1, 2))
    skew_norm = np.linalg.norm(matrix - transpose, axis=(1, 2))
    relative_skew = np.zeros_like(matrix_norm)
    reference_norm = float(np.max(matrix_norm))
    np.divide(
        skew_norm,
        matrix_norm,
        out=relative_skew,
        where=matrix_norm > np.finfo(float).eps,
    )

    # A local ratio is misleading when the whole matrix is nearly zero,
    # as is common for submerged-body damping outside its wave-radiating
    # frequency band. Warn only when the skew is significant both locally
    # and relative to the largest matrix norm in the solved frequency range.
    global_skew = np.zeros_like(skew_norm)
    if reference_norm > np.finfo(float).eps:
        global_skew = skew_norm / reference_norm
    warning_candidates = (
        (relative_skew > RECIPROCITY_WARNING_THRESHOLD)
        & (global_skew > RECIPROCITY_GLOBAL_WARNING_THRESHOLD)
    )

    if np.any(warning_candidates):
        worst = int(
            np.argmax(
                np.where(warning_candidates, relative_skew, 0.0)
            )
        )
        LOG.warning(
            "%s violates reciprocity by %.2f%% locally "
            "(global skew %.3g%%) at omega=%g rad/s before "
            "symmetrization; check or refine the hydrodynamic mesh",
            name,
            100.0 * relative_skew[worst],
            100.0 * global_skew[worst],
            omega[worst],
        )

    return 0.5 * (matrix + transpose)


def run(config_path: Path) -> Path:
    """Solve one JSON configuration and return its primary vessel MAT path.

    Geometry may be CSV offsets, ShipX MGF sections, WAMIT GDF or STL panels.
    JSON CG and inertia use MSS FSD; source MGF/GDF axes are converted internally.
    Check roll stability and any independent GM_T_m reference before the BEM
    solve. Export MSS FSD hydrodynamics, hull panels and hydrostatics for either
    a surface vessel or a fully submerged vehicle (MGF is surface-only).
    """
    import capytaine as cpt
    from capytaine.post_pro import rao
    import xarray as xr

    config_path = config_path.resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    submerged = config.get("submerged")
    if not isinstance(submerged, bool):
        raise ValueError("submerged must be true or false")
    if "total_damping" in config:
        raise ValueError(
            "Rename total_damping to viscous_damping and review the values: "
            "T_126/T_1236 specify total velocity-decay targets in seconds; "
            "delta_zeta entries specify additional restored-mode damping ratios"
        )
    damping_parameters = _viscous_damping_parameters(
        config.get("viscous_damping"),
        submerged=submerged,
    )
    # Fail before solving if the installed Capytaine lacks mean-drift support.
    from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force  # noqa: F401
    from .mean_drift import compute_mean_drift, kochin_angles, angular_difference
    theta = kochin_angles(config.get("mean_drift_theta_step_deg", 0.5))
    base = config_path.parent

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    from .hull_geometry import load_hull

    vertices, faces, geometry_source = load_hull(config, base, submerged)

    body_name = config.get("body_name", config_path.parent.name)
    if not isinstance(body_name, str) or not body_name.strip():
        raise ValueError("body_name must be a nonempty string")
    output_filename = config.get("output_filename", f"{body_name}.mat")
    if (
        not isinstance(output_filename, str)
        or Path(output_filename).name != output_filename
        or Path(output_filename).suffix.lower() != ".mat"
    ):
        raise ValueError("output_filename must be a .mat filename")

    output_dir = (base / config.get("output_dir", "results")).resolve()
    rho = _positive(
        config.get("water_density_kg_m3", 1025.0),
        "water density",
    )
    gravity = _positive(
        config.get("gravity_m_s2", 9.81),
        "gravity",
    )
    lpp = None if "Lpp_m" not in config else _positive(config["Lpp_m"], "Lpp_m")
    depth_value = config.get("water_depth_m")
    water_depth = (
        np.inf
        if depth_value is None
        else _positive(depth_value, "water depth")
    )

    center_of_mass_body = transform_positions(_vector3(
        config["center_of_mass_m"],
        "center_of_mass_m",
    ))
    if submerged:
        submergence_depth = _positive(
            config.get("submergence_depth_m"),
            "submergence_depth_m",
        )
        geometry_translation = np.array([0.0, 0.0, -submergence_depth])
        if np.max(vertices[:, 2] + geometry_translation[2]) >= 0.0:
            raise ValueError(
                "submergence_depth_m must place the complete hull below "
                "the free surface"
            )
    else:
        if "submergence_depth_m" in config:
            raise ValueError(
                "submergence_depth_m is only valid when submerged is true"
            )
        submergence_depth = 0.0
        geometry_translation = np.zeros(3)

    # Public inputs are MSS FSD; geometry and CG are now in solver axes.
    # A submerged case is translated only for the free-surface BEM solve.
    # Exported MSS centers and panel geometry remain body fixed.
    body_vertices = vertices.copy()
    vertices = vertices + geometry_translation
    center_of_mass = center_of_mass_body + geometry_translation
    rotation_center = center_of_mass

    # ------------------------------------------------------------------
    # Finite frequencies
    #
    # These frequencies are used for A(w), B(w), excitation-force RAOs,
    # and motion RAOs. The value 10 rad/s is reserved for the separately
    # computed infinite-frequency radiation result below and therefore
    # cannot be configured as a finite frequency.
    # ------------------------------------------------------------------

    omega = np.asarray(config["omega_rad_s"], dtype=float)
    if (
        omega.ndim != 1
        or len(omega) == 0
        or not np.all(np.isfinite(omega))
        or np.any(omega <= 0.0)
    ):
        raise ValueError(
            "omega_rad_s must be a nonempty list of positive finite "
            "angular frequencies"
        )
    if np.any(omega >= 10.0):
        raise ValueError(
            "omega_rad_s must contain only finite frequencies strictly "
            "below 10 rad/s; 10 rad/s is reserved for the "
            "infinite-frequency result"
        )
    omega = np.unique(omega)
    omega.sort()

    periods = 2.0 * np.pi / omega

    # ------------------------------------------------------------------
    # Wave headings
    #
    # Only the symmetric half-plane is solved. capytaine_vessel.py mirrors
    # these 19 headings to the full 36-heading directional set.
    # ------------------------------------------------------------------

    headings_rad = np.deg2rad(MSS_WAVE_DIRECTIONS_DEG)

    # ------------------------------------------------------------------
    # Capytaine body and hydrostatics
    # ------------------------------------------------------------------

    mesh = cpt.Mesh(
        vertices=vertices,
        faces=faces,
        name="hull",
    )
    # Surface vessels get an internal free-surface lid to suppress irregular
    # frequencies. A submerged vehicle has no waterplane to lid.
    if submerged:
        lid = None
    elif "max_panel_radius_m" in config:
        # Refinement preserves the hull surface. Generate the same interior
        # waterplane from its smaller source mesh to avoid the lid reader's
        # temporary point-by-panel arrays scaling with refined hull size.
        lid_config = dict(config)
        lid_radius = lid_config.pop("max_panel_radius_m")
        lid_vertices, lid_faces, _ = load_hull(lid_config, base, False)
        lid_hull = cpt.Mesh(vertices=lid_vertices, faces=lid_faces, name="lid source")
        from .hull_geometry import bounded_lid
        lid = bounded_lid(lid_hull, lid_radius)
    else:
        lid = mesh.generate_lid()

    mass = _positive(config["mass_kg"], "mass_kg")
    mass_matrix = _mass_matrix(config, base, mass)

    reflection = config.get("reflection_symmetry", False)
    if not isinstance(reflection, bool):
        raise ValueError("reflection_symmetry must be true or false")
    if reflection:
        from .geometry_format import geometry_format
        from .hull_geometry import reflection_mesh
        filename = config.get("geometry_file", config.get("offset_points_csv"))
        if geometry_format(base / filename) not in {".csv", ".mgf"}:
            raise ValueError("reflection_symmetry currently requires symmetric CSV/MGF sections")
        parity = np.array([1., -1., 1., -1., 1., -1.])
        if not np.isclose(center_of_mass[1], 0., atol=1e-10, rtol=0):
            raise ValueError("reflection_symmetry requires CG on the centerline")
        if not np.allclose(mass_matrix, parity[:, None]*mass_matrix*parity[None, :],
                           rtol=1e-12, atol=1e-10*np.max(np.abs(mass_matrix))):
            raise ValueError("reflection_symmetry requires reflection-symmetric inertia")
        mesh = reflection_mesh(mesh)
        if lid is not None:
            lid = reflection_mesh(lid)
        # Plane clipping can subdivide centerline-crossing panels. Export
        # the actual full mirrored solver mesh, rather than the input faces.
        body_vertices = mesh.vertices.copy() - geometry_translation
        faces = mesh.faces.copy()

    body = cpt.FloatingBody(
        mesh=mesh,
        lid_mesh=lid,
        dofs=cpt.rigid_body_dofs(rotation_center=rotation_center),
        center_of_mass=center_of_mass,
        mass=mass,
        name=body_name,
    )
    if tuple(body.dofs) != DOFS:
        raise ValueError(
            f"Expected six rigid-body DOFs in this order: {DOFS}"
        )
    body.inertia_matrix = body.add_dofs_labels_to_matrix(mass_matrix)

    immersed = body.immersed_part(water_depth=water_depth)
    displaced_mass = float(immersed.disp_mass(rho=rho))
    if not np.isfinite(displaced_mass) or displaced_mass <= 0.0:
        raise ValueError(
            "The hull does not enclose a positive submerged volume"
        )
    relative_mass_error = abs(mass - displaced_mass) / displaced_mass
    if relative_mass_error > 0.01:
        LOG.warning(
            "Specified mass differs from displaced mass by %.2f%%; "
            "check flotation equilibrium",
            100.0 * relative_mass_error,
        )

    if submerged:
        x_f_mesh = 0.0
        center_of_flotation_mss = None
    else:
        waterplane_area = float(immersed.waterplane_area)
        if not np.isfinite(waterplane_area) or waterplane_area <= 0.0:
            raise ValueError("The surface vessel has no positive waterplane area")
        waterplane_x = immersed.mesh.quadrature_points[0][:, :, 0]
        center_of_flotation_x = float(
            immersed.mesh.waterplane_integral(waterplane_x)
            / waterplane_area
        )
        # The mesh x axis points aft, whereas MSS x points forward.
        x_f_mesh = center_of_flotation_x - center_of_mass[0]
        center_of_flotation_mss = [
            -center_of_flotation_x,
            0.0,
            0.0,
        ]

    hydrostatic_stiffness = body.compute_hydrostatic_stiffness(
        rho=rho,
        g=gravity,
    )
    hydrostatic_matrix = np.asarray(
        hydrostatic_stiffness,
        dtype=float,
    ).copy()
    gm_t = float(hydrostatic_matrix[3, 3] / (rho * gravity * immersed.volume))
    _check_transverse_stability(config, gm_t)
    if submerged:
        # A freely submerged body has no restoring force or moment in surge,
        # sway, heave, or yaw.
        free_dofs = (0, 1, 2, 5)
        hydrostatic_matrix[list(free_dofs), :] = 0.0
        hydrostatic_matrix[:, list(free_dofs)] = 0.0

        # For a freely submerged constant-volume vehicle, gravity and
        # buoyancy provide independent roll and pitch restoring moments.
        # Roll-pitch hydrostatic coupling is identically zero.
        hydrostatic_matrix[3, 4] = 0.0
        hydrostatic_matrix[4, 3] = 0.0
    else:
        # Implement Eq. (4.28) in MSS FSD conventions. Capytaine supplies
        # the diagonal stiffness about CG in mesh axes. Remove its x_F shift
        # from pitch, then explicitly apply the equivalent CF-to-CG screw
        # transformation. The mesh-frame construction is converted to FSD
        # later, where x_F = -x_f_mesh.
        g33 = hydrostatic_matrix[2, 2]
        g44_cf = hydrostatic_matrix[3, 3]
        g55_cf = hydrostatic_matrix[4, 4] - g33 * x_f_mesh**2

        hydrostatic_matrix.fill(0.0)
        hydrostatic_matrix[2, 2] = g33
        hydrostatic_matrix[2, 4] = -g33 * x_f_mesh
        hydrostatic_matrix[3, 3] = g44_cf
        hydrostatic_matrix[4, 2] = -g33 * x_f_mesh
        hydrostatic_matrix[4, 4] = g55_cf + g33 * x_f_mesh**2

    body.hydrostatic_stiffness = body.add_dofs_labels_to_matrix(
        hydrostatic_matrix
    )
    immersed = body.immersed_part(water_depth=water_depth)

    # ------------------------------------------------------------------
    # Finite-frequency radiation and diffraction below 10 rad/s
    # ------------------------------------------------------------------

    test_matrix = xr.Dataset(
        coords={
            "omega": omega,
            "wave_direction": headings_rad,
            "radiating_dof": list(DOFS),
            "water_depth": [water_depth],
            "forward_speed": [0.0],
            "rho": [rho],
            "g": [gravity],
            "theta": theta,
        }
    )

    solver = cpt.BEMSolver()
    dataset = solver.fill_dataset(
        test_matrix,
        immersed,
        hydrostatics=True,
    )
    dataset["inertia_matrix"] = body.inertia_matrix
    dataset["hydrostatic_stiffness"] = body.hydrostatic_stiffness
    if "period" not in dataset.coords:
        dataset = dataset.assign_coords(period=("omega", periods))

    dataset.attrs["geometry_source"] = geometry_source
    dataset.attrs["coordinate_origin"] = (
        "body-fixed origin"
        if submerged
        else "midships, centerline, design waterline"
    )
    dataset.attrs["submerged"] = submerged
    dataset.attrs["submergence_depth_m"] = submergence_depth

    # ------------------------------------------------------------------
    # Limiting-frequency radiation
    #
    # These results add A(0), B(0), A(infinity), and B(infinity) to the
    # hydrodynamic matrices. The infinite-frequency result is exported at
    # the conventional practical endpoint of 10 rad/s. No limiting-frequency
    # RAOs are calculated.
    # ------------------------------------------------------------------

    if not np.isinf(water_depth):
        raise ValueError(
            "The limiting-frequency Capytaine radiation calculation "
            "requires infinite water depth"
        )

    limit_matrix = xr.Dataset(
        coords={
            "omega": [0.0, np.inf],
            "radiating_dof": list(DOFS),
            "water_depth": [water_depth],
            "forward_speed": [0.0],
            "rho": [rho],
            "g": [gravity],
        }
    )
    limit_dataset = solver.fill_dataset(
        limit_matrix,
        immersed,
        hydrostatics=False,
    )

    # ------------------------------------------------------------------
    # Convert all exported quantities to MSS FSD at CG
    # ------------------------------------------------------------------

    from .output_coordinates import (
        DOF_SIGNS,
        POSITION_SIGNS,
        to_mss_fsd,
    )

    output_dataset = to_mss_fsd(dataset)

    A0 = np.asarray(
        limit_dataset.added_mass.sel(omega=0.0).transpose(
            "influenced_dof",
            "radiating_dof",
        ).values,
        dtype=float,
    )
    B0 = np.asarray(
        limit_dataset.radiation_damping.sel(omega=0.0).transpose(
            "influenced_dof",
            "radiating_dof",
        ).values,
        dtype=float,
    )
    A0 = DOF_SIGNS[:, None] * A0 * DOF_SIGNS[None, :]
    B0 = DOF_SIGNS[:, None] * B0 * DOF_SIGNS[None, :]

    Ainf = np.asarray(
        limit_dataset.added_mass.sel(omega=np.inf).transpose(
            "influenced_dof",
            "radiating_dof",
        ).values,
        dtype=float,
    )
    Binf = np.asarray(
        limit_dataset.radiation_damping.sel(omega=np.inf).transpose(
            "influenced_dof",
            "radiating_dof",
        ).values,
        dtype=float,
    )
    Ainf = DOF_SIGNS[:, None] * Ainf * DOF_SIGNS[None, :]
    Binf = DOF_SIGNS[:, None] * Binf * DOF_SIGNS[None, :]

    # At zero forward speed, radiation added mass and damping satisfy the
    # reciprocity relation. Discretization errors can leave a skew component,
    # so diagnose it and use the symmetric coefficients consistently for
    # motion RAOs and vessel export.
    all_omega = np.concatenate(([0.0], omega, [np.inf]))
    added_mass = output_dataset.added_mass.transpose(
        "omega", "influenced_dof", "radiating_dof"
    )
    added_mass_values = _symmetrize_hydrodynamic_matrix(
        np.concatenate(
            (
                A0[None, :, :],
                np.asarray(added_mass.values, dtype=float),
                Ainf[None, :, :],
            ),
            axis=0,
        ),
        "Added-mass matrix",
        all_omega,
    )
    A0 = added_mass_values[0]
    Ainf = added_mass_values[-1]
    output_dataset["added_mass"] = added_mass.copy(
        data=added_mass_values[1:-1]
    )

    radiation_damping = output_dataset.radiation_damping.transpose(
        "omega", "influenced_dof", "radiating_dof"
    )
    damping_values = _symmetrize_hydrodynamic_matrix(
        np.concatenate(
            (
                B0[None, :, :],
                np.asarray(radiation_damping.values, dtype=float),
                Binf[None, :, :],
            ),
            axis=0,
        ),
        "Radiation-damping matrix",
        all_omega,
    )
    B0 = damping_values[0]
    Binf = damping_values[-1]
    output_dataset["radiation_damping"] = radiation_damping.copy(
        data=damping_values[1:-1]
    )

    # RAOs use potential-flow radiation damping only. The one viscous-damping
    # matrix is formed later by the local MATLAB processing helper from the
    # time constants and ratio increments exported under vessel.powerBased.
    output_dataset["motion_rao"] = rao(output_dataset)

    LOG.info("Computing diagonal mean-drift coefficients (surge, sway, yaw)")
    # Use exactly the exported, symmetrized potential-flow motion response.
    # Undo the FSD axis/heading conversion for the solver's Kochin functions.
    solver_headings = np.asarray(dataset.wave_direction)
    mss_headings = np.mod(np.pi - solver_headings, 2*np.pi)
    indices = [int(np.argmin(abs(output_dataset.wave_direction.values - beta)))
               for beta in mss_headings]
    motion_solver = dataset.excitation_force.rename(influenced_dof="radiating_dof").transpose(
        "omega", "wave_direction", "radiating_dof").copy(
        data=output_dataset.motion_rao.transpose("omega", "wave_direction", "radiating_dof").values[:, indices, :]
             * DOF_SIGNS)
    drift = compute_mean_drift(dataset, motion_solver, center_of_mass)
    angular_error = angular_difference(dataset, motion_solver, center_of_mass, drift)
    if angular_error is not None and np.any(angular_error > .02):
        LOG.warning("Mean-drift angular-grid norm differences (surge/sway/yaw): %s. "
                    "Reduce mean_drift_theta_step_deg and check convergence.",
                    ", ".join(f"{100*e:.2g}%" for e in angular_error))
    drift = drift * xr.DataArray(DOF_SIGNS[[0, 1, 5]], dims="drift_dof",
                                 coords={"drift_dof": drift.drift_dof})
    output_dataset["mean_drift_force"] = drift.assign_coords(
        wave_direction=mss_headings).sortby("wave_direction")

    hydrostatics = {
        "GM_T_m": gm_t,
        "volume_m3": float(immersed.volume),
        "displaced_mass_kg": displaced_mass,
        "body_mass_kg": mass,
        "center_of_buoyancy_m": (
            POSITION_SIGNS
            * (
                np.asarray(immersed.center_of_buoyancy)
                - geometry_translation
            )
        ).tolist(),
        "waterplane_area_m2": float(immersed.waterplane_area),
        "center_of_mass_m": (
            POSITION_SIGNS * center_of_mass_body
        ).tolist(),
        "rotation_center_m": (
            POSITION_SIGNS * center_of_mass_body
        ).tolist(),
        "inertia_matrix_SI": (
            DOF_SIGNS[:, None]
            * mass_matrix
            * DOF_SIGNS[None, :]
        ).tolist(),
        "hydrostatic_stiffness_SI": (
            DOF_SIGNS[:, None]
            * np.asarray(body.hydrostatic_stiffness)
            * DOF_SIGNS[None, :]
        ).tolist(),
        "coordinate_system": "MSS FSD: x forward, y starboard, z down",
        "geometry_coordinate_origin": "body-fixed origin",
        "matrix_reference": "CG",
        "submerged": submerged,
        "submergence_depth_m": submergence_depth,
        "mean_drift": {
            "method": "Capytaine far-field diagonal mean drift",
            "reference": "CG; MSS FSD",
            "normalization": "load per incident wave amplitude squared",
            "motion_damping": "potential radiation only",
            "theta_step_deg": float(np.rad2deg(theta[1] - theta[0])),
            "coarser_grid_relative_norm_difference": (
                None if angular_error is None else angular_error.tolist()),
        },
    }
    if center_of_flotation_mss is not None:
        hydrostatics["center_of_flotation_m"] = center_of_flotation_mss
    if "demihull_spacing_m" in config:
        hydrostatics["demihull_spacing_m"] = float(config["demihull_spacing_m"])

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "hydrostatics.json").open(
        "w",
        encoding="utf-8",
    ) as stream:
        json.dump(hydrostatics, stream, indent=2)

    from .capytaine_vessel import write_vessel

    write_vessel(
        output_dataset,
        hydrostatics,
        body_vertices,
        faces,
        body.name,
        rho,
        gravity,
        base / output_filename,
        zero_added_mass=A0,
        zero_radiation_damping=B0,
        infinite_added_mass=Ainf,
        infinite_radiation_damping=Binf,
        submerged=submerged,
        Lpp_m=lpp,
        supporting_output_dir=output_dir,
        **damping_parameters,
    )

    return base / output_filename
