"""Zero-speed, six-DOF Capytaine solver with MATLAB vessel export."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np


DOFS = ("Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw")
LOG = logging.getLogger(__name__)


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


def _mass_matrix(
    config: dict,
    base: Path,
    mass: float,
) -> np.ndarray:
    """Read the configured inertia matrix or build its diagonal form."""
    matrix_file = config.get("inertia_matrix_csv")
    if matrix_file:
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
    return matrix


def run(config_path: Path) -> Path:
    """Run Capytaine and write ``capytaineTestShip.mat``."""
    import capytaine as cpt
    from capytaine.post_pro import rao
    import xarray as xr

    config_path = config_path.resolve()
    with config_path.open(encoding="utf-8") as stream:
        config = json.load(stream)
    if "viscous_damping" in config:
        raise ValueError(
            "Rename viscous_damping to total_damping in the JSON configuration"
        )
    base = config_path.parent

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    shipx_archive = config.get("shipx_archive")
    if shipx_archive:
        from shipx_geometry import (
            panels_from_shipx_sections,
            read_shipx_sections,
            section_volume_m3,
        )

        archive_path = (base / shipx_archive).resolve()
        sections = read_shipx_sections(
            archive_path,
            float(config["draft_m"]),
        )
        vertices, faces = panels_from_shipx_sections(
            sections,
            int(config.get("shipx_samples_per_section", 20)),
        )
        LOG.info(
            "ShipX source sections: %d; offset-integrated volume: %.3f m3",
            len(sections),
            section_volume_m3(sections),
        )
    else:
        from .offset_geometry import (
            panels_from_sections,
            read_offset_sections,
        )

        offsets_path = (base / config["offset_points_csv"]).resolve()
        if not offsets_path.is_file():
            raise FileNotFoundError(
                f"Hull offset points not found: {offsets_path}"
            )
        sections = read_offset_sections(offsets_path)
        number_of_stations = config.get("number_of_stations")
        vertices, faces = panels_from_sections(
            sections,
            samples_per_section=int(config.get("samples_per_section", 20)),
            number_of_stations=(
                None
                if number_of_stations is None
                else int(number_of_stations)
            ),
        )

    if shipx_archive:
        raise ValueError(
            "The independent MSS vessel export requires offset_points_csv, "
            "not shipx_archive"
        )

    output_dir = (base / config.get("output_dir", "results")).resolve()
    rho = _positive(
        config.get("water_density_kg_m3", 1025.0),
        "water density",
    )
    gravity = _positive(
        config.get("gravity_m_s2", 9.81),
        "gravity",
    )
    depth_value = config.get("water_depth_m")
    water_depth = (
        np.inf
        if depth_value is None
        else _positive(depth_value, "water depth")
    )

    rotation_center = _vector3(
        config["rotation_center_m"],
        "rotation_center_m",
    )
    center_of_mass = _vector3(
        config["center_of_mass_m"],
        "center_of_mass_m",
    )
    if not np.allclose(rotation_center, center_of_mass):
        raise ValueError(
            "CG-referenced output requires rotation_center_m = "
            "center_of_mass_m"
        )

    # ------------------------------------------------------------------
    # Finite frequencies
    #
    # These frequencies are used for A(w), B(w), excitation-force RAOs,
    # and motion RAOs. A configured 10 rad/s value is reserved for the
    # separately computed infinite-frequency radiation result below.
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
    omega[np.isclose(omega, 10.0)] = 10.0
    if np.any(omega > 10.0):
        raise ValueError(
            "omega_rad_s must not contain frequencies above 10 rad/s"
        )

    # The MSS hydrodynamic endpoint at 10 rad/s contains the omega=infinity
    # radiation solution, so it must not also be solved as a finite-frequency
    # RAO point.
    omega = omega[omega < 10.0]
    if len(omega) == 0:
        raise ValueError(
            "omega_rad_s must contain at least one frequency below 10 rad/s"
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

    headings_deg = np.asarray(
        config["wave_directions_deg"],
        dtype=float,
    )
    expected_headings = np.arange(0.0, 181.0, 10.0)
    if headings_deg.shape != expected_headings.shape or not np.allclose(
        headings_deg,
        expected_headings,
    ):
        raise ValueError(
            "wave_directions_deg must be 0, 10, ..., 180 degrees"
        )
    headings_rad = np.deg2rad(headings_deg)

    # ------------------------------------------------------------------
    # Capytaine body and hydrostatics
    # ------------------------------------------------------------------

    mesh = cpt.Mesh(
        vertices=vertices,
        faces=faces,
        name="hull_from_offsets",
    )
    lid = (
        mesh.generate_lid()
        if config.get("generate_lid", True)
        else None
    )

    mass = _positive(config["mass_kg"], "mass_kg")
    mass_matrix = _mass_matrix(config, base, mass)

    body = cpt.FloatingBody(
        mesh=mesh,
        lid_mesh=lid,
        dofs=cpt.rigid_body_dofs(rotation_center=rotation_center),
        center_of_mass=center_of_mass,
        mass=mass,
        name=config.get("body_name", "capytaineTestShip"),
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
            "The offsets do not enclose a positive submerged volume"
        )
    relative_mass_error = abs(mass - displaced_mass) / displaced_mass
    if relative_mass_error > 0.01:
        LOG.warning(
            "Specified mass differs from displaced mass by %.2f%%; "
            "check flotation equilibrium",
            100.0 * relative_mass_error,
        )

    hydrostatic_stiffness = body.compute_hydrostatic_stiffness(
        rho=rho,
        g=gravity,
    )
    hydrostatic_matrix = np.asarray(
        hydrostatic_stiffness,
        dtype=float,
    ).copy()
    # A freely floating body has no hydrostatic restoring force or moment
    # in surge, sway, or yaw. Remove mesh-integration residuals in those
    # rows and columns before computing the motion RAOs.
    free_dofs = (0, 1, 5)
    hydrostatic_matrix[list(free_dofs), :] = 0.0
    hydrostatic_matrix[:, list(free_dofs)] = 0.0
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

    dataset.attrs["geometry_source"] = "independent offset points"

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

    # Bv is diagonal at the CG and has the same entries in both coordinate
    # frames. Compute it from the FSD coefficients that will be exported.
    from .viscous_damping import diagonal_viscous_damping

    viscous_damping = diagonal_viscous_damping(
        config.get("total_damping"),
        omega,
        np.asarray(output_dataset.inertia_matrix.values, dtype=float),
        A0,
        B0,
        np.asarray(output_dataset.added_mass.transpose(
            "omega", "influenced_dof", "radiating_dof"
        ).values, dtype=float),
        np.asarray(output_dataset.radiation_damping.transpose(
            "omega", "influenced_dof", "radiating_dof"
        ).values, dtype=float),
        np.asarray(output_dataset.hydrostatic_stiffness.values, dtype=float),
    )

    # Viscous damping changes the motion response, not the excitation force.
    # RAOs are still evaluated only at positive finite frequencies.
    dissipation = xr.DataArray(
        viscous_damping,
        dims=("influenced_dof", "radiating_dof"),
        coords={
            "influenced_dof": output_dataset.influenced_dof,
            "radiating_dof": output_dataset.radiating_dof,
        },
    )
    output_dataset["motion_rao"] = rao(output_dataset, dissipation=dissipation)

    hydrostatics = {
        "volume_m3": float(immersed.volume),
        "displaced_mass_kg": displaced_mass,
        "body_mass_kg": mass,
        "center_of_buoyancy_m": (
            POSITION_SIGNS
            * np.asarray(immersed.center_of_buoyancy)
        ).tolist(),
        "waterplane_area_m2": float(immersed.waterplane_area),
        "center_of_mass_m": (
            POSITION_SIGNS * np.asarray(center_of_mass)
        ).tolist(),
        "rotation_center_m": (
            POSITION_SIGNS * np.asarray(rotation_center)
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
        "matrix_reference": "CG",
    }

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
        vertices,
        faces,
        body.name,
        rho,
        gravity,
        output_dir / "capytaineTestShip.mat",
        zero_added_mass=A0,
        zero_radiation_damping=B0,
        infinite_added_mass=Ainf,
        infinite_radiation_damping=Binf,
        viscous_damping_matrix=viscous_damping,
    )

    return output_dir
