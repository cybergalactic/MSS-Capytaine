"""Build a wetted hull mesh from an independent station offset table."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def read_offset_sections(path: Path) -> list[tuple[float, np.ndarray]]:
    """Read x, z, half breadth in mesh axes: x aft, z up from waterline.

    Rows within a station run from keel to waterline. Repeated z values are
    allowed so a flat bottom can be represented without losing its corners.
    """
    sections: dict[float, list[tuple[float, float]]] = {}

    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        required = {"x_m", "z_m", "half_breadth_m"}

        if not required.issubset(reader.fieldnames or []):
            raise ValueError(
                f"{path} must have columns x_m,z_m,half_breadth_m"
            )

        for line, row in enumerate(reader, start=2):
            try:
                x, z, half_breadth = (
                    float(row[name])
                    for name in ("x_m", "z_m", "half_breadth_m")
                )
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid offset point on line {line}"
                ) from exc

            if (
                not np.all(np.isfinite((x, z, half_breadth)))
                or z > 1e-9
                or half_breadth < 0
            ):
                raise ValueError(
                    f"Line {line} requires finite x, z<=0, "
                    "half breadth>=0"
                )

            sections.setdefault(x, []).append((half_breadth, z))

    if len(sections) < 3:
        raise ValueError(
            "At least three longitudinal stations are required"
        )

    ordered = []

    for x in sorted(sections):
        profile = np.asarray(sections[x], dtype=float)

        if (
            len(profile) < 2
            or np.any(np.diff(profile[:, 1]) < -1e-9)
        ):
            raise ValueError(
                f"Station x={x} must contain "
                "keel-to-waterline ordered points"
            )

        if (
            profile[0, 1] >= 0
            or not np.isclose(profile[-1, 1], 0.0, atol=1e-9)
        ):
            raise ValueError(
                f"Station x={x} must start below and end at the waterline"
            )

        if not np.isclose(profile[0, 0], 0.0, atol=1e-9):
            raise ValueError(
                f"Station x={x} must start at the centerline keel"
            )

        if np.any(
            np.linalg.norm(np.diff(profile, axis=0), axis=1) <= 1e-10
        ):
            raise ValueError(
                f"Station x={x} has duplicate consecutive points"
            )

        ordered.append((x, profile))

    if not any(
        np.any(profile[:, 0] > 0)
        for _, profile in ordered
    ):
        raise ValueError("Offset table has no positive half breadths")

    return ordered


def panels_from_sections(
    sections: list[tuple[float, np.ndarray]],
    samples_per_section: int = 20,
    number_of_stations: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Connect arc-length-resampled half sections into a symmetric hull.

    ``samples_per_section`` controls the transverse resolution from keel
    to waterline on each half section.

    ``number_of_stations`` controls the longitudinal resolution. If omitted,
    the original station locations from the offset table are retained.
    """

    if samples_per_section < 3:
        raise ValueError(
            "samples_per_section must be at least 3"
        )

    # ---------------------------------------------------------------
    # Transverse interpolation
    # ---------------------------------------------------------------

    xs = np.asarray([x for x, _ in sections], dtype=float)
    profiles = []

    for _, profile in sections:
        along = np.r_[
            0.0,
            np.cumsum(
                np.linalg.norm(
                    np.diff(profile, axis=0),
                    axis=1,
                )
            ),
        ]

        target = np.linspace(
            0.0,
            along[-1],
            samples_per_section,
        )

        profiles.append(
            np.column_stack(
                (
                    np.interp(
                        target,
                        along,
                        profile[:, 0],
                    ),
                    np.interp(
                        target,
                        along,
                        profile[:, 1],
                    ),
                )
            )
        )

    profiles = np.asarray(profiles, dtype=float)

    # ---------------------------------------------------------------
    # Longitudinal interpolation
    # ---------------------------------------------------------------

    if number_of_stations is not None:
        if number_of_stations < len(xs):
            raise ValueError(
                "number_of_stations must be at least the number "
                "of stations in the offset table"
            )

        xs_new = np.linspace(
            xs[0],
            xs[-1],
            number_of_stations,
        )

        profiles_new = np.empty(
            (
                number_of_stations,
                samples_per_section,
                2,
            ),
            dtype=float,
        )

        for j in range(samples_per_section):
            profiles_new[:, j, 0] = np.interp(
                xs_new,
                xs,
                profiles[:, j, 0],
            )

            profiles_new[:, j, 1] = np.interp(
                xs_new,
                xs,
                profiles[:, j, 1],
            )

        xs = xs_new
        profiles = profiles_new

    # ---------------------------------------------------------------
    # Generate vertices and panels
    # ---------------------------------------------------------------

    vertices: list[tuple[float, float, float]] = []
    vertex_ids: dict[tuple[float, float, float], int] = {}
    faces: list[tuple[int, int, int, int]] = []

    def vertex(
        point: tuple[float, float, float],
    ) -> int:
        if point not in vertex_ids:
            vertex_ids[point] = len(vertices)
            vertices.append(point)

        return vertex_ids[point]

    def triangle(a, b, c) -> None:
        if (
            np.linalg.norm(
                np.cross(
                    np.subtract(b, a),
                    np.subtract(c, a),
                )
            )
            > 1e-10
        ):
            ia, ib, ic = (
                vertex(point)
                for point in (a, b, c)
            )

            faces.append((ia, ib, ic, ic))

    def quadrangle(a, b, c, d) -> None:
        triangle(a, b, c)
        triangle(a, c, d)

    # ---------------------------------------------------------------
    # Hull sides
    # ---------------------------------------------------------------

    for i in range(len(xs) - 1):
        for j in range(samples_per_section - 1):
            for sign in (1, -1):

                a = (
                    xs[i],
                    sign * profiles[i][j, 0],
                    profiles[i][j, 1],
                )

                b = (
                    xs[i],
                    sign * profiles[i][j + 1, 0],
                    profiles[i][j + 1, 1],
                )

                c = (
                    xs[i + 1],
                    sign * profiles[i + 1][j + 1, 0],
                    profiles[i + 1][j + 1, 1],
                )

                d = (
                    xs[i + 1],
                    sign * profiles[i + 1][j, 0],
                    profiles[i + 1][j, 1],
                )

                if sign == 1:
                    quadrangle(a, b, c, d)
                else:
                    triangle(a, c, b)
                    triangle(a, d, c)

    # ---------------------------------------------------------------
    # Close bow and stern
    # ---------------------------------------------------------------

    for i, end in (
        (0, -1),
        (len(xs) - 1, 1),
    ):
        for j in range(samples_per_section - 1):

            low_port = (
                xs[i],
                -profiles[i][j, 0],
                profiles[i][j, 1],
            )

            low_starboard = (
                xs[i],
                profiles[i][j, 0],
                profiles[i][j, 1],
            )

            high_starboard = (
                xs[i],
                profiles[i][j + 1, 0],
                profiles[i][j + 1, 1],
            )

            high_port = (
                xs[i],
                -profiles[i][j + 1, 0],
                profiles[i][j + 1, 1],
            )

            if end == -1:
                quadrangle(
                    low_port,
                    high_port,
                    high_starboard,
                    low_starboard,
                )
            else:
                quadrangle(
                    low_port,
                    low_starboard,
                    high_starboard,
                    high_port,
                )

    if not faces:
        raise ValueError(
            "No nondegenerate panels could be generated "
            "from the offset table"
        )

    return (
        np.asarray(vertices, dtype=float),
        np.asarray(faces, dtype=int),
    )