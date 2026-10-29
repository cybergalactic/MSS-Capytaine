"""Build a wetted hull mesh from an independent station offset table."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def refine_large_panels(vertices, faces, max_radius):
    """Bisect oversized triangles without moving the source hull surface.

    Reflections get identical refinement. Shared-edge midpoints are reused;
    neighbouring panels may meet at a subdivision of an unchanged edge.
    """
    if (isinstance(max_radius, bool) or not isinstance(max_radius, (int, float))
            or not np.isfinite(max_radius) or max_radius <= 0):
        raise ValueError("max_panel_radius_m must be positive and finite")
    points = list(map(tuple, vertices))
    midpoint_ids = {}
    pending = [tuple(map(int, face[:3])) for face in faces]
    refined = []
    while pending:
        triangle = pending.pop()
        p = np.asarray([points[i] for i in triangle])
        if np.linalg.norm(p - p.mean(axis=0), axis=1).max() <= max_radius:
            refined.append((*triangle, triangle[-1]))
            continue
        edges = ((0, 1, 2), (1, 2, 0), (2, 0, 1))
        a, b, c = max(edges, key=lambda edge: np.linalg.norm(p[edge[0]]-p[edge[1]]))
        ia, ib, ic = triangle[a], triangle[b], triangle[c]
        edge = tuple(sorted((ia, ib)))
        if edge not in midpoint_ids:
            midpoint_ids[edge] = len(points)
            points.append(tuple((p[a]+p[b])/2))
        midpoint = midpoint_ids[edge]
        pending.extend(((ia, midpoint, ic), (midpoint, ib, ic)))
    return np.asarray(points), np.asarray(refined, dtype=int)


def read_offset_sections(
    path: Path,
    *,
    submerged: bool = False,
) -> list[tuple[float, np.ndarray]]:
    """Read MSS FSD offsets, converting to internal x-aft, z-up sections.

    Input x points forward and z down. Rows run bottom to top (decreasing z).
    Surface-vessel sections end
    at the design waterline (z=0). Submerged sections are closed profiles with
    zero half breadth at both the bottom and top. Repeated z values are allowed
    so a flat bottom can be represented without losing its corners.
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
                or half_breadth < 0
                or (not submerged and z < -1e-9)
            ):
                raise ValueError(
                    f"Line {line} requires finite x and z, "
                    "half breadth>=0, and FSD z>=0 for a surface vessel"
                )

            sections.setdefault(-x, []).append((half_breadth, -z))

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
                f"Station x={-x} must contain bottom-to-top ordered points (decreasing FSD z)"
            )

        if submerged:
            if (
                profile[0, 1] >= profile[-1, 1]
                or not np.isclose(profile[-1, 0], 0.0, atol=1e-9)
            ):
                raise ValueError(
                    f"Submerged station x={-x} must run bottom-to-top "
                    "and end on the centerline"
                )
        elif (
            profile[0, 1] >= 0
            or not np.isclose(profile[-1, 1], 0.0, atol=1e-9)
        ):
            raise ValueError(
                f"Station x={-x} must start below and end at the waterline"
            )

        if not np.isclose(profile[0, 0], 0.0, atol=1e-9):
            raise ValueError(
                f"Station x={-x} must start at the centerline keel"
            )

        if np.any(
            np.linalg.norm(np.diff(profile, axis=0), axis=1) <= 1e-10
        ):
            raise ValueError(
                f"Station x={-x} has duplicate consecutive points"
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

    # A transom strip can span the entire beam. Subdivide it in both
    # directions to match the surrounding hull resolution, including
    # triangular strips whose lower width collapses to the centerline.
    cap_spacing = max(
        float(np.median(np.diff(xs))),
        float(np.median(np.linalg.norm(np.diff(profiles, axis=1), axis=2))),
    )

    def between(a, b, fraction):
        return tuple((1-fraction)*np.asarray(a) + fraction*np.asarray(b))

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

            width = max(2*profiles[i][j, 0], 2*profiles[i][j+1, 0])
            side_length = np.linalg.norm(np.subtract(high_port, low_port))
            columns = max(1, int(np.ceil(width / cap_spacing)))
            rows = max(1, int(np.ceil(side_length / cap_spacing)))
            for row in range(rows):
                left0 = between(low_port, high_port, row/rows)
                right0 = between(low_starboard, high_starboard, row/rows)
                left1 = between(low_port, high_port, (row+1)/rows)
                right1 = between(low_starboard, high_starboard, (row+1)/rows)
                for column in range(columns):
                    a = between(left0, right0, column/columns)
                    b = between(left1, right1, column/columns)
                    c = between(left1, right1, (column+1)/columns)
                    d = between(left0, right0, (column+1)/columns)
                    quadrangle(a, b, c, d) if end == -1 else quadrangle(a, d, c, b)

    if not faces:
        raise ValueError(
            "No nondegenerate panels could be generated "
            "from the offset table"
        )

    return (
        np.asarray(vertices, dtype=float),
        np.asarray(faces, dtype=int),
    )
