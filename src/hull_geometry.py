"""Generate panels from CSV/MGF sections or import GDF/STL panel meshes.

Recognize the filename extension, or inspect contents when no known extension
is present. Public coordinates use MSS FSD; solver-axis conversion is internal.
"""

from pathlib import Path

import numpy as np

from .geometry_format import geometry_format


def reflection_mesh(mesh):
    """Represent an xOz-symmetric hull or lid by one half and its mirror."""
    import capytaine as cpt
    half = mesh.clipped(origin=(0, 0, 0), normal=(0, 1, 0))
    return cpt.ReflectionSymmetricMesh(half, plane="xOz", name=mesh.name)


def bounded_lid(mesh, max_radius):
    """Generate Capytaine's interior lid with bounded temporary arrays.

    Classify candidate waterplane vertices in batches rather than allocating
    one array for every candidate vertex against every hull panel at once.
    """
    import capytaine as cpt
    from capytaine.meshes.predefined.rectangles import mesh_rectangle

    hull = mesh.clipped(origin=(0, 0, 0), normal=(0, 0, 1))
    if hull.nb_faces == 0:
        return cpt.Mesh(None, None, name="lid")
    low, high = hull.vertices.min(axis=0), hull.vertices.max(axis=0)
    span = high-low
    candidate = mesh_rectangle(
        size=(max(max_radius/2, 1.1*span[1]), max(max_radius/2, 1.1*span[0])),
        faces_max_radius=max_radius,
        center=((low[0]+high[0])/2, (low[1]+high[1])/2, 0), normal=(0, 0, -1),
    )
    panels = hull.vertices[hull.faces, :2]
    edges = panels[:, [1, 2, 3, 0], :]-panels
    inside = np.zeros(len(candidate.vertices), dtype=bool)
    for start in range(0, len(inside), 128):
        stop = min(start+128, len(inside))
        local = candidate.vertices[start:stop, None, None, :2]-panels[None]
        side = local[..., 0]*edges[None, ..., 1]-local[..., 1]*edges[None, ..., 0]
        above = np.all(side <= 0, axis=-1) | np.all(side >= 0, axis=-1)
        inside[start:stop] = above.sum(axis=-1) % 2 == 1
    faces = candidate.faces[np.all(inside[candidate.faces], axis=1)]
    if not len(faces):
        return cpt.Mesh(None, None, name="lid")
    return cpt.Mesh(candidate.vertices, faces, name="lid", auto_check=False)


def load_hull(config: dict, base: Path, submerged: bool):
    """Return panels in solver axes: x aft, y starboard, z up, in metres.

    CSV/STL use MSS FSD. MGF/GDF retain their standard file conventions.
    JSON centers always use MSS FSD; all conversions are internal.
    CSV demihull_spacing_m duplicates one locally symmetric hull at y=+/-s/2.
    Both components belong to the same rigid-body hydrodynamic solve.
    """
    if config.get("shipx_archive"):
        raise ValueError("Extract the MGF and select it with geometry_file")
    if "geometry_file" in config and "offset_points_csv" in config:
        raise ValueError("Provide only geometry_file or legacy offset_points_csv")
    filename = config.get("geometry_file", config.get("offset_points_csv"))
    if not isinstance(filename, str) or not filename.strip():
        raise ValueError("geometry_file must name a CSV, MGF, GDF or STL hull file")
    path = (base / filename).resolve()
    extension = geometry_format(path)
    if not path.is_file():
        raise FileNotFoundError(f"Hull geometry not found: {path}")
    spacing = config.get("demihull_spacing_m")
    if "demihull_spacing_m" in config:
        if extension != ".csv":
            raise ValueError("demihull_spacing_m is only supported for CSV demi-hull offsets")
        if (isinstance(spacing, bool) or not isinstance(spacing, (int, float))
                or not np.isfinite(spacing) or spacing <= 0):
            raise ValueError("demihull_spacing_m must be a positive finite centerline-to-centerline distance")

    if extension in {".csv", ".mgf"}:
        from .offset_geometry import panels_from_sections, read_offset_sections, refine_large_panels

        if extension == ".mgf":
            if submerged:
                raise ValueError("MGF import currently supports surface vessels only")
            from .shipx_geometry import read_mgf_sections
            try:
                draft = float(config["draft_m"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError("MGF geometry requires positive draft_m") from exc
            sections = read_mgf_sections(path, draft)
        else:
            sections = read_offset_sections(path, submerged=submerged)
        stations = config.get("number_of_stations")
        vertices, faces = panels_from_sections(
            sections,
            samples_per_section=int(config.get("samples_per_section", 20)),
            number_of_stations=None if stations is None else int(stations),
        )
        if "max_panel_radius_m" in config:
            vertices, faces = refine_large_panels(vertices, faces, config["max_panel_radius_m"])
        source = "ShipX MGF hull sections" if extension == ".mgf" else "independent offset points"
        if spacing is not None:
            beam = 2 * max(float(profile[:, 0].max()) for _, profile in sections)
            if spacing <= beam:
                raise ValueError("demihull_spacing_m must exceed the demi-hull beam so the hulls do not touch or overlap")
            # Build one complete locally symmetric hull, then translate two
            # copies. Keep the components separate and preserve outward winding.
            count = len(vertices)
            vertices = np.vstack((vertices + [0, -spacing/2, 0],
                                  vertices + [0, spacing/2, 0]))
            faces = np.vstack((faces, faces + count))
            source = "CSV demi-hull offsets (two identical rigidly connected hulls)"
        return vertices, faces, source

    import capytaine as cpt

    if extension == ".stl":
        from .stl_geometry import read_stl
        from .output_coordinates import transform_positions
        vertices, faces = read_stl(path)
        mesh = cpt.Mesh(vertices=transform_positions(vertices), faces=faces, name="STL hull")
        if not submerged:
            mesh = mesh.immersed_part()
        if mesh.nb_faces == 0:
            raise ValueError("STL has no immersed hull panels; check FSD axes and waterline")
        return mesh.vertices.copy(), mesh.faces.copy(), "STL hull mesh"

    # Capytaine expands the GDF ISX/ISY reflections. Flatten the collection
    # so the existing analysis and panel exports receive the entire hull.
    mesh = cpt.load_mesh(str(path), file_format="gdf").merged()
    # A proper 180-degree rotation about z preserves outward panel normals.
    vertices = np.asarray(mesh.vertices, dtype=float) * [-1.0, -1.0, 1.0]
    faces = np.asarray(mesh.faces, dtype=int).copy()
    return vertices, faces, "WAMIT GDF mesh"
