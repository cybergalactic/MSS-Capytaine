"""Reproduce both Otter loading cases from the same assumed demi-hull shape.

Run from the repository root with:
    python tests/otter/generate_offsets.py

This is an idealized hydrostatically fitted hull, not measured Otter lines.
The paper and manufacturer photographs establish the general hull type only.
"""

from pathlib import Path
import json
import re
import sys

import capytaine as cpt
import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.offset_geometry import panels_from_sections

BASE = ROOT / "vessels_capytaine/otter_80kg_csv"
RHO = 1025.0
MASS = 80.0
DRAFT = MASS / RHO / (2 * 0.4 * 0.25 * 2.0)
CG_BODY = np.array([0.1375, 0.0, -0.1375])
SAMPLES, STATIONS = 12, 48


def sections(parameters):
    """Assumed chined sections in FSD, keel first, on one local centerline."""
    bow_taper, chine_breadth, stern_depth = parameters
    # The long aft taper and rocker are assumptions fitted below, not
    # measurements extracted from the photographs or the cited paper.
    xs = np.unique(np.round(np.r_[np.linspace(-1, 1, 21), -0.1,
                                  1 - bow_taper, 0.2, 0.7], 12))
    heights = np.array([0, 0.15, 0.5, 0.8, 1])
    result = []
    for x in xs:
        width = np.interp(x, [-1, -0.1, 1-bow_taper, 1], [0.35, 1, 1, 0])
        depth = np.interp(x, [-1, 0.2, 0.7, 1], [stern_depth, 1, 1, 0.05])
        y = 0.125 * width * np.interp(heights, [0, 0.5, 1], [0, chine_breadth, 1])
        z = DRAFT * depth * (1-heights)
        result.append((x, np.column_stack((y, z))))
    return result


def section_mesh(source):
    internal = [(-x, profile * [1, -1]) for x, profile in source[::-1]]
    v, f = panels_from_sections(internal, SAMPLES, STATIONS)
    return cpt.Mesh(vertices=v, faces=f)


def residual(parameters):
    hull = section_mesh(sections(parameters))
    body = cpt.FloatingBody(mesh=hull)
    values = np.array([2*hull.volume, 2*hull.waterplane_area,
                       -body.center_of_buoyancy[0]])
    # Match displacement, MSS waterplane area and zero-trim longitudinal
    # balance for this load. These integral constraints do not identify
    # unique hull lines, nor reproduce all MSS hydrostatic coefficients.
    target = np.array([MASS/RHO, 2*0.75*2*0.25, CG_BODY[0]])
    return (values-target) / [MASS/RHO, target[1], 1.0]


def clipped_sections(source, waterline):
    """Cut the same physical hull at a lighter load's lower waterline.

    A tiny immersed end section avoids duplicate points at a zero-depth tip.
    The waterline is refitted using the actual meshed displaced volume.
    """
    xs = np.array([x for x, _ in source])
    profiles = np.array([p for _, p in source])
    keels = profiles[:, 0, 1]
    target = waterline+1e-5
    tips = []
    for i in range(len(xs)-1):
        if (keels[i]-target)*(keels[i+1]-target) < 0:
            tips.append(xs[i]+(xs[i+1]-xs[i])*(target-keels[i])/(keels[i+1]-keels[i]))
    stations = np.unique(np.r_[xs[keels >= target], tips])
    result = []
    for x in stations:
        p = np.array([[np.interp(x, xs, profiles[:, j, k]) for k in range(2)]
                      for j in range(profiles.shape[1])])
        wet = p[p[:, 1] > waterline]
        y = np.interp(waterline, p[::-1, 1], p[::-1, 0])
        wet = np.vstack((wet, [y, waterline]))
        wet[:, 1] -= waterline
        result.append((x, wet))
    return result


def pitch_rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def trimmed_sections(hull, rise, heave, pitch):
    """Slice the rigidly rotated original hull; retain the immersed half.

    The CSV remains a set of vertical sections in the level equilibrium FSD
    frame. Triangulation supplies intersections of the existing lofted mesh;
    the subsequent fit uses the final CSV mesh, including resampling error.
    """
    rotation = pitch_rotation(pitch)
    v = (hull.vertices * [-1, 1, -1] - [0, 0, rise]) @ rotation.T
    v += [0, 0, heave]
    triangles = v[np.vstack((hull.faces[:, :3], hull.faces[:, [0, 2, 3]]))]
    nonzero = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                     triangles[:, 2]-triangles[:, 0]), axis=1) > 1e-12
    triangles = triangles[nonzero]
    # Determine the wetted x extent, including free-surface intersections.
    wet = [v[v[:, 2] >= 0]]
    for j in range(3):
        a, b = triangles[:, j], triangles[:, (j+1) % 3]
        crossing = a[:, 2]*b[:, 2] < 0
        a, b = a[crossing], b[crossing]
        wet.append(a+(b-a)*(-a[:, 2]/(b[:, 2]-a[:, 2]))[:, None])
    wet = np.vstack(wet)
    # Avoid sliver sections at the wet tips: sub-millimetre profiles can
    # collapse differently after translation and Capytaine vertex merging.
    xs = np.linspace(wet[:, 0].min()+1e-3, wet[:, 0].max()-1e-3, STATIONS)
    result = []
    for x in xs:
        candidates = triangles[(triangles[:, :, 0].min(axis=1) <= x) &
                               (triangles[:, :, 0].max(axis=1) >= x)]
        points = []
        for triangle in candidates:
            intersections = []
            for j in range(3):
                a, b = triangle[j], triangle[(j+1) % 3]
                if (a[0]-x)*(b[0]-x) <= 0 and abs(b[0]-a[0]) > 1e-12:
                    intersections.append(a+(b-a)*((x-a[0])/(b[0]-a[0])))
            if len(intersections) < 2:
                continue
            a, b = intersections[0], intersections[-1]
            # Clip each intersection segment to starboard and below water.
            for dimension in (1, 2):
                if max(a[dimension], b[dimension]) < -1e-10:
                    break
                if a[dimension]*b[dimension] < 0:
                    crossing = a+(b-a)*(-a[dimension]/(b[dimension]-a[dimension]))
                    if a[dimension] < 0:
                        a = crossing
                    else:
                        b = crossing
            else:
                points.extend((a[1:], b[1:]))
        # Adjacent triangles calculate the same intersection independently.
        # Round before sorting so roundoff cannot reverse flat-edge ordering.
        p = np.unique(np.maximum(np.round(np.asarray(points), 10), 0), axis=0)
        p = p[np.lexsort((p[:, 0], -p[:, 1]))]
        p = p[np.r_[True, np.linalg.norm(np.diff(p, axis=0), axis=1) > 1e-9]]
        if p[0, 0] != 0 or p[-1, 1] != 0 or len(p) < 2:
            raise RuntimeError("Trimmed section is not a keel-to-waterline profile")
        result.append((x, p))
    return result


def write_case(base, source, payload, waterline_shift, pitch=0, heave=0):
    total = 55+payload
    rows = [[x, z, y] for x, profile in source for y, z in profile]
    np.savetxt(base / "offset_points.csv", rows, delimiter=",", fmt="%.12f",
               header="x_m,z_m,half_breadth_m", comments="")
    cg_body = np.array([0.2, 0, -0.2])*55/total
    rotation = pitch_rotation(pitch)
    cg = rotation @ (cg_body+[0, 0, waterline_shift])+[0, 0, heave]
    text = (base / "config.json").read_text()
    text = re.sub(r'("center_of_mass_m": )\[[^\]]+\]',
                  lambda m: m[1]+json.dumps(np.round(cg, 4).tolist()), text)

    # MSS radii already represent the complete twin-hull vessel.
    hull_cg = np.array([0.2, 0, -0.2])
    inertia = 55 * np.diag(np.array([0.4*1.08, 0.25*2, 0.25*2])**2)
    for mass, position in [(55, hull_cg), (payload, np.zeros(3))]:
        r = position-cg_body
        inertia += mass * (np.dot(r, r)*np.eye(3)-np.outer(r, r))
    matrix = np.zeros((6, 6))
    matrix[:3, :3] = total * np.eye(3)
    matrix[3:, 3:] = rotation @ inertia @ rotation.T
    matrix = np.round(matrix, 10)
    radii = np.sqrt(np.diag(matrix)[3:]/total)
    cross_radius = np.sign(matrix[3, 5])*np.sqrt(abs(matrix[3, 5])/total)
    for key, value in zip(("R44_m", "R55_m", "R66_m", "R46_m"),
                          [*radii, cross_radius]):
        text, count = re.subn(r'("'+key+r'":\s*)[^,\n]+',
                             lambda m: m[1]+str(round(float(value), 4)), text)
        if count != 1:
            raise ValueError(f"Expected one {key} in each Otter config")
    (base / "config.json").write_text(text)


def generate():
    fitted = least_squares(residual, [0.415, 0.85, 0.4],
                           bounds=([0.1, 0.1, 0.02], [0.9, 1, 0.99]),
                           gtol=1e-12, xtol=1e-12, ftol=1e-12)
    if not fitted.success or np.max(np.abs(residual(fitted.x))) > 1e-6:
        raise RuntimeError("Could not fit the assumed hydrostatic constraints")
    source = sections(fitted.x)
    rise = brentq(lambda h: 2*section_mesh(clipped_sections(source, h)).volume*RHO-55,
                  0, 0.08, xtol=1e-12)
    # The MSS body origin is the unloaded 55 kg design waterline. Define
    # that waterline on this fixed shape, rather than reshaping the hull
    # using a constant block coefficient independently for both loads.
    original = section_mesh(source)

    def trim_residual(state):
        heave, pitch = state
        light = section_mesh(trimmed_sections(original, rise, heave, pitch))
        cb = cpt.FloatingBody(mesh=light).center_of_buoyancy * [-1, 1, -1]
        cg = pitch_rotation(pitch) @ np.array([0.2, 0, -0.2])
        return [2*light.volume*RHO/55-1, cb[0]-cg[0]]

    equilibrium = least_squares(trim_residual, [0, -0.0045],
                                bounds=([-0.01, -0.05], [0.01, 0.05]),
                                diff_step=1e-3, gtol=1e-10, ftol=1e-10, xtol=1e-10)
    if not equilibrium.success or max(abs(np.asarray(trim_residual(equilibrium.x)))) > 1e-6:
        raise RuntimeError("Could not find the 55 kg heave/trim equilibrium")
    heave, pitch = equilibrium.x
    write_case(BASE, source, payload=25, waterline_shift=rise)
    write_case(BASE.with_name("otter_55kg"), trimmed_sections(original, rise, heave, pitch),
               payload=0, waterline_shift=0, pitch=pitch, heave=heave)
    print("Assumed shape parameters (bow taper, chine breadth, stern depth):", fitted.x)
    print(f"Untrimmed reference drafts: 55 kg = {DRAFT-rise:.9f} m; 80 kg = {DRAFT:.9f} m")
    print(f"Waterline separation on the fixed assumed hull: {rise:.9f} m")
    print(f"55 kg equilibrium: heave {heave:.9f} m, FSD pitch {np.rad2deg(pitch):.9f} deg")


if __name__ == "__main__":
    generate()
