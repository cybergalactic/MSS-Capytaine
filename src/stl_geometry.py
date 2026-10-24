"""Small ASCII/binary STL reader without optional mesh-library dependencies."""

from pathlib import Path
import re
import struct
import numpy as np


def read_stl(path: Path):
    """Return indexed triangles; source coordinates must be metres in MSS FSD."""
    data = path.read_bytes()
    count = struct.unpack_from("<I", data, 80)[0] if len(data) >= 84 else 0
    if count and len(data) == 84 + 50 * count:
        record = np.dtype([("normal", "<f4", 3), ("vertices", "<f4", (3, 3)),
                           ("attribute", "<u2")])
        triangles = np.frombuffer(data, record, count=count, offset=84)["vertices"].astype(float)
    else:
        try:
            text = data.decode("ascii")
            facets = re.findall(r"\bfacet\s+normal\b(.*?)\bendfacet\b", text, re.S | re.I)
            triangles = []
            for facet in facets:
                vertices = re.findall(r"\bvertex\s+([^\r\n]+)", facet, re.I)
                if len(vertices) != 3:
                    raise ValueError("Each ASCII STL facet must have three vertices")
                triangles.append([[float(x) for x in vertex.split()] for vertex in vertices])
            triangles = np.asarray(triangles, dtype=float)
        except (ValueError, UnicodeDecodeError) as exc:
            raise ValueError(f"Invalid STL file: {path}") from exc
    if triangles.ndim != 3 or triangles.shape[1:] != (3, 3) or not np.isfinite(triangles).all():
        raise ValueError("STL must contain finite triangular facets")
    vertices, inverse = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    faces = inverse.reshape(-1, 3)
    return vertices, np.column_stack((faces, faces[:, -1]))
