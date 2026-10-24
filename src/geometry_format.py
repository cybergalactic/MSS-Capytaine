"""Recognize hull geometry from its extension or unambiguous file structure."""
from pathlib import Path
import struct
import re


def geometry_format(path: Path) -> str:
    extension = path.suffix.lower()
    if extension in {".csv", ".mgf", ".gdf", ".stl"}:
        return extension
    if not path.is_file():
        raise ValueError(f"Cannot identify geometry format {extension!r}: {path}")
    size = path.stat().st_size
    with path.open('rb') as stream:
        prefix = stream.read(4096)
    if len(prefix) >= 84:
        count = struct.unpack_from('<I', prefix, 80)[0]
        if count > 0 and size == 84+50*count:
            return '.stl'
    text = prefix.decode('utf-8-sig', errors='replace')
    lines = text.splitlines()
    if lines and {'x_m', 'z_m', 'half_breadth_m'}.issubset(
            {x.strip() for x in lines[0].split(',')}):
        return '.csv'
    lowered = text.lower()
    if re.search(r'\bfacet\s+normal\b', lowered) and 'vertex' in lowered and 'endfacet' in lowered:
        return '.stl'
    try:
        if len(lines) >= 4:
            scales = list(map(float, lines[1].split()[:2]))
            flags = list(map(int, lines[2].split()[:2]))
            count = int(lines[3].split()[0])
            if len(scales) == 2 and min(scales) > 0 and len(flags) == 2 and set(flags) <= {0, 1} and count > 0:
                return '.gdf'
    except (ValueError, IndexError):
        pass
    try:
        if (len(lines) >= 9 and float(lines[4]) > 0 and int(lines[5]) >= 0
                and int(lines[7]) >= 2 and len(list(map(float, lines[8].split()))) == 2):
            float(lines[6])
            return '.mgf'
    except (ValueError, IndexError):
        pass
    raise ValueError(f"Cannot identify hull geometry in {path}; use CSV, MGF, GDF or STL")
