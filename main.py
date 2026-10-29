"""Run a zero-speed MSS-Capytaine case from CSV, MGF, GDF or STL geometry.

Examples (from the repository root):
    python main.py --list
    python main.py testShip_csv
    python main.py otter_80kg_csv
    python main.py supply_mgf
    python main.py tanker_gdf
    python main.py wigley_stl
    python main.py path/to/config.json
    python main.py supply_mgf --update-damping

CSV/MGF sections generate panels; GDF/STL provide existing panels. JSON
coordinates always use MSS FSD. See "How to configure the vessel JSON file.md".
"""

import argparse
import json
import logging
from pathlib import Path

from src.hydrodynamics import run
from src.damping import refresh_damping


ROOT = Path(__file__).resolve().parent
CATALOGUE = ROOT / "vessels_capytaine"


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "case",
        nargs="?",
        default="testShip",
        help=(
            "catalogue folder, unique vessel name, or path to a JSON configuration "
            "(default: testShip_csv)"
        ),
    )
    parser.add_argument(
        "--update-damping",
        action="store_true",
        help="update damping settings in an existing MAT file from JSON without a BEM solve",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="list available vessel-format folders and exit",
    )
    return parser.parse_args()


def _config_path(case: str) -> Path:
    path = Path(case)
    if path.suffix.lower() == ".json":
        return path
    direct = CATALOGUE / path / "config.json"
    if direct.is_file():
        return direct
    matches = [p for p in CATALOGUE.glob("*/config.json")
               if json.loads(p.read_text()).get("body_name") == case] if len(path.parts) == 1 else []
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValueError(f"Ambiguous vessel {case!r}; use its catalogue folder from --list")
    raise FileNotFoundError(f"Unknown catalogue case {case!r}; use --list or a JSON path")


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")
    arguments = _arguments()
    if arguments.list:
        cases = sorted(str(path.parent.relative_to(CATALOGUE)) for path in CATALOGUE.glob("*/config.json"))
        print("\n".join(cases))
    elif arguments.update_damping:
        print(f"Updated damping settings in {refresh_damping(_config_path(arguments.case))}")
    else:
        print(f"Wrote vessel to {run(_config_path(arguments.case))}")
