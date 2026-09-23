"""Run the synthetic offset-point Capytaine test ship case."""

import logging
from pathlib import Path

from src.hydrodynamics import run


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")
    config = Path(__file__).resolve().parent / "capytaineTestShip" / "config.json"
    print(f"Wrote capytaineTestShip results to {run(config)}")
