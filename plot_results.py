"""Plot the MSS vessel data written by main.py.

Examples:
    python plot_results.py
    python plot_results.py capytaineTestShip/results/capytaineTestShip.mat --heading 90
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


DOFS = ("Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw")
DEFAULT_RESULT = (
    Path(__file__).resolve().parent
    / "capytaineTestShip"
    / "results"
    / "capytaineTestShip.mat"
)


def load_vessel(path: Path) -> object:
    """Read the single MSS vessel structure from a MATLAB result file."""
    contents = loadmat(path, struct_as_record=False, squeeze_me=False)
    if "vessel" not in contents or contents["vessel"].size != 1:
        raise ValueError("The MATLAB file must contain one vessel structure")
    return contents["vessel"].item()


def coefficient(vessel: object, field: str, frequencies: np.ndarray) -> np.ndarray:
    values = np.asarray(getattr(vessel, field), dtype=float)
    expected = (6, 6, len(frequencies), 1)
    if values.shape != expected:
        raise ValueError(f"vessel.{field} has shape {values.shape}; expected {expected}")
    return values[:, :, :, 0]


def rao_data(
    vessel: object, field: str, headings: np.ndarray, frequencies: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return finite frequencies and frequency x heading x DOF RAO arrays."""
    rao = getattr(vessel, field)
    # loadmat retains a one-element array around nested MATLAB structures.
    if isinstance(rao, np.ndarray):
        if rao.size != 1:
            raise ValueError(f"vessel.{field} must contain one RAO structure")
        rao = rao.item()
    omega = np.asarray(rao.w, dtype=float).ravel()
    if omega.shape != frequencies[1:-1].shape or not np.allclose(
        omega, frequencies[1:-1]
    ):
        raise ValueError(f"vessel.{field}.w does not match finite frequencies")

    def unpack(cells: np.ndarray, name: str) -> np.ndarray:
        entries = np.asarray(cells, dtype=object).ravel()
        if len(entries) != 6:
            raise ValueError(f"vessel.{field}.{name} must contain six DOF cells")
        expected = (len(omega), len(headings), 1)
        arrays = [np.asarray(entry, dtype=float) for entry in entries]
        if any(array.shape != expected for array in arrays):
            raise ValueError(
                f"vessel.{field}.{name} cells must have shape {expected}"
            )
        return np.stack([array[:, :, 0] for array in arrays], axis=2)

    return omega, unpack(rao.amp, "amp"), unpack(rao.phase, "phase")


def plot_diagonal(
    frequencies: np.ndarray,
    matrix: np.ndarray,
    field: str,
    path: Path,
    show: bool,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for index, (dof, ax) in enumerate(zip(DOFS, axes.flat)):
        ax.plot(frequencies[:-1], matrix[index, index, :-1], marker="o", ms=3)
        # The 10 rad/s entry is the infinite-frequency solution, not A(10) or B(10).
        ax.plot(
            frequencies[-1], matrix[index, index, -1], marker="x",
            linestyle="none", color="tab:red", label="∞ (shown at 10 rad/s)",
        )
        ax.set(title=dof, xlabel="Angular frequency ω (rad/s)")
        if field == "A":
            ax.set_ylabel("kg" if index < 3 else "kg m²")
        else:
            ax.set_ylabel("kg/s" if index < 3 else "kg m²/s")
        ax.grid(True, alpha=0.3)
    axes.flat[0].legend(fontsize="small")
    title = "Diagonal added mass" if field == "A" else "Diagonal potential damping"
    fig.suptitle(f"Capytaine test ship: {title}")
    fig.savefig(path, dpi=180)
    if not show:
        plt.close(fig)


def plot_raos(
    omega: np.ndarray,
    amplitude: np.ndarray,
    phase: np.ndarray,
    heading_index: int,
    heading_deg: float,
    field: str,
    output_dir: Path,
    show: bool,
) -> list[Path]:
    label = "Motion" if field == "motionRAO" else "Excitation force"
    stem = "motion_rao" if field == "motionRAO" else "force_rao"
    order = np.argsort(omega)
    outputs = []
    for part in ("magnitude", "phase"):
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
        for index, (dof, ax) in enumerate(zip(DOFS, axes.flat)):
            magnitude = amplitude[:, heading_index, index][order]
            if part == "magnitude":
                plotted = magnitude
                if field == "motionRAO":
                    unit = "m/m" if index < 3 else "rad/m"
                else:
                    unit = "N/m" if index < 3 else "N m/m"
            else:
                plotted = np.rad2deg(
                    np.mod(phase[:, heading_index, index][order], 2 * np.pi)
                )
                cutoff = max(float(np.max(magnitude)) * 1e-8, np.finfo(float).tiny)
                plotted[magnitude <= cutoff] = np.nan
                unit = "degrees"
                ax.set_ylim(0, 360)
            ax.plot(omega[order], plotted, marker="o", ms=3)
            ax.set(title=dof, xlabel="Angular frequency ω (rad/s)", ylabel=unit)
            ax.grid(True, alpha=0.3)
        fig.suptitle(f"Capytaine test ship: {label} RAO {part}, {heading_deg:g}°")
        path = output_dir / f"{stem}_{part}_{heading_deg:g}deg.png"
        fig.savefig(path, dpi=180)
        if not show:
            plt.close(fig)
        outputs.append(path)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "result", nargs="?", type=Path, default=DEFAULT_RESULT,
        help="MSS vessel .mat file (default: capytaineTestShip result)",
    )
    parser.add_argument(
        "--heading", type=float, default=0.0,
        help="Wave heading in degrees (default: 0)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Directory for PNG files (default: plots/ beside the .mat file)",
    )
    parser.add_argument(
        "--show", action="store_true", help="Keep plot windows open",
    )
    args = parser.parse_args()

    result_path = args.result.expanduser().resolve()
    if not result_path.is_file():
        parser.error(f"Result file not found: {result_path}; run main.py first")
    if not np.isfinite(args.heading):
        parser.error("--heading must be a finite number")

    vessel = load_vessel(result_path)
    frequencies = np.asarray(vessel.freqs, dtype=float).ravel()
    headings = np.asarray(vessel.headings, dtype=float).ravel()
    if len(frequencies) < 3 or not np.isclose(frequencies[-1], 10.0):
        parser.error("Expected zero, finite, and 10 rad/s coefficient frequencies")
    if not len(headings):
        parser.error("The vessel has no wave headings")

    requested = np.deg2rad(args.heading % 360.0)
    matches = np.flatnonzero(np.isclose(headings, requested, rtol=0, atol=1e-8))
    if len(matches) != 1:
        choices = ", ".join(f"{value:g}" for value in np.rad2deg(headings))
        parser.error(f"Heading {args.heading:g}° is unavailable; choose: {choices}°")
    heading_index = int(matches[0])
    heading_deg = float(np.rad2deg(headings[heading_index]))

    output_dir = args.output_dir or result_path.parent / "plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for field, filename in (("A", "added_mass.png"), ("B", "radiation_damping.png")):
        path = output_dir / filename
        plot_diagonal(frequencies, coefficient(vessel, field, frequencies),
                      field, path, args.show)
        outputs.append(path)
    for field in ("motionRAO", "forceRAO"):
        omega, amplitude, phase = rao_data(vessel, field, headings, frequencies)
        outputs.extend(plot_raos(omega, amplitude, phase, heading_index,
                                 heading_deg, field, output_dir, args.show))

    print("Saved plots:")
    for path in outputs:
        print(path)
    if args.show:
        plt.show(block=True)


if __name__ == "__main__":
    main()
