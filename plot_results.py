"""Plot saved MSS FSD hydrodynamic results without rerunning the solver.

Examples:
    python plot_results.py
    python plot_results.py vessels_capytaine/testShip_csv/testShip.mat --heading 90
    python plot_results.py vessels_capytaine/wigley_stl/wigley.mat

Works with results from all four geometry formats. Saves six PNG files, plus
a signed mean-drift plot when vessel.drift is available, in
results/plots/ beside the result unless --output-dir is supplied. The JSON
output_dir customizes the supporting results folder. The 10 rad/s
coefficient marker represents infinite frequency; RAOs use finite frequencies.
MATLAB/Octave plotVesselGeometry shows source geometry and saved solver hull
panels together in two vertically stacked subplots; plotVesselHydrodynamics displays MSS analyses.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat


DOFS = ("Surge", "Sway", "Heave", "Roll", "Pitch", "Yaw")
DEFAULT_RESULT = (
    Path(__file__).resolve().parent
    / "vessels_capytaine"
    / "testShip_csv"
    / "testShip.mat"
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
    import matplotlib.pyplot as plt
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
    fig.suptitle(title)
    fig.savefig(path, dpi=180)
    if not show:
        plt.close(fig)


def drift_data(vessel, headings, frequencies):
    """Read Chapter 10 cells without squeezing singleton frequency/speed axes."""
    drift = vessel.drift.item()
    omega = np.asarray(drift.w, dtype=float).ravel()
    if not np.array_equal(omega, frequencies[1:-1]):
        raise ValueError("vessel.drift.w must match the finite RAO frequency grid")
    entries = np.asarray(drift.coeff, dtype=object).ravel()
    shape = (len(omega), len(headings), 1)
    if len(entries) != 3 or any(np.asarray(entry).shape != shape for entry in entries):
        raise ValueError(f"vessel.drift.coeff must contain three cells of shape {shape}")
    if any(np.iscomplexobj(entry) for entry in entries):
        raise ValueError("Mean-drift coefficients must be real and signed")
    values = np.stack([np.asarray(entry, dtype=float)[:, :, 0] for entry in entries], axis=2)
    if not np.all(np.isfinite(values)):
        raise ValueError("Mean-drift coefficients must be finite")
    return omega, values


def plot_drift(omega, values, headings, path, show):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)
    colors = plt.get_cmap('turbo')(np.linspace(.02, .98, len(headings)))
    degrees = np.rad2deg(headings)
    order = np.argsort(omega)
    for j, (label, ax) in enumerate(zip(('Surge', 'Sway', 'Yaw'), axes)):
        ax.axhline(0., color='0.55', linewidth=.7)
        for k, heading in enumerate(degrees):
            ax.plot(omega[order], values[order, k, j], color=colors[k],
                    linestyle='-' if heading < 180 else '--', linewidth=1.1,
                    label=f'{heading:g}°')
        ax.set(title=label, ylabel='N/m²' if j < 2 else 'N m/m²')
        ax.set_xlim(float(np.min(omega)), float(np.max(omega)))
        ax.grid(True, alpha=.3)
    axes[-1].set_xlabel('Angular frequency ω (rad/s)')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .01),
               ncol=min(9, len(headings)), fontsize=9, title='Propagation heading β')
    fig.subplots_adjust(left=.1, right=.97, bottom=.17, top=.92, hspace=.28)
    fig.suptitle('Signed mean-drift coefficients per squared wave amplitude, about CG')
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
    import matplotlib.pyplot as plt
    def phase_curve_with_gaps(
        x: np.ndarray,
        y: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Insert gaps at genuine 180-degree phase reversals."""
        x_plot: list[float] = []
        y_plot: list[float] = []
        for i, (x_value, y_value) in enumerate(zip(x, y)):
            if (
                i > 0
                and np.isfinite(y[i - 1])
                and np.isfinite(y_value)
                and abs(y_value - y[i - 1]) >= 180.0 - 1e-8
            ):
                x_plot.append(np.nan)
                y_plot.append(np.nan)
            x_plot.append(float(x_value))
            y_plot.append(float(y_value))
        return np.asarray(x_plot), np.asarray(y_plot)

    label = "Motion" if field == "motionRAO" else "Excitation force"
    stem = "motion_rao" if field == "motionRAO" else "force_rao"
    order = np.argsort(omega)
    omega_sorted = omega[order]
    omega_limits = (float(omega_sorted[0]), float(omega_sorted[-1]))
    outputs = []
    for part in ("magnitude", "phase"):
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
        for index, (dof, ax) in enumerate(zip(DOFS, axes.flat)):
            magnitude = amplitude[:, heading_index, index][order]
            zero_response = not np.any(
                magnitude
                > max(
                    float(np.max(magnitude)) * 1e-8,
                    np.finfo(float).tiny,
                )
            )
            if part == "magnitude":
                plotted = magnitude
                if field == "motionRAO":
                    unit = "m/m" if index < 3 else "rad/m"
                else:
                    unit = "N/m" if index < 3 else "N m/m"
            else:
                phase_degrees = np.rad2deg(
                    phase[:, heading_index, index][order]
                )
                # Principal phase eliminates artificial 0/360-degree jumps.
                plotted = (phase_degrees + 180.0) % 360.0 - 180.0
                cutoff = max(float(np.max(magnitude)) * 1e-6, np.finfo(float).tiny)
                plotted[magnitude <= cutoff] = np.nan
                unit = "degrees (principal value)"
                ax.set_ylim(-180, 180)
                ax.set_yticks([-180, -90, 0, 90, 180])
            if part == "phase":
                x_plot, y_plot = phase_curve_with_gaps(omega_sorted, plotted)
                ax.plot(x_plot, y_plot, marker="o", ms=3)
            else:
                ax.plot(omega_sorted, plotted, marker="o", ms=3)
            ax.set_xlim(*omega_limits)
            if zero_response:
                message = (
                    "Zero response\n(phase undefined)"
                    if part == "phase"
                    else "Zero response"
                )
                ax.text(
                    0.5,
                    0.5,
                    message,
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    color="0.35",
                )
            ax.set(title=dof, xlabel="Angular frequency ω (rad/s)", ylabel=unit)
            ax.grid(True, alpha=0.3)
        part_label = "principal phase" if part == "phase" else part
        fig.suptitle(f"{label} RAO {part_label}, {heading_deg:g}°")
        path = output_dir / f"{stem}_{part}_{heading_deg:g}deg.png"
        fig.savefig(path, dpi=180)
        if not show:
            plt.close(fig)
        outputs.append(path)
    return outputs


def main() -> None:
    import matplotlib.pyplot as plt
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "result", nargs="?", type=Path, default=DEFAULT_RESULT,
        help="MSS vessel .mat file (default: testShip_csv result)",
    )
    parser.add_argument(
        "--heading", type=float, default=0.0,
        help="MSS wave propagation heading in degrees on the 10-degree grid (default: 0)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Directory for PNG files (default: JSON output_dir/plots, or results/plots)",
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

    config_path = result_path.parent / 'config.json'
    config = json.loads(config_path.read_text()) if config_path.is_file() else {}
    output_dir = args.output_dir or result_path.parent / config.get('output_dir', 'results') / 'plots'
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
    if hasattr(vessel, 'drift'):
        omega, values = drift_data(vessel, headings, frequencies)
        path = output_dir / 'mean_drift.png'
        plot_drift(omega, values, headings, path, args.show)
        outputs.append(path)

    print("Saved plots:")
    for path in outputs:
        print(path)
    if args.show:
        plt.show(block=True)


if __name__ == "__main__":
    main()
