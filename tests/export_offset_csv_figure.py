"""Illustrate actual testShip CSV rows; run from any working directory."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


def main():
    root = Path(__file__).resolve().parents[1]
    source = root / "vessels_capytaine/testShip_csv/offset_points.csv"
    data = np.genfromtxt(source, delimiter=",", names=True)
    section = data[data["x_m"] == 0]
    waterline = np.sort(data[data["z_m"] == 0], order="x_m")
    blue, orange, gray = "#1763a6", "#c45d0a", "#909aa5"
    plt.rcParams.update({"font.size": 12, "axes.spines.top": False,
                         "axes.spines.right": False, "text.color": "#243344"})
    fig = plt.figure(figsize=(15, 8.6), facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 1.7], height_ratios=[1, 1],
                            left=0.045, right=0.97, bottom=0.14, top=0.86,
                            wspace=0.24, hspace=0.57)
    fig.suptitle("From CSV numbers to a hull section", x=0.045, y=0.97,
                 ha="left", fontsize=23, weight="bold")
    fig.text(0.045, 0.915, "Actual testShip offsets  |  Metres  |  MSS FSD: x forward, y starboard, z down",
             fontsize=13)

    rows = fig.add_subplot(grid[:, 0])
    rows.set_axis_off()
    rows.set_title("A. Five rows at x = 0 m", loc="left", pad=16, weight="bold", fontsize=15)
    rows.add_patch(FancyBboxPatch((0.12, 0.34), 0.86, 0.53,
                                 boxstyle="round,pad=0.018", facecolor="#f3f6fa",
                                 edgecolor="#dde3eb", transform=rows.transAxes))
    rows.text(0.17, 0.81, "x_m,z_m,half_breadth_m", family="monospace", fontsize=12,
              transform=rows.transAxes)
    for index, point in enumerate(section):
        height = 0.72 - 0.08 * index
        color = orange if index == 2 else blue
        rows.text(0.055, height, str(index + 1), color="white", ha="center", va="center",
                  weight="bold", fontsize=11, bbox=dict(boxstyle="circle,pad=0.33", fc=color, ec="none"),
                  transform=rows.transAxes)
        line = ",".join(f"{point[name]:g}" for name in data.dtype.names)
        rows.text(0.17, height, line, family="monospace", fontsize=16,
                  color=color, va="center", weight="bold" if index == 2 else "normal",
                  transform=rows.transAxes)
    rows.text(0.03, 0.26, "Numbers 1-5 identify the plotted points.\nDo not include these labels in the CSV.",
              fontsize=11, linespacing=1.6, transform=rows.transAxes)
    rows.text(0.03, 0.10, "Row 3:  0,1,3\nx = 0 m, z = 1 m, half breadth = 3 m\nPlotted location: (x, y, z) = (0, 3, 1) m",
              fontsize=12, color=orange, linespacing=1.8, transform=rows.transAxes)

    plan = fig.add_subplot(grid[0, 1])
    plan.set_title("B. x selects a station along the hull", loc="left", pad=15, weight="bold", fontsize=15)
    x, breadth = waterline["x_m"], waterline["half_breadth_m"]
    plan.fill_between(x, -breadth, breadth, color="#edf3f9")
    plan.plot(x, breadth, color=gray, lw=1.5)
    plan.plot(x, -breadth, color=gray, lw=1.5)
    plan.axhline(0, color=gray, ls="--", lw=0.9)
    plan.plot([0, 0], [-4, 4], color=orange, lw=3)
    plan.text(0.5, -4.3, "Selected station: x = 0", color=orange, fontsize=11)
    plan.set(xlim=(-11, 11), ylim=(5, -5.6), xlabel="x forward [m]", ylabel="y starboard [m]")
    plan.set_xticks([-10, -5, 0, 5, 10])
    plan.set_yticks([-4, 0, 4])
    plan.set_aspect("equal", adjustable="box")
    plan.text(-10.3, -4.5, "Waterline outline (top view)", fontsize=10, color=gray)

    profile = fig.add_subplot(grid[1, 1])
    profile.set_title("C. Half breadth becomes y; z is positive down", loc="left", pad=28,
                      weight="bold", fontsize=15)
    half, depth = section["half_breadth_m"], section["z_m"]
    profile.plot(-half, depth, "--", color=gray, lw=2)
    profile.plot(half, depth, color=blue, lw=2)
    profile.axhline(0, color="#738697", lw=1)
    profile.axvline(0, color=gray, ls=":", lw=1)
    profile.text(-4.55, -0.12, "Waterline: z = 0", fontsize=10)
    profile.text(-3.7, 1.5, "Port half\nmirrored automatically", color=gray, fontsize=11,
                  ha="center", linespacing=1.5)
    for index, (y, z) in enumerate(zip(half, depth)):
        color = orange if index == 2 else blue
        profile.scatter([y], [z], s=220, color=color, zorder=5)
        profile.text(y, z, str(index + 1), ha="center", va="center", color="white",
                      fontsize=10, weight="bold", zorder=6)
    profile.annotate("", xy=(3, 1), xytext=(0, 1),
                      arrowprops=dict(arrowstyle="<->", color=orange, lw=1.5))
    profile.text(1.5, 0.90, "half breadth = 3 m", ha="center", color=orange, fontsize=10)
    profile.text(0.16, 2.18, "1: keel", fontsize=10, color=blue)
    profile.text(3.3, -0.28, "5: waterline", fontsize=10, color=blue)
    profile.set(xlim=(-4.8, 4.8), ylim=(2.4, -0.45), xlabel="y starboard [m]", ylabel="z down [m]")
    profile.set_xticks([-4, -2, 0, 2, 4])
    profile.set_yticks([0, 1, 2])
    profile.set_aspect("equal", adjustable="box")
    fig.text(0.045, 0.055, "Keep station rows together, ordered from keel to waterline (1 to 5). Repeat for other x stations.",
             fontsize=12)
    output = root / "docs/images/offset-points-csv.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, facecolor="white")
    plt.close(fig)
    print(f"Exported {output}")


if __name__ == "__main__":
    main()
