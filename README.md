# Capytaine test ship

This project computes six-degree-of-freedom hydrodynamics for a synthetic hull
defined by offset points. The result is an MSS-style MATLAB vessel structure.

## Run

From the project directory, run:

```sh
python main.py
```

`main.py` reads `capytaineTestShip/config.json` and writes results to
`capytaineTestShip/results/`. NumPy, SciPy, Xarray, and Capytaine are required.
To load and plot the generated vessel in MATLAB, run
`matlab/dataCapytaine.m` with MSS on the MATLAB path.

The solver and export modules live in `src/`.

## Plot in Python

After running `main.py`, use Matplotlib to plot the saved MATLAB result:

```sh
python plot_results.py
python plot_results.py --heading 90 --show
```

The script reads `capytaineTestShip/results/capytaineTestShip.mat` and saves
six PNG figures in `capytaineTestShip/results/plots/`. It does not rerun the
hydrodynamics or change the result file. The coefficient plots show the
infinite-frequency values as separate markers at the `10 rad/s` label;
RAO plots use only finite frequencies.

## Inputs

`capytaineTestShip/offset_points.csv` contains the columns
`x_m,z_m,half_breadth_m`. The origin is at midships on the design waterline;
`x` points aft, `z` points upward, and half breadth is nonnegative. Each
section runs from keel to waterline (`z = 0`).

`capytaineTestShip/config.json` specifies the mesh resolution, mass, radii of
gyration (or an inertia-matrix CSV), center of mass, wave frequencies, and
headings. Set `rotation_center_m` equal to `center_of_mass_m` for the
CG-referenced export. The limiting-frequency calculations require infinite
water depth, which is the default when `water_depth_m` is absent or `null`.

`samples_per_section` controls resolution around each half section, while
`number_of_stations` controls resolution along the hull. Check Capytaine's
mesh-resolution warnings at the highest wave frequencies.

## Viscous damping

The `total_damping` object specifies reference-frequency targets for the combined damping
`D(ω) = B(ω) + Bv`, with one named value per DOF: time constants in seconds
for surge, sway, and yaw, and dimensionless damping ratios for heave, roll,
and pitch. The test-ship values are:

```json
"total_damping": {
  "surge_time_constant_s": 10,
  "sway_time_constant_s": 50,
  "heave_damping_ratio": 0.2,
  "roll_damping_ratio": 0.2,
  "pitch_damping_ratio": 0.3,
  "yaw_time_constant_s": 20
}
```

The exported `Bv` is diagonal at the center of gravity and constant across
all coefficient frequencies. Surge, sway, and yaw set the target total
damping to `(MRBii + Aii(0)) / Ti`, assuming `Bii(0) = 0`. For heave, roll,
and pitch, the target is `2 ζi ωn (MRBii + Aii(ωn))` at an estimated natural
frequency from `MRB + A(ω)` and `C`. The program subtracts interpolated
`Bii(ωn)` to obtain `Bvii`. If potential damping already exceeds the target,
that `Bv` entry is set to zero with a warning. If a natural frequency lies beyond
the finite calculation grid, extend `omega_rad_s` rather than extrapolating.

The motion RAOs include `B(ω) + Bv`; force RAOs are unchanged. The total
damping remains frequency-dependent because `B(ω)` varies. Set
`total_damping` to `null` to export zero `Bv`.

## Outputs

`capytaineTestShip/results/capytaineTestShip.mat` contains `MRB`, `A`, `B`,
`Bv`, `C`, force and motion RAOs, frequencies, and headings in MSS
forward-starboard-down axes at the center of gravity. The hydrostatic
restoring matrix `C` has entries only in the heave, roll, and pitch block.

The coefficient frequency grid includes zero frequency and the
infinite-frequency radiation solution, labeled `10 rad/s`. The force and
motion RAOs use only the positive finite frequencies below `10 rad/s`.
Headings solved from 0° to 180° are mirrored to a full 0° to 350° set.

The results directory also contains `hydrostatics.json` and the generated
hull panels in `.mat` and `.npz` formats. The hull and inertia values are
synthetic demonstration data; numerical accuracy depends on mesh resolution.
