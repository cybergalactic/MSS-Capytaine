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

The `viscous_damping` object specifies damping added to the potential-flow
radiation damping `B(ω)`. Its six nonnegative entries are viscous time
constants in seconds for surge, sway, and yaw, and additional dimensionless
damping ratios for heave, roll, and pitch. Surge, sway, and yaw have no
hydrostatic restoring; their time constants specify the intended diagonal
viscous damping, including yaw. A zero in any entry disables added damping
for that DOF. The test ship retains its surge/sway/yaw time constants while
defaulting to no added heave/roll/pitch damping:

```json
"viscous_damping": {
  "surge_viscous_time_constant_s": 10,
  "sway_viscous_time_constant_s": 50,
  "heave_additional_damping_ratio": 0,
  "roll_additional_damping_ratio": 0,
  "pitch_additional_damping_ratio": 0,
  "yaw_viscous_time_constant_s": 20
}
```

The exported `Bv` is diagonal at the center of gravity and constant across
all coefficient frequencies. For a nonzero time constant `Ti` in surge,
sway, or yaw, `Bvii = (MRBii + Aii(0)) / Ti`. For heave, roll, and pitch,
`Bvii = 2 ζv,i ωn,i (MRBii + Aii(ωn,i))`, where `ζv,i` is the
**additional** ratio and `ωn,i` is the estimated undamped natural frequency
from `MRB + A(ω)` and `C`. Nothing is subtracted from `B(ω)`.
A smaller nonzero time constant means more viscous damping. If a required
natural frequency lies beyond the finite calculation grid, extend
`omega_rad_s` rather than extrapolating.

The motion RAOs include `B(ω) + Bv`; force RAOs are unchanged. The total
damping remains frequency-dependent because `B(ω)` varies. Set all six
entries to zero, or `viscous_damping` to `null`, to export zero `Bv`.
The old `total_damping` settings are rejected: their values must be
reconsidered before using the additive model.

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
