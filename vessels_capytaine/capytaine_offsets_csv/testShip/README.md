# testShip

`testShip` is a synthetic surface monohull used to demonstrate and test the complete MSS-Capytaine workflow. It is not based on a particular vessel and should not be interpreted as a validated ship design. The case provides a simple, reproducible hull with smooth geometry, approximate mass properties, and enough asymmetry between the principal dimensions to exercise all parts of the Capytaine-to-MSS conversion.

## Hull geometry

The hull is a symmetric, double-ended monohull with the following principal dimensions:

| Quantity | Value |
| --- | ---: |
| Length | 20 m |
| Maximum beam | 8 m |
| Draft | 2 m |

The input geometry uses a body-fixed coordinate system whose origin is at midships on the design waterline. In the offset table, `x` points forward, `z` points down, and `half_breadth_m` is the nonnegative distance from the centerplane to one side of the hull.

The offsets are samples of the analytic half-breadth distribution

```text
half_breadth(x,z) = 4 [1 - (x/10)^2] [1 - (z/2)^2]
```

for `-10 <= x <= 10` m and `0 <= z <= 2` m. The two parabolic factors produce a maximum half-breadth of 4 m at midships on the waterline. The breadth decreases smoothly toward the keel, bow, and stern, reaching zero at the keel and at both ends. The resulting hull is symmetric about the centerplane and about midships.

[`offset_points.csv`](offset_points.csv) tabulates this equation at 11 longitudinal stations and five vertical levels. During mesh generation, each half-section is resampled uniformly in section arc length using 24 points. The sections are then interpolated to 60 longitudinal stations, mirrored to port and starboard, joined with triangular panels, and closed at the bow and stern. The current settings generate 2,774 vertices and 5,428 triangular hull panels.

The wetted hull remains open at the design waterplane in the exported panel file. For the surface-vessel boundary-element calculation, Capytaine generates an internal waterplane lid to suppress irregular-frequency artifacts. This lid is a numerical solver aid and is not included in `generated_hull_panels.mat`.

## Mass and hydrostatic assumptions

The case uses the synthetic properties in [`config.json`](config.json):

| Quantity | Value |
| --- | ---: |
| Mass | 141,954.8 kg |
| Center of gravity | `[0, 0, 1]` m in MSS FSD axes |
| Roll radius of gyration, `R44` | 2.4194 m |
| Pitch radius of gyration, `R55` | 5.8156 m |
| Yaw radius of gyration, `R66` | 6.2470 m |
| Water density | 1,025 kg/m^3 |
| Gravity | 9.81 m/s^2 |

The center of gravity is at midships, on the centerplane, and 1 m below the waterline. MSS FSD is always used; input CSV, JSON, exported panels, and hydrodynamic outputs all use FSD. The exported CG is `[0, 0, 1]` m. All exported rigid-body, hydrodynamic, and restoring matrices are referenced to this center of gravity.

Capytaine calculates the displaced volume, center of buoyancy, waterplane area, and hydrostatic restoring matrix from the generated mesh. The mass and radii of gyration are modeling inputs rather than measurements from a real ship. Their purpose is to produce a representative, approximately balanced demonstration case.

The exporter also stores the center of flotation as `vessel.main.CF` in MSS FSD coordinates. It constructs the surface-vessel restoring matrix using `x_F = CF(1) - CG(1)`, including the heave-pitch coupling and its associated pitch-stiffness shift. For this fore-aft-symmetric test hull, `CF(1) = CG(1)` and `x_F = 0`.

## Frequencies, headings, and damping inputs

The configured finite-frequency grid is nonuniform and extends from 0.01 to 4.5 rad/s. Zero-frequency and infinite-frequency radiation problems are solved separately. The infinite-frequency result is stored at the conventional 10 rad/s plotting location and is not treated as a finite-frequency solution.

Capytaine solves 19 wave propagation directions from 0 degrees to 180 degrees in 10-degree increments. The exporter mirrors these results to the full MSS 36-heading grid from 0 degrees to 350 degrees.

As a surface vessel, `testShip` uses total velocity-decay targets for surge, sway and yaw, and damping-ratio increments for restored heave, roll and pitch:

```json
"viscous_damping": {
  "T_126": [20, 8, 20],
  "delta_zeta_345": [0, 0.1, 0]
}
```

`T_126` is in seconds, ordered surge, sway, yaw. These are synthetic example estimates, checked below the computed potential-flow maxima at the 0.8 rad/s processing default and at 1.0 rad/s. MATLAB reports T_max and rejects infeasible targets rather than adding negative viscous damping. `delta_zeta_345` adds 10% of critical damping in roll, with no additional heave or pitch damping. See the [equations and migration instructions](../../../How%20to%20configure%20the%20vessel%20JSON%20file.md#viscous-damping-equations-and-examples).

Python stores these inputs in `vessel.powerBased`. The local MATLAB processing helper applies them while the native MSS surface-interface update is postponed. After damping-only JSON edits, run `python main.py testShip --update-damping`; no BEM rerun is needed. Submerged cases use T_1236 and delta_zeta_45.

## Generate and inspect the case

Run the hydrodynamic calculation from the repository root:

```sh
python main.py testShip
```

The generated MSS structure and panel geometry are written to:

```text
vessels_capytaine/capytaine_offsets_csv/testShip/results/testShip.mat
vessels_capytaine/capytaine_offsets_csv/testShip/results/generated_hull_panels.mat
```

Plot the generated coefficients and response amplitude operators with Python:

```sh
python plot_results.py vessels_capytaine/capytaine_offsets_csv/testShip/results/testShip.mat
python plot_results.py vessels_capytaine/capytaine_offsets_csv/testShip/results/testShip.mat --heading 90
```

With MSS on the MATLAB or GNU Octave path, load the generated vessel and use the vessel-independent functions to inspect its mesh and hydrodynamic data:

```matlab
addpath(genpath('/path/to/MSS'))
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('testShip')
plotVesselHydrodynamics('testShip');
```

The result files contain the rigid-body mass matrix, hydrostatic restoring, zero-, finite-, and infinite-frequency added mass and radiation damping, force and motion response amplitude operators, frequency and heading grids, and the inputs for the MSS power-based damping model.
