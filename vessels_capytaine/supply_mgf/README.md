# Supply vessel — MGF geometry

The workflow is **MGF geometry → Capytaine → MSS MAT file**. An existing MGF file is sufficient; no ShipX installation, licence or computed result files are required.

`supply.mgf` was copied unchanged from the single `shipx_geometry.mgf` inside `MSS/HYDRO/vessels_shipx/supply/Supply.shipx_ship`. This path records geometry provenance only. The importer supports upright, surface-piercing monohull half sections, clips them at the configured draft and converts coordinates internally to the Capytaine solver frame. JSON inputs and exported results use MSS FSD.

The loading is Lpp 82.8 m, draft 6 m, mass 6,362,208.5 kg, FSD CG approximately `[-5.3859196, 0, -1.3000002]` m, and CG radii of gyration `[6.72, 20.7, 20.7]` m. JSON contains only geometry, loading, frequency, meshing and damping inputs. Other supported vessel properties are computed from the hull and inertia. No ShipX coefficient tables, roll-tank damping or drift forces are appended, and no ShipX comparison workflow is included.

```sh
python main.py supply_mgf
```

Supply uses 36 positive finite frequencies from 0.1047198 to 3.141593 rad/s. Capytaine separately computes zero- and infinite-frequency radiation; infinity is stored at the MSS marker `10`, which is excluded from RAOs and finite-frequency integration.

The mesh uses 24 contour samples and 70 longitudinal stations, followed by local triangle subdivision to a 0.75 m maximum panel radius. Subdivision preserves the hull surface and increases runtime and memory use. The panel-radius criterion does not establish BEM convergence; inspect resolution and convergence for your application.

The included `supply.mat` has been regenerated with all 36 configured finite frequencies through 3.141593 rad/s, plus separately computed zero- and infinite-frequency radiation. The hull has 13,160 panels. Saved matrices and RAOs have been checked for matching grids, expected dimensions and finite values; this is an output-integrity check, not an experimental validation or mesh-convergence study.

Preview source geometry and inspect the generated result in MATLAB/Octave:

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('supply')
vessel = plotVesselHydrodynamics('supply');
```

Results are written to `supply.mat`, `results/hydrostatics.json` and the generated hull-panel files. The format suffix stays in the catalogue folder; the MAT filename and body name remain `supply`.

`T_126 = [6.4, 2.4, 3.3]` s, `delta_zeta_345 = [0, 0.1, 0]` and `nonlinear_456 = [5, 0, 0]` are example damping settings stored in JSON and the MAT file. At the default JONSWAP reference (peak 1.0 rad/s, gamma 3.3), T_max is approximately [23.2448, 2.49945, 4.51327] s. The sway target was reduced from 2.5 to 2.4 s to keep additional viscous damping nonnegative. Damping-only edits can be applied with `python main.py supply_mgf --update-damping` without repeating BEM. See the [JSON guide](../../How%20to%20configure%20the%20vessel%20JSON%20file.md).

The included saved MAT file contains zero-speed first-order hydrodynamics and diagonal mean-drift coefficients under `vessel.drift`. The drift coefficients are signed surge, sway and yaw loads per squared incident wave amplitude. Cross-frequency QTFs are not computed.
