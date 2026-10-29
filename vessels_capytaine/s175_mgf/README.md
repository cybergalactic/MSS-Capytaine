# S175 container ship — MGF geometry

The workflow is **MGF geometry → Capytaine → MSS MAT file**. An existing MGF file is sufficient; no ShipX installation, licence or computed result files are required.

`s175.mgf` was copied unchanged from the single `shipx_geometry.mgf` inside `MSS/HYDRO/vessels_shipx/s175/S175.shipx_ship`. This path records geometry provenance only. The importer supports upright, surface-piercing monohull half sections, clips them at the configured draft and converts coordinates internally to the Capytaine solver frame. JSON inputs and exported results use MSS FSD.

The loading is Lpp 175 m, draft 9.5 m, mass 24,609,620 kg, FSD CG approximately `[-2.5474801, 0, -0.0500002]` m, and CG radii of gyration approximately `[8.331, 42, 42]` m. JSON contains only geometry, loading, frequency, meshing and damping inputs. Other supported vessel properties are computed from the hull and inertia. No ShipX coefficient tables, roll-tank damping or drift forces are appended, and no ShipX comparison workflow is included.

```sh
python main.py s175_mgf
```

S175 uses 36 positive finite frequencies from 0.1047198 to 3.141593 rad/s. Capytaine separately computes zero- and infinite-frequency radiation; infinity is stored at the MSS marker `10`, which is excluded from RAOs and finite-frequency integration.

The mesh uses 22 contour samples and 180 longitudinal stations, followed by local triangle subdivision to a 0.75 m maximum panel radius. The more even panel spacing avoids excessive subdivision of long thin triangles. Native port-starboard reflection symmetry is enabled for the hull and internal lid. The full mirrored solver mesh is exported for plotting. The panel-radius criterion does not establish BEM convergence; inspect resolution and convergence for your application.

The full S175 solve is complete; `s175.mat` contains its first-order hydrodynamic data. Numerical consistency and MSS checks are recorded below.

Preview source geometry and inspect the generated result in MATLAB/Octave:

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('s175')
vessel = plotVesselHydrodynamics('s175');
```

Results are written to `s175.mat`, `results/hydrostatics.json` and the generated hull-panel files. The format suffix stays in the catalogue folder; the MAT filename and body name remain `s175`.

`T_126 = [11.9, 2, 1.8]` s, `delta_zeta_345 = [0, 0.2, 0]` and `nonlinear_456 = [5, 0, 0]` are example damping settings stored in JSON and the MAT file. At the default JONSWAP reference (peak 1.0 rad/s, gamma 3.3), T_max is approximately [55.1601, 2.03210, 2.05363] s. Damping-only edits can be applied with `python main.py s175_mgf --update-damping` without repeating BEM. See the [JSON guide](../../How%20to%20configure%20the%20vessel%20JSON%20file.md).

The included saved MAT file contains zero-speed first-order hydrodynamics and diagonal mean-drift coefficients under `vessel.drift`. The drift coefficients are signed surge, sway and yaw loads per squared incident wave amplitude. Cross-frequency QTFs are not computed.

## Memory and numerical checks

The optimized solver hull has 15,684 panels and the internal lid has 3,170 panels, compared with 38,662 hull panels in the initial configuration. Capytaine reports approximately 8 GB peak RAM for the symmetric hull. Accounting for the lid gives approximately 11.4 GB of dense matrix storage; allow additional headroom for temporary arrays and the operating system. The estimate is not a measured full-solve peak. 

The maximum hull-panel radius is 0.74983 m, below wavelength/8 at 3.141593 rad/s. Reflection preserves the generated hull volume to rounding. Small-hull regression solves check radiation at zero, finite and infinite frequency and diffraction against an expanded full-mesh solve. This checks the symmetry implementation; a mesh-convergence study has not been performed.

## Completed first-order checks (28 October 2026)

- JSON mass, CG, inertia, frequency grid and headings agree with the MAT file. All exported coefficient and RAO entries are finite; zero and infinity remain separate from the 36 finite RAO samples.
- Displaced mass differs from the specified 24,609,620 kg by -0.1505%. Computed transverse GM is 0.979997 m.
- Added mass and damping are symmetric, forbidden parity couplings are negligible, and effective mass is positive throughout the grid. Radiation damping is positive semidefinite within numerical roundoff (minimum normalized eigenvalue about -6.1e-10).
- Reconstructing the first-order impedance equation gives a maximum relative RAO residual of 4.2e-15.
- With the configured settings and the default 1.0 rad/s JONSWAP reference, time-constant feasibility checks pass. Natural heave/roll/pitch periods are 7.556, 18.014 and 7.090 s; damped periods are 7.631, 18.389 and 7.200 s, with damping ratios 0.140, 0.201 and 0.174.
- MSS first-order force and motion wave initialization passes with 100 frequencies and 24 directions. A 60-second heading-control simulation using the real MSS dynamics and RK4 remains finite (Hs 2 m, peak 1 rad/s, 0.05 s time step).

The raw motion RAOs include potential damping only and have a sharp roll resonance; additional MSS damping is not embedded in them. These checks establish numerical consistency and a limited runtime check, not experimental validation or mesh convergence.
