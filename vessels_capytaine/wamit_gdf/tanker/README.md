# MSS tanker GDF case

Run from the repository root with `python main.py wamit_gdf/tanker`.

`geometry_file` selects the supplied `tanker.gdf`. Capytaine's native GDF loader expands its port/starboard symmetry to 4750 panels. No hull panels are generated or resampled. The existing internal waterplane lid is still generated for the surface-vessel solve.

MSS FSD is always used. JSON centers, output panels, hydrodynamic results, and plots use x forward, y starboard, z down. Supply the original GDF unchanged; Python handles its native file convention and all solver-axis conversions internally.

The configured CG `[3.93, 0, -2.5]` is FSD relative to the waterline: 3.93 m forward and 2.5 m above the waterline, or 12.5 m above the keel at 10 m draft. This project's exported CG is the same FSD vector, whereas the legacy MSS tanker.main reports vertical centers as heights above the keel.

Inspect the original geometry before solving through the unified FSD preview (no solve, result file or MSS installation needed):

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('wamit_gdf/tanker')
```

`plotVesselGeometry` shows the source in the top subplot without a solve or MSS installation. After solving, the same call shows the saved Capytaine hull panels in the bottom subplot too. Both use MSS FSD.

Mass, radii of gyration, density and gravity follow the supplied MSS data. The 41 positive finite frequencies retain the MSS points through 1.8003 rad/s; higher finite frequencies were removed after the initial resolution warnings. Zero and the 10 rad/s marker for infinite frequency are added by the existing solver. Damping inputs use `T_126 = [45, 3, 4]` seconds for total surge, sway and yaw damping. These are example estimates, not calibrated tanker values or values derived from the former 5% rule. They pass the computed T_max checks at 0.8 and 1.0 rad/s. At 0.8 rad/s, T_max is approximately [51.9, 3.90, 5.51] s. `delta_zeta_345 = [0, 0.1, 0]` adds 10% of critical damping in roll, with no extra viscous heave or pitch damping. After damping-only JSON edits, use `python main.py wamit_gdf/tanker --update-damping` without a BEM rerun. See the [equations and migration instructions](../../../How%20to%20configure%20the%20vessel%20JSON%20file.md#viscous-damping-equations-and-examples). Volume, CB, GM, and dimensions are computed from the supplied geometry.

## Length and WAMIT validation

This example selects the **246 m** `tanker.gdf` copied unchanged from `MSS/HYDRO/vessels_wamit/tanker/tanker.gdf`. Its JSON sets `Lpp_m` to 246 m. The unused historical 252 m GDF copies have been removed from MSS; this case contains only the selected 246 m GDF.

MSS `fnames.wam` selects `tanker.gdf`; `tanker.out` and `wamitlog.txt` record that the 24 June 2008 WAMIT run used this 246 m file. Its vertex set matches the 252 m variant with x scaled by 246/252 to within 0.1 mm rounding. Breadth and draft are unchanged. The files do not explain why this longitudinal change was made.

The selected mesh spans about 256.864 m in x, with 46 m breadth and 10 m draft. Total mesh extent is not length between perpendiculars. `Lpp_m` sets the exported `vessel.main.Lpp` and block-coefficient reference length; `Lwl` retains the exporter’s mesh-extent estimate. Capytaine computes about 92337 m3 displacement versus WAMIT’s reported 92312 m3; non-planar-panel integration can differ. The original GDF panels trigger Capytaine warnings about non-coplanar and non-convex faces, so inspect accuracy warnings before accepting the results, particularly at the highest frequencies.

After solving, compare added mass, radiation damping, and complex force RAOs against the saved MSS WAMIT results at matching frequencies and headings. Put both sets in MSS FSD at CG and reconcile normalization, wave-heading and phase conventions. Legacy MSS vertical centers use keel heights; this project uses the waterline origin. The selected hull now matches the WAMIT geometry; the numerical solvers and internal lid construction may still produce differences.

## Completed comparison

The existing Capytaine solve has been compared with the saved MSS WAMIT results over 0.0424–1.8003 rad/s after shifting the reference quantities to CG. This is the project's primary BEM-to-BEM benchmark: WAMIT and Capytaine both model three-dimensional radiation and diffraction. The close agreement supports the numerical workflow and its MSS coordinate, reference-point and phase conversions. See the [comparison report and overlay plots](validation/README.md) and the [Python comparison test](../../../tests/compare_tanker_wamit.py). Diagonal added-mass relative norm errors are below 0.7%, radiation-damping errors are below 1.8%, and complex force-RAO errors pooled over all headings are below 0.7%. These are numerical comparisons over this band, not bounds on every sample or experimental validation. Higher frequencies remain subject to the mesh warnings. Reproduce with `python tests/compare_tanker_wamit.py`.
