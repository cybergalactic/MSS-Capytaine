# Supply vessel — ShipX MGF example

Run from the repository root:

```sh
python main.py shipx_mgf/supply
```

Preview the source in MATLAB/Octave before solving:

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('shipx_mgf/supply')
```

The top subplot shows source sections; the bottom shows saved Capytaine hull panels after solving, or a message before results are available. Both use FSD.

`supply.mgf` is copied byte-for-byte from `shipx_geometry.mgf` inside `MSS/HYDRO/vessels_shipx/supply/Supply.shipx_ship`. The archive contains one MGF run input, so the source selection is unambiguous. This is geometry, not a file of precomputed hydrodynamic coefficients.

The MGF holds 41 full-height half sections. Python clips them to the JSON `draft_m = 6`, drops dry sections, resamples the wetted contours and generates panels. `samples_per_section` and `number_of_stations` control resolution, as for CSV. End closures are subdivided across the transom to avoid oversized full-beam triangles. The maximum hull panel radius is about 1.47 m, below wavelength/8 at the configured 1.8 rad/s upper limit. Native ShipX coordinate conversion is internal; JSON CG and all public results use MSS FSD with the surface-vessel origin at the design waterline. The importer supports classic, upright monohull MGF sections; it does not apply trim, sinkage, or import the entire ShipX project.

The source Lpp is 82.8 m and breadth is 19.2 m. The loading condition follows the MSS ShipX supply example: mass 6,362,208.5 kg and FSD CG approximately `[-5.38592, 0, -1.30000]` m. Radii of gyration about CG are taken from `input.hyd`: 6.72 m, 20.7 m, 20.7 m. The stored MSS vessel has slightly different radius metadata; this case uses the documented source inputs. The default damping increments are assumptions, not calibrated supply-vessel viscous coefficients. See the [JSON guide](../../../How%20to%20configure%20the%20vessel%20JSON%20file.md#viscous-damping-equations-and-examples).

With the current 24-point sections and 70 stations, the mesh contains about 6,671 panels and displaces about 6,189 m3, versus approximately 6,207 m3 in the ShipX reference. Refine the sampling to check convergence before using this as an accuracy benchmark. New supply hydrodynamic results require a solve; the catalogue reorganization does not itself calculate them.

Inspect the generated results after solving; hydrodynamic post-processing requires MSS on the MATLAB/Octave path:

```matlab
plotVesselGeometry('shipx_mgf/supply')
addpath(genpath('/path/to/MSS'))
plotVesselHydrodynamics('shipx_mgf/supply');
```

Files are written to `results/supply.mat`, `results/hydrostatics.json`, and the FSD panel files in this case folder.

## Comparison with MSS ShipX

From the repository root, run:

```sh
python tests/compare_supply_shipx.py
```

Inspect the [Python comparison test](../../../tests/compare_supply_shipx.py) to see the frame transformations, direct sample matching and error calculation.

This uses the saved Capytaine solve and the reference at `MSS/HYDRO/vessels_shipx/supply/supply.mat` in the sibling checkout. Use `--shipx path/to/supply.mat` for another location. It does not rerun Capytaine or modify MSS. Plots, error tables and frame diagnostics are saved in [validation](validation/README.md).

The JSON uses the exact 33 ShipX frequency samples from 0.1047198 to 1.795196 rad/s. After changing from the older grid, rerun `python main.py shipx_mgf/supply` before comparing. The comparison selects matching reference samples directly and rejects grids requiring interpolation. The included results have been regenerated on this grid. Excitation-force amplitude norm differences span about 4.9–18.3% across the five compared DOFs; added mass and damping show larger differences. These are sample-norm differences, not maximum pointwise errors or experimental validation. It shifts the ShipX reference from the waterline above CG to the actual CG, using rigid-body mass coupling to identify the offset. The archived MSS reference is labelled CG despite retaining this vertical offset. After shifting, the rigid-body inertias agree.

This ShipX run uses ordinary strip theory, while Capytaine uses a 3D panel model. Surge is excluded from the error tables and plotted comparisons: MSS approximates surge added mass and damping by scaling sway data, which provides no meaningful surge benchmark. Reported force errors distinguish amplitude from complex phase differences; no wave-phase offset is fitted. Motion RAOs are excluded because the ShipX run includes viscous damping and roll tanks. Treat this as a diagnostic comparison; it does not establish experimental accuracy or mesh convergence.

## Velocity-decay time constants

`T_126 = [20, 2.5, 5]` seconds specifies total surge, sway and yaw damping targets. MATLAB displays the potential-flow T_max values and rejects infeasible requests; it does not clamp them. These are example estimates within the equivalent model limits, not sea-trial values transferred from Northern Clipper. At the 0.8 rad/s default, T_max is about [30.1, 4.29, 9.49] s. After damping-only edits, run `python main.py shipx_mgf/supply --update-damping` to refresh the MAT settings without solving again.
