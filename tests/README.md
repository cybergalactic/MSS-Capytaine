# Tests and numerical comparisons

Run from the MSS-Capytaine repository root. The comparison scripts read saved results; they do not run either hydrodynamic solver or change MSS files.

## Documentation figures

Run `python tests/export_offset_csv_figure.py` to regenerate the offset-guide illustration from the actual testShip CSV. With MATLAB, run `addpath('tests'); exportGeometryPreview()` to regenerate the README geometry overview. Neither runs the BEM solver.

## Optional tanker comparison with WAMIT

```sh
python tests/compare_tanker_wamit.py
```

This reads saved Capytaine and WAMIT results and does not run either solver. The default reference is in the sibling MSS checkout; use `--wamit path/to/tanker.mat` to select another file and `--max-omega` to reduce the finite band. See the [report](../vessels_capytaine/tanker_gdf/validation/README.md) and [script](compare_tanker_wamit.py). Comparison aligns MSS FSD axes, CG reference, zero speed and headings, with interpolation for rounded frequencies. Reported differences are relative norms, not pointwise error bounds or experimental validation.

## Python checks

The optional [Otter geometry development notes](otter/README.md) collect the source assumptions, loading/trim derivation and reproduction script. They are separate from the user-facing JSON/CSV cases; normal users need only `main.py`.

```sh
python -m unittest discover -s tests -v
```

Checks cover geometry readers, FSD transformations, CG/stability consistency, separated CSV demi-hulls, the Otter loading/inertia and lid-gap geometry, local triangle refinement, catalogue lookup and the MAT export contract. They do not require a full BEM solve or the sibling MSS checkout.

Damping checks reject the old `kappa_126` field, validate time constants, and verify that the metadata-only MAT refresh preserves numerical arrays and cells.

## MATLAB/Octave reader checks

```matlab
addpath('matlab')
addpath('tests')
testSourceGeometry();
testDemihullGeometry();
testSavedMesh();
```

These check source formats, saved panel data, catalogue lookup with MSS names on the path, and missing results, without graphics or a solver run. Check the two-subplot presentation interactively with `plotVesselGeometry('supply')`.

Also run `testDamping()` and `testStandalonePeriods()`. They check catalogue damping targets, infeasible inputs, migration of old MAT metadata, coupled natural periods, and exclusion of the 10 rad/s infinity marker. Neither requires an MSS installation, graphics or a BEM solve.

The default processing reference is JONSWAP with peak **1.0 rad/s** and gamma **3.3**. With an optional sibling MSS checkout, also run `testPowerBasedSpectrum()` to compare both averaging implementations against dense independent integration of the wave-generation spectrum. It checks explicit peak overrides, the fixed default despite stale cached peaks, infinity-marker exclusion, floating/submerged damping and the roll-only increment.

## MGF geometry workflow

MGF tests check geometry import, draft clipping, coordinates and hydrostatics. Supply and S175 are solved with `python main.py supply_mgf` and `python main.py s175_mgf`; no ShipX results, comparison scripts or commercial solver licences are needed. Full BEM solves are intentionally separate from the unit checks.

Symmetry regression checks use a small hull to compare symmetric and expanded BEM radiation at zero, finite and infinite frequency, and diffraction at three headings. They also reject incompatible symmetry inputs before solving. These small checks do not run either full-size vessel.

`testMSSDampingConfig()` checks that the sibling MSS configuration reads all catalogue damping inputs from the saved MAT files, while retaining the FPSO and semisub WAMIT presets. It does not replace installed MSS vessel files.

Mean-drift tests check an analytic normalization case, agreement with the diagonal of Capytaine's full directional output, CG moment transfer, a fixed-cylinder BEM case with translation and angular-grid checks, and real signed MAT cells including directional parity. The complete small-hull pipeline is exercised without a full catalogue solve. To generate and inspect the testShip export, run:

```sh
python main.py testShip
python tests/check_mean_drift.py
python plot_results.py
```

The saved-file checker verifies export consistency, not experimental accuracy or hull-mesh convergence. The saved testShip file passes its layout, finite-value and heading-symmetry checks at 25 finite frequencies and 36 headings. The 0.5°/1° Kochin grid comparison gives relative norm differences below 1e-12 for surge/sway and about 0.27% for yaw. The drift figure shows three signed line plots with one curve per heading. Plotting reads the saved MAT file without a BEM solve.

Optional checks with the sibling MSS installation on the MATLAB/Octave path:

```matlab
addpath(genpath(fullfile(fileparts(pwd),'MSS')))
addpath(fullfile(pwd,'tests'))
testWamitMeanDrift
testWamitImport
```

These checks run in temporary folders and verify WAMIT drift coefficients,
MAT-file damping inputs and imports without a new hydrodynamic solve.
