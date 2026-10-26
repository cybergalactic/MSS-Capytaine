# Tests and numerical comparisons

Run from the MSS-Capytaine repository root. The comparison scripts read saved results; they do not run either hydrodynamic solver or change MSS files.

## Documentation figures

Run `python tests/export_offset_csv_figure.py` to regenerate the offset-guide illustration from the actual testShip CSV. With MATLAB, run `addpath('tests'); exportGeometryPreview()` to regenerate the README geometry overview. Neither runs the BEM solver.

## Compare with MSS reference data

```sh
python tests/compare_tanker_wamit.py
python tests/compare_supply_shipx.py
```

The defaults use the sibling `MSS` checkout. Override reference locations with `--wamit path/to/tanker.mat` or `--shipx path/to/supply.mat`. Both scripts accept `--max-omega` to reduce the finite comparison band. The supply script also accepts `--output-dir` for a separate report directory.

Reports and overlay plots are kept with the cases:

- [Tanker versus WAMIT](../vessels_capytaine/wamit_gdf/tanker/validation/README.md).
- [Supply versus ShipX/VERES](../vessels_capytaine/shipx_mgf/supply/validation/README.md).

Inspect the scripts: [tanker/WAMIT](compare_tanker_wamit.py) and [supply/ShipX](compare_supply_shipx.py).

The comparisons align MSS FSD axes, actual matrix reference points, zero speed and headings. Tanker uses interpolation because its JSON frequency values are rounded. Supply uses the exact ShipX frequency grid up to 1.8 rad/s; no interpolation is permitted. Rerun an older supply solve after updating the JSON before comparing. Zero and infinity endpoints are excluded. The reported percentages are relative norm differences, not maximum pointwise errors or pass/fail tolerances. The supply reference uses strip theory and excludes surge from its error tables and plots because no meaningful strip-theory surge benchmark is available. Its report describes these limits and the uncertainty in comparing archived force phases.

## Python checks

The optional [Otter geometry development notes](otter/README.md) collect the source assumptions, loading/trim derivation and reproduction script. They are separate from the user-facing JSON/CSV cases; normal users need only `main.py`.

```sh
python -m unittest discover -s tests -v
```

Checks cover geometry readers, FSD transformations, CG/stability consistency, separated CSV demi-hulls, the Otter loading/inertia and lid-gap geometry, and comparison-frame transformations and direct frequency matching. They do not require a full BEM solve or the sibling MSS checkout.

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

With MSS on the path, also run `testDamping()`. It checks all six surface defaults at 0.8 and 1.0 rad/s, both Otter cases at 1.5 rad/s, the submerged LAUV case, the radiation-only boundary, zero radiation, combined errors for infeasible motions, and migration of old MAT metadata. Surface processing supports both older and current MSS damping interfaces and preserves the caller's warning settings. It requires no graphics or BEM solve.
