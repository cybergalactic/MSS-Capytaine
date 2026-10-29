# Capytaine vessel catalogue

Vessels are listed directly by name, followed by their geometry format. Each folder contains `config.json`, source geometry, a README and the main vessel MAT file. Supporting hydrostatics, generated panels and plots are stored under `results/` (or the configured `output_dir`).

| Folder | Geometry | MAT filename |
| --- | --- | --- |
| [testShip_csv](testShip_csv/) | Native CSV offsets | `testShip.mat` |
| [LAUV_marie_csv](LAUV_marie_csv/) | CSV offsets, submerged vehicle | `LAUV_marie.mat` |
| [otter_55kg_csv](otter_55kg_csv/) | Approximate twin hull, no payload | `otter_55kg.mat` |
| [otter_80kg_csv](otter_80kg_csv/) | Same hull, 25 kg payload | `otter_80kg.mat` |
| [supply_mgf](supply_mgf/) | Archived ShipX MGF sections | `supply.mat` |
| [s175_mgf](s175_mgf/) | Archived ShipX MGF sections | `s175.mat` |
| [tanker_gdf](tanker_gdf/) | WAMIT GDF panels | `tanker.mat` |
| [wigley_stl](wigley_stl/) | Public STL benchmark | `wigley.mat` |

The suffix identifies geometry, not the hydrodynamic solver. Every case is solved by Capytaine. Imported geometry requires no commercial solver licence. CSV remains the default native workflow.

```sh
python main.py --list
python main.py supply_mgf
python main.py s175_mgf
```

A unique JSON `body_name`, such as `supply` or `s175`, is also accepted. Use the full folder name when more than one geometry represents the same vessel. Direct JSON paths also work. `body_name` and the MAT filename retain the vessel identity without the format suffix.

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('s175')
plotVesselHydrodynamics('s175');
```

JSON contains geometry, loading, frequency, meshing and damping inputs. Derived vessel properties are computed; ShipX coefficient tables, roll-tank damping and drift forces are not appended. Current exports omit unused `vessel.main.S`, `Lwl` and `demihullSpacing`. Hull spacing remains in geometry configuration and the hydrostatics report. The remaining main fields support MSS dynamics, wave generation, loading checks or vessel displays.

All public geometry and results use MSS FSD. MGF/GDF keep their native file conventions and are converted internally. See the [JSON guide](../How%20to%20configure%20the%20vessel%20JSON%20file.md).

The MGF cases produce zero-speed data. The solver exports diagonal Capytaine mean drift in `vessel.drift`, included in all eight supplied vessel MAT files. The testShip case also includes a line plot of the signed coefficients at all 36 headings. Cross-frequency QTFs are not computed.

MATLAB/Octave post-processing uses the fixed JONSWAP reference with peak 1.0 rad/s and gamma 3.3 unless an explicit peak is supplied. Damping settings are example inputs, checked for nonnegative additional viscous damping at that reference. See the [release notes](../CHANGELOG.md) for the 1.2 folder and output-path migration.

Geometry readers supply geometry only; Capytaine solves every case. Additional geometry formats can follow the same interface. No ShipX result comparison is required.
