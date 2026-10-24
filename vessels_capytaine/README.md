# Capytaine vessel catalogue

Choose an input route, then a vessel. Each vessel folder contains `config.json`, its geometry, a README and, after solving, a `results/` folder.

| Group | Input and processing | Cases |
| --- | --- | --- |
| `capytaine_offsets_csv` | **Default workflow: free and open source.** Generate panels from your own CSV hull offsets; optional spacing duplicates a demi-hull. No commercial licence is required. | [testShip](capytaine_offsets_csv/testShip/), [LAUV_marie](capytaine_offsets_csv/LAUV_marie/), Otter at [55 kg](capytaine_offsets_csv/otter_55kg/) and [80 kg](capytaine_offsets_csv/otter_80kg/) — two loads of the same approximate twin hull |
| `mesh_stl` | Load an existing STL triangle mesh. | [wigley](mesh_stl/wigley/) — adapted public benchmark geometry |
| `shipx_mgf` | Generate panels from ShipX/VERES MGF hull sections. | [supply](shipx_mgf/supply/) |
| `wamit_gdf` | Load an existing WAMIT Geometric Data File panel mesh. | [tanker](wamit_gdf/tanker/) |

These group names sort alphabetically, placing the default CSV workflow first. The imported formats can also be solved with Capytaine without running the commercial software that originally produced the geometry.

From the repository root:

```sh
python main.py --list
python main.py                         # default: CSV testShip
python main.py capytaine_offsets_csv/testShip
python main.py capytaine_offsets_csv/LAUV_marie
python main.py capytaine_offsets_csv/otter_55kg
python main.py capytaine_offsets_csv/otter_80kg
python main.py wamit_gdf/tanker
python main.py shipx_mgf/supply
python main.py mesh_stl/wigley
```

A unique vessel name such as `tanker` is also accepted. If the same vessel appears in several groups, use `group/vessel` to select it explicitly. A direct path to a JSON configuration is accepted too.

Every JSON uses `geometry_file`; its extension selects the reader. If no known extension is present, the reader inspects the file contents. JSON centers and full inertia matrices, CSV offsets, STL vertices, exported panels, and plots use MSS FSD: x forward, y starboard, z down. GDF/MGF files retain their native conventions, which are converted internally. MGF additionally requires `draft_m`, measured above its baseline. See the root [JSON guide](../How%20to%20configure%20the%20vessel%20JSON%20file.md) for all four routes and the surface/submerged damping equations.

Inspect any source geometry before solving (MSS is not required):

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('capytaine_offsets_csv/testShip')
plotVesselGeometry('wamit_gdf/tanker')
plotVesselGeometry('shipx_mgf/supply')
plotVesselGeometry('mesh_stl/wigley')
```

Inspect results after solving (MSS is needed for hydrodynamic post-processing):

```matlab
plotVesselGeometry('capytaine_offsets_csv/testShip')
addpath(genpath('/path/to/MSS'))
plotVesselHydrodynamics('capytaine_offsets_csv/testShip');
```

The MATLAB functions accept either `group/vessel` or a unique vessel name. `plotVesselGeometry` shows source geometry at the top and the saved Capytaine hull mesh at the bottom. Before a solve, the mesh subplot shows a message; rerun Python after changing geometry or sampling to refresh the saved panels. The two plotting functions are public; readers and catalogue lookup live in `matlab/private/` and need no separate path setup. Folder name, JSON `body_name` and result filename agree within each vessel folder: `supply`, `supply`, `supply.mat`, for example. The group supplies the input-format context; a suffix on the vessel name is unnecessary.

To copy results into MSS later, use the matching vessel directory under `MSS/HYDRO/vessels_capytaine/`. Updating this catalogue does not update MSS.

For rigid symmetric twin hulls, use either a complete GDF/STL assembly or one locally symmetric CSV demi-hull with `demihull_spacing_m`. Always specify the total vessel mass, combined CG and full-vessel inertia. The Otter contours are assumed, fitted to simple hydrostatic constraints rather than measured hull lines. MGF multihull variants remain unsupported; see the [twin-hull guidance](../How%20to%20configure%20the%20vessel%20JSON%20file.md#twin-hulls-and-catamarans).
