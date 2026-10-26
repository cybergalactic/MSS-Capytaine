# How to configure the vessel JSON file

[Back to README](README.md) · [Generate offsets with AI](How%20to%20generate%20offsets%20points%20using%20AI.md)

A case consists of a JSON configuration and geometry in CSV, MGF, GDF or STL format. CSV and MGF describe hull sections from which panels are generated; GDF and STL describe existing surface panels. These files supply geometry, not computed hydrodynamic coefficients. All JSON coordinates use **MSS FSD: x forward, y starboard, z down**, with solver conversions handled internally.

## 1. Create a case folder

Choose an input group and create `vessels_capytaine/<group>/myVessel/`, using letters, numbers, underscores, or hyphens for the folder name. Use that same name for `body_name` and `output_filename` (with `.mat` appended) so the MATLAB catalogue helpers can locate the results.

Choose the closest starting point:

| Case | Geometry | Configuration |
| --- | --- | --- |
| Surface vessel | CSV offsets | [testShip/config.json](vessels_capytaine/capytaine_offsets_csv/testShip/config.json) |
| Submerged vehicle | CSV offsets | [LAUV_marie/config.json](vessels_capytaine/capytaine_offsets_csv/LAUV_marie/config.json) |
| Rigid twin hull, approximate Otter | One CSV demi-hull plus centerline spacing | [55 kg](vessels_capytaine/capytaine_offsets_csv/otter_55kg/config.json), [80 kg](vessels_capytaine/capytaine_offsets_csv/otter_80kg/config.json) |
| Existing hull panel mesh | WAMIT Geometric Data File (`.gdf`) | [tanker/config.json](vessels_capytaine/wamit_gdf/tanker/config.json) |
| ShipX hull sections | VERES MGF (`.mgf`) | [supply/config.json](vessels_capytaine/shipx_mgf/supply/config.json) |
| Existing triangle mesh | ASCII or binary STL (`.stl`) | [wigley/config.json](vessels_capytaine/mesh_stl/wigley/config.json) |

Copy the selected configuration into the new folder as `config.json`. Place the geometry file there too. Keep `config.json` inside the vessel folder, not in the format-group folder. File paths inside JSON are relative to that folder; the command can be run from the repository root.

If you need a CSV, follow [How to generate offset points using AI](How%20to%20generate%20offsets%20points%20using%20AI.md). The default group is `capytaine_offsets_csv`: the native, free, open-source workflow requiring no commercial licence. Other groups are `shipx_mgf`, `wamit_gdf`, and `mesh_stl`. A unique vessel name works in commands; `group/vessel` is explicit and also works if names repeat across groups. If a geometry file already exists in MGF, GDF or STL, skip AI offset generation.

## 2. Choose the vessel type

Use one complete case below. Surface vessels and submerged vehicles use different damping inputs; do not mix their arrays. The common field reference follows the examples, with [damping equations and worked examples](#viscous-damping-equations-and-examples) later in this guide.

### Case A: surface vessel

This complete example uses the synthetic testShip properties and a small frequency grid. Copy its [offset CSV](vessels_capytaine/capytaine_offsets_csv/testShip/offset_points.csv) into the same folder to try it. Replace geometry and physical properties when creating your own vessel; the example values are not general defaults.

```json
{
  "geometry_file": "offset_points.csv",
  "body_name": "myVessel",
  "output_filename": "myVessel.mat",
  "center_of_mass_m": [0.0, 0.0, 1.0],
  "mass_kg": 141954.8,
  "R44_m": 2.4194,
  "R55_m": 5.8156,
  "R66_m": 6.2470,
  "R46_m": 0.0,
  "samples_per_section": 24,
  "number_of_stations": 60,
  "viscous_damping": {
    "T_126": [20.0, 8.0, 20.0],
    "delta_zeta_345": [0.0, 0.1, 0.0]
  },
  "submerged": false,
  "water_density_kg_m3": 1025.0,
  "gravity_m_s2": 9.81,
  "omega_rad_s": [0.2, 0.5, 1.0],
  "output_dir": "results"
}
```

In this surface case, `T_126` specifies total velocity-decay time constants in surge, sway, and yaw: 20 s, 8 s, and 20 s. `delta_zeta_345` adds 10% of critical damping in roll, with no extra viscous term in heave or pitch. A surface vessel has hydrostatic heave restoring, so it uses a heave damping-ratio increment.

### Case B: submerged vehicle

This complete example uses the idealized LAUV Marie geometry and properties. Copy its [offset CSV](vessels_capytaine/capytaine_offsets_csv/LAUV_marie/offset_points.csv) into a new `vessels_capytaine/capytaine_offsets_csv/mySubmergedVehicle/` folder alongside this config. CG is 1 cm below the body origin in FSD. `submergence_depth_m` places that origin 1 m below the free surface without changing the body-fixed offsets.

```json
{
  "geometry_file": "offset_points.csv",
  "samples_per_section": 20,
  "number_of_stations": 50,
  "body_name": "mySubmergedVehicle",
  "output_filename": "mySubmergedVehicle.mat",
  "center_of_mass_m": [0, 0, 0.01],
  "mass_kg": 34.0,
  "R44_m": 0.053,
  "R55_m": 0.6218,
  "R66_m": 0.6218,
  "R46_m": 0.0,
  "viscous_damping": {
    "T_1236": [50.0, 5.0, 5.0, 5.0],
    "delta_zeta_45": [0.2, 0.2]
  },
  "submerged": true,
  "submergence_depth_m": 1.0,
  "water_density_kg_m3": 1025.0,
  "gravity_m_s2": 9.81,
  "omega_rad_s": [0.2, 0.5, 1.0],
  "output_dir": "results"
}
```

Here `T_1236` contains four time constants in **seconds**, ordered surge, sway, heave, yaw: 50 s, 5 s, 5 s, 5 s. There is no hydrostatic heave restoring for this freely submerged body, so heave also needs a time constant. `delta_zeta_45` adds 20% of critical damping in roll and pitch. This case supplies neither `T_126` nor `delta_zeta_345`.

### Common JSON rules and field reference

Use valid JSON: double quotes, lowercase `true`/`false`, no comments, and no trailing commas. Every length is in metres, mass in kilograms, angular frequency in rad/s, and full inertia entries in SI units.

| Field | What to enter |
| --- | --- |
| `geometry_file` | Required geometry path ending in `.csv`, `.mgf`, `.gdf` or `.stl`. |
| `draft_m` | MGF only: required positive design draft above the geometry baseline. It selects the wetted hull; it does not change JSON CG. |
| `demihull_spacing_m` | CSV only: optional positive centerline-to-centerline separation of two identical, locally symmetric hulls. Must exceed the single demi-hull beam. Omit for one hull. |
| `body_name` | Case name; defaults to the configuration folder name. |
| `output_filename` | A `.mat` filename; defaults to `<body_name>.mat`. |
| `center_of_mass_m` | Required FSD `[x, y, z]`, relative to the geometry origin. Also the rotation center. |
| `GM_T_m` | Optional independently known transverse GM in metres, used to check CG sign/origin before solving. Supplied for tanker and supply. |
| `GM_T_tolerance_m` | Optional absolute tolerance for that check; default is the larger of 0.05 m and 2% of the reference GM. |
| `mass_kg` | Required positive body mass. Compare it with displaced mass for equilibrium. |
| `R44_m`, `R55_m`, `R66_m` | Positive roll, pitch, yaw radii of gyration about CG. Required unless a full inertia matrix is supplied. |
| `R46_m` | Optional signed roll–yaw cross-radius in metres, default zero. Sets `MRB(4,6) = MRB(6,4) = mass_kg × R46_m × abs(R46_m)`. Supplied cases include it explicitly for consistency. |
| `inertia_matrix_csv` | Optional full 6-by-6 mass matrix at CG in MSS FSD and SI units; replaces the radii-based matrix. Omit `R46_m` when using this alternative. |
| `Lpp_m` | Optional positive length between perpendiculars. Defaults to mesh x extent; it does not rescale the hull. |
| `samples_per_section` | CSV/MGF only: integer at least 3, default 20; points along each half-section. |
| `number_of_stations` | CSV/MGF only: integer at least the number of input stations; omit to retain original wetted station positions. |
| `submerged` | Required boolean: `false` for a surface vessel, `true` for a submerged body. |
| `submergence_depth_m` | Required only when submerged: positive-down depth of the body-fixed origin. |
| `water_density_kg_m3` | Positive density; defaults to 1025.0. |
| `gravity_m_s2` | Positive gravity; defaults to 9.81. |
| `water_depth_m` | Omit or use `null` for infinite depth. Finite depth is currently unsupported by the complete limiting-frequency workflow. |
| `omega_rad_s` | Required nonempty list: every value must be positive, finite, and strictly below 10. Sorted and deduplicated internally. |
| `viscous_damping` | Required object: three surface or four submerged total velocity-decay time constants in seconds, plus restored-mode damping-ratio increments. See [equations and examples](#viscous-damping-equations-and-examples). |
| `output_dir` | Results folder relative to the config; defaults to `results`. |

For GDF/STL cases, select the mesh filename and omit the section-sampling settings. For MGF, select the MGF filename, set positive `draft_m` and configure section sampling. MGF import currently supports surface vessels only. For a submerged case, set `submerged` to `true`, supply `submergence_depth_m`, and replace the damping object with the submerged form shown below. No user setting selects solver axes.

## 3. Run and inspect

From the repository root:

```sh
python main.py capytaine_offsets_csv/myVessel
```

A direct path is also accepted: `python main.py path/to/config.json`. The results include `myVessel.mat`, `hydrostatics.json`, and the FSD panel files `generated_hull_panels.mat` and `generated_hull_panels.npz`.

Inspect the results in MATLAB or GNU Octave:

```matlab
addpath(genpath('/path/to/MSS'))
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('myVessel')
plotVesselHydrodynamics('myVessel');
```

For any input format, `plotVesselGeometry('myVessel')` previews the source at the top and the saved Capytaine hull mesh at the bottom, both in FSD. Before solving, the bottom subplot explains that no saved mesh is available. The mesh is read from `output_dir/generated_hull_panels.mat` (default `results/`); rerun Python after changing geometry or sampling to refresh it. Plotting itself does not generate a mesh and needs no MSS installation. Only `plotVesselHydrodynamics` requires MSS. These are the two public plotting functions; `matlab/private/` contains their support functions and needs no separate path setup. Check the mesh and all accuracy warnings before accepting hydrodynamic results. Changing mesh or physical inputs requires rerunning the solve; changing plotting options does not.

## Twin hulls and catamarans

For a rigid, port-starboard symmetric catamaran, use either one CSV demi-hull with `demihull_spacing_m`, or GDF/STL geometry containing both hulls at their actual centerline separation. The solver treats the complete assembly as one six-DOF rigid body and includes hydrodynamic interaction between its panels. Enter total mass, the combined FSD CG and inertia of the complete vessel about that CG; do not enter single-hull inertia or simply double its radii of gyration. Each hull's transverse offset contributes to the assembled roll and yaw inertia through the parallel-axis theorem.

Both wetted hulls need correct outward panel winding. STL must contain the complete port/starboard geometry in metres and MSS FSD; GDF symmetry flags may expand the native mesh as usual. Check hull spacing, the mesh between the hulls, and the generated internal free-surface lid. The open-water gap between the hulls must not be filled by artificial hull or lid panels. A coarse hull mesh may generate an incomplete or empty lid, so inspect/refine its resolution. The geometry preview shows exported hull panels, not the solver's internal lid.

For CSV, write the same `x_m,z_m,half_breadth_m` columns as for a monohull, but measure half breadth from the **local demi-hull centerline**. Do not add the hull separation to the CSV breadths. For example, these fields in an otherwise complete configuration build two hulls with centerlines at y = ±0.395 m:

```json
"geometry_file": "offset_points.csv",
"demihull_spacing_m": 0.79
```

The mesher first mirrors the local half-sections into one complete demi-hull, then places copies at y = ±spacing/2. It creates no connecting panels across the gap. Python and `plotVesselGeometry` both apply the spacing. CG and inertia in JSON always refer to the **whole vessel**, in the common FSD frame. Without the spacing field, the CSV continues to generate a single symmetric hull. The field is rejected for MGF, GDF and STL to avoid duplicating imported geometry. Spacing must exceed the demi-hull beam; touching or overlapping hulls are rejected.

The Otter cases use [55 kg without added payload](vessels_capytaine/capytaine_offsets_csv/otter_55kg/README.md) or [80 kg](vessels_capytaine/capytaine_offsets_csv/otter_80kg/README.md), MSS's 55 kg vessel plus a 25 kg payload at the body origin. Their assumed contours are fitted to displacement, waterplane area and longitudinal balance; they are not measured Otter hull lines. Both use the same physical hull; the 55 kg case solves its changed static heave and trim. Its CSV, CG and full inertia are rotated into the level equilibrium FSD frame together. The READMEs explain waterline origins, loading transforms and validation limits.

The MGF reader still supports a single monohull contour per station. ShipX/VERES supports multihulls, but this does not imply that every multihull MGF variant is handled by this importer. General asymmetric demi-hulls or different port/starboard hulls require complete panel geometry; the current heading mirroring also requires overall port/starboard symmetry.

Only the linear rigid-vessel workflow is covered. Independent hull motions, structural flexibility and empirical multihull resistance/maneuvering models are not added by importing a mesh.

## Geometry files

Use the same field for all four routes; the extension selects the reader:

| Example | Processing and coordinates |
| --- | --- |
| `"geometry_file": "offset_points.csv"` | Generate panels from your own CSV sections in MSS FSD. Default native workflow. |
| `"geometry_file": "supply.mgf"` | Read native ShipX/VERES sections; use `draft_m` to select the wetted hull, then generate panels. |
| `"geometry_file": "tanker.gdf"` | Native Capytaine GDF loading, including symmetry expansion. Leave WAMIT coordinates unchanged in the file. |
| `"geometry_file": "wigley.stl"` | Read ASCII/binary triangles in metres and MSS FSD, with outward winding and the full hull. |

When the extension is missing or unrecognized, Python and the preview reader inspect the file contents. Unidentifiable data produce an explicit error. Units and source axes are never inferred from geometry.

Paths are relative to `config.json`. `samples_per_section` and `number_of_stations` affect CSV/MGF only. GDF/STL reuse existing hull panels; refine those source files to change resolution. STL surface geometry above the waterline is clipped internally before the analysis; submerged STL bodies remain body-fixed until the configured depth translation is applied.

MGF contains full-height sections with heights relative to its baseline. `draft_m` is the design waterline height above that baseline; it is required because the geometry alone does not identify the loading condition. The importer preserves bulbous/flat-bottom contours and clips to an upright, untrimmed waterline. Dry stations are omitted. It does not import the entire `.shipx_ship` project archive; extract its MGF first. `supply.mgf` is already extracted in the example case.

STL carries no standard unit or axis declaration. Prepare the file in metres and MSS FSD at the intended origin. No automatic guessing or scaling is performed by the reader. For a surface vessel, submerged points have positive z and the waterline is z = 0. Existing mesh quality still limits useful frequencies. The [Wigley example](vessels_capytaine/mesh_stl/wigley/README.md) records the public source, preparation changes and licence.

The older `offset_points_csv` field is still an alias for `geometry_file`; provide only one of these two fields.

## Coordinates: MSS FSD

MSS FSD is always used; there is no coordinate-system option in JSON. JSON centers, offset CSV and STL coordinates, full inertia-matrix CSV inputs, exported panel files, hydrostatic centers, hydrodynamic matrices, RAOs, and plot axes use **x forward, y starboard, z down**. Python handles all solver sign changes internally. The configured center of mass is also the rotation center for the solve. The usual inputs are `mass_kg` and positive `R44_m`, `R55_m`, `R66_m`, with diagonal rotational entries `mass_kg × Rii_m²` about CG. Optional `R46_m` preserves roll–yaw cross-inertia without needing a separate matrix file:

```text
MRB(4,6) = MRB(6,4) = mass_kg × R46_m × abs(R46_m)
```

`R46_m` is a signed cross-radius, not a fourth positive principal radius. A negative value gives a negative matrix entry. Omission and explicit zero are equivalent. All catalogue examples include the field; only the Otter loads currently use nonzero values. For 80 kg, `R46_m ≈ 0.0927024811` gives the MSS cross-inertia of +0.6875 kg m². This is the actual matrix entry, not a product of inertia to which another minus sign should be applied. All gyration inputs are already in MSS FSD; do not change signs for Python.

For more general coupling, `inertia_matrix_csv` supplies a symmetric 6-by-6 rigid-body mass matrix about CG, ordered surge, sway, heave, roll, pitch, yaw in MSS FSD. It must contain six numeric rows without a header, be finite, symmetric and positive definite, and have its first three diagonal entries equal to `mass_kg`. Omit `R46_m` when selecting this full-matrix alternative. Optional `Lpp_m` specifies length between perpendiculars for exported `vessel.main.Lpp` and the block-coefficient reference length. Otherwise the exporter uses hull x extent. `Lwl` remains a mesh-extent estimate. The tanker sets `Lpp_m: 246` and selects the same 246 m GDF used by the MSS WAMIT run.

For a surface vessel, the origin is at midships on the design waterline: positive CG z is below the waterline and negative CG z is above it. The tanker uses `"center_of_mass_m": [3.93, 0, -2.5]`: 3.93 m forward and 2.5 m above the waterline. At 10 m draft, this is 12.5 m above the keel. This project's CG and CB are waterline-referenced FSD coordinates, not keel heights. For a submerged body, coordinates are relative to its body-fixed origin; `submergence_depth_m` separately specifies that origin's positive-down depth.

Offset CSV files contain `x_m,z_m,half_breadth_m`. Use FSD x and z and a nonnegative half breadth (distance to starboard). Rows within each station run from bottom to top, so z decreases; repeated z values represent flat bottom corners. Surface sections start at the centerline keel (z > 0) and end at the waterline (z = 0). Submerged sections run bottom to top with zero half breadth at both ends. Station groups may appear in any x order.

MGF/GDF files retain their standard source-file conventions. Supply the original file in metres without manually changing signs. The GDF loader handles symmetry; the MGF reader clips sections using `draft_m`. Both convert coordinates internally. The JSON CG is always MSS FSD regardless of geometry format. Inspect a GDF with `plotVesselGeometry('wamit_gdf/tanker')`: the top subplot previews the source before solving, and the bottom displays saved solver panels when available. No MSS installation is needed. Both triangles and quadrilaterals are supported. See [the tanker case](vessels_capytaine/wamit_gdf/tanker/README.md) for the supplied mesh's length variants and load-condition notes.

### CG and transverse-stability consistency check

Python checks roll stability before any BEM calculation. Non-positive calculated transverse GM stops the run and asks you to check FSD CG signs, vertical origin, geometry and loading condition.

When an independent GM is known, supply `GM_T_m`; for the tanker the MSS reference is 9.9823 m. Python compares it with GM calculated from the mesh and configured CG. A disagreement beyond `GM_T_tolerance_m` stops the solve and reports the approximate FSD zG implied by that hull and reference GM. The default tolerance is max(0.05 m, 2% of GM_T_m). Use a tolerance justified by the source accuracy and mesh convergence, not one chosen just to hide a sign error. The reference is a check, not a hydrostatic-stiffness override.

Positive GM alone does not prove that zG has the correct sign. A vessel can remain stable with either of two different CG heights. For an independent check, use GM from a trusted loading-condition calculation or measurement; do not copy the mesh's own calculated GM into JSON as if it were independent. Some demonstration cases have no independent reference, so they receive only the positive-GM check. Check these coordinate conventions when migrating older configurations.

## Frequencies, headings, and submerged depth

`vessels_capytaine/capytaine_offsets_csv/testShip/config.json` specifies the mesh resolution, mass, radii of gyration (or an inertia-matrix CSV), center of mass, and wave frequencies. The center of mass is also used as the rotation center, so all hydrodynamic matrices, forces, and motions are referenced to the CG. The solver always uses 19 wave directions from 0° to 180° in 10° increments; the MSS export workflow fixes these, so they are not configuration inputs. The limiting-frequency calculations assume infinite water depth, the default when `water_depth_m` is absent or `null`.

Set `submerged` to `false` for a surface vessel. The solver then generates an internal waterplane lid to suppress artifacts at irregular frequencies. Set it to `true` for a submerged vehicle such as an AUV; no lid is generated.

A submerged case also requires `submergence_depth_m`, the positive-down depth of the body-fixed origin used as the fixed equilibrium position for the free-surface solve. The complete hull must remain below the free surface. It is not an average of hydrodynamic coefficients over a depth range. Each value therefore defines a separate Capytaine model: small values give stronger free-surface interaction and wave excitation, while increasing the value progressively reduces wave excitation and makes the radiation coefficients approach their deeply submerged limits. Large values such as 20 m or 100 m are accepted, but much of the first-order wave excitation may then be extremely small; long, low-frequency waves penetrate more deeply than short, high-frequency waves. This translation affects the added mass, radiation damping, and RAOs but is removed from the exported MSS CG, CB, and panel coordinates.

For the `LAUV_marie` example, 1.0 m is the recommended reference depth when one generated model will be reused over a larger operating-depth range in MSS. The 0.15 m diameter hull remains fully submerged at this center depth, and choosing the shallowest intended operating depth avoids extrapolating the wave loads upward from a deeper model. In the MSS runtime model, the generated force RAOs are used without amplification when `d <= d_ref`, while greater depths may use the one-sided deep-water attenuation

```text
exp[-(omega^2/g) (d - d_ref)],    d > d_ref,
```

where `d_ref = vessel.main.submergenceDepth`. This runtime approximation attenuates the wave-excitation force only; it does not recompute the depth-dependent added mass or radiation damping. Users whose vehicle operates mainly at another depth can instead set `submergence_depth_m` to that representative depth and regenerate the complete Capytaine model. For submerged bodies, `vessel.main.T` is the body height rather than a surface-vessel draft. Use `output_filename` to give each catalog case its own `.mat` filename.

For CSV/MGF geometry, `samples_per_section` controls resolution around each half section, while `number_of_stations` controls resolution along the hull. Check Capytaine's mesh-resolution warnings at the highest wave frequencies. The workflow checks the reciprocal `A_ij`/`A_ji` and `B_ij`/`B_ji` coupling terms before symmetrizing the added-mass and radiation-damping matrices as required at zero speed. It reports a warning when the local reciprocity error exceeds 1%, and the skew also exceeds 0.1% of the largest matrix norm over the solved frequency range. This avoids misleading ratios where submerged-body radiation damping is numerically close to zero. For GDF/STL geometry, mesh resolution is fixed by the file. If a warning is reported, refine the mesh and rerun the calculation before accepting the exported model.

## Viscous damping: equations and examples

These inputs describe additional damping for the MSS maneuvering model. They are not all percentages: the surface-vessel inputs are dimensionless, while the submerged vehicle's four time constants are in **seconds**. The numbers in each field name identify the motions in its array:

| Index | Motion | Velocity |
| --- | --- | --- |
| 1 | Surge: forward/backward translation | u |
| 2 | Sway: sideways translation | v |
| 3 | Heave: vertical translation | w |
| 4 | Roll: rotation about x | p |
| 5 | Pitch: rotation about y | q |
| 6 | Yaw: rotation about z | r |

### What MSS computes

Capytaine supplies frequency-dependent added mass A(ω) and potential radiation damping B(ω). MSS `computeManeuveringModel` averages these with a normalized wave-spectrum weight to obtain constant matrices A_eq and B_eq. It then forms the effective inertia M and total linear damping D:

$$
M = M_{RB} + A_{eq}, \qquad D = B_{eq} + B_v.
$$

At zero forward speed, Capytaine supplies zero damping at the zero- and infinite-frequency limits. The exported 10 rad/s endpoint is a plotting marker for infinity, not a finite-frequency sample. A finite curve that is still rising at its highest computed frequency has not demonstrated convergence to that limit. Do not interpolate to the infinity marker or include it as a finite point in a damping integral; the unresolved tail needs its own check.

The weighting peak `omega_p` belongs to MSS post-processing, not the geometry JSON or the BEM solve. Select the wave-spectrum peak for the sea state being represented; do not move it to a vessel's radiation-damping peak to obtain larger equivalent damping. It does not change the saved frequency-dependent A/B data. The plotting helper accepts an optional peak frequency in rad/s:

```matlab
vessel = plotVesselHydrodynamics('otter_55kg', 1.5);
```

Omitting it uses 0.8 rad/s as an illustrative default, not a universal vessel constant. Use the local plotting helper to apply the current catalogue damping inputs; its returned vessel contains the processed model. At zero speed, vessels in the same sea state see the same wave-spectrum peak. Forward speed requires appropriate treatment of encounter frequencies and hydrodynamic coefficients; changing `omega_p` alone does not add forward-speed hydrodynamics. Small equivalent damping under low-frequency weighting does not imply small radiation damping at a higher roll resonance. No BEM rerun is needed to change this weighting.

Here B_v is the **additional diagonal viscous damping** calculated from the JSON settings. M includes added mass; it is not just the dry rigid-body mass. G denotes the hydrostatic restoring matrix. Translation damping coefficients have units N s/m; rotational damping coefficients have units N m s/rad. Positive damping opposes velocity: its force/moment is −Dν.

There are two kinds of motion to specify: unrestrained motions decay in velocity without a hydrostatic spring, while restored motions oscillate against hydrostatic stiffness. Their damping inputs therefore use different equations. The scalar equations below explain each diagonal setting; the full model can also contain coupling between motions.

### Surface vessels: three velocity-decay time constants

For `submerged: false`, use this JSON object (the testShip example):

```json
"viscous_damping": {
  "T_126": [20.0, 8.0, 20.0],
  "delta_zeta_345": [0.0, 0.1, 0.0]
}
```

`T_126 = [T₁, T₂, T₆]` follows **surge, sway, yaw**, in **seconds**. These are targets for total diagonal damping, including radiation, rather than viscous-only time constants. For one uncoupled velocity component:

$$
M_{ii}\dot{\nu}_i+D_{ii}\nu_i=\tau_i,\qquad
T_i=\frac{M_{ii}}{D_{ii}}.
$$

Without applied force, velocity falls to approximately 37% after Tᵢ. For a small force or moment step from rest, it reaches approximately 63% of its final value after Tᵢ. Use small-amplitude responses to estimate these values; quadratic drag, coupling and actuator dynamics affect larger or coupled responses. They are not position or heading return times.

The additional viscous damping is chosen so total diagonal damping meets that target:

$$
D_{ii}=\frac{M_{ii}}{T_i},\qquad
(B_v)_{ii}=\frac{M_{ii}}{T_i}-(B_{eq})_{ii},\qquad i\in\{1,2,6\}.
$$

For example, M₁₁ = 100 kg, B_eq,11 = 5 N s/m and T₁ = 5 s give D₁₁ = 20 N s/m and B_v,11 = 15 N s/m.

**Feasibility:** additional viscous damping must be nonnegative. For positive B_eq,ii, the potential-flow equivalent model permits:

$$
0<T_i\leq T_{i,max}=\frac{M_{ii}}{(B_{eq})_{ii}}.
$$

At Tᵢ = Tᵢ,max, radiation damping alone meets the requested time constant. If B_eq,ii = 0 there is no finite upper limit. A larger requested T would require cancelling radiation damping. MATLAB reports all offending motions with their requested and maximum values and stops; it never silently clamps inputs. A measured time constant outside this limit flags incompatibility with the selected equivalent model, not proof that the measurement is wrong. The maxima are recomputed whenever the post-processing weighting changes. Coupled modal decay times can differ from these diagonal targets.

The supplied surface examples use these **starting estimates**, checked against the saved results at the illustrative 0.8 rad/s default and at 1.0 rad/s:

| Case | T_126 [s]: surge, sway, yaw | Basis |
| --- | --- | --- |
| testShip | [20, 8, 20] | Synthetic example estimates. |
| Supply | [20, 2.5, 5] | Example estimates within this equivalent model's limits; not transferred from a different supply ship's sea trials. |
| Tanker | [45, 3, 4] | Example estimates within this equivalent model's limits; not calibrated tanker damping. |
| Wigley | [80, 8, 12] | Idealized STL example estimates. |
| Otter 80 kg | [5, 0.5, 1] | Starting values from the sea-trial-calibrated MSS Otter model; transfer to the approximate wave model still needs validation. |
| Otter 55 kg | [5, 0.5, 1] | Same starting values, without an independent 55 kg calibration. |

For example, supply T_max is approximately [30.1, 4.29, 9.49] s at the 0.8 rad/s processing default. The older comparison values [88, 45, 12] s cannot be used as total targets for this equivalent model at that setting. No value in this table is a universal vessel-class default.

### Restored motions: increments in damping ratio

For a surface vessel, `delta_zeta_345 = [Δζ₃, Δζ₄, Δζ₅]` follows the order **heave, roll, pitch**. For a submerged vehicle, `delta_zeta_45 = [Δζ₄, Δζ₅]` follows **roll, pitch**. Both use the same equation.

For an uncoupled displacement or angle ξᵢ with positive stiffness Gᵢᵢ:

$$
M_{ii}\ddot{\xi}_i + D_{ii}\dot{\xi}_i + G_{ii}\xi_i = 0,
$$

$$
D_{crit,i} = 2\sqrt{M_{ii}G_{ii}}, \qquad
\zeta_i = \frac{D_{ii}}{D_{crit,i}}.
$$

D_crit is the critical damping coefficient. ζ is a dimensionless damping ratio: values below 1 describe an underdamped scalar oscillator. The input Δζ specifies how much viscous damping to add as a fraction of D_crit:

$$
(B_v)_{ii} = 2\,\Delta\zeta_i\sqrt{M_{ii}G_{ii}},
\qquad
\zeta_{total,i} = \zeta_{potential,i} + \Delta\zeta_i.
$$

**An entry of 0.1 adds 10% of critical damping**, or 0.10 to the damping ratio. It does not multiply the existing damping by 1.1 and does not set the total damping ratio to 0.1. For example, if the potential damping ratio is 0.05, adding Δζ = 0.1 gives a total scalar damping ratio of 0.15. For an illustrative heave oscillator with M₃₃ = 1000 kg and G₃₃ = 4000 N/m, D_crit = 4000 N s/m, so Δζ₃ = 0.1 adds 400 N s/m.

Consequently, the surface example `[0.0, 0.1, 0.0]` adds no viscous damping in heave or pitch and adds 10% of critical damping in roll. Potential damping still acts in all three. For a submerged vehicle, `[0.2, 0.2]` adds 20% of critical damping in roll and pitch. These are scalar diagonal interpretations; the coupled system's modal damping ratios can differ.

### Submerged vehicles: four velocity-decay time constants

For `submerged: true`, replace the surface object with:

```json
"viscous_damping": {
  "T_1236": [50.0, 5.0, 5.0, 5.0],
  "delta_zeta_45": [0.2, 0.2]
}
```

`T_1236 = [T₁, T₂, T₃, T₆]` contains **four time constants in seconds**:

| Array entry | Motion | Example | Scalar interpretation |
| --- | --- | --- | --- |
| First: T₁ | Surge | 50 s | Forward velocity falls to about 37% after 50 s. |
| Second: T₂ | Sway | 5 s | Sideways velocity falls to about 37% after 5 s. |
| Third: T₃ | Heave | 5 s | Vertical velocity falls to about 37% after 5 s. |
| Fourth: T₆ | Yaw | 5 s | Yaw rate falls to about 37% after 5 s. |

Why four? A freely submerged, constant-volume vehicle in homogeneous water has no linear hydrostatic restoring in **surge, sway, heave, or yaw**. In particular, moving it a little vertically does not change its volume or buoyancy, so G₃₃ = 0. Heave is therefore specified with a velocity-decay time, whereas surface-vessel heave has waterplane restoring and uses Δζ₃. Roll and pitch can still have restoring from the separation of CG and CB; they use `delta_zeta_45`. This assumes no tether, spring, or other added restoring mechanism.

For one uncoupled, unforced velocity component νᵢ (u, v, w, or r):

$$
M_{ii}\dot{\nu}_i + D_{ii}\nu_i = 0,
\qquad
\nu_i(t) = \nu_i(0)e^{-t/T_i},
\qquad
T_i = \frac{M_{ii}}{D_{ii}}.
$$

After one time constant, about 37% of the initial velocity remains; after three, about 5% remains. These times describe velocity or yaw-rate decay, not time to return to the original position or heading. Smaller T means stronger damping for the same effective inertia.

MSS chooses B_v so the **total** diagonal damping has the target time constant:

$$
D_{ii} = \frac{M_{ii}}{T_i}, \qquad
(B_v)_{ii} = \frac{M_{ii}}{T_i} - (B_{eq})_{ii},
\qquad i \in \{1,2,3,6\}.
$$

For example, if M₁₁ = 100 kg, B_eq,11 = 5 N s/m, and T₁ = 5 s, then D₁₁ = 20 N s/m and MSS adds B_v,11 = 15 N s/m. Using T₁ = 10 s gives D₁₁ = 10 N s/m and adds only 5 N s/m.

The added viscous coefficient must be nonnegative. When B_eq,ii > 0:

$$
0 < T_i \leq \frac{M_{ii}}{(B_{eq})_{ii}}.
$$

In that example, the maximum allowed T₁ is 20 s. A request for 25 s would require negative viscous damping, so MSS rejects it. If B_eq,ii = 0, any positive finite T produces a positive viscous coefficient. Coupling can make the full model's modal decay times differ from these diagonal targets.

### Choosing values and applying them

Choose time constants using small force-step responses, free velocity decay, measured histories or another justified reference. Use seconds for total velocity-decay targets and damping-ratio increments for restored motions. All entries must be finite; T must be strictly positive and Δζ nonnegative. Use only the pair appropriate to the vessel type. The object name `viscous_damping` is retained, but its T entries target total damping while its Δζ entries add viscous damping to restored motions.

Python stores the inputs under `vessel.powerBased`; MATLAB checks the maxima when constructing B_v and D. The private helper used by `plotVesselHydrodynamics` supports both the older percentage-based MSS surface interface and the current T_126 interface. It applies the saved time constants and reports all infeasible targets together. Do not pass a T_126 vector directly to an older MSS function's percentage argument. Obtain the processed vessel through:

```matlab
vessel = plotVesselHydrodynamics('supply');
```

The helper prints specified time constants beside T_max. Its returned vessel contains M, D, G, and `powerBased.T_max_126` (or `T_max_1236` for submerged cases). The selected physical wave-spectrum weighting affects the equivalent matrices and their limits; it is not an extra JSON field.

After changing **only** damping settings, refresh an existing MAT file:

```sh
python main.py shipx_mgf/supply --update-damping
```

This replaces damping metadata without a BEM solve, preserves the raw hydrodynamic arrays, and clears derived damping/model caches if present. It does not apply changes to geometry, mass, CG, inertia or the frequency grid; those require an appropriate new solve/export. Repeat MATLAB processing after the refresh. The Python-exported motion RAOs use potential damping only, so the damping settings do not change those saved RAOs.

### Migrating older JSON and MAT files

Replace the old surface percentage field:

```json
"kappa_126": [0.05, 0.05, 0.05]
```

with positive total time constants in seconds, for example for this supply case:

```json
"T_126": [20.0, 2.5, 5.0]
```

Keep `delta_zeta_345`. **Do not copy the percentage numbers into T_126.** Python rejects `kappa_126` with migration instructions before meshing/solving; MATLAB likewise rejects stale MAT damping metadata. Run `--update-damping` after editing the JSON. Submerged `T_1236` and `delta_zeta_45` are unchanged. This damping migration needs no additional JSON version or coordinate-system field. JSON coordinates remain MSS FSD; check the [coordinate conventions](#coordinates-mss-fsd) when migrating older files.
