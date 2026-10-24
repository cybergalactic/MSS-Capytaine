# Wigley hull — public STL example

This case imports a prepared version of a public Wigley hull STL. The Wigley hull is an idealized benchmark hull, not a particular commercial ship.

```sh
python main.py mesh_stl/wigley
```

Preview the source in MATLAB/Octave before solving:

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('mesh_stl/wigley')
```

The top subplot shows the source STL; the bottom shows saved Capytaine hull panels after solving, or a message before results are available. Both use FSD.

## Source and licence

- Source: [Wigley hull STL on Wikimedia Commons](https://commons.wikimedia.org/wiki/File:Wigley_hull.stl).
- Credit: OpenFOAM Foundation; Wikimedia Commons upload by Ddcampayo. The file page identifies the hull as originating in OpenFOAM's geometry resources.
- Mesh licence: [Creative Commons Attribution-ShareAlike 4.0](https://creativecommons.org/licenses/by-sa/4.0/).
- This adapted `wigley.stl` is distributed under the same CC BY-SA 4.0 licence. This notice applies to the mesh adaptation, not to the application's code.
- `provenance.txt` records the source URL, hashes, modifications and panel count.

## Preparation and assumptions

The original file contains 12,896 triangles and uses upward-positive vertical coordinates. Preparation reversed its inward panel winding, clipped it to the source waterline, simplified one half to about 900 triangles using fast-simplification 0.2.0, and mirrored it for exact port/starboard symmetry. The decimated waterline boundary was pinned to z = 0. Coordinates were rotated to MSS FSD and uniformly scaled by 20, producing an approximately 20 m hull with about 2 m breadth and 1.25 m draft. The binary STL contains 1,798 panels. The simplified submerged volume is approximately 22.13 m3, about 0.1% below the original scaled wet surface. Simplification introduces geometric error; this is an import example, not a claim of experimental validation.

The configuration uses the mesh displacement for mass at 1025 kg/m3. CG, radii of gyration, and viscous damping are illustrative assumptions. The finite-frequency grid runs from 0.1 to 2.8 rad/s in 0.1 rad/s steps for smoother demonstration plots. The largest hull panel radius is about 0.909 m, which is below wavelength/8 at 2.8 rad/s (about 0.983 m). This heuristic check is not a substitute for a convergence study. The preparation libraries are not required to read or solve the shipped STL.

## Import convention

Unlike GDF and MGF, STL does not specify an axis convention or physical units. MSS-Capytaine requires STL coordinates in **metres and MSS FSD**. For a surface hull, put z = 0 at the design waterline, with the underwater hull at positive z. Supply the full port/starboard hull with outward-facing triangles. The importer loads both ASCII and binary STL, converts to solver coordinates internally, and clips off any above-water geometry for a surface vessel. It does not resample or automatically simplify the hull panels.

After solving, inspect the FSD output. Hydrodynamic post-processing requires MSS on the MATLAB/Octave path:

```matlab
plotVesselGeometry('mesh_stl/wigley')
addpath(genpath('/path/to/MSS'))
plotVesselHydrodynamics('mesh_stl/wigley');
```

## Velocity-decay time constants

`T_126 = [80, 8, 12]` seconds specifies total surge, sway and yaw damping targets. MATLAB displays the potential-flow T_max values and rejects infeasible requests; it does not clamp them. These are illustrative estimates, not measured Wigley damping coefficients. After damping-only edits, run `python main.py mesh_stl/wigley --update-damping` to refresh the MAT settings without solving again.
