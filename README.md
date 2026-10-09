# MSS-Capytaine

**Version 1.0 (2026)**

MSS-Capytaine is a Python add-on for the [Marine Systems Simulator (MSS)](https://github.com/cybergalactic/MSS). It uses the open-source [Capytaine](https://capytaine.org/) boundary-element solver to compute 6-DOF linear potential-flow hydrodynamics and exports the results as the standard MATLAB/Octave `vessel` structure used by MSS.

The project provides an open-source hydrodynamic-data workflow for MSS users without access to the commercial ShipX or WAMIT solvers. Capytaine performs the boundary-element calculations; MSS provides the MATLAB and GNU Octave functions for analysis, model reduction, nonlinear model refinement, plotting, and time-domain simulation.

The repository includes a synthetic surface monohull and an idealized submerged case named `LAUV_marie`. LAUV stands for Light Autonomous Underwater Vehicle; OceanScan commercializes the platform family, and Marie is NTNU AUR-Lab's vehicle. The cases demonstrate the complete Capytaine-to-MSS workflow rather than representing validated vessel designs.

Author: Thor I. Fossen

## Recommended citation

If you use MSS-Capytaine in your work, please cite the software as:

> Fossen, Thor I. (2026). *MSS-Capytaine* (Version 1.0) [Computer software]. GitHub. https://github.com/cybergalactic/MSS-Capytaine

BibTeX:

```bibtex
@software{Fossen2026_MSScapytaine,
  author  = {Fossen, Thor I.},
  title   = {MSS-Capytaine},
  version = {1.0},
  year    = {2026},
  url     = {https://github.com/cybergalactic/MSS-Capytaine}
}
```

## MSS Toolbox Integration

The version 1.0 workflow uses a named vessel folder containing an offset CSV file for the hull geometry and a JSON file for the vessel particulars and calculation settings. The folder name is also the vessel name used in the Python command. For example, a vessel named `myVessel` must have this catalog structure:

```text
vessels_capytaine/
  myVessel/
    config.json
    offset_points.csv
```

1. **Create the vessel catalogue and offset CSV with AI assistance.** Under `vessels_capytaine/`, create a folder named after the vessel, for example `vessels_capytaine/myVessel/`. Use a simple folder name without spaces because the same name will be entered on the command line. Copy the existing [`offset_points.csv`](vessels_capytaine/testShip/offset_points.csv) into the new folder and use it as the required template. Upload the template, along with the vessel datasheet, drawings, and relevant photographs, to an AI assistant such as ChatGPT, and ask it to produce a new CSV with the same columns and coordinate convention. Save the result as `offset_points.csv`. AI-generated geometry is a starting point and must be checked by the user.
2. **Enter the vessel particulars in JSON.** Copy the existing [`config.json`](vessels_capytaine/testShip/config.json) into the same folder. Set `body_name` to `myVessel`, `output_filename` to `myVessel.mat`, `offset_points_csv` to `offset_points.csv`, and `output_dir` to `results`. Update the principal dimensions, mass properties, center of mass, mesh resolution, frequency range, and damping inputs, then save the file as `config.json`.
3. **Run MSS-Capytaine and check its accuracy warnings.** From the repository root, type `python main.py myVessel`. This single command automatically generates the mesh, runs the Capytaine BEM calculations, converts the results to MSS forward-starboard-down (FSD) coordinates, and exports the MSS `vessel` structure to a MATLAB/Octave MAT-file. Before symmetrizing the hydrodynamic matrices, the program checks the reciprocal coupling terms in the added-mass matrix `A` and the potential damping matrix `B`. A reciprocity warning can indicate that the mesh is too sparse for the selected frequency range. Increase `samples_per_section` or `number_of_stations` in `config.json`, rerun the calculation, and repeat until the warning is resolved or otherwise understood.
4. **Inspect the generated mesh.** In MATLAB or Octave, add MSS and the MSS-Capytaine MATLAB folder to the path, then call `plotVesselMesh('myVessel')`. The function finds and loads the generated `.mat` and mesh files automatically.

   ```matlab
   addpath(genpath('/path/to/MSS'))
   addpath('/path/to/MSS-Capytaine/matlab')
   plotVesselMesh('myVessel')
   ```

   Check the hull shape, waterline, station spacing, panel resolution, and any Capytaine warnings. If the mesh is satisfactory, continue to the hydrodynamic-results inspection. If it is not satisfactory, return to the start of the workflow, revise the offset CSV or source information, and run the calculation again.
5. **Inspect the hydrodynamic results.** Use the same catalog name to plot the force RAOs, added mass, radiation damping, and viscous damping. This function also finds and loads the generated vessel file automatically:

   ```matlab
   plotVesselHydrodynamics('myVessel');
   ```

6. **Copy the MSS vessel file.** Copy the generated `vessels_capytaine/myVessel/results/myVessel.mat` file to a matching vessel directory under `MSS/HYDRO/vessels_capytaine/`.
7. **Simulate the generated model in MSS.** In MATLAB or GNU Octave, run [`SIMhydroVessel.m`](https://github.com/cybergalactic/MSS/blob/master/CRAFT/SIMhydroVessel.m) to exercise the generated 6-DOF model with waves and visualize the baseline response.
8. **Optionally refine and validate the model in MSS.** Treat the exported Capytaine model as a linear potential-flow baseline. Use experimental tests or measured time series to identify and calibrate the viscous damping matrix. Augment the model with the velocity-dependent rigid-body and added-mass Coriolis–centripetal matrices and the applicable nonlinear loads: ITTC quadratic surge drag, cross-flow drag, and quadratic roll, pitch, and yaw damping. When required, include force and moment models for actuators such as control surfaces, thrusters, and propellers. Compare simulated and measured responses and iterate until the model is adequate for its intended use.

The data flow is:

```text
Datasheets + drawings + photographs + CSV template
                         |
                         v
          AI-assisted geometry processing
       Generate offset data; review and verify
                         |
                         v
        vessels_capytaine/myVessel/
          |-- offset_points.csv
          `-- config.json
                         |
                         v
             python main.py myVessel
                         |
                         v
       Added mass (A) or potential damping (B)
                 reciprocity warning?
          |-- yes --> Refine mesh settings in config.json
          `-- no
               |
               v
  ------------- Python -> MATLAB/GNU Octave -------------
               |
               v
          plotVesselMesh('myVessel')
               |
               v
          Is the mesh visually acceptable?
          |-- no --> Return to top; revise inputs
          `-- yes
               |
               v
          plotVesselHydrodynamics('myVessel')
                                 |
                                 v
                    Copy results/myVessel.mat to MSS
                         |
                         v
          Run SIMhydroVessel.m in MSS
                         |
                         v
         Baseline 6-DOF simulation with waves
                         |
                         v
       Optional: refine the linear model in MSS
         |-- calibrate viscous damping from
         |   experimental tests or time series
         |-- add velocity-dependent Coriolis-
         |   centripetal terms and nonlinear drag
         `-- optionally add actuator force/moment
             models for control surfaces and propulsors
                         |
                  Validate against data
                         |
                         v
        Iterate calibration and validation
          until the model is fit for purpose
                         |
                         v
           Validated nonlinear MSS model
```

### MSS model refinement and validation

The exported mass, added-mass, radiation-damping, restoring, and wave-excitation data define a linear model about the selected equilibrium condition. They are not, by themselves, a validated nonlinear maneuvering model. Experimental data can come from model-basin tests, free-decay tests, captive tests, or full-scale trials; recorded input-output time series can also be used for parameter estimation. The calibration data should cover the operating conditions in which the model will be used.

MSS [`hydroVessel.m`](https://github.com/cybergalactic/MSS/blob/master/CRAFT/hydroVessel.m) shows how to extend this baseline with the velocity-dependent rigid-body and added-mass Coriolis–centripetal matrices, ITTC-1957 quadratic surge resistance, strip-theory cross-flow drag, and quadratic roll, pitch, and yaw damping. The supporting MSS functions include `rbody`, `m2c`, `XuuITTC`, and `crossFlowDrag`. Application models can additionally represent the forces and moments produced by control surfaces, thrusters, and propellers, together with relevant actuator dynamics and limits. These nonlinear terms and their coefficients are modeling assumptions; calibrate them against measurements and validate the complete model on data that were not used for calibration.

MSS-Capytaine converts the Capytaine results to the MSS conventions before export:

- Six degrees of freedom ordered as surge, sway, heave, roll, pitch, and yaw;
- Forward-starboard-down (FSD) body axes;
- Hydrodynamic matrices referenced to the CG;
- Wave headings expressed as MSS propagation directions on the fixed 10° grid;
- Zero vessel speed for the current Capytaine calculation; and
- A full 0°–350° directional set obtained by mirroring the symmetric 0°–180° solution.

The exported `vessel` structure contains:

| Field | Contents |
| --- | --- |
| `main` | Vessel particulars, mass properties, centers, and metacentric heights |
| `MRB` | Rigid-body mass matrix |
| `A` | Zero-, finite-, and infinite-frequency added mass |
| `B` | Potential damping |
| `C` | Hydrostatic restoring matrix |
| `forceRAO` | First-order wave-excitation force RAOs |
| `motionRAO` | First-order motion RAOs computed with potential damping `B` |
| `freqs` | Coefficient frequency grid |
| `headings` | Full directional grid |
| `velocities` | Vessel-speed grid, currently `[0]` |
| `powerBased` | Inputs used by MATLAB to compute the constant power-based `Bv` |

For surface vessels, `vessel.main.CF` stores the center of flotation in the same MSS FSD body coordinates as `vessel.main.CG` and `vessel.main.CB`. Thus `x_F = CF(1) - CG(1)`. The hydrostatic restoring matrix is constructed as `G35 = G53 = -G33*x_F` and includes the corresponding `G33*x_F^2` contribution in `G55`.

Pre-generated copies are included in the MSS [`HYDRO/vessels_capytaine`](https://github.com/cybergalactic/MSS/tree/master/HYDRO/vessels_capytaine) catalog. You need Python and Capytaine to regenerate the hydrodynamic data, but not to load the included `.mat` files in MSS.

## Requirements

The Python calculation requires:

- Python 3.11 or later;
- Capytaine 3.0 or later;
- NumPy, SciPy, and Xarray; and
- Matplotlib for the optional Python plots.

The current workflow was developed and tested with Python 3.11 and Capytaine 3.0.0. A basic virtual environment can be created with:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install capytaine numpy scipy xarray matplotlib
```

On Windows, activate the environment with `.venv\Scripts\activate` instead. MSS is not required to run the Python calculation, but it is required for the MATLAB/Octave post-processing workflow.

## Quick start

Run the default surface-vessel calculation from the repository root:

```sh
python main.py
```

`main.py` reads `vessels_capytaine/testShip/config.json` and writes the generated files to `vessels_capytaine/testShip/results/`.

List or select catalog vehicles by name:

```sh
python main.py --list
python main.py testShip
python main.py LAUV_marie
```

A direct path to another configuration is also accepted:

```sh
python main.py path/to/config.json
```

Run the submerged LAUV Marie case with:

```sh
python main.py LAUV_marie
```

This writes `vessels_capytaine/LAUV_marie/results/LAUV_marie.mat`. The case uses NTNU's published 2.15 m length and 34 kg mass for Marie, the published 0.15 m OceanScan LAUV-family diameter, and a documented idealized tapered-cylinder hull. Its finite-frequency grid extends to 9.5 rad/s; 10 rad/s remains the plotting location for the separately computed infinite-frequency result.

To inspect either result using MSS, add MSS and the MSS-Capytaine MATLAB directory to the MATLAB or GNU Octave path, then pass the catalogue name to the two vessel-independent plotting functions. Each function automatically locates and loads the generated vessel data.

For `testShip`:

```matlab
addpath(genpath('/path/to/MSS'))
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselMesh('testShip')
plotVesselHydrodynamics('testShip');
```

For `LAUV_marie`:

```matlab
plotVesselMesh('LAUV_marie')
plotVesselHydrodynamics('LAUV_marie');
```

`plotVesselMesh` loads `myVessel.mat` and the matching generated panel file from the named catalog folder. `plotVesselHydrodynamics` loads the same vessel structure, calls `computeManeuveringModel` using the damping inputs exported from `config.json`, then uses `plotTF`, `plotABC`, and `plotBv` to inspect the MSS hydrodynamic data. It also calls `vesselPeriods` for surface vessels; fully submerged vehicles have no hydrostatic heave stiffness, so that calculation is skipped automatically.

## Repository layout

```text
main.py                         Run the example calculation
plot_results.py                 Plot the exported result using Python
src/                            Mesh, solver, coordinate, damping, and export code
vessels_capytaine/
  testShip/config.json          Example configuration
  testShip/offset_points.csv
                                Example hull offsets
  testShip/README.md            Geometry and modeling assumptions
  testShip/results/             Generated hydrodynamic data
  LAUV_marie/config.json        Submerged LAUV Marie-inspired configuration
  LAUV_marie/offset_points.csv  Idealized closed-body offsets
  LAUV_marie/README.md          Sources and modeling assumptions
matlab/plotVesselMesh.m         Plot the generated mesh for any catalog vessel
matlab/plotVesselHydrodynamics.m
                                Plot MSS hydrodynamic data for any vessel
```

## Plot in Python

After running `main.py`, use Matplotlib to plot the saved MATLAB result:

```sh
python plot_results.py
python plot_results.py --heading 90 --show
```

The script reads `vessels_capytaine/testShip/results/testShip.mat` and saves six PNG figures in `vessels_capytaine/testShip/results/plots/`. It does not rerun the hydrodynamics or change the result file. The coefficient plots show the infinite-frequency values as separate markers at the `10 rad/s` label; RAO plots use only finite frequencies.

To plot `LAUV_marie`, pass its result file explicitly:

```sh
python plot_results.py vessels_capytaine/LAUV_marie/results/LAUV_marie.mat
python plot_results.py vessels_capytaine/LAUV_marie/results/LAUV_marie.mat --heading 90
```

These figures are saved in `vessels_capytaine/LAUV_marie/results/plots/`. All six RAO panels use the same frequency limits. At symmetry headings such as 0° and 180°, the sway, roll, and yaw responses are identically zero for a port-starboard symmetric body; their phase is undefined and is labeled as such instead of being plotted on an arbitrary autoscaled frequency axis.

RAO phase figures use the principal interval from -180° to 180°. This maps equivalent 0° and 360° values to the same phase and removes artificial two-pi jumps. Curves are broken at genuine 180° phase reversals, and phase is hidden where the RAO magnitude is negligible, since phase is undefined at a zero response.

## Inputs

Offset CSV files contain the columns `x_m,z_m,half_breadth_m`; `x` points aft, `z` points upward, and half breadth is nonnegative. For a surface vessel, the origin is at midships on the design waterline, and each section runs from keel to `z = 0`. For a submerged body, each section is a closed bottom-to-top profile with zero half breadth at both ends, and the origin is body-fixed.

`vessels_capytaine/testShip/config.json` specifies the mesh resolution, mass, radii of gyration (or an inertia-matrix CSV), center of mass, and wave frequencies. The center of mass is also used as the rotation center, so all hydrodynamic matrices, forces, and motions are referenced to the CG. The solver always uses 19 wave directions from 0° to 180° in 10° increments; the MSS export workflow fixes these, so they are not configuration inputs. The limiting-frequency calculations assume infinite water depth, the default when `water_depth_m` is absent or `null`.

Set `submerged` to `false` for a surface vessel. The solver then generates an internal waterplane lid to suppress artifacts at irregular frequencies. Set it to `true` for a submerged vehicle such as an AUV; no lid is generated.

A submerged case also requires `submergence_depth_m`, the mean operating depth of the body-fixed origin used as the fixed equilibrium position for the free-surface solve. It is not an average of hydrodynamic coefficients over a depth range. This translation affects the added mass, radiation damping, and RAOs but is removed from the exported MSS CG, CB, and panel coordinates. For submerged bodies, `vessel.main.T` is the body height rather than a surface-vessel draft. Use `output_filename` to give each catalog case its own `.mat` filename.

`samples_per_section` controls resolution around each half section, while `number_of_stations` controls resolution along the hull. Check Capytaine's mesh-resolution warnings at the highest wave frequencies. The workflow checks the reciprocal `A_ij`/`A_ji` and `B_ij`/`B_ji` coupling terms before symmetrizing the added-mass and radiation-damping matrices as required at zero speed. It reports a warning when the local reciprocity error exceeds 1%, and the skew also exceeds 0.1% of the largest matrix norm over the solved frequency range. This avoids misleading ratios where submerged-body radiation damping is numerically close to zero. If a warning is reported, refine the mesh and rerun the calculation before accepting the exported model.

## Viscous damping correction

The `viscous_damping` object contains the vessel-type-specific inputs used by MSS `computeManeuveringModel`. They are stored unchanged in `vessel.powerBased`; Python does not construct a second viscous-damping matrix.

For a surface vessel, `kappa_126` gives dimensionless relative damping increments for surge, sway, and yaw. `delta_zeta_345` gives additional damping ratios for heave, roll, and pitch:

```json
"viscous_damping": {
  "kappa_126": [0.05, 0.05, 0.05],
  "delta_zeta_345": [0, 0.1, 0]
}
```

For a submerged vehicle, heave is unrestrained along with surge, sway, and yaw. `T_1236` therefore specifies four positive target time constants in seconds for those modes, while `delta_zeta_45` gives damping-ratio increments for the restored roll and pitch modes:

```json
"viscous_damping": {
  "T_1236": [50, 5, 5, 5],
  "delta_zeta_45": [0.2, 0.2]
}
```

The stored values document the selected damping assumptions and match the defaults in `computeManeuveringModel`. After loading the exported vessel, call `computeManeuveringModel(vessel, omega_p)`. It computes `A_eq`, `B_eq`, the single constant diagonal `powerBased.Bv`, and `D = B_eq + Bv`. For a surface vessel, `Bvii = kappa_i B_eq,ii` in DOFs 1, 2, and 6, and `Bvii = 2 delta_zeta_i sqrt(Mii Gii)` in DOFs 3, 4, and 5. For a submerged vehicle, the target damping in DOFs 1, 2, 3, and 6 is `Mii / T_i`; the viscous contribution is `Bvii = Mii / T_i - B_eq,ii`. The requested time constants must not require a negative viscous contribution. Roll and pitch use the same damping-ratio formula as the restored modes of a surface vessel.

Capytaine motion RAOs use potential-flow damping only. It does not export the top-level frequency-dependent `vessel.Bv`.

## Outputs

`vessels_capytaine/testShip/results/testShip.mat` contains `MRB`, `A`, `B`, `C`, the power-based damping inputs, force and motion RAOs, frequencies, and headings in MSS forward-starboard-down axes at the center of gravity. For a surface vessel, the hydrostatic restoring matrix `C` is confined to the heave, roll, and pitch block. For a fully submerged vehicle, heave waterplane stiffness is zero; `LAUV_marie` has restoring stiffness only in roll and pitch from its assumed vertical CG–CB separation. The exported `vessel.main` structure records the `submerged` flag and nominal `submergenceDepth` used by the free-surface solve.

The coefficient frequency grid includes zero frequency and the infinite-frequency radiation solution, labeled `10 rad/s`. The force and motion RAOs use only the positive finite frequencies below `10 rad/s`. The fixed 19 headings from 0° to 180° in 10° increments are mirrored to the full 36-heading set from 0° to 350° required by MSS.

Every value configured in `omega_rad_s` must be strictly below 10 rad/s. The solver rejects values greater than or equal to 10 because the exporter always reserves and appends 10 rad/s for the separately computed infinite-frequency
radiation result.

The results directory also contains `hydrostatics.json` and the generated hull panels in `.mat` and `.npz` formats. The hull and inertia values are synthetic demonstration data; numerical accuracy depends on mesh resolution.

## Current scope and limitations

- The calculation is restricted to zero forward speed.
- The exported RAOs are first-order; second-order wave-drift loads are not computed.
- Directional mirroring assumes a port-starboard symmetric monohull.
- Zero- and infinite-frequency radiation calculations currently require infinite water depth.
- Mesh convergence and frequency-range convergence must be checked for each new hull.

## Related projects

- [MSS](https://github.com/cybergalactic/MSS) — MATLAB and GNU Octave marine systems simulation library.
- [Capytaine](https://github.com/capytaine/capytaine) — Python boundary-element solver for linear potential-flow hydrodynamics.
