# Otter geometry development notes

The user-facing cases need only their JSON and offset CSV. This folder keeps the optional generator, derivation, sources and equilibrium details together.

Run from the repository root to reproduce both cases:

```sh
python tests/otter/generate_offsets.py
```

This performs geometry fitting only, not a BEM solve. Both configurations use `mass_kg`, `R44_m`, `R55_m`, `R66_m` and signed `R46_m`. The cross term is `mass_kg × R46_m × abs(R46_m)`; omission defaults to zero. The generator rounds JSON CG coordinates and radii to four decimal places. The more detailed values below describe the fitted construction before JSON rounding; they do not imply measured hull accuracy. Rotational inertia changes by less than 0.01% under this rounding. The case offset CSV retains numerical precision to preserve the constructed mesh, without claiming physical accuracy.

The current additional roll damping ratio is 0.2, guided by MSS `otter.m` as a starting estimate. Under the default JONSWAP reference (peak 1.0 rad/s, gamma 3.3), the equivalent roll radiation damping is about 0.1551 and 0.1766 N m s/rad for 55 kg and 80 kg. The total uncoupled equivalent-model roll ratio is approximately 0.2005. It remains an assumed damping setting, not experimental calibration. An explicit weighting peak in `plotVesselHydrodynamics(case, omega_p)` selects another post-processing reference without a BEM solve.

Both JSON cases use total surge/sway/yaw targets `T_126 = [5, 0.5, 1]` s. The user identifies these as sea-trial-calibrated values in MSS `otter.m`. They are starting values for the approximate wave model; transferring them to its new inertia/reference point and to the 55 kg loading still needs validation. They replace the uncalibrated 5% radiation-damping correction. The local MATLAB helper computes additional viscous damping as M_ii/T_i minus B_eq,ii, reports T_max, and rejects negative increments. Damping-only edits can be applied with `python main.py <case> --update-damping`.

## 80 kg derivation

This case demonstrates CSV demi-hull meshing for a rigid catamaran using the **MSS default 80 kg loading: 55 kg vessel plus a 25 kg point payload at the MSS body origin**. It is an idealized Otter-like underwater hull, not a measured or manufacturer-approved reconstruction. Its wave coefficients have not been validated against an Otter experiment or reference solver. The companion [55 kg case](../../vessels_capytaine/otter_55kg_csv/README.md) uses the same physical hull at a lower waterline, with no added payload. The names identify loading conditions, not different vessel designs.

## Sources and assumptions

The [MSS Otter model](https://github.com/cybergalactic/MSS/blob/master/CRAFT/USV/models/otter.m) supplies length 2 m, overall breadth 1.08 m, demi-hull breadth 0.25 m, demi-hull centerlines at y = ±0.395 m, unloaded mass and CG, hull radii of gyration, and the default payload. Its block coefficient 0.4 and waterplane coefficient 0.75 are approximate model inputs, not measured offsets.

[Bhagavathi, Kufoalor and Hasan (2023), *Digital Twin-Driven Fault Diagnosis for Autonomous Surface Vehicles*](https://doi.org/10.1109/ACCESS.2023.3268711) describes the Otter as weighing 55 kg "in total" and illustrates its general arrangement. It does not establish a separate 25 kg payload or identify that 55 kg explicitly as an unloaded mass. The 80 kg condition here comes from MSS and the selected payload, not from that experiment. [Manufacturer photographs and brochures](https://www.maritimerobotics.com/otter) also show the twin-hull form. These sources do not establish measured underwater station coordinates. Newer brochure weights describe a different specification and are not substituted for the MSS mass in this case.

The CSV uses assumed symmetric, chined sections with a narrow bow, an aft taper and keel rocker. Three shape parameters are fitted to:

- Total displacement 80/1025 = 0.07804878049 m³.
- Total waterplane area 2 × 0.75 × 2 × 0.25 = 0.75 m².
- Longitudinal buoyancy center x = 0.1375 m, matching the loaded CG so the assumed geometry has zero static trim.

Fitting these integrals does not identify the real hull shape. In particular, the fitted aft taper and shallow stern rocker are assumptions needed by this simple shape family. Different plausible contours can have different wave coefficients. This is a reproducible starting example for replacing with measured offsets, rather than a validated Otter wave model.

The model draft is T = 80/(1025 × 2 × 0.4 × 0.25 × 2) = 0.195122 m. The wetted hull span is 0.79 + 0.25 = 1.04 m; MSS's 1.08 m is the overall vessel breadth, including the structure. The deck, thrusters, appendages, above-water structure and viscous flow are not represented by the hull mesh.

## Coordinates and loading

The CSV describes **one demi-hull about its own centerline**. Half breadth is measured from that local centerline. `demihull_spacing_m: 0.79` builds two complete hulls at y = ±0.395 m and solves them together as one six-DOF body, including hydrodynamic interaction between hulls. No hull panels bridge the gap.

All public coordinates are MSS FSD. The CSV origin is at midships on the assumed **80 kg** design waterline. MSS's maneuvering model instead keeps its body origin at the unloaded design waterline. The generator finds that 55 kg waterline by clipping the same fixed hull and matching displacement:

```
55 kg draft on this fixed hull:       0.1614833434 m
80 kg draft on this fixed hull:       0.1951219512 m
Waterline separation:                 0.0336386078 m
Combined CG from MSS body origin:     [0.1375, 0, -0.1375] m
Combined CG from loaded CSV origin:   [0.1375, 0, -0.1038613922] m
```

MSS's constant-block-coefficient draft formula gives 0.134146 m at 55 kg. It does not describe the change of immersed shape of this fixed hull. Applying that formula independently to reshape both cases would create two different physical hulls. Here, the 80 kg draft anchors the approximate geometry; the 55 kg waterline is found from its actual displaced volume. These are not measured drafts or a reproduction of the nonlinear MSS simulation's steady heave/trim. For the lighter case, the generator then solves heave and pitch equilibrium on the same fixed hull, with the CG and inertia rotated consistently. The resulting static pitch is about -0.289° (bow up); its final maximum draft is about 0.164732 m. The 55 kg figure in the table above is the untrimmed reference used to define the common physical body origin.

The JSON mass and gyration radii reconstruct the complete 6-by-6 rigid-body mass matrix at the combined CG, in FSD. The rotational block is

```
[10.95182   0        0.6875
  0        15.125    0
  0.6875    0       14.4375]  kg m²
```

It combines MSS hull inertia with the 25 kg point payload using the parallel-axis theorem. MSS hull radii already describe the complete vessel: do not add the hull-separation contribution a second time. A different mass or payload position requires new CG, inertia and equilibrium geometry.

The JSON viscous-damping fractions are illustrative surface-vessel settings; they are not fitted Otter measurements or copies of MSS maneuvering damping. The zero-speed linear wave model does not replace the nonlinear MSS Otter propulsion and maneuvering model. No independent `GM_T_m` measurement is available for this approximate case, so that optional reference is omitted.

## Preview and solve

Both Otter loading cases use the same **76 finite frequencies**, from 0.2 to 9.5 rad/s. The grid uses 0.25 rad/s spacing from 2 to 4.5 rad/s, 0.1 rad/s spacing from 4.5 to 9 rad/s, and 0.05 rad/s spacing from 7.5 to 8.5 rad/s. The original low-frequency samples and all 24 previous samples are retained; 9.25 rad/s is also included. This resolves the changes around 5–6 rad/s and the sharp features near 8 rad/s more closely in both loads.

Rerun both cases if their saved results use the earlier 24-point grid. Extra frequencies improve sampling, rather than proving that a sharp feature is physical. Check its position and magnitude with a finer hull/lid mesh before accepting it as converged. The highest frequency and hull sampling settings are unchanged by this grid update. JSON lists remain compact, eight values per line.

Preview the sections before running Python:

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('otter_80kg')
```

Run from the repository root:

```sh
python main.py otter_80kg_csv
```

Then repeat the geometry plot to see the saved mesh below the source sections. Inspect mesh resolution, the two separate hulls, flotation and the internal free-surface lid; refine sampling and compare coefficients before relying on the results. Hydrodynamic plotting works without MSS on the MATLAB/Octave path:

```matlab
plotVesselHydrodynamics('otter_80kg')
```

`generate_offsets.py` reproduces **both cases'** assumed CSVs, CGs and inertia from the documented constraints using 12 transverse samples and 48 stations:

```sh
python tests/otter/generate_offsets.py
```

Run it from the repository root if regenerating these original examples. Changing those mesh settings can change the discretized displacement slightly; changing the loading requires re-equilibrating the same source geometry.

## Choosing a loading condition

Public sources report different masses and often omit a loading breakdown:

| Source | Reported figure | Meaning and limitation |
| --- | --- | --- |
| [Bhagavathi, Kufoalor and Hasan (2023)](https://doi.org/10.1109/ACCESS.2023.3268711), Section III-B | 55 kg total | No separate payload mass or trial weighing breakdown is established. |
| [IMR field report (2022)](https://www.hi.no/en/hi/nettrapporter/toktrapport-en-2022-9), Section 2.1.2 | 65 kg | Describes Akvaplan-niva's Otter with EK80 echosounders; does not separately weigh the installed payload. |
| [SeaBee/NIVA infrastructure report (2025)](https://seabee.no/wp-content/uploads/2025/07/8120-2025-SeaBees-Drone-and-Sensor-Infrastructure.pdf), p. 29 | 85 kg with batteries; 30 kg maximum payload | Operator inventory/performance specification; not a measured 30 kg mission payload. |
| [Current manufacturer specification](https://www.maritimerobotics.com/otter) | 62 kg dry; 30 kg payload capacity | A different specification from the historical 55 kg vessel; capacity is not typical carried mass. |

The sources do not establish a universal Otter operating mass. A 25 kg payload is an illustrative MSS loading choice, not a demonstrated typical research payload. These cases retain the historical MSS 55 kg base mass and offer 55 kg and 80 kg alternatives without changing MSS defaults or claiming to reproduce either experiment. Choose the case at simulation initialization. Both exported matrices are at their respective CGs; use at a fixed body origin requires each case's CG and waterline-origin shift. The equilibrium records below specify the exact pose and point transform from the common unloaded, level body frame.

See the [JSON guide](../../How%20to%20configure%20the%20vessel%20JSON%20file.md#twin-hulls-and-catamarans).

## 55 kg derivation

This is the **MSS 55 kg base vessel with no added payload**. It is paired with [otter_80kg](../../vessels_capytaine/otter_80kg_csv/README.md), which adds a 25 kg point payload at the MSS body origin. Both use the same assumed physical demi-hull shape, not two independently resized hull designs.

The companion case README records the MSS dimensions, published sources, geometry assumptions and payload discussion. These are idealized Otter-like contours, not measured hull lines. Hasan and co-authors' 2023 paper reports 55 kg total, but does not establish the loading breakdown needed to reproduce its experiment exactly. This case uses MSS's own base-vessel CG and inertia.

## Loading and coordinates

| Quantity | Value |
| --- | --- |
| Total mass | 55 kg |
| Additional payload | 0 kg |
| Displacement | 55/1025 = 0.05365853659 m³ |
| FSD CG in the level equilibrium frame | Approximately [0.201005, 0, -0.198963] m |
| Demi-hull centerline spacing | 0.79 m |
| Maximum draft of the trimmed wetted geometry | Approximately 0.164732 m |
| Waterplane area | Approximately 0.681372 m² |
| Principal CG rotational inertias | [10.26432, 13.75, 13.75] kg m²; matrix rotated with the hull |

The shared generator lowers the waterline on the fixed 80 kg hull, then solves heave and pitch together so the final CSV mesh displaces 55 kg with its longitudinal buoyancy center below CG. It rigidly rotates the original hull and reslices its immersed surface into vertical CSV stations. It rotates CG and the inertia matrix into the same level equilibrium frame. Both CSVs use MSS FSD and their own design waterline at z = 0. Half breadth remains local to one demi-hull; JSON spacing places the two hulls at y = ±0.395 m. Total-vessel inertia is supplied at CG.

The 55 kg draft differs from MSS's constant-block-coefficient estimate of 0.134146 m. A fixed physical hull's immersed shape changes with waterline; its block coefficient need not remain 0.4 at both loads. The unloaded waterline of this assumed geometry defines the MSS body-origin reference for both cases.

Removing the payload changes the physical equilibrium trim. Relative to the unloaded, level body frame, the fitted FSD pitch is approximately -0.28851° (bow up), with a heave translation of +0.000026907 m. The untrimmed waterline used to define that body origin is 0.033638608 m below the 80 kg waterline in the fixed hull coordinates. The generator balances volume and pitch using the final resampled mesh; it does not force both loading cases to have the same attitude. Waves are computed about the respective equilibrium shapes. This is hydrostatic consistency of an assumed shape, not measured Otter trim or experimental wave validation. Check mesh convergence and calibrate the geometry before relying on the coefficients. Viscous damping is illustrative. The very ends of the trimmed hull are cut back by 1 mm to avoid numerically collapsed tip panels; the final displacement fit includes that approximation. The exact pose and mapping from the common body frame are in the equilibrium records below.

## Preview and solve

Both loads use the same **76-point finite frequency grid**, from 0.2 to 9.5 rad/s, retaining all 24 original samples. Spacing is 0.25 rad/s from 2 to 4.5 rad/s, 0.1 rad/s from 4.5 to 9 rad/s, and 0.05 rad/s from 7.5 to 8.5 rad/s to resolve the sharp features near 8 rad/s. See the companion [frequency-grid guidance](../../vessels_capytaine/otter_80kg_csv/README.md#preview-and-solve). Rerun results made with the previous grid. More frequency samples resolve features; mesh convergence must still establish their numerical reliability.

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('otter_55kg')
```

Run from the repository root:

```sh
python main.py otter_55kg_csv
```

Then repeat the geometry plot to compare source sections with saved panels. Use `plotVesselHydrodynamics('otter_55kg')` with the MSS-Capytaine MATLAB directory on the path. The output is `otter_55kg.mat`; the 80 kg case writes a separate file.

To reproduce both original geometries and loading properties:

```sh
python tests/otter/generate_offsets.py
```

See the provenance records below. Select a loading case at simulation initialization. Forces/matrices are exported at each load's CG; a fixed MSS body-origin model must account for that reference change, rather than changing only the scalar mass.

## otter_55kg equilibrium record

```text
Equilibrium frame for the assumed Otter geometry

Total mass: 55 kg
FSD rotation about y: -0.288510208457 degrees
Pre-rotation z offset: 0 m
Post-rotation heave: 2.69071786816e-05 m (positive down)

Point mapping from the common, unloaded level MSS body frame:
  p_equilibrium = R_y(pitch) * (p_body + [0,0,z_offset]) + [0,0,heave]
  R_y = [cos(pitch),0,sin(pitch); 0,1,0; -sin(pitch),0,cos(pitch)]

The CSV, configured CG and inertia already use the equilibrium frame.
Do not apply this rotation or translation again in the Python solve.
This records a fitted equilibrium of assumed contours, not measured Otter trim.
```

## otter_80kg equilibrium record

```text
Equilibrium frame for the assumed Otter geometry

Total mass: 80 kg
FSD rotation about y: 0 degrees
Pre-rotation z offset: 0.0336386078497 m
Post-rotation heave: 0 m (positive down)

Point mapping from the common, unloaded level MSS body frame:
  p_equilibrium = R_y(pitch) * (p_body + [0,0,z_offset]) + [0,0,heave]
  R_y = [cos(pitch),0,sin(pitch); 0,1,0; -sin(pitch),0,cos(pitch)]

The CSV, configured CG and inertia already use the equilibrium frame.
Do not apply this rotation or translation again in the Python solve.
This records a fitted equilibrium of assumed contours, not measured Otter trim.
```

## Provenance

```text
Otter example provenance

These CSV offsets are an original idealized construction for MSS-Capytaine.
They are not downloaded CAD, measured hull lines, or digitized manufacturer
offsets. No paper figures or manufacturer images are redistributed.

Loading and principal model parameters:
  MSS CRAFT/USV/models/otter.m, inspected 23 October 2026
  https://github.com/cybergalactic/MSS/blob/master/CRAFT/USV/models/otter.m
  55 kg unloaded vessel plus 25 kg point payload at the MSS body origin.

General arrangement and historical vessel description:
  R. Bhagavathi, D. K. M. Kufoalor and A. Hasan (2023).
  Digital Twin-Driven Fault Diagnosis for Autonomous Surface Vehicles.
  IEEE Access 11, 41096-41104.
  https://doi.org/10.1109/ACCESS.2023.3268711

Visual context; current specifications are not used to replace MSS loading:
  https://www.maritimerobotics.com/otter
  https://maritimerobotics.ams3.digitaloceanspaces.com/statamic/files/mr-otter-usv-brochure-a3-digital.pdf

Reproduction:
  generate_offsets.py constructs assumed chined demi-hull sections, fitting
  total displaced volume, waterplane area and longitudinal buoyancy center.
  The local demi-hull is duplicated at 0.79 m centerline spacing by the mesher.
  The 55 kg case is cut from the same physical hull and fitted to its static
  heave/trim equilibrium; its CG and inertia are rotated consistently.
  The equilibrium records above document the reference transform for each loading condition.
  The case READMEs record assumptions, coordinate origins and limitations.

Status: approximate demonstration geometry; no experimental wave validation.
```
