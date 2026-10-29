# LAUV_marie

LAUV means Light Autonomous Underwater Vehicle. The LAUV platform family was developed by the University of Porto's LSTS and is commercialized by OceanScan. Marie is NTNU AUR-Lab's particular OceanScan LAUV. This `LAUV_marie` case is an idealized submerged hydrodynamic model inspired by [NTNU's LAUV Marie](https://www.ntnu.edu/aur-lab/lauv-marie) and the [OceanScan LAUV family](https://lsts.fe.up.pt/systems/53).

NTNU publishes a length of 2.15 m and a mass in air of 34 kg for Marie. The 0.15 m diameter is the published diameter of the OceanScan LAUV family.

The exact LAUV Marie hull offsets, mass distribution, centers of gravity and buoyancy, and appendage geometry are not public. This case therefore uses a generic axisymmetric body: a 1.45 m cylindrical midbody with equal 0.35 m ellipsoidal tapers. Fins, antennas, sonar heads, and the propeller are omitted. The transverse radius of gyration is based on a uniform cylinder, and the center of gravity is assumed to be 0.01 m below the center of buoyancy.

MSS FSD is always used. JSON centers, offset CSV coordinates, exported panels, and plots use x forward, y starboard, z down. `center_of_mass_m: [0, 0, 0.01]` places CG 1 cm below the origin. The body-fixed offset geometry is centered on the origin. The recommended Capytaine reference position for this example is a center depth of 1.0 m. The 0.15 m diameter hull is then fully submerged, with its upper surface about 0.925 m below the free surface. This relatively shallow reference is useful when one generated model will be reused over a larger operating-depth range: MSS can attenuate its wave-excitation force for deeper operation without having to amplify a model generated at 5 m.

`submergence_depth_m` always defines one fixed equilibrium depth, not a depth-averaged hydrodynamic model. Setting it to another value makes MSS-Capytaine recompute the added mass, radiation damping, and RAOs at that depth. Values between 1 m and 20 m therefore represent progressively deeper models with weaker wave excitation and less free-surface influence. A value such as 100 m is also accepted, but most first-order excitation from ordinary wave frequencies will be extremely small; very long, low-frequency waves penetrate more deeply.

When this 1 m model is used by the MSS depth-attenuation model, the generated force RAOs are unchanged for `d <= 1 m`. For `d > 1 m`, each frequency component is multiplied by

```text
exp[-(omega^2/g) (d - 1 m)].
```

This is a one-sided approximation: wave loads are never amplified above the reference position. It changes the wave-excitation force only and does not recompute added mass or radiation damping as the vehicle changes depth. Users operating mainly at another depth can generate a dedicated Capytaine model at that representative depth instead. The equilibrium translation is not included in the MSS body-fixed CG, CB, or panel coordinates.

The finite-frequency grid runs from 0.1 to 9.5 rad/s. It resolves the low-frequency free-surface interaction and becomes coarser at higher frequencies as the coefficients converge toward the separately calculated infinite-frequency solution. Zero-frequency radiation coefficients are calculated separately, so near-zero finite points are unnecessary. The infinity result is stored at the conventional 10 rad/s plotting position and is not a finite-frequency calculation.

## Viscous damping inputs

Because the fully submerged vehicle has no hydrostatic heave stiffness, surge, sway, heave, and yaw are treated as aperiodic modes by MSS. The configuration specifies a positive target time constant for each of these four modes and damping-ratio increments for the restored roll and pitch modes:

```json
"viscous_damping": {
  "T_1236": [50.0, 5.0, 5.0, 5.0],
  "delta_zeta_45": [0.2, 0.2]
}
```

`T_1236` is ordered surge, sway, heave, yaw, in **seconds**. The four entries set uncoupled velocity-decay targets: after 50 s in surge or 5 s in the other three motions, about 37% of an initial velocity remains. They do not describe position-return times. Heave needs a time constant because this freely submerged body has no hydrostatic heave stiffness. `delta_zeta_45 = [0.2, 0.2]` adds 20% of critical damping in roll and pitch. See the [damping equations and examples](../../How%20to%20configure%20the%20vessel%20JSON%20file.md#viscous-damping-equations-and-examples), including the constraint that a requested time constant cannot require negative viscous damping.

The exporter stores these values as `vessel.powerBased.T_1236` and `vessel.powerBased.delta_zeta_45`. MSS `computeManeuveringModel` uses them to construct the constant viscous damping matrix. For each aperiodic mode, smaller target time constants imply greater damping.

Run this case from the repository root:

    python main.py LAUV_marie

The generated MSS structure is written to vessels_capytaine/LAUV_marie_csv/LAUV_marie.mat.

Plot the generated coefficients and RAOs from the repository root with:

    python plot_results.py vessels_capytaine/LAUV_marie_csv/LAUV_marie.mat

Select another wave heading, for example 90 degrees, with:

    python plot_results.py vessels_capytaine/LAUV_marie_csv/LAUV_marie.mat --heading 90

The figures are saved in vessels_capytaine/LAUV_marie_csv/results/plots/.

With MSS on the MATLAB or GNU Octave path, load the generated vessel and use the vessel-independent functions to inspect its mesh and hydrodynamic data:

```matlab
addpath('/path/to/MSS-Capytaine/matlab')
plotVesselGeometry('LAUV_marie')
plotVesselHydrodynamics('LAUV_marie');
```

The hydrodynamic plotting function detects that LAUV Marie is submerged and therefore skips the surface-vessel heave-period calculation automatically.
