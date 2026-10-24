# Otter at 80 kg

Approximate Otter twin-hull example: the MSS 55 kg base vessel plus 25 kg additional payload. The companion [otter_55kg](../otter_55kg/README.md) case uses the same underlying physical hull at a different loading condition. These are assumed contours, not measured Otter lines or a validated wave model.

- `config.json`: total mass, FSD CG, gyration radii, hull spacing and frequencies.
- `offset_points.csv`: one demi-hull, with breadths measured from its local centerline.
- `results/`: created by the Python solve.

The JSON spacing of 0.79 m places the two hulls at y = ±0.395 m. All public coordinates use MSS FSD, relative to this case's equilibrium waterline. Its static trim is level; geometry, CG and inertia already include that attitude. The familiar gyration radii are supplemented by signed `R46_m` for roll–yaw cross-inertia. Omit it or set it to zero for diagonal inertia. The two loads share 76 finite frequencies over 0.2–9.5 rad/s, with closer sampling around 5–6 and 8 rad/s. The additional roll damping ratio is 0.2, a starting estimate guided by the MSS Otter model rather than a measured decay test. JSON CG and radii are rounded to four decimals for readability; the assumed hull has no claimed sub-millimetre physical accuracy.

Run from the repository root:

```sh
python main.py capytaine_offsets_csv/otter_80kg
```

Preview the sections before solving, then repeat to see the saved mesh:

```matlab
plotVesselGeometry('otter_80kg')
plotVesselHydrodynamics('otter_80kg')       % default weighting peak: 0.8 rad/s
plotVesselHydrodynamics('otter_80kg', 1.5)  % choose another weighting peak
```

Hydrodynamic plotting requires MSS on the MATLAB/Octave path. Check mesh convergence before relying on the coefficients.

Optional [geometry derivation, sources and reproduction script](../../../tests/otter/README.md) are kept with development files; users do not need to run the generator.

## Velocity-decay time constants

`T_126 = [5, 0.5, 1]` seconds specifies total surge, sway and yaw damping targets. MATLAB displays the potential-flow T_max values and rejects infeasible requests; it does not clamp them. The starting values come from the MSS Otter sea-trial calibration. Their transfer to the approximate wave model still needs validation. After damping-only edits, run `python main.py capytaine_offsets_csv/otter_80kg --update-damping` to refresh the MAT settings without solving again.
