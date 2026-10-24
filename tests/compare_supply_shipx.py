"""Compare saved shipx_mgf/supply results with MSS ShipX/VERES; no BEM solve.

Run from the repository root: python tests/compare_supply_shipx.py
Use --shipx for a reference outside the sibling MSS checkout. Reports and plots
are saved under vessels_capytaine/shipx_mgf/supply/validation/ by default.
Frequency samples must match ShipX directly; no interpolation is performed.
"""
from pathlib import Path
import argparse
import csv
import json

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / 'vessels_capytaine/shipx_mgf/supply'
DOFS = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')
COMPARISON_DOFS = (1, 2, 3, 4, 5)  # Strip theory does not provide a surge benchmark.


def norm_error(actual, reference):
    """Relative norm error in percent; undefined for a zero reference."""
    scale = np.linalg.norm(reference)
    return float(100 * np.linalg.norm(actual-reference) / scale) if scale else None


def zero_speed(vessel):
    indices = np.flatnonzero(np.isclose(np.atleast_1d(vessel['velocities']), 0, atol=1e-12))
    if len(indices) != 1:
        raise ValueError('The vessel must contain exactly one zero-speed dataset')
    return int(indices[0])


def coefficient(vessel, key, speed):
    values = np.asarray(vessel[key], dtype=float)
    if values.ndim == 4:
        values = values[:, :, :, speed]
    elif speed != 0:
        raise ValueError('Missing speed dimension in coefficient data')
    if values.shape != (6, 6, len(vessel['freqs'])) or not np.isfinite(values).all():
        raise ValueError(f'Invalid vessel.{key} dimensions or values')
    return values


def complex_force(vessel, speed):
    values = []
    rao = vessel['forceRAO']
    for amp, phase in zip(rao['amp'], rao['phase']):
        amp, phase = np.asarray(amp), np.asarray(phase)
        if amp.ndim == 3:
            amp, phase = amp[:, :, speed], phase[:, :, speed]
        elif speed != 0:
            raise ValueError('Missing speed dimension in force RAOs')
        values.append(amp * np.exp(1j*phase))
    result = np.stack(values, axis=-1)
    if result.shape != (len(rao['w']), len(vessel['headings']), 6) or not np.isfinite(result).all():
        raise ValueError('Invalid force RAO dimensions or values')
    return result


def reference_to_cg(vessel):
    """Derive the reference-to-CG displacement from rigid-body mass coupling.

    At a point P, MRB[:3,3:] = -m*skew(r_P_to_CG). This avoids trusting a
    legacy CG label when the run actually uses the waterline above CG.
    nu_P = H*nu_CG, so M_CG = H.T*M_P*H and force_CG = H.T*force_P.
    """
    mass = float(vessel['main']['m'])
    matrix = np.asarray(vessel['MRB'], dtype=float)
    if mass <= 0 or matrix.shape != (6, 6) or not np.isfinite(matrix).all():
        raise ValueError('Invalid reference rigid-body mass matrix')
    np.testing.assert_allclose(matrix[:3, :3], mass*np.eye(3), atol=mass*1e-7)
    r = np.array([matrix[1, 5], matrix[2, 3], matrix[0, 4]]) / mass
    x, y, z = r
    skew = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    np.testing.assert_allclose(matrix[:3, 3:], -mass*skew, atol=mass*1e-6)
    np.testing.assert_allclose(matrix[3:, :3], (-mass*skew).T, atol=mass*1e-6)
    transform = np.eye(6)
    transform[:3, 3:] = skew
    return r, transform


def matching_samples(omega, source_omega, values):
    """Select existing reference samples; never interpolate or extrapolate."""
    source_omega = np.asarray(source_omega, dtype=float)
    if not np.isfinite(source_omega).all() or np.any(np.diff(source_omega) <= 0):
        raise ValueError('Reference frequencies must be finite and increasing')
    matches = np.isclose(np.asarray(omega)[:, None], source_omega[None, :], rtol=0, atol=1e-12)
    if np.any(matches.sum(axis=1) != 1):
        raise ValueError('Frequency grids do not match; rerun python main.py shipx_mgf/supply with the current JSON. Interpolation is disabled.')
    return np.asarray(values)[np.argmax(matches, axis=1)]


def compare(c, s, max_omega):
    """Align zero-speed reference frames and grids, then compute diagnostics."""
    if not np.isfinite(max_omega) or max_omega <= 0:
        raise ValueError('max_omega must be positive and finite')
    if c.get('hydrodynamic_reference') != 'CG' or c.get('hydrodynamic_axes') != 'MSS FSD':
        raise ValueError('Expected Capytaine export at CG in MSS FSD')
    if s.get('hydrodynamic_axes') not in ('FSD', 'MSS FSD'):
        raise ValueError('Expected an MSS ShipX reference in FSD axes')
    np.testing.assert_allclose(c['headings'], s['headings'], atol=1e-10)
    for key in ('CG', 'm', 'Lpp', 'B', 'T', 'rho', 'g'):
        np.testing.assert_allclose(c['main'][key], s['main'][key], rtol=1e-5, atol=1e-5)
    cs, ss = zero_speed(c), zero_speed(s)
    r, transform = reference_to_cg(s)
    # Capytaine endpoints are zero and the separately computed infinity marker.
    cw = np.asarray(c['freqs'])[1:-1]
    np.testing.assert_allclose(cw, c['forceRAO']['w'])
    if c['freqs'][0] != 0 or c['freqs'][-1] != 10 or np.any(np.diff(cw) <= 0):
        raise ValueError('Unexpected Capytaine coefficient grid')
    sw = np.asarray(s['freqs'])
    fw = np.asarray(s['forceRAO']['w'])
    lower, upper = max(sw[0], fw[0]), min(sw[-1], fw[-1], max_omega)
    selected = np.flatnonzero((cw >= lower) & (cw <= upper))
    if len(selected) < 2:
        raise ValueError('At least two overlapping finite frequencies are required')
    omega = cw[selected]
    actual = {k: coefficient(c, k, cs)[:, :, 1:-1][:, :, selected] for k in ('A', 'B', 'C')}
    reference = {}
    for k in actual:
        shifted = np.einsum('ji,jkn,kl->iln', transform, coefficient(s, k, ss), transform)
        reference[k] = matching_samples(omega, sw, shifted.transpose(2, 0, 1)).transpose(1, 2, 0)
    cf = complex_force(c, cs)[selected]
    sf = matching_samples(omega, fw, np.einsum('ji,nhj->nhi', transform, complex_force(s, ss)))
    rows = []
    for i in COMPARISON_DOFS:
        rows.append(dict(DOF=DOFS[i],
            A_error_percent=norm_error(actual['A'][i, i], reference['A'][i, i]),
            B_error_percent=norm_error(actual['B'][i, i], reference['B'][i, i]),
            force_complex_error_percent=norm_error(cf[:, :, i], sf[:, :, i]),
            force_amplitude_error_percent=norm_error(abs(cf[:, :, i]), abs(sf[:, :, i]))))
    hydro = []
    for key in ('m', 'nabla', 'GM_T', 'GM_L'):
        cv, sv = float(c['main'][key]), float(s['main'][key])
        hydro.append(dict(quantity=key, capytaine=cv, shipx=sv, difference_percent=100*(cv-sv)/sv))
    common_gm = []
    for matrices, vessel in ((actual, c), (reference, s)):
        matrix = matrices['C'][:, :, 0]
        main = vessel['main']
        common_gm.append(float((matrix[4,4]-matrix[2,4]**2/matrix[2,2]) /
                               (main['rho']*main['g']*main['nabla'])))
    hydro.append(dict(quantity='GM_L (common stiffness definition)', capytaine=common_gm[0],
                      shipx=common_gm[1], difference_percent=100*(common_gm[0]-common_gm[1])/common_gm[1]))
    stiffness = {f'C{i+1}{j+1}': float(100*(actual['C'][i,j,0]-reference['C'][i,j,0])/reference['C'][i,j,0])
                 for i,j in ((2,2), (2,4), (3,3), (4,4))}
    diagnostics = dict(band_rad_s=[float(omega[0]), float(omega[-1])], samples=len(omega),
        headings=len(c['headings']), frequency_alignment='direct matching samples; no interpolation',
        CG_FSD_m=np.asarray(c['main']['CG']).tolist(),
        reference_label=s.get('hydrodynamic_reference'), reference_shift_to_CG_FSD_m=r.tolist(),
        reference_point_FSD_m=(np.asarray(s['main']['CG'])-r).tolist(),
        H_reference_from_CG=transform.tolist(),
        shifted_MRB_error_percent=norm_error(c['MRB'], transform.T@s['MRB']@transform),
        excluded_DOFs={'Surge':'Ordinary strip theory; MSS estimates A11 and B11 from sway. No surge validation is performed.'},
        errors=rows, hydrostatics=hydro, stiffness_entry_error_percent=stiffness,
        excluded_capytaine_frequencies_rad_s=cw[(cw<lower)|(cw>upper)].tolist())
    return diagnostics, omega, actual, reference, cf, sf


def save_plots(output, omega, actual, reference, cf, sf, headings):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def exclude_surge(ax):
        ax.axis('off')
        ax.set_title('Surge excluded')
        ax.text(.5, .5, 'Strip theory provides no\nmeaningful surge benchmark',
                transform=ax.transAxes, ha='center', va='center')
    for key, label, unit in (('A','Added mass','kg'), ('B','Radiation damping','kg/s')):
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
        exclude_surge(axes.flat[0])
        for i in COMPARISON_DOFS:
            ax = axes.flat[i]
            ax.plot(omega, reference[key][i,i], 'k-', label='ShipX/VERES, shifted to CG')
            ax.plot(omega, actual[key][i,i], 'o--', ms=3, label='Capytaine at CG')
            ax.set(title=f'{DOFS[i]}: {key}{i+1}{i+1}', xlabel='Angular frequency (rad/s)',
                   ylabel=unit if i<3 else ('kg m²' if key=='A' else 'kg m²/s'))
            ax.grid(alpha=.25)
        axes.flat[1].legend(fontsize=8)
        fig.suptitle(f'Supply MGF: {label}, zero speed, MSS FSD at CG')
        fig.savefig(output/('added_mass.png' if key=='A' else 'radiation_damping.png'), dpi=160)
        plt.close(fig)
    for heading in (0, 90, 180):
        j = int(np.flatnonzero(np.isclose(np.rad2deg(headings), heading))[0])
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
        exclude_surge(axes.flat[0])
        for i in COMPARISON_DOFS:
            ax = axes.flat[i]
            ax.plot(omega, abs(sf[:,j,i]), 'k-', label='ShipX/VERES, shifted to CG')
            ax.plot(omega, abs(cf[:,j,i]), 'o--', ms=3, label='Capytaine at CG')
            ax.set(title=DOFS[i], xlabel='Angular frequency (rad/s)', ylabel='N/m' if i<3 else 'N m/m')
            ax.grid(alpha=.25)
            if heading in (0, 180) and i in (1, 3, 5):
                ax.set_ylim(0, 1)
                ax.text(.5, .5, 'Symmetry-zero response', transform=ax.transAxes, ha='center')
        axes.flat[1].legend(fontsize=8)
        fig.suptitle(f'Supply MGF: force RAO amplitudes at {heading}°, MSS FSD at CG')
        fig.savefig(output/f'force_rao_{heading}deg.png', dpi=160)
        plt.close(fig)


def write_report(output, diagnostics, source, shipx):
    d = diagnostics
    with (output/'errors.csv').open('w', newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=d['errors'][0].keys())
        writer.writeheader(); writer.writerows(d['errors'])
    (output/'comparison.json').write_text(json.dumps(d, indent=2, allow_nan=False)+'\n')
    def percent(value):
        return 'n/a' if value is None else f'{value:.2f}%'
    r = d['reference_shift_to_CG_FSD_m']
    lines = ['# Supply MGF versus MSS ShipX/VERES', '',
        f'Finite comparison band: **{d["band_rad_s"][0]:g}–{d["band_rad_s"][1]:g} rad/s**, {d["samples"]} samples, {d["headings"]} headings at zero speed.', '',
        '## Alignment', '',
        'The two datasets use the same supply-vessel MGF, design draft, mass and CG. Both are compared in SI units and MSS FSD.',
        f'Capytaine input: `{source.name}` in this case’s `results/` folder. ShipX reference: `{shipx.name}` from the MSS ShipX supply case.',
        'The shipped ShipX `input.out` selects ordinary strip theory (IMETH=1), a motion reference vertically at the waterline (IMTNCOORD=1), viscous damping and roll tanks.',
        f'The saved MSS reference label is `{d["reference_label"]}`, but its rigid-body translation/rotation coupling implies a displacement to the actual CG of [{r[0]:.6g}, {r[1]:.6g}, {r[2]:.6g}] m.',
        'For the archived supply reference, the point is longitudinally at CG and vertically at the waterline, as recorded in `input.re7` and `input.re8`. The diagnostic JSON records the point inferred from the actual saved mass matrix, so a future export already shifted to CG is handled without a second shift.',
        'For r from that point to CG, H = [[I, skew(r)], [0, I]]. Use Hᵀ A H, Hᵀ B H and Hᵀ C H for matrices, and Hᵀ F for excitation forces.',
        f'The shifted rigid-body mass matrix agrees with the configured Capytaine CG inertia to {d["shifted_MRB_error_percent"]:.6g}% in relative norm.', '',
        'The script derives the displacement from mass coupling; it does not assume that an older reference-point label is correct. It rejects mismatched loading conditions or axes.',
        'Select the zero-speed entry, including only its first populated speed slot in the legacy RAO cells. The same 36 MSS propagation headings are used without reversing them again.',
        '**No frequency interpolation is used.** The supply JSON uses the exact ShipX frequency samples up to 1.8 rad/s. Reference coefficients and complex forces are selected directly at those samples; a mismatched grid produces an error.',
        'Zero and the 10 rad/s infinity marker are excluded; ShipX has no matching limiting solutions. No phase correction is fitted.',
        'The archived ShipX run specifies a 5 m regular-wave amplitude. This does not change linear A, B or excitation transfer functions normalized per metre of wave amplitude. It can affect equivalent viscous roll damping in motion calculations, which are excluded here.', '',
        '## Relative errors', '',
        '**Surge is excluded from all error tables and plotted comparisons.** Ordinary strip theory does not provide a meaningful surge benchmark; MSS estimates A11 and B11 from sway. Comparisons cover sway, heave, roll, pitch and yaw only.', '',
        'Each percentage is 100 × ||Capytaine − ShipX||₂ / ||ShipX||₂ over this band, not a maximum pointwise error. Matrix rows use diagonal entries; force errors pool all 36 headings for each DOF.',
        'These are unweighted norms over the nonuniform native frequency samples, not frequency-integrated errors. Percentages from a different grid are not directly comparable because the sample weighting changes.',
        'Complex forces are reconstructed as amplitude × exp(i × phase) using the stored MSS phases. Their errors include phase differences; no spatial wave-phase origin adjustment is assumed or fitted. Amplitude errors are shown separately.', '',
        '| DOF | Added mass | Radiation damping | Complex force RAO | Force amplitude |',
        '| --- | ---: | ---: | ---: | ---: |']
    for row in d['errors']:
        lines.append('| '+row['DOF']+' | '+' | '.join(percent(row[k]) for k in list(row)[1:])+' |')
    lines += ['', '## Hydrostatic checks', '',
        '| Quantity | Capytaine | ShipX | Difference |', '| --- | ---: | ---: | ---: |']
    for row in d['hydrostatics']:
        lines.append(f'| {row["quantity"]} | {row["capytaine"]:.6g} | {row["shipx"]:.6g} | {row["difference_percent"]:+.3f}% |')
    lines += ['', 'Stiffness-entry differences after shifting to CG: '+', '.join(
        f'{key}: {value:+.3f}%' for key,value in d['stiffness_entry_error_percent'].items())+'.', '',
        'The stored ShipX GM_L scalar and the value derived from its restoring matrix use different conventions. The common-definition row uses (C55 − C35²/C33)/(rho g volume) for both datasets; the stored-scalar difference alone is not a measure of geometry error.', '',
        '## Interpretation and limits', '',
        '- Capytaine solves a three-dimensional panel model; this ShipX run uses ordinary strip theory. Differences are not a solver pass/fail criterion and cannot be attributed to frequency limits alone.',
        '- MSS `ABCtransform` estimates ShipX A11 and B11 by scaling sway coefficients. These estimates are excluded from validation. Full six-DOF frame transformations are still applied before selecting the five compared DOFs, so the remaining forces and moments use the correct reference point. Strip theory also omits some longitudinal coupling terms.',
        '- The MGF is resampled and closed for Capytaine. Its volume is about 0.28% smaller; geometry discretization and strip-theory assumptions can affect agreement. A short hydrostatics check doubling samples and stations (24/70 to 48/140) changed volume by about 0.04% and common-definition GM_L by about 0.04%; this does not establish convergence of A and B.',
        '- Complex force differences use the archived MSS phases. They should not be interpreted as pure hydrodynamic-model error without independently confirming the incident-wave phase origin in both calculations.',
        '- Motion RAOs are not compared: the ShipX run includes viscous damping and roll tanks, while exported Capytaine motion RAOs use potential damping only. Second-order drift is not computed by this Capytaine workflow.',
        '- No experimental data or acceptance tolerance is used. These are comparisons of saved numerical datasets. Neither solver is rerun.', '',
        '## Plots', '', '![Added mass](added_mass.png)', '', '![Radiation damping](radiation_damping.png)', '']
    for heading in (0,90,180):
        lines += [f'![Force RAOs at {heading} degrees](force_rao_{heading}deg.png)', '']
    lines += ['Inspect the [Python comparison test](../../../../tests/compare_supply_shipx.py).', '',
        'Reproduce from the repository root with `python tests/compare_supply_shipx.py`.',
        'Use `--shipx path/to/supply.mat` for another MSS checkout and `--max-omega` to reduce the band. `--output-dir` selects a separate report directory.', '']
    (output/'README.md').write_text('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--shipx',type=Path,default=ROOT.parent/'MSS/HYDRO/vessels_shipx/supply/supply.mat',
                        help='MSS ShipX supply.mat reference (default: sibling MSS checkout)')
    parser.add_argument('--max-omega',type=float,default=1.8,help='upper finite-frequency bound in rad/s (default: 1.8)')
    parser.add_argument('--output-dir',type=Path,default=CASE/'validation',help='report directory (default: supply/validation)')
    args=parser.parse_args()
    source=CASE/'results/supply.mat'
    if not source.is_file() or not args.shipx.is_file():
        parser.error('Both saved supply.mat files must exist; run the Capytaine supply case first and select the MSS ShipX reference')
    c=loadmat(source,simplify_cells=True)['vessel']; s=loadmat(args.shipx,simplify_cells=True)['vessel']
    try:
        diagnostics,omega,actual,reference,cf,sf=compare(c,s,args.max_omega)
    except (ValueError,AssertionError) as exc:
        parser.error(str(exc))
    # Check the documented motion point against the native run header when present.
    raw=args.shipx.parent/'input.re7'
    if raw.is_file():
        lines=raw.read_text().splitlines()
        if len(lines)>17 and len(lines[10].split())==4:
            point=np.array(list(map(float,lines[17].split())))
            expected=np.array([-point[3],0,float(s['main']['T'])-point[4]])
            diagnostics['native_run_point_FSD_m']=expected.tolist()
            diagnostics['native_point_matches_saved_reference']=bool(np.allclose(
                diagnostics['reference_point_FSD_m'],expected,atol=1e-5))
    diagnostics['source_files']={'capytaine':str(source.resolve()),'shipx':str(args.shipx.resolve())}
    output=args.output_dir.resolve(); output.mkdir(parents=True,exist_ok=True)
    save_plots(output,omega,actual,reference,cf,sf,np.asarray(c['headings']))
    write_report(output,diagnostics,source,args.shipx)
    print('DOF       Added mass   Damping   Complex force   Force amplitude')
    for row in diagnostics['errors']:
        print(f'{row["DOF"]:7s} '+ ' '.join(f'{row[k]:12.2f}%' if row[k] is not None else '          n/a' for k in list(row)[1:]))
    print(f'\nReport: {output/"README.md"}')


if __name__=='__main__':
    main()
