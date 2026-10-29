"""Compare saved tanker_gdf results with MSS WAMIT without a solve.

Shift the WAMIT reference to the same CG in MSS FSD and compare added mass,
radiation damping and complex force RAOs over overlapping finite frequencies.
Write plots and error reports to vessels_capytaine/tanker_gdf/validation.
The default reference is in the sibling MSS checkout; use --wamit to override.
"""
from pathlib import Path
import argparse
import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
DOFS = ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']


def complex_force(vessel):
    values = []
    for amp, phase in zip(vessel['forceRAO']['amp'], vessel['forceRAO']['phase']):
        amp, phase = np.asarray(amp), np.asarray(phase)
        # Legacy WAMIT cells allocate extra unused speed entries; use speed zero.
        if amp.ndim == 3:
            amp, phase = amp[:, :, 0], phase[:, :, 0]
        values.append(amp * np.exp(1j * phase))
    return np.stack(values, axis=-1)


def interpolate(x, xp, values):
    shape = values.shape
    flat = values.reshape(shape[0], -1)
    return np.column_stack([np.interp(x, xp, col) for col in flat.T]).reshape(
        (len(x),) + shape[1:]
    )


def error_percent(actual, reference):
    denominator = np.linalg.norm(reference)
    return float(100 * np.linalg.norm(actual - reference) / denominator)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wamit', type=Path, default=ROOT.parent / 'MSS/HYDRO/vessels_wamit/tanker.mat',
                        help='reference tanker.mat (default: sibling MSS/HYDRO/vessels_wamit/tanker.mat)')
    parser.add_argument('--max-omega', type=float, default=1.8003,
                        help='upper finite-frequency comparison bound (rad/s)')
    args = parser.parse_args()
    if not np.isfinite(args.max_omega) or args.max_omega <= 0:
        parser.error('--max-omega must be positive and finite')
    source = ROOT / 'vessels_capytaine/tanker_gdf/tanker.mat'
    output = source.parent / 'validation'
    output.mkdir(exist_ok=True)
    c = loadmat(source, simplify_cells=True)['vessel']
    w = loadmat(args.wamit, simplify_cells=True)['vessel']
    np.testing.assert_allclose(c['headings'], w['headings'])
    assert c['hydrodynamic_reference'] == 'CG'
    assert float(np.asarray(w['velocities']).flat[0]) == 0
    # MSS legacy CG z is keel height. WAMIT's origin is the design waterline.
    r = np.array(w['main']['CG'], dtype=float)
    r[2] = w['main']['T'] - r[2]
    np.testing.assert_allclose(r, c['main']['CG'])
    x, y, z = r
    skew = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    # nu_origin = H nu_CG; conservation of work gives force_CG = H.T force_origin.
    H = np.eye(6)
    H[:3, 3:] = skew
    shifted = {k: np.einsum('ji,jkn,kl->iln', H, w[k], H) for k in ('A', 'B', 'C')}
    cw, ww = c['freqs'][1:-1], w['freqs'][1:-1]
    # Avoid endpoint extrapolation. Keep the 1.8003 point as the nominal 1.8 limit.
    indices = np.flatnonzero((cw >= ww[0]) & (cw <= ww[-1]) & (cw <= args.max_omega))
    if len(indices) < 2:
        parser.error('At least two overlapping frequency samples are required')
    omega = cw[indices]
    cf = complex_force(c)
    wf = np.einsum('ji,nhj->nhi', H, complex_force(w))
    wforce = interpolate(cw, w['forceRAO']['w'], wf)
    reference = {k: interpolate(cw, ww, a[:, :, 1:-1].transpose(2, 0, 1)).transpose(1, 2, 0)
                 for k, a in shifted.items() if k != 'C'}
    rows = []
    for i, name in enumerate(DOFS):
        rows.append(dict(DOF=name,
            A_error_percent=error_percent(c['A'][i, i, 1:-1][indices], reference['A'][i, i, indices]),
            B_error_percent=error_percent(c['B'][i, i, 1:-1][indices], reference['B'][i, i, indices]),
            force_complex_error_percent=error_percent(cf[indices, :, i], wforce[indices, :, i]),
            force_amplitude_error_percent=error_percent(abs(cf[indices, :, i]), abs(wforce[indices, :, i]))))
    with (output / 'errors.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    plt.rcParams.update({'font.size': 10})
    for k, title, filename in [('A', 'Added mass', 'added_mass.png'), ('B', 'Radiation damping', 'radiation_damping.png')]:
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
        for i, ax in enumerate(axes.flat):
            ax.plot(omega, reference[k][i, i, indices], color='black', label='WAMIT, shifted to CG')
            ax.plot(omega, c[k][i, i, 1:-1][indices], 'o--', ms=3, color='#1976b9', label='Capytaine at CG')
            ax.set(title=f'{DOFS[i]}: {k}{i+1}{i+1}', xlabel='Angular frequency (rad/s)')
            ax.set_ylabel(('kg' if i < 3 else 'kg m²') + ('/s' if k == 'B' else ''))
            ax.grid(alpha=.25)
        axes.flat[0].legend(fontsize=8)
        fig.suptitle(f'Tanker GDF: {title}, MSS FSD at CG')
        fig.savefig(output / filename, dpi=160); plt.close(fig)
    for heading in (0, 90, 180):
        j = int(np.argmin(abs(np.rad2deg(c['headings']) - heading)))
        fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
        for i, ax in enumerate(axes.flat):
            ax.plot(omega, abs(wforce[indices, j, i]), color='black', label='WAMIT, shifted to CG')
            ax.plot(omega, abs(cf[indices, j, i]), 'o--', ms=3, color='#1976b9', label='Capytaine at CG')
            ax.set(title=DOFS[i], xlabel='Angular frequency (rad/s)', ylabel='N/m' if i < 3 else 'N m/m')
            ax.grid(alpha=.25)
            if heading in (0, 180) and i in (1, 3, 5):
                # Avoid magnifying symmetry-zero numerical remnants.
                ax.set_ylim(0, 1)
                ax.text(.5, .5, 'Symmetry-zero response', transform=ax.transAxes, ha='center')
        axes.flat[0].legend(fontsize=8)
        fig.suptitle(f'Tanker GDF: force RAO amplitudes at {heading}°, MSS FSD at CG')
        fig.savefig(output / f'force_rao_{heading}deg.png', dpi=160); plt.close(fig)
    # All coupling coefficients, normalized by matching diagonal scales to handle
    # reference coefficients that vanish by symmetry without division by zero.
    coupling_stats = {}
    for k in ('A', 'B'):
        ref = reference[k][:, :, indices]
        difference = c[k][:, :, 1:-1][:, :, indices] - ref
        norms = np.linalg.norm(ref, axis=2)
        diagonal = np.diag(norms)
        normalized = np.linalg.norm(difference, axis=2) / np.sqrt(diagonal[:, None] * diagonal[None, :])
        coupling_stats[k] = (100 * normalized).tolist()
    hydro_rows = []
    for k in ('m', 'nabla', 'GM_T', 'GM_L'):
        cv, wv = c['main'][k], w['main'][k]
        hydro_rows.append((k, float(cv), float(wv), float(100 * (cv - wv) / wv)))
    cg_c = c['C'][:, :, 0]
    cg_w = shifted['C'][:, :, 0]
    gm_l_common = (cg_w[4, 4] - cg_w[2, 4]**2 / cg_w[2, 2]) / (w['main']['rho'] * w['main']['g'] * w['main']['nabla'])
    hydro_rows.append(('GM_L, common stiffness definition', float(c['main']['GM_L']), float(gm_l_common),
        float(100 * (c['main']['GM_L'] - gm_l_common) / gm_l_common)))
    stiffness_errors = {f'C{i+1}{j+1}': float(100 * (cg_c[i,j] - cg_w[i,j]) / cg_w[i,j])
        for i,j in ((2,2), (2,4), (3,3), (4,4))}
    md = ['# Tanker GDF versus MSS WAMIT' , '',
        f'Finite comparison band: **{omega[0]:.4f}–{omega[-1]:.4f} rad/s**, {len(omega)} samples.', '',
        'The selected GDF matches the 246 m MSS WAMIT geometry. Both datasets are in SI units and MSS FSD.',
        'WAMIT data are shifted from the waterline origin to CG [3.93, 0, -2.5] m.',
        'For r from origin to CG, H = [[I, skew(r)], [0, I]] maps CG velocities to origin velocities.',
        'Matrices are compared using Hᵀ A H and Hᵀ B H; forces use Hᵀ F. No phase offset is fitted.', '',
        'WAMIT complex forces are reconstructed as amplitude × exp(i × phase); the Capytaine export already uses the MSS time convention.',
        'Both datasets have the same 36 headings at zero speed. The first speed entry is selected from legacy WAMIT RAO cells.',
        'WAMIT values are linearly interpolated in angular frequency, using complex real/imaginary parts for forces.',
        'The first Capytaine point (0.0383 rad/s) is below WAMIT’s 0.03834016 and is excluded to avoid extrapolation.', '',
        '## Relative errors', '',
        'Each number is 100 × ||Capytaine − WAMIT||₂ / ||WAMIT||₂ over the band; it is a relative RMS/norm error, not the worst pointwise percentage.',
        'A and B use diagonal entries. Force errors pool all 36 headings for each DOF; complex errors include phase disagreement.', '',
        '| DOF | Added mass | Radiation damping | Complex force RAO | Force amplitude |',
        '| --- | ---: | ---: | ---: | ---: |']
    for row in rows:
        md.append('| '+row['DOF']+' | '+' | '.join(f'{row[k]:.2f}%' for k in list(row)[1:])+' |')
    md += ['', '## Hydrostatic checks', '', '| Quantity | Capytaine | Legacy MSS WAMIT | Difference |', '| --- | ---: | ---: | ---: |']
    md += [f'| {k} | {cv:.6g} | {wv:.6g} | {diff:+.3f}% |' for k, cv, wv, diff in hydro_rows]
    md += ['', 'The legacy GM_L scalar differs substantially, but the hydrostatic stiffness matrices agree closely after shifting the WAMIT matrix to CG.',
        'The common definition uses (C55_CG − C35_CG²/C33) / (rho g volume) on both sides. The legacy stored scalar is retained in the table for transparency.',
        'Stiffness-entry differences after the CG shift: ' + ', '.join(f'{name}: {value:+.3f}%' for name, value in stiffness_errors.items()) + '.',
        'Treat legacy GM_L as a separate stored-metadata/derivation issue; it is not evidence of a matching discrepancy in the BEM radiation coefficients.', '',

        '## Limits of this comparison', '',
        '- These are comparisons between two numerical solutions of the same hull, not experimental validation.',
        '- Frequencies above the selected limit are excluded from this summary. The earlier full-range run flagged mesh resolution from 1.9574 rad/s and radiation reciprocity at 5.3702 rad/s; the current configuration excludes that band.',
        '- The infinite-frequency marker at 10 rad/s is a separate limiting solution and is not interpolated as a finite frequency.',
        '- Motion RAOs are not summarized: the saved WAMIT MRB has origin coupling, and its origin inertias do not yield exactly the same CG inertia as the Capytaine diagonal MRB after transformation. The WAMIT run also specifies additional model terms in its FRC file. These must be reconciled before comparing motion RAOs.',
        '- Near-zero force components at symmetry headings are displayed as zero-scale responses rather than assigned misleading relative errors.', '',
        '## Plots', '',
        '![Added mass](added_mass.png)', '', '![Radiation damping](radiation_damping.png)', '']
    for heading in (0, 90, 180):
        md += [f'![Force RAOs at {heading} degrees](force_rao_{heading}deg.png)', '']
    md += ['Inspect the [Python comparison test](../../../tests/compare_tanker_wamit.py).', '',
        'Reproduce from the repository root with `python tests/compare_tanker_wamit.py`.',
        'Use `--wamit path/to/tanker.mat` if MSS is not in the sibling folder; `--max-omega` changes the finite comparison limit.', '']
    (output / 'README.md').write_text('\n'.join(md))
    diagnostics = {'band_rad_s': [float(omega[0]), float(omega[-1])], 'samples': len(omega),
        'reference_shift_origin_to_CG_FSD_m': r.tolist(), 'H_origin_from_CG': H.tolist(),
        'errors': rows, 'stiffness_entry_error_percent': stiffness_errors, 'GM_L_wamit_common_definition_m': float(gm_l_common),
        'coupling_error_percent_diagonal_normalized': coupling_stats,
        'infinite_added_mass_diagonal_error_percent':
            (100 * (np.diag(c['A'][:, :, -1]) - np.diag(shifted['A'][:, :, -1])) / np.diag(shifted['A'][:, :, -1])).tolist()}
    (output / 'comparison.json').write_text(json.dumps(diagnostics, indent=2)+'\n')
    print('\n'.join(md[:20])); print(f'\nReport: {output / "README.md"}')


if __name__ == '__main__':
    main()
