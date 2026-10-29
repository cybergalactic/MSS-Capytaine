"""Check a regenerated Chapter 10 mean-drift MAT file without running BEM.

Usage: python tests/check_mean_drift.py [path/to/vessel.mat]
This checks the export contract, not experimental accuracy or mesh convergence.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from plot_results import drift_data, load_vessel  # noqa: E402


def check(path):
    vessel = load_vessel(path)
    if not hasattr(vessel, 'drift'):
        raise ValueError('No vessel.drift data: regenerate this case with python main.py <case>')
    headings = np.asarray(vessel.headings).ravel()
    frequencies = np.asarray(vessel.freqs).ravel()
    omega, values = drift_data(vessel,headings,frequencies)
    for field in ('forceRAO','motionRAO'):
        if not np.array_equal(omega,np.asarray(getattr(vessel,field).item().w).ravel()):
            raise ValueError(f'Mean-drift frequencies do not match {field}.w')
    if not np.array_equal(np.asarray(vessel.velocities).ravel(),[0.]):
        raise ValueError('This export requires zero speed')
    if np.any(omega<=0) or np.any(omega>=10) or np.any(np.diff(omega)<=0):
        raise ValueError('Mean-drift frequencies must be sorted, finite and exclude limit markers')
    for k,beta in enumerate(headings):
        mirrored = np.mod(-beta,2*np.pi)
        if np.isclose(mirrored,2*np.pi): mirrored=0.
        matches = np.flatnonzero(np.isclose(headings,mirrored,rtol=0,atol=1e-10))
        if len(matches)!=1: raise ValueError('Missing mirrored mean-drift heading')
        np.testing.assert_allclose(values[:,k,:],values[:,matches[0],:]*[1,-1,-1],
                                   rtol=1e-10,atol=1e-10)
    return dict(status='passed',finite_samples=len(omega),headings=len(headings),
                finite_band_rad_s=[float(omega[0]),float(omega[-1])],
                physical_dofs=[1,2,6],reference='CG',axes='MSS FSD',
                normalization='per incident wave amplitude squared',
                coefficient_min=values.min(axis=(0,1)).tolist(),
                coefficient_max=values.max(axis=(0,1)).tolist())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result',nargs='?',type=Path,
        default=ROOT/'vessels_capytaine/testShip_csv/testShip.mat')
    args=parser.parse_args()
    print(json.dumps(check(args.result),indent=2))
