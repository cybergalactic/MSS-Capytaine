"""Checks for reference-point shifts and finite-frequency ShipX alignment."""
import unittest

import numpy as np

from compare_supply_shipx import compare, matching_samples, reference_to_cg


class ShipXComparisonTest(unittest.TestCase):
    def test_waterline_reference_despite_cg_label(self):
        mass = 1000.0
        r = np.array([2.0, 0.0, -1.3])
        x, y, z = r
        h = np.eye(6)
        h[:3, 3:] = [[0, -z, y], [z, 0, -x], [-y, x, 0]]
        at_cg = np.diag([mass]*3+[2000., 5000., 6000.])
        inverse = np.linalg.inv(h)
        vessel = dict(main=dict(m=mass), MRB=inverse.T@at_cg@inverse,
                      hydrodynamic_reference='CG')
        displacement, transform = reference_to_cg(vessel)
        np.testing.assert_allclose(displacement, r)
        np.testing.assert_allclose(transform.T@vessel['MRB']@transform, at_cg, atol=1e-10)
        # Generalized force transformation must preserve instantaneous power.
        velocity_cg = np.array([1, 2, 3, .1, .2, .3])
        force_reference = np.array([10, 20, 30, 40, 50, 60])
        self.assertAlmostEqual(force_reference@(transform@velocity_cg),
                               (transform.T@force_reference)@velocity_cg)

    def test_direct_complex_samples_reject_interpolation(self):
        np.testing.assert_array_equal(matching_samples(np.array([1., 2.]), [0., 1., 2.],
                                      np.array([1+2j, 3+4j, 5+6j])), [3+4j, 5+6j])
        for omega in (-.01, .5):
            with self.assertRaisesRegex(ValueError, 'Interpolation is disabled'):
                matching_samples(np.array([omega]), [0., 1.], np.array([1., 2.]))

    def test_known_equal_data_on_different_frequency_grids(self):
        def vessel(frequencies, capytaine):
            omega = np.asarray(frequencies)
            force_omega = omega[1:-1] if capytaine else omega
            base = np.diag(np.arange(2., 8.))
            base[2, 4] = base[4, 2] = .2
            force = (1+force_omega[:, None]) * np.ones((1, 2)) * (1+1j)
            return dict(main=dict(CG=np.array([0., 0., -1.]), m=1000., Lpp=10., B=2., T=1.,
                                  rho=1000., g=9.81, nabla=1., GM_T=1., GM_L=10.),
                        MRB=np.diag([1000.]*3+[2000., 5000., 6000.]),
                        freqs=omega, headings=np.array([0., np.pi/2]), velocities=0,
                        A=base[:, :, None]*(1+omega), B=base[:, :, None]*(1+omega),
                        C=base[:, :, None]*np.ones(len(omega)),
                        forceRAO=dict(w=force_omega, amp=[abs(force)*(i+1) for i in range(6)],
                                      phase=[np.angle(force) for _ in range(6)]),
                        hydrodynamic_reference='CG',
                        hydrodynamic_axes='MSS FSD' if capytaine else 'FSD')
        c = vessel([0., .05, .2, .5, 1.5, 1.8, 10.], True)
        s = vessel([.1, .2, .5, 1., 1.5, 1.8, 2.], False)
        # Arbitrary surge disagreement must not contribute to validation errors.
        c['A'][0, 0] *= 20
        c['B'][0, 0] *= 30
        c['forceRAO']['amp'][0] *= 40
        diagnostics, omega, *_ = compare(c, s, 1.8)
        np.testing.assert_allclose(omega, [.2, .5, 1.5, 1.8])
        self.assertEqual(diagnostics['excluded_capytaine_frequencies_rad_s'], [.05])
        self.assertEqual([row['DOF'] for row in diagnostics['errors']],
                         ['Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'])
        self.assertIn('Surge', diagnostics['excluded_DOFs'])
        for row in diagnostics['errors']:
            for key, value in row.items():
                if key != 'DOF': self.assertLess(value, 1e-10)


if __name__ == '__main__':
    unittest.main()
