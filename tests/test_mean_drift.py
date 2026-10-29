"""Check mean-drift normalization, CG reference, headings and quadrature."""
import unittest

import capytaine as cpt
from capytaine.post_pro.mean_drift_force import far_field_mean_drift_force
import numpy as np
import xarray as xr

from src.mean_drift import compute_mean_drift, kochin_angles, angular_difference


class MeanDriftTest(unittest.TestCase):
    def test_theta_grid_and_invalid_settings(self):
        for step in (1., .5, 2.):
            theta = kochin_angles(step)
            self.assertLess(theta[0], 0)
            self.assertGreater(theta[-1], 2*np.pi)
            self.assertEqual(theta[1], 0)
            self.assertEqual(theta[-2], 2*np.pi)
        for invalid in (True, 0, -1, 7, 11, np.nan, np.inf, 'bad'):
            with self.subTest(step=invalid), self.assertRaises(ValueError):
                kochin_angles(invalid)

    def synthetic(self):
        omega = np.array([.5, 1.])
        beta = np.array([0., np.pi/2, np.pi])
        theta = kochin_angles()
        ds = xr.Dataset(coords=dict(omega=omega, wave_direction=beta, theta=theta,
                                    radiating_dof=['Heave'], forward_speed=0.,
                                    water_depth=np.inf, rho=1025., g=9.81,
                                    wavenumber=('omega', omega**2/9.81)))
        ds['kochin_diffraction'] = (('omega','wave_direction','theta'),
                                    np.ones((2,3,len(theta)), dtype=complex))
        ds['kochin_radiation'] = (('omega','radiating_dof','theta'),
                                  np.zeros((2,1,len(theta)), dtype=complex))
        motion = xr.DataArray(np.zeros((2,3,1)), dims=('omega','wave_direction','radiating_dof'),
                              coords={k:ds.coords[k] for k in ('omega','wave_direction','radiating_dof')})
        return ds, motion

    def test_analytic_normalization_and_CG_transfer(self):
        # H_total=i: the integral terms vanish and the forward scattering
        # term is exactly -2*pi*rho*omega*[cos(beta),sin(beta)].
        ds, motion = self.synthetic()
        cg = np.array([2., -3., 0.])
        result = compute_mean_drift(ds,motion,cg).values
        scale = -2*np.pi*float(ds.rho)*ds.omega.values[:,None]
        expected_x = scale*np.cos(ds.wave_direction.values)
        expected_y = scale*np.sin(ds.wave_direction.values)
        np.testing.assert_allclose(result[:,:,0],expected_x,atol=1e-10)
        np.testing.assert_allclose(result[:,:,1],expected_y,atol=1e-10)
        np.testing.assert_allclose(result[:,:,2],-cg[0]*expected_y+cg[1]*expected_x,atol=1e-10)
        np.testing.assert_allclose(angular_difference(ds,motion,cg,
            compute_mean_drift(ds,motion,cg)),0.,atol=1e-12)

    def test_diagonal_matches_full_directional_Capytaine_result(self):
        ds, motion = self.synthetic()
        ds['kochin_diffraction'] *= np.exp(1j*ds.theta)*(.2+ds.wave_direction)
        ds['kochin_radiation'] += .1*np.cos(ds.theta)
        motion = motion.astype(complex)+.2+.3j
        full = far_field_mean_drift_force(motion,ds)
        diagonal = compute_mean_drift(ds,motion,[0,0,0])
        for j,key in enumerate(('drift_force_surge','drift_force_sway','drift_force_yaw')):
            expected = np.diagonal(full[key].values,axis1=1,axis2=2).real
            np.testing.assert_allclose(diagonal.values[:,:,j],expected,rtol=1e-12,atol=1e-10)

    def test_fixed_cylinder_symmetry_translation_and_angular_convergence(self):
        mesh = cpt.mesh_vertical_cylinder(length=1.,radius=.5,center=(0,0,-.5),
                                          resolution=(2,24,4)).immersed_part()
        headings = np.array([0., np.pi/2, np.pi])
        solver = cpt.BEMSolver()
        loads = []
        for center in ([0.,0.,-.5], [1.3,-.7,-.5]):
            shifted = mesh.translated([center[0],center[1],0.])
            body = cpt.FloatingBody(mesh=shifted, dofs=cpt.rigid_body_dofs(rotation_center=center),
                                    center_of_mass=center)
            ds = solver.fill_dataset(xr.Dataset(coords=dict(omega=[2.],wave_direction=headings,
                theta=kochin_angles(.5),radiating_dof=['Heave'])),body,hydrostatics=False,
                progress_bar=False)
            # Fixed body: radiation is solved but its motion multiplier is zero.
            motion = xr.DataArray(np.zeros((1,3,1)),dims=('omega','wave_direction','radiating_dof'),
                coords={k:ds.coords[k] for k in ('omega','wave_direction','radiating_dof')})
            fine = compute_mean_drift(ds,motion,center).values[0]
            # Keep the seam's neighbor on both sides on the coarser grid.
            coarse_ds = ds.isel(theta=np.arange(1,ds.sizes['theta']-1,2))
            # Add periodic neighbors, retaining values calculated at the same angles.
            coarse_ds = xr.concat([ds.isel(theta=[-4]).assign_coords(theta=[-np.pi/180]),
                                   coarse_ds,
                                   ds.isel(theta=[3]).assign_coords(theta=[2*np.pi+np.pi/180])],
                                  dim='theta',data_vars='minimal',coords='minimal',compat='equals')
            coarse = compute_mean_drift(coarse_ds,motion,center).values[0]
            np.testing.assert_allclose(coarse,fine,rtol=2e-3,atol=.03)
            along = fine[:,0]*np.cos(headings)+fine[:,1]*np.sin(headings)
            across = -fine[:,0]*np.sin(headings)+fine[:,1]*np.cos(headings)
            self.assertTrue(np.all(along>0))
            np.testing.assert_allclose(across,0.,atol=.03)
            np.testing.assert_allclose(fine[:,2],0.,atol=.03)
            loads.append(fine)
        # Moving the geometry and the CG together must preserve CG loads.
        np.testing.assert_allclose(loads[0],loads[1],rtol=2e-3,atol=.03)


if __name__ == '__main__':
    unittest.main()
