"""Check symmetric BEM against the identical expanded mesh on a small hull."""
import json
from pathlib import Path
import tempfile
import unittest

import capytaine as cpt
import numpy as np
from scipy.io import loadmat

from src.hull_geometry import load_hull, bounded_lid, reflection_mesh
from src.hydrodynamics import run


class ReflectionTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.base = Path(self.directory.name)
        rows = ['x_m,z_m,half_breadth_m']
        for x in (-1, 0, 1):
            rows.extend([f'{x},1,0', f'{x},1,0.5', f'{x},0,0.5'])
        (self.base/'box.csv').write_text('\n'.join(rows)+'\n')
        self.config = dict(geometry_file='box.csv', samples_per_section=4,
                           submerged=False, center_of_mass_m=[0., 0., .5],
                           mass_kg=2000., R44_m=.5, R55_m=.7, R66_m=.7,
                           omega_rad_s=[.8], reflection_symmetry=True,
                           viscous_damping=dict(T_126=[5, 2, 2], delta_zeta_345=[0, .1, 0]))

    def test_radiation_limits_and_diffraction_match_expanded_mesh(self):
        v, f, _ = load_hull(self.config, self.base, False)
        full_mesh = cpt.Mesh(vertices=v, faces=f)
        mesh = reflection_mesh(full_mesh)
        lid = reflection_mesh(bounded_lid(full_mesh, .3))
        cg = [0., 0., -.5]
        symmetric = cpt.FloatingBody(mesh=mesh, lid_mesh=lid,
            dofs=cpt.rigid_body_dofs(rotation_center=cg), center_of_mass=cg, mass=2000.)
        expanded = cpt.FloatingBody(mesh=mesh.merged(), lid_mesh=lid.merged(),
            dofs=cpt.rigid_body_dofs(rotation_center=cg), center_of_mass=cg, mass=2000.)
        self.assertIsInstance(symmetric.mesh_including_lid, cpt.ReflectionSymmetricMesh)
        self.assertAlmostEqual(symmetric.mesh.volume, expanded.mesh.volume)
        self.assertAlmostEqual(full_mesh.volume, symmetric.mesh.volume)
        solvers = [cpt.BEMSolver(), cpt.BEMSolver()]
        for omega in (0., .8, np.inf):
            matrices = [dict(added_masses=[], radiation_dampings=[]) for _ in range(2)]
            for dof in symmetric.dofs:
                results = [solver.solve(cpt.RadiationProblem(body=body, omega=omega,
                    radiating_dof=dof, rho=1025., g=9.81), keep_details=False)
                    for body, solver in zip((symmetric, expanded), solvers)]
                for key in ('added_masses', 'radiation_dampings'):
                    for index, result in enumerate(results):
                        matrices[index][key].append(list(getattr(result, key).values()))
            for key in ('added_masses', 'radiation_dampings'):
                a, b = [np.asarray(matrix[key]) for matrix in matrices]
                # Compare full coefficient matrices rather than relative
                # errors in near-zero individual modes or couplings.
                self.assertLessEqual(np.linalg.norm(a-b), 1e-3*max(np.linalg.norm(b), 1.))
                parity = np.array([1., -1., 1., -1., 1., -1.])
                np.testing.assert_allclose(a[parity[:, None]*parity[None, :] < 0], 0., atol=1e-7)
        for heading in (0., np.pi/3, np.pi):
            results = [solver.solve(cpt.DiffractionProblem(body=body, omega=.8,
                wave_direction=heading, rho=1025., g=9.81), keep_details=False)
                for body, solver in zip((symmetric, expanded), solvers)]
            a, b = [np.asarray(list(result.forces.values())) for result in results]
            self.assertLessEqual(np.linalg.norm(a-b), 1e-3*max(np.linalg.norm(b), 1.))
            if heading in (0., np.pi):
                np.testing.assert_allclose(a[[1, 3, 5]], 0., atol=1e-7)
        solver = solvers[0]
        symmetric_estimate = solver.engine.compute_ram_estimation(cpt.RadiationProblem(
            body=symmetric, omega=.8, radiating_dof='Heave'))
        expanded_estimate = solver.engine.compute_ram_estimation(cpt.RadiationProblem(
            body=expanded, omega=.8, radiating_dof='Heave'))
        self.assertLess(symmetric_estimate, expanded_estimate)

    def test_reject_invalid_symmetry_settings_before_solving(self):
        for edits, message in ((dict(reflection_symmetry='true'), 'true or false'),
                               (dict(center_of_mass_m=[0, .1, .5]), 'centerline')):
            config = dict(self.config, **edits)
            path = self.base/'config.json'
            path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, message):
                run(path)

    def test_symmetric_pipeline_exports_the_full_solver_mesh(self):
        config = dict(self.config, body_name='box', mass_kg=2050., output_dir='support_files')
        path = self.base/'config.json'
        path.write_text(json.dumps(config))
        output = run(path)
        self.assertEqual(output, (self.base/'box.mat').resolve())
        vessel = loadmat(output, simplify_cells=True)['vessel']
        panels = loadmat(self.base/'support_files/generated_hull_panels.mat')
        self.assertTrue((self.base/'support_files/hydrostatics.json').is_file())
        self.assertFalse((self.base/'support_files/box.mat').exists())
        self.assertEqual(vessel['hydrodynamic_source'], 'Capytaine')
        self.assertTrue(np.isfinite(vessel['A']).all())
        np.testing.assert_array_equal(vessel['drift']['w'], [.8])
        self.assertEqual(len(vessel['drift']['coeff']), 3)
        for coefficient in vessel['drift']['coeff']:
            self.assertEqual(coefficient.shape, (36,)) # singleton frequency/speed squeezed by scipy
            self.assertTrue(np.isfinite(coefficient).all())
        np.testing.assert_array_equal(vessel['freqs'], [0., .8, 10.])
        np.testing.assert_allclose(panels['vertices'].min(axis=0), [-1., -.5, 0.])
        np.testing.assert_allclose(panels['vertices'].max(axis=0), [1., .5, 1.])
        self.assertTrue(np.all(panels['faces'] >= 1))


if __name__ == '__main__':
    unittest.main()
