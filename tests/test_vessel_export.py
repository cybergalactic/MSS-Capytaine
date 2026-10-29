"""Check the exported MSS contract without a hydrodynamic solve."""
from pathlib import Path
import json
import tempfile
import unittest

import numpy as np
from scipy.io import loadmat
import xarray as xr

from src.capytaine_vessel import DOFS, write_vessel
from plot_results import coefficient, drift_data, load_vessel, rao_data


class VesselExportTest(unittest.TestCase):
    def test_saved_catalogue_files_match_python_reader_contract(self):
        catalogue = Path(__file__).resolve().parents[1]/'vessels_capytaine'
        for config_path in sorted(catalogue.glob('*/config.json')):
            with self.subTest(case=config_path.parent.name):
                config = json.loads(config_path.read_text())
                vessel = load_vessel(config_path.parent/config['output_filename'])
                frequencies = np.asarray(vessel.freqs).ravel()
                headings = np.asarray(vessel.headings).ravel()
                for name in ('A', 'B', 'C'):
                    values = coefficient(vessel, name, frequencies)
                    self.assertTrue(np.isfinite(values).all())
                for name in ('forceRAO', 'motionRAO'):
                    _, amplitude, phase = rao_data(vessel, name, headings, frequencies)
                    self.assertTrue(np.isfinite(amplitude).all())
                    self.assertTrue(np.isfinite(phase).all())
                _, drift = drift_data(vessel, headings, frequencies)
                self.assertTrue(np.isfinite(drift).all())

    def test_basic_inputs_produce_required_fields_and_separate_limits(self):
        omega = np.array([.3, .8, 1.5])
        headings = np.array([0., np.pi/2, np.pi])
        mass = np.diag([100., 100., 100., 400., 900., 1600.])
        stiffness = np.diag([0., 0., 100., 200., 300., 0.])
        ds = xr.Dataset(coords={'omega': omega, 'wave_direction': headings,
                                'influenced_dof': list(DOFS), 'radiating_dof': list(DOFS)})
        ds['inertia_matrix'] = (('influenced_dof', 'radiating_dof'), mass)
        ds['hydrostatic_stiffness'] = (('influenced_dof', 'radiating_dof'), stiffness)
        for key in ('added_mass', 'radiation_damping'):
            ds[key] = (('omega', 'influenced_dof', 'radiating_dof'),
                       np.repeat(np.eye(6)[None], len(omega), axis=0))
        for key, dof in (('excitation_force', 'influenced_dof'), ('motion_rao', 'radiating_dof')):
            ds[key] = (('omega', 'wave_direction', dof), np.ones((3, 3, 6), dtype=complex))
        ds.attrs.update(matrix_reference='CG', coordinate_system='MSS FSD')
        ds['mean_drift_force'] = (('omega','wave_direction','drift_dof'),
            np.tile([[1.,0.,0.],[2.,-3.,4.],[-1.,0.,0.]], (3,1,1)))
        ds = ds.assign_coords(drift_dof=['Surge','Sway','Yaw'])
        hydro = dict(body_mass_kg=100., volume_m3=.1, center_of_mass_m=[0., 0., 0.],
                     center_of_buoyancy_m=[0., 0., 1.], center_of_flotation_m=[0., 0., 0.],
                     demihull_spacing_m=2.)
        vertices = np.array([[-2., -1., -1.], [2., -1., -1.], [0., 1., 0.]])
        faces = np.array([[0, 1, 2, 2]])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'basic.mat'
            write_vessel(ds, hydro, vertices, faces, 'basic', 1025., 9.81, path,
                         zero_added_mass=2*np.eye(6), zero_radiation_damping=np.zeros((6, 6)),
                         infinite_added_mass=3*np.eye(6), infinite_radiation_damping=np.zeros((6, 6)),
                         T_126=np.array([10., 2., 3.]), delta_zeta_345=np.array([0., .1, 0.]),
                         nonlinear_456=np.array([5., 0., 0.]))
            vessel = loadmat(path, simplify_cells=True)['vessel']
            raw = loadmat(path)['vessel'][0,0]['drift'][0,0]
            self.assertTrue((path.parent/'results/generated_hull_panels.mat').is_file())
            self.assertTrue((path.parent/'results/generated_hull_panels.npz').is_file())
            self.assertFalse((path.parent/'generated_hull_panels.mat').exists())
        np.testing.assert_allclose([vessel['main'][k] for k in ('k44', 'k55', 'k66')], [2., 3., 4.])
        for key in ('S', 'Lwl', 'demihullSpacing'):
            self.assertNotIn(key, vessel['main'])
        for key in ('roll', 'driftfrc', 'hydrodynamic_axes',
                    'viscous_damping_model', 'second_order_drift'):
            self.assertNotIn(key, vessel)
        self.assertEqual(vessel['hydrodynamic_source'], 'Capytaine')
        np.testing.assert_array_equal(vessel['powerBased']['nonlinear_456'], [5.,0.,0.])
        for key in ('m', 'CG', 'CB', 'CF', 'GM_T', 'GM_L', 'C_B', 'nabla', 'Lpp', 'B', 'T'):
            self.assertIn(key, vessel['main'])
        np.testing.assert_array_equal(vessel['freqs'], [0., .3, .8, 1.5, 10.])
        np.testing.assert_array_equal(vessel['forceRAO']['w'], omega)
        np.testing.assert_array_equal(vessel['motionRAO']['w'], omega)
        np.testing.assert_array_equal(vessel['A'][:, :, 0], 2*np.eye(6))
        np.testing.assert_array_equal(vessel['A'][:, :, -1], 3*np.eye(6))
        np.testing.assert_array_equal(vessel['drift']['w'], omega)
        self.assertEqual(raw['coeff'].shape, (1,3))
        for cell in raw['coeff'].flat:
            self.assertEqual(cell.shape, (3,4,1))
        np.testing.assert_array_equal(vessel['drift']['coeff'][0], [[1,2,-1,2]]*3)
        np.testing.assert_array_equal(vessel['drift']['coeff'][1], [[0,-3,0,3]]*3)
        np.testing.assert_array_equal(vessel['drift']['coeff'][2], [[0,4,0,-4]]*3)


if __name__ == '__main__':
    unittest.main()
