"""Check damping migration and MAT refresh without running a BEM solve."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
from scipy.io import loadmat

from src.damping import refresh_damping, validate_damping

ROOT = Path(__file__).resolve().parents[1]


def assert_mat_equal(test, actual, expected):
    test.assertEqual(actual.shape, expected.shape)
    test.assertEqual(actual.dtype, expected.dtype)
    if actual.dtype.names:
        for key in actual.dtype.names:
            assert_mat_equal(test, actual[key], expected[key])
    elif actual.dtype == object:
        for a, b in zip(actual.flat, expected.flat):
            assert_mat_equal(test, a, b)
    else:
        np.testing.assert_array_equal(actual, expected)


class DampingTest(unittest.TestCase):
    def test_legacy_percentages_have_migration_error(self):
        with self.assertRaisesRegex(ValueError, 'Replace.*T_126.*seconds'):
            validate_damping({'kappa_126': [.05]*3, 'delta_zeta_345': [0,.1,0]})

    def test_time_constants_and_vessel_type(self):
        valid = {'T_126': [20,2.5,5], 'delta_zeta_345': [0,.1,0]}
        np.testing.assert_array_equal(validate_damping(valid)['T_126'], [20,2.5,5])
        for invalid in ([0,1,1], [-1,1,1], [float('nan'),1,1], [1,1],
                        [True,1,1], ['20',1,1], [[1,1,1]], [float('inf'),1,1]):
            with self.subTest(values=invalid), self.assertRaises(ValueError):
                validate_damping(dict(valid,T_126=invalid))
        with self.assertRaises(ValueError):
            validate_damping(valid,submerged=True)
        with self.assertRaises(ValueError):
            validate_damping(dict(valid,delta_zeta_345=[0,-.1,0]))
        with self.assertRaises(ValueError):
            validate_damping(dict(valid,T_1236=[1]*4))
        np.testing.assert_array_equal(validate_damping(valid)['nonlinear_456'], [0,0,0])
        for invalid in ([-1,0,0], [1,2], [True,0,0], [float('inf'),0,0]):
            with self.subTest(nonlinear=invalid), self.assertRaises(ValueError):
                validate_damping(dict(valid,nonlinear_456=invalid))

    def test_all_catalogue_inputs(self):
        for path in (ROOT/'vessels_capytaine').glob('*/config.json'):
            with self.subTest(case=path.parent.name):
                config = json.loads(path.read_text())
                validate_damping(config['viscous_damping'],config['submerged'])

    def test_refresh_preserves_hydrodynamics_and_cells(self):
        source = ROOT/'vessels_capytaine/supply_mgf/supply.mat'
        before = loadmat(source)['vessel']
        with tempfile.TemporaryDirectory() as folder:
            base = Path(folder)
            shutil.copy2(source,base/'supply.mat')
            config = {'submerged': False,'output_dir': 'support_files', 'output_filename': 'supply.mat',
                      'viscous_damping': {'T_126': [20,2.5,5],'delta_zeta_345':[0,.1,0],
                                         'nonlinear_456':[5,0,0]}}
            (base/'config.json').write_text(json.dumps(config))
            result = refresh_damping(base/'config.json')
            after = loadmat(result)['vessel']
            for key in before.dtype.names:
                if key != 'powerBased':
                    assert_mat_equal(self,after[key],before[key])
            p = loadmat(result,simplify_cells=True)['vessel']['powerBased']
            self.assertNotIn('kappa_126',p)
            np.testing.assert_array_equal(p['T_126'],[20,2.5,5])
            np.testing.assert_array_equal(p['nonlinear_456'],[5,0,0])
            config['viscous_damping'] = {'kappa_126':[.05]*3}
            (base/'config.json').write_text(json.dumps(config))
            saved = result.read_bytes()
            with self.assertRaisesRegex(ValueError,'Legacy'):
                refresh_damping(base/'config.json')
            self.assertEqual(saved,result.read_bytes())


if __name__ == '__main__':
    unittest.main()
