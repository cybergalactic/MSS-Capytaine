"""Check reference GM_T catches CG sign errors before any BEM solve."""
import json
from pathlib import Path
import unittest

import capytaine as cpt
import numpy as np

from src.hull_geometry import load_hull
from src.hydrodynamics import _check_transverse_stability

ROOT = Path(__file__).resolve().parents[1]


class StabilityTest(unittest.TestCase):
    def test_reference_and_error_message(self):
        c = {'center_of_mass_m': [3.93, 0, -2.5], 'GM_T_m': 9.9823}
        _check_transverse_stability(c, 9.98)
        c['center_of_mass_m'][2] = 2.5
        with self.assertRaisesRegex(ValueError, 'zG should be approximately -2.5'):
            _check_transverse_stability(c, 14.9823)
        with self.assertRaisesRegex(ValueError, 'not positive'):
            _check_transverse_stability(c, -1)
        with self.assertRaises(ValueError):
            _check_transverse_stability({**c, 'GM_T_tolerance_m': 0}, 9.98)

    def test_positive_gm_without_reference_does_not_prove_sign(self):
        # Both signs may be stable: no invented universal CG sign restriction.
        for z in (-2.5, 2.5):
            _check_transverse_stability({'center_of_mass_m': [0, 0, z]}, 10)

    def test_real_tanker_and_supply_signs(self):
        for case in ('tanker_gdf', 'supply_mgf'):
            with self.subTest(case=case):
                base = ROOT / 'vessels_capytaine' / case
                config = json.loads((base / 'config.json').read_text())
                v, f, _ = load_hull(config, base, False)
                mesh = cpt.Mesh(vertices=v, faces=f)
                for reverse in (False, True):
                    c = dict(config)
                    c['center_of_mass_m'] = list(config['center_of_mass_m'])
                    if reverse:
                        c['center_of_mass_m'][2] *= -1
                    cg = np.asarray(c['center_of_mass_m']) * [-1, 1, -1]
                    body = cpt.FloatingBody(mesh=mesh, mass=c['mass_kg'],
                        center_of_mass=cg, dofs=cpt.rigid_body_dofs(rotation_center=cg))
                    stiffness = np.asarray(body.compute_hydrostatic_stiffness(
                        rho=c['water_density_kg_m3'], g=c['gravity_m_s2']))
                    gm = stiffness[3, 3] / (c['water_density_kg_m3']*c['gravity_m_s2']*mesh.volume)
                    if reverse:
                        with self.assertRaisesRegex(ValueError, 'GM_T consistency check failed'):
                            _check_transverse_stability(c, gm)
                    else:
                        _check_transverse_stability(c, gm)


if __name__ == '__main__':
    unittest.main()
