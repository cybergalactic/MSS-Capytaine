"""MGF/STL geometry, flat catalogue lookup and hydrostatics checks."""
import json
from pathlib import Path
import struct
import tempfile
import unittest

import capytaine as cpt
import numpy as np

from main import _config_path
from src.hull_geometry import load_hull, bounded_lid
from src.shipx_geometry import read_mgf_sections
from src.stl_geometry import read_stl
from src.geometry_format import geometry_format
from src.offset_geometry import refine_large_panels

ROOT = Path(__file__).resolve().parents[1]


class FormatsTest(unittest.TestCase):
    def test_bounded_lid_matches_native_lid(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            rows = ['x_m,z_m,half_breadth_m']
            for x in (-1, 0, 1):
                rows.extend([f'{x},1,0', f'{x},1,0.5', f'{x},0,0.5'])
            (base/'box.csv').write_text('\n'.join(rows)+'\n')
            v, f, _ = load_hull({'geometry_file': 'box.csv', 'samples_per_section': 4}, base, False)
            mesh = cpt.Mesh(vertices=v, faces=f)
            native = mesh.generate_lid(faces_max_radius=.2)
            batched = bounded_lid(mesh, .2)
            np.testing.assert_array_equal(batched.vertices, native.vertices)
            np.testing.assert_array_equal(batched.faces, native.faces)

    def test_local_refinement_preserves_surface(self):
        vertices = np.array([[0., 0., -1.], [4., 0., -1.], [0., 2., -1.]])
        faces = np.array([[0, 1, 2, 2]])
        v, f = refine_large_panels(vertices, faces, .4)
        triangles = v[f[:, :3]]
        radius = np.linalg.norm(triangles-triangles.mean(axis=1)[:, None], axis=2)
        self.assertLessEqual(radius.max(), .4)
        areas = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
        self.assertTrue(np.all(areas[:, 2] > 0))
        self.assertAlmostEqual(np.linalg.norm(areas, axis=1).sum()/2, 4.)
        np.testing.assert_allclose(v[:, 2], -1.)
        for invalid in (True, '0.4', None, 0, -1, float('nan')):
            with self.assertRaises(ValueError):
                refine_large_panels(vertices, faces, invalid)

    def test_missing_extensions(self):
        for case in ('testShip_csv', 'supply_mgf',
                     'tanker_gdf', 'wigley_stl'):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as directory:
                base = ROOT / 'vessels_capytaine' / case
                c = json.loads((base / 'config.json').read_text())
                source = base / c['geometry_file']
                path = Path(directory) / 'hull_geometry'
                path.write_bytes(source.read_bytes())
                self.assertEqual(geometry_format(path), source.suffix)
                original_v, original_f, _ = load_hull(c, base, c['submerged'])
                c['geometry_file'] = str(path)
                v, f, _ = load_hull(c, base, c['submerged'])
                np.testing.assert_array_equal(v, original_v)
                np.testing.assert_array_equal(f, original_f)

    def test_folder_and_body_name_lookup(self):
        for case in ('testShip_csv', 'tanker_gdf',
                     'supply_mgf', 'wigley_stl'):
            config = json.loads(_config_path(case).read_text())
            self.assertEqual(_config_path(case), _config_path(config['body_name']))
        self.assertEqual(_config_path('custom.json'), Path('custom.json'))

    def test_mgf_clips_draft_preserves_flat_keel(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'box.mgf'
            text = 'Box\n\n\n\n4\n'
            for i, x in enumerate((-2, 0, 2), 1):
                text += f'{i}\n{x}\n4\n0 2\n1 2\n1 0\n0 0\n'
            path.write_text(text)
            sections = read_mgf_sections(path, 1.)
            np.testing.assert_allclose(sections[0][1], [[0, -1], [1, -1], [1, 0]])
            v, f, _ = load_hull({'geometry_file': path.name, 'draft_m': 1.,
                                 'samples_per_section': 3}, path.parent, False)
            self.assertGreater(cpt.Mesh(vertices=v, faces=f).volume, 0)
            for draft in (0, -1, float('nan')):
                with self.assertRaises(ValueError):
                    read_mgf_sections(path, draft)
            with self.assertRaisesRegex(ValueError, 'draft_m'):
                load_hull({'geometry_file': path.name}, path.parent, False)
            with self.assertRaisesRegex(ValueError, 'surface vessels'):
                load_hull({'geometry_file': path.name, 'draft_m': 1}, path.parent, True)
            path.write_text('Bad\n\n\n\n4\n1\n0\n4\n0 0\n')
            with self.assertRaises(ValueError):
                read_mgf_sections(path, 1)

    def test_ascii_binary_stl_agree_in_fsd(self):
        vertices = np.array([[0., 0, 0], [2., 0, 0], [0, 1., 0], [0, 0, 1.]])
        triangles = vertices[[[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'hull.stl'
            text = 'solid hull\n'
            for triangle in triangles:
                text += 'facet normal 0 0 0\nouter loop\n'
                text += ''.join('vertex '+' '.join(map(str, v))+'\n' for v in triangle)
                text += 'endloop\nendfacet\n'
            path.write_text(text+'endsolid hull\n')
            a, fa = read_stl(path)
            # Binary headers may begin with 'solid'; length/count detection wins.
            data = bytearray(b'solid binary'.ljust(80, b' '))
            data.extend(struct.pack('<I', len(triangles)))
            for triangle in triangles:
                data.extend(struct.pack('<12fH', 0, 0, 0, *triangle.ravel(), 0))
            path.write_bytes(data)
            b, fb = read_stl(path)
            np.testing.assert_array_equal(a, b)
            np.testing.assert_array_equal(fa, fb)
            v, f, _ = load_hull({'geometry_file': path.name}, path.parent, False)
            self.assertLessEqual(v[:, 0].max(), 0)
            self.assertLessEqual(v[:, 2].max(), 0)
            self.assertAlmostEqual(cpt.Mesh(vertices=v, faces=f).volume, 1/3)
            path.write_bytes(data[:-1])
            with self.assertRaises(ValueError):
                read_stl(path)

    def test_mgf_and_wigley_hydrostatics(self):
        for case in ('supply_mgf', 's175_mgf', 'wigley_stl'):
            with self.subTest(case=case):
                base = ROOT / 'vessels_capytaine' / case
                c = json.loads((base / 'config.json').read_text())
                v, f, _ = load_hull(c, base, False)
                mesh = cpt.Mesh(vertices=v, faces=f)
                if 'max_panel_radius_m' in c:
                    coarse = dict(c)
                    radius = coarse.pop('max_panel_radius_m')
                    cv, cf, _ = load_hull(coarse, base, False)
                    lid = bounded_lid(cpt.Mesh(vertices=cv, faces=cf), radius)
                else:
                    lid = mesh.generate_lid()
                largest_radius = max(mesh.faces_radiuses.max(), lid.faces_radiuses.max())
                wavelength_eighth = 2*np.pi*c['gravity_m_s2'] / max(c['omega_rad_s'])**2 / 8
                self.assertLess(largest_radius, wavelength_eighth)
                cg = np.asarray(c['center_of_mass_m']) * [-1, 1, -1]
                body = cpt.FloatingBody(mesh=mesh,
                    dofs=cpt.rigid_body_dofs(rotation_center=cg),
                    center_of_mass=cg, mass=c['mass_kg'])
                self.assertLess(abs(mesh.volume*c['water_density_kg_m3']/c['mass_kg']-1), .01)
                stiffness = np.asarray(body.compute_hydrostatic_stiffness(rho=1025, g=9.81))
                self.assertTrue(np.isfinite(stiffness).all())
                self.assertGreater(stiffness[3, 3], 0)
                self.assertGreater(stiffness[4, 4], 0)


if __name__ == '__main__':
    unittest.main()
