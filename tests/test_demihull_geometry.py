"""Check separated CSV hulls, common loading and the approximate Otter case."""
import json
from pathlib import Path
import tempfile
import unittest

import capytaine as cpt
import numpy as np

from src.hull_geometry import load_hull
from src.hydrodynamics import _mass_matrix
from src.output_coordinates import transform_matrix

ROOT = Path(__file__).resolve().parents[1]


class DemihullTest(unittest.TestCase):
    def test_two_separate_copies_preserve_volume_and_panel_winding(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            rows = ["x_m,z_m,half_breadth_m"]
            for x in (-1, 0, 1):
                rows.extend([f"{x},0.2,0", f"{x},0.2,0.1", f"{x},0,0.1"])
            (base / "box.csv").write_text("\n".join(rows)+"\n")
            config = {"geometry_file": "box.csv", "samples_per_section": 4}
            single_v, single_f, _ = load_hull(config, base, False)
            config["demihull_spacing_m"] = 0.8
            v, f, _ = load_hull(config, base, False)
            count, panels = len(single_v), len(single_f)
            np.testing.assert_allclose(v[:count], single_v+[0, -0.4, 0])
            np.testing.assert_allclose(v[count:], single_v+[0, 0.4, 0])
            np.testing.assert_array_equal(f[:panels], single_f)
            np.testing.assert_array_equal(f[panels:], single_f+count)
            self.assertTrue(np.all(f[:panels] < count))
            self.assertTrue(np.all(f[panels:] >= count))
            single = cpt.Mesh(vertices=single_v, faces=single_f)
            twin = cpt.Mesh(vertices=v, faces=f)
            self.assertGreater(single.volume, 0)
            self.assertAlmostEqual(twin.volume, 2*single.volume)
            np.testing.assert_allclose(twin.faces_normals[:panels], single.faces_normals)
            self.assertLess(v[:count, 1].max(), 0)
            self.assertGreater(v[count:, 1].min(), 0)
            for spacing in (True, None, "0.8", 0, -1, float("nan"), float("inf")):
                with self.subTest(spacing=spacing), self.assertRaisesRegex(ValueError, "positive finite"):
                    load_hull(dict(config, demihull_spacing_m=spacing), base, False)
            for spacing in (0.2, 0.1):
                with self.assertRaisesRegex(ValueError, "touch or overlap"):
                    load_hull(dict(config, demihull_spacing_m=spacing), base, False)

    def test_spacing_is_csv_only(self):
        for case in ("tanker_gdf", "supply_mgf", "wigley_stl"):
            base = ROOT / "vessels_capytaine" / case
            config = json.loads((base / "config.json").read_text())
            config["demihull_spacing_m"] = 0.79
            with self.subTest(case=case), self.assertRaisesRegex(ValueError, "only supported for CSV"):
                load_hull(config, base, False)

    def test_otter_loading_and_lid_gap(self):
        base = ROOT / "vessels_capytaine/otter_80kg_csv"
        config = json.loads((base / "config.json").read_text())
        v, f, _ = load_hull(config, base, False)
        mesh = cpt.Mesh(vertices=v, faces=f)
        cg = np.asarray(config["center_of_mass_m"]) * [-1, 1, -1]
        body = cpt.FloatingBody(mesh=mesh, center_of_mass=cg, mass=config["mass_kg"],
                                dofs=cpt.rigid_body_dofs(rotation_center=cg))
        self.assertEqual(config["mass_kg"], 80)
        self.assertAlmostEqual(mesh.volume*1025, 80, places=5)
        self.assertAlmostEqual(mesh.waterplane_area, 0.75, places=6)
        self.assertAlmostEqual(body.center_of_buoyancy[0], cg[0], places=6)
        self.assertAlmostEqual(body.center_of_buoyancy[1], 0, places=6)
        stiffness = body.compute_hydrostatic_stiffness(rho=1025, g=9.81).values
        self.assertGreater(stiffness[3, 3], 0)
        self.assertGreater(stiffness[4, 4]-stiffness[2, 4]**2/stiffness[2, 2], 0)
        # The fixed hull's unloaded waterline supplies the origin shift.
        self.assertAlmostEqual(config["center_of_mass_m"][2], -0.10386139215, delta=5.1e-5)
        inertia = transform_matrix(_mass_matrix(config, base, 80))
        reference = np.array([[10.95182, 0, 0.6875], [0, 15.125, 0], [0.6875, 0, 14.4375]])
        # Four-decimal JSON radii retain the MSS inertia within 0.01%.
        self.assertLess(np.linalg.norm(inertia[3:, 3:]-reference)/np.linalg.norm(reference), 1e-4)
        lid = mesh.generate_lid()
        self.assertGreater(lid.nb_faces, 0)
        inner_edge = (config["demihull_spacing_m"]-0.25)/2
        self.assertTrue(np.all(np.abs(lid.vertices[:, 1]) >= inner_edge-1e-9))
        self.assertTrue(np.all(np.ptp(lid.vertices[lid.faces, 1], axis=1) < 0.25))
        wavelength_eighth = 2*np.pi*9.81/max(config["omega_rad_s"])**2/8
        self.assertLess(max(mesh.faces_radiuses.max(), lid.faces_radiuses.max()), wavelength_eighth)

    def test_light_otter_static_trim_and_rotated_inertia(self):
        base = ROOT / "vessels_capytaine/otter_55kg_csv"
        config = json.loads((base / "config.json").read_text())
        v, f, _ = load_hull(config, base, False)
        mesh = cpt.Mesh(vertices=v, faces=f)
        cg = np.array(config["center_of_mass_m"])*[-1, 1, -1]
        body = cpt.FloatingBody(mesh=mesh, mass=55, center_of_mass=cg,
                                dofs=cpt.rigid_body_dofs(rotation_center=cg))
        self.assertEqual(config["mass_kg"], 55)
        self.assertAlmostEqual(mesh.volume*1025, 55, places=5)
        self.assertAlmostEqual(body.center_of_buoyancy[0], cg[0], delta=5.1e-5)
        inertia = transform_matrix(_mass_matrix(config, base, 55))[3:, 3:]
        np.testing.assert_allclose(np.linalg.eigvalsh(inertia), [10.26432, 13.75, 13.75], atol=2e-4)
        pitch = .5*np.arctan2(2*inertia[0, 2], inertia[2, 2]-inertia[0, 0])
        self.assertAlmostEqual(np.rad2deg(pitch), -0.28851, delta=0.005)
        self.assertAlmostEqual(config["center_of_mass_m"][0], .2*(np.cos(pitch)-np.sin(pitch)), delta=5.1e-5)
        C = body.compute_hydrostatic_stiffness(rho=1025, g=9.81).values
        self.assertGreater(C[3, 3], 0)
        self.assertGreater(C[4, 4]-C[2, 4]**2/C[2, 2], 0)
        lid = mesh.generate_lid()
        self.assertTrue(np.all(np.abs(lid.vertices[:, 1]) >= .27-1e-9))


if __name__ == "__main__":
    unittest.main()
