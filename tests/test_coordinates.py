"""Fast convention checks; no BEM solves are run."""

import json
from pathlib import Path
import tempfile
import unittest

import capytaine as cpt
import numpy as np
from scipy.io import loadmat

from src.hull_geometry import load_hull
from src.hydrodynamics import _mass_matrix
from src.output_coordinates import (
    transform_matrix, transform_positions,
)

ROOT = Path(__file__).resolve().parents[1]


class CoordinatesTest(unittest.TestCase):
    def test_cg_and_coupled_inertia_roundtrip(self):
        np.testing.assert_array_equal(
            transform_positions([3.93, 0, -2.5]), [-3.93, 0, 2.5]
        )
        # Cross terms exercise signs hidden by the diagonal catalogue inputs.
        mass = np.diag([100., 100., 100., 40., 50., 60.])
        mass[0, 4] = mass[4, 0] = 2.
        mass[1, 3] = mass[3, 1] = 3.
        mass[3, 5] = mass[5, 3] = 4.
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            np.savetxt(base / "inertia.csv", mass, delimiter=",")
            solver = _mass_matrix({"inertia_matrix_csv": "inertia.csv"}, base, 100.)
            with self.assertRaisesRegex(ValueError, "not inertia_matrix_csv"):
                _mass_matrix({"R46_m": 0.1, "inertia_matrix_csv": "inertia.csv"}, base, 100.)
        self.assertEqual(solver[0, 4], -2.)
        self.assertEqual(solver[1, 3], -3.)
        self.assertEqual(solver[3, 5], 4.)
        np.testing.assert_array_equal(transform_matrix(solver), mass)

    def test_optional_signed_cross_radius(self):
        config = {"R44_m": 0.4, "R55_m": 0.5, "R66_m": 0.6}
        diagonal = transform_matrix(_mass_matrix(config, Path('.'), 100.))
        np.testing.assert_allclose(diagonal, np.diag([100., 100., 100., 16., 25., 36.]))
        np.testing.assert_array_equal(_mass_matrix(config, Path('.'), 100.),
                                     _mass_matrix(dict(config, R46_m=0), Path('.'), 100.))
        for radius in (-0.2, 0.2):
            matrix = transform_matrix(_mass_matrix(dict(config, R46_m=radius), Path('.'), 100.))
            expected = diagonal.copy()
            expected[3, 5] = expected[5, 3] = np.copysign(4., radius)
            np.testing.assert_allclose(matrix, expected)
        for invalid in (True, "0.1", None, [], float('inf'), float('nan'), 1.0):
            with self.subTest(radius=invalid), self.assertRaises(ValueError):
                _mass_matrix(dict(config, R46_m=invalid), Path('.'), 100.)

    def test_asymmetric_csv_and_legacy_key(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base / "hull.CSV").write_text(
                "x_m,z_m,half_breadth_m\n"
                "-2,2,0\n-2,0,1\n0,3,0\n0,0,2\n4,1,0\n4,0,0.5\n"
            )
            v, f, _ = load_hull({"geometry_file": "hull.CSV"}, base, False)
            np.testing.assert_allclose(v.min(axis=0), [-4, -2, -3])
            np.testing.assert_allclose(v.max(axis=0), [2, 2, 0])
            v2, f2, _ = load_hull({"offset_points_csv": "hull.CSV"}, base, False)
            np.testing.assert_array_equal(v, v2)
            np.testing.assert_array_equal(f, f2)
            self.assertGreater(cpt.Mesh(vertices=v, faces=f).volume, 0)
            (base / "hull.CSV").write_text(
                "x_m,z_m,half_breadth_m\n0,-1,0\n"
            )
            with self.assertRaisesRegex(ValueError, "FSD z>=0"):
                load_hull({"geometry_file": "hull.CSV"}, base, False)

    def test_catalogue_config_and_panels(self):
        for case in ("capytaine_offsets_csv/testShip", "capytaine_offsets_csv/LAUV_marie", "wamit_gdf/tanker"):
            with self.subTest(case=case):
                base = ROOT / "vessels_capytaine" / case
                config = json.loads((base / "config.json").read_text())
                self.assertNotIn("coordinate_system", config)
                v, f, _ = load_hull(config, base, config["submerged"])
                self.assertTrue(np.isfinite(v).all())
                self.assertGreater(cpt.Mesh(vertices=v, faces=f).volume, 0)
                if case == "wamit_gdf/tanker":
                    original = cpt.load_mesh(str(base / "tanker.gdf")).merged()
                    self.assertEqual(len(f), 4750)
                    # The composed file->solver->public transform is WAMIT->FSD.
                    np.testing.assert_allclose(
                        transform_positions(v), original.vertices * [1, -1, -1]
                    )
                    np.testing.assert_array_equal(f, original.faces)
                else:
                    with np.load(base / "results/generated_hull_panels.npz") as panel:
                        self.assertEqual(str(panel["coordinate_system"]), "MSS_FSD")
                        np.testing.assert_allclose(panel["vertices"], transform_positions(v))
                    mat = loadmat(base / "results/generated_hull_panels.mat")
                    np.testing.assert_allclose(mat["vertices"], transform_positions(v))
                    np.testing.assert_array_equal(mat["faces"], f + 1)

    def test_geometry_errors(self):
        for config in ({}, {"geometry_file": "hull.obj"},
                       {"geometry_file": "a.gdf", "offset_points_csv": "b.csv"}):
            with self.assertRaises(ValueError):
                load_hull(config, ROOT, False)


if __name__ == "__main__":
    unittest.main()
