"""Tests for scripts/check_structure.py (Task 1)."""
import os
import unittest

import helpers
import herlib
from check_structure import check_dir, check_geometry


class Geometry(unittest.TestCase):
    def test_template_site_found_by_geometry(self):
        s = helpers.structure()
        m, ligs = herlib.find_site(s)
        self.assertEqual(s["elements"][m], "Ir")
        self.assertEqual(sorted(s["elements"][j] for j in ligs), ["C", "C", "N", "S"])

    def test_minimum_image(self):
        s = helpers.structure()
        s["frac"] = [[0.01, 0.0, 0.5], [0.99, 0.0, 0.5]]
        s["elements"] = ["C", "C"]
        self.assertAlmostEqual(herlib.distance(s, 0, 1), 0.02 * 9.84, places=6)

    def test_clean_h_on_metal_passes(self):
        rep, flags, warns = check_geometry(helpers.structure([(0, 0, 1.6)]))
        self.assertEqual(flags, [])
        self.assertEqual(rep["coordination"], 4)
        self.assertEqual(len(rep["ligands"]), 4)
        # mixed C/N/S ligands: spread > 0.10 gives a Jahn–Teller warning, never a flag
        self.assertTrue(any("Jahn–Teller" in w for w in warns))

    def test_h2_formed(self):
        _, flags, _ = check_geometry(helpers.structure([(0, 0, 1.6), (0.75, 0, 1.6)]))
        self.assertTrue(any("dihydrogen" in f for f in flags))

    def test_h_in_vacuum(self):
        _, flags, _ = check_geometry(helpers.structure([(0, 0, 5.0)]))
        self.assertTrue(any("not bound" in f for f in flags))
        self.assertTrue(any("vacuum" in f for f in flags))

    def test_metal_out_of_plane(self):
        rep, flags, _ = check_geometry(helpers.structure(metal_shift=(0, 0, 1.2)))
        self.assertGreater(rep["out_of_plane"], 0.5)
        self.assertTrue(any("out of the ligand plane" in f for f in flags))

    def test_broken_bond(self):
        s = helpers.structure()
        m, ligs = herlib.find_site(s)
        far = ligs[-1]
        vec = herlib.mic_vector(s, m, far)
        rep, flags, _ = check_geometry(helpers.structure(ligand_moves={far: [1.0 * x for x in vec]}))
        self.assertEqual(rep["coordination"], 3)
        self.assertTrue(any("3-coordinated" in f for f in flags))

    def test_relaxed_template_out_of_plane(self):
        # the relaxed Ir-CCNS@G template: Ir 0.605 Å from its least-squares ligand plane
        rep, flags, _ = check_geometry(helpers.structure(flat=False))
        self.assertAlmostEqual(rep["out_of_plane"], 0.605, places=3)
        self.assertEqual(len(flags), 1)
        self.assertEqual(check_geometry(helpers.structure(flat=False), {"max_out_of_plane": 0.7})[1], [])

    def test_thresholds_are_options(self):
        _, _, warns = check_geometry(helpers.structure([(0, 0, 1.6)]), {"jt_spread": 1.0})
        self.assertFalse(any("Jahn–Teller" in w for w in warns))


class Dirs(unittest.TestCase):
    def test_converged_relaxation(self):
        d = helpers.relax_job(helpers.structure([(0, 0, 1.6)]), mag=1.234)
        r = check_dir(d, use_scheduler=False)
        self.assertEqual(r["status"], "pass")
        self.assertAlmostEqual(r["mag"], 1.234)

    def test_unconverged_not_judged(self):
        d = helpers.relax_job(helpers.structure(), extra={"log": "still going\n"})
        self.assertEqual(check_dir(d, use_scheduler=False)["status"], "not_judged")

    def test_frequency_imaginary_modes(self):
        d = helpers.freq_job(helpers.structure([(0, 0, 1.6)]), 0.3, 0.01, imaginary=(52.3,))
        r = check_dir(d, use_scheduler=False)
        self.assertEqual(r["status"], "flag")
        self.assertEqual([x["cm-1"] for x in r["imaginary"]], [52.3])
        r = check_dir(d, {"imag_tol": 60.0}, use_scheduler=False)
        self.assertEqual(r["status"], "pass")

    def test_missing_contcar(self):
        d = helpers.relax_job(helpers.structure())
        os.remove(os.path.join(d, "CONTCAR"))
        self.assertEqual(check_dir(d, use_scheduler=False)["status"], "flag")


if __name__ == "__main__":
    unittest.main()
