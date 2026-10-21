"""Regression tests: python -m unittest discover -s tests"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from her_analysis import analyse  # noqa: E402
from parse_candidates import canonical, parse  # noqa: E402
from check_jobs import check  # noqa: E402
import tempfile  # noqa: E402


class WorkedCheck(unittest.TestCase):
    def test_cr_psss(self):
        # references/dft_her.md §7: ΔG_H* = -0.352; 2H* levels -0.397 (C), -0.401 (M); G_* = G_H2 = 0
        r = analyse(0.0, -0.352, -0.397, -0.401, 0.0)
        self.assertAlmostEqual(r["eta_VH"], 0.352)
        self.assertAlmostEqual(r["eta_VTC"], 0.397)
        self.assertAlmostEqual(r["eta_VTM"], 0.401)
        self.assertEqual(r["pathway"], "VH")
        self.assertAlmostEqual(r["eta_limiting"], 0.352)


class Canonical(unittest.TestCase):
    def test_known_names_are_canonical(self):
        for lab in ("CCNS", "PSSS", "OPPS"):
            self.assertEqual(canonical(lab), lab)

    def test_equivalent_orderings(self):
        self.assertEqual({canonical(x) for x in ("PONS", "SNOP", "OPSN", "NSPO")}, {"NSPO"})


class Parse(unittest.TestCase):
    def test_missing_dg_flagged(self):
        rows = parse("Top\n配置 1:\n  Z=44, 配体=[15, 8, 7, 16]  (PONS)\n")
        self.assertEqual(rows[0]["label_check"], "MISSING dG")
        self.assertEqual(rows[0]["canonical"], "Ru-NSPO@G")


def job(files):
    d = tempfile.mkdtemp()
    for name, text in files.items():
        with open(os.path.join(d, name), "w") as fh:
            fh.write(text)
    return d


class CheckJobs(unittest.TestCase):
    def test_converged(self):
        d = job({"INCAR": "IBRION = 2\nNSW = 500\n",
                 "log": "  1 F= -1\n reached required accuracy - stopping structural energy minimisation\n"})
        self.assertEqual(check(d, use_scheduler=False)["status"], "converged")

    def test_zbrent(self):
        d = job({"INCAR": "IBRION = 2\nNSW = 500\n", "log": "ZBRENT: fatal error in bracketing\n"})
        r = check(d, use_scheduler=False)
        self.assertEqual((r["status"], r["error"]), ("failed", "zbrent"))

    def test_nsw_reached(self):
        d = job({"INCAR": "IBRION = 2\nNSW = 3\n", "log": "x\n",
                 "OSZICAR": "DAV: 1\n   3 F= -.1E+03 E0= -.1E+03  d E =-.1E-02  mag=  1.0\n"})
        r = check(d, use_scheduler=False)
        self.assertEqual(r["status"], "unconverged")
        self.assertIn("CONTCAR", r["fix"])

    def test_frequency_done(self):
        d = job({"INCAR": "IBRION = 5  # finite differences\nNSW = 1\n",
                 "OUTCAR": " General timing and accounting informations for this job:\n"})
        self.assertEqual(check(d, use_scheduler=False)["status"], "converged")


if __name__ == "__main__":
    unittest.main()
