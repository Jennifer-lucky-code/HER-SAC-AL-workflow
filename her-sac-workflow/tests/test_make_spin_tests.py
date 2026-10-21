"""Tests for scripts/make_spin_tests.py (Task 6)."""
import argparse
import json
import os
import unittest

import helpers
import herlib
from make_spin_tests import setup, summarise


def args(**kw):
    a = dict(max=4, seed=7, dmin=0.02, dmax=0.04, nelect=None, force=False, script="sub.sh", log="log",
             no_scheduler=True, json=None)
    a.update(kw)
    return argparse.Namespace(**a)


def converged(nelect=130, s=None):
    s = s or helpers.structure([(0, 0, 1.6)])
    return helpers.relax_job(s, -200.0, mag=2.0, extra={"OUTCAR": f"   NELECT =     {nelect:.4f}    total number\n",
                                                        "KPOINTS": "k\n", "POTCAR": "p\n", "sub.sh": "#!/bin/bash\n"})


class SpinTests(unittest.TestCase):
    def test_even_parity_folders_and_incar(self):
        d = converged(130)
        made = setup(d, args())
        self.assertEqual([os.path.basename(x) for x in made], ["spin_N0", "spin_N2", "spin_N4"])
        orig = herlib.parse_incar(os.path.join(d, "INCAR"))
        new = herlib.parse_incar(os.path.join(d, "spin_N2", "INCAR"))
        self.assertEqual((new.pop("NUPDOWN"), new.pop("ISYM")), ("2", "0"))
        self.assertEqual(new, orig)
        for f in ("KPOINTS", "POTCAR", "sub.sh"):
            self.assertTrue(os.path.isfile(os.path.join(d, "spin_N4", f)))

    def test_odd_parity(self):
        d = converged(131)
        self.assertEqual([os.path.basename(x) for x in setup(d, args(max=5))], ["spin_N1", "spin_N3", "spin_N5"])

    def test_displacements_recorded_and_bounded(self):
        d = converged()
        setup(d, args())
        with open(os.path.join(d, "spin_tests.json")) as fh:
            man = json.load(fh)[0]
        self.assertEqual(man["seed"], 7)
        self.assertEqual(len(man["displacements"]), 5)
        self.assertEqual(man["displacements"][0]["atom"][:2], "Ir")
        a = herlib.read_structure(os.path.join(d, "CONTCAR"))
        b = herlib.read_structure(os.path.join(d, "spin_N0", "POSCAR"))
        moved = [i for i in range(len(a["elements"]))
                 if any(abs(x - y) > 1e-7 for x, y in zip(a["cart"][i], b["cart"][i]))]
        self.assertEqual(len(moved), 5)
        for i in moved:
            dist = sum((x - y) ** 2 for x, y in zip(a["cart"][i], b["cart"][i])) ** 0.5
            self.assertTrue(0.02 - 1e-6 <= dist <= 0.04 + 1e-6, dist)
        self.assertNotIn("H", [a["elements"][i] for i in moved])
        # same seed -> same displacements
        d2 = converged()
        setup(d2, args())
        with open(os.path.join(d2, "spin_tests.json")) as fh:
            self.assertEqual(json.load(fh)[0]["displacements"], man["displacements"])

    def test_no_overwrite(self):
        d = converged()
        setup(d, args())
        self.assertEqual(setup(d, args()), [])

    def test_d0_metal_skipped(self):
        s = helpers.structure()
        s["elements"] = [("Y" if e == "Ir" else e) for e in s["elements"]]
        d = converged(s=s)
        with self.assertRaises(SystemExit):
            setup(d, args())
        self.assertTrue(setup(d, args(force=True)))

    def test_unconverged_refused(self):
        d = converged()
        with open(os.path.join(d, "log"), "w") as fh:
            fh.write("running\n")
        with self.assertRaises(SystemExit):
            setup(d, args())

    def test_summarise(self):
        d = converged()
        setup(d, args(max=2))
        for k, e0 in ((0, -199.9), (2, -200.05)):
            sub = os.path.join(d, f"spin_N{k}")
            with open(os.path.join(sub, "POSCAR")) as fh:
                pos = fh.read()
            for name, text in (("CONTCAR", pos), ("log", helpers.CONVERGED_LOG), ("OSZICAR", helpers.oszicar(e0, k))):
                with open(os.path.join(sub, name), "w") as fh:
                    fh.write(text)
        rows = summarise(d, args())
        self.assertEqual([r["NUPDOWN"] for r in rows], [0, 2])
        self.assertAlmostEqual(rows[0]["dE"], 0.1)
        self.assertAlmostEqual(rows[1]["dE"], -0.05)
        self.assertEqual(len(rows[0]["bond_change"]), 4)
        self.assertTrue(all(abs(v) < 0.1 for v in rows[0]["bond_change"].values()))


if __name__ == "__main__":
    unittest.main()
