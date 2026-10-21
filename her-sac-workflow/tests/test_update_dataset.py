"""Tests for scripts/update_dataset.py (Task 3)."""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

import helpers
import herlib

SCRIPTS = os.path.join(helpers.ROOT, "scripts")
G_H2 = -6.8901
NAME = "Ru-PONS@G"


LIGAND_C = herlib.find_site(helpers.structure())[1][0]


def project(h2c_offsets=((0, 0, 1.6),), h2c_on=((LIGAND_C, (0, 0, 1.1)),)):
    work = tempfile.mkdtemp()
    cand = os.path.join(work, NAME)
    s1 = helpers.structure([(0, 0, 1.6)])
    s2c = helpers.structure(h2c_offsets, h_on=h2c_on)
    helpers.relax_job(helpers.structure(), -200.0, mag=2.0, parent=cand, name="slab")
    for state, s, e, zpe in (("H1_metal", s1, -203.8, 0.2), ("H2_C", s2c, -207.5, 0.4),
                             ("H2_M", helpers.structure([(0.6, 0, 1.5), (-0.6, 0, 1.5)]), -207.4, 0.4)):
        helpers.relax_job(s, e, mag=1.0, parent=cand, name=state)
        helpers.freq_job(s, zpe, 0.01, parent=cand, name="freq_" + state)
    run(work, "collect_energies.py", NAME, "--no-scheduler")
    run(work, "her_analysis.py", "energies.csv")
    return work, cand


def run(work, script, *args, check=True):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args], cwd=work, check=check,
                          capture_output=True, text=True)


def dataset(work):
    with open(os.path.join(work, "dft_dataset.csv"), newline="") as fh:
        return list(csv.DictReader(fh))


class UpdateDataset(unittest.TestCase):
    def test_add_then_unchanged_then_supersede(self):
        work, cand = project()
        herlib.save_state({"jobs": {f"{NAME}/H2_C": {"retries": 2, "status": "done"}}},
                          os.path.join(work, "loop_state.json"))
        run(work, "update_dataset.py", "her_results.csv", "--no-scheduler", "--date", "2026-01-01")
        rows = dataset(work)
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual((r["canonical"], r["TM"], r["ligands"]), ("Ru-NSPO@G", "Ru", "PONS"))
        self.assertEqual((r["H_site"], r["retries"], r["mag_slab"], r["incar_diff"]), ("H1_metal", "2", "2.000", ""))
        self.assertEqual(r["eta_screen"] != "", True)
        self.assertIn(f"{NAME}/freq_H2_M", r["job_paths"].split(";"))
        log = herlib.read_log(os.path.join(work, "workflow_log.jsonl"))
        self.assertEqual([e["event"] for e in log], ["added"])
        self.assertEqual(herlib.load_state(os.path.join(work, "loop_state.json"))["added"], [NAME])

        out = run(work, "update_dataset.py", "her_results.csv", "--no-scheduler").stdout
        self.assertIn("unchanged", out)
        self.assertEqual(len(dataset(work)), 1)

        # same canonical key, different value -> new row, old row kept
        with open(os.path.join(work, "her_results.csv"), newline="") as fh:
            res = list(csv.DictReader(fh))
        res[0]["name"] = "Ru-SNOP@G"          # symmetry-equivalent ordering
        os.rename(cand, os.path.join(work, "Ru-SNOP@G"))
        res[0]["eta_VH"] = "0.9999"
        with open(os.path.join(work, "her_results.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(res[0]))
            w.writeheader()
            w.writerows(res)
        run(work, "update_dataset.py", "her_results.csv", "--no-scheduler")
        rows = dataset(work)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["eta_VH"], r["eta_VH"])
        self.assertIn("supersedes row 1", rows[1]["notes"])

    def test_structure_flag_refused_unless_allowed(self):
        work, _ = project(h2c_offsets=((0, 0, 1.6), (0.7, 0, 1.6)), h2c_on=())     # H–H 0.7 Å in H2_C
        out = run(work, "update_dataset.py", "her_results.csv", "--no-scheduler").stdout
        self.assertIn("refused", out)
        self.assertIn("dihydrogen", out)
        self.assertFalse(os.path.exists(os.path.join(work, "dft_dataset.csv")))
        run(work, "update_dataset.py", "her_results.csv", "--no-scheduler", "--allow-flagged")
        self.assertIn("dihydrogen", dataset(work)[0]["notes"])

    def test_incar_diff_recorded(self):
        work, cand = project()
        herlib.set_incar(os.path.join(cand, "H2_C", "INCAR"), {"IBRION": 1, "NPAR": 4})
        run(work, "update_dataset.py", "her_results.csv", "--no-scheduler")
        self.assertEqual(dataset(work)[0]["incar_diff"], "H2_C:IBRION=1 (template 2)")

    def test_unconverged_refused_even_if_allowed(self):
        work, cand = project()
        with open(os.path.join(cand, "H2_M", "log"), "w") as fh:
            fh.write("killed\n")
        out = run(work, "update_dataset.py", "her_results.csv", "--no-scheduler", "--allow-flagged").stdout
        self.assertIn("refused", out)
        self.assertIn("H2_M: job not converged", out)

    def test_extra_columns(self):
        work, _ = project()
        with open(os.path.join(work, "features.csv"), "w") as fh:
            fh.write(f"name,Ef,d_band_centre\n{NAME},-1.23,-2.5\n")
        run(work, "update_dataset.py", "her_results.csv", "--no-scheduler", "--extra", "features.csv")
        r = dataset(work)[0]
        self.assertEqual((r["Ef"], r["d_band_centre"]), ("-1.2300", "-2.5000"))

    def test_foreign_header_not_touched(self):
        work, _ = project()
        with open(os.path.join(work, "dft_dataset.csv"), "w") as fh:
            fh.write("a,b\n1,2\n")
        p = run(work, "update_dataset.py", "her_results.csv", "--no-scheduler", check=False)
        self.assertNotEqual(p.returncode, 0)
        with open(os.path.join(work, "dft_dataset.csv")) as fh:
            self.assertEqual(fh.read(), "a,b\n1,2\n")


if __name__ == "__main__":
    unittest.main()
