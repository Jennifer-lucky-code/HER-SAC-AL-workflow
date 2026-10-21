"""Tests for scripts/collect_energies.py (Task 2), including an end-to-end run through her_analysis.py."""
import csv
import os
import subprocess
import sys
import tempfile
import unittest

import helpers
import herlib
from collect_energies import collect

SCRIPTS = os.path.join(helpers.ROOT, "scripts")
G_H2 = -6.8901


def candidate(states, name="Cr-PSSS@G"):
    """states: {folder: (E0, ZPE, TS) or (E0,) for slab}; freq folders are made when ZPE/TS are given."""
    root = tempfile.mkdtemp()
    cand = os.path.join(root, name)
    s = helpers.structure([(0, 0, 1.6)])
    for state, vals in states.items():
        helpers.relax_job(s, vals[0], parent=cand, name=state)
        if len(vals) == 3:
            helpers.freq_job(s, vals[1], vals[2], parent=cand, name="freq_" + state)
    return cand


class Settings(unittest.TestCase):
    def test_g_h2_single_source(self):
        self.assertEqual(herlib.g_h2(), G_H2)

    def test_vaspkit_lines_from_vasp_settings(self):
        text = ("Zero-point energy E_ZPE   :      13.814 kcal/mol    0.599026 eV\n"
                "Entropy contribution T*S  :   15147.343 J/(mol)     0.156991 eV\n")
        self.assertEqual(herlib.parse_vaspkit_501(text), {"ZPE": 0.599026, "TS": 0.156991})


class Collect(unittest.TestCase):
    def test_lowest_g_site_and_missing_config(self):
        cand = candidate({"slab": (-200.0,), "H1_metal": (-203.9, 0.20, 0.01), "H1_lig2": (-204.1, 0.30, 0.01),
                          "H2_C": (-207.5, 0.40, 0.02)})
        row, detail = collect(cand, G_H2, use_scheduler=False)
        self.assertEqual(row["H_site"], "H1_lig2")               
        self.assertAlmostEqual(row["G_H"], -204.1 + 0.30 - 0.01)
        self.assertAlmostEqual(row["G_2H_C"], -207.5 + 0.40 - 0.02)
        self.assertIsNone(row["G_2H_M"])
        self.assertIn("H2_M: no folder", row["notes"])
        self.assertEqual(row["G_H2"], G_H2)

    def test_unconverged_is_blank_not_estimated(self):
        cand = candidate({"slab": (-200.0,), "H1_metal": (-203.9, 0.20, 0.01)})
        with open(os.path.join(cand, "slab", "log"), "w") as fh:
            fh.write("still running\n")
        row, _ = collect(cand, G_H2, use_scheduler=False)
        self.assertIsNone(row["E_slab"])
        self.assertIn("slab: unconverged", row["notes"])

    def test_missing_freq_gives_note_and_site_warning(self):
        cand = candidate({"slab": (-200.0,), "H1_metal": (-203.9, 0.20, 0.01), "H1_lig1": (-205.0,)})
        row, _ = collect(cand, G_H2, use_scheduler=False)
        self.assertEqual(row["H_site"], "H1_metal")
        self.assertIn("freq_H1_lig1: no folder", row["notes"])
        self.assertIn("H1_lig1 has a lower E than H1_metal", row["notes"])


class EndToEnd(unittest.TestCase):
    def test_cr_psss_worked_check(self):
        # E values chosen so that dG_H* = -0.352 and 2H* levels -0.397 (C), -0.401 (M) with G_H2 = -6.8901
        E_s = -200.0
        G_H = E_s + 0.5 * G_H2 - 0.352
        G_C, G_M = E_s + G_H2 - 0.397, E_s + G_H2 - 0.401
        cand = candidate({"slab": (E_s,), "H1_metal": (G_H - 0.2 + 0.01, 0.2, 0.01),
                          "H2_C": (G_C - 0.4 + 0.02, 0.4, 0.02), "H2_M": (G_M - 0.4 + 0.02, 0.4, 0.02)})
        work = tempfile.mkdtemp()
        run = lambda *args: subprocess.run([sys.executable, *args], cwd=work, check=True,
                                           capture_output=True, text=True)
        run(os.path.join(SCRIPTS, "collect_energies.py"), cand, "--no-scheduler")
        run(os.path.join(SCRIPTS, "her_analysis.py"), "energies.csv")
        with open(os.path.join(work, "her_results.csv")) as fh:
            r = next(csv.DictReader(fh))
        self.assertEqual((r["eta_VH"], r["eta_VTC"], r["eta_VTM"]), ("0.3520", "0.3970", "0.4010"))
        self.assertEqual(r["pathway"], "VH")
        self.assertEqual(r["H_site"], "H1_metal")


if __name__ == "__main__":
    unittest.main()
