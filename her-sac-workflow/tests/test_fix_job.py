"""Tests for scripts/fix_job.py and the structured actions in check_jobs.ERRORS (Task 4)."""
import os
import tempfile
import unittest

import helpers
import herlib
from check_jobs import ERRORS
from fix_job import PROTECTED_INCAR, fix

SUB = "#!/bin/bash\n#PBS -l ncpus=48\n#PBS -l mem=190GB\nmpirun vasp_std > log\n"


def failed_job(log_text, **extra):
    s = helpers.structure()
    moved = helpers.structure(metal_shift=(0, 0, 0.05))
    files = {"INCAR": helpers.RELAX_INCAR, "POSCAR": helpers.poscar_text(s), "CONTCAR": helpers.poscar_text(moved),
             "log": log_text, "sub.sh": SUB, "WAVECAR": "x", "JOBID": "111.gadi-pbs\n"}
    files.update(extra)
    return helpers.make_dir(files)


class Fix(unittest.TestCase):
    def setUp(self):
        root = tempfile.mkdtemp()
        self.state_path = os.path.join(root, "loop_state.json")
        self.opts = {"log": "log", "use_scheduler": False, "state": self.state_path,
                     "event_log": os.path.join(root, "workflow_log.jsonl"), "max_retries": 3, "submit": False,
                     "scheduler": "auto", "script": "sub.sh", "submit_cmd": None, "dry_run": False}
        self.state = herlib.load_state(self.state_path)

    def read(self, d, name):
        with open(os.path.join(d, name)) as fh:
            return fh.read()

    def events(self):
        return [(e["event"], e["error_type"]) for e in herlib.read_log(self.opts["event_log"])]

    def entry(self, d):
        return self.state["jobs"][herlib.job_key(d, self.state_path)]

    def test_zbrent_then_else_then_submit(self):
        d = failed_job("ZBRENT: fatal error in bracketing\n")
        r = fix(d, dict(self.opts, submit=True, submit_cmd="echo 4242.gadi-pbs"), self.state)
        self.assertEqual(r["result"], "fixed and resubmitted")
        self.assertEqual(herlib.parse_incar(os.path.join(d, "INCAR"))["IBRION"], "1")
        self.assertEqual(herlib.parse_incar(os.path.join(d, "INCAR-1"))["IBRION"], "2")
        self.assertIn("ZBRENT", self.read(d, "log-1"))
        self.assertFalse(os.path.exists(os.path.join(d, "log")))
        self.assertEqual(self.read(d, "POSCAR"), self.read(d, "CONTCAR"))
        self.assertEqual(self.read(d, "JOBID").strip(), "4242.gadi-pbs")
        self.assertEqual((self.entry(d)["retries"], self.entry(d)["status"]), (1, "queued"))
        self.assertEqual(self.events(), [("fixed", "zbrent"), ("submitted", "zbrent")])

        with open(os.path.join(d, "log"), "w") as fh:
            fh.write("ZBRENT: fatal error in bracketing\n")
        r = fix(d, self.opts, self.state)
        inc = herlib.parse_incar(os.path.join(d, "INCAR"))
        self.assertEqual((inc["IBRION"], inc["POTIM"]), ("1", "0.2"))
        self.assertEqual(r["result"], "fixed (resubmit with --submit)")
        self.assertEqual(self.entry(d)["status"], "ready")
        self.assertTrue(os.path.exists(os.path.join(d, "log-2")))

    def test_retry_limit_flags(self):
        d = failed_job("ZBRENT: fatal error in bracketing\n")
        herlib.job_entry(self.state, herlib.job_key(d, self.state_path))["retries"] = 3
        r = fix(d, self.opts, self.state)
        self.assertEqual(r["result"], "flagged")
        self.assertIn("retry limit", r["change"])
        self.assertEqual(self.entry(d)["status"], "flagged")
        self.assertFalse(os.path.exists(os.path.join(d, "log-4")))
        # a flagged job is skipped afterwards
        self.assertTrue(fix(d, self.opts, self.state)["result"].startswith("skipped"))

    def test_crash_never_auto_fixed(self):
        d = failed_job("forrtl: severe (174): SIGSEGV, segmentation fault occurred\n")
        before = self.read(d, "INCAR")
        r = fix(d, self.opts, self.state)
        self.assertEqual(r["result"], "flagged")
        self.assertIn("needs judgement", r["change"])
        self.assertEqual(self.read(d, "INCAR"), before)
        self.assertFalse(os.path.exists(os.path.join(d, "log-1")))
        self.assertEqual(self.events(), [("flagged", "crash")])

    def test_unknown_never_auto_fixed(self):
        d = helpers.make_dir({"INCAR": helpers.RELAX_INCAR})
        self.assertEqual(fix(d, self.opts, self.state)["result"], "flagged")

    def test_potcar_protected(self):
        d = failed_job("ERROR: number of potentials on File POTCAR incompatible with number\n")
        r = fix(d, self.opts, self.state)
        self.assertEqual(r["result"], "flagged")
        self.assertIn("POTCAR", r["change"])

    def test_memory_scaled(self):
        d = failed_job("=>> PBS: job killed: mem 200GB exceeded limit 190GB\n", log="")
        with open(os.path.join(d, "gadi.o111"), "w") as fh:
            fh.write("Out Of Memory\n")
        r = fix(d, self.opts, self.state)
        self.assertEqual(r["error"], "memory")
        self.assertIn("#PBS -l mem=285GB", self.read(d, "sub.sh"))

    def test_nsw_reached_restarts_from_contcar(self):
        d = failed_job("x\n", OSZICAR=helpers.oszicar(-200.0, step=500))
        r = fix(d, self.opts, self.state)
        self.assertEqual(r["change"], "cp CONTCAR POSCAR")
        self.assertEqual(self.read(d, "POSCAR"), self.read(d, "CONTCAR"))

    def test_brmix_recurrence_uses_else(self):
        d = failed_job("BRMIX: very serious problems\n")
        self.assertEqual(fix(d, self.opts, self.state)["change"], "rm WAVECAR")
        # not submitted yet -> a second run must not fix it again
        self.assertEqual(fix(d, self.opts, self.state)["result"], "skipped: fixed, waiting for submission")
        self.entry(d)["status"] = "queued"
        with open(os.path.join(d, "log"), "w") as fh:
            fh.write("BRMIX: very serious problems\n")
        self.assertEqual(fix(d, self.opts, self.state)["change"], "ISYM=0")

    def test_cell_change_refused(self):
        s = helpers.structure()
        s["lattice"][0][0] += 0.2
        d = failed_job("BRIONS problems: POTIM should be increased\n", CONTCAR=helpers.poscar_text(s))
        r = fix(d, self.opts, self.state)
        self.assertEqual(r["result"], "flagged")
        self.assertIn("cell", r["change"])

    def test_dry_run_changes_nothing(self):
        d = failed_job("ZBRENT: fatal error in bracketing\n")
        r = fix(d, dict(self.opts, dry_run=True), self.state)
        self.assertEqual(r["result"], "would fix")
        self.assertEqual(r["change"], "IBRION=1; cp CONTCAR POSCAR")
        self.assertTrue(os.path.exists(os.path.join(d, "log")))
        self.assertEqual(self.events(), [])

    def test_converged_left_alone(self):
        d = failed_job(helpers.CONVERGED_LOG)
        self.assertEqual(fix(d, self.opts, self.state)["result"], "nothing to fix")


class Table(unittest.TestCase):
    def test_every_error_has_a_structured_action(self):
        for pat, name, text, action in ERRORS:
            self.assertIsInstance(action, dict, name)
            for a in (action, action.get("else", {})):
                self.assertFalse(set(a.get("incar", {})) & PROTECTED_INCAR, name)


if __name__ == "__main__":
    unittest.main()
