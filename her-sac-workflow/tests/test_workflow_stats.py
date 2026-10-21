"""Tests for the structured log/state (Task 5): workflow_stats.py, check_jobs --state/--event-log, submit_job.py."""
import json
import os
import subprocess
import sys
import tempfile
import unittest

import helpers
import herlib
from workflow_stats import stats

SCRIPTS = os.path.join(helpers.ROOT, "scripts")

LOG = [
    ("A", "submitted", ""), ("A", "checked", ""), ("A", "converged", ""),
    ("B", "submitted", ""), ("B", "fixed", "zbrent"), ("B", "submitted", "zbrent"), ("B", "converged", ""),
    ("C", "submitted", ""), ("C", "fixed", "zbrent"), ("C", "fixed", "brmix"), ("C", "flagged", "brmix"),
    ("D", "submitted", ""),
    ("cand", "added", ""),
]
STATE = {"jobs": {"A": {"status": "done", "retries": 0}, "B": {"status": "done", "retries": 1},
                  "C": {"status": "flagged", "retries": 2}, "D": {"status": "queued", "retries": 0}}}


def run(cwd, script, *args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args], cwd=cwd, check=True,
                          capture_output=True, text=True).stdout


class Stats(unittest.TestCase):
    def test_synthetic_log(self):
        events = [{"job_dir": j, "event": e, "error_type": t} for j, e, t in LOG]
        s = stats(STATE, events)
        self.assertEqual(s["total_jobs"], 4)
        self.assertEqual(s["converged"], 2)
        self.assertEqual(s["converged_first_attempt"], 1)
        self.assertEqual(s["auto_fixed_jobs"], 2)
        self.assertEqual(s["fixes_per_error_type"], {"brmix": 1, "zbrent": 2})
        self.assertEqual(s["jobs_fixed_per_error_type"], {"brmix": 1, "zbrent": 2})
        self.assertEqual(s["total_retries"], 3)
        self.assertEqual((s["flagged_now"], s["flagged_ever"]), (1, 1))
        self.assertEqual(s["running_or_queued"], 1)
        self.assertEqual(s["pct_without_manual_of_all"], 50.0)
        self.assertEqual(s["pct_without_manual_of_finished"], 66.7)
        self.assertEqual(s["rows_added"], 1)

    def test_cli_on_files(self):
        work = tempfile.mkdtemp()
        herlib.save_state(STATE, os.path.join(work, "loop_state.json"))
        log = os.path.join(work, "workflow_log.jsonl")
        for j, e, t in LOG:
            herlib.log_event(log, e, j, t)
        out = run(work, "workflow_stats.py", "--json", "s.json")
        self.assertIn("50.0 % of all jobs", out)
        with open(os.path.join(work, "s.json")) as fh:
            self.assertEqual(json.load(fh)["total_retries"], 3)

    def test_unknown_event_rejected(self):
        with self.assertRaises(ValueError):
            herlib.log_event(os.path.join(tempfile.mkdtemp(), "l.jsonl"), "guessed")


class Record(unittest.TestCase):
    def test_check_jobs_records_converged_once(self):
        work = tempfile.mkdtemp()
        helpers.make_dir({"INCAR": helpers.RELAX_INCAR, "log": helpers.CONVERGED_LOG}, work, "job1")
        helpers.make_dir({"INCAR": helpers.RELAX_INCAR, "log": "ZBRENT: fatal error\n"}, work, "job2")
        for _ in range(2):
            run(work, "check_jobs.py", "job1", "job2", "--no-scheduler", "--state", "loop_state.json",
                "--event-log", "workflow_log.jsonl")
        ev = [(e["job_dir"], e["event"]) for e in herlib.read_log(os.path.join(work, "workflow_log.jsonl"))]
        self.assertEqual(ev.count(("job1", "converged")), 1)
        self.assertEqual(ev.count(("job2", "checked")), 2)
        st = herlib.load_state(os.path.join(work, "loop_state.json"))["jobs"]
        self.assertEqual((st["job1"]["status"], st["job2"]["status"]), ("done", "failed"))


class Submit(unittest.TestCase):
    def job(self, potcar_elems):
        work = tempfile.mkdtemp()
        potcar = "".join(f"   TITEL  = PAW_PBE {e} 08Apr2002\n" for e in potcar_elems)
        helpers.make_dir({"INCAR": helpers.RELAX_INCAR, "KPOINTS": "k\n", "sub.sh": "#!/bin/bash\n",
                          "POSCAR": helpers.poscar_text(helpers.structure()), "POTCAR": potcar}, work, "slab")
        return work

    def test_submit_registers_and_logs(self):
        work = self.job(["C", "N", "S", "Ir"])
        out = run(work, "submit_job.py", "slab", "--submit-cmd", "echo 777.gadi-pbs")
        self.assertIn("SUBMITTED", out)
        st = herlib.load_state(os.path.join(work, "loop_state.json"))["jobs"]["slab"]
        self.assertEqual((st["jobid"], st["status"]), ("777.gadi-pbs", "queued"))
        self.assertEqual([e["event"] for e in herlib.read_log(os.path.join(work, "workflow_log.jsonl"))],
                         ["submitted"])
        self.assertIn("SKIP", run(work, "submit_job.py", "slab", "--submit-cmd", "echo 778"))

    def test_potcar_order_mismatch_refused(self):
        work = self.job(["C", "S", "N", "Ir"])
        self.assertIn("REFUSE", run(work, "submit_job.py", "slab", "--submit-cmd", "echo 1"))
        self.assertFalse(os.path.exists(os.path.join(work, "slab", "JOBID")))


if __name__ == "__main__":
    unittest.main()
