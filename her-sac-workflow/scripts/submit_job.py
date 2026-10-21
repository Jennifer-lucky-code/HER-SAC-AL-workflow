#!/usr/bin/env python3
"""Submit job folders, write each JOBID, register the jobs in loop_state.json and log 'submitted' events.

    python submit_job.py <job_dir> [<job_dir> ...] [--scheduler auto|pbs|slurm] [--script sub.sh]
                         [--state loop_state.json] [--event-log workflow_log.jsonl] [--force] [--dry-run]

Before submitting it checks that INCAR, KPOINTS, POSCAR, POTCAR and the job script exist and that the POTCAR
TITEL elements match POSCAR line 6 in order (references/vasp_settings.md); a folder that fails is not submitted.
Jobs the state lists as queued or running are not resubmitted, and flagged jobs are refused, unless --force.
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import herlib  # noqa: E402


def potcar_elements(path):
    with open(path, errors="replace") as fh:
        return [m.group(1).split("_")[0] for m in re.finditer(r"TITEL\s*=\s*\S+\s+(\S+)", fh.read())]


def preflight(d, script):
    missing = [f for f in ("INCAR", "KPOINTS", "POSCAR", "POTCAR", script) if not os.path.isfile(os.path.join(d, f))]
    if missing:
        return f"missing {', '.join(missing)}"
    with open(os.path.join(d, "POSCAR")) as fh:
        species = fh.read().splitlines()[5].split()
    pot = potcar_elements(os.path.join(d, "POTCAR"))
    if pot != species:
        return f"POTCAR order {' '.join(pot)} does not match POSCAR {' '.join(species)}"
    return ""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--scheduler", default="auto", choices=["auto", "pbs", "slurm"])
    ap.add_argument("--script", default="sub.sh", help="job script in each folder")
    ap.add_argument("--submit-cmd", help=argparse.SUPPRESS)  # testing: replaces qsub/sbatch
    ap.add_argument("--state", default=herlib.STATE_FILE)
    ap.add_argument("--event-log", default=herlib.LOG_FILE)
    ap.add_argument("--force", action="store_true", help="submit even if queued/running/flagged in the state")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    state = herlib.load_state(a.state)
    for d in a.dirs:
        key = herlib.job_key(d, a.state)
        entry = herlib.job_entry(state, key)
        if not a.force and entry.get("jobid") and entry.get("status") in ("queued", "running"):
            print(f"SKIP   {d}: already {entry['status']} (JOBID {entry['jobid']})")
            continue
        if not a.force and entry.get("status") == "flagged":
            print(f"SKIP   {d}: flagged ({entry.get('flag_reason', '')})")
            continue
        problem = preflight(d, a.script)
        if problem:
            print(f"REFUSE {d}: {problem}")
            continue
        if a.dry_run:
            print(f"READY  {d}")
            continue
        try:
            jobid = herlib.submit(d, a.scheduler, a.script, a.submit_cmd)
        except RuntimeError as e:
            print(f"FAILED {d}: {e}")
            continue
        entry.update(jobid=jobid, status="queued")
        herlib.log_event(a.event_log, "submitted", key, "", "", jobid, retry=entry.get("retries", 0))
        print(f"SUBMITTED {d}  JOBID {jobid}")
    if not a.dry_run:
        herlib.save_state(state, a.state)


if __name__ == "__main__":
    main()
