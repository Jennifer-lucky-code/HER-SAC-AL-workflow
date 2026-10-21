#!/usr/bin/env python3
"""Classify VASP job folders: converged, still running/queued, failed (with the error), or stopped unconverged.

Run on the HPC login node every 5 h while a cycle has jobs out (see SKILL.md step 6):

    python check_jobs.py <job_dir> [<job_dir> ...] [--log log] [--no-scheduler] [--json status.json]
                         [--state loop_state.json --event-log workflow_log.jsonl]

Per folder it reads only small pieces: the tail of the job log, the tail of OSZICAR/OUTCAR, the scheduler
output file, INCAR (IBRION/NSW) and JOBID (scheduler ID written at submission: `qsub sub.sh > JOBID`).

Status:
    converged    relaxation: 'reached required accuracy' in the last two log lines;
                 frequency (IBRION 5-8) / static (NSW 0): OUTCAR finished normally
    running      the scheduler still lists the job (queued, held or running)
    failed       a known VASP or scheduler error was found -> 'error' and 'fix' say what to change
    unconverged  job ended without convergence and without a known error (e.g. NSW reached)
    unknown      nothing to judge from (no log, no JOBID); report to the user
With --state/--event-log it also records the check: one 'checked' event per job, a 'converged' event when a job
first converges, and the job status in loop_state.json (converged -> done, running -> running,
otherwise failed; 'flagged' and 'ready' (fixed, waiting for submission) are kept).
Each failed/unconverged result also carries 'action', the machine-readable fix (see ERRORS below).
The script only reports; scripts/fix_job.py applies the action following references/dft_her.md section 6.
"""
import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import herlib  # noqa: E402

CONVERGED = "reached required accuracy - stopping structural energy minimisation"

# (pattern, short name, fix, action). First match wins, so specific patterns come first.
# 'action' is the machine-readable form of 'fix', applied by fix_job.py (one action per retry):
#   incar              {KEY: value} set in INCAR (None removes the key)
#   contcar_to_poscar  cp CONTCAR POSCAR before resubmitting (relaxations only)
#   remove             files deleted before resubmitting
#   sub                {"mem_factor": f} scale the memory request in the job script
#   nbands_factor      NBANDS -> ceil(f x current NBANDS)
#   magmom_fit         rewrite a uniform MAGMOM to the POSCAR atom count
#   protected          the fix would touch an energy-reference setting -> flag, never change
#   manual             needs judgement -> never auto-fixed, reported and flagged
#   else               the action to use instead when the INCAR part of the action is already in place,
#                      or (actions without an INCAR part) when the same error comes back
ERRORS = [
    (r"walltime .* exceeded limit|DUE TO TIME LIMIT", "walltime",
     "cp CONTCAR POSCAR and resubmit (continuation, not a settings change); ask for longer walltime if repeated",
     {"contcar_to_poscar": True, "else": {"manual": "walltime exceeded again: ask for a longer walltime"}}),
    (r"exceeded memory|Out Of Memory|oom-kill|killed by signal 9", "memory",
     "resubmit with more memory (mem / nodes); settings unchanged",
     {"sub": {"mem_factor": 1.5}}),
    (r"number of potentials on File POTCAR incompatible|POTCAR.*not.*match", "potcar",
     "regenerate POTCAR with vaspkit 102 and check TITEL order against POSCAR line 6",
     {"protected": "POTCAR"}),
    (r"Error reading item 'MAGMOM'|MAGMOM.*number of items", "magmom",
     "rewrite MAGMOM so its count equals the POSCAR atom count",
     {"magmom_fit": True}),
    (r"ZBRENT: fatal error", "zbrent",
     "cp CONTCAR POSCAR; set IBRION = 1 (if it recurs: POTIM = 0.2)",
     {"incar": {"IBRION": 1}, "contcar_to_poscar": True,
      "else": {"incar": {"POTIM": 0.2}, "contcar_to_poscar": True}}),
    (r"BRIONS problems: POTIM should be increased", "brions",
     "cp CONTCAR POSCAR and resubmit",
     {"contcar_to_poscar": True}),
    (r"RSPHER|REAL_OPTLAY|REAL_OPT: internal ERROR", "lreal",
     "set LREAL = .FALSE.",
     {"incar": {"LREAL": ".FALSE."}}),
    (r"PRICEL|IBZKPT|inverse of rotation matrix was not found|SGRCON|POSMAP", "symmetry",
     "set ISYM = 0",
     {"incar": {"ISYM": 0}}),
    (r"Call to ZHEGV failed|ZHEEV|Sub-Space-Matrix is not hermitian|EDDDAV|EDWAV", "eigensolver",
     "set ALGO = All (if already All: ALGO = Normal); remove WAVECAR",
     {"incar": {"ALGO": "All"}, "remove": ["WAVECAR"],
      "else": {"incar": {"ALGO": "Normal"}, "remove": ["WAVECAR"]}}),
    (r"BRMIX: very serious problems", "brmix",
     "remove WAVECAR/CHGCAR and resubmit; if it recurs set ISYM = 0",
     {"remove": ["WAVECAR", "CHGCAR"], "else": {"incar": {"ISYM": 0}}}),
    (r"TOO FEW BANDS|NBANDS.*too small", "nbands",
     "increase NBANDS by ~20 %",
     {"nbands_factor": 1.2}),
    (r"segmentation fault|forrtl: severe|BAD TERMINATION|MPI_ABORT", "crash",
     "check the log lines just above the crash; if no VASP error is shown, resubmit once unchanged",
     {"manual": "crash: read the log lines above it"}),
]


def tail(path, n=60):
    if not os.path.isfile(path):
        return []
    with open(path, "rb") as fh:
        fh.seek(0, 2)
        fh.seek(max(0, fh.tell() - 200 * n))
        return fh.read().decode("utf-8", "replace").splitlines()[-n:]


def incar(d):
    vals = {}
    for line in tail(os.path.join(d, "INCAR"), 500):
        m = re.match(r"\s*([A-Za-z_]+)\s*=\s*([^\s#!(]+)", line)
        if m:
            vals[m.group(1).upper()] = m.group(2)
    return vals


def scheduler_state(jobid):
    """'Q'/'R'/'H'/... while the scheduler lists the job, None if it is gone or unknown."""
    if not jobid:
        return None
    if shutil.which("qstat"):
        r = subprocess.run(["qstat", "-f", jobid], capture_output=True, text=True)
        m = re.search(r"job_state\s*=\s*(\w)", r.stdout)
        return m.group(1) if r.returncode == 0 and m and m.group(1) not in "FX" else None
    if shutil.which("squeue"):
        r = subprocess.run(["squeue", "-h", "-j", jobid, "-o", "%t"], capture_output=True, text=True)
        s = r.stdout.strip()
        return s or None
    return None


def check(d, logname="log", use_scheduler=True):
    inc = incar(d)
    ibrion, nsw = int(inc.get("IBRION", 0)), int(inc.get("NSW", 0))
    jobid = ""
    if os.path.isfile(os.path.join(d, "JOBID")):
        with open(os.path.join(d, "JOBID")) as fh:
            jobid = fh.read().strip().split(".")[0]
    log = tail(os.path.join(d, logname))
    sched_out = []
    for pat in ("*.o*", "*.e*", "slurm-*.out"):
        for f in glob.glob(os.path.join(d, pat)):
            if not jobid or jobid in f:
                sched_out += tail(f, 40)
    res = {"dir": d, "jobid": jobid, "status": "unknown", "error": "", "fix": "", "detail": "",
           "action": {"manual": "nothing to judge from (no log, no JOBID)"}}

    if 5 <= ibrion <= 8 or nsw == 0:
        outcar = tail(os.path.join(d, "OUTCAR"), 30)
        if any("General timing and accounting" in l for l in outcar):
            res.update(status="converged", action={})
            return res
    elif any(CONVERGED in l for l in log[-2:]):
        res.update(status="converged", action={})
        return res

    state = scheduler_state(jobid) if use_scheduler else None
    if state:
        res.update(status="running", detail=f"scheduler state {state}", action={})
        return res

    text = "\n".join(log + sched_out)
    for pat, name, fix, action in ERRORS:
        m = re.search(pat, text, re.I)
        if m:
            line = next((l for l in text.splitlines() if re.search(pat, l, re.I)), m.group(0))
            res.update(status="failed", error=name, fix=fix, detail=line.strip()[:160], action=action)
            return res

    osz = tail(os.path.join(d, "OSZICAR"), 5)
    steps = [l for l in osz if re.match(r"\s*\d+\s+F=", l)]
    if steps:
        last = int(steps[-1].split()[0])
        nsw_reached = bool(nsw and last >= nsw)
        res.update(status="unconverged", detail=f"{last} ionic steps (NSW = {nsw})",
                   fix="cp CONTCAR POSCAR and resubmit" if nsw_reached
                   else "job ended early without a known error: read the log tail and the scheduler output",
                   action={"contcar_to_poscar": True} if nsw_reached
                   else {"manual": "ended early without a known error"})
    elif log or sched_out:
        res.update(status="unconverged", detail="no ionic step finished",
                   fix="read the log tail and the scheduler output",
                   action={"manual": "no ionic step finished"})
    return res


def record(results, state_path, event_log):
    """Write job statuses to loop_state.json and 'checked'/'converged' events to workflow_log.jsonl."""
    state = herlib.load_state(state_path)
    for r in results:
        key = herlib.job_key(r["dir"], state_path)
        entry = herlib.job_entry(state, key)
        prev = entry.get("status")
        if r["status"] == "converged":
            entry["status"] = "done"
        elif r["status"] == "running":
            entry["status"] = "running"
        elif prev not in ("flagged", "ready", "done"):
            entry["status"] = "failed"
        if r["jobid"]:
            entry["jobid"] = entry.get("jobid") or r["jobid"]
        if event_log:
            herlib.log_event(event_log, "checked", key, r["error"], r["status"], r["jobid"])
            if r["status"] == "converged" and prev != "done":
                herlib.log_event(event_log, "converged", key, "", "", r["jobid"],
                                 retries=entry.get("retries", 0))
    herlib.save_state(state, state_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--log", default="log", help="VASP stdout file written by the job script")
    ap.add_argument("--no-scheduler", action="store_true", help="do not query qstat/squeue")
    ap.add_argument("--json", help="also write the results to this file")
    ap.add_argument("--state", help="record job statuses in this loop_state.json")
    ap.add_argument("--event-log", help="append 'checked'/'converged' events to this workflow_log.jsonl")
    a = ap.parse_args()
    results = [check(d, a.log, not a.no_scheduler) for d in a.dirs]
    if a.state or a.event_log:
        record(results, a.state or herlib.STATE_FILE, a.event_log)
    for r in results:
        extra = f"  [{r['error']}] {r['detail']}  -> {r['fix']}" if r["status"] == "failed" else \
            f"  {r['detail']}" + (f"  -> {r['fix']}" if r["fix"] else "")
        print(f"{r['status']:<12} {r['dir']}{extra}")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(results, fh, indent=1)


if __name__ == "__main__":
    main()
