#!/usr/bin/env python3
"""Apply the fix that check_jobs.py names for a failed job, exactly as in references/dft_her.md §6.

    python fix_job.py <job_dir> [<job_dir> ...] [--submit] [--scheduler auto|pbs|slurm] [--script sub.sh]
                      [--max-retries 3] [--state loop_state.json] [--event-log workflow_log.jsonl]
                      [--log log] [--no-scheduler] [--dry-run] [--json fixes.json]

Per failed or unconverged job (converged and running jobs are left alone):
    1. back up: log -> log-N (moved), INCAR -> INCAR-N (copied), N = retry number
    2. apply exactly one action, the structured 'action' of check_jobs.ERRORS ('else' when the first one is
       already in place or the same error came back)
    3. cp CONTCAR POSCAR when the action says so (relaxations only; refused if CONTCAR has another cell)
    4. resubmit only with --submit, and 5. write the new JOBID
    6. retries + 1 in loop_state.json, one 'fixed' (and 'submitted') event in workflow_log.jsonl
A job is marked 'flagged' (event 'flagged') and not touched when:
    - it already has --max-retries retries;
    - the action needs judgement ('manual': crash, unknown, ...);
    - the fix would change a protected setting: ENCUT, PREC, IVDW, ISPIN, POTCAR, KPOINTS or the cell;
    - no fix is left (the action and its 'else' are both already in place).
Jobs that are already flagged are skipped until the user resets their status; jobs fixed without --submit
('ready') are skipped until submit_job.py submits them.
"""
import argparse
import json
import math
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import herlib  # noqa: E402
from check_jobs import check  # noqa: E402

PROTECTED_INCAR = {"ENCUT", "PREC", "IVDW", "ISPIN"}
PROTECTED_FILES = {"POTCAR", "KPOINTS"}


class Flag(Exception):
    """The job cannot be fixed automatically; the message is the reason."""


def incar_in_place(incar, changes):
    return all(k in incar and herlib.norm_value(incar[k]) == herlib.norm_value(v) for k, v in changes.items())


def choose(action, incar, entry, error):
    """The single action to apply now (primary or its 'else')."""
    primary = {k: v for k, v in action.items() if k != "else"}
    if "incar" in primary:
        use_else = incar_in_place(incar, primary["incar"])
    else:
        use_else = bool(entry.get("last_error")) and entry.get("last_error") == error
    if not use_else:
        return primary
    if "else" not in action:
        if "incar" in primary:
            raise Flag("fix already in place and no further fix in the table")
        return primary
    alt = action["else"]
    if "incar" in alt and incar_in_place(incar, alt["incar"]):
        raise Flag("fix and its alternative are both already in place")
    return alt


def magmom_fit(d, incar):
    poscar = herlib.read_structure(os.path.join(d, "POSCAR"))
    n = len(poscar["elements"])
    vals = []
    for tok in incar.get("MAGMOM", "").split():
        if "*" in tok:
            k, v = tok.split("*", 1)
            vals += [v] * int(k)
        else:
            vals.append(tok)
    if not vals or len({herlib.norm_value(v) for v in vals}) != 1:
        raise Flag("MAGMOM is not uniform; rewriting it needs judgement")
    return f"{n}*{vals[0]}"


def current_nbands(d, incar):
    if "NBANDS" in incar:
        return int(incar["NBANDS"])
    path = os.path.join(d, "OUTCAR")
    if os.path.isfile(path):
        with open(path, errors="replace") as fh:
            for line in fh:
                m = re.search(r"NBANDS\s*=\s*(\d+)", line)
                if m:
                    return int(m.group(1))
    raise Flag("NBANDS not found in INCAR or OUTCAR")


def scale_memory(script_path, factor):
    """New text of the job script with its memory request scaled (PBS -l mem= or SBATCH --mem=)."""
    if not os.path.isfile(script_path):
        raise Flag(f"no job script {os.path.basename(script_path)} to raise the memory in")
    with open(script_path) as fh:
        text = fh.read()
    rx = re.compile(r"^(#PBS\s+-l\s+mem=|#SBATCH\s+--mem(?:-per-cpu)?=)(\d+(?:\.\d+)?)(\S*)", re.M)
    if not rx.search(text):
        raise Flag("no memory request found in the job script")
    return rx.sub(lambda m: f"{m.group(1)}{math.ceil(float(m.group(2)) * factor)}{m.group(3)}", text)


def same_cell(a, b, tol=1e-4):
    return all(abs(a["lattice"][i][j] - b["lattice"][i][j]) < tol for i in range(3) for j in range(3))


def plan(d, r, entry, opts):
    """(list of change descriptions, callable that applies them). Raises Flag."""
    action = r.get("action") or {}
    if not action:
        raise Flag(f"no action for status {r['status']}")
    incar = herlib.parse_incar(os.path.join(d, "INCAR"))
    act = choose(action, incar, entry, r["error"] or r["status"])
    if "manual" in act:
        raise Flag(f"needs judgement: {act['manual']}")
    if "protected" in act:
        raise Flag(f"fix would change {act['protected']} (protected)")
    changes = dict(act.get("incar", {}))
    if act.get("magmom_fit"):
        changes["MAGMOM"] = magmom_fit(d, incar)
    if act.get("nbands_factor"):
        changes["NBANDS"] = math.ceil(current_nbands(d, incar) * act["nbands_factor"])
    bad = (set(changes) & PROTECTED_INCAR) | (set(act.get("remove", [])) & PROTECTED_FILES)
    if bad:
        raise Flag(f"fix would change {', '.join(sorted(bad))} (protected)")
    desc = [f"{k}={v}" for k, v in changes.items()]
    new_script = None
    if "sub" in act:
        script = os.path.join(d, opts["script"])
        new_script = scale_memory(script, act["sub"]["mem_factor"])
        desc.append(f"memory x{act['sub']['mem_factor']} in {opts['script']}")
    if act.get("contcar_to_poscar"):
        if herlib.is_frequency(incar) or incar.get("NSW", "1").strip() == "0":
            raise Flag("cp CONTCAR POSCAR asked for a frequency/static job")
        cpath = os.path.join(d, "CONTCAR")
        if not os.path.isfile(cpath) or os.path.getsize(cpath) == 0:
            raise Flag("CONTCAR missing or empty")
        try:
            cont, pos = herlib.read_structure(cpath), herlib.read_structure(os.path.join(d, "POSCAR"))
        except (OSError, ValueError, IndexError) as e:
            raise Flag(f"CONTCAR/POSCAR unreadable: {e}")
        if len(cont["elements"]) != len(pos["elements"]) or not same_cell(cont, pos):
            raise Flag("CONTCAR has a different cell or atom count than POSCAR (cell is protected)")
        desc.append("cp CONTCAR POSCAR")
    removes = [f for f in act.get("remove", []) if os.path.exists(os.path.join(d, f))]
    desc += [f"rm {f}" for f in removes]
    if not desc:
        desc.append("resubmit unchanged")

    def apply():
        if changes:
            herlib.set_incar(os.path.join(d, "INCAR"), changes)
        if new_script is not None:
            with open(os.path.join(d, opts["script"]), "w") as fh:
                fh.write(new_script)
        if act.get("contcar_to_poscar"):
            shutil.copyfile(os.path.join(d, "CONTCAR"), os.path.join(d, "POSCAR"))
        for f in removes:
            os.remove(os.path.join(d, f))
    return desc, apply


def fix(d, opts, state):
    r = check(d, opts["log"], opts["use_scheduler"])
    key = herlib.job_key(d, opts["state"])
    entry = herlib.job_entry(state, key)
    out = {"dir": d, "status": r["status"], "error": r["error"], "result": "", "change": "", "jobid": ""}
    if r["status"] in ("converged", "running"):
        out["result"] = "nothing to fix"
        return out
    if entry.get("status") == "flagged":
        out["result"] = "skipped: already flagged (" + entry.get("flag_reason", "") + ")"
        return out
    if entry.get("status") == "ready":
        out["result"] = "skipped: fixed, waiting for submission"
        return out
    error = r["error"] or r["status"]
    try:
        if entry.get("retries", 0) >= opts["max_retries"]:
            raise Flag(f"retry limit reached ({entry.get('retries', 0)} retries)")
        desc, apply = plan(d, r, entry, opts)
        n = entry.get("retries", 0) + 1
        backups = [(opts["log"], f"{opts['log']}-{n}"), ("INCAR", f"INCAR-{n}")]
        for src, dst in backups:
            if os.path.exists(os.path.join(d, dst)):
                raise Flag(f"backup {dst} already exists")
    except Flag as e:
        out.update(result="flagged", change=str(e))
        if not opts["dry_run"]:
            entry.update(status="flagged", flag_reason=str(e), last_error=error)
            herlib.log_event(opts["event_log"], "flagged", key, error, str(e), r["jobid"])
        return out
    out["change"] = "; ".join(desc)
    if opts["dry_run"]:
        out["result"] = "would fix"
        return out
    if os.path.exists(os.path.join(d, opts["log"])):
        os.rename(os.path.join(d, opts["log"]), os.path.join(d, f"{opts['log']}-{n}"))
    shutil.copyfile(os.path.join(d, "INCAR"), os.path.join(d, f"INCAR-{n}"))
    apply()
    entry["retries"] = n
    entry["last_error"] = error
    entry["status"] = "ready"
    jobid, submit_error = "", ""
    if opts["submit"]:
        try:
            jobid = herlib.submit(d, opts["scheduler"], opts["script"], opts.get("submit_cmd"))
            entry.update(status="queued", jobid=jobid)
        except RuntimeError as e:
            submit_error = str(e)
    entry.setdefault("history", []).append({"time": herlib.now(), "error": error, "change": out["change"],
                                            "jobid": jobid})
    herlib.log_event(opts["event_log"], "fixed", key, error, out["change"], jobid, retry=n)
    if jobid:
        herlib.log_event(opts["event_log"], "submitted", key, error, "", jobid, retry=n)
    out.update(result="fixed" + (" and resubmitted" if jobid else f"; NOT submitted: {submit_error}"
                                 if submit_error else " (resubmit with --submit)"), jobid=jobid)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--submit", action="store_true", help="resubmit after fixing")
    ap.add_argument("--scheduler", default="auto", choices=["auto", "pbs", "slurm"])
    ap.add_argument("--script", default="sub.sh", help="job script in each folder")
    ap.add_argument("--submit-cmd", help=argparse.SUPPRESS)  # testing: replaces qsub/sbatch
    ap.add_argument("--max-retries", type=int, default=3)
    ap.add_argument("--state", default=herlib.STATE_FILE)
    ap.add_argument("--event-log", default=herlib.LOG_FILE)
    ap.add_argument("--log", default="log", help="VASP stdout file written by the job script")
    ap.add_argument("--no-scheduler", action="store_true", help="do not query qstat/squeue")
    ap.add_argument("--dry-run", action="store_true", help="show what would be done; change nothing")
    ap.add_argument("--json", help="also write the results to this file")
    a = ap.parse_args()
    opts = {"log": a.log, "use_scheduler": not a.no_scheduler, "state": a.state, "event_log": a.event_log,
            "max_retries": a.max_retries, "submit": a.submit, "scheduler": a.scheduler, "script": a.script,
            "submit_cmd": a.submit_cmd, "dry_run": a.dry_run}
    state = herlib.load_state(a.state)
    results = [fix(d, opts, state) for d in a.dirs]
    if not a.dry_run:
        herlib.save_state(state, a.state)
    for r in results:
        print(f"{r['result']:<34} {r['dir']}  [{r['error'] or r['status']}] {r['change']}"
              + (f"  JOBID {r['jobid']}" if r["jobid"] else ""))
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(results, fh, indent=1)


if __name__ == "__main__":
    main()
