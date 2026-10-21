#!/usr/bin/env python3
"""Run statistics from workflow_log.jsonl and loop_state.json (the numbers for the reviewer response).

    python workflow_stats.py [--state loop_state.json] [--event-log workflow_log.jsonl] [--json stats.json]

Jobs = every job in the state plus every job directory that appears in the log ('added' events, which name a
candidate folder, are counted separately as rows added). Per job:
    converged first attempt  converged, with no 'fixed' and no 'flagged' event
    auto-fixed               at least one 'fixed' event; per error type both the number of fixes and of jobs
    retries                  retry count in the state (or the number of 'fixed' events if larger)
    flagged                  flagged now (state) or at any time (log)
    completed without manual intervention
                             converged (state 'done' or a 'converged' event) and never flagged;
                             reported as a % of all jobs and of finished jobs (converged + currently flagged)
"""
import argparse
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import herlib  # noqa: E402


def stats(state, events):
    jobs = collections.defaultdict(lambda: {"events": [], "retries": 0, "status": ""})
    added = [e for e in events if e["event"] == "added"]
    for e in events:
        if e["event"] != "added":
            jobs[e["job_dir"]]["events"].append(e)
    for key, entry in state.get("jobs", {}).items():
        jobs[key]["retries"] = entry.get("retries", 0)
        jobs[key]["status"] = entry.get("status", "")
    fixes_per_type, jobs_per_type = collections.Counter(), collections.Counter()
    out = {"total_jobs": len(jobs), "converged": 0, "converged_first_attempt": 0, "auto_fixed_jobs": 0,
           "total_retries": 0, "flagged_now": 0, "flagged_ever": 0, "running_or_queued": 0,
           "completed_without_manual": 0}
    for key, j in jobs.items():
        kinds = [e["event"] for e in j["events"]]
        fixed = [e for e in j["events"] if e["event"] == "fixed"]
        converged = j["status"] == "done" or "converged" in kinds
        flagged_ever = j["status"] == "flagged" or "flagged" in kinds
        out["converged"] += converged
        out["converged_first_attempt"] += converged and not fixed and not flagged_ever
        out["auto_fixed_jobs"] += bool(fixed)
        out["total_retries"] += max(j["retries"], len(fixed))
        out["flagged_now"] += j["status"] == "flagged"
        out["flagged_ever"] += flagged_ever
        out["running_or_queued"] += j["status"] in ("queued", "running", "ready")
        out["completed_without_manual"] += converged and not flagged_ever
        for e in fixed:
            fixes_per_type[e["error_type"] or "?"] += 1
        for t in {e["error_type"] or "?" for e in fixed}:
            jobs_per_type[t] += 1
    finished = out["converged"] + out["flagged_now"]
    pct = lambda n, d: round(100.0 * n / d, 1) if d else None
    out["pct_without_manual_of_all"] = pct(out["completed_without_manual"], out["total_jobs"])
    out["pct_without_manual_of_finished"] = pct(out["completed_without_manual"], finished)
    out["fixes_per_error_type"] = dict(sorted(fixes_per_type.items()))
    out["jobs_fixed_per_error_type"] = dict(sorted(jobs_per_type.items()))
    out["rows_added"] = len(added)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--state", default=herlib.STATE_FILE)
    ap.add_argument("--event-log", default=herlib.LOG_FILE)
    ap.add_argument("--json", help="also write the statistics to this file")
    a = ap.parse_args()
    if not os.path.isfile(a.state) and not os.path.isfile(a.event_log):
        raise SystemExit(f"neither {a.state} nor {a.event_log} exists")
    s = stats(herlib.load_state(a.state), herlib.read_log(a.event_log))
    fmt = lambda v: "n/a" if v is None else f"{v:.1f} %"
    print(f"Total jobs                         {s['total_jobs']}")
    print(f"Converged                          {s['converged']}")
    print(f"  on the first attempt             {s['converged_first_attempt']}")
    print(f"Auto-fixed jobs                    {s['auto_fixed_jobs']}")
    for t, n in s["fixes_per_error_type"].items():
        print(f"  {t:<32} {n} fix(es) on {s['jobs_fixed_per_error_type'][t]} job(s)")
    print(f"Total retries                      {s['total_retries']}")
    print(f"Flagged (now / ever)               {s['flagged_now']} / {s['flagged_ever']}")
    print(f"Still queued/running               {s['running_or_queued']}")
    print(f"Completed without manual action    {s['completed_without_manual']}  "
          f"({fmt(s['pct_without_manual_of_all'])} of all jobs, "
          f"{fmt(s['pct_without_manual_of_finished'])} of finished jobs)")
    print(f"Dataset rows added                 {s['rows_added']}")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(s, fh, indent=1)


if __name__ == "__main__":
    main()
