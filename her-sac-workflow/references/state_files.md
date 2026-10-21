# State and log files

Both files live in the project root (the folder that holds the `<name>/` candidate folders). Scripts write
them; the agent reads them but does not edit them by hand, except to reset a `flagged` job after the user has
dealt with it (set its `status` back to `failed` and say so in the reply).

## `loop_state.json`

One JSON object, rewritten atomically by the scripts.

| Key | Type | Written by | Meaning |
| --- | --- | --- | --- |
| `cycle` | int | agent | active-learning cycle number |
| `status` | string | agent | loop step: `waiting_for_log`, `picking`, `running_dft`, `processing`, `updating` |
| `dataset_size` | int | `update_dataset.py` | rows in the DFT dataset after the last update |
| `candidates` | list of names | agent | candidates picked this cycle |
| `jobs` | object | scripts | one entry per job folder, keyed by its path relative to the project root |
| `jobs.<path>.jobid` | string | `submit_job.py`, `fix_job.py --submit` | scheduler ID of the latest submission |
| `jobs.<path>.status` | string | `submit_job.py`, `check_jobs.py --state`, `fix_job.py` | `queued`, `running`, `done`, `failed`, `ready` (fixed, waiting for submission), `flagged` |
| `jobs.<path>.retries` | int | `fix_job.py` | automatic fixes applied so far (at most 3, then `flagged`) |
| `jobs.<path>.last_error` | string | `fix_job.py` | error type of the last fix (`check_jobs.ERRORS` name, or `unconverged`/`unknown`) |
| `jobs.<path>.flag_reason` | string | `fix_job.py` | why the job was flagged |
| `jobs.<path>.history` | list | `fix_job.py` | one `{time, error, change, jobid}` per fix |
| `added` | list of names | `update_dataset.py` | candidates whose rows were added to the dataset |

Example:

```json
{
 "cycle": 3,
 "status": "running_dft",
 "dataset_size": 412,
 "candidates": ["Ru-PONS@G", "Cr-PSSS@G"],
 "jobs": {
  "Ru-PONS@G/slab": {"jobid": "131234567.gadi-pbs", "status": "done", "retries": 0, "history": []},
  "Ru-PONS@G/H2_C": {"jobid": "131234890.gadi-pbs", "status": "queued", "retries": 1,
                     "last_error": "zbrent",
                     "history": [{"time": "2026-10-21T09:12:03+11:00", "error": "zbrent",
                                  "change": "IBRION=1; cp CONTCAR POSCAR", "jobid": "131234890.gadi-pbs"}]},
  "Cr-PSSS@G/H1_metal": {"jobid": "131233001.gadi-pbs", "status": "flagged", "retries": 0,
                         "last_error": "crash", "flag_reason": "needs judgement: crash: read the log lines above it",
                         "history": []}
 },
 "added": ["Ir-CCNS@G"]
}
```

## `workflow_log.jsonl`

Append-only, one JSON object per line and per event. It replaces the free-text `workflow_log.md`.

| Field | Meaning |
| --- | --- |
| `time` | ISO 8601 with time zone |
| `job_dir` | job path relative to the project root (for `added`: the candidate folder) |
| `event` | `submitted`, `checked`, `fixed`, `converged`, `flagged` or `added` |
| `error_type` | `check_jobs.ERRORS` name, `unconverged` or `unknown`; empty when not applicable |
| `change` | what was changed (`fixed`), the check result (`checked`), the reason (`flagged`), notes (`added`) |
| `jobid` | scheduler ID, when there is one |
| extra keys | `retry` (`fixed`, `submitted`), `retries` (`converged`), `name` and `canonical` (`added`) |

| Event | Written by |
| --- | --- |
| `submitted` | `submit_job.py`, `fix_job.py --submit` |
| `checked`, `converged` | `check_jobs.py --state loop_state.json --event-log workflow_log.jsonl` |
| `fixed`, `flagged` | `fix_job.py` |
| `added` | `update_dataset.py` |

Example lines:

```json
{"time": "2026-10-21T04:00:11+11:00", "job_dir": "Ru-PONS@G/H2_C", "event": "checked", "error_type": "zbrent", "change": "failed", "jobid": "131234567"}
{"time": "2026-10-21T04:00:12+11:00", "job_dir": "Ru-PONS@G/H2_C", "event": "fixed", "error_type": "zbrent", "change": "IBRION=1; cp CONTCAR POSCAR", "jobid": "131234890.gadi-pbs", "retry": 1}
{"time": "2026-10-21T04:00:12+11:00", "job_dir": "Ru-PONS@G/H2_C", "event": "submitted", "error_type": "zbrent", "change": "", "jobid": "131234890.gadi-pbs", "retry": 1}
{"time": "2026-10-21T09:00:05+11:00", "job_dir": "Ru-PONS@G/H2_C", "event": "converged", "error_type": "", "change": "", "jobid": "131234890", "retries": 1}
```

`python scripts/workflow_stats.py` summarises both files: total jobs, converged on the first attempt,
auto-fixed (per error type), total retries, flagged, and the share completed without manual intervention.
