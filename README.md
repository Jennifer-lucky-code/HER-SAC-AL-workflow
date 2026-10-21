# HER-SAC-workflow

An active-learning workflow for hydrogen evolution reaction (HER) screening of single-atom catalysts (SACs), developed by Jennifer-lucky-code.

## Workflow

- Candidate-log parsing
- POSCAR generation and VASP input preparation
- HPC submission and job tracking, with rule-based automatic fixes
- Structure checks, energy collection and HER free-energy analysis
- Dataset updates with provenance, and run statistics

See [her-sac-workflow](her-sac-workflow/) for the workflow files.

Before use, fill in your actual HPC project IDs and storage paths in `her-sac-workflow/references/vasp_settings.md`.

## What is code and what is the agent

The agent (an LLM following `her-sac-workflow/SKILL.md`) only orchestrates. Every number and every decision on a
job comes from VASP output or from a deterministic script in `her-sac-workflow/scripts/`:

| Step | Script | Output |
| --- | --- | --- |
| Parse the screening log | `parse_candidates.py` | `candidates.csv` (with the C2v canonical key) |
| Build structures | `build_poscar.py`, `make_freq_poscar.py` | POSCARs |
| Submit | `submit_job.py` | JOBID, `loop_state.json`, `submitted` event |
| Check jobs every 5 h | `check_jobs.py --state --event-log` | status per job, structured fix |
| Fix failed jobs | `fix_job.py --submit` | one rule-based change per retry, flag after 3 retries |
| Check structures | `check_structure.py` | pass / flag with reasons |
| Collect energies | `collect_energies.py` | `energies.csv` |
| Free energies, mechanism | `her_analysis.py` | `her_results.csv`, diagrams |
| Update the dataset | `update_dataset.py` | appended rows with job paths, date, retries, INCAR differences |
| Run statistics | `workflow_stats.py` | first-attempt convergence, fixes per error type, retries, flags |
| Optional spin-state test | `make_spin_tests.py` | fixed-moment runs and their ΔE / M–X changes |

All actions are logged in `workflow_log.jsonl` (one JSON object per event); the job state is in `loop_state.json`.
Both are described in `her-sac-workflow/references/state_files.md`.

Requirements: Python 3 (standard library); `matplotlib` only for `her_analysis.py --plot`
(`pip install -r her-sac-workflow/requirements.txt`). Run the regression tests with
`cd her-sac-workflow && python -m unittest discover -s tests`.

## Open items

These values are still `TODO` in `her-sac-workflow/references/` and must come from the project, not be guessed:

- Stopping criterion of the active-learning loop (`references/screening.md`).
- E_TM, the bulk energy per metal atom used in E_f (`references/vasp_settings.md`, reference energies).
- AIMD settings: ensemble, temperature, time step, total time, thermostat (`references/vasp_settings.md`).
- Setonix job-script template and modules (`references/vasp_settings.md`, HPC).
- Supercell size (`references/vasp_settings.md`, shared settings).
- QUT Aquarius scheduler and details (`references/vasp_settings.md`, HPC).
