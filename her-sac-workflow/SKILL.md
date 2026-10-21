---
name: her-sac-workflow
description: Closed-loop, active-learning HER screening workflow for single-atom catalysts (SACs) such as TMX4@G / M-N4/C. Covers the full loop - asking the user to run their private screening program, parsing the candidate log they provide (metal Z, ligands, predicted ΔG ± σ), POSCAR generation by element substitution, VASP input assembly, HPC submission and job tracking, failure handling, result processing (H*/2H* free energies, Volmer–Heyrovsky vs Volmer–Tafel mechanism, overpotential, stability, feature extraction) and dataset update before retraining. Use this skill whenever the user asks to run, continue, check or analyse an active-learning cycle, select or submit catalyst candidates, set up or parse HER calculations, compute ΔG_H*, ΔG_2H*, overpotentials or free-energy diagrams, or update an HER/ML dataset, even if they do not say "skill" or "workflow".
---

# Active-learning HER workflow for single-atom catalysts

Reusable version of the TMX4@G workflow.

Note: Before use, fill in your actual HPC project IDs and storage paths as described in `references/vasp_settings.md`.

**Ground rule:** the agent orchestrates; it never invents numbers. Structures come from scripts, energies from
VASP, descriptors and free energies from scripted post-processing, predictions from the user's screening program.
If a value is missing, report it instead of estimating it.

## The loop

```
 DFT-labelled dataset ──▶ (1) user runs screening ──▶ (2) parse log ──▶ (3) user picks candidates
        ▲                                                                           │
        │                                                                           ▼
 (8) update dataset ◀── (7) result processing ◀── (6) track jobs ◀── (5) HPC ◀── (4) POSCAR + inputs
```

Keep a state file `loop_state.json` in the project root: cycle number, loop step, dataset size, candidates
picked, each job's path, scheduler ID and status (`queued / running / done / failed / ready / flagged`), retry
count, and which candidates were added. The scripts write the job part and log every action to
`workflow_log.jsonl`; the agent sets cycle, step and candidates (schemas in `references/state_files.md`).
Always read it first when resuming, so a cycle is never repeated or skipped.

## Steps

**1. Ask the user to run the screening program.** The agent does not run it. When a new cycle starts
(dataset updated with the last results), remind the user to run their screening program to pick new
structures, then stop and wait. Set `loop_state.json` to `waiting_for_log`. Continue only when the user
provides the log file.

**2. Parse the log the user provided.** `python scripts/parse_candidates.py <log>` → `candidates.csv`
(name such as `Ru-PONS@G`, metal Z, ligand atomic numbers, ΔG ± σ). Report flagged blocks.

**3. Pick candidates.** Show the table and ask the user which ones to run this cycle. Skip anything already
calculated.

**4. Build inputs.** `python scripts/build_poscar.py candidates.csv --template assets/POSCAR_template --sites assets/sites.json --only <names>`
(element substitution on the project template, see `references/screening.md`), then copy the INCAR from `references/vasp_settings.md` and generate KPOINTS + POTCAR with
`echo -e "102\n1\n0.03" | vaspkit` (`submit_job.py` checks the POSCAR/POTCAR element order). Job set per
candidate, in the folders named in `references/dft_her.md` §1: `slab/`, `H1_metal/` and `H1_lig<k>/`, `H2_C/`
and `H2_M/`, and `freq_<state>/`. Relaxations fix no atoms; only frequency (ZPE/TS) jobs fix the substrate
(`scripts/make_freq_poscar.py`). Details in `references/dft_her.md`.

**5. Submit.** Submit to the HPC scheduler (Gadi: PBS; Setonix: SLURM). Ask the user before the first
submission of a session. Submit with `python scripts/submit_job.py <job dirs>`: it checks the input files and the
POSCAR/POTCAR element order, runs `qsub sub.sh` (SLURM: `sbatch --parsable sub.sh`), writes each folder's
`JOBID`, registers the job in `loop_state.json` and logs a `submitted` event. Then set up the 5-hourly check (step 6).

**6. Check every 5 h and fix failures automatically.** While a cycle has jobs out, check them every 5 hours:
schedule a recurring check with whatever scheduler this agent has (e.g. `/loop 5h`, a cron task or a
Routine) and cancel it once every job of the cycle is `done` or `flagged`. Each check runs two scripts and
nothing else decides a fix:
1. `python scripts/check_jobs.py <job dirs> --state loop_state.json --event-log workflow_log.jsonl --json job_status.json`
   (logs `checked`/`converged`, sets `done` / `running` / `failed`).
2. `python scripts/fix_job.py <failed or unconverged dirs> --submit` applies the one structured action
   `check_jobs` names (backup, change, `cp CONTCAR POSCAR` when the action says so, resubmit, new JOBID,
   retry + 1, `fixed` event). It flags instead of fixing when the action needs judgement (`crash`, `unknown`),
   would touch ENCUT, PREC, IVDW, ISPIN, POTCAR, KPOINTS or the cell, or the job already had 3 retries.
Then: start the jobs that depend on newly `done` ones (e.g. frequency jobs, via `submit_job.py`); tell the
user about every `flagged` job with its reason; give a one-line summary per check (counts of
done/running/fixed/flagged). Never edit an INCAR by hand to fix a job. Every action is in
`workflow_log.jsonl` (one JSON object per event; schema in `references/state_files.md`).

**7. Process results.** For each finished candidate (every number comes from a script, never typed in):
- structure checks of `references/dft_her.md` §6: `python scripts/check_structure.py <name>/*/ --json structure_check.json`
  (pass / flag with reasons; a Jahn–Teller spread is only a note). Report every flag.
- energies: `python scripts/collect_energies.py <names> [--run-vaspkit]` → `energies.csv` (E0 of converged
  relaxations, ZPE/TS from `freq_*/vaspkit_501.txt`, lowest-G H* site, G_H₂ from `vasp_settings.md`; missing
  values stay blank with a note)
- free energies, VH/VT mechanism and limiting potential: `python scripts/her_analysis.py energies.csv --plot`
- formation energy E_f (and AIMD only for shortlisted candidates)
- feature extraction for the ML dataset: ε_d (VASPKIT), Bader charge transfer, d_M−X, d_M−H, plus the
  screening target η_HER (Eq. 4) kept separate from the mechanism values.
- optional, only when the user asks: spin-state test for shortlisted open-shell candidates
  (`scripts/make_spin_tests.py`, `references/dft_her.md` §5b).

**8. Update the dataset and loop.** `python scripts/update_dataset.py her_results.csv --structure structure_check.json [--extra features.csv]`
appends one row per validated candidate (deduplicated by the canonical key, never overwriting a row; rows with
a structure flag are refused unless the user agrees to `--allow-flagged`). Report refused rows. Run
`python scripts/workflow_stats.py` for the cycle's job statistics, then return to step 1. Stop when the
stopping criterion in `references/screening.md` is met, and report the cycle summary.

## Reference files

- `references/screening.md`: the user-run screening step, log format, POSCAR building, loop control.
- `references/dft_her.md`: structures, VASP jobs, checks, free energies, VH/VT mechanism, stability.
- `references/vasp_settings.md`: INCAR/KPOINTS settings, reference energies, HPC details.
- `references/state_files.md`: `loop_state.json` and `workflow_log.jsonl` schemas with examples.
- `scripts/parse_candidates.py`: screening log → `candidates.csv`.
- `scripts/build_poscar.py`: candidates → POSCARs from the template and site map.
- `scripts/submit_job.py`: input check, submission, JOBID, state and `submitted` event.
- `scripts/check_jobs.py`: job folders → converged / running / failed (error + fix + structured action) / unconverged.
- `scripts/fix_job.py`: applies the structured action of a failed job (backup, one change, resubmit, retry count, flag).
- `scripts/make_freq_poscar.py`: relaxed CONTCAR → frequency POSCAR (substrate fixed, H free).
- `scripts/check_structure.py`: converged job → pass / flag (coordination, out-of-plane, H, magnetic moment, imaginary modes).
- `scripts/collect_energies.py`: candidate folders → `energies.csv`.
- `scripts/her_analysis.py`: ΔG_H*, ΔG_2H* (C and metal "horns" configurations), η_VH, η_VT, pathway, diagrams.
- `scripts/update_dataset.py`: `her_results.csv` + checks → dataset rows with provenance.
- `scripts/workflow_stats.py`: run statistics from the log and state.
- `scripts/make_spin_tests.py`: optional fixed-moment spin-state test.
- `scripts/herlib.py`: shared helpers (geometry, file parsing, settings, state/log); not run directly.
- `assets/POSCAR_template`, `assets/sites.json`: TMX4@G template structure and metal/ligand atom indices
  (template indices only; the other scripts find atoms by geometry).

## Reporting to the user

At the end of each cycle give a short table: candidate, predicted ΔG ± σ, DFT
limiting potential, preferred pathway, E_f, flags.
List failed, flagged or missing items separately. Round to 2 decimals only for display.
