# DFT and HER analysis steps

Read this for loop steps 4–8 (structure generation, VASP jobs, checks, free energies, mechanism, stability).

## 1. Model and naming

- Model: a transition-metal (TM) atom coordinated by four ligand atoms X ∈ {C, N, O, P, S} in a graphene supercell.
- Name each configuration `TM-X1X2X3X4@G` (e.g. `Ir-CCNS@G`), using the canonical (rotation-reduced) ligand order.
- Symmetry reduction: map every ligand arrangement to its canonical rotationally invariant form before running
  DFT, so equivalent structures are calculated once. 28 TMs × 5⁴ arrangements = 17,500 nominal → 4,900 unique.
- Keep one folder per configuration, `<name>/`, with exactly these job sub-folders (the scripts rely on the names):

  | Folder | Job |
  | --- | --- |
  | `slab/` | bare catalyst `*` (start from the `<name>/POSCAR` written by `build_poscar.py`) |
  | `H1_metal/` | `H*` on top of the metal |
  | `H1_lig<k>/` | `H*` on ligand X<k>, k = 1–4 in the X1..X4 order of the name (only inequivalent sites) |
  | `H2_C/` | `2H*-C`: Tafel, second H on a ligand atom (any X, not only C; the `_C` name is kept for the scripts) |
  | `H2_M/` | `2H*-M` ("horns"): both H on the metal |
  | `freq_<state>/` | frequency job of an adsorbed state, e.g. `freq_H1_metal/`, `freq_H1_lig2/`, `freq_H2_C/`, `freq_H2_M/`; keep the vaspkit 501 output there as `vaspkit_501.txt` |

  `python scripts/collect_energies.py <name> [...]` reads this layout and writes `energies.csv`.

## 2–4. VASP calculations

Settings used in the TMX4@G project are in `references/vasp_settings.md`. Read it before writing any INCAR.

- **Bare catalyst (`*`)**: relax, spin-polarised.
- **First H (`H*`)**: try H on top of the metal and on each inequivalent ligand site; keep the lowest-energy
  converged structure and record its site.
- **Second H**: the second H reacts either by the Heyrovsky or by the Tafel mechanism.
  - **Heyrovsky** (H* + H⁺ + e⁻ → * + H₂): no extra calculation; its step follows from ΔG_H* alone
    (ΔG_Heyrovsky = −ΔG_H*, §7).
  - **Tafel** (2H* → * + H₂): needs a `2H*` structure. Start from the best `H*` and always compute BOTH cases:
    - **2H*-C**: the second H adsorbs on a ligand atom X (C, N, O, P or S) next to the metal, i.e. metal–ligand
      cooperative recombination. The label "C" is only the folder/column name, not the element.
    - **2H*-M ("horns")**: the second H also adsorbs on the metal, both H on the same metal atom pointing up like
      two horns (dihydride).
    Each gives its own Volmer–Tafel route (VT-C, VT-M); keep both, do not keep only the lower-energy one.
- After relaxation, check the H geometry (section 6): an H–H distance below ~0.9 Å means H₂ formed or a
  dihydrogen complex; label it, do not silently treat it as 2H*.

## 5. Frequencies (ZPE and TS)

Fixed atoms are used ONLY here. Structure relaxations (`*`, `H*`, `2H*`) fix no atoms.
Make the frequency POSCAR from the relaxed CONTCAR with `python scripts/make_freq_poscar.py CONTCAR`
(all substrate atoms F F F, all H atoms T T T).


- Run a frequency job for every adsorbed state, with only the H atoms free and the slab fixed.
- G = E_DFT + E_ZPE − TS, T = 298.15 K (ZPE/TS via VASPKIT 501, E_DFT = OSZICAR `E0`; commands in
  `references/vasp_settings.md`). The bare slab uses G_* = E_*.
- G_H₂ for gas-phase H₂: use the project reference value in `references/vasp_settings.md`; never mix references
  from different setups.
- Imaginary modes on H: report them; the structure may not be a minimum (`check_structure.py` lists them).
- `python scripts/collect_energies.py <name> [...]` assembles E, ZPE, TS and G from the folders of §1 into
  `energies.csv`; save each vaspkit 501 output as `freq_<state>/vaspkit_501.txt` (or use `--run-vaspkit`).

## 5b. Optional: spin-state test (shortlisted open-shell candidates only)

Not part of the standard loop; run it only when the user asks, for shortlisted candidates whose converged
state is open-shell. Skip Y and the other formally d⁰ (group-3) metals.

1. `python scripts/make_spin_tests.py <name>/<state>` on a converged folder (e.g. `Cr-PSSS@G/slab`) writes
   `spin_N<k>/` sub-folders, one per parity-allowed fixed total moment (`NUPDOWN` = 0 or 1 up to `--max`,
   step 2). The metal and its four ligands are displaced randomly by 0.02–0.04 Å (seeded; displacements are
   recorded in `spin_tests.json`), `ISYM = 0` is set, every other INCAR key is unchanged.
2. Submit the sub-folders like any other job (`scripts/submit_job.py`) and check them with the 5-hourly run.
3. `python scripts/make_spin_tests.py <name>/<state> --summarise` tabulates ΔE against the original solution
   and the M–X bond-length changes. Report the table; do not replace the dataset values with a lower-energy
   spin state without the user's decision.

## 6. Checks before any result enters the dataset

A job passes only if all hold:

- Electronic SCF converged in every ionic step, and the final ionic step reached EDIFFG (forces < 0.01 eV/Å).
  The job log must end with `reached required accuracy - stopping structural energy minimisation`
  (last or second-to-last line).
- No error messages in the output or scheduler log; the job finished normally.
- H stayed bound: metal–H or ligand–H distance within a sensible bonding range; no H drifted into vacuum.
- Final magnetic moment recorded; note large changes between `*`, `H*` and `2H*`.
- Structure integrity: TM still four-coordinated, no large out-of-plane distortion. Flag symmetry breaking
  (Jahn–Teller-like distortion) in the log rather than discarding the job.

The first two are `scripts/check_jobs.py`; the rest is `scripts/check_structure.py <job dirs> --json
structure_check.json`, which judges only converged jobs and reports `pass` or `flag` with reasons. It finds the
metal (the single TM atom) and its ligands (4 nearest non-H neighbours) by geometry, with minimum-image
distances, and checks: coordination number 4 (covalent-radius bond criterion), the 4 M–X bond lengths and
their spread (spread > 0.10 Å → note "possible Jahn–Teller distortion", never a failure), metal displacement
from the ligand plane (flag > 0.5 Å), every H bound and not in the vacuum, H–H < 0.9 Å → "H₂ / dihydrogen
formed", final `mag=` from OSZICAR, and imaginary (`f/i=`) modes of frequency jobs. Every threshold is a
command-line option; the defaults are listed at the top of the script.

Failed jobs are fixed automatically (no need to ask the user) by scripts, never by hand. The 5-hourly
`scripts/check_jobs.py --state loop_state.json --event-log workflow_log.jsonl` run finds them; every row of its
`ERRORS` table (walltime, memory, POTCAR, MAGMOM, ZBRENT, BRIONS, LREAL, symmetry, eigensolver, BRMIX, NBANDS,
crashes) carries a structured action. `scripts/fix_job.py <dir> --submit` then applies exactly this procedure:

1. Back up before changing anything: `log` → `log-N`, `INCAR` → `INCAR-N` (N = retry number).
2. Apply the one action the table names for the error (its `else` alternative when the first one is already in
   place or the same error came back). When the action says so (relaxations that had made progress),
   `cp CONTCAR POSCAR` first so the restart keeps it.
3. Never change anything that shifts the energy reference to make a job run: ENCUT, PREC, IVDW, ISPIN,
   POTCAR set, KPOINTS, cell. If the only fix would need that, the job is flagged instead.
4. Resubmit (`--submit`), write the new JOBID, add the retry to `loop_state.json`.
5. After 3 retries of the same job, stop and flag it for the user with the errors seen. Actions that need
   judgement (`crash`, `unknown`) are never auto-fixed: the job is flagged and reported.

Every action (what failed, what was changed, new job ID) is one JSON line in `workflow_log.jsonl` (events
`submitted / checked / fixed / converged / flagged / added`, schema in `references/state_files.md`);
`scripts/workflow_stats.py` turns it into run statistics. `update_dataset.py` writes any INCAR key that differs
from the template into the dataset row, so it is clear which jobs deviate from the template INCAR.

## 7. Free energies, mechanism and overpotential

Run `python scripts/her_analysis.py energies.csv --plot` on the `energies.csv` from `collect_energies.py`
(input columns are documented at the top of the script; rows with a blank E_slab or G_H are skipped and reported).
It computes, with the reference state * + 2H⁺ + 2e⁻ = 0 at U = 0 V:

- ΔG_H* = G_H* − G_* − ½G_H₂ (first Volmer)
- ΔG_2H* = G_2H* − G_H* − ½G_H₂ (second Volmer, a step value, not the cumulative 2H* level)
- ΔG_Tafel = −(ΔG_H* + ΔG_2H*) (chemical step, no electron transfer)
- ΔG_Heyrovsky = −ΔG_H*

Mechanism: only uphill (positive) steps cost energy.

- η_VH = largest positive step of {ΔG_H*, −ΔG_H*} = |ΔG_H*|
- η_VT = largest positive step of {ΔG_H*, ΔG_2H*, ΔG_Tafel}, computed separately for VT-C and VT-M
- Preferred pathway = the smallest of η_VH, η_VT-C, η_VT-M; report its value as the limiting potential.
- If one 2H* configuration is missing, say so next to the result.
- If the largest VT step is the Tafel step, say so: it is chemical and cannot be lowered by applied potential.

Keep the project's screening descriptor separate and clearly named: η_HER = max(|ΔG_H*|, |ΔG_2H*|) (Eq. 4 of
the TMX4@G paper, using the more stable 2H* configuration). Never mix screening-descriptor values and mechanism values in one column or one figure axis.

Worked check (Cr-PSSS@G, from the project data): ΔG_H* = −0.352; 2H* levels −0.397 (C) and −0.401 (M)
→ η_VH = 0.352, η_VT-C = 0.397, η_VT-M = 0.401 → Volmer–Heyrovsky, 0.35 V. Use this to test any new script version.

## 8. Stability

- E_f = E_TMX4@G − E_X4@G − E_TM (E_TM per atom of the bulk metal). E_f < 0 → thermodynamically stable.
- For shortlisted candidates only: AIMD (see `references/vasp_settings.md`) and check the TM stays anchored.

## 9. Dataset row (DFT part)

Append one row per configuration to the dataset with: name, TM, ligands, E_f, ΔG_H*, H* site, ΔG_2H* (C and M),
η_VH, η_VT-C, η_VT-M, preferred pathway, η_HER (screening), magnetic moment, d-band centre (VASPKIT), job paths, date.
Never overwrite existing rows; if a value changes, add a new row with a note and keep the old one.

`python scripts/update_dataset.py her_results.csv --structure structure_check.json [--extra features.csv]` does
this. It adds the canonical key, the magnetic moments of `*`, `H*`, `2H*-C` and `2H*-M`, the retry count from
`loop_state.json` and every INCAR key that differs from the templates in `references/vasp_settings.md`
(SYSTEM, MAGMOM and NPAR are allowed to differ). It refuses a candidate whose jobs are not converged, and one with
a structure flag unless `--allow-flagged` is given (the flags then go into `notes`). E_f and the d-band centre
come from the optional `--extra` CSV (columns `name, Ef, d_band_centre`) and stay blank otherwise.

