# Screening program, candidate log and POSCAR building

## The screening step (run by the user)

The user's screening program is private and is run by the user, not the agent.

- At the start of each cycle, after the dataset has been updated with the last cycle's DFT results, remind the
  user: "The dataset is updated. Please run your screening program to pick new structures and send me the log."
- Then stop and wait. Do not look for, open, run or describe the program, and do not guess candidates.
- When the user provides the log file, continue with parsing (below).

## Log format

The program writes blocks like this (values below are an example only):

```
======================================================================
Top 10 最优配置 (ΔG ≈ 0)
======================================================================

配置 14:
  Z=44, 配体=[15, 8, 7, 16]  (PONS)
  ΔG = -0.0244 ± 0.5133
```

- `Z` = atomic number of the transition metal (44 → Ru).
- `配体` = atomic numbers of the four ligands X1–X4 in the program's order (15, 8, 7, 16 → P, O, N, S).
- Label in brackets = ligand symbols in the same order; name the configuration `Ru-PONS@G`.
- `ΔG = value ± σ` = predicted value and its uncertainty.

`python scripts/parse_candidates.py <log>` turns the log into `candidates.csv` (section, config ID, name,
canonical, Z, ligands, ΔG, σ). `canonical` is the symmetry key: the lexicographically smallest ligand string
over the four C2v operations of the site rectangle (e.g. `Ru-PONS@G` and `Ru-SNOP@G` both → `Ru-NSPO@G`;
Ir-CCNS, Cr-PSSS and Y-OPPS are already canonical). It flags blocks where the bracket label does not match the ligand numbers, and blocks with a
missing ΔG line; report those to the user instead of guessing.

## Choosing what to run

Show the parsed table to the user and ask which candidates to run this cycle. Skip any configuration that is
already in the dataset or already has a job folder; compare by the `canonical` column, not by `name`,
so symmetry-equivalent orderings are caught.

## Building POSCARs

`python scripts/build_poscar.py candidates.csv --template <POSCAR_template> --sites sites.json [--only NAME ...]`

- Template: `assets/POSCAR_template` (from the relaxed Ir-CCNS@G structure with the two H removed;
  31 atoms, 9.84 × 8.52 Å cell, c = 15 Å; no fixed atoms).
- Site map: `assets/sites.json` = `{"metal": 31, "ligands": [27, 28, 29, 30]}` (1-based).
- Ligand convention (confirmed against two calculated structures, Ir-CCNS@G and Y-OPPS@G): the four sites form
  a rectangle read like a 2 × 2 grid, X1 X2 / X3 X4. X1–X4 and X2–X3 are the diagonals
  (atoms 27–30 and 28–29). This matches the C2v transform in the paper (s3 = X1 − X2 − X3 + X4 is the
  C2-invariant combination). It is NOT a cyclic order around the metal.
- The template is a relaxed Ir structure, so for other metals it is only a starting geometry; every new POSCAR
  must be fully relaxed.
- The script swaps elements, regroups atoms by element for the POSCAR, keeps selective-dynamics flags,
  and never overwrites an existing POSCAR. Generate KPOINTS and POTCAR from the new POSCAR with
  `echo -e "102\n1\n0.03" | vaspkit` so the element order matches, then check it against the POTCAR `TITEL` lines.

## Loop control

- One cycle = user runs screening and sends the log → parse log → user picks candidates → DFT → results into the dataset.
- Stopping criterion: TODO (user to decide).
- Per cycle, log: cycle number, dataset size, candidates picked, DFT results, flags.
