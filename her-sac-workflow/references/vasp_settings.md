# VASP settings (TMX4@G project)

Values marked TODO were not in the manuscript. Fill them in from the project's own input files before first use;
until then, ask the user rather than guessing.

## Shared settings (from the manuscript methods)

| Setting | Value |
| --- | --- |
| Code | VASP 5.4.4, PAW, GGA-PBE |
| Spin | ISPIN = 2 (spin-polarised) |
| ENCUT | 500 eV |
| k-points | Monkhorst-Pack, from VASPKIT k-spacing 0.03 (gives 4 × 4 × 1 for the template cell; see KPOINTS below) |
| EDIFF (relaxation) | 1E-5 |
| EDIFF (static) | 1E-6 |
| EDIFFG | -0.01 (forces < 0.01 eV/Å) |
| Pre/post-processing | VASPKIT (d-band centre, PDOS) |
| Vacuum | 15 Å (user-stated; not yet written in the manuscript) |
| Supercell size | TODO |
| vdW correction | DFT-D3 (IVDW = 11, zero damping) |
| Dipole correction (LDIPOL/IDIPOL) | not used (not in the project INCAR) |
| Initial MAGMOM | not set in template (VASP default 1 μB per atom); set explicitly if spin states matter |
| POTCAR | generated with VASPKIT together with KPOINTS (its recommended PAW_PBE set; see KPOINTS below) |
| Solvation (VASPsol) | not used (not in the project INCAR) |

## Relaxation INCAR (project template)

Use this as is for `*`, `H*` and `2H*` relaxations. Only change SYSTEM, MAGMOM (if needed) and NPAR (machine-dependent).

```
# Collinear magnetic calculation
ISPIN      = 2
LASPH      = .TRUE.
GGA_COMPAT = .FALSE.
VOSKOWN    = 1
LMAXMIX    = 4        # d elements; use 6 for f elements
AMIX       = 0.2
BMIX       = 0.0001
AMIX_MAG   = 0.4
BMIX_MAG   = 0.0001

# DFT-D3 correction
IVDW   = 11           # DFT-D3, zero damping

# Global parameters
ISTART = 0            # start from scratch
LREAL  = Auto
ENCUT  = 500
PREC   = Accurate
LWAVE  = .TRUE.
LCHARG = .FALSE.      # set .TRUE. when Bader charges are needed
ADDGRID = .TRUE.
NPAR   = 8            # depends on node count

# Electronic relaxation
ISMEAR = 0
SIGMA  = 0.05
NELM   = 90
NELMIN = 6
EDIFF  = 1E-05        # 1E-06 for static single-point runs

# Ionic relaxation
NSW    = 500
IBRION = 2
ISIF   = 2
EDIFFG = -1E-02

LORBIT = 12
```

Static single-point (for ε_d, PDOS, Bader): same file with NSW = 0, IBRION = -1, EDIFF = 1E-06,
LCHARG = .TRUE., plus LAECHG = .TRUE. for Bader.

## KPOINTS (project template)

```
K-Spacing Value to Generate K-Mesh: 0.030
0
Monkhorst-Pack
   4   4   1
0.0  0.0  0.0
```

Production runs use a plain Monkhorst-Pack grid (confirmed by the user). The manuscript methods say
"Γ-centred Monkhorst–Pack"; that wording is wrong for the actual runs, keep new jobs Monkhorst-Pack.

Generate KPOINTS and POTCAR in every job folder (after the final POSCAR is in place) with:

```bash
echo -e "102\n1\n0.03" | vaspkit
```

(102 = KPOINTS + POTCAR from k-spacing, 1 = Monkhorst-Pack, 0.03 = k-spacing. Plain numbers, no prefix.)
Then check POSCAR–POTCAR consistency: the element order on POSCAR line 6 must match the order of the
`TITEL` lines in POTCAR. If it does not, report it to the user; do not submit.

## Frequency INCAR (project template)

Fix all slab atoms with selective dynamics (F F F); leave only the H atoms free (T T T).
Use `scripts/make_freq_poscar.py`. This is the only job type with fixed atoms; relaxations fix nothing.

```
# Collinear magnetic calculation
ISPIN  = 2
LPLANE = .TRUE.
LSCALU = .FALSE.
NSIM   = 4

# DFT-D3 correction
IVDW   = 11

# Global parameters
ISTART = 0
LREAL  = Auto
ENCUT  = 500
PREC   = Normal
LWAVE  = .FALSE.
LCHARG = .FALSE.
ADDGRID = .TRUE.
ALGO   = Fast

# Electronic relaxation
ISMEAR = 0
SIGMA  = 0.05
NELM   = 90
NELMIN = 6
EDIFF  = 1E-06        # tighter for frequencies
EDIFFG = -1E-02

# Frequency calculation
IBRION = 5            # finite differences
NFREE  = 2            # central differences
POTIM  = 0.015        # displacement in Å
NSW    = 1

LORBIT = 12
```

ZPE and TS at 298.15 K, run inside the finished frequency folder:

```bash
printf "501\n298.15\n" | vaspkit | grep -E "Zero-point energy E_ZPE|Entropy contribution T\*S"
```

```text
Zero-point energy E_ZPE   :      13.814 kcal/mol    0.599026 eV
Entropy contribution T*S  :   15147.343 J/(mol)     0.156991 eV
```

ZPE and TS are the second-to-last column (eV) of these two lines. G = E + ZPE − TS for every adsorbed
state, with E from the matching relaxation (below).

## Convergence and energies

- A relaxation counts as converged only if the last or second-to-last line of the job log (`log`, from the
  job script below) contains `reached required accuracy - stopping structural energy minimisation`.
  Do not use its CONTCAR for frequency jobs or results otherwise.
- E (E_DFT) = the `E0` value on the last line of that job's `OSZICAR`.
- Whenever MAGMOM is set explicitly, its total count must equal the POSCAR atom count (lines 6–7),
  otherwise VASP stops.

## Reference energies

| Quantity | Value (eV) | Note |
| --- | --- | --- |
| G_H₂ | -6.8901 | From the project spreadsheet (H2+surface − surface). G = E + ZPE − TS, same as every other step (confirmed by the user). |
| E_TM (bulk, per atom) | TODO | One value per metal |

## AIMD (shortlisted candidates)

TODO: ensemble, temperature, time step, total time, thermostat. Pass criterion: TM remains anchored and
coordination is retained.

## HPC

Note: Before use, replace `{PROJECT_ID}` with your actual HPC allocation/project ID, `{STORAGE_PROJECT_ID}` with your storage project ID, and `{STORAGE_PATHS}` with the storage paths required by your HPC system.

| System | Scheduler | Project | Notes |
| --- | --- | --- | --- |
| NCI Gadi | PBS Pro | {PROJECT_ID} (storage {STORAGE_PROJECT_ID} + {PROJECT_ID}) | template below |
| Pawsey Setonix | SLURM | {PROJECT_ID} / {STORAGE_PROJECT_ID} | TODO: job script template, modules |
| QUT Aquarius | TODO | | |

### Gadi job script (project template)

Fill every `{PLACEHOLDER}` before submitting; ask the user for any value not given for the job.
`{NCPUS}` must be a multiple of the cores per node of `{QUEUE}`, and NPAR a divisor of `{NCPUS}`.

```bash
#!/bin/bash
#PBS -P {PROJECT_ID}
#PBS -N {JOBNAME}
#PBS -q {QUEUE}
#PBS -l ncpus={NCPUS}
#PBS -l walltime={WALLTIME}
#PBS -l mem={MEM}
#PBS -l software=vasp
#PBS -l wd
#PBS -j oe
#PBS -l storage={STORAGE_PATHS}

module purge
module load vasp/5.4.4

mpirun -np {NCPUS} vasp_std > log
```
