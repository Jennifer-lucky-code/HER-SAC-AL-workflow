#!/usr/bin/env python3
"""Optional spin-state test for shortlisted open-shell candidates (references/dft_her.md §5b).

Set up, from a converged job folder:
    python make_spin_tests.py <job_dir> [--max 6] [--seed 0] [--dmin 0.02] [--dmax 0.04] [--nelect N] [--force]
Summarise the finished runs:
    python make_spin_tests.py <job_dir> --summarise [--json spin_summary.json]

Set-up writes one sub-folder <job_dir>/spin_N<k>/ per parity-allowed fixed total moment k = NUPDOWN
(0, 2, 4, ... for an even electron count, 1, 3, 5, ... for an odd one, up to --max). Each holds:
    POSCAR   the converged CONTCAR with the metal and its 4 ligands (found by geometry) displaced by a random
             vector of length dmin-dmax Å; the same displacements, from --seed, in every sub-folder
    INCAR    the job's INCAR with NUPDOWN = k and ISYM = 0; every other key unchanged
    KPOINTS, POTCAR, the job script: copied unchanged
The seed and displacements are recorded in <job_dir>/spin_tests.json. Existing sub-folders are never
overwritten. Group-3 metals (Sc, Y, La; formally d0) are skipped unless --force.
The electron count comes from NELECT in the job's OUTCAR, or --nelect.

--summarise tabulates, for each spin_N<k> run that converged: E0, ΔE = E0 - E0(original), the final magnetic
moment and the change of each M-X bond length relative to the original CONTCAR.
"""
import argparse
import glob
import json
import math
import os
import random
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import herlib  # noqa: E402
from check_jobs import check  # noqa: E402

D0_METALS = {"Sc", "Y", "La"}
COPY = ("KPOINTS", "POTCAR")
MANIFEST = "spin_tests.json"


def random_vector(rng, dmin, dmax):
    while True:
        v = [rng.gauss(0, 1) for _ in range(3)]
        n = math.sqrt(sum(x * x for x in v))
        if n > 1e-8:
            break
    length = rng.uniform(dmin, dmax)
    return [x / n * length for x in v]


def setup(d, a):
    r = check(d, a.log, not a.no_scheduler)
    if r["status"] != "converged":
        raise SystemExit(f"{d}: job is {r['status']}, not converged")
    s = herlib.read_structure(os.path.join(d, "CONTCAR"))
    m, ligs = herlib.find_site(s)
    metal = s["elements"][m]
    if metal in D0_METALS and not a.force:
        raise SystemExit(f"{d}: {metal} is formally d0; spin test skipped (use --force to override)")
    nelect = a.nelect if a.nelect is not None else herlib.outcar_nelect(os.path.join(d, "OUTCAR"))
    if nelect is None:
        raise SystemExit(f"{d}: NELECT not found in OUTCAR; give --nelect")
    if abs(nelect - round(nelect)) > 1e-6:
        raise SystemExit(f"{d}: NELECT = {nelect} is not an integer; parity undefined")
    last = herlib.oszicar_last(os.path.join(d, "OSZICAR"))
    if last and last["mag"] is not None and abs(last["mag"]) < 0.1:
        print(f"note: {d} converged closed-shell (mag = {last['mag']:.3f}); spin test is meant for open-shell cases")
    moments = list(range(int(round(nelect)) % 2, a.max + 1, 2))

    rng = random.Random(a.seed)
    disp = []
    inv = herlib.inverse(s["lattice"])
    for i in [m] + ligs:
        v = random_vector(rng, a.dmin, a.dmax)
        dfrac = herlib.mat_vec_t(inv, v)
        s["frac"][i] = [s["frac"][i][k] + dfrac[k] for k in range(3)]
        disp.append({"atom": f"{s['elements'][i]}{i + 1}", "vector": [round(x, 6) for x in v],
                     "length": round(math.sqrt(sum(x * x for x in v)), 6)})
    script = a.script
    made = []
    for k in moments:
        sub = os.path.join(d, f"spin_N{k}")
        if os.path.exists(sub):
            print(f"KEEP {sub} (exists, not overwritten)")
            continue
        os.makedirs(sub)
        herlib.write_structure(s, os.path.join(sub, "POSCAR"), f"{s['comment']} | spin test NUPDOWN={k}")
        shutil.copyfile(os.path.join(d, "INCAR"), os.path.join(sub, "INCAR"))
        herlib.set_incar(os.path.join(sub, "INCAR"), {"NUPDOWN": k, "ISYM": 0})
        for f in COPY + (script,):
            if os.path.isfile(os.path.join(d, f)):
                shutil.copyfile(os.path.join(d, f), os.path.join(sub, f))
        made.append(sub)
        print(f"WROTE {sub}  (NUPDOWN = {k})")
    manifest = {"source": d, "metal": metal, "nelect": nelect, "moments": moments, "seed": a.seed,
                "dmin": a.dmin, "dmax": a.dmax, "displacements": disp, "created": herlib.now()}
    if made:
        path = os.path.join(d, MANIFEST)
        old = []
        if os.path.isfile(path):
            with open(path) as fh:
                old = json.load(fh)
        with open(path, "w") as fh:
            json.dump((old if isinstance(old, list) else [old]) + [manifest], fh, indent=1)
    return made


def bonds(path):
    s = herlib.read_structure(path)
    m, ligs = herlib.find_site(s)
    return {f"{s['elements'][j]}{j + 1}": herlib.distance(s, m, j) for j in ligs}


def summarise(d, a):
    ref = herlib.oszicar_last(os.path.join(d, "OSZICAR"))
    if ref is None:
        raise SystemExit(f"{d}: no OSZICAR")
    ref_bonds = bonds(os.path.join(d, "CONTCAR"))
    rows = []
    subs = []
    for p in glob.glob(os.path.join(d, "spin_N*")):
        m = re.fullmatch(r"spin_N(\d+)", os.path.basename(p))
        if m and os.path.isdir(p):
            subs.append((int(m.group(1)), p))
    for k, sub in sorted(subs):
        st = check(sub, a.log, not a.no_scheduler)["status"]
        row = {"NUPDOWN": k, "status": st, "E0": None, "dE": None, "mag": None, "bond_change": {}}
        if st == "converged":
            last = herlib.oszicar_last(os.path.join(sub, "OSZICAR"))
            if last is not None:
                row.update(E0=last["E0"], dE=last["E0"] - ref["E0"], mag=last["mag"])
            b = bonds(os.path.join(sub, "CONTCAR"))
            row["bond_change"] = {k2: round(b[k2] - v, 4) if k2 in b else None for k2, v in ref_bonds.items()}
        rows.append(row)
    print(f"original: E0 = {ref['E0']:.5f} eV  mag = {ref['mag']}")
    print(f"{'NUPDOWN':>7} {'status':<12} {'dE (eV)':>9} {'mag':>7}  M-X change (Å)")
    for r in rows:
        dE = "" if r["dE"] is None else f"{r['dE']:+.4f}"
        mag = "" if r["mag"] is None else f"{r['mag']:.3f}"
        bc = " ".join(f"{k}:{v:+.3f}" for k, v in r["bond_change"].items() if v is not None)
        print(f"{r['NUPDOWN']:>7} {r['status']:<12} {dE:>9} {mag:>7}  {bc}")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump({"source": d, "E0": ref["E0"], "mag": ref["mag"], "bonds": ref_bonds, "runs": rows}, fh,
                      indent=1)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dir", help="converged job folder")
    ap.add_argument("--summarise", action="store_true", help="tabulate finished spin_N* runs")
    ap.add_argument("--max", type=int, default=6, help="largest NUPDOWN (default 6)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--dmin", type=float, default=0.02, help="smallest displacement in Å (default 0.02)")
    ap.add_argument("--dmax", type=float, default=0.04, help="largest displacement in Å (default 0.04)")
    ap.add_argument("--nelect", type=float, help="electron count if OUTCAR has no NELECT")
    ap.add_argument("--force", action="store_true", help="also for group-3 (formally d0) metals")
    ap.add_argument("--script", default="sub.sh", help="job script to copy")
    ap.add_argument("--log", default="log", help="VASP stdout file written by the job script")
    ap.add_argument("--no-scheduler", action="store_true", help="do not query qstat/squeue")
    ap.add_argument("--json", help="--summarise: also write the table to this file")
    a = ap.parse_args()
    if a.summarise:
        summarise(a.dir, a)
    else:
        setup(a.dir, a)


if __name__ == "__main__":
    main()
