#!/usr/bin/env python3
"""Collect energies from finished job folders into energies.csv (the input of her_analysis.py).

No number is typed by hand: E comes from OSZICAR, ZPE/TS from the saved vaspkit 501 output, G_H2 from the
reference-energy table in references/vasp_settings.md.

    python collect_energies.py <candidate_dir> [<candidate_dir> ...] [-o energies.csv] [--json energies.json]
                               [--run-vaspkit] [--no-scheduler] [--log log]

Folder convention per candidate (references/dft_her.md §1), all optional except slab/:
    slab/                      bare catalyst *            G_* = E_*
    H1_metal/, H1_lig<k>/      H* on the metal or on ligand X<k> (k = 1-4, the X1..X4 order of the name)
    H2_C/, H2_M/               2H*: second H on an adjacent C / both H on the metal ("horns")
    freq_<state>/              frequency job of each adsorbed state, e.g. freq_H1_metal/, freq_H2_C/,
                               holding vaspkit_501.txt = output of `printf "501\\n298.15\\n" | vaspkit`

Rules:
    E   = E0 on the last OSZICAR line, only if check_jobs reports the relaxation converged.
    G   = E + ZPE - TS for adsorbed states; ZPE/TS only from a converged frequency job.
    H*  = the H1_* site with the lowest G; the site is written to H_site.
    A missing or unconverged value leaves the cell blank and adds a note. Nothing is estimated.
--run-vaspkit runs vaspkit 501 in a converged frequency folder without vaspkit_501.txt (if vaspkit is on PATH)
and saves its output there.
"""
import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import herlib  # noqa: E402
from check_jobs import check as job_check  # noqa: E402

VASPKIT_FILE = "vaspkit_501.txt"
FIELDS = ["name", "E_slab", "G_H", "G_2H_C", "G_2H_M", "G_H2", "H_site", "notes"]


def energy(d, notes, logname, use_scheduler):
    """E0 of a converged relaxation, else None (with a note)."""
    tag = os.path.basename(os.path.normpath(d))
    if not os.path.isdir(d):
        notes.append(f"{tag}: no folder")
        return None
    st = job_check(d, logname, use_scheduler)["status"]
    if st != "converged":
        notes.append(f"{tag}: {st}")
        return None
    last = herlib.oszicar_last(os.path.join(d, "OSZICAR"))
    if last is None:
        notes.append(f"{tag}: no ionic step in OSZICAR")
        return None
    return last["E0"]


def zpe_ts(d, notes, logname, use_scheduler, run_vaspkit=False):
    tag = os.path.basename(os.path.normpath(d))
    if not os.path.isdir(d):
        notes.append(f"{tag}: no folder")
        return None
    st = job_check(d, logname, use_scheduler)["status"]
    if st != "converged":
        notes.append(f"{tag}: {st}")
        return None
    path = os.path.join(d, VASPKIT_FILE)
    if not os.path.isfile(path) and run_vaspkit:
        if shutil.which("vaspkit"):
            r = subprocess.run(["vaspkit"], input=f"501\n{herlib.TEMPERATURE}\n", capture_output=True,
                               text=True, cwd=d)
            with open(path, "w") as fh:
                fh.write(r.stdout)
        else:
            notes.append(f"{tag}: vaspkit not on PATH")
    if not os.path.isfile(path):
        notes.append(f"{tag}: no {VASPKIT_FILE}")
        return None
    with open(path, errors="replace") as fh:
        v = herlib.parse_vaspkit_501(fh.read())
    if v["ZPE"] is None or v["TS"] is None:
        notes.append(f"{tag}: ZPE/TS lines missing in {VASPKIT_FILE}")
        return None
    return v


def free_energy(cand, state, notes, logname, use_scheduler, run_vaspkit):
    """{'E', 'ZPE', 'TS', 'G'} of an adsorbed state; missing parts are None."""
    E = energy(os.path.join(cand, state), notes, logname, use_scheduler)
    v = zpe_ts(os.path.join(cand, "freq_" + state), notes, logname, use_scheduler, run_vaspkit) \
        if E is not None else None
    out = {"E": E, "ZPE": v and v["ZPE"], "TS": v and v["TS"], "G": None}
    if E is not None and v is not None:
        out["G"] = E + v["ZPE"] - v["TS"]
    return out


def collect(cand, g_h2, logname="log", use_scheduler=True, run_vaspkit=False):
    name = os.path.basename(os.path.normpath(cand))
    notes = []
    detail = {"slab": {"E": energy(os.path.join(cand, "slab"), notes, logname, use_scheduler)}}
    sites = sorted(os.path.basename(p) for p in glob.glob(os.path.join(cand, "H1_*")) if os.path.isdir(p))
    if not sites:
        notes.append("no H1_* folder")
    for s in sites:
        detail[s] = free_energy(cand, s, notes, logname, use_scheduler, run_vaspkit)
    for s in ("H2_C", "H2_M"):
        detail[s] = free_energy(cand, s, notes, logname, use_scheduler, run_vaspkit)

    done = {s: detail[s] for s in sites if detail[s]["G"] is not None}
    best = min(done, key=lambda s: done[s]["G"]) if done else None
    if best:
        for s in sites:
            E = detail[s]["E"]
            if s not in done and E is not None and E < done[best]["E"]:
                notes.append(f"{s} has a lower E than {best} but no G; H* site may change")
    row = {"name": name, "E_slab": detail["slab"]["E"], "G_H": done[best]["G"] if best else None,
           "G_2H_C": detail["H2_C"]["G"], "G_2H_M": detail["H2_M"]["G"], "G_H2": g_h2,
           "H_site": best or "", "notes": "; ".join(notes)}
    if g_h2 is None:
        row["notes"] = "; ".join(notes + ["G_H2 is TODO in vasp_settings.md"])
    return row, detail


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("candidates", nargs="+", help="candidate folders, e.g. Ru-PONS@G")
    ap.add_argument("-o", "--out", default="energies.csv")
    ap.add_argument("--json", help="also write per-state E, ZPE, TS, G to this file")
    ap.add_argument("--run-vaspkit", action="store_true", help="run vaspkit 501 where vaspkit_501.txt is missing")
    ap.add_argument("--log", default="log", help="VASP stdout file written by the job script")
    ap.add_argument("--no-scheduler", action="store_true", help="do not query qstat/squeue")
    ap.add_argument("--settings", default=herlib.SETTINGS_MD, help="vasp_settings.md with the G_H2 value")
    a = ap.parse_args()
    g_h2 = herlib.g_h2(a.settings)
    rows, details = [], {}
    for c in a.candidates:
        row, detail = collect(c, g_h2, a.log, not a.no_scheduler, a.run_vaspkit)
        rows.append(row)
        details[row["name"]] = detail
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    for r in rows:
        print(f"{r['name']:<14} H*={r['H_site'] or '-':<10} {r['notes'] or 'complete'}")
    print(f"Saved {a.out} ({len(rows)} rows)")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(details, fh, indent=1)


if __name__ == "__main__":
    main()
