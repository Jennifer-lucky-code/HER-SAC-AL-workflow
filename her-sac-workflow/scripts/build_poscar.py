#!/usr/bin/env python3
"""Build POSCARs for candidates by element substitution on a template TMX4@G structure.

Template: a relaxed or ideal POSCAR (VASP 5 format) of the graphene divacancy with a metal at the centre and
four coordinating atoms. Site map (JSON) gives 1-based atom indices in that template:
    {"metal": 61, "ligands": [52, 53, 70, 71]}
The ligand order MUST follow the same convention as the screening program's 配体 list [X1, X2, X3, X4];
confirm it with the user once and keep the site map with the project.

Usage:
    python build_poscar.py candidates.csv --template POSCAR_template --sites sites.json [--only Ru-PONS@G ...]
Writes <name>/POSCAR for each candidate. Existing POSCARs are never overwritten.
"""
import argparse
import csv
import json
import os

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se "
           "Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy "
           "Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg").split()


def read_poscar(path):
    lines = open(path).read().splitlines()
    species, counts = lines[5].split(), [int(x) for x in lines[6].split()]
    i = 7
    selective = lines[i].strip().lower().startswith("s")
    if selective:
        i += 1
    coord_type = lines[i].strip()
    i += 1
    elems = [s for s, n in zip(species, counts) for _ in range(n)]
    atoms = []
    for k in range(len(elems)):
        parts = lines[i + k].split()
        atoms.append({"el": elems[k], "xyz": parts[:3], "flags": parts[3:6] if selective else []})
    return {"head": lines[:5], "selective": selective, "coord_type": coord_type, "atoms": atoms}


def write_poscar(p, path, comment):
    order = []
    for a in p["atoms"]:
        if a["el"] not in order:
            order.append(a["el"])
    grouped = [a for el in order for a in p["atoms"] if a["el"] == el]
    out = [comment] + p["head"][1:5]
    out.append("  ".join(order))
    out.append("  ".join(str(sum(a["el"] == el for a in p["atoms"])) for el in order))
    if p["selective"]:
        out.append("Selective dynamics")
    out.append(p["coord_type"])
    for a in grouped:
        out.append("  ".join(a["xyz"] + a["flags"]))
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidates")
    ap.add_argument("--template", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args()
    sites = json.load(open(a.sites))
    with open(a.candidates, encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r.get("Z")]
    if a.only:
        rows = [r for r in rows if r["name"] in a.only]
        missing = set(a.only) - {r["name"] for r in rows}
        for name in sorted(missing):
            print(f"NOT FOUND {name} (not in {a.candidates})")
    for r in rows:
        p = read_poscar(a.template)
        lig = [int(x) for x in r["ligands"].split()]
        if len(lig) != len(sites["ligands"]):
            print(f"SKIP {r['name']}: {len(lig)} ligands vs {len(sites['ligands'])} sites")
            continue
        p["atoms"][sites["metal"] - 1]["el"] = SYMBOLS[int(r["Z"]) - 1]
        for idx, z in zip(sites["ligands"], lig):
            p["atoms"][idx - 1]["el"] = SYMBOLS[z - 1]
        os.makedirs(r["name"], exist_ok=True)
        out = os.path.join(r["name"], "POSCAR")
        if os.path.exists(out):
            print(f"KEEP {out} (exists, not overwritten)")
            continue
        write_poscar(p, out, r["name"])
        print(f"WROTE {out}")


if __name__ == "__main__":
    main()
