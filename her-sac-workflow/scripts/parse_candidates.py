#!/usr/bin/env python3
"""Parse the screening program's log into a candidate table.

The screening program itself is private: run it as instructed, never open, copy or describe its source.
This script only reads its log. Expected blocks (any section headers are kept as 'section'):

    Top 10 最优配置 (ΔG ≈ 0)
    配置 14:
      Z=44, 配体=[15, 8, 7, 16]  (PONS)
      ΔG = -0.0244 ± 0.5133

Usage:
    python parse_candidates.py screening.log [-o candidates.csv]
"""
import argparse
import csv
import re

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se "
           "Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy "
           "Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg").split()


def symbol(z):
    return SYMBOLS[z - 1]


# C2v of the X1 X2 / X3 X4 site rectangle (see references/screening.md): identity, C2, two mirrors.
C2V = ((0, 1, 2, 3), (3, 2, 1, 0), (1, 0, 3, 2), (2, 3, 0, 1))


def canonical(label):
    """Symmetry-equivalent ligand strings map to one key (lexicographic minimum over C2v)."""
    return min("".join(label[i] for i in perm) for perm in C2V)


RE_CONF = re.compile(r"配置\s*(\d+)\s*[:：]")
RE_ZL = re.compile(r"Z\s*=\s*(\d+)\s*,\s*配体\s*=\s*\[([^\]]+)\]\s*(?:\(\s*([A-Za-z]+)\s*\))?")
RE_DG = re.compile(r"ΔG\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*±\s*(\d*\.?\d+(?:[eE][-+]?\d+)?)")


def parse(text):
    rows, cur, section = [], None, ""
    for line in text.splitlines():
        s = line.strip()
        if not s or set(s) <= {"="}:
            continue
        m = RE_CONF.search(s)
        if m:
            cur = {"section": section, "config_id": m.group(1)}
            rows.append(cur)
            continue
        m = RE_ZL.search(s)
        if m and cur is not None:
            z = int(m.group(1))
            lig = [int(x) for x in m.group(2).split(",")]
            label = "".join(symbol(x) for x in lig)
            cur.update(Z=z, metal=symbol(z), ligands=" ".join(map(str, lig)),
                       label=label, name=f"{symbol(z)}-{label}@G",
                       canonical=f"{symbol(z)}-{canonical(label)}@G" if len(lig) == 4 else "",
                       label_check="ok" if (m.group(3) or label) == label else f"MISMATCH log={m.group(3)}")
            continue
        m = RE_DG.search(s)
        if m and cur is not None:
            cur.update(dG=float(m.group(1)), sigma=float(m.group(2)))
            continue
        if cur is None or "Z" in cur:
            section = s  # any other non-empty line is treated as a section header
            cur = None
    for r in rows:
        r.setdefault("dG", "")
        r.setdefault("sigma", "")
        if "Z" not in r:
            r["label_check"] = "INCOMPLETE block"
        elif r["dG"] == "":
            r["label_check"] = ("MISSING dG; " + r["label_check"]) if r["label_check"] != "ok" else "MISSING dG"
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("log")
    p.add_argument("-o", "--out", default="candidates.csv")
    a = p.parse_args()
    rows = parse(open(a.log, encoding="utf-8").read())
    fields = ["section", "config_id", "name", "canonical", "Z", "metal", "ligands", "label", "dG", "sigma", "label_check"]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        dg = "missing" if r["dG"] == "" else f"{r['dG']:+.4f} ± {r['sigma']:.4f}"
        print(f"[{r['section'][:30]}] #{r['config_id']:<4} {r.get('name', '?'):<14} ΔG = {dg}  {r['label_check']}")
    print(f"Saved {a.out} ({len(rows)} candidates)")


if __name__ == "__main__":
    main()
