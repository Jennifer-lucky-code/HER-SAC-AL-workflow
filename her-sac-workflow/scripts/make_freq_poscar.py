#!/usr/bin/env python3
"""Make a frequency (ZPE/TS) POSCAR from a relaxed CONTCAR: substrate fixed, adsorbates free.

Relaxations use no fixed atoms. Only frequency jobs fix atoms: every substrate atom gets F F F and the
adsorbate atoms (default: all H) get T T T.

Usage:
    python make_freq_poscar.py CONTCAR [-o POSCAR_freq] [--free H] [--free-idx 32 33]
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("contcar")
    ap.add_argument("-o", "--out", default="POSCAR_freq")
    ap.add_argument("--free", nargs="*", default=["H"], help="elements left free")
    ap.add_argument("--free-idx", nargs="*", type=int, default=[], help="extra 1-based atom indices left free")
    a = ap.parse_args()
    L = open(a.contcar).read().splitlines()
    species, counts = L[5].split(), [int(x) for x in L[6].split()]
    elems = [s for s, n in zip(species, counts) for _ in range(n)]
    i = 7 + (1 if L[7].strip().lower().startswith("s") else 0)
    coord_type, start = L[i].strip(), i + 1
    out = L[:7] + ["Selective dynamics", coord_type]
    n_free = 0
    for k, el in enumerate(elems):
        xyz = L[start + k].split()[:3]
        free = el in a.free or (k + 1) in a.free_idx
        n_free += free
        out.append("  ".join(xyz + (["T", "T", "T"] if free else ["F", "F", "F"])))
    open(a.out, "w").write("\n".join(out) + "\n")
    print(f"Wrote {a.out}: {n_free} free atoms, {len(elems) - n_free} fixed")
    if n_free == 0:
        print("WARNING: no free atoms; check --free / --free-idx")


if __name__ == "__main__":
    main()
