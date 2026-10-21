#!/usr/bin/env python3
"""Structure checks of references/dft_her.md §6 for finished VASP job folders: pass or flag (with reasons).

A structure is judged only when check_jobs.check() reports the job as converged; otherwise the result is
'not_judged'. Reads CONTCAR and OSZICAR, plus OUTCAR for frequency jobs (IBRION 5-8).

    python check_structure.py <job_dir> [<job_dir> ...] [--json structure_check.json] [--no-scheduler]
                              [--bond-factor 1.15] [--jt-spread 0.10] [--max-out-of-plane 0.5]
                              [--h-bond-factor 1.30] [--h-max-height 4.0] [--hh-min 0.9] [--imag-tol 0]

Atoms are found by geometry (minimum-image distances), never by template index:
metal = the single transition-metal atom; ligands = its 4 nearest non-H neighbours.

Checks (every threshold is an option; defaults in brackets):
    coordination   the number of non-H atoms bonded to the metal must be 4. Bonded = distance
                   <= bond_factor x (r_cov(M) + r_cov(X)), Cordero covalent radii [1.15].
                   Reports the 4 M-X bond lengths and their spread (max - min).
    Jahn-Teller    spread > jt_spread [0.10 Å] -> warning "possible Jahn–Teller distortion". A warning only:
                   it never changes the status. Mixed-ligand sites have an intrinsic spread, so the spread over
                   same-element ligands is reported too.
    out of plane   distance of the metal from the least-squares plane of its 4 ligands > max_out_of_plane [0.5 Å]
    H bound        every H has a non-H neighbour within h_bond_factor x (r_cov(H) + r_cov(X)) [1.30]
    H in vacuum    every H lies within h_max_height [4.0 Å] of the mean height of the slab (non-H atoms)
    H-H            any H-H distance < hh_min [0.9 Å] -> "H₂ / dihydrogen formed"
    magnetic moment  final value from the last mag= in OSZICAR (reported, not judged)
    imaginary modes  frequency jobs: every f/i= mode in OUTCAR above imag_tol [0 cm-1] is listed and flagged

Status: pass | flag | not_judged. Exit code 0 unless a folder could not be read at all.
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import herlib  # noqa: E402
from check_jobs import check as job_check  # noqa: E402

DEFAULTS = {"bond_factor": 1.15, "jt_spread": 0.10, "max_out_of_plane": 0.5, "h_bond_factor": 1.30,
            "h_max_height": 4.0, "hh_min": 0.9, "imag_tol": 0.0}


def _radius(el):
    if el not in herlib.COVALENT_RADII:
        raise ValueError(f"no covalent radius for {el}; add it to herlib.COVALENT_RADII")
    return herlib.COVALENT_RADII[el]


def plane_distance(points, p):
    """Distance of p from the least-squares plane z = a x + b y + c through points (slab normal ~ z)."""
    n = len(points)
    sx = sum(q[0] for q in points); sy = sum(q[1] for q in points); sz = sum(q[2] for q in points)
    sxx = sum(q[0] ** 2 for q in points); syy = sum(q[1] ** 2 for q in points)
    sxy = sum(q[0] * q[1] for q in points)
    sxz = sum(q[0] * q[2] for q in points); syz = sum(q[1] * q[2] for q in points)
    A = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, n]]
    inv = herlib.inverse(A)
    a, b, c = (sum(inv[i][j] * v for j, v in enumerate((sxz, syz, sz))) for i in range(3))
    return (p[2] - (a * p[0] + b * p[1] + c)) / math.sqrt(1 + a * a + b * b)


def check_geometry(s, opts=None):
    """Geometry checks on a structure dict (herlib.read_structure). Returns (report, flags, warnings)."""
    o = dict(DEFAULTS, **(opts or {}))
    rep, flags, warns = {}, [], []
    el = s["elements"]
    try:
        m, ligs = herlib.find_site(s)
    except ValueError as e:
        return rep, [str(e)], warns
    rep["metal"] = f"{el[m]}{m + 1}"
    bonds = [herlib.distance(s, m, j) for j in ligs]
    rep["ligands"] = [{"atom": f"{el[j]}{j + 1}", "d": round(d, 4)} for j, d in zip(ligs, bonds)]
    rep["bond_spread"] = round(max(bonds) - min(bonds), 4)
    same = {}
    for j, d in zip(ligs, bonds):
        same.setdefault(el[j], []).append(d)
    rep["bond_spread_same_element"] = round(max((max(v) - min(v) for v in same.values()), default=0.0), 4)

    bonded = [(d, j) for d, j in herlib.neighbours(s, m, exclude_h=True)
              if d <= o["bond_factor"] * (_radius(el[m]) + _radius(el[j]))]
    rep["coordination"] = len(bonded)
    if len(bonded) != 4:
        flags.append(f"metal is {len(bonded)}-coordinated (expected 4; bond_factor {o['bond_factor']})")
    if rep["bond_spread"] > o["jt_spread"]:
        warns.append(f"possible Jahn–Teller distortion (M–X spread {rep['bond_spread']:.3f} Å > {o['jt_spread']})")

    origin = s["cart"][m]
    lig_pts = [[origin[k] + v[k] for k in range(3)] for v in (herlib.mic_vector(s, m, j) for j in ligs)]
    rep["out_of_plane"] = round(abs(plane_distance(lig_pts, origin)), 4)
    if rep["out_of_plane"] > o["max_out_of_plane"]:
        flags.append(f"metal {rep['out_of_plane']:.3f} Å out of the ligand plane (> {o['max_out_of_plane']})")

    hs = [i for i, e in enumerate(el) if e == "H"]
    heavy = [i for i, e in enumerate(el) if e != "H"]
    c = s["lattice"][2][2]
    zs = [s["frac"][i][2] for i in heavy]
    # mean slab height along c, unwrapped around the first heavy atom
    zmean = sum(zs[0] + ((z - zs[0]) - round(z - zs[0])) for z in zs) / len(zs)
    rep["H"] = []
    for i in hs:
        d, j = herlib.neighbours(s, i, exclude_h=True)[0]
        dz = s["frac"][i][2] - zmean
        height = abs(dz - round(dz)) * c
        rep["H"].append({"atom": f"H{i + 1}", "nearest": f"{el[j]}{j + 1}", "d": round(d, 4),
                         "height": round(height, 3)})
        if d > o["h_bond_factor"] * (_radius("H") + _radius(el[j])):
            flags.append(f"H{i + 1} not bound (nearest non-H {el[j]}{j + 1} at {d:.3f} Å)")
        if height > o["h_max_height"]:
            flags.append(f"H{i + 1} drifted into the vacuum ({height:.2f} Å from the slab plane)")
    hh = [(herlib.distance(s, a, b), a, b) for k, a in enumerate(hs) for b in hs[k + 1:]]
    rep["HH_min"] = round(min(hh)[0], 4) if hh else None
    for d, a, b in hh:
        if d < o["hh_min"]:
            flags.append(f"H₂ / dihydrogen formed (H{a + 1}–H{b + 1} {d:.3f} Å < {o['hh_min']})")
    return rep, flags, warns


def check_dir(d, opts=None, logname="log", use_scheduler=True):
    o = dict(DEFAULTS, **(opts or {}))
    res = {"dir": d, "status": "not_judged", "job_status": "", "flags": [], "warnings": [], "mag": None,
           "imaginary": None}
    job = job_check(d, logname, use_scheduler)
    res["job_status"] = job["status"]
    if job["status"] != "converged":
        res["flags"].append(f"job not converged ({job['status']}); structure not judged")
        return res
    contcar = os.path.join(d, "CONTCAR")
    if not os.path.isfile(contcar) or os.path.getsize(contcar) == 0:
        res.update(status="flag", flags=["CONTCAR missing or empty"])
        return res
    rep, flags, warns = check_geometry(herlib.read_structure(contcar), o)
    res.update(rep)
    res["flags"] += flags
    res["warnings"] += warns
    last = herlib.oszicar_last(os.path.join(d, "OSZICAR"))
    res["mag"] = None if last is None else last["mag"]
    if res["mag"] is None:
        res["warnings"].append("no mag= in OSZICAR")
    if herlib.is_frequency(herlib.parse_incar(os.path.join(d, "INCAR"))):
        modes = herlib.outcar_imaginary(os.path.join(d, "OUTCAR"))
        if modes is None:
            res["flags"].append("frequency job without OUTCAR")
        else:
            res["imaginary"] = [x for x in modes if x["cm-1"] > o["imag_tol"]]
            if res["imaginary"]:
                res["flags"].append("imaginary modes: " + ", ".join(f"{x['cm-1']:.1f}i cm-1"
                                                                     for x in res["imaginary"]))
    res["status"] = "flag" if res["flags"] else "pass"
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--json", help="also write the results to this file")
    ap.add_argument("--log", default="log", help="VASP stdout file written by the job script")
    ap.add_argument("--no-scheduler", action="store_true", help="do not query qstat/squeue")
    for k, v in DEFAULTS.items():
        ap.add_argument("--" + k.replace("_", "-"), type=float, default=v, help=f"default {v}")
    a = ap.parse_args()
    opts = {k: getattr(a, k) for k in DEFAULTS}
    results = [check_dir(d, opts, a.log, not a.no_scheduler) for d in a.dirs]
    for r in results:
        bonds = " ".join(f"{x['atom']}:{x['d']:.3f}" for x in r.get("ligands", []))
        mag = "" if r["mag"] is None else f"  mag={r['mag']:.3f}"
        print(f"{r['status']:<10} {r['dir']}  {r.get('metal', '')} {bonds}{mag}")
        for f in r["flags"]:
            print(f"    FLAG  {f}")
        for w in r["warnings"]:
            print(f"    note  {w}")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(results, fh, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
