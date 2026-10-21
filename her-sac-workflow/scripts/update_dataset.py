#!/usr/bin/env python3
"""Append validated results to the DFT dataset (loop step 8, references/dft_her.md §9). Never overwrites a row.

    python update_dataset.py her_results.csv [--dataset dft_dataset.csv] [--root .] [--extra features.csv]
                             [--structure structure_check.json] [--state loop_state.json]
                             [--event-log workflow_log.jsonl] [--allow-flagged] [--date YYYY-MM-DD]
                             [--no-scheduler] [--dry-run] [--json report.json]

For each row of her_results.csv (from her_analysis.py), the jobs that produced it are
<root>/<name>/{slab, <H_site>, H2_C, H2_M} and their freq_* folders (those that exist).
Every one of them must pass check_structure.py: results are taken from --structure (the --json output of
check_structure.py) or, if a job is not in that file, computed here with the default thresholds.

    - a job that is not converged or has no check -> row refused
    - a job with a structure flag -> row refused, unless --allow-flagged (the flags go into notes)
    - deduplication by the canonical (C2v) key of parse_candidates.py: an existing row with the same values is
      left alone; different values -> a new row whose notes say which row it supersedes
    - provenance: job paths (relative to the folder of loop_state.json), date, retry count (loop_state.json), INCAR keys that differ from the templates in
      references/vasp_settings.md (SYSTEM, MAGMOM, NPAR ignored)
--extra is an optional CSV keyed by name with Ef and/or d_band_centre (blank = not available).
Each added row is logged as an 'added' event and its name is appended to "added" in loop_state.json.
"""
import argparse
import csv
import datetime
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import herlib  # noqa: E402
from check_structure import check_dir  # noqa: E402
from parse_candidates import canonical  # noqa: E402

FIELDS = ["name", "canonical", "TM", "ligands", "Ef", "dG_H", "H_site", "dG_2H_C", "dG_2H_M", "eta_VH",
          "eta_VTC", "eta_VTM", "pathway", "limiting_step", "eta_limiting", "eta_screen", "mag_slab", "mag_H",
          "mag_2H_C", "mag_2H_M", "d_band_centre", "job_paths", "date", "retries", "incar_diff", "notes"]
# columns compared to decide whether an existing row already holds the same result
VALUES = ["Ef", "dG_H", "H_site", "dG_2H_C", "dG_2H_M", "eta_VH", "eta_VTC", "eta_VTM", "pathway",
          "limiting_step", "eta_limiting", "eta_screen", "mag_slab", "mag_H", "mag_2H_C", "mag_2H_M",
          "d_band_centre"]
FROM_RESULTS = {"dG_H": "dG_H", "H_site": "H_site", "dG_2H_C": "dG_2H_VTC", "dG_2H_M": "dG_2H_VTM",
                "eta_VH": "eta_VH", "eta_VTC": "eta_VTC", "eta_VTM": "eta_VTM", "pathway": "pathway",
                "limiting_step": "limiting_step", "eta_limiting": "eta_limiting", "eta_screen": "eta_screen",
                "Ef": "Ef"}


def split_name(name):
    """'Ru-PONS@G' -> ('Ru', 'PONS')."""
    tm, rest = name.split("-", 1)
    return tm, rest.split("@", 1)[0]


def _fmt(v, nd=4):
    if v is None or v == "":
        return ""
    try:
        return f"{float(v):.{nd}f}"
    except ValueError:
        return str(v)


def job_dirs(cand, h_site):
    out = []
    for state in ["slab", h_site, "H2_C", "H2_M"]:
        if state and os.path.isdir(os.path.join(cand, state)):
            out.append(state)
            if state != "slab" and os.path.isdir(os.path.join(cand, "freq_" + state)):
                out.append("freq_" + state)
    return out


def build_row(res, root, checks, state, state_path, date, allow_flagged, opts):
    """(row or None, reasons) for one her_results row."""
    name = res["name"]
    cand = os.path.join(root, name)
    if not os.path.isdir(cand):
        return None, [f"no folder {cand}"]
    tm, label = split_name(name)
    if not res.get("H_site"):
        return None, ["no H_site in her_results.csv (run collect_energies.py)"]
    jobs = job_dirs(cand, res["H_site"])
    if "slab" not in jobs or res["H_site"] not in jobs:
        return None, ["slab/ or the H* folder is missing"]
    refuse, notes, mags, diffs, retries = [], [], {}, [], 0
    for j in jobs:
        d = os.path.join(cand, j)
        r = checks.get(os.path.abspath(d)) or check_dir(d, None, opts["logname"], opts["use_scheduler"])
        if r["status"] == "not_judged":
            refuse.append(f"{j}: " + "; ".join(r["flags"]))
        elif r["status"] == "flag":
            (notes if allow_flagged else refuse).append(f"{j}: structure flag: " + "; ".join(r["flags"]))
        elif r["status"] != "pass":
            refuse.append(f"{j}: no structure check")
        mags[j] = r.get("mag")
        inc = herlib.parse_incar(os.path.join(d, "INCAR"))
        kind = "frequency" if j.startswith("freq_") else "relaxation"
        diffs += [f"{j}:{x}" for x in herlib.incar_diff(inc, herlib.incar_template(kind, opts.get("settings")))]
        retries += state["jobs"].get(herlib.job_key(d, state_path), {}).get("retries", 0)
    if refuse:
        return None, refuse
    row = {k: _fmt(res.get(src)) if k not in ("H_site", "pathway", "limiting_step") else res.get(src, "")
           for k, src in FROM_RESULTS.items()}
    row.update(name=name, canonical=f"{tm}-{canonical(label)}@G" if len(label) == 4 else "", TM=tm,
               ligands=label, mag_slab=_fmt(mags.get("slab"), 3), mag_H=_fmt(mags.get(res["H_site"]), 3),
               mag_2H_C=_fmt(mags.get("H2_C"), 3), mag_2H_M=_fmt(mags.get("H2_M"), 3),
               job_paths=";".join(herlib.job_key(os.path.join(cand, j), state_path) for j in jobs), date=date, retries=str(retries),
               incar_diff="; ".join(diffs), d_band_centre="")
    if res.get("notes"):
        notes.insert(0, res["notes"])
    row["notes"] = "; ".join(notes)
    return row, []


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("results", help="her_results.csv from her_analysis.py")
    ap.add_argument("--dataset", default="dft_dataset.csv")
    ap.add_argument("--root", default=".", help="folder that holds the <name>/ candidate folders")
    ap.add_argument("--extra", help="CSV keyed by name with Ef and/or d_band_centre")
    ap.add_argument("--structure", help="--json output of check_structure.py")
    ap.add_argument("--state", default=herlib.STATE_FILE)
    ap.add_argument("--event-log", default=herlib.LOG_FILE)
    ap.add_argument("--allow-flagged", action="store_true", help="add rows with structure flags (noted)")
    ap.add_argument("--date", default=datetime.date.today().isoformat())
    ap.add_argument("--log", default="log", help="VASP stdout file written by the job script")
    ap.add_argument("--no-scheduler", action="store_true", help="do not query qstat/squeue")
    ap.add_argument("--settings", default=herlib.SETTINGS_MD)
    ap.add_argument("--dry-run", action="store_true", help="report only; write nothing")
    ap.add_argument("--json", help="write the per-row decisions to this file")
    a = ap.parse_args()

    checks = {}
    if a.structure:
        with open(a.structure) as fh:
            checks = {os.path.abspath(r["dir"]): r for r in json.load(fh)}
    extra = {}
    if a.extra:
        with open(a.extra, newline="", encoding="utf-8") as fh:
            extra = {r["name"]: r for r in csv.DictReader(fh)}
    existing = []
    if os.path.isfile(a.dataset):
        with open(a.dataset, newline="", encoding="utf-8") as fh:
            rd = csv.DictReader(fh)
            if rd.fieldnames and rd.fieldnames != FIELDS:
                raise SystemExit(f"{a.dataset}: header differs from the dataset columns; not touching it.\n"
                                 f"expected: {','.join(FIELDS)}")
            existing = list(rd)
    state = herlib.load_state(a.state)
    opts = {"logname": a.log, "use_scheduler": not a.no_scheduler}

    with open(a.results, newline="", encoding="utf-8") as fh:
        results = list(csv.DictReader(fh))
    report, new = [], []
    for res in results:
        ex = extra.get(res["name"], {})
        if ex.get("Ef"):
            res["Ef"] = ex["Ef"]
        row, why = build_row(res, a.root, checks, state, a.state, a.date, a.allow_flagged,
                             dict(opts, settings=a.settings))
        if row is None:
            report.append({"name": res["name"], "action": "refused", "reasons": why})
            continue
        row["d_band_centre"] = _fmt(ex.get("d_band_centre"))
        same = [(n, r) for n, r in enumerate(existing + new, 1) if r["canonical"] == row["canonical"]]
        if same:
            n, last = same[-1]
            if all(last.get(k, "") == row[k] for k in VALUES):
                report.append({"name": res["name"], "action": "unchanged", "reasons": [f"same as row {n}"]})
                continue
            row["notes"] = "; ".join(x for x in (f"supersedes row {n} ({last['name']}, {last['date']})",
                                                 row["notes"]) if x)
        new.append(row)
        report.append({"name": res["name"], "action": "added", "reasons": [row["notes"]] if row["notes"] else []})

    for r in report:
        print(f"{r['action']:<10} {r['name']:<14} {'; '.join(r['reasons'])}")
    if not a.dry_run and new:
        write_header = not os.path.isfile(a.dataset) or os.path.getsize(a.dataset) == 0
        with open(a.dataset, "a", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            if write_header:
                w.writeheader()
            w.writerows(new)
        for row in new:
            herlib.log_event(a.event_log, "added", job_dir=herlib.job_key(os.path.join(a.root, row["name"]), a.state),
                             change=row["notes"], name=row["name"], canonical=row["canonical"])
            state.setdefault("added", []).append(row["name"])
        state["dataset_size"] = len(existing) + len(new)
        herlib.save_state(state, a.state)
        print(f"Appended {len(new)} row(s) to {a.dataset}")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(report, fh, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
