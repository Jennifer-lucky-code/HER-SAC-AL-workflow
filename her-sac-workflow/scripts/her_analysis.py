#!/usr/bin/env python3
"""HER free-energy and mechanism analysis for single-atom catalysts.

Input CSV columns (one row per configuration; energies in eV):
    name, E_slab, G_H, G_2H_C, G_2H_M, G_H2
      E_slab : E_DFT of the bare catalyst (G_* = E_*)
      G_H    : G of H* (E_DFT + ZPE - TS), lowest-energy site
      G_2H_C : G of 2H* with the second H on a ligand atom X, any element (blank if not computed)
      G_2H_M : G of 2H* with both H on the same metal atom, "horns" (blank if not computed)
      G_H2   : G of gas-phase H2 (same reference for every row)
Each 2H* configuration gives its own Volmer-Tafel route (VT-C, VT-M). The VT route with the
smaller largest-uphill step is compared with VH.
Optional columns kept as-is in the output: H_site, H2_site, Ef, notes
Rows with a blank E_slab, G_H or G_H2 (e.g. from collect_energies.py when a job is missing) are skipped
and reported, never filled in.

Reference state: * + 2H+ + 2e- = 0 eV at U = 0 V (CHE model).

Usage:
    python her_analysis.py energies.csv [-o results.csv] [--plot]
"""
import argparse
import csv
import math


def vt_route(dG_H, G_H, G_2H, G_H2):
    dG_2H = G_2H - G_H - 0.5 * G_H2                       # second Volmer (step value)
    dG_T = -(dG_H + dG_2H)                                # Tafel, chemical step
    steps = {"Volmer-1": dG_H, "Volmer-2": dG_2H, "Tafel": dG_T}
    step, eta = max(steps.items(), key=lambda kv: kv[1])
    return {"dG_2H": dG_2H, "level_2H": dG_H + dG_2H, "dG_Tafel": dG_T,
            "eta": max(eta, 0.0), "step": step}


def analyse(E_slab, G_H, G_2H_C, G_2H_M, G_H2):
    dG_H = G_H - E_slab - 0.5 * G_H2                      # first Volmer
    out = {"dG_H": dG_H, "dG_Heyrovsky": -dG_H, "eta_VH": abs(dG_H)}
    routes = {}
    for tag, G2 in (("C", G_2H_C), ("M", G_2H_M)):
        r = vt_route(dG_H, G_H, G2, G_H2) if G2 is not None else None
        routes[tag] = r
        for k in ("dG_2H", "level_2H", "dG_Tafel", "eta"):
            out[f"{k}_VT{tag}" if k != "eta" else f"eta_VT{tag}"] = None if r is None else r[k]
    done = {t: r for t, r in routes.items() if r is not None}
    vh_step = "Volmer" if dG_H > 0 else "Heyrovsky"
    if not done:
        out.update(pathway="VH (2H* missing)", limiting_step=vh_step, eta_limiting=out["eta_VH"],
                   eta_screen=None)
        return out
    best_tag, best = min(done.items(), key=lambda kv: kv[1]["eta"])
    if best["eta"] < out["eta_VH"]:
        out.update(pathway=f"VT-{best_tag}", limiting_step=best["step"], eta_limiting=best["eta"])
    else:
        out.update(pathway="VH", limiting_step=vh_step, eta_limiting=out["eta_VH"])
    # Eq. 4 screening descriptor, using the more stable 2H* configuration
    stable = min(done.values(), key=lambda r: r["level_2H"])
    out["eta_screen"] = max(abs(dG_H), abs(stable["dG_2H"]))
    if len(done) < 2:
        out["pathway"] += " (one 2H* config missing)"
    return out


def _f(x):
    x = "" if x is None else x.strip()
    return "" if x == "" else float(x)


def plot(rows, path="her_diagrams.png"):
    import matplotlib.pyplot as plt
    n = len(rows)
    cols = min(3, n)
    nrows = math.ceil(n / cols)
    fig, axes = plt.subplots(nrows, cols, figsize=(4 * cols, 3.4 * nrows), squeeze=False)
    colors = {"C": "#4FA8E0", "M": "#1F5FA8"}
    labels = {"C": "VT (adjacent C)", "M": "VT (metal, horns)"}
    for ax, r in zip(axes.flat, rows):
        dH = r["dG_H"]

        def bar(x, y, c):
            ax.hlines(y, x - 0.3, x + 0.3, colors=c, lw=3)
        bar(0, 0, "k"); bar(3, 0, "k")
        for tag in ("C", "M"):
            L2 = r.get(f"level_2H_VT{tag}")
            if L2 is None:
                continue
            c = colors[tag]
            ax.plot([1.3, 1.7], [dH, L2], "--", c=c)
            ax.plot([2.3, 2.7], [L2, 0], "--", c=c)
            bar(2, L2, c)
            ax.text(2.35, L2, f"{L2:.2f}", ha="left", va="bottom" if tag == "C" else "top",
                    color=c, fontsize=9)
            ax.plot([], [], c=c, label=labels[tag])
        ax.plot([0.3, 0.7], [0, dH], "--", c="#F28C38")
        ax.plot([1.3, 2.7], [dH, 0], "--", c="#F28C38")
        bar(1, dH, "#F28C38")
        ax.text(1, dH, f"{dH:.2f}", ha="center", va="top", color="#C0601A", fontsize=9)
        ax.plot([], [], c="#F28C38", label="VH")
        ax.set_title(f"{r['name']}  ({r['pathway']})", fontsize=10)
        ax.set_xticks([]); ax.set_ylabel("ΔG (eV)")
        ax.legend(frameon=False, fontsize=7)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=600)
    print(f"Saved {path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv")
    p.add_argument("-o", "--out", default="her_results.csv")
    p.add_argument("--plot", action="store_true")
    a = p.parse_args()

    rows = []
    with open(a.csv, newline="") as fh:
        for raw in csv.DictReader(fh):
            blank = [k for k in ("E_slab", "G_H", "G_H2") if _f(raw.get(k)) == ""]
            if blank:
                print(f"SKIP {raw['name']}: blank {', '.join(blank)}  {raw.get('notes', '')}".rstrip())
                continue
            gc, gm = _f(raw.get("G_2H_C")), _f(raw.get("G_2H_M"))
            res = analyse(float(raw["E_slab"]), float(raw["G_H"]),
                          None if gc == "" else gc, None if gm == "" else gm,
                          float(raw["G_H2"]))
            rows.append({"name": raw["name"],
                         **{k: raw[k] for k in ("H_site", "H2_site", "Ef", "notes") if k in raw},
                         **res})

    if not rows:
        raise SystemExit(f"No data rows in {a.csv}")
    fields = list(rows[0].keys())
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"Saved {a.out}")
    fmt = lambda v: "n/a" if v is None else f"{v:.3f}"
    for r in rows:
        print(f"{r['name']:<12} dG_H={r['dG_H']:+.3f}  eta_VH={r['eta_VH']:.3f}  "
              f"eta_VT-C={fmt(r['eta_VTC'])}  eta_VT-M={fmt(r['eta_VTM'])}  "
              f"-> {r['pathway']} ({r['limiting_step']}, {r['eta_limiting']:.3f} V)")
    if a.plot:
        plot(rows)


if __name__ == "__main__":
    main()
