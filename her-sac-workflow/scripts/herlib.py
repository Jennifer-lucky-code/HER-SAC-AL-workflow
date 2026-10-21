"""Shared helpers for the HER-SAC workflow scripts (standard library only; not a CLI).

- POSCAR/CONTCAR reading with minimum-image geometry, and finding the metal site by geometry
  (metal = the single transition-metal atom, ligands = its 4 nearest non-H neighbours). Atom indices in
  assets/sites.json refer to the template only: build_poscar.py regroups atoms by element.
- INCAR, OSZICAR, OUTCAR and vaspkit 501 parsing.
- The INCAR templates and reference energies, read from references/vasp_settings.md, which stays the one place
  where they are written down.
- loop_state.json and workflow_log.jsonl (schemas in references/state_files.md).
"""
import datetime
import itertools
import json
import math
import os
import re
import shlex
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
SETTINGS_MD = os.path.join(HERE, "..", "references", "vasp_settings.md")
TEMPERATURE = 298.15  # K, ZPE/TS temperature used with vaspkit 501 (references/vasp_settings.md)

TRANSITION_METALS = set(("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
                         "La Hf Ta W Re Os Ir Pt Au Hg").split())

# Covalent radii in Å (Cordero et al., Dalton Trans. 2008, 2832; C sp2; Mn/Fe/Co high-spin values).
# Used only for geometric bond detection in check_structure.py, never in an energy.
COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.73, "N": 0.71, "O": 0.66, "F": 0.57, "Si": 1.11, "P": 1.07, "S": 1.05,
    "Cl": 1.02, "Se": 1.20,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.61, "Fe": 1.52, "Co": 1.50, "Ni": 1.24,
    "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46,
    "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "La": 2.07, "Hf": 1.75, "Ta": 1.70, "W": 1.62,
    "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32,
}


# --------------------------------------------------------------------------- structures

def read_structure(path):
    """POSCAR/CONTCAR (VASP 5) -> dict(lattice, elements, frac, cart, selective, flags, coord_type, comment)."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    scale = float(lines[1].split()[0])
    lat = [[scale * float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    species, counts = lines[5].split(), [int(x) for x in lines[6].split()]
    elems = [s for s, n in zip(species, counts) for _ in range(n)]
    i = 7
    selective = lines[i].strip().lower().startswith("s")
    if selective:
        i += 1
    coord_type = lines[i].strip()
    cartesian = coord_type[:1].lower() in ("c", "k")
    i += 1
    frac, flags = [], []
    inv = inverse(lat)
    for k in range(len(elems)):
        parts = lines[i + k].split()
        xyz = [float(x) for x in parts[:3]]
        if cartesian:
            xyz = mat_vec_t(inv, [scale * x for x in xyz])
        frac.append(xyz)
        flags.append(parts[3:6] if selective else [])
    return {"comment": lines[0], "lattice": lat, "elements": elems, "frac": frac,
            "cart": [to_cart(lat, f) for f in frac], "selective": selective, "flags": flags,
            "cartesian": cartesian}


def write_structure(s, path, comment=None):
    """Write a structure from read_structure() (same element order, same coordinate type)."""
    order = []
    for e in s["elements"]:
        if not order or order[-1] != e:
            order.append(e)
    counts, prev = [], None
    for e in s["elements"]:
        if e != prev:
            counts.append(0)
            prev = e
        counts[-1] += 1
    out = [comment or s["comment"], "1.0"]
    out += ["  ".join(f"{x:.10f}" for x in row) for row in s["lattice"]]
    out += ["  ".join(order), "  ".join(map(str, counts))]
    if s["selective"]:
        out.append("Selective dynamics")
    out.append("Cartesian" if s["cartesian"] else "Direct")
    for f, fl in zip(s["frac"], s["flags"]):
        xyz = to_cart(s["lattice"], f) if s["cartesian"] else f
        out.append("  ".join(f"{x:.10f}" for x in xyz) + ("  " + "  ".join(fl) if fl else ""))
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def to_cart(lat, f):
    return [sum(f[k] * lat[k][j] for k in range(3)) for j in range(3)]


def mat_vec_t(inv, v):
    """Fractional coordinates of Cartesian v, given inv = inverse(lattice) (rows = lattice vectors)."""
    return [sum(v[j] * inv[j][k] for j in range(3)) for k in range(3)]


def inverse(m):
    a, b, c = m[0], m[1], m[2]
    det = (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
           + a[2] * (b[0] * c[1] - b[1] * c[0]))
    if abs(det) < 1e-12:
        raise ValueError("singular lattice")
    cof = [[b[1] * c[2] - b[2] * c[1], -(a[1] * c[2] - a[2] * c[1]), a[1] * b[2] - a[2] * b[1]],
           [-(b[0] * c[2] - b[2] * c[0]), a[0] * c[2] - a[2] * c[0], -(a[0] * b[2] - a[2] * b[0])],
           [b[0] * c[1] - b[1] * c[0], -(a[0] * c[1] - a[1] * c[0]), a[0] * b[1] - a[1] * b[0]]]
    return [[cof[i][j] / det for j in range(3)] for i in range(3)]


def mic_vector(s, i, j):
    """Shortest Cartesian vector from atom i to atom j over periodic images (minimum image)."""
    d = [s["frac"][j][k] - s["frac"][i][k] for k in range(3)]
    d = [x - round(x) for x in d]
    best = None
    for shift in itertools.product((-1, 0, 1), repeat=3):
        v = to_cart(s["lattice"], [d[k] + shift[k] for k in range(3)])
        n = math.sqrt(sum(x * x for x in v))
        if best is None or n < best[0]:
            best = (n, v)
    return best[1]


def distance(s, i, j):
    return math.sqrt(sum(x * x for x in mic_vector(s, i, j)))


def find_site(s):
    """(metal index, [4 ligand indices nearest first]) by geometry; raises ValueError if not exactly one TM."""
    tms = [i for i, e in enumerate(s["elements"]) if e in TRANSITION_METALS]
    if len(tms) != 1:
        raise ValueError(f"expected exactly one transition-metal atom, found {len(tms)}")
    m = tms[0]
    others = sorted((distance(s, m, j), j) for j, e in enumerate(s["elements"]) if j != m and e != "H")
    return m, [j for _, j in others[:4]]


def neighbours(s, i, exclude_h=False):
    """[(distance, index)] of all other atoms, nearest first."""
    return sorted((distance(s, i, j), j) for j, e in enumerate(s["elements"])
                  if j != i and not (exclude_h and e == "H"))


# --------------------------------------------------------------------------- VASP files

def parse_incar(path):
    """INCAR -> {KEY: value string} (comments after # or ! removed; ';' separates statements)."""
    vals = {}
    if not os.path.isfile(path):
        return vals
    with open(path) as fh:
        text = fh.read()
    for line in text.splitlines():
        line = re.split(r"[#!]", line, maxsplit=1)[0]
        for stmt in line.split(";"):
            if "=" in stmt:
                k, v = stmt.split("=", 1)
                k = k.strip().upper()
                if re.fullmatch(r"[A-Z_][A-Z0-9_]*", k):
                    vals[k] = " ".join(v.split())
    return vals


def norm_value(v):
    """Normalise an INCAR value for comparison: 1E-05 == 0.00001, .TRUE. == T, case-insensitive."""
    out = []
    for tok in str(v).split():
        t = tok.upper().strip()
        if t in (".TRUE.", "T", "TRUE"):
            out.append("T")
            continue
        if t in (".FALSE.", "F", "FALSE"):
            out.append("F")
            continue
        try:
            out.append(repr(float(t.replace("D", "E"))))
        except ValueError:
            out.append(t)
    return " ".join(out)


def set_incar(path, changes):
    """Apply {KEY: value or None (= remove)} to an INCAR in place, keeping every other line as it was."""
    lines = []
    if os.path.isfile(path):
        with open(path) as fh:
            lines = fh.read().splitlines()
    done = set()
    out = []
    for line in lines:
        m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
        key = m.group(1).upper() if m else None
        if key in changes and ";" not in line:
            done.add(key)
            if changes[key] is not None:
                out.append(f"{key} = {changes[key]}")
            continue
        out.append(line)
    for k, v in changes.items():
        if k not in done and v is not None:
            out.append(f"{k} = {v}")
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def is_frequency(incar):
    try:
        return 5 <= int(incar.get("IBRION", 0)) <= 8
    except ValueError:
        return False


def oszicar_last(path):
    """Last ionic-step line of OSZICAR -> {'step', 'F', 'E0', 'mag'} (mag None if not printed); None if absent."""
    if not os.path.isfile(path):
        return None
    last = None
    with open(path, errors="replace") as fh:
        text = fh.read()
    for line in text.splitlines():
        if re.match(r"\s*\d+\s+F=", line):
            last = line
    if last is None:
        return None
    num = r"([-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?)"
    get = lambda key: re.search(key + r"\s*=\s*" + num, last)
    m_mag = get("mag")
    return {"step": int(last.split()[0]), "F": float(get("F").group(1)), "E0": float(get("E0").group(1)),
            "mag": float(m_mag.group(1)) if m_mag else None}


def outcar_imaginary(path):
    """Imaginary modes from OUTCAR: [{'mode', 'THz', 'cm-1', 'meV'}] (last frequency block only)."""
    if not os.path.isfile(path):
        return None
    blocks, cur = [], None
    rx = re.compile(r"^\s*(\d+)\s+f\s*(/i)?\s*=\s*([\d.]+)\s+THz\s+[\d.]+\s+2PiTHz\s+([\d.]+)\s+cm-1\s+([\d.]+)\s+meV")
    with open(path, errors="replace") as fh:
        text = fh.read()
    for line in text.splitlines():
        if ("Eigenvectors and eigenvalues of the dynamical matrix" in line
                or "Eigenvectors after division by SQRT(mass)" in line):
            cur = []
            blocks.append(cur)
            continue
        m = rx.match(line)
        if m:
            if cur is None:
                cur = []
                blocks.append(cur)
            if m.group(2):
                cur.append({"mode": int(m.group(1)), "THz": float(m.group(3)), "cm-1": float(m.group(4)),
                            "meV": float(m.group(5))})
    return blocks[-1] if blocks else []


def outcar_nelect(path):
    if not os.path.isfile(path):
        return None
    with open(path, errors="replace") as fh:
        for line in fh:
            m = re.search(r"NELECT\s*=\s*([\d.]+)", line)
            if m:
                return float(m.group(1))
    return None


def parse_vaspkit_501(text):
    """vaspkit 501 output -> {'ZPE': eV, 'TS': eV} (second-to-last column of the two lines, see vasp_settings.md).

    Missing lines give None for that value.
    """
    out = {"ZPE": None, "TS": None}
    for line in text.splitlines():
        for key, label in (("ZPE", "Zero-point energy E_ZPE"), ("TS", "Entropy contribution T*S")):
            if label in line:
                cols = line.split()
                out[key] = float(cols[-2])
    return out


# --------------------------------------------------------------------------- references/vasp_settings.md

def _settings_text(path=None):
    with open(path or SETTINGS_MD, encoding="utf-8") as fh:
        return fh.read()


def reference_energy(name, path=None):
    """Value of a row in the 'Reference energies' table of vasp_settings.md (e.g. 'G_H₂'); None if TODO."""
    text = _settings_text(path)
    sec = text.split("## Reference energies", 1)[1].split("\n## ", 1)[0]
    for line in sec.splitlines():
        cols = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cols) >= 2 and cols[0] == name:
            try:
                return float(cols[1])
            except ValueError:
                return None
    raise KeyError(f"{name} not in the reference-energy table of vasp_settings.md")


def g_h2(path=None):
    return reference_energy("G_H₂", path)


def incar_template(kind, path=None):
    """INCAR template ('relaxation' or 'frequency') from vasp_settings.md as {KEY: value}."""
    heading = {"relaxation": "## Relaxation INCAR", "frequency": "## Frequency INCAR"}[kind]
    sec = _settings_text(path).split(heading, 1)[1]
    block = sec.split("```", 2)[1]
    vals = {}
    for line in block.splitlines():
        line = re.split(r"[#!]", line, maxsplit=1)[0]
        if "=" in line:
            k, v = line.split("=", 1)
            vals[k.strip().upper()] = " ".join(v.split())
    return vals


# Keys the relaxation template allows to change per job (vasp_settings.md: SYSTEM, MAGMOM, NPAR).
INCAR_FREE_KEYS = {"SYSTEM", "MAGMOM", "NPAR"}


def incar_diff(job_incar, template):
    """['KEY=job (template T)'] for every key that differs from the template, ignoring INCAR_FREE_KEYS."""
    out = []
    for k in sorted(set(job_incar) | set(template)):
        if k in INCAR_FREE_KEYS:
            continue
        a, b = job_incar.get(k), template.get(k)
        if a is None or b is None or norm_value(a) != norm_value(b):
            out.append(f"{k}={'<unset>' if a is None else a} (template {'<unset>' if b is None else b})")
    return out


# --------------------------------------------------------------------------- state and log

STATE_FILE = "loop_state.json"
LOG_FILE = "workflow_log.jsonl"
EVENTS = ("submitted", "checked", "fixed", "converged", "flagged", "added")


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def load_state(path=STATE_FILE):
    if os.path.isfile(path):
        with open(path) as fh:
            state = json.load(fh)
    else:
        state = {}
    state.setdefault("jobs", {})
    return state


def save_state(state, path=STATE_FILE):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(state, fh, indent=1, ensure_ascii=False)
    os.replace(tmp, path)


def job_key(job_dir, state_path=STATE_FILE):
    """Job path relative to the folder that holds loop_state.json (the project root)."""
    root = os.path.dirname(os.path.abspath(state_path))
    return os.path.relpath(os.path.abspath(job_dir), root).replace(os.sep, "/")


def job_entry(state, key):
    return state["jobs"].setdefault(key, {"jobid": "", "status": "queued", "retries": 0, "history": []})


def log_event(path, event, job_dir="", error_type="", change="", jobid="", **extra):
    if event not in EVENTS:
        raise ValueError(f"unknown event {event!r}")
    rec = {"time": now(), "job_dir": job_dir, "event": event, "error_type": error_type,
           "change": change, "jobid": jobid}
    rec.update(extra)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
    return rec


def read_log(path=LOG_FILE):
    out = []
    if not os.path.isfile(path):
        return out
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    for n, line in enumerate(lines, 1):
        line = line.strip()
        if line:
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                raise ValueError(f"{path}:{n}: not valid JSON")
    return out


# --------------------------------------------------------------------------- submission

def submit(job_dir, scheduler="auto", script="sub.sh", cmd=None):
    """Submit job_dir/script, write the scheduler ID to job_dir/JOBID and return it.

    scheduler: 'pbs' (qsub), 'slurm' (sbatch --parsable) or 'auto' (whichever is on PATH).
    cmd overrides the command (split with shlex; the script name is not appended), e.g. for tests.
    """
    if cmd:
        argv = shlex.split(cmd)
    else:
        if scheduler == "auto":
            scheduler = "pbs" if shutil.which("qsub") else "slurm" if shutil.which("sbatch") else ""
        if scheduler == "pbs":
            argv = ["qsub", script]
        elif scheduler == "slurm":
            argv = ["sbatch", "--parsable", script]
        else:
            raise RuntimeError("no scheduler found (qsub/sbatch); use --scheduler or run on the login node")
    r = subprocess.run(argv, cwd=job_dir, capture_output=True, text=True)
    jobid = r.stdout.strip().splitlines()[-1].strip() if r.stdout.strip() else ""
    if r.returncode != 0 or not jobid:
        raise RuntimeError(f"submission failed in {job_dir}: {r.stderr.strip() or r.stdout.strip()}")
    with open(os.path.join(job_dir, "JOBID"), "w") as fh:
        fh.write(jobid + "\n")
    return jobid
