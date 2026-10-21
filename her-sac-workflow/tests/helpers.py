"""Synthetic fixtures for the tests (no real VASP output is needed)."""
import copy
import os
import sys
import tempfile

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import herlib  # noqa: E402

TEMPLATE = os.path.join(ROOT, "assets", "POSCAR_template")
CONVERGED_LOG = "  5 F= -.2E+03\n reached required accuracy - stopping structural energy minimisation\n"
RELAX_INCAR = "".join(f"{k} = {v}\n" for k, v in herlib.incar_template("relaxation").items())
FREQ_INCAR = "".join(f"{k} = {v}\n" for k, v in herlib.incar_template("frequency").items())
OUTCAR_DONE = " General timing and accounting informations for this job:\n"
VASPKIT_501 = ("Zero-point energy E_ZPE   :      13.814 kcal/mol    {zpe:.6f} eV\n"
               "Entropy contribution T*S  :   15147.343 J/(mol)     {ts:.6f} eV\n")


def oszicar(e0, mag=1.0, step=5):
    return (f"DAV:   1    -0.1E+03   -0.1E+03   -0.1E+00  1000   0.1E+00\n"
            f"  {step} F= {e0 - 0.001:.8E} E0= {e0:.8E}  d E =-.1E-05  mag=     {mag:.4f}\n")


def structure(h_offsets=(), metal_shift=(0.0, 0.0, 0.0), ligand_moves=None, flat=True, h_on=()):
    """Template structure (C N S Ir) with H atoms placed at Cartesian offsets from the metal.

    h_on: extra H atoms as (atom index, Cartesian offset) pairs, e.g. on a ligand C.

    flat=True first moves the metal into the plane of its ligands; in the relaxed template the Ir sits
    0.61 Å out of that plane.
    """
    s = copy.deepcopy(herlib.read_structure(TEMPLATE))
    if flat:
        from check_structure import check_geometry
        dz = check_geometry(s, {})[0]["out_of_plane"]
        metal_shift = (metal_shift[0], metal_shift[1], metal_shift[2] - dz)
    m, _ = herlib.find_site(s)
    inv = herlib.inverse(s["lattice"])
    cart = [s["cart"][m][k] + metal_shift[k] for k in range(3)]
    s["frac"][m] = herlib.mat_vec_t(inv, cart)
    for idx, vec in (ligand_moves or {}).items():
        c = [s["cart"][idx][k] + vec[k] for k in range(3)]
        s["frac"][idx] = herlib.mat_vec_t(inv, c)
    s["cart"] = [herlib.to_cart(s["lattice"], f) for f in s["frac"]]
    for at, off in [(m, o) for o in h_offsets] + list(h_on):
        c = [s["cart"][at][k] + off[k] for k in range(3)]
        s["elements"].append("H")
        s["frac"].append(herlib.mat_vec_t(inv, c))
        s["flags"].append([])
    s["cart"] = [herlib.to_cart(s["lattice"], f) for f in s["frac"]]
    return s


def poscar_text(s):
    d = tempfile.mkdtemp()
    p = os.path.join(d, "POSCAR")
    herlib.write_structure(s, p)
    with open(p) as fh:
        return fh.read()


def make_dir(files, parent=None, name=None):
    d = tempfile.mkdtemp() if parent is None else os.path.join(parent, name)
    os.makedirs(d, exist_ok=True)
    for fname, text in files.items():
        with open(os.path.join(d, fname), "w") as fh:
            fh.write(text)
    return d


def relax_job(s, e0=-200.0, mag=1.0, parent=None, name=None, extra=None):
    files = {"INCAR": RELAX_INCAR, "log": CONVERGED_LOG, "OSZICAR": oszicar(e0, mag),
             "CONTCAR": poscar_text(s), "POSCAR": poscar_text(s)}
    files.update(extra or {})
    return make_dir(files, parent, name)


def freq_job(s, zpe, ts, parent=None, name=None, imaginary=()):
    outcar = " Eigenvectors and eigenvalues of the dynamical matrix\n"
    outcar += "   1 f  =   99.000000 THz   622.000000 2PiTHz 3302.000000 cm-1   409.000000 meV\n"
    for k, cm in enumerate(imaginary, 2):
        outcar += f"   {k} f/i=    1.000000 THz     6.000000 2PiTHz {cm:10.6f} cm-1     4.000000 meV\n"
    outcar += OUTCAR_DONE
    return make_dir({"INCAR": FREQ_INCAR, "OUTCAR": outcar, "OSZICAR": oszicar(-1.0),
                     "CONTCAR": poscar_text(s), "vaspkit_501.txt": VASPKIT_501.format(zpe=zpe, ts=ts)},
                    parent, name)
