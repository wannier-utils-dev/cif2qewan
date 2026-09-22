"""Characterization tests for the physical content of the generated inputs.

``test_examples.py`` pins the exact text. These tests record what the text
means, so that a future implementation can be checked for physical
equivalence even when the formatting changes: which spin settings each flag
produces, how the k meshes relate to each other, that QE and Wannier90 see
the same cell and atoms, and how the band counts follow from the
pseudopotential table.
"""

import re

import numpy as np
import pytest

from cif2qewan.qe.pseudopotential import PseudopotentialTable
from conftest import PACKAGE

# --------------------------------------------------------------------------
# Small parsers for QE and Wannier90 inputs
# --------------------------------------------------------------------------


def namelist(text):
    """``key = value`` assignments of a QE input, keys in lower case."""
    values = {}
    for line in text.splitlines():
        match = re.match(r"\s*([A-Za-z_0-9()]+)\s*=\s*(.*?)\s*$", line)
        if match:
            values[match.group(1).lower()] = match.group(2)
    return values


def card(text, name):
    """The non-empty lines of a QE card, and its option, e.g. ('crystal', [...])."""
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if line.startswith(name):
            option = line[len(name) :].strip(" {}")
            body = []
            for entry in lines[i + 1 :]:
                if not entry.strip() or re.match(r"^[A-Z_]+( |$|\{)", entry):
                    break
                body.append(entry)
            return option, body
    raise AssertionError(f"no {name} card")


def cell_from_qe(text):
    """Lattice vectors in angstrom from CELL_PARAMETERS {alat} and A."""
    values = namelist(text)
    alat = float(values["a"])
    option, rows = card(text, "CELL_PARAMETERS")
    assert option == "alat"
    return np.array([[float(x) for x in row.split()] for row in rows]) * alat


def atoms_from_qe(text):
    option, rows = card(text, "ATOMIC_POSITIONS")
    assert option == "crystal"
    return [(r.split()[0], np.array([float(x) for x in r.split()[1:4]])) for r in rows]


def species_from_qe(text):
    _, rows = card(text, "ATOMIC_SPECIES")
    return {row.split()[0]: row.split()[2] for row in rows}


def win_block(text, name):
    """The lines between ``begin name`` and ``end name`` in a .win file."""
    match = re.search(rf"begin {name}\n(.*?)end {name}", text, re.S)
    assert match, f"no {name} block"
    return [line for line in match.group(1).splitlines() if line.strip()]


def win_value(text, key):
    match = re.search(rf"^{key}\s*[=:]\s*(.*?)\s*$", text, re.M)
    assert match, f"no {key} in .win"
    return match.group(1)


def read(workdir, name):
    return (workdir / name).read_text()


# --------------------------------------------------------------------------
# Spin settings
# --------------------------------------------------------------------------


def test_scf_spin_settings_follow_flags(generated):
    """scf.in: --mag is a collinear spin-polarized run, --so alone is spinor."""
    case, workdir = generated
    scf = namelist(read(workdir, "scf.in"))

    if case.mag:
        assert scf["nspin"] == "2"
        assert scf["starting_magnetization(1)"] == "3.0"
        assert "noncolin" not in scf and "lspinorb" not in scf
    elif case.so:
        assert scf["lspinorb"] == ".true."
        assert scf["noncolin"] == ".true."
        assert "nspin" not in scf
    else:
        for key in ("nspin", "noncolin", "lspinorb", "starting_magnetization(1)"):
            assert key not in scf


@pytest.mark.parametrize("name", ["nscf.in", "check_wannier/nscf.in", "band/nscf.in"])
def test_nscf_spin_settings_follow_flags(generated, name):
    """From nscf.in onwards --mag switches to a noncollinear run with lforcet."""
    case, workdir = generated
    nscf = namelist(read(workdir, name))

    if case.mag:
        assert nscf["noncolin"] == ".true."
        assert nscf["lforcet"] == ".true."
        assert nscf["angle1"] == "0" and nscf["angle2"] == "0"
        assert nscf["lspinorb"] == (".true." if case.so else ".false.")
        assert nscf["starting_magnetization(1)"] == "3.0"
        assert "nspin" not in nscf
    elif case.so:
        assert nscf["lspinorb"] == ".true."
        assert nscf["noncolin"] == ".true."
        assert "lforcet" not in nscf
    else:
        for key in ("nspin", "noncolin", "lspinorb", "lforcet"):
            assert key not in nscf


def test_spinors_in_win_iff_spinor_run(generated):
    case, workdir = generated
    win = read(workdir, "pwscf.win")
    assert ("spinors = .true." in win) == (case.so or case.mag)


def test_pseudopotentials_are_relativistic_only_with_so(generated):
    """--so selects the fully relativistic pseudopotentials from nscf.in on.

    With --mag the collinear scf.in keeps the scalar-relativistic ones.
    """
    case, workdir = generated

    def is_rel(filename):
        return ".rel-" in filename or "_fr" in filename

    scf_rel = [is_rel(f) for f in species_from_qe(read(workdir, "scf.in")).values()]
    nscf_rel = [is_rel(f) for f in species_from_qe(read(workdir, "nscf.in")).values()]
    assert all(r == (case.so and not case.mag) for r in scf_rel)
    assert all(r == case.so for r in nscf_rel)


# --------------------------------------------------------------------------
# Calculation types, convergence and k meshes
# --------------------------------------------------------------------------


def test_calculation_types_and_convergence(generated):
    case, workdir = generated
    scf = namelist(read(workdir, "scf.in"))
    nscf = namelist(read(workdir, "nscf.in"))
    check = namelist(read(workdir, "check_wannier/nscf.in"))
    band = namelist(read(workdir, "band/nscf.in"))

    assert scf["calculation"] == "'scf'"
    assert scf["conv_thr"] == "1.0d-8"
    assert nscf["calculation"] == "'nscf'"
    assert nscf["conv_thr"] == "1.d-10"
    assert nscf["nosym"] == ".true."
    assert check["calculation"] == "'nscf'"
    assert check["conv_thr"] == "1.d-8"
    assert check["verbosity"] == "'high'"
    assert "nosym" not in check
    assert band["calculation"] == "'bands'"
    for values in (scf, nscf, check, band):
        assert values["occupations"] == "'smearing'"
        assert values["smearing"] == "'m-p'"
        assert values["degauss"] == "0.010"


def test_k_meshes_are_consistent(generated):
    """nscf mesh = scf mesh clamped to [4, 8]; .win and check_wannier reuse it."""
    case, workdir = generated

    option, rows = card(read(workdir, "cif_scf.in"), "K_POINTS")
    assert option == "automatic"
    scf_mesh = [int(x) for x in rows[0].split()[:3]]
    option, rows = card(read(workdir, "scf.in"), "K_POINTS")
    assert [int(x) for x in rows[0].split()[:3]] == scf_mesh

    expected = [min(max(n, 4), 8) for n in scf_mesh]

    option, rows = card(read(workdir, "nscf.in"), "K_POINTS")
    assert option == "crystal"
    assert int(rows[0]) == np.prod(expected)
    assert len(rows) == np.prod(expected) + 1
    weights = [float(r.split()[3]) for r in rows[1:]]
    assert sum(weights) == pytest.approx(1.0)

    win = read(workdir, "pwscf.win")
    assert [int(x) for x in win_value(win, "mp_grid").split()] == expected
    assert len(win_block(win, "kpoints")) == np.prod(expected)
    for qe_row, win_row in zip(rows[1:], win_block(win, "kpoints")):
        qe_k = [float(x) for x in qe_row.split()[:3]]
        win_k = [float(x) for x in win_row.split()[:3]]
        assert qe_k == pytest.approx(win_k)

    option, rows = card(read(workdir, "check_wannier/nscf.in"), "K_POINTS")
    assert option == "automatic"
    assert rows[0].split() == [str(n) for n in expected] + ["1", "1", "1"]


def test_band_path_is_shared_by_qe_and_win(generated):
    """band/nscf.in and the .win kpoint_path describe the same k path.

    In the QE crystal_b list a point with zero divisions marks a jump in the
    path; the .win file lists one segment per consecutive pair of points
    that is not separated by such a jump.
    """
    case, workdir = generated

    option, rows = card(read(workdir, "band/nscf.in"), "K_POINTS")
    assert option == "crystal_b"
    assert int(rows[0]) == len(rows) - 1
    qe_points = [[float(x) for x in r.split()[:3]] for r in rows[1:]]
    qe_ndiv = [int(r.split()[3]) for r in rows[1:]]
    qe_labels = [r.split("!")[1].strip() for r in rows[1:]]
    segments = [i for i in range(len(qe_points) - 1) if qe_ndiv[i] > 0]

    path = win_block(read(workdir, "pwscf.win"), "kpoint_path")
    assert len(path) == len(segments)
    for i, segment in zip(segments, path):
        fields = segment.split()
        assert fields[0] == qe_labels[i] and fields[4] == qe_labels[i + 1]
        assert [float(x) for x in fields[1:4]] == pytest.approx(qe_points[i])
        assert [float(x) for x in fields[5:8]] == pytest.approx(qe_points[i + 1])


# --------------------------------------------------------------------------
# Structure: QE and Wannier90 must describe the same crystal
# --------------------------------------------------------------------------


def test_qe_and_win_describe_the_same_cell(generated):
    case, workdir = generated
    qe_cell = cell_from_qe(read(workdir, "scf.in"))
    win = read(workdir, "pwscf.win")
    block = win_block(win, "unit_cell_cart")
    assert block[0].strip() == "ang"
    win_cell = np.array([[float(x) for x in row.split()] for row in block[1:]])

    assert np.linalg.det(qe_cell) > 0
    assert win_cell == pytest.approx(qe_cell, abs=1e-6)
    assert cell_from_qe(read(workdir, "cif_scf.in")) == pytest.approx(qe_cell)


def test_qe_and_win_describe_the_same_atoms(generated):
    case, workdir = generated
    qe_atoms = atoms_from_qe(read(workdir, "scf.in"))
    assert len(qe_atoms) == int(namelist(read(workdir, "scf.in"))["nat"])

    rows = win_block(read(workdir, "pwscf.win"), "atoms_frac")
    win_atoms = [
        (r.split()[0], np.array([float(x) for x in r.split()[1:4]])) for r in rows
    ]
    assert len(win_atoms) == len(qe_atoms)
    for (qe_el, qe_pos), (win_el, win_pos) in zip(qe_atoms, win_atoms):
        assert qe_el == win_el
        assert (qe_pos - win_pos) % 1.0 == pytest.approx(0.0, abs=1e-8)


# --------------------------------------------------------------------------
# Wannier function and band counts
# --------------------------------------------------------------------------

ORBITAL_SIZE = {"s": 1, "p": 3, "d": 5, "f": 7}


def test_wannier_counts_follow_the_pseudopotential_table(generated):
    case, workdir = generated
    table = PseudopotentialTable.from_csv(PACKAGE / case.pp_csv)
    win = read(workdir, "pwscf.win")
    species = species_from_qe(read(workdir, "scf.in"))
    atoms = atoms_from_qe(read(workdir, "scf.in"))

    # projections: one "El:s,p,d" line per species, orbitals from the table
    projections = {}
    for line in win_block(win, "projections"):
        element, orbitals = line.split(":")
        projections[element] = orbitals.split(",")
    assert set(projections) == set(species)
    for element, orbitals in projections.items():
        assert "".join(orbitals) == table.lookup(element).orbitals

    num_wann = sum(ORBITAL_SIZE[o] for el, _ in atoms for o in projections[el])
    nexclude = sum(table.lookup(el).nexclude for el, _ in atoms)
    factor = case.spin_factor

    assert int(win_value(win, "num_wann")) == num_wann * factor
    assert int(win_value(win, "num_bands")) == 3 * num_wann * factor
    assert win_value(win, "exclude_bands") == f"1-{nexclude * factor}"

    nbnd = int(namelist(read(workdir, "nscf.in"))["nbnd"])
    assert nbnd == nexclude * factor + 3 * num_wann * factor

    pp = namelist(read(workdir, "band/pp.in"))
    assert int(pp["kband(1)"]) == nexclude * factor + 1
    assert int(pp["kband(2)"]) == nexclude * factor + num_wann * factor


def test_pw2wan_and_win_share_the_seedname(generated):
    case, workdir = generated
    pw2wan = namelist(read(workdir, "pw2wan.in"))
    scf = namelist(read(workdir, "scf.in"))
    assert pw2wan["prefix"] == scf["prefix"] == "'pwscf'"
    assert pw2wan["seedname"] == "'pwscf'"
    assert pw2wan["outdir"] == scf["outdir"]
    assert pw2wan["write_unk"] == ".true."
    assert pw2wan["spin_component"] == "'none'"
