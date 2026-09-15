"""Tests for the QE renderer.

The models below are written by hand from the values in
``examples/PSLibrary/Fe_nonmag``; the rendered text must match the reference
files byte for byte.
"""

import itertools

import numpy as np
import pytest

from cif2qewan.exceptions import InputModelError
from cif2qewan.qe.model import (
    AtomicPosition,
    CellParameters,
    KPathPoint,
    KPointsAutomatic,
    KPointsList,
    KPointsPath,
    Namelist,
    NamelistInput,
    PwInput,
    RawValue,
    Species,
)
from cif2qewan.qe.writer import (
    format_value,
    render,
    render_kpoints,
    render_namelist,
    render_namelist_input,
    render_pw_input,
)
from conftest import EXAMPLES

FE = EXAMPLES / "PSLibrary" / "Fe_nonmag"
ALAT = 2.86304
CELL_ALAT = ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5))


def control(calculation, extra=None):
    entries = {
        "calculation": calculation,
        "restart_mode": "from_scratch",
        "prefix": "pwscf",
        "tstress": True,
        "tprnfor": True,
        "pseudo_dir": "/path/to/pslibrary",
        "outdir": "./work",
        "wf_collect": True,
        "disk_io": "low",
    }
    entries.update(extra or {})
    return Namelist("control", entries)


def system(first=None):
    entries = dict(first or {})
    entries.update(
        {
            "ibrav": 0,
            "A": RawValue(f"{ALAT:10.5f}"),
            "nat": 1,
            "ntyp": 1,
            "ecutwfc": 64.0,
            "ecutrho": 782.0,
            "occupations": "smearing",
            "smearing": "m-p",
            "degauss": RawValue("0.010"),
        }
    )
    return Namelist("system", entries)


def electrons(conv_thr):
    return Namelist(
        "electrons",
        {"mixing_mode": "plain", "mixing_beta": 0.1, "conv_thr": RawValue(conv_thr)},
    )


def pw_input(ctrl, sys_, elec, kpoints):
    return PwInput(
        control=ctrl,
        system=sys_,
        electrons=elec,
        species=(Species("Fe", 55.845, "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"),),
        positions=(AtomicPosition("Fe", (0.0, 0.0, 0.0)),),
        kpoints=kpoints,
        cell=CellParameters("alat", CELL_ALAT),
    )


def reference(name):
    return (FE / name).read_text()


# --------------------------------------------------------------------------


def test_format_value():
    assert format_value(True) == ".true." and format_value(False) == ".false."
    assert format_value(31) == "31"
    assert format_value(64.0) == "64.0" and format_value(0.1) == "0.1"
    assert format_value("m-p") == "'m-p'"
    assert format_value(RawValue("1.0d-8")) == "1.0d-8"
    assert format_value(RawValue("")) == ""
    with pytest.raises(InputModelError):
        format_value("it's")
    with pytest.raises(InputModelError):
        format_value([1, 2])


def test_render_namelist_keeps_order_and_indent():
    text = render_namelist(
        Namelist("plot", {"iflag": 3, "fileout": ".xsf"}), indent=" "
    )
    assert text == "&plot\n iflag = 3\n fileout = '.xsf'\n/\n"
    assert render_namelist(Namelist("empty")) == "&empty\n/\n"


def test_render_kpoints_variants():
    assert render_kpoints(KPointsAutomatic((8, 8, 8), (1, 1, 1))) == (
        "K_POINTS {automatic}\n8 8 8  1 1 1\n"
    )
    text = render_kpoints(KPointsList("crystal", ((0, 0, 0, 0.5), (0.5, 0, 0, 0.5))))
    assert text.splitlines()[1] == "2"
    assert (
        text.splitlines()[2]
        == "   0.0000000000    0.0000000000    0.0000000000    0.5000000000"
    )
    path = KPointsPath(
        "crystal_b",
        (KPathPoint((0, 0, 0), 20, "G"), KPathPoint((0.5, -0.5, 0.5), 0, "H")),
    )
    assert render_kpoints(path).splitlines()[2] == (
        "   0.0000000000    0.0000000000    0.0000000000    20    !  G"
    )


def test_scf_input_matches_the_reference():
    pw = pw_input(
        control("scf"), system(), electrons("1.0d-8"), KPointsAutomatic((21, 21, 21))
    )
    assert render_pw_input(pw) == reference("scf.in")
    assert render(pw) == reference("scf.in")


def test_nscf_input_matches_the_reference():
    n = 8
    points = tuple(
        (i / n, j / n, k / n, 1.0 / n**3)
        for i, j, k in itertools.product(range(n), range(n), range(n))
    )
    pw = pw_input(
        control("nscf"),
        system({"nosym": True, "nbnd": 31}),
        electrons("1.d-10"),
        KPointsList("crystal", points),
    )
    assert render_pw_input(pw) == reference("nscf.in")


def test_check_wannier_input_matches_the_reference():
    pw = pw_input(
        control("nscf", {"verbosity": "high"}),
        system({"nbnd": 17}),
        electrons("1.d-8"),
        KPointsAutomatic((8, 8, 8), (1, 1, 1)),
    )
    assert render_pw_input(pw) == reference("check_wannier/nscf.in")


def test_bands_input_matches_the_reference():
    vertices = [
        ((0, 0, 0), 20, "G"),
        ((0.5, -0.5, 0.5), 14, "H"),
        ((0, 0, 0.5), 14, "N"),
        ((0, 0, 0), 17, "G"),
        ((0.25, 0.25, 0.25), 17, "P"),
        ((0.5, -0.5, 0.5), 0, "H"),
        ((0.25, 0.25, 0.25), 10, "P"),
        ((0, 0, 0.5), 10, "N"),
    ]
    path = KPointsPath("crystal_b", tuple(KPathPoint(*v) for v in vertices))
    # band/nscf.in inherits the check_wannier state of 0.2.x: verbosity,
    # nbnd = nexclude + int(1.5 num_wann), conv_thr 1.d-8 and no nosym.
    pw = pw_input(
        control("bands", {"verbosity": "high"}),
        system({"nbnd": 17}),
        electrons("1.d-8"),
        path,
    )
    assert render_pw_input(pw) == reference("band/nscf.in")


def test_namelist_inputs_match_the_references():
    pw2wan = NamelistInput(
        (
            Namelist(
                "inputpp",
                {
                    "outdir": "./work",
                    "prefix": "pwscf",
                    "seedname": "pwscf",
                    "spin_component": "none",
                    "write_mmn": True,
                    "write_amn": True,
                    "write_unk": RawValue(".true."),
                },
            ),
        )
    )
    assert render_namelist_input(pw2wan) == reference("pw2wan.in")

    bands = NamelistInput(
        (
            Namelist(
                "bands",
                {"prefix": "pwscf", "outdir": "./work/", "filband": "bands.out"},
            ),
        )
    )
    assert render(bands) == reference("band/band.in")

    proj = NamelistInput(
        (
            Namelist(
                "projwfc",
                {
                    "prefix": "pwscf",
                    "outdir": "./work",
                    "kresolveddos": False,
                    "degauss": RawValue("0.010"),
                    "Emax": RawValue(""),
                    "Emin": RawValue(""),
                },
            ),
        )
    )
    assert render(proj) == reference("band/proj.in")

    pp = NamelistInput(
        (
            Namelist(
                "inputpp",
                {
                    "prefix": "pwscf",
                    "outdir": "./work",
                    "filplot": "wf_pp",
                    "plot_num": 7,
                    "kpoint": 1,
                    "kband(1)": 5,
                    "kband(2)": 13,
                    "lsign": RawValue(".TRUE."),
                },
            ),
            Namelist("plot", {"iflag": 3, "output_format": 5, "fileout": ".xsf"}),
        )
    )
    assert render(pp) == reference("band/pp.in")


def test_render_rejects_unknown_models():
    with pytest.raises(InputModelError):
        render(np.zeros(3))


def test_cell_rows_use_fifteen_decimals():
    pw = pw_input(
        control("scf"), system(), electrons("1.0d-8"), KPointsAutomatic((21, 21, 21))
    )
    lines = render_pw_input(pw).splitlines()
    i = lines.index("CELL_PARAMETERS {alat}")
    assert (
        lines[i + 1] == " -0.500000000000000   0.500000000000000   0.500000000000000 "
    )
    j = lines.index("ATOMIC_POSITIONS {crystal}")
    assert (
        lines[j + 1]
        == "Fe   0.000000000000000   0.000000000000000   0.000000000000000 "
    )
