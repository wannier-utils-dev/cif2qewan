"""Tests for the Wannier90 renderer (DEVELOPMENT_PLAN.md Step 5).

The rendered ``.win`` for the Fe_nonmag and Fe_so examples must contain
the same blocks (projections, cell, atoms, k mesh, k path) and the same
keyword values as the reference files. The layout of the keyword lines is
normalized (see ``cif2qewan.wannier90.writer``), so the files are compared
block by block rather than byte by byte.
"""

import itertools
import re

import numpy as np
import pytest

from cif2qewan.exceptions import InputModelError
from cif2qewan.wannier90.model import AtomFrac, KPathSegment, Projection, Wannier90Input
from cif2qewan.wannier90.writer import HEADER, format_parameter, render_win
from conftest import EXAMPLES

ALAT = 2.86304
BCC = np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]) * ALAT
PATH = [
    ("G", (0, 0, 0), "H", (0.5, -0.5, 0.5)),
    ("H", (0.5, -0.5, 0.5), "N", (0, 0, 0.5)),
    ("N", (0, 0, 0.5), "G", (0, 0, 0)),
    ("G", (0, 0, 0), "P", (0.25, 0.25, 0.25)),
    ("P", (0.25, 0.25, 0.25), "H", (0.5, -0.5, 0.5)),
    ("P", (0.25, 0.25, 0.25), "N", (0, 0, 0.5)),
]
PARAMETERS = {
    "dis_num_iter": 200,
    "num_iter": 0,
    "dis_froz_max": -200,
    "dis_froz_min": -200,
    "bands_plot": True,
    "write_hr": True,
    "write_tb": True,
    "fermi_surface_plot": True,
    "wannier_plot": True,
}


def fe_win(spinors):
    n = 8
    factor = 2 if spinors else 1
    return Wannier90Input(
        num_wann=9 * factor,
        num_bands=27 * factor,
        unit_cell_cart=BCC,
        atoms_frac=(AtomFrac("Fe", (0.0, 0.0, 0.0)),),
        mp_grid=(n, n, n),
        kpoints=tuple(
            (i / n, j / n, k / n) for i, j, k in itertools.product(range(n), repeat=3)
        ),
        projections=(Projection("Fe", ("s", "p", "d")),),
        spinors=spinors,
        exclude_bands=(1, 4 * factor),
        kpoint_path=tuple(KPathSegment(*segment) for segment in PATH),
        parameters=dict(PARAMETERS),
    )


def blocks(text):
    """name -> list of lines for every begin/end block."""
    return {
        name: body.strip("\n").splitlines()
        for name, body in re.findall(r"begin (\w+)\n(.*?)end \1", text, re.S)
    }


def keywords(text):
    """key -> value for every ``key = value`` / ``key: value`` line outside blocks."""
    outside = re.sub(r"begin (\w+)\n.*?end \1\n", "", text, flags=re.S)
    return dict(re.findall(r"^(\w+)\s*[=:]\s*(.*?)\s*$", outside, re.M))


@pytest.mark.parametrize("example,spinors", [("Fe_nonmag", False), ("Fe_so", True)])
def test_rendered_win_matches_the_reference_content(example, spinors):
    reference = (EXAMPLES / "PSLibrary" / example / "pwscf.win").read_text()
    rendered = render_win(fe_win(spinors))

    assert rendered.startswith(HEADER + "\n")
    assert blocks(rendered) == blocks(reference)
    assert keywords(rendered) == keywords(reference)
    assert ("spinors = .true." in rendered) == spinors
    assert rendered.endswith("end kpoint_path\n")


def test_layout_is_fixed():
    text = render_win(fe_win(True))
    head = text.splitlines()[:6]
    assert head == [
        HEADER,
        "num_bands = 54",
        "num_wann = 18",
        "exclude_bands = 1-8",
        "",
        "spinors = .true.",
    ]
    order = [
        text.index("dis_num_iter"),
        text.index("begin projections"),
        text.index("begin unit_cell_cart"),
        text.index("begin atoms_frac"),
        text.index("mp_grid:"),
        text.index("begin kpoints"),
        text.index("begin kpoint_path"),
    ]
    assert order == sorted(order)
    assert "\n\n\n" not in text  # never more than one blank line


def test_optional_parts_are_omitted():
    win = fe_win(False)
    win.exclude_bands = None
    win.kpoint_path = ()
    win.parameters = {}
    text = render_win(win)
    assert "exclude_bands" not in text and "spinors" not in text
    assert "kpoint_path" not in text
    assert text.endswith("end kpoints\n")
    assert text.splitlines()[3] == ""  # blank line after the counts, then projections
    assert text.splitlines()[4] == "begin projections"


def test_numeric_formats_match_the_reference():
    lines = render_win(fe_win(False)).splitlines()
    i = lines.index("begin unit_cell_cart")
    assert lines[i + 1] == "ang"
    assert lines[i + 2] == "  -1.4315200    1.4315200    1.4315200"
    j = lines.index("begin atoms_frac")
    assert (
        lines[j + 1]
        == "Fe   0.000000000000000   0.000000000000000   0.000000000000000 "
    )
    k = lines.index("begin kpoints")
    assert lines[k + 2] == "   0.0000000000    0.0000000000    0.1250000000"
    p = lines.index("begin kpoint_path")
    assert lines[p + 1] == (
        "G   0.0000000000   0.0000000000   0.0000000000  "
        "H   0.5000000000  -0.5000000000   0.5000000000"
    )


def test_format_parameter():
    assert format_parameter(True) == ".true."
    assert format_parameter(False) == ".false."
    assert format_parameter(200) == "200"
    assert format_parameter(1.5) == "1.5"
    assert format_parameter("xsf") == "xsf"
    with pytest.raises(InputModelError):
        format_parameter([1])
