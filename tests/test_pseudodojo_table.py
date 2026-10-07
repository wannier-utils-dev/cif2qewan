"""Tests of ``tools/pseudodojo_table.py`` on a tiny fake PseudoDojo set.

The real sets are not shipped; the rules (cutoff from the "high" hint,
nexclude carried over from a reference table and corrected for semicore
shells) are checked on hand-made files.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

from conftest import REPO_ROOT

SCRIPT = REPO_ROOT / "tools" / "pseudodojo_table.py"


@pytest.fixture(scope="module")
def tool():
    spec = importlib.util.spec_from_file_location("pseudodojo_table", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def write_upf(path, element, z, core, valence):
    lines = [f"{n} {l} {f}" for n, l, f in core + valence]
    path.write_text(
        '<UPF version="2.0.1">\n<PP_INFO>\n<PP_INPUTFILE>\n'
        "# ATOM AND REFERENCE CONFIGURATION\n"
        "# atsym  z   nc   nv     iexc    psfile\n"
        f"{element} {z:.2f}    {len(core)}    {len(valence)}       4      upf\n"
        "#\n#   n    l    f        energy (Ha)\n" + "\n".join(lines) + "\n"
        "#\n# PSEUDOPOTENTIAL AND OPTIMIZATION\n# lmax\n2\n</PP_INPUTFILE>\n"
    )


def write_djrepo(path, high=None):
    data = (
        {
            "hints": {
                "low": {"ecut": high - 8},
                "normal": {"ecut": high - 4},
                "high": {"ecut": high},
            }
        }
        if high
        else {}
    )
    path.write_text(json.dumps(data))


def make_set(root, name, fe_valence, hints):
    upf = root / f"{name}_upf"
    djrepo = root / f"{name}_djrepo"
    upf.mkdir()
    djrepo.mkdir()
    write_upf(upf / "Fe.upf", "Fe", 26, [(1, 0, 2), (2, 0, 2), (2, 1, 6)], fe_valence)
    write_upf(upf / "I.upf", "I", 53, [(1, 0, 2)], [(4, 2, 10), (5, 0, 2), (5, 1, 5)])
    write_djrepo(djrepo / "Fe.djrepo", hints.get("Fe"))
    write_djrepo(djrepo / "I.djrepo", hints.get("I"))
    return upf, djrepo


def test_valence_configuration_and_rules(tmp_path, tool):
    fe_old = [(3, 0, 2), (3, 1, 6), (3, 2, 6), (4, 0, 2)]
    ref_upf, ref_djrepo = make_set(tmp_path, "ref", fe_old, {"Fe": 53.0, "I": 41.0})
    assert tool.valence_configuration(ref_upf / "Fe.upf") == (
        (3, "s"),
        (3, "p"),
        (3, "d"),
        (4, "s"),
    )
    reference = tmp_path / "ref.csv"
    reference.write_text(
        "atom,pp_file_name,nexclude,orbitals,ecutwfc,ecutrho\n"
        "Fe,Fe_sr,4,dsp,106.0,424.0\nI,I_sr,6,p,82.0,328.0\n"
    )
    # the new set drops the 3s shell of Fe and has no hints for I
    new_upf, new_djrepo = make_set(tmp_path, "new", fe_old[1:], {"Fe": 50.0})
    output = tmp_path / "new.csv"
    status = tool.main(
        [
            "--upf",
            str(new_upf),
            "--djrepo",
            str(new_djrepo),
            "--fallback-djrepo",
            str(ref_djrepo),
            "--reference",
            str(reference),
            "--reference-upf",
            str(ref_upf),
            "--output",
            str(output),
        ]
    )
    assert status == 0
    table = tool.read_table(output)
    assert len(table) == len(tool.ELEMENTS)
    assert table["Fe"] == {
        "atom": "Fe",
        "pp_file_name": "Fe_sr",
        "nexclude": "3",
        "orbitals": "dsp",
        "ecutwfc": "100.0",
        "ecutrho": "400.0",
    }
    assert table["I"]["ecutwfc"] == "82.0" and table["I"]["nexclude"] == "6"
    assert table["Cu"] == {
        "atom": "Cu",
        "pp_file_name": "Cu_sr",
        "nexclude": "",
        "orbitals": "",
        "ecutwfc": "",
        "ecutrho": "",
    }
    # --check against the file just written finds no differences
    assert (
        tool.main(
            [
                "--upf",
                str(new_upf),
                "--djrepo",
                str(new_djrepo),
                "--fallback-djrepo",
                str(ref_djrepo),
                "--reference",
                str(reference),
                "--reference-upf",
                str(ref_upf),
                "--output",
                str(output),
                "--check",
            ]
        )
        == 0
    )
