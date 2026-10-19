import pathlib
import re

"""Regression tests: regenerate the reference examples and compare them.

The files under ``examples/`` are the reference output of cif2qewan. These
tests run the generator on the same inputs and require every generated file
to be identical, so that a change which alters a generated Quantum ESPRESSO
or Wannier90 input file cannot pass unnoticed.
"""

from conftest import generated_files


def test_example_is_reproduced(generated):
    """Every generated file must match the reference example exactly."""
    case, workdir = generated

    expected_files = generated_files(case.reference)
    produced_files = generated_files(workdir)
    assert produced_files == expected_files

    for name in sorted(expected_files):
        expected = (case.reference / name).read_text()
        produced = (workdir / name).read_text()
        assert produced == expected, f"{case.example}/{name} differs from the reference"


def test_nbnd_relations(generated):
    """nbnd follows the documented relation to num_wann and nexclude.

    nscf uses nexclude + 3*num_wann and the shifted check_wannier mesh uses
    nexclude + int(1.5*num_wann); with --so or --mag both are doubled.
    """
    case, workdir = generated

    def nbnd_of(path):
        for line in (workdir / path).read_text().splitlines():
            if "nbnd" in line:
                return int(line.split("=")[1])
        raise AssertionError(f"no nbnd in {path}")

    num_wann = 0
    nexclude = 0
    for line in (workdir / "pwscf.win").read_text().splitlines():
        if line.startswith("num_wann"):
            num_wann = int(line.split("=")[1]) // case.spin_factor
        if line.startswith("exclude_bands"):
            nexclude = int(line.split("-")[1]) // case.spin_factor
    assert num_wann > 0

    factor = case.spin_factor
    assert nbnd_of("nscf.in") == (nexclude + num_wann * 3) * factor
    assert nbnd_of("check_wannier/nscf.in") == (nexclude + int(num_wann * 1.5)) * factor


# --------------------------------------------------------------------------
# 0.2.0 outputs (tests/fixtures/v0.2.0): the generated inputs differ from the
# 0.2.0 ones only in the K_POINTS spelling of scf.in and the layout of
# pwscf.win; every other file is byte for byte the 0.2.0 output, which the
# reference examples pin above.
# --------------------------------------------------------------------------

V020 = pathlib.Path(__file__).parent / "fixtures" / "v0.2.0"


def _normalized(text):
    """Lines with collapsed whitespace, without blank lines, sorted."""
    return sorted(
        re.sub(r"\s+", " ", line).strip() for line in text.splitlines() if line.strip()
    )


def test_scf_in_matches_0_2_0_up_to_the_kpoints_spelling(generated):
    case, workdir = generated
    old = (V020 / case.example / "scf.in").read_text()
    new = (workdir / "scf.in").read_text()
    assert "K_POINTS automatic\n" in old and "K_POINTS {automatic}\n" in new
    assert old.replace("K_POINTS automatic\n", "K_POINTS {automatic}\n") == new


def test_win_matches_0_2_0_up_to_layout(generated):
    case, workdir = generated
    old = (V020 / case.example / "pwscf.win").read_text()
    new = (workdir / "pwscf.win").read_text()
    assert old != new  # the layout did change ...
    assert _normalized(old) == _normalized(new)  # ... and nothing else
