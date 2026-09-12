"""Regression tests: regenerate the reference examples and compare them.

The files under ``examples/`` are the reference output of cif2qewan. These
tests run the generator on the same inputs and require every generated file
to be identical, so that a change which alters a generated Quantum ESPRESSO
or Wannier90 input file cannot pass unnoticed.

``cif2cell`` is not needed: ``qe_wannier_in`` reuses an existing
``cif_scf.in``, and the reference examples ship one.
"""

import os
import re
import shutil
import sys

import pytest

from conftest import EXAMPLES

pytest.importorskip("seekpath", reason="the band k-path needs seekpath")
pytest.importorskip("pymatgen", reason="the band k-path needs pymatgen")

from cif2qewan.cif2qewan import main  # noqa: E402

# (example directory, pseudopotential table) pairs.
EXAMPLE_CASES = [
    ("PSLibrary/Fe", "pp_psl_rrkj.csv"),
    ("PseudoDojo/Fe", "nc-sr-05_pbe_standard_upf.csv"),
]

# The example inputs; everything else in the directory is generated.
INPUT_FILES = {"README.md", "cif2qewan.toml", "mp-13_Fe.cif", "cif_scf.in"}


def generate(example, pp_csv, workdir, repo_root):
    """Regenerate one example inside workdir and return the reference path."""
    reference = EXAMPLES / example

    shutil.copy(reference / "cif_scf.in", workdir)
    shutil.copy(reference / "mp-13_Fe.cif", workdir)

    # Use the example's own configuration, with pp_list_path pointed at the
    # table in this checkout instead of the documented placeholder.
    toml = (reference / "cif2qewan.toml").read_text()
    toml = re.sub(
        r'pp_list_path = ".*"',
        f'pp_list_path = "{repo_root / "cif2qewan" / pp_csv}"',
        toml,
    )
    (workdir / "cif2qewan.toml").write_text(toml)

    # main() reads sys.argv and writes into the current working directory.
    cwd = os.getcwd()
    argv = sys.argv
    os.chdir(workdir)
    sys.argv = ["cif2qewan", "mp-13_Fe.cif", "cif2qewan.toml", "--so", "--mag"]
    try:
        main()
    finally:
        os.chdir(cwd)
        sys.argv = argv

    return reference


def generated_files(root):
    """Every generated file below root, as paths relative to it."""
    return {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path.name not in INPUT_FILES
    }


@pytest.mark.parametrize("example,pp_csv", EXAMPLE_CASES)
def test_example_is_reproduced(example, pp_csv, tmp_path, repo_root):
    """Every generated file must match the reference example exactly."""
    reference = generate(example, pp_csv, tmp_path, repo_root)

    expected_files = generated_files(reference)
    produced_files = generated_files(tmp_path)
    assert produced_files == expected_files

    for name in sorted(expected_files):
        expected = (reference / name).read_text()
        produced = (tmp_path / name).read_text()
        assert produced == expected, f"{example}/{name} differs from the reference"


@pytest.mark.parametrize("example,pp_csv", EXAMPLE_CASES)
def test_nbnd_relations(example, pp_csv, tmp_path, repo_root):
    """nbnd follows the documented relation to num_wann and nexclude.

    With --so/--mag the counts are doubled: nscf uses nexclude + 3*num_wann
    and the shifted check_wannier mesh uses nexclude + 1.5*num_wann.
    """
    generate(example, pp_csv, tmp_path, repo_root)

    def nbnd_of(path):
        for line in (tmp_path / path).read_text().splitlines():
            if "nbnd" in line:
                return int(line.split("=")[1])
        raise AssertionError(f"no nbnd in {path}")

    num_wann = 0
    nexclude = 0
    for line in (tmp_path / "pwscf.win").read_text().splitlines():
        if line.startswith("num_wann"):
            num_wann = int(line.split("=")[1]) // 2  # spinor doubling
        if line.startswith("exclude_bands"):
            nexclude = int(line.split("-")[1]) // 2
    assert num_wann > 0

    assert nbnd_of("nscf.in") == (nexclude + num_wann * 3) * 2
    assert nbnd_of("check_wannier/nscf.in") == (nexclude + int(num_wann * 1.5)) * 2
