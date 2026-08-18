"""Regression tests: regenerate the reference examples and compare them.

The files under ``examples/`` are the reference output of cif2qewan. These
tests run the generator on the same inputs and require the result to be
identical, so that a change which alters a generated Quantum ESPRESSO or
Wannier90 input file cannot pass unnoticed.

``cif2cell`` is not needed: ``qe_wannier_in`` reuses an existing
``cif_scf.in``, and the reference examples ship one.
"""

import os
import shutil

import pytest

from conftest import EXAMPLES, read_pseudo_dir

pytest.importorskip("seekpath", reason="the band k-path needs seekpath")
pytest.importorskip("pymatgen", reason="the band k-path needs pymatgen")

from cif2qewan.cif2qewan import qe_wannier_in  # noqa: E402

# (example directory, pseudopotential table) pairs.
EXAMPLE_CASES = [
    ("PSLibrary/Fe", "pp_psl_rrkj.csv"),
    ("PseudoDojo/Fe", "nc-sr-05_pbe_standard_upf.csv"),
]

# Files that the current code reproduces byte for byte.
#
# band/nscf.in is deliberately absent: the reference files predate the
# adaptive band_nks k-point counts and still carry a uniform 20 points per
# segment. They need to be regenerated before this file can be compared.
COMPARED_FILES = [
    "scf.in",
    "nscf.in",
    "pw2wan.in",
    "pwscf.win",
    "band/band.in",
    "check_wannier/nscf.in",
]

TOML_TEMPLATE = """\
cif2cell_path = "cif2cell"
pseudo_dir = "{pseudo_dir}"
pp_list_path = "{pp_list_path}"
scf_k_resolution = 0.15
degauss = 0.01

[pw2wan]
write_unk = ".true."
"""


def generate(example, pp_csv, workdir, repo_root):
    """Regenerate one example inside workdir and return its path."""
    reference = EXAMPLES / example

    shutil.copy(reference / "cif_scf.in", workdir)
    shutil.copy(reference / "mp-13_Fe.cif", workdir)

    toml_file = workdir / "cif2qewan.toml"
    toml_file.write_text(
        TOML_TEMPLATE.format(
            pseudo_dir=read_pseudo_dir(reference / "scf.in"),
            pp_list_path=repo_root / "cif2qewan" / pp_csv,
        )
    )

    # The generator writes into the current working directory.
    cwd = os.getcwd()
    os.chdir(workdir)
    try:
        qe_wan = qe_wannier_in("mp-13_Fe.cif", str(toml_file), True, True)
        qe_wan.write_pwscf_in("scf.in")
        qe_wan.convert2nscf()
        qe_wan.write_pwscf_in("nscf.in")
        qe_wan.calc_bands_seekpath()
        qe_wan.write_pw2wan("pw2wan.in")
        qe_wan.write_wannier("pwscf.win")

        os.makedirs("check_wannier", exist_ok=True)
        qe_wan.shift_k_nscf()
        qe_wan.write_pwscf_in("check_wannier/nscf.in")

        os.makedirs("band", exist_ok=True)
        qe_wan.convert2band()
        qe_wan.write_pwscf_in("band/nscf.in")
        qe_wan.write_band_in("band/band.in")
    finally:
        os.chdir(cwd)

    return reference


@pytest.mark.parametrize("example,pp_csv", EXAMPLE_CASES)
def test_example_is_reproduced(example, pp_csv, tmp_path, repo_root):
    """The generated inputs must match the reference examples exactly."""
    reference = generate(example, pp_csv, tmp_path, repo_root)

    for name in COMPARED_FILES:
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

    def num_wann_of():
        for line in (tmp_path / "pwscf.win").read_text().splitlines():
            if line.startswith("num_wann"):
                return int(line.split("=")[1])
        raise AssertionError("no num_wann in pwscf.win")

    num_wann = num_wann_of() // 2  # spinor doubling
    nexclude = 0
    for line in (tmp_path / "pwscf.win").read_text().splitlines():
        if line.startswith("exclude_bands"):
            nexclude = int(line.split("-")[1]) // 2

    assert nbnd_of("nscf.in") == (nexclude + num_wann * 3) * 2
    assert nbnd_of("check_wannier/nscf.in") == (nexclude + int(num_wann * 1.5)) * 2
