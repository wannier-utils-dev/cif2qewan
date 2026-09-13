"""Shared fixtures for the cif2qewan tests.

The reference examples under ``examples/`` are the characterization baseline
of the input generator. ``generate`` re-runs the generator on one example in a
scratch directory; the ``generated`` fixture does this once per example for
the whole test session so that several test modules can inspect the result.

``cif2cell`` is not needed: ``qe_wannier_in`` reuses an existing
``cif_scf.in``, and every reference example ships one.
"""

import os
import pathlib
import re
import toml
import shutil
import sys
from typing import NamedTuple

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
EXAMPLES = REPO_ROOT / "examples"
PACKAGE = REPO_ROOT / "cif2qewan"

# The example inputs; everything else in an example directory is generated.
INPUT_FILES = {"README.md", "cif2qewan.toml", "mp-13_Fe.cif", "cif_scf.in"}


class ExampleCase(NamedTuple):
    """One reference example: directory, pseudopotential table and CLI flags."""

    example: str
    pp_csv: str
    flags: tuple

    @property
    def reference(self):
        return EXAMPLES / self.example

    @property
    def so(self):
        return "--so" in self.flags

    @property
    def mag(self):
        return "--mag" in self.flags

    @property
    def spin_factor(self):
        """Band and Wannier-function counts double for spinor calculations."""
        return 2 if (self.so or self.mag) else 1


EXAMPLE_CASES = [
    ExampleCase("PSLibrary/Fe", "pp_psl_rrkj.csv", ("--so", "--mag")),
    ExampleCase("PSLibrary/Fe_nonmag", "pp_psl_rrkj.csv", ()),
    ExampleCase("PSLibrary/Fe_so", "pp_psl_rrkj.csv", ("--so",)),
    ExampleCase("PseudoDojo/Fe", "nc-sr-05_pbe_standard_upf.csv", ("--so", "--mag")),
]


def write_example_toml(case, workdir):
    """Copy the example's TOML with pp_list_path pointed at this checkout."""
    toml = (case.reference / "cif2qewan.toml").read_text()
    toml = re.sub(
        r'pp_list_path = ".*"',
        f'pp_list_path = "{PACKAGE / case.pp_csv}"',
        toml,
    )
    (workdir / "cif2qewan.toml").write_text(toml)


def generate(case, workdir):
    """Regenerate one example inside workdir with the current code."""
    pytest.importorskip("seekpath", reason="the band k-path needs seekpath")
    pytest.importorskip("pymatgen", reason="the band k-path needs pymatgen")
    from cif2qewan.cif2qewan import main

    shutil.copy(case.reference / "cif_scf.in", workdir)
    shutil.copy(case.reference / "mp-13_Fe.cif", workdir)
    write_example_toml(case, workdir)

    # main() reads sys.argv and writes into the current working directory.
    cwd = os.getcwd()
    argv = sys.argv
    os.chdir(workdir)
    sys.argv = ["cif2qewan", "mp-13_Fe.cif", "cif2qewan.toml", *case.flags]
    try:
        main()
    finally:
        os.chdir(cwd)
        sys.argv = argv


def generated_files(root):
    """Every generated file below root, as paths relative to it."""
    return {
        str(path.relative_to(root))
        for path in root.rglob("*")
        if path.is_file() and path.name not in INPUT_FILES
    }


@pytest.fixture
def repo_root():
    """Path to the repository root."""
    return REPO_ROOT


@pytest.fixture(scope="session", params=EXAMPLE_CASES, ids=lambda c: c.example)
def generated(request, tmp_path_factory):
    """(case, workdir) for each reference example, regenerated once per session."""
    case = request.param
    workdir = tmp_path_factory.mktemp(case.example.replace("/", "_"))
    generate(case, workdir)
    return case, workdir


def example_config(case, so=None, mag=None):
    """The example's Config with pp_list_path pointed at this checkout."""
    from cif2qewan.config import Config

    data = toml.load(case.reference / "cif2qewan.toml")
    data["pp_list_path"] = str(PACKAGE / case.pp_csv)
    return Config.from_dict(
        data,
        so=case.so if so is None else so,
        mag=case.mag if mag is None else mag,
    )


def win_blocks(text):
    """name -> list of lines for every begin/end block of a .win file."""
    return {
        name: body.strip("\n").splitlines()
        for name, body in re.findall(r"begin (\w+)\n(.*?)end \1", text, re.S)
    }


def win_keywords(text):
    """key -> value for every ``key = value`` / ``key: value`` line outside blocks."""
    outside = re.sub(r"begin (\w+)\n.*?end \1\n", "", text, flags=re.S)
    return dict(re.findall(r"^(\w+)\s*[=:]\s*(.*?)\s*$", outside, re.M))


def normalize_reference(path, text):
    """Reference text with the documented 0.3 formatting normalizations applied."""
    if path == "scf.in":
        text = text.replace("K_POINTS automatic", "K_POINTS {automatic}")
    return text
