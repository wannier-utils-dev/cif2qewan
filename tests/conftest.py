"""Shared fixtures for the cif2qewan tests.

The reference examples under ``examples/`` are the characterization baseline
of the input generator. ``generate`` re-runs the generator on one example in a
scratch directory; the ``generated`` fixture does this once per example for
the whole test session so that several test modules can inspect the result.

``cif2cell`` is not needed: the shipped ``cif_scf.in`` of every example is
passed to the CLI with ``--cif2cell-output``.
"""

import pathlib
import shutil
from typing import NamedTuple

import pytest
import toml

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


def write_example_toml(case, workdir, **extra):
    """Copy the example's TOML into workdir, with ``key = "value"`` lines added.

    The extra keys go before the tables. The examples select their
    pseudopotential table by its bundled file name (or use the default), so
    nothing needs to be rewritten.
    """
    text = (case.reference / "cif2qewan.toml").read_text()
    text = "".join(f'{key} = "{value}"\n' for key, value in extra.items()) + text
    (workdir / "cif2qewan.toml").write_text(text)


def generate(case, workdir):
    """Regenerate one example inside workdir with the current code.

    The shipped cif2cell output is passed explicitly, so cif2cell itself is
    not needed; the CIF is copied only to mirror the example directory.
    """
    pytest.importorskip("seekpath", reason="the band k-path needs seekpath")
    pytest.importorskip("pymatgen", reason="the band k-path needs pymatgen")
    from cif2qewan.cli import main

    shutil.copy(case.reference / "cif_scf.in", workdir)
    shutil.copy(case.reference / "mp-13_Fe.cif", workdir)
    write_example_toml(case, workdir)

    status = main(
        [
            "mp-13_Fe.cif",
            str(workdir / "cif2qewan.toml"),
            *case.flags,
            "--cif2cell-output",
            str(workdir / "cif_scf.in"),
            "--output-dir",
            str(workdir),
        ]
    )
    assert status == 0


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
    """The example's Config; its table is the bundled ``case.pp_csv``."""
    from cif2qewan.config import Config

    data = toml.load(case.reference / "cif2qewan.toml")
    return Config.from_dict(
        data,
        so=case.so if so is None else so,
        mag=case.mag if mag is None else mag,
    )
