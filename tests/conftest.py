"""Shared fixtures for the cif2qewan tests."""

import pathlib
import re

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
EXAMPLES = REPO_ROOT / "examples"


@pytest.fixture
def repo_root():
    """Path to the repository root."""
    return REPO_ROOT


def read_pseudo_dir(scf_in):
    """Read the pseudo_dir that a reference scf.in was generated with.

    The reference inputs contain the absolute pseudopotential directory of the
    machine they were produced on. Reusing it keeps the comparison exact
    without hard-coding the path in the test.
    """
    text = pathlib.Path(scf_in).read_text()
    match = re.search(r"pseudo_dir\s*=\s*'([^']*)'", text)
    assert match is not None, f"no pseudo_dir in {scf_in}"
    return match.group(1)
