"""Shared fixtures for the cif2qewan tests."""

import pathlib

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
EXAMPLES = REPO_ROOT / "examples"


@pytest.fixture
def repo_root():
    """Path to the repository root."""
    return REPO_ROOT
