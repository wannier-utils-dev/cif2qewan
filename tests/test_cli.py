"""Characterization tests for the command-line behaviour of ``cif2qewan``.

They record the exit status and the failure mode of the 0.2.x CLI for
representative invalid inputs. Some of these behaviours (the bare
tracebacks, the silent reuse of ``cif_scf.in``) are limitations that the
0.3.0 refactoring intends to change; when that happens, update the test in
the same PR and record the change in CHANGELOG.md.
"""

import shutil
import subprocess
import sys

import pytest

from conftest import EXAMPLES, EXAMPLE_CASES, generated_files, write_example_toml

FE = EXAMPLE_CASES[0]  # PSLibrary/Fe, --so --mag


def run_cli(args, cwd):
    return subprocess.run(
        [sys.executable, "-m", "cif2qewan.cif2qewan", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def test_help_exits_zero(tmp_path):
    result = run_cli(["--help"], tmp_path)
    assert result.returncode == 0
    assert "cif2qewan" in result.stdout
    assert "--so" in result.stdout and "--mag" in result.stdout


def test_missing_arguments_exit_nonzero(tmp_path):
    result = run_cli([], tmp_path)
    assert result.returncode != 0
    assert "Usage" in result.stdout + result.stderr


def test_missing_toml_exits_nonzero(tmp_path):
    result = run_cli(["structure.cif", "missing.toml"], tmp_path)
    assert result.returncode != 0
    assert "missing.toml" in result.stderr
    assert generated_files(tmp_path) == set()


def test_missing_cif_without_cif2cell_exits_nonzero(tmp_path):
    """With no cif_scf.in and an unusable cif2cell path nothing is generated."""
    write_example_toml(FE, tmp_path)
    result = run_cli(["missing.cif", "cif2qewan.toml"], tmp_path)
    assert result.returncode != 0
    assert "cif_scf.in" in result.stderr
    assert generated_files(tmp_path) == set()


def test_unsupported_element_exits_nonzero(tmp_path):
    """An element without a pseudopotential in the table aborts the run."""
    pytest.importorskip("pymatgen")
    cif_scf = (FE.reference / "cif_scf.in").read_text().replace("Fe", "Fr")
    (tmp_path / "cif_scf.in").write_text(cif_scf)
    write_example_toml(FE, tmp_path)

    result = run_cli(["structure.cif", "cif2qewan.toml"], tmp_path)
    assert result.returncode != 0
    assert "scf.in" not in generated_files(tmp_path)


@pytest.mark.parametrize("flags", [(), ("--so",), ("--mag",), ("--so", "--mag")])
def test_valid_run_exits_zero_and_writes_all_files(tmp_path, flags):
    pytest.importorskip("seekpath")
    pytest.importorskip("pymatgen")
    shutil.copy(FE.reference / "cif_scf.in", tmp_path)
    shutil.copy(FE.reference / "mp-13_Fe.cif", tmp_path)
    write_example_toml(FE, tmp_path)

    result = run_cli(["mp-13_Fe.cif", "cif2qewan.toml", *flags], tmp_path)
    assert result.returncode == 0, result.stderr
    assert generated_files(tmp_path) >= {
        "scf.in",
        "nscf.in",
        "pw2wan.in",
        "pwscf.win",
        "check_wannier/nscf.in",
        "band/nscf.in",
        "band/band.in",
        "band/proj.in",
        "band/pp.in",
    }


def test_existing_cif_scf_in_is_reused_without_the_cif(tmp_path):
    """0.2.x reuses cif_scf.in from the working directory unconditionally.

    The CIF file itself is not read when cif_scf.in exists. This is the
    behaviour the example tests rely on; DEVELOPMENT_PLAN.md Step 3 replaces
    it with an explicit way to pass a pre-generated cif2cell output.
    """
    pytest.importorskip("seekpath")
    pytest.importorskip("pymatgen")
    shutil.copy(FE.reference / "cif_scf.in", tmp_path)
    write_example_toml(FE, tmp_path)

    result = run_cli(["does_not_exist.cif", "cif2qewan.toml"], tmp_path)
    assert result.returncode == 0, result.stderr
    assert "scf.in" in generated_files(tmp_path)
