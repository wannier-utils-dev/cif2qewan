"""Tests for the ``cif2qewan`` command line.

They cover the exit status and the error reporting for representative
invalid inputs, the ``--cif2cell-output`` compatibility path used by the
examples, the ``--reader pymatgen`` path, and the deprecated
``python -m cif2qewan.cif2qewan`` entry point.
"""

import shutil
import subprocess
import sys

import pytest

from conftest import EXAMPLE_CASES, generated_files, write_example_toml

FE = EXAMPLE_CASES[0]  # PSLibrary/Fe, --so --mag
ALL_OUTPUTS = {
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


def run_cli(args, cwd, module="cif2qewan.cli"):
    return subprocess.run(
        [sys.executable, "-m", module, *args],
        cwd=cwd,
        capture_output=True,
        text=True,
    )


def test_help_and_version_exit_zero(tmp_path):
    result = run_cli(["--help"], tmp_path)
    assert result.returncode == 0
    assert "--so" in result.stdout and "--mag" in result.stdout
    assert "--cif2cell-output" in result.stdout and "--reader" in result.stdout
    assert run_cli(["--version"], tmp_path).stdout.startswith("cif2qewan ")


def test_missing_arguments_exit_nonzero(tmp_path):
    result = run_cli([], tmp_path)
    assert result.returncode != 0
    assert "usage" in (result.stdout + result.stderr).lower()


def test_missing_toml_is_reported_without_a_traceback(tmp_path):
    result = run_cli(["structure.cif", "missing.toml"], tmp_path)
    assert result.returncode == 1
    assert "cif2qewan: error:" in result.stderr and "missing.toml" in result.stderr
    assert "Traceback" not in result.stderr
    assert generated_files(tmp_path) == set()


def test_unusable_cif2cell_is_reported(tmp_path):
    """The configured cif2cell path does not exist: no output, clear message."""
    write_example_toml(FE, tmp_path, cif2cell_path="/path/to/cif2cell")
    (tmp_path / "structure.cif").write_text((FE.reference / "mp-13_Fe.cif").read_text())
    result = run_cli(["structure.cif", "cif2qewan.toml"], tmp_path)
    assert result.returncode == 1
    assert (
        "cannot run cif2cell" in result.stderr and "/path/to/cif2cell" in result.stderr
    )
    assert "pip install cif2cell" in result.stderr
    assert "--reader pymatgen" in result.stderr
    assert "Traceback" not in result.stderr
    assert generated_files(tmp_path) == {"structure.cif"}  # only the input we wrote


def test_existing_cif_scf_in_is_not_reused_implicitly(tmp_path):
    """0.2.x silently reused cif_scf.in from the working directory; 0.3 does not."""
    shutil.copy(FE.reference / "cif_scf.in", tmp_path)
    write_example_toml(FE, tmp_path)
    result = run_cli(["does_not_exist.cif", "cif2qewan.toml"], tmp_path)
    assert result.returncode == 1
    assert "CIF file not found" in result.stderr
    assert generated_files(tmp_path) == set()


def test_unsupported_element_is_reported(tmp_path):
    cif_scf = (FE.reference / "cif_scf.in").read_text().replace("Fe", "Fr")
    (tmp_path / "cif_scf.in").write_text(cif_scf)
    write_example_toml(FE, tmp_path)
    result = run_cli(
        ["structure.cif", "cif2qewan.toml", "--cif2cell-output", "cif_scf.in"], tmp_path
    )
    assert result.returncode == 1
    assert "Fr has no pseudopotential" in result.stderr
    assert generated_files(tmp_path) == set()


@pytest.mark.parametrize("flags", [(), ("--so",), ("--mag",), ("--so", "--mag")])
def test_cif2cell_output_path_writes_all_files(tmp_path, flags):
    pytest.importorskip("seekpath")
    pytest.importorskip("pymatgen")
    shutil.copy(FE.reference / "cif_scf.in", tmp_path)
    write_example_toml(FE, tmp_path)

    result = run_cli(
        ["mp-13_Fe.cif", "cif2qewan.toml", *flags, "--cif2cell-output", "cif_scf.in"],
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("cif2qewan: wrote scf.in, ")
    assert generated_files(tmp_path) == ALL_OUTPUTS
    # the explicitly given cif2cell output is left untouched
    assert (tmp_path / "cif_scf.in").read_text() == (
        FE.reference / "cif_scf.in"
    ).read_text()


def test_output_dir_and_copied_cif2cell_output(tmp_path):
    pytest.importorskip("seekpath")
    pytest.importorskip("pymatgen")
    write_example_toml(FE, tmp_path)
    out = tmp_path / "run"
    result = run_cli(
        [
            "x.cif",
            "cif2qewan.toml",
            "--cif2cell-output",
            str(FE.reference / "cif_scf.in"),
            "--output-dir",
            str(out),
        ],
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert generated_files(out) == ALL_OUTPUTS
    assert (out / "cif_scf.in").read_text() == (FE.reference / "cif_scf.in").read_text()
    assert generated_files(tmp_path) == {"run/" + f for f in ALL_OUTPUTS}


def test_pymatgen_reader_path(tmp_path):
    pytest.importorskip("seekpath")
    pytest.importorskip("pymatgen")
    shutil.copy(FE.reference / "mp-13_Fe.cif", tmp_path)
    write_example_toml(FE, tmp_path)
    result = run_cli(
        ["mp-13_Fe.cif", "cif2qewan.toml", "--reader", "pymatgen"], tmp_path
    )
    assert result.returncode == 0, result.stderr
    assert generated_files(tmp_path) == ALL_OUTPUTS  # no cif_scf.in on this path
    scf = (tmp_path / "scf.in").read_text()
    assert "21 21 21  0 0 0" in scf
    assert "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF" in scf


def test_old_module_entry_point_still_works(tmp_path):
    pytest.importorskip("seekpath")
    pytest.importorskip("pymatgen")
    shutil.copy(FE.reference / "cif_scf.in", tmp_path)
    write_example_toml(FE, tmp_path)
    result = run_cli(
        ["x.cif", "cif2qewan.toml", "--cif2cell-output", "cif_scf.in"],
        tmp_path,
        module="cif2qewan.cif2qewan",
    )
    assert result.returncode == 0, result.stderr
    assert generated_files(tmp_path) == ALL_OUTPUTS


def test_verbose_reports_the_decisions_and_special_paths_work(tmp_path):
    pytest.importorskip("seekpath")
    pytest.importorskip("pymatgen")
    workdir = tmp_path / "Fe (bcc) ü"
    workdir.mkdir()
    shutil.copy(FE.reference / "cif_scf.in", workdir / "cif scf.in")
    write_example_toml(FE, workdir)
    result = run_cli(
        [
            "x.cif",
            str(workdir / "cif2qewan.toml"),
            "--cif2cell-output",
            str(workdir / "cif scf.in"),
            "--output-dir",
            str(workdir / "out dir"),
            "--verbose",
        ],
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert "cif2qewan: INFO:" in result.stderr
    assert "scf mesh 21x21x21, nscf mesh 8x8x8" in result.stderr
    assert generated_files(workdir / "out dir") == ALL_OUTPUTS
    # without --verbose the decisions are not reported
    quiet = run_cli(
        [
            "x.cif",
            str(workdir / "cif2qewan.toml"),
            "--cif2cell-output",
            str(workdir / "cif scf.in"),
            "--output-dir",
            str(workdir / "quiet"),
        ],
        tmp_path,
    )
    assert quiet.returncode == 0 and "INFO" not in quiet.stderr
