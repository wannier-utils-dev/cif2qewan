"""Tests for the output boundary and the analysis-script error handling (Step 8)."""

import subprocess
import sys

import pytest

from cif2qewan.exceptions import Cif2qewanError, DataFileError
from cif2qewan.wannier_conv import Hamiltonian
from cif2qewan.workflow.output import write_rendered


def test_write_rendered_creates_subdirectories_and_overwrites(tmp_path):
    out = tmp_path / "run (1) ü"
    written = write_rendered({"scf.in": "a\n", "band/nscf.in": "b\n"}, out)
    assert [p.relative_to(out).as_posix() for p in written] == [
        "scf.in",
        "band/nscf.in",
    ]
    assert (out / "band" / "nscf.in").read_text() == "b\n"
    write_rendered({"scf.in": "c\n"}, out)
    assert (out / "scf.in").read_text() == "c\n"


def test_write_rendered_refuses_paths_outside_the_directory(tmp_path):
    with pytest.raises(ValueError, match="outside"):
        write_rendered({"../escape.in": "x\n"}, tmp_path / "out")
    assert not (tmp_path / "escape.in").exists()


def test_hamiltonian_reader_raises_a_typed_error(tmp_path):
    with pytest.raises(DataFileError, match="cannot read the Wannier90 Hamiltonian"):
        Hamiltonian(str(tmp_path / "missing_hr.dat"))
    bad = tmp_path / "bad_hr.dat"
    bad.write_text("header\nnot a number\n")
    with pytest.raises(DataFileError):
        Hamiltonian(str(bad))
    assert issubclass(DataFileError, Cif2qewanError)


def run_module(module, args, cwd):
    return subprocess.run(
        [sys.executable, "-m", module, *args], cwd=cwd, capture_output=True, text=True
    )


def test_band_comp_reports_missing_files(tmp_path):
    result = run_module(
        "cif2qewan.band_comp", ["-o", str(tmp_path / "missing")], tmp_path
    )
    assert result.returncode == 1
    assert result.stderr.startswith("band_comp: error:")
    assert "Traceback" not in result.stderr


def test_band_comp_reports_missing_scf_after_reading_band_files(tmp_path):
    work = tmp_path / "wannier"
    work.mkdir()
    band = tmp_path / "band"
    band.mkdir()
    data = "0.0 1.0\n1.0 2.0\n"
    (work / "pwscf_band.dat").write_text(data)
    (band / "bands.out.gnu").write_text(data)

    result = run_module("cif2qewan.band_comp", ["-o", str(tmp_path)], work)
    assert result.returncode == 1
    assert result.stderr.startswith("band_comp: error:")
    assert "scf.out" in result.stderr
    assert "Traceback" not in result.stderr


def test_wannier_conv_reports_missing_files(tmp_path):
    result = run_module(
        "cif2qewan.wannier_conv",
        ["-e", "5.0", "-o", str(tmp_path), "-i", "missing.out"],
        tmp_path,
    )
    assert result.returncode == 1
    assert result.stderr.startswith("wannier_conv: error:")
    assert "Traceback" not in result.stderr
