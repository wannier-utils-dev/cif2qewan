"""Unit tests for the parsers of the analysis scripts and the pseudopotential table."""

import pytest

from cif2qewan.band_comp import get_ef_from_scfout, get_froz_max
from cif2qewan.qe.pseudopotential import PseudopotentialTable
from cif2qewan.wannier_conv import get_nexclude


def test_fermi_energy_is_the_last_one(tmp_path):
    """scf.out can hold several Fermi lines; the converged one is the last.

    Regression test: a premature break made this return the first value.
    """
    scfout = tmp_path / "scf.out"
    scfout.write_text(
        "     the Fermi energy is    11.1234 ev\n"
        "     some other output\n"
        "     the Fermi energy is    12.5678 ev\n"
    )

    assert get_ef_from_scfout(str(scfout)) == pytest.approx(12.5678)


def test_froz_max_is_the_last_one(tmp_path):
    """The last dis_froz_max assignment in pwscf.win is the effective one.

    Regression test for the same premature break.
    """
    win = tmp_path / "pwscf.win"
    win.write_text("dis_froz_max = -200\ndis_froz_max = 13.5\n")

    assert get_froz_max(str(win)) == pytest.approx(13.5)


def test_froz_max_default_when_absent(tmp_path):
    win = tmp_path / "pwscf.win"
    win.write_text("num_wann = 18\n")

    assert get_froz_max(str(win)) == pytest.approx(-200.0)


def test_get_nexclude(tmp_path):
    win = tmp_path / "pwscf.win"
    win.write_text("num_wann = 18\nexclude_bands = 1-8\n")

    assert get_nexclude(str(win)) == 8


def test_get_nexclude_defaults_to_zero(tmp_path):
    win = tmp_path / "pwscf.win"
    win.write_text("num_wann = 18\n")

    assert get_nexclude(str(win)) == 0


def test_pseudo_table_is_read(repo_root):
    """The shipped PSLibrary table maps an element to its pseudopotential."""
    table = PseudopotentialTable.from_csv(repo_root / "cif2qewan" / "pp_psl_rrkj.csv")
    entry = table.lookup("Fe")

    assert entry.file_name == "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"
    assert entry.nexclude == 4
    assert entry.orbitals == "spd"
    assert entry.num_wann == 1 + 3 + 5
    assert entry.ecutwfc == pytest.approx(64.0)
    assert entry.ecutrho == pytest.approx(782.0)


def test_num_wann_counts_each_orbital_once(repo_root):
    """s/p/d/f contribute 1/3/5/7 Wannier functions."""
    table = PseudopotentialTable.from_csv(repo_root / "cif2qewan" / "pp_psl_rrkj.csv")

    assert table.lookup("O").num_wann == 3  # p
    assert table.lookup("Li").num_wann == 1  # s
