"""Unit tests for the small parsers used by the analysis scripts."""

import pytest

from cif2qewan.band_comp import get_ef_from_scfout, get_froz_max
from cif2qewan.cif2qewan import pseudo_list
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
    pslist = pseudo_list(
        "/nonexistent", str(repo_root / "cif2qewan" / "pp_psl_rrkj.csv")
    )
    pp_file, nexclude, orbitals, num_wann, ecutwfc, ecutrho = pslist.pseudo("Fe")

    assert pp_file == "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"
    assert nexclude == 4
    assert orbitals == "spd"
    assert num_wann == 1 + 3 + 5
    assert ecutwfc == pytest.approx(64.0)
    assert ecutrho == pytest.approx(782.0)


def test_num_wann_counts_each_orbital_once(repo_root):
    """s/p/d/f contribute 1/3/5/7 Wannier functions."""
    pslist = pseudo_list(
        "/nonexistent", str(repo_root / "cif2qewan" / "pp_psl_rrkj.csv")
    )

    assert pslist.pseudo("O")[3] == 3  # p
    assert pslist.pseudo("Li")[3] == 1  # s
