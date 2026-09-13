"""Tests for MCIF magnetism (DEVELOPMENT_PLAN.md Step 9).

Covers the moment-aware pymatgen reader, the tolerance-based site
classification, the collinear / noncollinear decision, the QE
starting_magnetization / angle1 / angle2 values, and the resulting QE and
Wannier90 inputs for a noncollinear (Mn3Sn) and a collinear
(antiferromagnetic bcc Fe) structure.
"""

import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

from cif2qewan.config import Config
from cif2qewan.exceptions import StructureError
from cif2qewan.qe.model import RawValue
from cif2qewan.structure.magnetism import (
    COLLINEAR,
    NONCOLLINEAR,
    NONMAGNETIC,
    classify_magnetic_sites,
    magnetic_order,
    qe_magnetization,
)
from cif2qewan.structure.model import AtomicSite, MagneticMoment, NormalizedStructure
from cif2qewan.structure.readers import PymatgenReader
from cif2qewan.workflow.builder import SpinPolicy, build_plan
from cif2qewan.workflow.render import render_plan
from conftest import EXAMPLE_CASES, PACKAGE, example_config, write_example_toml

pytest.importorskip("pymatgen")

FIXTURES = EXAMPLE_CASES[0].reference.parent.parent.parent / "tests" / "fixtures"
MN3SN = FIXTURES / "Mn3Sn.mcif"
A_FE = 2.8630355


def up(m=2.2):
    return MagneticMoment((0.0, 0.0, m))


def cubic(sites):
    return NormalizedStructure(np.eye(3) * A_FE, sites)


def afm_fe():
    return cubic(
        (AtomicSite("Fe", (0, 0, 0), up()), AtomicSite("Fe", (0.5, 0.5, 0.5), up(-2.2)))
    )


def fm_fe():
    return cubic(
        (AtomicSite("Fe", (0, 0, 0), up()), AtomicSite("Fe", (0.5, 0.5, 0.5), up()))
    )


def write_mcif(structure, path):
    from pymatgen.core import Lattice, Structure
    from pymatgen.electronic_structure.core import Magmom
    from pymatgen.io.cif import CifWriter

    pmg = Structure(
        Lattice(structure.lattice_matrix),
        [s.element for s in structure.sites],
        [s.frac_coords for s in structure.sites],
        site_properties={
            "magmom": [Magmom(list(s.magmom.vector)) for s in structure.sites]
        },
    )
    CifWriter(pmg, write_magmoms=True).write_file(str(path))


# --------------------------------------------------------------------------
# order and classification
# --------------------------------------------------------------------------


def test_magnetic_order():
    assert magnetic_order(cubic((AtomicSite("Fe", (0, 0, 0)),))) == NONMAGNETIC
    assert (
        magnetic_order(cubic((AtomicSite("Fe", (0, 0, 0), MagneticMoment.zero()),)))
        == NONMAGNETIC
    )
    assert magnetic_order(afm_fe()) == COLLINEAR
    assert magnetic_order(fm_fe()) == COLLINEAR
    tilted = cubic(
        (
            AtomicSite("Fe", (0, 0, 0), up()),
            AtomicSite("Fe", (0.5, 0.5, 0.5), MagneticMoment.from_polar(2.2, 4.0, 0.0)),
        )
    )
    assert magnetic_order(tilted) == COLLINEAR  # within the 5 degree tolerance
    assert magnetic_order(tilted, angle_tol_deg=1.0) == NONCOLLINEAR
    canted = cubic(
        (
            AtomicSite("Fe", (0, 0, 0), up()),
            AtomicSite("Fe", (0.5, 0.5, 0.5), MagneticMoment((2.2, 0, 0))),
        )
    )
    assert magnetic_order(canted) == NONCOLLINEAR


def test_classification_splits_by_moment_with_tolerance():
    relabelled, species = classify_magnetic_sites(afm_fe())
    assert [s.label for s in relabelled.sites] == ["Fe1", "Fe2"]
    assert [(s.label, s.element, s.count) for s in species] == [
        ("Fe1", "Fe", 1),
        ("Fe2", "Fe", 1),
    ]
    assert species[1].moment.vector == pytest.approx((0.0, 0.0, -2.2))

    # moments equal within the tolerance stay one species and keep the plain label
    nearly = cubic(
        (
            AtomicSite("Fe", (0, 0, 0), up(2.2)),
            AtomicSite("Fe", (0.5, 0.5, 0.5), up(2.2004)),
        )
    )
    relabelled, species = classify_magnetic_sites(nearly)
    assert [s.label for s in relabelled.sites] == ["Fe", "Fe"]
    assert len(species) == 1 and species[0].count == 2

    # a site without a moment counts as zero moment and is its own species
    mixed = cubic(
        (AtomicSite("Fe", (0, 0, 0), up()), AtomicSite("Fe", (0.5, 0.5, 0.5)))
    )
    relabelled, species = classify_magnetic_sites(mixed)
    assert [s.label for s in species] == ["Fe1", "Fe2"]
    assert species[1].moment.is_zero()

    # existing labels are refined, not replaced
    labelled = cubic(
        (
            AtomicSite("Fe", (0, 0, 0), up(), label="Fea"),
            AtomicSite("Fe", (0.5, 0.5, 0.5), up(-2.2), label="Fea"),
            AtomicSite("Fe", (0.5, 0.5, 0.0), up(), label="Feb"),
        )
    )
    relabelled, species = classify_magnetic_sites(labelled)
    assert [s.label for s in species] == ["Fea1", "Fea2", "Feb"]


def test_classification_rejects_label_clashes():
    clashing = cubic(
        (
            AtomicSite("Fe", (0, 0, 0), up()),
            AtomicSite("Fe", (0.5, 0.5, 0.5), up(-2.2)),
            AtomicSite("Fe", (0.5, 0.5, 0.0), up(), label="Fe1"),
        )
    )
    with pytest.raises(StructureError, match="clashes"):
        classify_magnetic_sites(clashing)


def test_qe_magnetization_values():
    _, species = classify_magnetic_sites(afm_fe())
    mag = qe_magnetization(species, COLLINEAR)
    assert mag["Fe1"].starting_magnetization == pytest.approx(1.0)
    assert mag["Fe2"].starting_magnetization == pytest.approx(-1.0)
    assert (mag["Fe1"].angle1, mag["Fe1"].angle2) == (0.0, 0.0)
    assert (mag["Fe2"].angle1, mag["Fe2"].angle2) == (
        0.0,
        0.0,
    )  # common axis, sign carries the direction

    canted = cubic(
        (
            AtomicSite("Fe", (0, 0, 0), MagneticMoment((0, 0, 3.0))),
            AtomicSite("Fe", (0.5, 0.5, 0.5), MagneticMoment((1.5, 0, 0))),
            AtomicSite("Sn", (0.5, 0, 0)),
        )
    )
    _, species = classify_magnetic_sites(canted)
    mag = qe_magnetization(species, NONCOLLINEAR)
    assert mag["Fe1"].starting_magnetization == pytest.approx(1.0)
    assert mag["Fe2"].starting_magnetization == pytest.approx(0.5)
    assert (mag["Fe2"].angle1, mag["Fe2"].angle2) == pytest.approx((90.0, 0.0))
    assert mag["Sn"].is_zero and (mag["Sn"].angle1, mag["Sn"].angle2) == (0.0, 0.0)

    assert all(m.is_zero for m in qe_magnetization(species, NONMAGNETIC).values())


# --------------------------------------------------------------------------
# reader
# --------------------------------------------------------------------------


def test_reader_keeps_mn3sn_moments_and_frame():
    structure = PymatgenReader().read(MN3SN)
    assert structure.num_sites == 8
    assert structure.elements == ("Mn", "Sn")
    assert structure.lattice[0] == pytest.approx(
        (5.665, 0.0, 0.0), abs=1e-6
    )  # input frame kept
    moments = [s.magmom for s in structure.sites if s.element == "Mn"]
    assert all(m.magnitude == pytest.approx(3.0, abs=1e-3) for m in moments)
    assert all(abs(m.vector[2]) < 1e-6 for m in moments)  # in-plane triangular order
    assert all(s.magmom.is_zero() for s in structure.sites if s.element == "Sn")
    assert magnetic_order(structure) == NONCOLLINEAR


def test_reader_does_not_fold_an_antiferromagnet(tmp_path):
    afm_path = tmp_path / "afm.mcif"
    write_mcif(afm_fe(), afm_path)
    structure = PymatgenReader().read(afm_path)
    assert structure.num_sites == 2  # the chemical primitive cell would have one site
    assert sorted(round(s.magmom.vector[2], 3) for s in structure.sites) == [-2.2, 2.2]
    assert structure.volume == pytest.approx(A_FE**3)


def test_reader_reduces_a_ferromagnet_and_keeps_its_moment(tmp_path):
    fm_path = tmp_path / "fm.mcif"
    write_mcif(fm_fe(), fm_path)
    structure = PymatgenReader().read(fm_path)
    assert structure.num_sites == 1
    assert structure.volume == pytest.approx(A_FE**3 / 2)
    assert structure.sites[0].magmom.vector == pytest.approx((0.0, 0.0, 2.2))


# --------------------------------------------------------------------------
# builder
# --------------------------------------------------------------------------


def psl_config(so=False, mag=False):
    return Config.from_dict(
        {
            "cif2cell_path": "c",
            "pseudo_dir": "/pp",
            "pp_list_path": str(PACKAGE / "pp_psl_rrkj.csv"),
            "scf_k_resolution": 0.15,
            "degauss": 0.01,
            "pw2wan": {"write_unk": ".true."},
        },
        so=so,
        mag=mag,
    )


def test_spin_policy_matrix():
    assert not SpinPolicy(NONMAGNETIC, {}, so=False, mag=False).spinor
    assert SpinPolicy(NONMAGNETIC, {}, so=True, mag=False).spinor
    assert SpinPolicy(NONMAGNETIC, {}, so=False, mag=True).collinear_scf
    assert not SpinPolicy(NONMAGNETIC, {}, so=True, mag=False).collinear_scf
    _, species = classify_magnetic_sites(afm_fe())
    collinear = SpinPolicy(
        COLLINEAR, qe_magnetization(species, COLLINEAR), so=True, mag=False
    )
    assert collinear.spinor and collinear.collinear_scf
    assert not collinear.scf_relativistic and collinear.nscf_relativistic
    assert collinear.scf(["Fe1", "Fe2"])["nspin"] == 2
    assert collinear.nscf(["Fe1", "Fe2"])["lforcet"] is True


def test_mn3sn_plan_is_noncollinear_from_the_scf(tmp_path):
    structure = PymatgenReader().read(MN3SN)
    plan = build_plan(structure, psl_config(so=True))

    labels = [s.label for s in plan.scf.species]
    assert labels == ["Mn1", "Mn2", "Mn3", "Sn"]
    assert [p.label for p in plan.scf.positions] == [
        "Mn1",
        "Mn2",
        "Mn1",
        "Mn2",
        "Mn3",
        "Mn3",
        "Sn",
        "Sn",
    ]
    assert all(s.pseudo_file.startswith("Mn.rel-") for s in plan.scf.species[:3])

    scf = plan.scf.system.entries
    assert scf["noncolin"] is True and scf["lspinorb"] is True and "nspin" not in scf
    assert scf["starting_magnetization(1)"] == RawValue("1.0000")
    assert scf["starting_magnetization(4)"] == RawValue("0.0000")
    assert (scf["angle1(1)"], scf["angle2(1)"]) == (
        RawValue("90.0000"),
        RawValue("60.0000"),
    )
    assert (scf["angle1(2)"], scf["angle2(2)"]) == (
        RawValue("90.0000"),
        RawValue("180.0000"),
    )
    assert (scf["angle1(3)"], scf["angle2(3)"]) == (
        RawValue("90.0000"),
        RawValue("-60.0000"),
    )
    assert "angle1(4)" not in scf  # Sn carries no moment
    nscf = plan.nscf.system.entries
    assert "lforcet" not in nscf and nscf["noncolin"] is True
    assert nscf["angle2(3)"] == RawValue("-60.0000")

    win = plan.wannier90
    assert win.spinors
    assert (win.num_wann, win.num_bands, win.exclude_bands) == (124, 372, (1, 68))
    assert [p.site for p in win.projections] == ["Mn1", "Mn2", "Mn3", "Sn"]
    assert [a.label for a in win.atoms_frac] == [p.label for p in plan.scf.positions]

    # the moments survive as direction and relative magnitude in every QE input
    files = render_plan(plan)
    for name in ("scf.in", "nscf.in", "check_wannier/nscf.in", "band/nscf.in"):
        assert "  angle2(3) = -60.0000\n" in files[name], name
    assert "Mn1:s,p,d\nMn2:s,p,d\nMn3:s,p,d\nSn:s,p\n" in files["pwscf.win"]


def test_afm_plan_is_collinear_scf_then_noncollinear_nscf():
    plan = build_plan(afm_fe(), psl_config(so=False))
    scf = plan.scf.system.entries
    assert scf["nspin"] == 2 and "noncolin" not in scf
    assert scf["starting_magnetization(1)"] == RawValue("1.0000")
    assert scf["starting_magnetization(2)"] == RawValue("-1.0000")
    assert "angle1(1)" not in scf
    nscf = plan.nscf.system.entries
    assert (
        nscf["noncolin"] is True
        and nscf["lforcet"] is True
        and nscf["lspinorb"] is False
    )
    assert nscf["starting_magnetization(2)"] == RawValue("-1.0000")
    assert (nscf["angle1(1)"], nscf["angle2(1)"]) == (
        RawValue("0.0000"),
        RawValue("0.0000"),
    )
    assert plan.wannier90.spinors and plan.wannier90.num_wann == 2 * 9 * 2
    assert (
        plan.scf.species[0].pseudo_file == "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"
    )  # scalar rel.


def test_structure_moments_take_precedence_over_mag_flag():
    with_flag = build_plan(afm_fe(), psl_config(mag=True))
    without = build_plan(afm_fe(), psl_config(mag=False))
    assert with_flag.scf.system == without.scf.system
    assert with_flag.nscf.system == without.nscf.system


def test_legacy_mag_examples_are_unchanged():
    """The ferromagnetic --mag policy for structures without moments is untouched."""
    case = EXAMPLE_CASES[0]
    from cif2qewan.structure.readers import Cif2cellReader
    from cif2qewan.workflow.builder import plan_from_cif2cell_output

    plan = plan_from_cif2cell_output(
        Cif2cellReader.read_output(case.reference / "cif_scf.in"), example_config(case)
    )
    for path, text in render_plan(plan).items():
        assert text == (case.reference / path).read_text(), path


def test_cli_reads_an_mcif_with_the_pymatgen_reader(tmp_path):
    pytest.importorskip("seekpath")
    shutil.copy(MN3SN, tmp_path / "Mn3Sn.mcif")
    write_example_toml(EXAMPLE_CASES[0], tmp_path)
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "cif2qewan.cli",
            "Mn3Sn.mcif",
            "cif2qewan.toml",
            "--so",
            "--reader",
            "pymatgen",
            "-v",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (
        "noncollinear magnetic order from the structure; species Mn1(2), Mn2(2), Mn3(2), Sn(2)"
        in result.stderr
    )
    scf = (tmp_path / "scf.in").read_text()
    assert "  ntyp = 4\n" in scf and "  angle2(2) = 180.0000\n" in scf
    assert re.search(r"^  Mn3 .*Mn\.rel-pbe", scf, re.M)
