"""Unit tests for the typed internal models (DEVELOPMENT_PLAN.md Step 2).

The models carry units and validation only; nothing here runs the input
generator or touches the file system.
"""

import math

import numpy as np
import pytest

from cif2qewan.exceptions import (
    Cif2qewanError,
    ExternalCommandError,
    InputModelError,
    StructureError,
)
from cif2qewan.qe.model import (
    AtomicPosition,
    CellParameters,
    KPathPoint,
    KPointsAutomatic,
    KPointsList,
    KPointsPath,
    Namelist,
    NamelistInput,
    PwInput,
    RawValue,
    Species,
)
from cif2qewan.structure.model import AtomicSite, MagneticMoment, NormalizedStructure
from cif2qewan.wannier90.model import AtomFrac, KPathSegment, Projection, Wannier90Input
from cif2qewan.workflow.model import CalculationPlan

# bcc Fe, primitive cell, as in examples/PSLibrary/Fe/cif_scf.in
ALAT = 2.86304
BCC = np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]) * ALAT


def fe_structure(magmom=None):
    return NormalizedStructure(BCC, (AtomicSite("Fe", (0.0, 0.0, 0.0), magmom),))


# --------------------------------------------------------------------------
# exceptions
# --------------------------------------------------------------------------


def test_exceptions_share_a_base_class():
    for exc in (StructureError, InputModelError, ExternalCommandError):
        assert issubclass(exc, Cif2qewanError)


def test_external_command_error_message():
    err = ExternalCommandError(["cif2cell", "x.cif"], returncode=2, stderr="boom\n")
    assert err.command == ("cif2cell", "x.cif")
    assert "exit status 2" in str(err) and "boom" in str(err)


# --------------------------------------------------------------------------
# MagneticMoment
# --------------------------------------------------------------------------


def test_moment_magnitude_and_direction():
    m = MagneticMoment((0.0, 0.0, 3.0))
    assert m.magnitude == pytest.approx(3.0)
    assert m.direction() == pytest.approx([0.0, 0.0, 1.0])
    assert not m.is_zero()
    assert MagneticMoment.zero().is_zero()


def test_moment_polar_roundtrip():
    """QE angle1 (from z) and angle2 (from x in the xy plane), in degrees."""
    m = MagneticMoment.from_polar(2.0, 90.0, 30.0)
    assert m.vector == pytest.approx((2.0 * math.cos(math.radians(30)), 1.0, 0.0))
    angle1, angle2 = m.angles()
    assert angle1 == pytest.approx(90.0)
    assert angle2 == pytest.approx(30.0)

    along_z = MagneticMoment.from_polar(1.5, 0.0, 123.0)
    assert along_z.vector == pytest.approx((0.0, 0.0, 1.5))
    assert along_z.angles() == pytest.approx((0.0, 0.0))

    assert MagneticMoment.from_polar(1.0, 180.0, 0.0).angles()[0] == pytest.approx(
        180.0
    )
    assert MagneticMoment.zero().angles() == (0.0, 0.0)


def test_moment_equivalence_uses_a_tolerance():
    a = MagneticMoment((0.0, 0.0, 3.0))
    b = MagneticMoment((0.0, 0.0004, 3.0002))
    assert a.is_equivalent(b)
    assert not a.is_equivalent(MagneticMoment((0.0, 0.0, -3.0)))
    assert not a.is_equivalent(b, tol=1e-6)


def test_zero_moment_has_no_direction():
    with pytest.raises(StructureError):
        MagneticMoment.zero().direction()


@pytest.mark.parametrize("bad", [(1.0, 2.0), (1.0, "x", 3.0), (float("nan"), 0, 0)])
def test_moment_rejects_bad_vectors(bad):
    with pytest.raises(StructureError):
        MagneticMoment(bad)


# --------------------------------------------------------------------------
# AtomicSite
# --------------------------------------------------------------------------


def test_site_defaults_and_wrapping():
    site = AtomicSite("Fe", (1.25, -0.5, 0.0))
    assert site.label == "Fe"
    assert site.occupancy == 1.0
    assert not site.is_magnetic
    assert site.wrapped().frac_coords == pytest.approx((0.25, 0.5, 0.0))
    assert AtomicSite("Fe", (0, 0, 0), MagneticMoment((0, 0, 2.2))).is_magnetic
    assert not AtomicSite("Fe", (0, 0, 0), MagneticMoment.zero()).is_magnetic


def test_site_rejects_unknown_element_and_bad_occupancy():
    with pytest.raises(StructureError):
        AtomicSite("Xx", (0, 0, 0))
    with pytest.raises(StructureError):
        AtomicSite("Fe", (0, 0, 0), occupancy=0.0)
    with pytest.raises(StructureError):
        AtomicSite("Fe", (0, 0, 0), occupancy=1.5)
    with pytest.raises(StructureError):
        AtomicSite("Fe", (0, 0, 0), label="  ")


# --------------------------------------------------------------------------
# NormalizedStructure
# --------------------------------------------------------------------------


def test_structure_volume_and_cartesian_coordinates():
    s = fe_structure()
    assert s.volume == pytest.approx(ALAT**3 / 2.0)
    assert s.lattice_matrix == pytest.approx(BCC)
    assert s.num_sites == 1
    assert s.elements == ("Fe",)
    assert s.species_labels == ("Fe",)
    assert s.frac_to_cart((1.0, 0.0, 0.0)) == pytest.approx(BCC[0])

    two = NormalizedStructure(
        np.eye(3) * 4.0,
        (AtomicSite("Na", (0, 0, 0)), AtomicSite("Cl", (0.5, 0.5, 0.5))),
    )
    assert two.cart_coords()[1] == pytest.approx([2.0, 2.0, 2.0])
    assert two.elements == ("Na", "Cl")


def test_reciprocal_lattice_includes_two_pi():
    s = fe_structure()
    product = s.lattice_matrix @ s.reciprocal_lattice().T
    assert product == pytest.approx(2.0 * np.pi * np.eye(3))


def test_structure_is_immutable_and_comparable():
    a = fe_structure()
    b = fe_structure()
    assert a == b
    with pytest.raises(Exception):
        a.sites = ()


def test_structure_magnetism_flag():
    assert not fe_structure().is_magnetic
    assert fe_structure(MagneticMoment((0, 0, 2.2))).is_magnetic


def test_structure_rejects_degenerate_lattice_and_bad_sites():
    with pytest.raises(StructureError):
        NormalizedStructure(np.zeros((3, 3)), (AtomicSite("Fe", (0, 0, 0)),))
    with pytest.raises(StructureError):
        NormalizedStructure([[1, 0, 0], [0, 1, 0]], (AtomicSite("Fe", (0, 0, 0)),))
    with pytest.raises(StructureError):
        NormalizedStructure(np.eye(3), ())
    with pytest.raises(StructureError):
        NormalizedStructure(
            np.eye(3),
            (
                AtomicSite("Fe", (0, 0, 0), label="M"),
                AtomicSite("Co", (0.5, 0.5, 0.5), label="M"),
            ),
        )


def test_with_sites_keeps_the_lattice():
    s = fe_structure().with_sites(
        [AtomicSite("Fe", (0, 0, 0), MagneticMoment((0, 0, 2)))]
    )
    assert s.lattice_matrix == pytest.approx(BCC)
    assert s.is_magnetic


# --------------------------------------------------------------------------
# QE model
# --------------------------------------------------------------------------


def fe_pw_input(calculation="scf", kpoints=None, **system_extra):
    system = {"ibrav": 0, "A": ALAT, "nat": 1, "ntyp": 1, "ecutwfc": 64.0}
    system.update(system_extra)
    return PwInput(
        control=Namelist("control", {"calculation": calculation, "prefix": "pwscf"}),
        system=Namelist("system", system),
        electrons=Namelist("electrons", {"conv_thr": RawValue("1.0d-8")}),
        species=(Species("Fe", 55.845, "Fe.pbe-spn-rrkjus_psl.0.2.1.UPF"),),
        positions=(AtomicPosition("Fe", (0, 0, 0)),),
        kpoints=kpoints or KPointsAutomatic((21, 21, 21)),
        cell=CellParameters("alat", BCC / ALAT),
    )


def test_namelist_validation():
    n = Namelist("system", {"nbnd": 31, "nosym": True, "conv_thr": RawValue("1.d-10")})
    assert n.get("nbnd") == 31 and n.get("missing") is None
    with pytest.raises(InputModelError):
        Namelist("", {})
    with pytest.raises(InputModelError):
        Namelist("system", {"": 1})
    with pytest.raises(InputModelError):
        Namelist("system", {"x": [1, 2]})
    with pytest.raises(InputModelError):
        Namelist("system", {"x": float("inf")})
    with pytest.raises(InputModelError):
        RawValue(" ")


def test_pw_input_accepts_a_consistent_scf():
    pw = fe_pw_input()
    assert pw.calculation == "scf"
    assert pw.nat == 1 and pw.ntyp == 1
    assert pw.cell.matrix() == pytest.approx(BCC / ALAT)


def test_pw_input_checks_counts_and_labels():
    with pytest.raises(InputModelError, match="nat"):
        fe_pw_input(nat=2)
    with pytest.raises(InputModelError, match="ntyp"):
        fe_pw_input(ntyp=2)
    with pytest.raises(InputModelError, match="ATOMIC_SPECIES"):
        PwInput(
            control=Namelist("control"),
            system=Namelist("system"),
            electrons=Namelist("electrons"),
            species=(Species("Fe", 55.845, "Fe.UPF"),),
            positions=(AtomicPosition("Co", (0, 0, 0)),),
            kpoints=KPointsAutomatic((4, 4, 4)),
        )
    with pytest.raises(InputModelError, match="CELL_PARAMETERS"):
        PwInput(
            control=Namelist("control"),
            system=Namelist("system", {"ibrav": 0}),
            electrons=Namelist("electrons"),
            species=(Species("Fe", 55.845, "Fe.UPF"),),
            positions=(AtomicPosition("Fe", (0, 0, 0)),),
            kpoints=KPointsAutomatic((4, 4, 4)),
        )
    with pytest.raises(InputModelError, match="&control"):
        PwInput(
            control=Namelist("system"),
            system=Namelist("system"),
            electrons=Namelist("electrons"),
            species=(Species("Fe", 55.845, "Fe.UPF"),),
            positions=(AtomicPosition("Fe", (0, 0, 0)),),
            kpoints=KPointsAutomatic((4, 4, 4)),
        )


def test_species_and_position_validation():
    with pytest.raises(InputModelError):
        Species("Fe", 0.0, "Fe.UPF")
    with pytest.raises(InputModelError):
        Species("Fe", 55.8, "")
    with pytest.raises(InputModelError):
        AtomicPosition("Fe", (0, 0))
    with pytest.raises(InputModelError):
        CellParameters("alat", np.zeros((3, 3)))
    with pytest.raises(InputModelError):
        CellParameters("nm", np.eye(3))


def test_kpoints_models():
    auto = KPointsAutomatic((8, 8, 8), (1, 1, 1))
    assert auto.num_points == 512
    with pytest.raises(InputModelError):
        KPointsAutomatic((0, 8, 8))
    with pytest.raises(InputModelError):
        KPointsAutomatic((8, 8, 8), (2, 0, 0))

    explicit = KPointsList("crystal", ((0, 0, 0, 0.5), (0.5, 0, 0, 0.5)))
    assert explicit.num_points == 2
    with pytest.raises(InputModelError):
        KPointsList("crystal", ((0, 0, 0),))
    with pytest.raises(InputModelError):
        KPointsList("crystal", ((0, 0, 0, -1.0),))

    path = KPointsPath(
        "crystal_b",
        (
            KPathPoint((0, 0, 0), 20, "G"),
            KPathPoint((0.5, -0.5, 0.5), 0, "H"),
            KPathPoint((0.25, 0.25, 0.25), 10, "P"),
            KPathPoint((0, 0, 0.5), 10, "N"),
        ),
    )
    assert path.num_points == 4
    assert [(a.label, b.label) for a, b in path.segments()] == [("G", "H"), ("P", "N")]
    with pytest.raises(InputModelError):
        KPointsPath("crystal_b", (KPathPoint((0, 0, 0), 1),))
    with pytest.raises(InputModelError):
        KPathPoint((0, 0, 0), -1)


def test_namelist_input():
    inp = NamelistInput(
        (Namelist("inputpp", {"plot_num": 7}), Namelist("plot", {"iflag": 3}))
    )
    assert inp.namelist("plot").get("iflag") == 3
    with pytest.raises(KeyError):
        inp.namelist("missing")
    with pytest.raises(InputModelError):
        NamelistInput(())
    with pytest.raises(InputModelError):
        NamelistInput((Namelist("a"), Namelist("a")))


# --------------------------------------------------------------------------
# Wannier90 model
# --------------------------------------------------------------------------


def mesh_kpoints(n):
    return tuple(
        (i / n, j / n, k / n) for i in range(n) for j in range(n) for k in range(n)
    )


def fe_win(spinors=False, num_wann=None, **extra):
    factor = 2 if spinors else 1
    kwargs = dict(
        num_wann=9 * factor if num_wann is None else num_wann,
        num_bands=27 * factor,
        unit_cell_cart=BCC,
        atoms_frac=(AtomFrac("Fe", (0, 0, 0)),),
        mp_grid=(2, 2, 2),
        kpoints=mesh_kpoints(2),
        projections=(Projection("Fe", ("s", "p", "d")),),
        spinors=spinors,
        exclude_bands=(1, 4 * factor),
        parameters={"dis_num_iter": 200, "write_hr": True},
    )
    kwargs.update(extra)
    return Wannier90Input(**kwargs)


def test_projection_sizes():
    assert Projection("Fe", ("s", "p", "d")).num_wann_per_site == 9
    assert Projection("C", ("sp3",)).num_wann_per_site == 4
    with pytest.raises(InputModelError):
        Projection("Fe", ("g",))
    with pytest.raises(InputModelError):
        Projection("Fe", ())


def test_win_counts_are_consistent():
    win = fe_win()
    assert win.projected_num_wann() == 9
    assert win.num_exclude == 4
    spinor = fe_win(spinors=True)
    assert spinor.projected_num_wann() == 18
    assert spinor.num_exclude == 8


def test_win_projection_count_must_match_num_wann():
    with pytest.raises(InputModelError, match="Wannier functions"):
        fe_win(num_wann=10, num_bands=30)


def test_win_unresolvable_projection_site_is_not_checked():
    win = fe_win(projections=(Projection("f=0,0,0", ("s", "p", "d")),))
    assert win.projected_num_wann() is None


def test_win_validation():
    with pytest.raises(InputModelError, match="mp_grid"):
        fe_win(kpoints=mesh_kpoints(3))
    with pytest.raises(InputModelError, match="num_bands"):
        fe_win(num_bands=5)
    with pytest.raises(InputModelError, match="even"):
        fe_win(spinors=True, num_wann=9, num_bands=27)
    with pytest.raises(InputModelError, match="exclude_bands"):
        fe_win(exclude_bands=(4, 1))
    with pytest.raises(InputModelError, match="field"):
        fe_win(parameters={"num_wann": 9})
    with pytest.raises(InputModelError):
        fe_win(unit_cell_cart=np.zeros((3, 3)))
    with pytest.raises(InputModelError):
        KPathSegment("", (0, 0, 0), "H", (0.5, -0.5, 0.5))


# --------------------------------------------------------------------------
# CalculationPlan
# --------------------------------------------------------------------------


def fe_plan(**override):
    n = 2
    explicit = KPointsList("crystal", tuple(k + (1.0 / n**3,) for k in mesh_kpoints(n)))
    path = KPointsPath(
        "crystal_b",
        (KPathPoint((0, 0, 0), 20, "G"), KPathPoint((0.5, -0.5, 0.5), 10, "H")),
    )
    inputs = dict(
        structure=fe_structure(),
        scf=fe_pw_input(),
        nscf=fe_pw_input("nscf", explicit, nbnd=31, nosym=True),
        check_wannier=fe_pw_input(
            "nscf", KPointsAutomatic((n, n, n), (1, 1, 1)), nbnd=17
        ),
        bands_nscf=fe_pw_input("bands", path, nbnd=17),
        bands=NamelistInput((Namelist("bands", {"prefix": "pwscf"}),)),
        pw2wan=NamelistInput((Namelist("inputpp", {"seedname": "pwscf"}),)),
        projwfc=NamelistInput((Namelist("projwfc", {"prefix": "pwscf"}),)),
        pp=NamelistInput(
            (Namelist("inputpp", {"plot_num": 7}), Namelist("plot", {"iflag": 3}))
        ),
        wannier90=fe_win(),
    )
    inputs.update(override)
    return CalculationPlan(**inputs)


def test_plan_layout_matches_the_generated_files():
    plan = fe_plan()
    assert list(plan.files()) == [
        "scf.in",
        "nscf.in",
        "pw2wan.in",
        "pwscf.win",
        "check_wannier/nscf.in",
        "band/nscf.in",
        "band/band.in",
        "band/proj.in",
        "band/pp.in",
    ]
    assert plan.files()["pwscf.win"] is plan.wannier90


def test_plan_checks_calculation_types_and_meshes():
    with pytest.raises(InputModelError, match="calculation = 'scf'"):
        fe_plan(scf=fe_pw_input("nscf"))
    with pytest.raises(InputModelError, match="mp_grid"):
        fe_plan(
            check_wannier=fe_pw_input("nscf", KPointsAutomatic((4, 4, 4), (1, 1, 1)))
        )
    with pytest.raises(InputModelError, match="same k points"):
        fe_plan(nscf=fe_pw_input("nscf", KPointsAutomatic((2, 2, 2))))
    with pytest.raises(InputModelError, match="K_POINTS path"):
        fe_plan(bands_nscf=fe_pw_input("bands"))
