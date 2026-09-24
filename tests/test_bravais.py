"""Tests of the ``ibrav /= 0`` representation (:mod:`cif2qewan.qe.bravais`).

Every ibrav value that cif2qewan writes is exercised with a structure built
by pymatgen from a space group, scrambled by a rotation and a change of
primitive basis, and re-expressed in QE's lattice; the result must be the
same crystal (pymatgen's StructureMatcher, which allows a rigid rotation)
and its lattice must be exactly QE's ``latgen`` lattice for the parameters
written into the input. None of this needs QE itself.
"""

import numpy as np
import pytest

from cif2qewan.exceptions import StructureError
from cif2qewan.qe.bravais import (
    IBRAV_OF_LATTICE,
    PARAMETERS_OF_IBRAV,
    BravaisLattice,
    bravais_structure,
    crystal_system,
    ibrav_of,
    qe_lattice,
)
from cif2qewan.structure.model import AtomicSite, MagneticMoment, NormalizedStructure
from cif2qewan.structure.readers import PymatgenReader
from conftest import REPO_ROOT

pytest.importorskip("pymatgen")
pytest.importorskip("spglib")

from pymatgen.analysis.structure_matcher import StructureMatcher  # noqa: E402
from pymatgen.core import Lattice, Structure  # noqa: E402
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # noqa: E402

MN3SN = REPO_ROOT / "tests" / "fixtures" / "Mn3Sn.mcif"

# (expected ibrav, space group, lattice, species, coordinates of the
# inequivalent sites); the space group is the one spglib should find.
CASES = [
    (1, "Pm-3m", Lattice.cubic(3.3), ["Po"], [[0, 0, 0]]),
    (2, "Fm-3m", Lattice.cubic(3.6), ["Cu"], [[0, 0, 0]]),
    (-3, "Im-3m", Lattice.cubic(2.87), ["Fe"], [[0, 0, 0]]),
    (4, "P6_3/mmc", Lattice.hexagonal(3.2, 5.2), ["Mg"], [[1 / 3, 2 / 3, 1 / 4]]),
    (
        4,
        "P-3m1",
        Lattice.hexagonal(4.24, 6.84),
        ["Cd", "I"],
        [[0, 0, 0], [1 / 3, 2 / 3, 0.25]],
    ),
    (5, "R-3m", Lattice.hexagonal(3.76, 10.55), ["As"], [[0, 0, 0.227]]),
    (6, "P4/mmm", Lattice.tetragonal(3.0, 4.2), ["Po"], [[0, 0, 0]]),
    (7, "I4/mmm", Lattice.tetragonal(3.25, 4.95), ["In"], [[0, 0, 0]]),
    (8, "Pmmm", Lattice.orthorhombic(3.0, 4.0, 5.0), ["Po"], [[0, 0, 0]]),
    (9, "Cmcm", Lattice.orthorhombic(2.85, 5.87, 4.96), ["U"], [[0, 0.1025, 0.25]]),
    (
        91,
        "Amm2",
        Lattice.orthorhombic(3.0, 4.0, 5.0),
        ["Na", "Cl"],
        [[0, 0, 0], [0.5, 0, 0.3]],
    ),
    (10, "Fmmm", Lattice.orthorhombic(3.0, 4.0, 5.0), ["Po"], [[0, 0, 0]]),
    (11, "Immm", Lattice.orthorhombic(3.0, 4.0, 5.0), ["Po"], [[0, 0, 0]]),
    (
        -12,
        "P2_1/m",
        Lattice.monoclinic(4.0, 5.0, 6.0, 100.0),
        ["S", "Se"],
        [[0.1, 0.25, 0.3], [0.6, 0.25, 0.85]],
    ),
    (-13, "C2/m", Lattice.monoclinic(5.0, 3.0, 6.0, 105.0), ["S"], [[0.2, 0, 0.3]]),
    (
        14,
        "P-1",
        Lattice.from_parameters(4.0, 5.0, 6.0, 80.0, 85.0, 95.0),
        ["S"],
        [[0.1, 0.2, 0.3]],
    ),
]


def to_pymatgen(structure: NormalizedStructure) -> Structure:
    return Structure(
        structure.lattice_matrix,
        [site.element for site in structure.sites],
        structure.frac_coords(),
    )


def scrambled_primitive(spacegroup, lattice, species, coords, seed=0):
    """A primitive cell of the structure in a random orientation and basis."""
    conventional = Structure.from_spacegroup(spacegroup, lattice, species, coords)
    primitive = SpacegroupAnalyzer(conventional, symprec=1e-3).find_primitive()
    rng = np.random.default_rng(seed)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(rotation) < 0:
        rotation[:, 0] *= -1
    basis = np.array([[1, 1, 0], [0, 1, 0], [1, 0, 1]])  # det +1
    new_lattice = basis @ primitive.lattice.matrix @ rotation
    frac = (primitive.frac_coords @ np.linalg.inv(basis)) % 1.0
    sites = [
        AtomicSite(str(site.specie.symbol), tuple(f))
        for site, f in zip(primitive, frac)
    ]
    return NormalizedStructure(new_lattice, sites)


def assert_same_crystal(old: NormalizedStructure, new: NormalizedStructure):
    matcher = StructureMatcher(
        ltol=1e-4,
        stol=1e-4,
        angle_tol=0.01,
        primitive_cell=False,
        scale=False,
        attempt_supercell=False,
    )
    assert matcher.fit(to_pymatgen(old), to_pymatgen(new))
    assert new.volume == pytest.approx(old.volume, rel=1e-8)
    assert [s.label for s in new.sites] == [s.label for s in old.sites]


@pytest.mark.parametrize(
    ("ibrav", "spacegroup", "lattice", "species", "coords"),
    CASES,
    ids=[f"ibrav{c[0]}_{c[1]}" for c in CASES],
)
def test_every_ibrav_reproduces_the_crystal(
    ibrav, spacegroup, lattice, species, coords
):
    structure = scrambled_primitive(spacegroup, lattice, species, coords)
    new, bravais = bravais_structure(structure)
    assert bravais.ibrav == ibrav
    assert bravais.symbol.replace("_", "") == spacegroup.replace("_", "")
    assert tuple(bravais.parameters) == PARAMETERS_OF_IBRAV[ibrav]
    assert new.lattice_matrix == pytest.approx(qe_lattice(ibrav, bravais.parameters))
    assert np.linalg.det(new.lattice_matrix) > 0
    assert np.all(new.frac_coords() >= 0) and np.all(new.frac_coords() < 1)
    assert_same_crystal(structure, new)


@pytest.mark.parametrize("shear", [4, 17])
def test_primitive_cell_with_large_integer_shear(shear):
    """An arbitrary primitive basis must not be limited to small coefficients."""
    old = NormalizedStructure(
        np.array([[3.0, 0, 0], [3.0 * shear, 3.0, 0], [0, 0, 3.0]]),
        (AtomicSite("Po", (0.123, 0.234, 0.345)),),
    )
    new, bravais = bravais_structure(old)
    assert bravais.ibrav == 1
    assert_same_crystal(old, new)


def test_qe_lattice_definitions():
    """Lengths, angles and centering vectors of QE's latgen lattices."""
    p = {"A": 3.0, "B": 4.0, "C": 5.0, "cosBC": 0.1, "cosAC": -0.2, "cosAB": 0.3}
    for ibrav, names in PARAMETERS_OF_IBRAV.items():
        params = {k: p[k] for k in names}
        if ibrav == 5:
            params["cosAB"] = 0.4
        v = qe_lattice(ibrav, params)
        assert np.linalg.det(v) > 0, ibrav
    a1, a2, a3 = qe_lattice(4, {"A": 3.0, "C": 5.0})
    assert np.linalg.norm(a1) == np.linalg.norm(a2) == 3.0
    assert np.dot(a1, a2) / 9.0 == pytest.approx(-0.5)  # 120 degrees
    assert list(a3) == [0, 0, 5.0]
    v = qe_lattice(5, {"A": 3.0, "cosAB": 0.4})
    assert np.linalg.norm(v, axis=1) == pytest.approx([3.0, 3.0, 3.0])
    assert (v @ v.T)[np.triu_indices(3, 1)] / 9.0 == pytest.approx([0.4] * 3)
    a1, a2, a3 = qe_lattice(9, {"A": 3.0, "B": 4.0, "C": 5.0})
    assert list(a1 - a2) == [3.0, 0.0, 0.0] and list(a1 + a2) == [0.0, 4.0, 0.0]
    a1, a2, a3 = qe_lattice(91, {"A": 3.0, "B": 4.0, "C": 5.0})
    assert list(a2 + a3) == [0.0, 4.0, 0.0] and list(a3 - a2) == [0.0, 0.0, 5.0]
    a1, a2, a3 = qe_lattice(-13, {"A": 3.0, "B": 4.0, "C": 5.0, "cosAC": -0.2})
    assert list(a1 + a2) == [0.0, 4.0, 0.0]
    assert np.dot(a1 - a2, a3) / (3.0 * 5.0) == pytest.approx(-0.2)
    v = qe_lattice(
        14, {"A": 3.0, "B": 4.0, "C": 5.0, "cosBC": 0.1, "cosAC": -0.2, "cosAB": 0.3}
    )
    g = v @ v.T
    assert np.sqrt(np.diag(g)) == pytest.approx([3.0, 4.0, 5.0])
    assert g[1, 2] / 20 == pytest.approx(0.1)
    assert g[0, 2] / 15 == pytest.approx(-0.2)
    assert g[0, 1] / 12 == pytest.approx(0.3)


def test_ibrav_table_covers_the_230_space_groups():
    """The (system, centering) table agrees with the space-group list of mcif2qewan."""
    from pymatgen.symmetry.groups import SpaceGroup

    reference = {
        1: [195, 198, 200, 201, 205, 207, 208, 212, 213, 215, 218, 221, 222, 223, 224],
        2: [196, 202, 203, 209, 210, 216, 219, 225, 226, 227, 228],
        -3: [197, 199, 204, 206, 211, 214, 217, 220, 229, 230],
        4: [143, 144, 145, 147]
        + list(range(149, 155))
        + list(range(156, 160))
        + list(range(162, 166))
        + list(range(168, 195)),
        5: [146, 148, 155, 160, 161, 166, 167],
        6: list(range(75, 79))
        + [81]
        + list(range(83, 87))
        + list(range(89, 97))
        + list(range(99, 107))
        + list(range(111, 119))
        + list(range(123, 139)),
        7: [
            79,
            80,
            82,
            87,
            88,
            97,
            98,
            107,
            108,
            109,
            110,
            119,
            120,
            121,
            122,
            139,
            140,
            141,
            142,
        ],
        8: [16, 17, 18, 19] + list(range(25, 35)) + list(range(47, 63)),
        9: [20, 21, 35, 36, 37, 63, 64, 65, 66, 67, 68],
        91: [38, 39, 40, 41],
        10: [22, 42, 43, 69, 70],
        11: [23, 24, 44, 45, 46, 71, 72, 73, 74],
        -12: [3, 4, 6, 7, 10, 11, 13, 14],
        -13: [5, 8, 9, 12, 15],
        14: [1, 2],
    }
    expected = {n: ibrav for ibrav, numbers in reference.items() for n in numbers}
    assert sorted(expected) == list(range(1, 231))
    for number in range(1, 231):
        symbol = SpaceGroup.from_int_number(number).symbol
        assert ibrav_of(number, symbol) == expected[number], (number, symbol)
    assert crystal_system(1) == "triclinic" and crystal_system(230) == "cubic"
    assert set(IBRAV_OF_LATTICE.values()) == set(PARAMETERS_OF_IBRAV)


def test_bravais_lattice_validation():
    with pytest.raises(StructureError, match="not supported"):
        BravaisLattice(3, {"A": 1.0}, 229, "Im-3m")
    with pytest.raises(StructureError, match="takes"):
        BravaisLattice(4, {"A": 1.0}, 194, "P6_3/mmc")
    with pytest.raises(StructureError, match="not supported"):
        qe_lattice(12, {"A": 1.0, "B": 1.0, "C": 1.0, "cosAB": 0.0})
    with pytest.raises(StructureError, match="right-handed"):
        bravais_structure(
            NormalizedStructure(-np.eye(3) * 3.0, (AtomicSite("Po", (0, 0, 0)),))
        )


def test_moments_rotate_with_the_crystal():
    """Mn3Sn (noncollinear, hexagonal): ibrav 4 and moments carried along."""
    import cif2qewan.qe.bravais as bravais_module

    structure = PymatgenReader().read(MN3SN)
    new, bravais = bravais_structure(structure)
    assert bravais.ibrav == 4 and bravais.space_group == 194
    assert_same_crystal(structure, new)
    old_m = np.array([s.magmom.as_array() for s in structure.sites])
    new_m = np.array([s.magmom.as_array() for s in new.sites])
    # a rotation preserves magnitudes and mutual angles ...
    assert np.linalg.norm(new_m, axis=1) == pytest.approx(np.linalg.norm(old_m, axis=1))
    assert new_m @ new_m.T == pytest.approx(old_m @ old_m.T, abs=1e-8)
    # ... and the moments follow the same rotation as the atoms: the in-plane
    # kagome moments stay perpendicular to the hexagonal axis, which is z in QE
    assert new_m[:, 2] == pytest.approx(0.0, abs=1e-6)
    matrix, rotation = bravais_module._basis_change(
        structure.lattice_matrix, new.lattice_matrix
    )
    assert new_m == pytest.approx(old_m @ rotation, abs=1e-8)
    frac = structure.frac_coords() @ np.linalg.inv(matrix)
    delta = new.frac_coords() - frac
    assert delta == pytest.approx(np.rint(delta), abs=1e-8)


def test_antiferromagnet_keeps_its_magnetic_cell():
    """Two Fe with opposite moments in a cubic cell: simple cubic, not bcc."""
    a = 2.87
    up, down = MagneticMoment((0, 0, 2.2)), MagneticMoment((0, 0, -2.2))
    afm = NormalizedStructure(
        np.eye(3) * a,
        (AtomicSite("Fe", (0, 0, 0), up), AtomicSite("Fe", (0.5, 0.5, 0.5), down)),
    )
    new, bravais = bravais_structure(afm)
    assert bravais.ibrav == 1 and bravais.space_group == 221
    assert new.num_sites == 2
    assert [s.magmom.vector[2] for s in new.sites] == pytest.approx([2.2, -2.2])
    # without the moments the same cell is a conventional bcc cell: not primitive
    plain = afm.with_sites(AtomicSite(s.element, s.frac_coords) for s in afm.sites)
    with pytest.raises(StructureError, match="not a primitive cell"):
        bravais_structure(plain)


def test_non_primitive_cell_is_rejected():
    conventional = NormalizedStructure(
        np.eye(3) * 2.87,
        (AtomicSite("Fe", (0, 0, 0)), AtomicSite("Fe", (0.5, 0.5, 0.5))),
    )
    with pytest.raises(StructureError, match="not a primitive cell"):
        bravais_structure(conventional)


def test_lattice_mismatch_is_reported(monkeypatch):
    """If the QE lattice cannot reproduce the structure's lattice, fail clearly."""
    import cif2qewan.qe.bravais as bravais

    structure = scrambled_primitive(*CASES[2][1:])
    monkeypatch.setattr(
        bravais, "_symmetry", lambda s, symprec: (229, "Im-3m", np.eye(3) * 2.0)
    )
    with pytest.raises(StructureError, match="does not match"):
        bravais_structure(structure)
