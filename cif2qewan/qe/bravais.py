"""The ``ibrav /= 0`` representation of a structure for Quantum ESPRESSO.

By default cif2qewan writes ``ibrav = 0`` and the lattice vectors of the
structure in ``CELL_PARAMETERS``. With ``use_ibrav = true`` in the
configuration the cell is written the way ``pw.x`` classifies it instead:
the Bravais-lattice index ``ibrav`` and the lattice parameters ``A``,
``B``, ``C``, ``cosAB``, ``cosAC``, ``cosBC`` (angstrom and cosines), and no
``CELL_PARAMETERS`` card. ``pw.x`` then builds the lattice vectors itself
from these numbers, so the structure has to be re-expressed in exactly
those vectors: :func:`bravais_structure` returns the structure with new
fractional coordinates, rigidly rotated (moments included) into QE's
orientation of the lattice, together with the :class:`BravaisLattice`.

The lattice type is found with spglib from the space group of the
structure: the number gives the crystal system, the first letter of the
Hermann-Mauguin symbol in the standard setting the centering, and spglib's
standardized conventional cell the lattice parameters. Sites are
distinguished by species label (what ``pw.x`` sees as species); the
magnetic moments are ignored as long as the cell is a primitive cell for
that description, because moments that only break point symmetry leave the
Bravais lattice unchanged (Mn3Sn stays hexagonal). When the cell is a
supercell of the chemical primitive cell (an antiferromagnet), sites are
distinguished by moment as well, so that the magnetic cell keeps its own
lattice type instead of being folded. The vectors follow ``Modules/latgen.f90`` of Quantum
ESPRESSO 7.x. ``ibrav = -13`` changed in QE 6.4.1; the inputs generated
here assume the current definition (as does the ``BEWARE`` note pw.x
prints). ``-3``, ``-12`` and ``-13`` are used rather than ``3``, ``12`` and
``13`` because they match the standard settings (unique axis b for the
monoclinic groups).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np

from cif2qewan.exceptions import StructureError
from cif2qewan.structure.magnetism import classify_magnetic_sites
from cif2qewan.structure.model import AtomicSite, MagneticMoment, NormalizedStructure

logger = logging.getLogger(__name__)

#: Relative deviation allowed between the metric of the structure's lattice
#: and the metric of the QE lattice built from the standardized parameters.
LATTICE_TOLERANCE = 1.0e-3

_CRYSTAL_SYSTEMS = (
    (2, "triclinic"),
    (15, "monoclinic"),
    (74, "orthorhombic"),
    (142, "tetragonal"),
    (167, "trigonal"),
    (194, "hexagonal"),
    (230, "cubic"),
)

#: (crystal system, centering letter of the standard setting) -> ibrav.
IBRAV_OF_LATTICE = {
    ("cubic", "P"): 1,
    ("cubic", "F"): 2,
    ("cubic", "I"): -3,
    ("hexagonal", "P"): 4,
    ("trigonal", "P"): 4,
    ("trigonal", "R"): 5,
    ("tetragonal", "P"): 6,
    ("tetragonal", "I"): 7,
    ("orthorhombic", "P"): 8,
    ("orthorhombic", "C"): 9,
    ("orthorhombic", "A"): 91,
    ("orthorhombic", "F"): 10,
    ("orthorhombic", "I"): 11,
    ("monoclinic", "P"): -12,
    ("monoclinic", "C"): -13,
    ("triclinic", "P"): 14,
}

#: The &system parameters each ibrav takes, in the order they are written.
PARAMETERS_OF_IBRAV = {
    1: ("A",),
    2: ("A",),
    -3: ("A",),
    4: ("A", "C"),
    5: ("A", "cosAB"),
    6: ("A", "C"),
    7: ("A", "C"),
    8: ("A", "B", "C"),
    9: ("A", "B", "C"),
    91: ("A", "B", "C"),
    10: ("A", "B", "C"),
    11: ("A", "B", "C"),
    -12: ("A", "B", "C", "cosAC"),
    -13: ("A", "B", "C", "cosAC"),
    14: ("A", "B", "C", "cosBC", "cosAC", "cosAB"),
}


def crystal_system(space_group: int) -> str:
    """Crystal system of a space-group number (1-230)."""
    if not 1 <= int(space_group) <= 230:
        raise StructureError(f"invalid space-group number {space_group}")
    for last, system in _CRYSTAL_SYSTEMS:
        if space_group <= last:
            return system
    raise AssertionError("unreachable")


def ibrav_of(space_group: int, symbol: str) -> int:
    """``ibrav`` for a space group given by number and standard HM symbol."""
    letter = symbol.strip()[:1]
    try:
        return IBRAV_OF_LATTICE[(crystal_system(space_group), letter)]
    except KeyError:
        raise StructureError(
            f"no ibrav for space group {space_group} ({symbol}): "
            f"unexpected centering {letter!r}"
        ) from None


def qe_lattice(ibrav: int, parameters: Dict[str, float]) -> np.ndarray:
    """QE's lattice vectors (rows, angstrom) for ``ibrav`` and ``A``, ``B``, ...

    Transcribed from ``latgen_lib`` in ``Modules/latgen.f90`` (QE 7.x); the
    angles enter as the cosines ``cosAB`` (gamma), ``cosAC`` (beta) and
    ``cosBC`` (alpha) exactly as ``pw.x`` reads them.
    """
    a = float(parameters["A"])
    b = float(parameters.get("B", 0.0))
    c = float(parameters.get("C", 0.0))
    cosab = float(parameters.get("cosAB", 0.0))
    cosac = float(parameters.get("cosAC", 0.0))
    cosbc = float(parameters.get("cosBC", 0.0))
    if ibrav == 1:
        rows = [[a, 0, 0], [0, a, 0], [0, 0, a]]
    elif ibrav == 2:
        rows = [[-a / 2, 0, a / 2], [0, a / 2, a / 2], [-a / 2, a / 2, 0]]
    elif ibrav == -3:
        rows = [[-a / 2, a / 2, a / 2], [a / 2, -a / 2, a / 2], [a / 2, a / 2, -a / 2]]
    elif ibrav == 4:
        rows = [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, c]]
    elif ibrav == 5:
        tx = math.sqrt((1 - cosab) / 2)
        ty = math.sqrt((1 - cosab) / 6)
        tz = math.sqrt((1 + 2 * cosab) / 3)
        rows = [
            [a * tx, -a * ty, a * tz],
            [0, 2 * a * ty, a * tz],
            [-a * tx, -a * ty, a * tz],
        ]
    elif ibrav == 6:
        rows = [[a, 0, 0], [0, a, 0], [0, 0, c]]
    elif ibrav == 7:
        rows = [[a / 2, -a / 2, c / 2], [a / 2, a / 2, c / 2], [-a / 2, -a / 2, c / 2]]
    elif ibrav == 8:
        rows = [[a, 0, 0], [0, b, 0], [0, 0, c]]
    elif ibrav == 9:
        rows = [[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [0, 0, c]]
    elif ibrav == 91:
        rows = [[a, 0, 0], [0, b / 2, -c / 2], [0, b / 2, c / 2]]
    elif ibrav == 10:
        rows = [[a / 2, 0, c / 2], [a / 2, b / 2, 0], [0, b / 2, c / 2]]
    elif ibrav == 11:
        rows = [[a / 2, b / 2, c / 2], [-a / 2, b / 2, c / 2], [-a / 2, -b / 2, c / 2]]
    elif ibrav == -12:
        sinac = math.sqrt(1 - cosac**2)
        rows = [[a, 0, 0], [0, b, 0], [c * cosac, 0, c * sinac]]
    elif ibrav == -13:
        sinac = math.sqrt(1 - cosac**2)
        rows = [[a / 2, b / 2, 0], [-a / 2, b / 2, 0], [c * cosac, 0, c * sinac]]
    elif ibrav == 14:
        sinab = math.sqrt(1 - cosab**2)
        term = 1 + 2 * cosbc * cosac * cosab - cosbc**2 - cosac**2 - cosab**2
        if term < 0:
            raise StructureError("triclinic cell angles are inconsistent")
        rows = [
            [a, 0, 0],
            [b * cosab, b * sinab, 0],
            [
                c * cosac,
                c * (cosbc - cosac * cosab) / sinab,
                c * math.sqrt(term / (1 - cosab**2)),
            ],
        ]
    else:
        raise StructureError(f"ibrav = {ibrav} is not supported")
    return np.array(rows, dtype=float)


@dataclass(frozen=True)
class BravaisLattice:
    """QE's description of a Bravais lattice: ``ibrav`` and its parameters.

    Attributes
    ----------
    ibrav : int
        The QE lattice index.
    parameters : dict
        ``A``, ``B``, ``C`` in angstrom and ``cosAB``, ``cosAC``, ``cosBC``
        as needed by ``ibrav``, in the order ``pw.x`` documents them.
    space_group : int
        Space-group number of the structure (sites distinguished by species
        and moment).
    symbol : str
        Hermann-Mauguin symbol of the standard setting.
    """

    ibrav: int
    parameters: Dict[str, float]
    space_group: int
    symbol: str

    def __post_init__(self) -> None:
        expected = PARAMETERS_OF_IBRAV.get(self.ibrav)
        if expected is None:
            raise StructureError(f"ibrav = {self.ibrav} is not supported")
        if tuple(self.parameters) != expected:
            raise StructureError(
                f"ibrav = {self.ibrav} takes {expected}, got {tuple(self.parameters)}"
            )
        object.__setattr__(
            self, "parameters", {k: float(v) for k, v in self.parameters.items()}
        )

    @property
    def lattice_matrix(self) -> np.ndarray:
        """The lattice vectors (rows, angstrom) exactly as ``pw.x`` builds them."""
        return qe_lattice(self.ibrav, self.parameters)

    def describe(self) -> str:
        values = ", ".join(f"{k} = {v:.6g}" for k, v in self.parameters.items())
        return f"ibrav = {self.ibrav} ({self.symbol}, No. {self.space_group}; {values})"


def _lattice_parameters(lattice: np.ndarray) -> Tuple[float, ...]:
    """(a, b, c, alpha, beta, gamma) of lattice vectors given as rows; angles in degrees."""
    a, b, c = (float(np.linalg.norm(v)) for v in lattice)

    def angle(u, v):
        return math.degrees(
            math.acos(
                max(
                    -1.0,
                    min(
                        1.0,
                        float(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v)),
                    ),
                )
            )
        )

    return (
        a,
        b,
        c,
        angle(lattice[1], lattice[2]),
        angle(lattice[0], lattice[2]),
        angle(lattice[0], lattice[1]),
    )


def _parameters(ibrav: int, conventional: np.ndarray) -> Dict[str, float]:
    """QE's ``A``, ``B``, ... from spglib's standardized conventional cell."""
    a, b, c, alpha, beta, gamma = _lattice_parameters(conventional)
    if ibrav == 5:
        # spglib standardizes R lattices on hexagonal axes; QE wants the
        # rhombohedral edge and angle. Accept rhombohedral axes as well.
        if (
            abs(gamma - 120.0) < 1.0
            and abs(alpha - 90.0) < 1.0
            and abs(beta - 90.0) < 1.0
        ):
            a_rh = math.sqrt(a**2 / 3 + c**2 / 9)
            cos_rh = (2 * c**2 - 3 * a**2) / (2 * c**2 + 6 * a**2)
        else:
            a_rh, cos_rh = a, math.cos(math.radians(alpha))
        return {"A": a_rh, "cosAB": cos_rh}
    values = {
        "A": a,
        "B": b,
        "C": c,
        "cosBC": math.cos(math.radians(alpha)),
        "cosAC": math.cos(math.radians(beta)),
        "cosAB": math.cos(math.radians(gamma)),
    }
    return {key: values[key] for key in PARAMETERS_OF_IBRAV[ibrav]}


def _dataset_value(dataset, key: str):
    return getattr(dataset, key) if hasattr(dataset, key) else dataset[key]


def _spglib_dataset(structure: NormalizedStructure, numbers, symprec: float):
    """spglib's dataset and the number of atoms in the primitive cell."""
    import spglib

    cell = (
        structure.lattice_matrix.tolist(),
        structure.frac_coords().tolist(),
        numbers,
    )
    try:
        dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
        primitive = spglib.find_primitive(cell, symprec=symprec)
    except Exception as exc:  # spglib raises its own error types
        raise StructureError(f"spglib failed on the structure: {exc}") from exc
    if dataset is None or primitive is None:
        raise StructureError(
            f"spglib could not determine the space group (symprec = {symprec})"
        )
    return dataset, len(primitive[1])


def _symmetry(structure: NormalizedStructure, symprec: float):
    """(space-group number, HM symbol, standardized conventional lattice) from spglib.

    Sites are told apart by species label first; when the cell is then not
    primitive (a magnetic supercell), by label and moment.
    """
    labels = structure.species_labels
    numbers = [labels.index(site.label) + 1 for site in structure.sites]
    dataset, primitive_atoms = _spglib_dataset(structure, numbers, symprec)
    if primitive_atoms != structure.num_sites:
        _, species = classify_magnetic_sites(structure)
        for kind, group in enumerate(species, start=1):
            for index in group.site_indices:
                numbers[index] = kind
        dataset, primitive_atoms = _spglib_dataset(structure, numbers, symprec)
        if primitive_atoms != structure.num_sites:
            raise StructureError(
                f"the cell with {structure.num_sites} sites is not a primitive cell "
                f"(spglib finds one with {primitive_atoms}); ibrav /= 0 needs a "
                "primitive cell"
            )
        logger.info(
            "magnetic supercell: sites distinguished by moment for the symmetry"
        )
    return (
        int(_dataset_value(dataset, "number")),
        str(_dataset_value(dataset, "international")),
        np.array(_dataset_value(dataset, "std_lattice"), dtype=float),
    )


def _vectors_of_length(
    old: np.ndarray, length_squared: float, tolerance: float
) -> list[np.ndarray]:
    """Integer lattice vectors with the requested squared length.

    QR makes ``|m @ old|**2`` a sum of three squares. Enumerating the
    coefficients from the last square to the first bounds each coefficient
    by the remaining length, even for an arbitrarily sheared input basis.
    """
    _, upper = np.linalg.qr(old.T)
    coefficients = np.zeros(3, dtype=int)
    matches: list[np.ndarray] = []
    limit = length_squared + tolerance

    def search(index: int, used: float) -> None:
        if index < 0:
            if abs(used - length_squared) <= tolerance:
                matches.append(coefficients.copy())
            return
        remaining = limit - used
        if remaining < 0:
            return
        offset = float(np.dot(upper[index, index + 1 :], coefficients[index + 1 :]))
        diagonal = float(upper[index, index])
        center = -offset / diagonal
        radius = math.sqrt(remaining) / abs(diagonal)
        # Widen the floating-point bounds slightly so a boundary vector is
        # still checked against the exact tolerance at the recursion leaf.
        first = math.ceil(center - radius - 1.0e-10)
        last = math.floor(center + radius + 1.0e-10)
        for value in range(first, last + 1):
            coefficients[index] = value
            search(index - 1, used + (diagonal * value + offset) ** 2)

    search(2, 0.0)
    # Preserve the former lexicographic tie order for symmetry-equivalent
    # bases, so the chosen Cartesian frame stays stable.
    matches.sort(key=tuple)
    return matches


def _basis_change(old: np.ndarray, new: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Integer ``M`` (det +1) and rotation ``R`` with ``new = M @ old @ R``.

    ``old`` and ``new`` are two primitive bases (rows) of the same lattice
    in different Cartesian orientations. Every row of ``new`` is an integer
    combination of the rows of ``old`` up to the rotation, so ``M`` is found
    by matching the metric tensors and ``R`` follows.
    """
    g_old = old @ old.T
    g_new = new @ new.T
    scale = float(np.max(np.abs(np.diag(g_new))))
    tolerance = 2 * LATTICE_TOLERANCE * scale

    candidates = [
        _vectors_of_length(old, float(g_new[i, i]), tolerance) for i in range(3)
    ]

    best = None
    for m1 in candidates[0]:
        for m2 in candidates[1]:
            if abs(m1 @ g_old @ m2 - g_new[0, 1]) > tolerance:
                continue
            for m3 in candidates[2]:
                if (
                    abs(m1 @ g_old @ m3 - g_new[0, 2]) > tolerance
                    or abs(m2 @ g_old @ m3 - g_new[1, 2]) > tolerance
                ):
                    continue
                matrix = np.array([m1, m2, m3], dtype=int)
                if int(round(np.linalg.det(matrix))) != 1:
                    continue
                deviation = (
                    float(np.max(np.abs(matrix @ g_old @ matrix.T - g_new))) / scale
                )
                if best is None or deviation < best[0]:
                    best = (deviation, matrix)
    if best is None or best[0] > LATTICE_TOLERANCE:
        found = f" (closest metric deviation {best[0]:.2e})" if best else ""
        raise StructureError(
            "the lattice of the structure does not match the QE lattice built from "
            f"the standardized cell parameters{found}"
        )
    matrix = best[1]
    approx = np.linalg.inv(matrix @ old) @ new
    u, _, vt = np.linalg.svd(approx)
    rotation = u @ vt  # the orthogonal factor: exact rotation, no strain
    return matrix, rotation


def bravais_structure(
    structure: NormalizedStructure, symprec: float = 1.0e-3
) -> Tuple[NormalizedStructure, BravaisLattice]:
    """Re-express ``structure`` in QE's ``ibrav`` lattice.

    Parameters
    ----------
    structure : NormalizedStructure
        A primitive cell (of the magnetic structure, if it has moments);
        sites may carry moments and labels, both are preserved.
    symprec : float
        Symmetry tolerance in angstrom passed to spglib.

    Returns
    -------
    NormalizedStructure, BravaisLattice
        The structure with QE's lattice vectors, the same sites (labels,
        elements, order) at their new fractional coordinates, and the
        moments rotated with the crystal; and the ``ibrav`` description.

    Raises
    ------
    StructureError
        If the lattice is left-handed, spglib fails, or the QE lattice does
        not reproduce the structure's lattice within ``LATTICE_TOLERANCE``.
    """
    old = structure.lattice_matrix
    if np.linalg.det(old) <= 0.0:
        raise StructureError("the lattice vectors must be right-handed for ibrav /= 0")

    number, symbol, conventional = _symmetry(structure, symprec)
    ibrav = ibrav_of(number, symbol)
    bravais = BravaisLattice(ibrav, _parameters(ibrav, conventional), number, symbol)
    new = bravais.lattice_matrix
    logger.info("Bravais lattice: %s", bravais.describe())

    matrix, rotation = _basis_change(old, new)
    inverse = np.rint(np.linalg.inv(matrix)).astype(int)
    frac = structure.frac_coords() @ inverse
    frac -= np.floor(frac + 1.0e-8)
    frac[np.abs(frac) < 1.0e-8] = 0.0

    sites = []
    for site, coords in zip(structure.sites, frac):
        moment = site.magmom
        if moment is not None:
            moment = MagneticMoment(tuple(moment.as_array() @ rotation))
        sites.append(
            AtomicSite(site.element, tuple(coords), moment, site.label, site.occupancy)
        )
    return NormalizedStructure(new, sites), bravais


__all__ = [
    "IBRAV_OF_LATTICE",
    "LATTICE_TOLERANCE",
    "PARAMETERS_OF_IBRAV",
    "BravaisLattice",
    "bravais_structure",
    "crystal_system",
    "ibrav_of",
    "qe_lattice",
]
