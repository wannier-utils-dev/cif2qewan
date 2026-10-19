"""Code-independent description of a crystal structure.

Units and conventions:

- lattice vectors in angstrom, one vector per row;
- atomic positions as fractional coordinates of the lattice vectors;
- magnetic moments as Cartesian vectors in Bohr magneton, in the same
  Cartesian frame as the lattice vectors.

Nothing here knows about Quantum ESPRESSO or Wannier90; in particular the
QE lattice classification ``ibrav`` belongs to the QE layer.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np

from cif2qewan.exceptions import StructureError

Vector3 = Tuple[float, float, float]
Matrix3 = Tuple[Vector3, Vector3, Vector3]

#: Element symbols in order of atomic number (index + 1 = Z).
ELEMENTS = tuple("""
    H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni
    Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe
    Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg
    Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg
    Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
    """.split())
ELEMENT_SYMBOLS = frozenset(ELEMENTS)


def atomic_number(symbol: str) -> int:
    """Atomic number of an element symbol; raises StructureError if unknown."""
    try:
        return ELEMENTS.index(symbol) + 1
    except ValueError:
        raise StructureError(f"unknown element symbol {symbol!r}") from None


#: Below this magnitude (in Bohr magneton) a moment is treated as zero.
MOMENT_TOLERANCE = 1.0e-3


def _vector3(values: Iterable[float], what: str) -> Vector3:
    """Convert ``values`` to a tuple of three finite floats or raise StructureError."""
    try:
        vector = tuple(float(x) for x in values)
    except (TypeError, ValueError):
        raise StructureError(f"{what} must be three numbers, got {values!r}") from None
    if len(vector) != 3 or not all(math.isfinite(x) for x in vector):
        raise StructureError(f"{what} must be three finite numbers, got {values!r}")
    return vector  # type: ignore[return-value]


@dataclass(frozen=True)
class MagneticMoment:
    """Magnetic moment of one site as a Cartesian vector in Bohr magneton.

    Parameters
    ----------
    vector : sequence of three floats
        (m_x, m_y, m_z) in mu_B, in the Cartesian frame of the lattice vectors.
    """

    vector: Vector3

    def __post_init__(self) -> None:
        object.__setattr__(self, "vector", _vector3(self.vector, "magnetic moment"))

    @classmethod
    def zero(cls) -> "MagneticMoment":
        return cls((0.0, 0.0, 0.0))

    @classmethod
    def from_polar(
        cls, magnitude: float, angle1: float, angle2: float
    ) -> "MagneticMoment":
        """Build a moment from its magnitude and the Quantum ESPRESSO angles.

        Parameters
        ----------
        magnitude : float
            Length of the moment in mu_B (may be negative to flip it).
        angle1 : float
            Polar angle from the z axis, in degrees (QE ``angle1``).
        angle2 : float
            Azimuthal angle in the xy plane from the x axis, in degrees
            (QE ``angle2``).
        """
        theta = math.radians(angle1)
        phi = math.radians(angle2)
        return cls(
            (
                magnitude * math.sin(theta) * math.cos(phi),
                magnitude * math.sin(theta) * math.sin(phi),
                magnitude * math.cos(theta),
            )
        )

    def as_array(self) -> np.ndarray:
        return np.array(self.vector, dtype=float)

    @property
    def magnitude(self) -> float:
        """Length of the moment in mu_B."""
        return float(np.linalg.norm(self.vector))

    def is_zero(self, tol: float = MOMENT_TOLERANCE) -> bool:
        return self.magnitude <= tol

    def direction(self) -> np.ndarray:
        """Unit vector along the moment; raises StructureError for a zero moment."""
        if self.is_zero():
            raise StructureError("a zero magnetic moment has no direction")
        return self.as_array() / self.magnitude

    def angles(self) -> Tuple[float, float]:
        """(angle1, angle2) in degrees, QE convention; (0, 0) for a zero moment."""
        if self.is_zero():
            return 0.0, 0.0
        x, y, z = self.direction()
        angle1 = math.degrees(math.acos(max(-1.0, min(1.0, z))))
        angle2 = math.degrees(math.atan2(y, x)) if math.hypot(x, y) > 1e-12 else 0.0
        return angle1, angle2

    def is_equivalent(
        self, other: "MagneticMoment", tol: float = MOMENT_TOLERANCE
    ) -> bool:
        """True when the two vectors agree component-wise within ``tol`` mu_B."""
        return bool(np.max(np.abs(self.as_array() - other.as_array())) <= tol)


@dataclass(frozen=True)
class AtomicSite:
    """One atom in the unit cell.

    Parameters
    ----------
    element : str
        Element symbol.
    frac_coords : sequence of three floats
        Fractional coordinates with respect to the lattice vectors.
    magmom : MagneticMoment, optional
        Magnetic moment; None for a non-magnetic description.
    label : str, optional
        Species label used to distinguish inequivalent sites of the same
        element (for example ``Fe1``). Defaults to the element symbol.
    occupancy : float
        Site occupancy in (0, 1].
    """

    element: str
    frac_coords: Vector3
    magmom: Optional[MagneticMoment] = None
    label: Optional[str] = None
    occupancy: float = 1.0

    def __post_init__(self) -> None:
        if self.element not in ELEMENT_SYMBOLS:
            raise StructureError(f"unknown element symbol {self.element!r}")
        object.__setattr__(
            self, "frac_coords", _vector3(self.frac_coords, "fractional coordinates")
        )
        if self.label is None:
            object.__setattr__(self, "label", self.element)
        elif not str(self.label).strip():
            raise StructureError("site label must not be empty")
        if not (0.0 < float(self.occupancy) <= 1.0):
            raise StructureError(f"occupancy must be in (0, 1], got {self.occupancy}")
        object.__setattr__(self, "occupancy", float(self.occupancy))

    @property
    def is_magnetic(self) -> bool:
        return self.magmom is not None and not self.magmom.is_zero()

    def wrapped(self) -> "AtomicSite":
        """The same site with fractional coordinates moved into [0, 1)."""
        coords = tuple(x - math.floor(x) for x in self.frac_coords)
        return AtomicSite(self.element, coords, self.magmom, self.label, self.occupancy)


@dataclass(frozen=True)
class NormalizedStructure:
    """A periodic crystal structure in the units of this package.

    Parameters
    ----------
    lattice : 3x3 array-like
        Lattice vectors in angstrom, one per row.
    sites : sequence of AtomicSite
        At least one site. Sites sharing a label must share the element.
    """

    lattice: Matrix3
    sites: Tuple[AtomicSite, ...]

    def __post_init__(self) -> None:
        try:
            matrix = np.array(self.lattice, dtype=float)
        except (TypeError, ValueError):
            raise StructureError("lattice must be a 3x3 array of numbers") from None
        if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
            raise StructureError("lattice must be a 3x3 array of finite numbers")
        if abs(np.linalg.det(matrix)) < 1.0e-8:
            raise StructureError("lattice vectors are linearly dependent")
        object.__setattr__(
            self, "lattice", tuple(tuple(float(x) for x in row) for row in matrix)
        )

        sites = tuple(self.sites)
        if not sites:
            raise StructureError("a structure needs at least one site")
        if not all(isinstance(site, AtomicSite) for site in sites):
            raise StructureError("sites must be AtomicSite instances")
        element_of_label = {}
        for site in sites:
            previous = element_of_label.setdefault(site.label, site.element)
            if previous != site.element:
                raise StructureError(
                    f"label {site.label!r} is used for both {previous} and {site.element}"
                )
        object.__setattr__(self, "sites", sites)

    @property
    def lattice_matrix(self) -> np.ndarray:
        """Lattice vectors in angstrom as a (3, 3) array, one per row."""
        return np.array(self.lattice, dtype=float)

    @property
    def volume(self) -> float:
        """Cell volume in cubic angstrom."""
        return float(abs(np.linalg.det(self.lattice_matrix)))

    @property
    def num_sites(self) -> int:
        return len(self.sites)

    @property
    def elements(self) -> Tuple[str, ...]:
        """Distinct elements in order of first appearance."""
        return tuple(dict.fromkeys(site.element for site in self.sites))

    @property
    def species_labels(self) -> Tuple[str, ...]:
        """Distinct site labels in order of first appearance."""
        return tuple(dict.fromkeys(site.label for site in self.sites))  # type: ignore

    @property
    def is_magnetic(self) -> bool:
        return any(site.is_magnetic for site in self.sites)

    def frac_coords(self) -> np.ndarray:
        """Fractional coordinates as an (n, 3) array."""
        return np.array([site.frac_coords for site in self.sites], dtype=float)

    def cart_coords(self) -> np.ndarray:
        """Cartesian coordinates in angstrom as an (n, 3) array."""
        return self.frac_coords() @ self.lattice_matrix

    def frac_to_cart(self, frac: Sequence[float]) -> np.ndarray:
        return np.asarray(frac, dtype=float) @ self.lattice_matrix

    def reciprocal_lattice(self) -> np.ndarray:
        """Reciprocal lattice vectors (rows) in 1/angstrom, including the 2*pi."""
        return 2.0 * np.pi * np.linalg.inv(self.lattice_matrix).T

    def with_sites(self, sites: Iterable[AtomicSite]) -> "NormalizedStructure":
        return NormalizedStructure(self.lattice, tuple(sites))
