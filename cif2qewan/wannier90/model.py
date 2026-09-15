"""Typed model of a Wannier90 ``.win`` input.

The model holds the cell, atoms, k mesh, projections and pass-through
parameters; :mod:`cif2qewan.wannier90.writer` renders it. The number
of Wannier functions implied by the projections is checked against
``num_wann`` when every projection site can be resolved to atoms.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple, Union

import numpy as np

from cif2qewan.exceptions import InputModelError

Vector3 = Tuple[float, float, float]

#: Wannier functions per projection orbital (Wannier90 user guide, sec. 3).
ORBITAL_SIZES = {
    "s": 1,
    "p": 3,
    "d": 5,
    "f": 7,
    "sp": 2,
    "sp2": 3,
    "sp3": 4,
    "sp3d": 5,
    "sp3d2": 6,
    "px": 1,
    "py": 1,
    "pz": 1,
    "dxy": 1,
    "dxz": 1,
    "dyz": 1,
    "dx2-y2": 1,
    "dz2": 1,
}

ParameterValue = Union[bool, int, float, str]


def _vector3(values, what: str) -> Vector3:
    vector = tuple(float(x) for x in values)
    if len(vector) != 3 or not all(math.isfinite(x) for x in vector):
        raise InputModelError(f"{what} must be three finite numbers, got {values!r}")
    return vector  # type: ignore[return-value]


@dataclass(frozen=True)
class Projection:
    """One ``projections`` line: ``site:orbital,orbital,...``.

    ``site`` is a species label (``Fe``), a site label, or any other
    Wannier90 site specification such as ``f=0.5,0.5,0.5``.
    """

    site: str
    orbitals: Tuple[str, ...]

    def __post_init__(self) -> None:
        if not str(self.site).strip():
            raise InputModelError("projection site must not be empty")
        orbitals = tuple(str(o).strip() for o in self.orbitals)
        if not orbitals:
            raise InputModelError(f"projection on {self.site} has no orbitals")
        for orbital in orbitals:
            if orbital not in ORBITAL_SIZES:
                raise InputModelError(
                    f"unknown projection orbital {orbital!r} on {self.site}"
                )
        object.__setattr__(self, "orbitals", orbitals)

    @property
    def num_wann_per_site(self) -> int:
        """Wannier functions this projection gives on one site (without spin)."""
        return sum(ORBITAL_SIZES[o] for o in self.orbitals)


@dataclass(frozen=True)
class AtomFrac:
    """One ``atoms_frac`` row."""

    label: str
    coords: Vector3

    def __post_init__(self) -> None:
        if not str(self.label).strip():
            raise InputModelError("atom label must not be empty")
        object.__setattr__(self, "coords", _vector3(self.coords, "atom position"))


@dataclass(frozen=True)
class KPathSegment:
    """One ``kpoint_path`` line: two labelled points in fractional coordinates."""

    start_label: str
    start: Vector3
    end_label: str
    end: Vector3

    def __post_init__(self) -> None:
        for label in (self.start_label, self.end_label):
            if not str(label).strip():
                raise InputModelError("k-path labels must not be empty")
        object.__setattr__(self, "start", _vector3(self.start, "k-path point"))
        object.__setattr__(self, "end", _vector3(self.end, "k-path point"))


@dataclass
class Wannier90Input:
    """A Wannier90 ``.win`` file.

    Parameters
    ----------
    num_wann, num_bands : int
        Wannier functions and bands passed to Wannier90 (spin included).
    unit_cell_cart : 3x3 array-like
        Lattice vectors in angstrom, one per row.
    atoms_frac : sequence of AtomFrac
    mp_grid : (int, int, int)
    kpoints : sequence of three-vectors
        Fractional coordinates of the full mesh; ``len == prod(mp_grid)``.
    projections : sequence of Projection
    spinors : bool
    exclude_bands : (first, last), optional
        Inclusive 1-based range of excluded bands.
    kpoint_path : sequence of KPathSegment
    parameters : dict
        Pass-through keywords (``dis_num_iter``, ``write_hr``, ...). They
        must not repeat a keyword that is a field of this model.
    """

    num_wann: int
    num_bands: int
    unit_cell_cart: Tuple[Vector3, Vector3, Vector3]
    atoms_frac: Tuple[AtomFrac, ...]
    mp_grid: Tuple[int, int, int]
    kpoints: Tuple[Vector3, ...]
    projections: Tuple[Projection, ...]
    spinors: bool = False
    exclude_bands: Optional[Tuple[int, int]] = None
    kpoint_path: Tuple[KPathSegment, ...] = ()
    parameters: Dict[str, ParameterValue] = field(default_factory=dict)

    RESERVED = frozenset(
        {
            "num_wann",
            "num_bands",
            "unit_cell_cart",
            "atoms_frac",
            "atoms_cart",
            "mp_grid",
            "kpoints",
            "projections",
            "spinors",
            "exclude_bands",
            "kpoint_path",
        }
    )

    def __post_init__(self) -> None:
        self.num_wann = int(self.num_wann)
        self.num_bands = int(self.num_bands)
        if self.num_wann < 1:
            raise InputModelError(f"num_wann must be positive, got {self.num_wann}")
        if self.num_bands < self.num_wann:
            raise InputModelError(
                f"num_bands = {self.num_bands} is smaller than num_wann = {self.num_wann}"
            )
        if self.spinors and self.num_wann % 2:
            raise InputModelError("spinors = .true. needs an even num_wann")

        cell = np.array(self.unit_cell_cart, dtype=float)
        if cell.shape != (3, 3) or not np.all(np.isfinite(cell)):
            raise InputModelError("unit_cell_cart needs three finite vectors")
        if abs(np.linalg.det(cell)) < 1.0e-12:
            raise InputModelError("unit_cell_cart vectors are linearly dependent")
        self.unit_cell_cart = tuple(tuple(float(x) for x in row) for row in cell)

        self.atoms_frac = tuple(self.atoms_frac)
        if not self.atoms_frac or not all(
            isinstance(a, AtomFrac) for a in self.atoms_frac
        ):
            raise InputModelError("atoms_frac needs at least one AtomFrac")

        mp_grid = tuple(int(n) for n in self.mp_grid)
        if len(mp_grid) != 3 or any(n < 1 for n in mp_grid):
            raise InputModelError(
                f"mp_grid must be three positive integers, got {self.mp_grid}"
            )
        self.mp_grid = mp_grid  # type: ignore[assignment]
        self.kpoints = tuple(_vector3(k, "k point") for k in self.kpoints)
        if len(self.kpoints) != int(np.prod(mp_grid)):
            raise InputModelError(
                f"{len(self.kpoints)} kpoints do not fill mp_grid {mp_grid}"
            )

        self.projections = tuple(self.projections)
        if not all(isinstance(p, Projection) for p in self.projections):
            raise InputModelError("projections must be Projection instances")
        implied = self.projected_num_wann()
        if implied is not None and implied != self.num_wann:
            raise InputModelError(
                f"projections give {implied} Wannier functions but num_wann = {self.num_wann}"
            )

        if self.exclude_bands is not None:
            first, last = (int(x) for x in self.exclude_bands)
            if first < 1 or last < first:
                raise InputModelError(
                    f"invalid exclude_bands range {self.exclude_bands}"
                )
            self.exclude_bands = (first, last)

        self.kpoint_path = tuple(self.kpoint_path)
        if not all(isinstance(s, KPathSegment) for s in self.kpoint_path):
            raise InputModelError("kpoint_path must be KPathSegment instances")

        for key in self.parameters:
            if key.lower() in self.RESERVED:
                raise InputModelError(f"{key} must be set as a field, not a parameter")

    @property
    def num_exclude(self) -> int:
        if self.exclude_bands is None:
            return 0
        return self.exclude_bands[1] - self.exclude_bands[0] + 1

    def projected_num_wann(self) -> Optional[int]:
        """Wannier functions implied by the projections, or None if unknown.

        Projections whose site is an atom label are multiplied by the number
        of such atoms and, for spinors, by two. If any site cannot be matched
        to ``atoms_frac`` (for example ``f=...`` or ``random``) the count is
        unknown and None is returned.
        """
        total = 0
        for projection in self.projections:
            count = sum(1 for atom in self.atoms_frac if atom.label == projection.site)
            if count == 0:
                return None
            total += count * projection.num_wann_per_site
        return total * (2 if self.spinors else 1)
