"""Typed model of Quantum ESPRESSO input files.

The model mirrors the structure of a ``pw.x`` input: Fortran namelists
followed by cards. It stores typed values only; turning them into text is
the job of the writer (DEVELOPMENT_PLAN.md Step 4), and choosing the
physical parameters is the job of the workflow builder (Step 6).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple, Union

import numpy as np

from cif2qewan.exceptions import InputModelError

Vector3 = Tuple[float, float, float]


@dataclass(frozen=True)
class RawValue:
    """A namelist value kept exactly as written, e.g. ``1.0d-8`` or ``'m-p'``.

    An empty literal renders as ``key = `` with nothing after the equals
    sign; the 0.2.x ``proj.in`` writes ``Emax`` and ``Emin`` that way.
    """

    literal: str

    def __post_init__(self) -> None:
        if not isinstance(self.literal, str):
            raise InputModelError(
                f"a raw namelist value must be a string, got {self.literal!r}"
            )


NamelistValue = Union[bool, int, float, str, RawValue]


@dataclass
class Namelist:
    """One Fortran namelist: an ordered mapping from variable to typed value.

    Parameters
    ----------
    name : str
        Namelist name without the ampersand, e.g. ``control``.
    entries : dict
        Variable name (may include an index, ``starting_magnetization(1)``)
        to value. ``bool`` renders as ``.true.``/``.false.``, ``str`` is
        quoted, numbers use their Python repr, ``RawValue`` is copied.
    """

    name: str
    entries: Dict[str, NamelistValue] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not str(self.name).strip():
            raise InputModelError("namelist name must not be empty")
        for key, value in self.entries.items():
            if not str(key).strip():
                raise InputModelError(f"namelist {self.name}: empty variable name")
            if not isinstance(value, (bool, int, float, str, RawValue)):
                raise InputModelError(
                    f"namelist {self.name}: unsupported value {value!r} for {key}"
                )
            if isinstance(value, float) and not math.isfinite(value):
                raise InputModelError(f"namelist {self.name}: {key} is not finite")

    def get(self, key: str, default: Optional[NamelistValue] = None):
        return self.entries.get(key, default)


@dataclass(frozen=True)
class Species:
    """One ATOMIC_SPECIES row: label, atomic mass (amu) and pseudopotential file."""

    label: str
    mass: float
    pseudo_file: str

    def __post_init__(self) -> None:
        if not str(self.label).strip():
            raise InputModelError("species label must not be empty")
        if not (float(self.mass) > 0.0):
            raise InputModelError(f"species {self.label}: mass must be positive")
        if not str(self.pseudo_file).strip():
            raise InputModelError(f"species {self.label}: pseudopotential file missing")


@dataclass(frozen=True)
class CellParameters:
    """CELL_PARAMETERS card: three lattice vectors (rows) in the given unit.

    ``option`` is ``alat``, ``bohr`` or ``angstrom``. With ``alat`` the
    vectors are in units of ``A`` or ``celldm(1)`` of the &system namelist.
    """

    option: str
    vectors: Tuple[Vector3, Vector3, Vector3]

    def __post_init__(self) -> None:
        if self.option not in ("alat", "bohr", "angstrom"):
            raise InputModelError(f"unknown CELL_PARAMETERS option {self.option!r}")
        matrix = np.array(self.vectors, dtype=float)
        if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)):
            raise InputModelError("CELL_PARAMETERS needs three finite vectors")
        if abs(np.linalg.det(matrix)) < 1.0e-12:
            raise InputModelError("CELL_PARAMETERS vectors are linearly dependent")
        object.__setattr__(
            self, "vectors", tuple(tuple(float(x) for x in row) for row in matrix)
        )

    def matrix(self) -> np.ndarray:
        return np.array(self.vectors, dtype=float)


@dataclass(frozen=True)
class AtomicPosition:
    """One ATOMIC_POSITIONS row."""

    label: str
    coords: Vector3

    def __post_init__(self) -> None:
        if not str(self.label).strip():
            raise InputModelError("atomic position label must not be empty")
        coords = tuple(float(x) for x in self.coords)
        if len(coords) != 3 or not all(math.isfinite(x) for x in coords):
            raise InputModelError(
                f"position of {self.label} must be three finite numbers"
            )
        object.__setattr__(self, "coords", coords)


@dataclass(frozen=True)
class KPointsAutomatic:
    """K_POINTS automatic: Monkhorst-Pack mesh and (0/1) shifts."""

    mesh: Tuple[int, int, int]
    shift: Tuple[int, int, int] = (0, 0, 0)

    def __post_init__(self) -> None:
        mesh = tuple(int(n) for n in self.mesh)
        shift = tuple(int(s) for s in self.shift)
        if len(mesh) != 3 or any(n < 1 for n in mesh):
            raise InputModelError(
                f"k mesh must be three positive integers, got {self.mesh}"
            )
        if len(shift) != 3 or any(s not in (0, 1) for s in shift):
            raise InputModelError(f"k shift must be three 0/1 values, got {self.shift}")
        object.__setattr__(self, "mesh", mesh)
        object.__setattr__(self, "shift", shift)

    @property
    def num_points(self) -> int:
        return int(np.prod(self.mesh))


@dataclass(frozen=True)
class KPointsList:
    """K_POINTS crystal / tpiba: explicit points with weights."""

    option: str
    points: Tuple[Tuple[float, float, float, float], ...]

    def __post_init__(self) -> None:
        if self.option not in ("crystal", "tpiba"):
            raise InputModelError(f"unknown K_POINTS list option {self.option!r}")
        points = tuple(tuple(float(x) for x in p) for p in self.points)
        if not points or any(len(p) != 4 for p in points):
            raise InputModelError("each explicit k point needs kx, ky, kz and a weight")
        if any(p[3] < 0.0 for p in points):
            raise InputModelError("k-point weights must not be negative")
        object.__setattr__(self, "points", points)

    @property
    def num_points(self) -> int:
        return len(self.points)


@dataclass(frozen=True)
class KPathPoint:
    """One vertex of a K_POINTS *_b path: coordinates, divisions to the next, label."""

    coords: Vector3
    ndiv: int
    label: str = ""

    def __post_init__(self) -> None:
        coords = tuple(float(x) for x in self.coords)
        if len(coords) != 3 or not all(math.isfinite(x) for x in coords):
            raise InputModelError("a k-path vertex needs three finite coordinates")
        if int(self.ndiv) < 0:
            raise InputModelError("k-path divisions must not be negative")
        object.__setattr__(self, "coords", coords)
        object.__setattr__(self, "ndiv", int(self.ndiv))


@dataclass(frozen=True)
class KPointsPath:
    """K_POINTS crystal_b / tpiba_b: a band-structure path.

    A vertex with ``ndiv == 0`` marks a jump: no line is drawn from it to the
    next vertex.
    """

    option: str
    points: Tuple[KPathPoint, ...]

    def __post_init__(self) -> None:
        if self.option not in ("crystal_b", "tpiba_b"):
            raise InputModelError(f"unknown K_POINTS path option {self.option!r}")
        points = tuple(self.points)
        if len(points) < 2 or not all(isinstance(p, KPathPoint) for p in points):
            raise InputModelError("a k path needs at least two KPathPoint vertices")
        object.__setattr__(self, "points", points)

    @property
    def num_points(self) -> int:
        return len(self.points)

    def segments(self) -> Tuple[Tuple[KPathPoint, KPathPoint], ...]:
        """Consecutive vertex pairs that are connected by a line."""
        return tuple(
            (a, b) for a, b in zip(self.points[:-1], self.points[1:]) if a.ndiv > 0
        )


KPoints = Union[KPointsAutomatic, KPointsList, KPointsPath]


@dataclass
class PwInput:
    """A complete ``pw.x`` input.

    ``system`` may carry ``nat`` and ``ntyp``; when present they must match
    the cards. ``ibrav = 0`` requires ``cell``.
    """

    control: Namelist
    system: Namelist
    electrons: Namelist
    species: Tuple[Species, ...]
    positions: Tuple[AtomicPosition, ...]
    kpoints: KPoints
    cell: Optional[CellParameters] = None
    positions_option: str = "crystal"

    def __post_init__(self) -> None:
        for namelist, name in (
            (self.control, "control"),
            (self.system, "system"),
            (self.electrons, "electrons"),
        ):
            if not isinstance(namelist, Namelist) or namelist.name != name:
                raise InputModelError(f"expected a &{name} namelist, got {namelist!r}")
        self.species = tuple(self.species)
        self.positions = tuple(self.positions)
        if not self.species:
            raise InputModelError("ATOMIC_SPECIES must not be empty")
        if not self.positions:
            raise InputModelError("ATOMIC_POSITIONS must not be empty")
        if self.positions_option not in ("alat", "bohr", "angstrom", "crystal"):
            raise InputModelError(
                f"unknown ATOMIC_POSITIONS option {self.positions_option!r}"
            )
        if not isinstance(self.kpoints, (KPointsAutomatic, KPointsList, KPointsPath)):
            raise InputModelError(f"unsupported K_POINTS {self.kpoints!r}")

        labels = [s.label for s in self.species]
        if len(set(labels)) != len(labels):
            raise InputModelError("duplicate species labels in ATOMIC_SPECIES")
        for position in self.positions:
            if position.label not in labels:
                raise InputModelError(
                    f"position label {position.label!r} has no ATOMIC_SPECIES entry"
                )

        nat = self.system.get("nat")
        if nat is not None and int(nat) != len(self.positions):
            raise InputModelError(f"nat = {nat} but {len(self.positions)} positions")
        ntyp = self.system.get("ntyp")
        if ntyp is not None and int(ntyp) != len(self.species):
            raise InputModelError(f"ntyp = {ntyp} but {len(self.species)} species")
        ibrav = self.system.get("ibrav")
        if ibrav is not None and int(ibrav) == 0 and self.cell is None:
            raise InputModelError("ibrav = 0 requires CELL_PARAMETERS")
        if self.cell is not None and self.cell.option == "alat":
            if self.system.get("A") is None and self.system.get("celldm(1)") is None:
                raise InputModelError("CELL_PARAMETERS alat requires A or celldm(1)")

    @property
    def calculation(self) -> str:
        value = self.control.get("calculation", "scf")
        return value.literal.strip("'\"") if isinstance(value, RawValue) else str(value)

    @property
    def nat(self) -> int:
        return len(self.positions)

    @property
    def ntyp(self) -> int:
        return len(self.species)


@dataclass
class NamelistInput:
    """An input made of namelists only (``bands.x``, ``pw2wannier90.x``, ...)."""

    namelists: Tuple[Namelist, ...]

    def __post_init__(self) -> None:
        self.namelists = tuple(self.namelists)
        if not self.namelists or not all(
            isinstance(n, Namelist) for n in self.namelists
        ):
            raise InputModelError("NamelistInput needs at least one Namelist")
        names = [n.name for n in self.namelists]
        if len(set(names)) != len(names):
            raise InputModelError(f"duplicate namelists in input: {names}")

    def namelist(self, name: str) -> Namelist:
        for namelist in self.namelists:
            if namelist.name == name:
                return namelist
        raise KeyError(name)
