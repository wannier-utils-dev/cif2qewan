"""Classification of magnetic structures and their QE representation.

Sites are grouped into magnetic species by element (or existing label) and
by the equivalence of their Cartesian moments within a tolerance; exact
floating-point equality is never used. The order of a structure is
``nonmagnetic``, ``collinear`` (all moments parallel or antiparallel within
an angular tolerance) or ``noncollinear``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from cif2qewan.exceptions import StructureError
from cif2qewan.structure.model import (
    MOMENT_TOLERANCE,
    AtomicSite,
    MagneticMoment,
    NormalizedStructure,
)

NONMAGNETIC = "nonmagnetic"
COLLINEAR = "collinear"
NONCOLLINEAR = "noncollinear"

#: Moments whose directions differ by less than this are treated as (anti)parallel.
COLLINEAR_ANGLE_TOLERANCE_DEG = 5.0


def magnetic_order(
    structure: NormalizedStructure,
    tol: float = MOMENT_TOLERANCE,
    angle_tol_deg: float = COLLINEAR_ANGLE_TOLERANCE_DEG,
) -> str:
    """``nonmagnetic``, ``collinear`` or ``noncollinear``."""
    directions = [
        site.magmom.direction()
        for site in structure.sites
        if site.magmom is not None and not site.magmom.is_zero(tol)
    ]
    if not directions:
        return NONMAGNETIC
    reference = directions[0]
    sin_tol = math.sin(math.radians(angle_tol_deg))
    for direction in directions[1:]:
        if np.linalg.norm(np.cross(reference, direction)) > sin_tol:
            return NONCOLLINEAR
    return COLLINEAR


@dataclass(frozen=True)
class MagneticSpecies:
    """A group of sites with the same element and equivalent moments."""

    label: str
    element: str
    moment: MagneticMoment
    site_indices: Tuple[int, ...]

    @property
    def count(self) -> int:
        return len(self.site_indices)


def classify_magnetic_sites(
    structure: NormalizedStructure, tol: float = MOMENT_TOLERANCE
) -> Tuple[NormalizedStructure, List[MagneticSpecies]]:
    """Split the species by magnetic moment and return the relabelled structure.

    Sites sharing a label whose moments are equivalent within ``tol`` form
    one species. A label that splits into several groups gets a running
    index appended (``Mn`` -> ``Mn1``, ``Mn2``, ...); labels that do not
    split are kept. Sites without a moment count as zero moment. The order
    of species follows the first appearance of each group.
    """
    groups: Dict[str, List[Tuple[MagneticMoment, List[int]]]] = {}
    for index, site in enumerate(structure.sites):
        moment = site.magmom if site.magmom is not None else MagneticMoment.zero()
        buckets = groups.setdefault(site.label, [])  # type: ignore[arg-type]
        for representative, members in buckets:
            if representative.is_equivalent(moment, tol):
                members.append(index)
                break
        else:
            buckets.append((moment, [index]))

    label_of_site: Dict[int, str] = {}
    species: List[Tuple[int, MagneticSpecies]] = []
    for base, buckets in groups.items():
        split = len(buckets) > 1
        for number, (moment, members) in enumerate(buckets, start=1):
            label = f"{base}{number}" if split else base
            if split and any(label == other for other in groups if other != base):
                raise StructureError(
                    f"magnetic species label {label!r} clashes with an existing label"
                )
            element = structure.sites[members[0]].element
            species.append(
                (members[0], MagneticSpecies(label, element, moment, tuple(members)))
            )
            for index in members:
                label_of_site[index] = label

    species.sort(key=lambda item: item[0])
    sites = tuple(
        AtomicSite(
            site.element,
            site.frac_coords,
            site.magmom,
            label_of_site[i],
            site.occupancy,
        )
        for i, site in enumerate(structure.sites)
    )
    return structure.with_sites(sites), [item[1] for item in species]


@dataclass(frozen=True)
class QEMagnetization:
    """``starting_magnetization``, ``angle1`` and ``angle2`` of one species."""

    starting_magnetization: float
    angle1: float
    angle2: float

    @property
    def is_zero(self) -> bool:
        return self.starting_magnetization == 0.0


def qe_magnetization(
    species: Sequence[MagneticSpecies], order: str, tol: float = MOMENT_TOLERANCE
) -> Dict[str, QEMagnetization]:
    """QE starting moments per species label.

    The magnitude is relative to the largest moment (``1.0`` for the largest
    species, ``0.0`` for non-magnetic ones). For a collinear order the sign
    gives the direction along the common axis and the angles are those of
    the axis (the direction of the first magnetic species); for a
    noncollinear order every species carries its own polar angles.
    """
    magnitudes = [s.moment.magnitude for s in species]
    largest = max(magnitudes) if magnitudes else 0.0
    if order == NONMAGNETIC or largest <= tol:
        return {s.label: QEMagnetization(0.0, 0.0, 0.0) for s in species}

    axis = None
    if order == COLLINEAR:
        first = next(s for s in species if not s.moment.is_zero(tol))
        axis = first.moment

    result = {}
    for s in species:
        if s.moment.is_zero(tol):
            result[s.label] = QEMagnetization(0.0, 0.0, 0.0)
            continue
        relative = s.moment.magnitude / largest
        if axis is not None:
            sign = (
                1.0
                if float(np.dot(axis.direction(), s.moment.direction())) >= 0.0
                else -1.0
            )
            angle1, angle2 = axis.angles()
            result[s.label] = QEMagnetization(sign * relative, angle1, angle2)
        else:
            angle1, angle2 = s.moment.angles()
            result[s.label] = QEMagnetization(relative, angle1, angle2)
    return result
