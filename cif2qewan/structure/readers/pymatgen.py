"""Read CIF, MCIF and the other formats pymatgen understands."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np

from cif2qewan.exceptions import StructureError
from cif2qewan.structure.model import AtomicSite, MagneticMoment, NormalizedStructure
from cif2qewan.structure.readers.base import PathLike

logger = logging.getLogger(__name__)


class PymatgenReader:
    """Read a structure with pymatgen and normalize it.

    Parameters
    ----------
    primitive : bool
        Reduce to the primitive cell (as cif2cell does by default). The
        Cartesian frame of the input is kept.
    symprec : float
        Symmetry tolerance in angstrom used for the primitive-cell search.

    Notes
    -----
    A ``magmom`` site property (MagCIF) is kept as a Cartesian moment in
    Bohr magneton, taken from pymatgen's global frame. The classification
    of magnetic sites happens later, in :mod:`cif2qewan.structure.magnetism`.
    """

    def __init__(self, primitive: bool = True, symprec: float = 1.0e-3) -> None:
        self.primitive = primitive
        self.symprec = symprec

    def read(self, path: PathLike) -> NormalizedStructure:
        from pymatgen.core import Structure

        path = Path(path)
        if not path.is_file():
            raise StructureError(f"structure file not found: {path}")
        try:
            structure = Structure.from_file(str(path))
        except Exception as exc:  # pymatgen raises many different types
            raise StructureError(f"pymatgen could not read {path}: {exc}") from exc
        if self.primitive:
            structure = self._primitive(structure)
        return self.from_pymatgen(structure)

    def _primitive(self, structure):
        """Primitive cell that preserves the magnetic structure, else the input.

        pymatgen's symmetry search ignores site properties, so the chemical
        primitive cell can fold sites with different moments onto each other
        (an antiferromagnet). The candidate cell is accepted only when every
        original site maps onto a candidate site with an equivalent moment.
        """
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        magmoms = structure.site_properties.get("magmom")
        chemical = (
            structure.copy(site_properties={"magmom": None})
            if magmoms is not None
            else structure
        )
        try:
            candidate = SpacegroupAnalyzer(
                chemical, symprec=self.symprec
            ).find_primitive()
        except Exception as exc:
            raise StructureError(f"primitive cell search failed: {exc}") from exc
        if candidate is None:
            return structure
        if len(candidate) == len(structure):
            return structure  # already primitive: keep the input's lattice vectors
        if magmoms is None:
            return candidate

        moments = [_magnetic_moment(m) for m in magmoms]
        assigned = _fold_moments(structure, moments, candidate, self.symprec)
        if assigned is None:
            logger.info(
                "keeping the %d-site cell: the primitive cell would fold different moments",
                len(structure),
            )
            return structure
        candidate = candidate.copy(
            site_properties={"magmom": [m.vector for m in assigned]}
        )
        return candidate

    @staticmethod
    def from_pymatgen(structure) -> NormalizedStructure:
        """Convert a ``pymatgen.core.Structure`` to a NormalizedStructure."""
        magmoms = structure.site_properties.get("magmom")
        sites = []
        for index, site in enumerate(structure):
            element, occupancy = _single_species(site)
            magmom = _magnetic_moment(magmoms[index]) if magmoms is not None else None
            sites.append(
                AtomicSite(
                    element=element,
                    frac_coords=tuple(float(x) for x in site.frac_coords),
                    magmom=magmom,
                    occupancy=occupancy,
                )
            )
        return NormalizedStructure(
            np.array(structure.lattice.matrix, dtype=float), sites
        )


def _single_species(site):
    """(element, occupancy) of a fully occupied site with one species."""
    species = list(site.species.items())
    if len(species) != 1:
        raise StructureError(
            f"site at {tuple(site.frac_coords)} is shared by several species "
            f"({site.species_string}); mixed occupancy is not supported"
        )
    specie, occupancy = species[0]
    occupancy = float(occupancy)
    if abs(occupancy - 1.0) > 1.0e-8:
        raise StructureError(
            f"site at {tuple(site.frac_coords)} has occupancy {occupancy}; "
            "partial occupancy is not supported"
        )
    symbol = (
        getattr(specie, "symbol", None) or getattr(specie, "element", specie).symbol
    )
    return str(symbol), occupancy


def _magnetic_moment(value) -> Optional[MagneticMoment]:
    """Cartesian moment from a pymatgen Magmom, a 3-vector or a collinear scalar."""
    if value is None:
        return None
    from pymatgen.electronic_structure.core import Magmom

    try:
        vector: Sequence[float] = Magmom(value).global_moment
    except Exception as exc:
        raise StructureError(f"cannot interpret magnetic moment {value!r}") from exc
    return MagneticMoment(tuple(float(x) for x in vector))


def _fold_moments(
    structure, moments, candidate, symprec: float
) -> Optional[List[MagneticMoment]]:
    """Moments of the candidate sites, or None if the folding is inconsistent."""
    inverse = np.linalg.inv(candidate.lattice.matrix)
    candidate_frac = np.array([site.frac_coords for site in candidate])
    assigned: List[Optional[MagneticMoment]] = [None] * len(candidate)
    for site, moment in zip(structure, moments):
        frac = site.coords @ inverse
        delta = candidate_frac - frac
        delta -= np.rint(delta)
        cart = np.abs(delta @ candidate.lattice.matrix)
        matches = [
            j
            for j in range(len(candidate))
            if np.max(cart[j]) <= max(symprec, 1e-6)
            and candidate[j].species == site.species
        ]
        if len(matches) != 1:
            return None
        j = matches[0]
        moment = moment if moment is not None else MagneticMoment.zero()
        if assigned[j] is None:
            assigned[j] = moment
        elif not assigned[j].is_equivalent(moment):
            return None
    if any(m is None for m in assigned):
        return None
    return assigned  # type: ignore[return-value]
