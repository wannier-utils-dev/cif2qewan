"""Read CIF, MCIF and the other formats pymatgen understands."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from cif2qewan.exceptions import StructureError
from cif2qewan.structure.model import AtomicSite, MagneticMoment, NormalizedStructure
from cif2qewan.structure.readers.base import PathLike


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
    A ``magmom`` site property is kept as a Cartesian moment in Bohr
    magneton, taken from pymatgen's global frame. Reading of MagCIF files
    and the classification of magnetic sites are the subject of
    DEVELOPMENT_PLAN.md Step 9; until then the moments are passed through
    as pymatgen provides them.
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
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        try:
            primitive = SpacegroupAnalyzer(
                structure, symprec=self.symprec
            ).find_primitive()
        except Exception as exc:
            raise StructureError(f"primitive cell search failed: {exc}") from exc
        return primitive if primitive is not None else structure

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
    """(element, occupancy) of a site occupied by one species."""
    species = list(site.species.items())
    if len(species) != 1:
        raise StructureError(
            f"site at {tuple(site.frac_coords)} is shared by several species "
            f"({site.species_string}); mixed occupancy is not supported"
        )
    specie, occupancy = species[0]
    symbol = (
        getattr(specie, "symbol", None) or getattr(specie, "element", specie).symbol
    )
    return str(symbol), float(occupancy)


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
