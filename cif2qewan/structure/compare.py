"""Physical comparison of two structures.

Two structures are equivalent when they describe the same crystal in the
same Cartesian frame: the lattices generate the same set of translations
(possibly with different primitive vectors), every site of one has a
counterpart in the other with the same element, occupancy and magnetic
moment at the same position up to a lattice translation.
"""

from __future__ import annotations

from typing import List

import numpy as np

from cif2qewan.structure.model import MOMENT_TOLERANCE, NormalizedStructure


def compare_structures(
    a: NormalizedStructure,
    b: NormalizedStructure,
    length_tol: float = 1.0e-4,
    frac_tol: float = 1.0e-5,
    moment_tol: float = MOMENT_TOLERANCE,
) -> List[str]:
    """Return the list of differences between ``a`` and ``b``; empty when equivalent.

    Parameters
    ----------
    length_tol : float
        Tolerance in angstrom on the lattice relation.
    frac_tol : float
        Tolerance on fractional coordinates (in the cell of ``a``).
    moment_tol : float
        Tolerance in Bohr magneton on magnetic moments.
    """
    differences: List[str] = []

    la = a.lattice_matrix
    lb = b.lattice_matrix
    if abs(a.volume - b.volume) > length_tol * max(a.volume, b.volume) ** (2.0 / 3.0):
        differences.append(f"volumes differ: {a.volume:.6f} vs {b.volume:.6f} A^3")
        return differences

    # b's vectors expressed in a's: must be an integer unimodular matrix
    transform = lb @ np.linalg.inv(la)
    rounded = np.rint(transform)
    if np.max(np.abs(transform - rounded) @ np.abs(la)) > length_tol:
        differences.append("lattices are not the same set of translations")
        return differences
    if abs(abs(np.linalg.det(rounded)) - 1.0) > 1.0e-8:
        differences.append("one cell is a supercell of the other")
        return differences

    if a.num_sites != b.num_sites:
        differences.append(f"site counts differ: {a.num_sites} vs {b.num_sites}")
        return differences

    frac_b_in_a = b.cart_coords() @ np.linalg.inv(la)
    unmatched = list(range(a.num_sites))
    for j, site_b in enumerate(b.sites):
        found = None
        for i in unmatched:
            site_a = a.sites[i]
            if site_a.element != site_b.element:
                continue
            if abs(site_a.occupancy - site_b.occupancy) > 1.0e-6:
                continue
            delta = np.asarray(site_a.frac_coords) - frac_b_in_a[j]
            delta -= np.rint(delta)
            if np.max(np.abs(delta)) > frac_tol:
                continue
            if not _moments_equivalent(site_a, site_b, moment_tol):
                continue
            found = i
            break
        if found is None:
            differences.append(
                f"no counterpart for {site_b.element} at {tuple(round(x, 6) for x in site_b.frac_coords)}"
            )
        else:
            unmatched.remove(found)
    return differences


def _moments_equivalent(site_a, site_b, tol: float) -> bool:
    ma, mb = site_a.magmom, site_b.magmom
    zero_a = ma is None or ma.is_zero(tol)
    zero_b = mb is None or mb.is_zero(tol)
    if zero_a or zero_b:
        return zero_a and zero_b
    return ma.is_equivalent(mb, tol)


def structures_equivalent(
    a: NormalizedStructure, b: NormalizedStructure, **tolerances
) -> bool:
    return not compare_structures(a, b, **tolerances)
