"""k-point policies: SCF and NSCF meshes and the band-structure path."""

from __future__ import annotations

import itertools
import logging
import math
from typing import List, Sequence, Tuple

import numpy as np

from cif2qewan.qe.model import KPathPoint, KPointsList, KPointsPath
from cif2qewan.structure.model import NormalizedStructure, atomic_number
from cif2qewan.wannier90.model import KPathSegment

Mesh = Tuple[int, int, int]

logger = logging.getLogger(__name__)

#: Bounds of the NSCF (Wannier90) mesh per direction.
NSCF_MESH_MIN = 4
NSCF_MESH_MAX = 8

#: Divisions of the shortest band-path segment; longer segments scale with length.
BAND_PATH_BASE_DIVISIONS = 10

#: Divisions of every segment of the fallback path used without seekpath.
FALLBACK_PATH_DIVISIONS = 20

#: Tolerance on the segment-length ratio when truncating to an integer.
_RATIO_TOLERANCE = 1.0e-6


def scf_mesh(structure: NormalizedStructure, k_resolution: float) -> Mesh:
    """Monkhorst-Pack mesh with spacing at most ``k_resolution`` (1/angstrom).

    Uses cif2cell's rule: ``round(|b_i| / k_resolution)`` with the reciprocal
    vectors including 2*pi, and at least one point per direction.
    """
    lengths = np.linalg.norm(structure.reciprocal_lattice(), axis=1)
    return tuple(max(1, int(round(float(b) / k_resolution))) for b in lengths)  # type: ignore


def nscf_mesh(mesh: Sequence[int]) -> Mesh:
    """The SCF mesh clamped to [NSCF_MESH_MIN, NSCF_MESH_MAX] per direction."""
    return tuple(min(max(int(n), NSCF_MESH_MIN), NSCF_MESH_MAX) for n in mesh)  # type: ignore


def mesh_points(mesh: Sequence[int]) -> Tuple[Tuple[float, float, float], ...]:
    """The full Gamma-centred mesh in fractional coordinates, kz fastest."""
    n1, n2, n3 = (int(n) for n in mesh)
    return tuple(
        (i / n1, j / n2, k / n3)
        for i, j, k in itertools.product(range(n1), range(n2), range(n3))
    )


def mesh_kpoints_list(mesh: Sequence[int]) -> KPointsList:
    """``K_POINTS {crystal}`` listing the full mesh with equal weights."""
    points = mesh_points(mesh)
    weight = 1.0 / len(points)
    return KPointsList("crystal", tuple(point + (weight,) for point in points))


def fallback_band_path() -> KPointsPath:
    """R-G-X-M-G of a simple cubic cell, used when seekpath is unavailable."""
    vertices = [
        ("R", (0.5, 0.5, 0.5)),
        ("G", (0.0, 0.0, 0.0)),
        ("X", (0.5, 0.0, 0.0)),
        ("M", (0.5, 0.5, 0.0)),
        ("G", (0.0, 0.0, 0.0)),
    ]
    return KPointsPath(
        "crystal_b",
        tuple(
            KPathPoint(coords, FALLBACK_PATH_DIVISIONS, label)
            for label, coords in vertices
        ),
    )


def _seekpath():
    """The seekpath module if it is importable and recent enough, else None."""
    try:
        import seekpath
    except ImportError:
        return None
    version = tuple(int(x) for x in seekpath.__version__.split(".")[:2] if x.isdigit())
    if version < (2, 1):
        return None
    return seekpath


def band_path(structure: NormalizedStructure) -> Tuple[KPointsPath, bool]:
    """Band-structure path through the high-symmetry points of ``structure``.

    Returns the path and a flag telling whether seekpath was used. Segment
    divisions are ``BAND_PATH_BASE_DIVISIONS`` times the segment length
    relative to the shortest segment, truncated to an integer (with a small
    tolerance against floating-point noise); a vertex with 0 divisions marks
    a jump in the path. The last vertex carries the base count, as in 0.2.x.
    """
    seekpath = _seekpath()
    if seekpath is None:
        logger.warning(
            "seekpath >= 2.1 is not available; using the simple-cubic R-G-X-M-G band path"
        )
        return fallback_band_path(), False

    cell = structure.lattice_matrix.tolist()
    positions = structure.frac_coords().tolist()
    numbers = [atomic_number(site.element) for site in structure.sites]
    kpath = seekpath.getpaths.get_explicit_k_path_orig_cell([cell, positions, numbers])
    return _path_from_seekpath(kpath), True


def _path_from_seekpath(kpath) -> KPointsPath:
    labels: List[str] = list(kpath["explicit_kpoints_labels"])
    rel = np.asarray(kpath["explicit_kpoints_rel"], dtype=float)
    absolute = np.asarray(kpath["explicit_kpoints_abs"], dtype=float)

    vertices: List[Tuple[str, np.ndarray, bool]] = []  # (label, k, jump before next)
    indices = [i for i, label in enumerate(labels) if label != ""]
    for n, i in enumerate(indices):
        label = labels[i].replace("GAMMA", "G").replace("SIGMA", "S")
        # two labelled points in a row mark a discontinuity: the previous
        # vertex ends a line and this one starts the next
        if n > 0 and indices[n - 1] == i - 1:
            vertices[-1] = (vertices[-1][0], vertices[-1][1], True)
        vertices.append((label, rel[i], False))

    distances = []
    for n in range(len(indices) - 1):
        if vertices[n][2]:
            distances.append(0.0)
        else:
            distances.append(
                float(np.linalg.norm(absolute[indices[n + 1]] - absolute[indices[n]]))
            )
    nonzero = [d for d in distances if d > 0.0]
    shortest = min(nonzero) if nonzero else 1.0

    points = []
    for n, (label, k, _) in enumerate(vertices):
        if n < len(distances):
            ratio = BAND_PATH_BASE_DIVISIONS * distances[n] / shortest
            ndiv = int(math.floor(ratio + _RATIO_TOLERANCE))
        else:
            ndiv = BAND_PATH_BASE_DIVISIONS
        points.append(KPathPoint(tuple(float(x) for x in k), ndiv, label))
    return KPointsPath("crystal_b", tuple(points))


def win_kpoint_path(path: KPointsPath) -> Tuple[KPathSegment, ...]:
    """The Wannier90 ``kpoint_path`` segments of a QE band path (no jumps)."""
    return tuple(
        KPathSegment(a.label, a.coords, b.label, b.coords) for a, b in path.segments()
    )


def atomic_mass(element: str) -> float:
    """Standard atomic mass in amu (from pymatgen)."""
    from pymatgen.core.periodic_table import Element

    return float(Element(element).atomic_mass)


__all__ = [
    "Mesh",
    "atomic_mass",
    "band_path",
    "fallback_band_path",
    "mesh_kpoints_list",
    "mesh_points",
    "nscf_mesh",
    "scf_mesh",
    "win_kpoint_path",
]
