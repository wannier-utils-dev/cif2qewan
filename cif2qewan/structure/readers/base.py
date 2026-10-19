"""The reader boundary: anything with ``read(path) -> NormalizedStructure``."""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, Union, runtime_checkable

from cif2qewan.structure.model import NormalizedStructure

PathLike = Union[str, Path]


@runtime_checkable
class StructureReader(Protocol):
    """Reads one structure file and returns it in the units of this package.

    A reader must not write Quantum ESPRESSO input files, change the current
    working directory, or leave files next to the input.
    """

    def read(self, path: PathLike) -> NormalizedStructure:
        """Return the structure stored in ``path``; raise StructureError on failure."""
