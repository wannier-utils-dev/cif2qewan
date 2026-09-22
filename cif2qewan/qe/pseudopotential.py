"""The pseudopotential / projection table (``pp_psl_rrkj.csv`` and friends).

Columns: ``atom, pp_file_name, nexclude, orbitals, ecutwfc, ecutrho``. An
empty ``pp_file_name`` marks an element without a supported
pseudopotential; empty ``nexclude``/``ecutwfc``/``ecutrho`` count as 0 and
an empty ``orbitals`` as no projection.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

from cif2qewan.exceptions import PseudopotentialError

PathLike = Union[str, Path]

#: Wannier functions per projection letter in the ``orbitals`` column.
ORBITAL_SIZES = {"s": 1, "p": 3, "d": 5, "f": 7}


@dataclass(frozen=True)
class PseudopotentialEntry:
    """One row of the table for a supported element."""

    element: str
    file_name: str
    nexclude: int
    orbitals: str
    ecutwfc: float
    ecutrho: float

    def __post_init__(self) -> None:
        for letter in self.orbitals:
            if letter not in ORBITAL_SIZES:
                raise PseudopotentialError(
                    f"{self.element}: unknown orbital {letter!r} in table entry {self.orbitals!r}"
                )
        if self.nexclude < 0:
            raise PseudopotentialError(f"{self.element}: nexclude must not be negative")

    @property
    def num_wann(self) -> int:
        """Wannier functions per atom (without spin)."""
        return sum(ORBITAL_SIZES[letter] for letter in self.orbitals)

    @property
    def projection_orbitals(self) -> Tuple[str, ...]:
        return tuple(self.orbitals)

    def relativistic_file_name(self) -> str:
        """The fully relativistic file: ``.pbe`` -> ``.rel-pbe`` (PSLibrary), ``_sr`` -> ``_fr`` (PseudoDojo)."""
        return self.file_name.replace(".pbe", ".rel-pbe").replace("_sr", "_fr")

    def pseudo_file(self, relativistic: bool) -> str:
        return self.relativistic_file_name() if relativistic else self.file_name


class PseudopotentialTable:
    """Look up pseudopotential and projection information by element."""

    def __init__(
        self, entries: Dict[str, Optional[PseudopotentialEntry]], source: str = ""
    ):
        self._entries = dict(entries)
        self.source = source

    @classmethod
    def from_csv(cls, path: PathLike) -> "PseudopotentialTable":
        path = Path(path)
        if not path.is_file():
            raise PseudopotentialError(f"pseudopotential table not found: {path}")
        entries: Dict[str, Optional[PseudopotentialEntry]] = {}
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            required = {
                "atom",
                "pp_file_name",
                "nexclude",
                "orbitals",
                "ecutwfc",
                "ecutrho",
            }
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise PseudopotentialError(
                    f"{path}: expected columns {sorted(required)}, got {reader.fieldnames}"
                )
            for row in reader:
                element = (row["atom"] or "").strip()
                if not element:
                    continue
                file_name = (row["pp_file_name"] or "").strip()
                if not file_name:
                    entries[element] = None
                    continue
                try:
                    entries[element] = PseudopotentialEntry(
                        element=element,
                        file_name=f"{file_name}.UPF",
                        nexclude=int(float(row["nexclude"] or 0)),
                        orbitals=(row["orbitals"] or "").strip(),
                        ecutwfc=float(row["ecutwfc"] or 0.0),
                        ecutrho=float(row["ecutrho"] or 0.0),
                    )
                except ValueError as exc:
                    raise PseudopotentialError(
                        f"{path}: bad row for {element}: {exc}"
                    ) from exc
        return cls(entries, source=str(path))

    def __contains__(self, element: str) -> bool:
        return self._entries.get(element) is not None

    def lookup(self, element: str) -> PseudopotentialEntry:
        """The entry for ``element``; raises PseudopotentialError if unsupported."""
        if element not in self._entries:
            raise PseudopotentialError(
                f"{element} is not listed in the pseudopotential table {self.source}"
            )
        entry = self._entries[element]
        if entry is None:
            raise PseudopotentialError(
                f"{element} has no pseudopotential in the table {self.source}"
            )
        return entry
