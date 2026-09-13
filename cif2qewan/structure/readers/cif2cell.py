"""Run cif2cell, or parse its output, to obtain a structure.

This is the compatibility path of the 0.2.x generator: cif2cell writes a
``pw.x`` input (``cif_scf.in``) from which the cell, the atoms and the
automatic k mesh are taken. The command is executed in a fresh temporary
directory with an argument list, and any failure raises
:class:`ExternalCommandError` instead of falling back to an old file.
"""

from __future__ import annotations

import logging
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np

from cif2qewan.exceptions import ExternalCommandError, StructureError
from cif2qewan.structure.model import AtomicSite, NormalizedStructure
from cif2qewan.structure.readers.base import PathLike

BOHR_IN_ANGSTROM = 0.529177210903

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Cif2cellOutput:
    """The parts of a cif2cell ``pw.x`` input that cif2qewan uses.

    Attributes
    ----------
    structure : NormalizedStructure
        Cell in angstrom and fractional positions.
    alat : float
        The lattice parameter ``A`` in angstrom.
    cell_alat : 3x3 tuple
        CELL_PARAMETERS in units of ``alat``, exactly as written by cif2cell.
    masses : dict
        Species label -> atomic mass from ATOMIC_SPECIES.
    kmesh : (int, int, int) or None
        The automatic k mesh cif2cell chose for the requested resolution.
    text : str
        The complete file, for callers that need other details.
    """

    structure: NormalizedStructure
    alat: float
    cell_alat: Tuple[Tuple[float, float, float], ...]
    masses: Dict[str, float]
    kmesh: Optional[Tuple[int, int, int]]
    text: str


def _card(lines: Sequence[str], name: str) -> Tuple[str, list]:
    """(option, body lines) of the QE card ``name``; body ends at a blank or a new card."""
    for i, line in enumerate(lines):
        if line.strip().upper().startswith(name):
            option = line.strip()[len(name) :].strip(" {}()").lower()
            body = []
            for entry in lines[i + 1 :]:
                if not entry.strip() or entry.lstrip().startswith("#"):
                    break
                if re.match(r"^\s*[A-Z][A-Z_]+(\s|$|\{|\()", entry):
                    break  # next card (CELL_PARAMETERS, ATOMIC_SPECIES, ...)
                body.append(entry)
            return option, body
    raise StructureError(f"cif2cell output has no {name} card")


def parse_cif2cell_output(text: str) -> Cif2cellOutput:
    """Parse the ``pw.x`` input written by ``cif2cell -p pwscf``."""
    lines = text.splitlines()

    def namelist_value(key: str) -> Optional[str]:
        pattern = re.compile(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", re.I)
        for line in lines:
            match = pattern.match(line)
            if match:
                return match.group(1).rstrip(",")
        return None

    nat = namelist_value("nat")
    ntyp = namelist_value("ntyp")
    if nat is None or ntyp is None:
        raise StructureError(
            "cif2cell output lacks nat/ntyp; it is not a usable pw.x input"
        )
    ibrav = namelist_value("ibrav")
    if ibrav is not None and int(ibrav) != 0:
        raise StructureError(
            f"cif2cell output uses ibrav = {ibrav}; only ibrav = 0 is supported"
        )

    option, rows = _card(lines, "CELL_PARAMETERS")
    try:
        cell = np.array(
            [[float(x) for x in row.split()[:3]] for row in rows], dtype=float
        )
    except ValueError:
        raise StructureError(
            "cannot parse CELL_PARAMETERS in cif2cell output"
        ) from None
    if cell.shape != (3, 3):
        raise StructureError("CELL_PARAMETERS in cif2cell output must have three rows")

    alat_value = namelist_value("A")
    if option == "alat":
        if alat_value is None:
            raise StructureError("CELL_PARAMETERS {alat} without A in cif2cell output")
        alat = float(alat_value)
        lattice = cell * alat
    elif option == "angstrom":
        alat = (
            float(alat_value)
            if alat_value is not None
            else float(np.linalg.norm(cell[0]))
        )
        lattice = cell
    elif option == "bohr":
        lattice = cell * BOHR_IN_ANGSTROM
        alat = (
            float(alat_value)
            if alat_value is not None
            else float(np.linalg.norm(lattice[0]))
        )
    else:
        raise StructureError(
            f"unsupported CELL_PARAMETERS option {option!r} in cif2cell output"
        )

    _, rows = _card(lines, "ATOMIC_SPECIES")
    masses: Dict[str, float] = {}
    for row in rows:
        fields = row.split()
        if len(fields) < 2:
            raise StructureError(f"cannot parse ATOMIC_SPECIES row {row!r}")
        masses[fields[0]] = float(fields[1])
    if len(masses) != int(ntyp):
        raise StructureError(f"ntyp = {ntyp} but {len(masses)} ATOMIC_SPECIES rows")

    option, rows = _card(lines, "ATOMIC_POSITIONS")
    if option != "crystal":
        raise StructureError(
            f"ATOMIC_POSITIONS {option!r} in cif2cell output; only crystal is supported"
        )
    sites = []
    for row in rows:
        fields = row.split()
        if len(fields) < 4:
            raise StructureError(f"cannot parse ATOMIC_POSITIONS row {row!r}")
        label = fields[0]
        if label not in masses:
            raise StructureError(
                f"position label {label!r} has no ATOMIC_SPECIES entry"
            )
        element = re.match(r"[A-Z][a-z]?", label)
        if element is None:
            raise StructureError(f"cannot derive an element from label {label!r}")
        sites.append(
            AtomicSite(
                element=element.group(0),
                frac_coords=tuple(float(x) for x in fields[1:4]),
                label=label,
            )
        )
    if len(sites) != int(nat):
        raise StructureError(f"nat = {nat} but {len(sites)} ATOMIC_POSITIONS rows")

    kmesh = None
    try:
        option, rows = _card(lines, "K_POINTS")
    except StructureError:
        option, rows = "", []
    if option == "automatic" and rows:
        kmesh = tuple(int(x) for x in rows[0].split()[:3])  # type: ignore[assignment]

    structure = NormalizedStructure(lattice, sites)
    return Cif2cellOutput(
        structure=structure,
        alat=alat,
        cell_alat=tuple(tuple(float(x) for x in row) for row in lattice / alat),
        masses=masses,
        kmesh=kmesh,
        text=text,
    )


class Cif2cellReader:
    """Obtain a structure through cif2cell.

    Parameters
    ----------
    cif2cell_path : str or Path
        The cif2cell executable.
    k_resolution : float
        ``--k-resolution`` in 1/angstrom; also determines ``kmesh``.
    print_digits : int
        ``--print-digits`` for the coordinates.
    timeout : float, optional
        Seconds to wait for cif2cell.
    """

    output_name = "cif_scf.in"

    def __init__(
        self,
        cif2cell_path: PathLike,
        k_resolution: float = 0.15,
        print_digits: int = 10,
        timeout: Optional[float] = 600.0,
    ) -> None:
        self.cif2cell_path = str(cif2cell_path)
        self.k_resolution = float(k_resolution)
        self.print_digits = int(print_digits)
        self.timeout = timeout

    def command(self, cif_path: Path, output: Path) -> list:
        return [
            self.cif2cell_path,
            "-p",
            "pwscf",
            "--setup-all",
            f"--k-resolution={self.k_resolution:.3f}",
            f"--print-digits={self.print_digits}",
            "-o",
            str(output),
            str(cif_path),
        ]

    def run(self, cif_path: PathLike) -> Cif2cellOutput:
        """Run cif2cell on ``cif_path`` in a temporary directory and parse the result."""
        cif_path = Path(cif_path).resolve()
        if not cif_path.is_file():
            raise StructureError(f"CIF file not found: {cif_path}")
        with tempfile.TemporaryDirectory(prefix="cif2qewan-cif2cell-") as tmp:
            output = Path(tmp) / self.output_name
            command = self.command(cif_path, output)
            logger.info("running %s", " ".join(command))
            try:
                completed = subprocess.run(
                    command,
                    cwd=tmp,
                    capture_output=True,
                    text=True,
                    timeout=self.timeout,
                )
            except OSError as exc:
                raise ExternalCommandError(
                    command,
                    message=f"cannot run cif2cell {self.cif2cell_path!r}: {exc}",
                ) from exc
            except subprocess.TimeoutExpired as exc:
                raise ExternalCommandError(
                    command, message=f"cif2cell did not finish within {self.timeout} s"
                ) from exc
            if completed.stderr:
                logger.debug("cif2cell stderr:\n%s", completed.stderr.rstrip())
            if completed.returncode != 0:
                raise ExternalCommandError(
                    command, completed.returncode, completed.stderr
                )
            if not output.is_file() or output.stat().st_size == 0:
                raise ExternalCommandError(
                    command,
                    completed.returncode,
                    completed.stderr,
                    message=f"cif2cell exited normally but wrote no {self.output_name}",
                )
            text = output.read_text()
        return parse_cif2cell_output(text)

    def read(self, path: PathLike) -> NormalizedStructure:
        return self.run(path).structure

    @staticmethod
    def read_output(path: PathLike) -> Cif2cellOutput:
        """Parse an existing cif2cell output file given explicitly by the caller."""
        path = Path(path)
        if not path.is_file():
            raise StructureError(f"cif2cell output file not found: {path}")
        return parse_cif2cell_output(path.read_text())
