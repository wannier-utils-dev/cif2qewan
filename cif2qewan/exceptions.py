"""Exceptions raised by cif2qewan.

Every error that cif2qewan raises on purpose derives from
:class:`Cif2qewanError`, so that callers can catch failures of the tool
without catching programming errors.
"""

from __future__ import annotations

from typing import Optional, Sequence


class Cif2qewanError(Exception):
    """Base class for errors raised by cif2qewan."""


class ConfigError(Cif2qewanError):
    """A configuration value (TOML file, command-line option) is missing or invalid."""


class StructureError(Cif2qewanError):
    """A crystal structure is malformed or physically inconsistent."""


class InputModelError(Cif2qewanError):
    """A Quantum ESPRESSO or Wannier90 input model is internally inconsistent."""


class PseudopotentialError(Cif2qewanError):
    """No usable pseudopotential or projection information for an element."""


class DataFileError(Cif2qewanError):
    """A QE or Wannier90 output file (scf.out, *_hr.dat, ...) is missing or malformed."""


class ExternalCommandError(Cif2qewanError):
    """An external program (cif2cell, ...) failed or produced no output.

    Parameters
    ----------
    command : sequence of str
        The argument list that was executed.
    returncode : int, optional
        Exit status of the program, if it ran.
    stderr : str, optional
        Captured standard error of the program.
    message : str, optional
        Additional explanation; by default one is built from the other fields.
    """

    def __init__(
        self,
        command: Sequence[str],
        returncode: Optional[int] = None,
        stderr: Optional[str] = None,
        message: Optional[str] = None,
    ) -> None:
        self.command = tuple(command)
        self.returncode = returncode
        self.stderr = stderr
        if message is None:
            message = f"command {' '.join(self.command)!r} failed"
            if returncode is not None:
                message += f" with exit status {returncode}"
            if stderr:
                message += f":\n{stderr.rstrip()}"
        super().__init__(message)
