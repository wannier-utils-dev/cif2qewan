"""The set of inputs that cif2qewan generates for one structure."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Union

from cif2qewan.exceptions import InputModelError
from cif2qewan.qe.model import (
    KPointsAutomatic,
    KPointsList,
    KPointsPath,
    NamelistInput,
    PwInput,
)
from cif2qewan.structure.model import NormalizedStructure
from cif2qewan.wannier90.model import Wannier90Input

RenderableInput = Union[PwInput, NamelistInput, Wannier90Input]


@dataclass
class CalculationPlan:
    """Everything needed to run the QE + Wannier90 workflow for one structure.

    The plan is produced by :mod:`cif2qewan.workflow.builder`
    and consumed by the writers; it does not touch the file system itself.
    ``files()`` gives the intended job layout.
    """

    structure: NormalizedStructure
    scf: PwInput
    nscf: PwInput
    check_wannier: PwInput
    bands_nscf: PwInput
    bands: NamelistInput
    pw2wan: NamelistInput
    projwfc: NamelistInput
    pp: NamelistInput
    wannier90: Wannier90Input

    def __post_init__(self) -> None:
        if not isinstance(self.structure, NormalizedStructure):
            raise InputModelError("CalculationPlan needs a NormalizedStructure")
        if self.scf.calculation != "scf":
            raise InputModelError("scf input must have calculation = 'scf'")
        for name in ("nscf", "check_wannier"):
            if getattr(self, name).calculation != "nscf":
                raise InputModelError(f"{name} input must have calculation = 'nscf'")
        if self.bands_nscf.calculation != "bands":
            raise InputModelError("bands_nscf input must have calculation = 'bands'")
        if not isinstance(self.bands_nscf.kpoints, KPointsPath):
            raise InputModelError("bands_nscf input must use a K_POINTS path")

        mesh = self.wannier90.mp_grid
        nscf_k = self.nscf.kpoints
        if not isinstance(nscf_k, KPointsList) or nscf_k.num_points != len(
            self.wannier90.kpoints
        ):
            raise InputModelError(
                "nscf must list the same k points as the Wannier90 mp_grid"
            )
        check_k = self.check_wannier.kpoints
        if not isinstance(check_k, KPointsAutomatic) or check_k.mesh != mesh:
            raise InputModelError("check_wannier must use the Wannier90 mp_grid")

        nat = self.structure.num_sites
        for name in ("scf", "nscf", "check_wannier", "bands_nscf"):
            if getattr(self, name).nat != nat:
                raise InputModelError(
                    f"{name} has {getattr(self, name).nat} atoms, structure has {nat}"
                )
        if len(self.wannier90.atoms_frac) != nat:
            raise InputModelError("Wannier90 atoms_frac does not match the structure")

    def files(self) -> Dict[str, RenderableInput]:
        """Relative output path -> input model, in the order they are written."""
        return {
            "scf.in": self.scf,
            "nscf.in": self.nscf,
            "pw2wan.in": self.pw2wan,
            "pwscf.win": self.wannier90,
            "check_wannier/nscf.in": self.check_wannier,
            "band/nscf.in": self.bands_nscf,
            "band/band.in": self.bands,
            "band/proj.in": self.projwfc,
            "band/pp.in": self.pp,
        }
