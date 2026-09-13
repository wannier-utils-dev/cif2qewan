"""Build the :class:`CalculationPlan` for one structure.

This module holds the scientific policy of cif2qewan: which pseudopotential
and cutoffs to use, how many bands and Wannier functions, which k meshes,
and how the ``--so`` / ``--mag`` options translate into QE settings. The
numbers reproduce the 0.2.x generator; the conventions are listed in
CLAUDE.md under "Architecture notes" and must not change silently.

Magnetism policy (``--mag``): a ferromagnetic starting guess of
``STARTING_MAGNETIZATION`` on every species. The SCF is run collinear
(``nspin = 2``, scalar-relativistic pseudopotentials); from the NSCF on the
run is noncollinear with ``lforcet = .true.`` and the moment along z, with
the fully relativistic pseudopotentials when ``--so`` is also given.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Tuple

import numpy as np

from cif2qewan.config import Config
from cif2qewan.exceptions import InputModelError
from cif2qewan.qe.kpoints import (
    Mesh,
    atomic_mass,
    band_path,
    mesh_kpoints_list,
    mesh_points,
    nscf_mesh,
    scf_mesh,
    win_kpoint_path,
)
from cif2qewan.qe.model import (
    AtomicPosition,
    CellParameters,
    KPoints,
    KPointsAutomatic,
    Namelist,
    NamelistInput,
    PwInput,
    RawValue,
    Species,
)
from cif2qewan.qe.pseudopotential import PseudopotentialEntry, PseudopotentialTable
from cif2qewan.structure.model import NormalizedStructure
from cif2qewan.structure.readers.cif2cell import Cif2cellOutput
from cif2qewan.wannier90.model import AtomFrac, Projection, Wannier90Input
from cif2qewan.workflow.model import CalculationPlan

logger = logging.getLogger(__name__)

PREFIX = "pwscf"
OUTDIR = "./work"
SEEDNAME = "pwscf"

#: Initial magnetization per species for ``--mag`` (QE starting_magnetization).
STARTING_MAGNETIZATION = 3.0

#: nbnd = nexclude + NSCF_BANDS_PER_WANN * num_wann for the Wannier90 NSCF run.
NSCF_BANDS_PER_WANN = 3
#: nbnd = nexclude + int(CHECK_BANDS_PER_WANN * num_wann) for the check_wannier run.
CHECK_BANDS_PER_WANN = 1.5

CONV_THR_SCF = "1.0d-8"
CONV_THR_NSCF = "1.d-10"
CONV_THR_CHECK = "1.d-8"

WANNIER90_PARAMETERS = {
    "dis_num_iter": 200,
    "num_iter": 0,
    "dis_froz_max": -200,
    "dis_froz_min": -200,
    "bands_plot": True,
    "write_hr": True,
    "write_tb": True,
    "fermi_surface_plot": True,
    "wannier_plot": True,
}


@dataclass(frozen=True)
class StructureHints:
    """Values taken from the structure source instead of being recomputed.

    ``alat`` and ``masses`` fix the QE representation of the cell and the
    ATOMIC_SPECIES rows; ``scf_kmesh`` overrides the mesh derived from the
    k resolution. The cif2cell reader supplies all three so that the inputs
    match the 0.2.x generator exactly; the pymatgen reader supplies none.
    """

    alat: Optional[float] = None
    masses: Optional[Mapping[str, float]] = None
    scf_kmesh: Optional[Mesh] = None

    @classmethod
    def from_cif2cell(cls, output: Cif2cellOutput) -> "StructureHints":
        return cls(alat=output.alat, masses=output.masses, scf_kmesh=output.kmesh)


@dataclass(frozen=True)
class WannierCounts:
    """Band and Wannier-function counts derived from the pseudopotential table."""

    num_wann: int  # per spin channel
    nexclude: int  # per spin channel
    spinor: bool

    @property
    def factor(self) -> int:
        return 2 if self.spinor else 1

    @property
    def num_wann_total(self) -> int:
        return self.num_wann * self.factor

    @property
    def nexclude_total(self) -> int:
        return self.nexclude * self.factor

    @property
    def num_bands(self) -> int:
        return NSCF_BANDS_PER_WANN * self.num_wann_total

    @property
    def nbnd_nscf(self) -> int:
        return (self.nexclude + NSCF_BANDS_PER_WANN * self.num_wann) * self.factor

    @property
    def nbnd_check(self) -> int:
        return (self.nexclude + int(CHECK_BANDS_PER_WANN * self.num_wann)) * self.factor


class WorkflowBuilder:
    """Turn a structure and a configuration into a CalculationPlan."""

    def __init__(
        self, config: Config, table: Optional[PseudopotentialTable] = None
    ) -> None:
        self.config = config
        self.table = (
            table
            if table is not None
            else PseudopotentialTable.from_csv(config.pp_list_path)
        )

    # -- public -----------------------------------------------------------

    def build(
        self, structure: NormalizedStructure, hints: Optional[StructureHints] = None
    ) -> CalculationPlan:
        hints = hints or StructureHints()
        config = self.config

        labels = self._species_labels(structure, hints)
        entries = {
            label: self.table.lookup(self._element_of(structure, label))
            for label in labels
        }
        masses = self._masses(structure, labels, hints)
        ecutwfc, ecutrho = self._cutoffs(entries.values())
        counts = self._counts(structure, entries)

        alat = (
            hints.alat
            if hints.alat is not None
            else float(np.linalg.norm(structure.lattice_matrix[0]))
        )
        cell = CellParameters("alat", structure.lattice_matrix / alat)
        positions = tuple(
            AtomicPosition(site.label, site.frac_coords) for site in structure.sites
        )

        kmesh_scf = (
            hints.scf_kmesh
            if hints.scf_kmesh is not None
            else scf_mesh(structure, config.scf_k_resolution)
        )
        kmesh_nscf = nscf_mesh(kmesh_scf)
        path, _ = band_path(structure)
        logger.info(
            "species %s; ecutwfc %g Ry, ecutrho %g Ry; num_wann %d, nexclude %d (per spin); "
            "scf mesh %s, nscf mesh %s, %d band-path vertices",
            ", ".join(f"{label}:{entries[label].file_name}" for label in labels),
            ecutwfc,
            ecutrho,
            counts.num_wann,
            counts.nexclude,
            "x".join(map(str, kmesh_scf)),
            "x".join(map(str, kmesh_nscf)),
            path.num_points,
        )

        def species(relativistic: bool) -> Tuple[Species, ...]:
            return tuple(
                Species(label, masses[label], entries[label].pseudo_file(relativistic))
                for label in labels
            )

        def pw_input(
            control, system, conv_thr, relativistic, kpoints: KPoints
        ) -> PwInput:
            return PwInput(
                control=control,
                system=system,
                electrons=self._electrons(conv_thr),
                species=species(relativistic),
                positions=positions,
                kpoints=kpoints,
                cell=cell,
            )

        base_system = self._base_system(
            alat, structure.num_sites, len(labels), ecutwfc, ecutrho
        )

        scf = pw_input(
            self._control("scf"),
            self._system(base_system, first={}, spin=self._scf_spin(len(labels))),
            CONV_THR_SCF,
            relativistic=config.so and not config.mag,
            kpoints=KPointsAutomatic(kmesh_scf),
        )
        nscf = pw_input(
            self._control("nscf"),
            self._system(
                base_system,
                first={"nosym": True, "nbnd": counts.nbnd_nscf},
                spin=self._nscf_spin(len(labels)),
            ),
            CONV_THR_NSCF,
            relativistic=config.so,
            kpoints=mesh_kpoints_list(kmesh_nscf),
        )
        check_wannier = pw_input(
            self._control("nscf", verbosity="high"),
            self._system(
                base_system,
                first={"nbnd": counts.nbnd_check},
                spin=self._nscf_spin(len(labels)),
            ),
            CONV_THR_CHECK,
            relativistic=config.so,
            kpoints=KPointsAutomatic(kmesh_nscf, (1, 1, 1)),
        )
        # The bands run inherits the check_wannier settings (0.2.x behaviour).
        bands_nscf = pw_input(
            self._control("bands", verbosity="high"),
            self._system(
                base_system,
                first={"nbnd": counts.nbnd_check},
                spin=self._nscf_spin(len(labels)),
            ),
            CONV_THR_CHECK,
            relativistic=config.so,
            kpoints=path,
        )

        wannier90 = Wannier90Input(
            num_wann=counts.num_wann_total,
            num_bands=counts.num_bands,
            unit_cell_cart=structure.lattice_matrix,
            atoms_frac=tuple(
                AtomFrac(site.label, site.frac_coords) for site in structure.sites
            ),
            mp_grid=kmesh_nscf,
            kpoints=mesh_points(kmesh_nscf),
            projections=tuple(
                Projection(label, entries[label].projection_orbitals)
                for label in labels
                if entries[label].orbitals
            ),
            spinors=counts.spinor,
            exclude_bands=(1, counts.nexclude_total) if counts.nexclude > 0 else None,
            kpoint_path=win_kpoint_path(path),
            parameters=dict(WANNIER90_PARAMETERS),
        )

        return CalculationPlan(
            structure=structure,
            scf=scf,
            nscf=nscf,
            check_wannier=check_wannier,
            bands_nscf=bands_nscf,
            bands=self._bands_input(),
            pw2wan=self._pw2wan_input(),
            projwfc=self._projwfc_input(),
            pp=self._pp_input(counts),
            wannier90=wannier90,
        )

    # -- species, cutoffs, counts -----------------------------------------

    @staticmethod
    def _species_labels(
        structure: NormalizedStructure, hints: StructureHints
    ) -> Tuple[str, ...]:
        labels = structure.species_labels
        if hints.masses is not None:
            listed = tuple(hints.masses)
            if set(listed) != set(labels):
                raise InputModelError(
                    f"species {sorted(listed)} from the structure source do not match the sites {sorted(labels)}"
                )
            return listed  # keep the source's ATOMIC_SPECIES order
        return labels

    @staticmethod
    def _element_of(structure: NormalizedStructure, label: str) -> str:
        for site in structure.sites:
            if site.label == label:
                return site.element
        raise InputModelError(f"no site with label {label!r}")

    def _masses(self, structure, labels, hints) -> Dict[str, float]:
        if hints.masses is not None:
            return {label: float(hints.masses[label]) for label in labels}
        return {
            label: atomic_mass(self._element_of(structure, label)) for label in labels
        }

    @staticmethod
    def _cutoffs(entries) -> Tuple[float, float]:
        ecutwfc = max(entry.ecutwfc for entry in entries)
        ecutrho = max(entry.ecutrho for entry in entries)
        if ecutrho < 4 * ecutwfc:
            warnings.warn(
                f"ecut_rho should be bigger than 4*ecut_wfc, but {ecutrho} < {4 * ecutwfc}."
            )
        return ecutwfc, ecutrho

    def _counts(
        self,
        structure: NormalizedStructure,
        entries: Mapping[str, PseudopotentialEntry],
    ) -> WannierCounts:
        num_wann = sum(entries[site.label].num_wann for site in structure.sites)
        nexclude = sum(entries[site.label].nexclude for site in structure.sites)
        if num_wann == 0:
            raise InputModelError(
                "no Wannier projections: every species has an empty orbitals entry"
            )
        return WannierCounts(num_wann, nexclude, self.config.spinor)

    # -- namelists ----------------------------------------------------------

    def _control(self, calculation: str, **extra) -> Namelist:
        entries = {
            "calculation": calculation,
            "restart_mode": "from_scratch",
            "prefix": PREFIX,
            "tstress": True,
            "tprnfor": True,
            "pseudo_dir": self.config.pseudo_dir,
            "outdir": OUTDIR,
            "wf_collect": True,
            "disk_io": "low",
        }
        entries.update(extra)
        return Namelist("control", entries)

    def _base_system(self, alat, nat, ntyp, ecutwfc, ecutrho) -> Dict:
        return {
            "ibrav": 0,
            "A": RawValue(f"{alat:10.5f}"),
            "nat": nat,
            "ntyp": ntyp,
            "ecutwfc": ecutwfc,
            "ecutrho": ecutrho,
            "occupations": "smearing",
            "smearing": "m-p",
            "degauss": RawValue(f"{self.config.degauss:.3f}"),
        }

    @staticmethod
    def _system(base: Mapping, first: Mapping, spin: Mapping) -> Namelist:
        entries = dict(first)
        entries.update(base)
        entries.update(spin)
        return Namelist("system", entries)

    def _scf_spin(self, ntyp: int) -> Dict:
        if self.config.mag:
            spin: Dict = {"nspin": 2}
            spin.update(self._starting_magnetization(ntyp))
            return spin
        if self.config.so:
            return {"lspinorb": True, "noncolin": True}
        return {}

    def _nscf_spin(self, ntyp: int) -> Dict:
        if self.config.mag:
            spin: Dict = {
                "lspinorb": bool(self.config.so),
                "noncolin": True,
                "lforcet": True,
                "angle1": 0,
                "angle2": 0,
            }
            spin.update(self._starting_magnetization(ntyp))
            return spin
        if self.config.so:
            return {"lspinorb": True, "noncolin": True}
        return {}

    @staticmethod
    def _starting_magnetization(ntyp: int) -> Dict:
        return {
            f"starting_magnetization({i + 1})": STARTING_MAGNETIZATION
            for i in range(ntyp)
        }

    @staticmethod
    def _electrons(conv_thr: str) -> Namelist:
        return Namelist(
            "electrons",
            {
                "mixing_mode": "plain",
                "mixing_beta": 0.1,
                "conv_thr": RawValue(conv_thr),
            },
        )

    # -- auxiliary inputs -----------------------------------------------------

    @staticmethod
    def _bands_input() -> NamelistInput:
        return NamelistInput(
            (
                Namelist(
                    "bands",
                    {"prefix": PREFIX, "outdir": OUTDIR + "/", "filband": "bands.out"},
                ),
            )
        )

    def _pw2wan_input(self) -> NamelistInput:
        options = self.config.pw2wan
        entries: Dict = {
            "outdir": OUTDIR,
            "prefix": PREFIX,
            "seedname": SEEDNAME,
            "spin_component": "none",
        }
        if "wannier_plot_supercell" in options:
            entries["wannier_plot_supercell"] = RawValue(
                str(options["wannier_plot_supercell"])
            )
        entries["write_mmn"] = True
        entries["write_amn"] = True
        entries["write_unk"] = RawValue(str(options["write_unk"]))
        return NamelistInput((Namelist("inputpp", entries),))

    def _projwfc_input(self) -> NamelistInput:
        return NamelistInput(
            (
                Namelist(
                    "projwfc",
                    {
                        "prefix": PREFIX,
                        "outdir": OUTDIR,
                        "kresolveddos": False,
                        "degauss": RawValue(f"{self.config.degauss:.3f}"),
                        "Emax": RawValue(""),
                        "Emin": RawValue(""),
                    },
                ),
            )
        )

    @staticmethod
    def _pp_input(counts: WannierCounts) -> NamelistInput:
        return NamelistInput(
            (
                Namelist(
                    "inputpp",
                    {
                        "prefix": PREFIX,
                        "outdir": OUTDIR,
                        "filplot": "wf_pp",
                        "plot_num": 7,
                        "kpoint": 1,
                        "kband(1)": counts.nexclude_total + 1,
                        "kband(2)": counts.nexclude_total + counts.num_wann_total,
                        "lsign": RawValue(".TRUE."),
                    },
                ),
                Namelist("plot", {"iflag": 3, "output_format": 5, "fileout": ".xsf"}),
            )
        )


def build_plan(
    structure: NormalizedStructure,
    config: Config,
    hints: Optional[StructureHints] = None,
    table: Optional[PseudopotentialTable] = None,
) -> CalculationPlan:
    """Convenience wrapper around :class:`WorkflowBuilder`."""
    return WorkflowBuilder(config, table).build(structure, hints)


def plan_from_cif2cell_output(
    output: Cif2cellOutput, config: Config
) -> CalculationPlan:
    """The plan for a cif2cell output, reproducing the 0.2.x inputs."""
    return build_plan(output.structure, config, StructureHints.from_cif2cell(output))


__all__ = [
    "STARTING_MAGNETIZATION",
    "StructureHints",
    "WannierCounts",
    "WorkflowBuilder",
    "build_plan",
    "plan_from_cif2cell_output",
]
