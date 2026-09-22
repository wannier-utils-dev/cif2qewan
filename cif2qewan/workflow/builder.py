"""Build the :class:`CalculationPlan` for one structure.

This module holds the scientific policy of cif2qewan: which pseudopotential
and cutoffs to use, how many bands and Wannier functions, which k meshes,
and how the ``--so`` / ``--mag`` options and the magnetic moments of the
structure translate into QE settings. The numbers reproduce the 0.2.x
generator and must not change silently; changing any of them changes the
results of published workflows:

- SCF mesh: cif2cell's mesh, or ``round(|b_i| / scf_k_resolution)`` (at
  least 1) for the pymatgen reader; NSCF mesh: each dimension clamped to
  ``[4, 8]``.
- ``nbnd``: ``nexclude + 3 * num_wann`` (nscf), ``nexclude + int(1.5 *
  num_wann)`` (check_wannier and bands), doubled for spinors.
- ``conv_thr``: ``1.0d-8`` (scf), ``1.d-10`` (nscf), ``1.d-8`` (check_wannier,
  bands); ``ecutwfc`` / ``ecutrho``: maximum over the species from the table.
- Wannier90: ``num_bands = 3 * num_wann``, ``exclude_bands = 1-nexclude``,
  ``dis_num_iter = 200``, ``num_iter = 0``, ``dis_froz_max = dis_froz_min = -200``.

Magnetism policy:

- ``--mag`` without moments in the structure: a ferromagnetic starting guess
  of ``STARTING_MAGNETIZATION`` on every species. The SCF is run collinear
  (``nspin = 2``, scalar-relativistic pseudopotentials); from the NSCF on
  the run is noncollinear with ``lforcet = .true.`` and the moment along z,
  with the fully relativistic pseudopotentials when ``--so`` is also given.
- Moments in the structure (MCIF): species are split by moment
  (:mod:`cif2qewan.structure.magnetism`) and ``starting_magnetization`` is
  the moment relative to the largest one. A collinear order follows the
  two-step scheme above with the sign and the common axis of the moments; a
  noncollinear order is run noncollinear from the SCF on with ``angle1`` /
  ``angle2`` per species. Either way the Wannier functions are spinors.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Tuple

import numpy as np

from cif2qewan.config import Config
from cif2qewan.exceptions import InputModelError, StructureError
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
from cif2qewan.structure.magnetism import (
    COLLINEAR,
    NONMAGNETIC,
    QEMagnetization,
    classify_magnetic_sites,
    magnetic_order,
    qe_magnetization,
)
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


class SpinPolicy:
    """The &system spin settings of the SCF and NSCF runs.

    Parameters
    ----------
    order : str
        Magnetic order of the structure (``nonmagnetic`` when it has no moments).
    magnetization : mapping
        Species label -> :class:`QEMagnetization` from the structure's moments.
    so, mag : bool
        The command-line options.
    """

    def __init__(
        self,
        order: str,
        magnetization: Mapping[str, QEMagnetization],
        so: bool,
        mag: bool,
    ) -> None:
        self.order = order
        self.magnetization = dict(magnetization)
        self.so = bool(so)
        self.mag = bool(mag)

    @property
    def from_structure(self) -> bool:
        return self.order != NONMAGNETIC

    @property
    def spinor(self) -> bool:
        return self.so or self.mag or self.from_structure

    @property
    def collinear_scf(self) -> bool:
        """The SCF is a collinear nspin = 2 run (ferromagnetic guess or collinear order)."""
        if self.from_structure:
            return self.order == COLLINEAR
        return self.mag

    @property
    def scf_relativistic(self) -> bool:
        return self.so and not self.collinear_scf

    @property
    def nscf_relativistic(self) -> bool:
        return self.so

    def scf(self, labels) -> Dict:
        if self.from_structure:
            if self.order == COLLINEAR:
                spin: Dict = {"nspin": 2}
                spin.update(self._magnitudes(labels))
                return spin
            spin = {"lspinorb": self.so, "noncolin": True}
            spin.update(self._magnitudes(labels))
            spin.update(self._angles(labels))
            return spin
        if self.mag:
            spin = {"nspin": 2}
            spin.update(self._ferromagnetic(labels))
            return spin
        if self.so:
            return {"lspinorb": True, "noncolin": True}
        return {}

    def nscf(self, labels) -> Dict:
        if self.from_structure:
            spin: Dict = {"lspinorb": self.so, "noncolin": True}
            if self.order == COLLINEAR:
                spin["lforcet"] = True
            spin.update(self._magnitudes(labels))
            spin.update(self._angles(labels))
            return spin
        if self.mag:
            spin = {
                "lspinorb": self.so,
                "noncolin": True,
                "lforcet": True,
                "angle1": 0,
                "angle2": 0,
            }
            spin.update(self._ferromagnetic(labels))
            return spin
        if self.so:
            return {"lspinorb": True, "noncolin": True}
        return {}

    @staticmethod
    def _ferromagnetic(labels) -> Dict:
        return {
            f"starting_magnetization({i + 1})": STARTING_MAGNETIZATION
            for i in range(len(labels))
        }

    def _magnitudes(self, labels) -> Dict:
        return {
            f"starting_magnetization({i + 1})": RawValue(
                f"{self.magnetization[label].starting_magnetization:.4f}"
            )
            for i, label in enumerate(labels)
        }

    def _angles(self, labels) -> Dict:
        angles: Dict = {}
        for i, label in enumerate(labels):
            m = self.magnetization[label]
            if m.is_zero:
                continue
            angles[f"angle1({i + 1})"] = RawValue(f"{m.angle1:.4f}")
            angles[f"angle2({i + 1})"] = RawValue(f"{m.angle2:.4f}")
        return angles


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

        partial = next(
            (site for site in structure.sites if abs(site.occupancy - 1.0) > 1.0e-8),
            None,
        )
        if partial is not None:
            raise StructureError(
                f"site {partial.label!r} at {partial.frac_coords} has occupancy "
                f"{partial.occupancy}; partial occupancy is not supported"
            )

        order = magnetic_order(structure)
        if order != NONMAGNETIC:
            structure, magnetic_species = classify_magnetic_sites(structure)
            magnetization = qe_magnetization(magnetic_species, order)
            logger.info(
                "%s magnetic order from the structure; species %s",
                order,
                ", ".join(f"{sp.label}({sp.count})" for sp in magnetic_species),
            )
        else:
            magnetization = {}
        spin = SpinPolicy(order, magnetization, so=config.so, mag=config.mag)

        labels = self._species_labels(structure, hints)
        entries = {
            label: self.table.lookup(self._element_of(structure, label))
            for label in labels
        }
        masses = self._masses(structure, labels, hints)
        ecutwfc, ecutrho = self._cutoffs(entries.values())
        counts = self._counts(structure, entries, spin.spinor)

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
            self._system(base_system, first={}, spin=spin.scf(labels)),
            CONV_THR_SCF,
            relativistic=spin.scf_relativistic,
            kpoints=KPointsAutomatic(kmesh_scf),
        )
        nscf = pw_input(
            self._control("nscf"),
            self._system(
                base_system,
                first={"nosym": True, "nbnd": counts.nbnd_nscf},
                spin=spin.nscf(labels),
            ),
            CONV_THR_NSCF,
            relativistic=spin.nscf_relativistic,
            kpoints=mesh_kpoints_list(kmesh_nscf),
        )
        check_wannier = pw_input(
            self._control("nscf", verbosity="high"),
            self._system(
                base_system, first={"nbnd": counts.nbnd_check}, spin=spin.nscf(labels)
            ),
            CONV_THR_CHECK,
            relativistic=spin.nscf_relativistic,
            kpoints=KPointsAutomatic(kmesh_nscf, (1, 1, 1)),
        )
        # The bands run inherits the check_wannier settings (0.2.x behaviour).
        bands_nscf = pw_input(
            self._control("bands", verbosity="high"),
            self._system(
                base_system, first={"nbnd": counts.nbnd_check}, spin=spin.nscf(labels)
            ),
            CONV_THR_CHECK,
            relativistic=spin.nscf_relativistic,
            kpoints=path,
        )

        wannier90_parameters = dict(WANNIER90_PARAMETERS)
        # wannier_plot reads UNK files written by pw2wannier90.x.  Keeping it
        # enabled while write_unk is false makes Wannier90 abort after the
        # Hamiltonian and band files have otherwise been generated.
        wannier90_parameters["wannier_plot"] = config.write_unk

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
            parameters=wannier90_parameters,
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
            if set(listed) == set(labels):
                return listed  # keep the source's ATOMIC_SPECIES order
            elements = {
                WorkflowBuilder._element_of(structure, label) for label in labels
            }
            if not elements <= set(listed):
                raise InputModelError(
                    f"species {sorted(listed)} from the structure source do not match "
                    f"the sites {sorted(labels)}"
                )
        return labels

    @staticmethod
    def _element_of(structure: NormalizedStructure, label: str) -> str:
        for site in structure.sites:
            if site.label == label:
                return site.element
        raise InputModelError(f"no site with label {label!r}")

    def _masses(self, structure, labels, hints) -> Dict[str, float]:
        masses = {}
        for label in labels:
            element = self._element_of(structure, label)
            if hints.masses is not None and label in hints.masses:
                masses[label] = float(hints.masses[label])
            elif hints.masses is not None and element in hints.masses:
                masses[label] = float(hints.masses[element])
            else:
                masses[label] = atomic_mass(element)
        return masses

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
        spinor: bool,
    ) -> WannierCounts:
        num_wann = sum(entries[site.label].num_wann for site in structure.sites)
        nexclude = sum(entries[site.label].nexclude for site in structure.sites)
        if num_wann == 0:
            raise InputModelError(
                "no Wannier projections: every species has an empty orbitals entry"
            )
        return WannierCounts(num_wann, nexclude, spinor)

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
        entries["write_unk"] = self.config.write_unk
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
    "SpinPolicy",
    "StructureHints",
    "WannierCounts",
    "WorkflowBuilder",
    "build_plan",
    "plan_from_cif2cell_output",
]
