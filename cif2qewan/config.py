"""User configuration of cif2qewan: the TOML file plus the CLI options."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Union

import toml

from cif2qewan.exceptions import ConfigError
from cif2qewan.qe.pseudopotential import DEFAULT_TABLE, resolve_table_path

PathLike = Union[str, Path]

#: ``cif2cell_path`` when the TOML file gives none: the command on PATH,
#: which is where ``pip install cif2cell`` puts it.
DEFAULT_CIF2CELL = "cif2cell"


@dataclass(frozen=True)
class Config:
    """Everything the workflow builder needs besides the structure.

    Attributes
    ----------
    pseudo_dir : str
        ``pseudo_dir`` written into the QE inputs.
    cif2cell_path : str
        The cif2cell executable (used by the cif2cell reader only); by
        default the ``cif2cell`` command on PATH.
    pp_list_path : str
        CSV table of pseudopotentials and projections. The bare file name of
        a bundled table selects the installed copy
        (:func:`cif2qewan.qe.pseudopotential.resolve_table_path`); the
        attribute holds the resolved path. Default: ``pp_psl_rrkj.csv``.
    scf_k_resolution : float
        Target k-space resolution of the SCF mesh in 1/angstrom.
    degauss : float
        Smearing width in Ry.
    pw2wan : dict
        Options copied into ``pw2wan.in``: ``write_unk`` (required) and
        optionally ``wannier_plot_supercell``. ``write_unk`` accepts a TOML
        boolean or the legacy ``.true.`` / ``.false.`` strings.
    so, mag : bool
        The ``--so`` and ``--mag`` command-line options.
    use_ibrav : bool
        Write the cell as QE's ``ibrav`` and ``A``, ``B``, ``C``, ``cosAB``,
        ... instead of ``ibrav = 0`` with ``CELL_PARAMETERS``
        (:mod:`cif2qewan.qe.bravais`). Default false.
    """

    pseudo_dir: str
    scf_k_resolution: float
    degauss: float
    cif2cell_path: str = DEFAULT_CIF2CELL
    pp_list_path: str = DEFAULT_TABLE
    pw2wan: Dict[str, Any] = field(default_factory=dict)
    so: bool = False
    mag: bool = False
    use_ibrav: bool = False

    def __post_init__(self) -> None:
        for name in ("cif2cell_path", "pseudo_dir", "pp_list_path"):
            if not str(getattr(self, name)).strip():
                raise ConfigError(f"{name} must not be empty")
        object.__setattr__(
            self, "pp_list_path", str(resolve_table_path(self.pp_list_path))
        )
        if not (float(self.scf_k_resolution) > 0.0):
            raise ConfigError(
                f"scf_k_resolution must be positive, got {self.scf_k_resolution}"
            )
        if not (float(self.degauss) > 0.0):
            raise ConfigError(f"degauss must be positive, got {self.degauss}")
        if "write_unk" not in self.pw2wan:
            raise ConfigError("the [pw2wan] table needs write_unk")
        # Keep accepting the Fortran-style strings used by the legacy TOML
        # files, but reject values that would produce an invalid QE namelist.
        self._logical_option(self.pw2wan["write_unk"], "pw2wan.write_unk")

    @staticmethod
    def _logical_option(value: Any, name: str) -> bool:
        if isinstance(value, bool):
            return value
        normalized = str(value).strip().lower()
        if normalized in {".true.", "true"}:
            return True
        if normalized in {".false.", "false"}:
            return False
        raise ConfigError(f"{name} must be a boolean or .true./.false., got {value!r}")

    @property
    def write_unk(self) -> bool:
        """Whether ``pw2wannier90.x`` writes the UNK files needed for WF plots."""
        return self._logical_option(self.pw2wan["write_unk"], "pw2wan.write_unk")

    @property
    def spinor(self) -> bool:
        """True when the calculation uses two-component spinors (``--so`` or ``--mag``)."""
        return self.so or self.mag

    @classmethod
    def from_dict(
        cls, data: Dict[str, Any], so: bool = False, mag: bool = False
    ) -> "Config":
        required = ("pseudo_dir", "scf_k_resolution", "degauss")
        missing = [key for key in required if key not in data]
        if missing:
            raise ConfigError(f"missing configuration keys: {', '.join(missing)}")
        pw2wan = data.get("pw2wan")
        if not isinstance(pw2wan, dict):
            raise ConfigError("the configuration needs a [pw2wan] table")
        use_ibrav = data.get("use_ibrav", False)
        if not isinstance(use_ibrav, bool):
            raise ConfigError(f"use_ibrav must be true or false, got {use_ibrav!r}")
        try:
            return cls(
                pseudo_dir=str(data["pseudo_dir"]),
                scf_k_resolution=float(data["scf_k_resolution"]),
                degauss=float(data["degauss"]),
                cif2cell_path=str(data.get("cif2cell_path", DEFAULT_CIF2CELL)),
                pp_list_path=str(data.get("pp_list_path", DEFAULT_TABLE)),
                pw2wan=dict(pw2wan),
                so=bool(so),
                mag=bool(mag),
                use_ibrav=use_ibrav,
            )
        except (TypeError, ValueError) as exc:
            raise ConfigError(f"invalid configuration value: {exc}") from exc

    @classmethod
    def from_toml(cls, path: PathLike, so: bool = False, mag: bool = False) -> "Config":
        path = Path(path)
        if not path.is_file():
            raise ConfigError(f"configuration file not found: {path}")
        try:
            data = toml.load(str(path))
        except (toml.TomlDecodeError, OSError) as exc:
            raise ConfigError(f"cannot read {path}: {exc}") from exc
        return cls.from_dict(data, so=so, mag=mag)
