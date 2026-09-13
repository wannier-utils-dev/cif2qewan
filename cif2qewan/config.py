"""User configuration of cif2qewan: the TOML file plus the CLI options."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Union

import toml

from cif2qewan.exceptions import ConfigError

PathLike = Union[str, Path]


@dataclass(frozen=True)
class Config:
    """Everything the workflow builder needs besides the structure.

    Attributes
    ----------
    cif2cell_path : str
        The cif2cell executable (used by the cif2cell reader only).
    pseudo_dir : str
        ``pseudo_dir`` written into the QE inputs.
    pp_list_path : str
        CSV table of pseudopotentials and projections.
    scf_k_resolution : float
        Target k-space resolution of the SCF mesh in 1/angstrom.
    degauss : float
        Smearing width in Ry.
    pw2wan : dict
        Options copied into ``pw2wan.in``: ``write_unk`` (required) and
        optionally ``wannier_plot_supercell``.
    so, mag : bool
        The ``--so`` and ``--mag`` command-line options.
    """

    cif2cell_path: str
    pseudo_dir: str
    pp_list_path: str
    scf_k_resolution: float
    degauss: float
    pw2wan: Dict[str, Any] = field(default_factory=dict)
    so: bool = False
    mag: bool = False

    def __post_init__(self) -> None:
        for name in ("cif2cell_path", "pseudo_dir", "pp_list_path"):
            if not str(getattr(self, name)).strip():
                raise ConfigError(f"{name} must not be empty")
        if not (float(self.scf_k_resolution) > 0.0):
            raise ConfigError(
                f"scf_k_resolution must be positive, got {self.scf_k_resolution}"
            )
        if not (float(self.degauss) > 0.0):
            raise ConfigError(f"degauss must be positive, got {self.degauss}")
        if "write_unk" not in self.pw2wan:
            raise ConfigError("the [pw2wan] table needs write_unk")

    @property
    def spinor(self) -> bool:
        """True when the calculation uses two-component spinors (``--so`` or ``--mag``)."""
        return self.so or self.mag

    @classmethod
    def from_dict(
        cls, data: Dict[str, Any], so: bool = False, mag: bool = False
    ) -> "Config":
        required = (
            "cif2cell_path",
            "pseudo_dir",
            "pp_list_path",
            "scf_k_resolution",
            "degauss",
        )
        missing = [key for key in required if key not in data]
        if missing:
            raise ConfigError(f"missing configuration keys: {', '.join(missing)}")
        pw2wan = data.get("pw2wan")
        if not isinstance(pw2wan, dict):
            raise ConfigError("the configuration needs a [pw2wan] table")
        try:
            return cls(
                cif2cell_path=str(data["cif2cell_path"]),
                pseudo_dir=str(data["pseudo_dir"]),
                pp_list_path=str(data["pp_list_path"]),
                scf_k_resolution=float(data["scf_k_resolution"]),
                degauss=float(data["degauss"]),
                pw2wan=dict(pw2wan),
                so=bool(so),
                mag=bool(mag),
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
