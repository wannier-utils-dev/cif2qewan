# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- Reference examples `examples/PSLibrary/Fe_nonmag` (no options) and
  `examples/PSLibrary/Fe_so` (`--so` only), generated from the same Fe CIF as
  `examples/PSLibrary/Fe`.
- Characterization tests for the 0.2.x behaviour (DEVELOPMENT_PLAN.md
  Step 1): the physical content of the generated inputs (spin settings per
  option, k-mesh relations, cell and atoms shared by QE and Wannier90, band
  and Wannier-function counts) and the CLI exit status for invalid inputs.
- Typed internal models (DEVELOPMENT_PLAN.md Step 2), not yet used by the
  CLI: `cif2qewan.structure.model` (`NormalizedStructure`, `AtomicSite`,
  Cartesian `MagneticMoment` in Bohr magneton), `cif2qewan.qe.model`
  (namelists, cards and `PwInput`), `cif2qewan.wannier90.model`
  (`Wannier90Input`), `cif2qewan.workflow.model` (`CalculationPlan`) and
  the exception hierarchy in `cif2qewan.exceptions`.
- Structure readers (DEVELOPMENT_PLAN.md Step 3), not yet used by the CLI:
  `cif2qewan.structure.readers` with the `StructureReader` protocol, a
  `PymatgenReader` for CIF/MCIF and the other pymatgen formats, and a
  `Cif2cellReader` that runs cif2cell with `subprocess` in a fresh temporary
  directory (or parses an explicitly given cif2cell output) and raises
  `ExternalCommandError` on failure. `cif2qewan.structure.compare` checks
  two structures for physical equivalence (same lattice, same sites up to
  lattice translations, same magnetic moments) within tolerances.
- QE renderer (DEVELOPMENT_PLAN.md Step 4), not yet used by the CLI:
  `cif2qewan.qe.writer` turns `PwInput` and `NamelistInput` models into
  text with the 0.2.x formatting. It always writes the `K_POINTS` option in
  braces (`K_POINTS {automatic}`); 0.2.x copied `K_POINTS automatic` from
  cif2cell for scf.in. Both spellings are accepted by pw.x.
- Wannier90 renderer (DEVELOPMENT_PLAN.md Step 5), not yet used by the
  CLI: `cif2qewan.wannier90.writer.render_win` turns a `Wannier90Input`
  into `.win` text with the 0.2.x numeric formats and a fixed layout
  (counts, `spinors`, pass-through parameters, then the blocks).

## [0.2.0] - 2026-09-12

### Added
- NumPy-style docstrings, comments and type hints throughout the code.
- Packaging: the scripts moved into a `cif2qewan` package that can be
  installed with `pip install .`, with `cif2qewan`, `band_comp` and
  `wannier_conv` console scripts.
- Tests: the reference examples under `examples/` are regenerated and compared
  byte for byte, plus unit tests for the scf.out/pwscf.win parsers.
- A CI workflow running those tests on Python 3.9 and 3.12.
- A `.gitignore` for Python and Quantum ESPRESSO / Wannier90 working files.
- Expanded README with installation, configuration and workflow instructions.

### Changed
- The scripts are now invoked as `cif2qewan ...` (or
  `python -m cif2qewan.cif2qewan ...`) instead of `python cif2qewan.py ...`.
- The pseudopotential CSV tables moved from the repository root into the
  `cif2qewan` package, so that they are installed with it. `pp_list_path` in
  `cif2qewan.toml` is an absolute path chosen by the user and is unaffected.
- `band_comp` and `wannier_conv` now require `-o`; both already failed
  without it.
- Python 3.9 or later is required, because of the new type hints.

### Fixed
- `cif2cell_scf_in()` used `sys.argv[1]` instead of its `cif_file` argument.
- `read_set_system()` now raises instead of continuing with an empty system
  when ntyp/nat cannot be read from the cif2cell output.
- `submit_all.sh` started with a Python docstring, which bash tried to run as
  a command.
- `matplotlib.use("Agg")` is called before `pyplot` is imported again, so
  plotting works on machines without a display.

- The license changed from GPLv2 to MIT, with the agreement of the copyright
  holders.

### Not changed
- The numerical behaviour of the generated inputs is unchanged: both
  reference examples are reproduced exactly.

## [0.1.0] - 2026-09-12

### Added
- Historical release of the original script-based version from the former
  `master` branch, preserved for existing users before the packaging and
  documentation changes in 0.2.0.
