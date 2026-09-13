# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Changed
- The `cif2qewan` command is now implemented by `cif2qewan.cli` on top of
  the reader / workflow builder / renderer modules (DEVELOPMENT_PLAN.md
  Step 7). The command name and the `--so` / `--mag` options are unchanged.
- Generated inputs: scf.in writes `K_POINTS {automatic}` (0.2.x copied
  `K_POINTS automatic` from cif2cell); pwscf.win lists the keywords in a
  fixed layout (`num_wann = ...` with a single space, `spinors` right after
  the counts, all pass-through parameters in one block before
  `projections`). All numbers and all other files are identical to 0.2.x;
  the reference examples are regenerated.
- An existing `cif_scf.in` in the working directory is no longer reused
  implicitly. cif2cell runs in a fresh temporary directory with an argument
  list (no shell), its exit status is checked, and its output is written to
  `cif_scf.in` in the output directory as before. To reuse a cif2cell
  output, pass it explicitly with `--cif2cell-output FILE`.
- Errors in the configuration, the structure, the pseudopotential table or
  cif2cell are reported as `cif2qewan: error: ...` with exit status 1
  instead of a traceback.
- `cif2qewan` no longer uses `docopt`; the dependency stays until the
  deprecated module is removed.

- Error handling (DEVELOPMENT_PLAN.md Step 8): `band_comp` and
  `wannier_conv` report missing or unreadable files as `<prog>: error: ...`
  with exit status 1 instead of a traceback; the Wannier90 Hamiltonian
  reader raises `DataFileError` instead of printing and continuing with
  empty data. `cif2qewan -v/--verbose` logs the decisions taken (species and
  pseudopotentials, cutoffs, band counts, k meshes, the cif2cell command);
  the fallback to the simple-cubic band path without seekpath is reported
  as a warning (0.2.x printed it, the first 0.3 development versions were
  silent). The deprecated 0.2.x module is unchanged and still uses
  `os.system`; it is removed in Step 10.

### Deprecated
- `cif2qewan.cif2qewan.qe_wannier_in` and `cif2qewan.cif2qewan.main`
  (the 0.2.x implementation). `main` delegates to the new CLI;
  `qe_wannier_in` still works but warns. `python -m cif2qewan.cif2qewan`
  keeps working. Both will be removed in a later release.

### Added
- `cif2qewan --reader pymatgen` reads the structure with pymatgen instead
  of cif2cell (primitive cell in the CIF's Cartesian frame; the SCF mesh
  follows `scf_k_resolution` with cif2cell's rule). `--cif2cell-output FILE`
  reuses a cif2cell output, `--output-dir DIR` selects the target directory,
  `--version` prints the version.
- Reference examples `examples/PSLibrary/Fe_nonmag` (no options) and
  `examples/PSLibrary/Fe_so` (`--so` only), generated from the same Fe CIF as
  `examples/PSLibrary/Fe`.
- Characterization tests for the 0.2.x behaviour (DEVELOPMENT_PLAN.md
  Step 1): the physical content of the generated inputs (spin settings per
  option, k-mesh relations, cell and atoms shared by QE and Wannier90, band
  and Wannier-function counts) and the CLI exit status for invalid inputs.
- Typed internal models (DEVELOPMENT_PLAN.md Step 2): `cif2qewan.structure.model` (`NormalizedStructure`, `AtomicSite`,
  Cartesian `MagneticMoment` in Bohr magneton), `cif2qewan.qe.model`
  (namelists, cards and `PwInput`), `cif2qewan.wannier90.model`
  (`Wannier90Input`), `cif2qewan.workflow.model` (`CalculationPlan`) and
  the exception hierarchy in `cif2qewan.exceptions`.
- Structure readers (DEVELOPMENT_PLAN.md Step 3):
  `cif2qewan.structure.readers` with the `StructureReader` protocol, a
  `PymatgenReader` for CIF/MCIF and the other pymatgen formats, and a
  `Cif2cellReader` that runs cif2cell with `subprocess` in a fresh temporary
  directory (or parses an explicitly given cif2cell output) and raises
  `ExternalCommandError` on failure. `cif2qewan.structure.compare` checks
  two structures for physical equivalence (same lattice, same sites up to
  lattice translations, same magnetic moments) within tolerances.
- QE renderer (DEVELOPMENT_PLAN.md Step 4):
  `cif2qewan.qe.writer` turns `PwInput` and `NamelistInput` models into
  text with the 0.2.x formatting. It always writes the `K_POINTS` option in
  braces (`K_POINTS {automatic}`); 0.2.x copied `K_POINTS automatic` from
  cif2cell for scf.in. Both spellings are accepted by pw.x.
- Wannier90 renderer (DEVELOPMENT_PLAN.md Step 5): `cif2qewan.wannier90.writer.render_win` turns a `Wannier90Input`
  into `.win` text with the 0.2.x numeric formats and a fixed layout
  (counts, `spinors`, pass-through parameters, then the blocks).
- Workflow layer (DEVELOPMENT_PLAN.md Step 6):
  `cif2qewan.config.Config` (TOML + `--so`/`--mag`),
  `cif2qewan.qe.pseudopotential.PseudopotentialTable` (the CSV tables, read
  with the standard library), `cif2qewan.qe.kpoints` (SCF mesh from the
  k resolution with cif2cell's rule, the [4, 8] NSCF mesh, the seekpath band
  path with the simple-cubic fallback) and
  `cif2qewan.workflow.builder.WorkflowBuilder`, which turns a structure and
  a `Config` into a `CalculationPlan` with the 0.2.x conventions, including
  the two-step `--mag` policy. `cif2qewan.workflow.render.render_plan`
  renders every input of a plan. From the cif2cell output of each reference
  example the new code reproduces all QE inputs byte for byte (after the
  `K_POINTS` normalization) and pwscf.win up to layout.
- The band-path divisions are now `floor(10 * L / L_min + 1e-6)` instead of
  `int(10 * L / L_min)`, so that a ratio that is an integer up to
  floating-point noise is not truncated to the integer below. Identical for
  the reference examples.

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
