# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

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

### Not changed
- The license is still GPLv2. Relicensing to MIT is proposed separately.
- The numerical behaviour of the generated inputs is unchanged: both
  reference examples are reproduced exactly.
