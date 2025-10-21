# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Comprehensive documentation and comments for all Python modules
- NumPy-style docstrings for all functions and classes
- Type hints throughout the codebase
- MIT license
- Detailed README.md with installation and usage instructions
- Python package configuration (pyproject.toml, setup.py)
- Command-line entry points for all tools
- Automated workflow script (submit_all.sh) with detailed comments
- Band structure comparison tool (band_comp.py)
- Wannier90 convergence checker (wannier_conv.py)
- Materials Project CIF downloader (get_cif.py)
- Configuration file support (cif2qewan.toml)
- Git ignore file for Python and QE/Wannier90 projects

### Changed
- Improved code organization and readability
- Enhanced error handling and user feedback
- Better command-line interface with help text
- More robust file parsing and data handling

### Fixed
- Improved import organization
- Consistent coding style throughout the project
- Better error messages and warnings

## [1.0.0] - 2024-01-01

### Added
- Initial release of cif2qewan
- CIF to Quantum ESPRESSO input file generation
- CIF to Wannier90 input file generation
- Support for spin-orbit coupling calculations
- Support for magnetic calculations
- Band structure comparison functionality
- Wannier90 convergence checking
- Materials Project integration
- Configuration file support
- Comprehensive documentation
