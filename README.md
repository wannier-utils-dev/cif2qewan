# cif2qewan

A comprehensive Python toolkit for generating Quantum ESPRESSO and Wannier90 input files from CIF (Crystallographic Information File) structures, with automated workflow management and band structure analysis capabilities.

[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Usage](#usage)
- [Workflow](#workflow)
- [Band Structure Analysis](#band-structure-analysis)
- [Convergence Checking](#convergence-checking)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [References](#references)

## Features

- **Automated Input Generation**: Generate Quantum ESPRESSO and Wannier90 input files from CIF structures
- **Band Structure Comparison**: Compare DFT and Wannier90 band structures with publication-ready plots
- **Convergence Checking**: Automated Wannier90 convergence analysis
- **Workflow Automation**: Complete end-to-end workflow from CIF to analysis
- **Spin-Orbit Coupling Support**: Handle SOC calculations and magnetic systems
- **Flexible Configuration**: TOML-based configuration system

## Installation

### Prerequisites

- Python 3.9 or higher
- Quantum ESPRESSO (QE)
- Wannier90
- cif2cell
- Required Python packages (see below)

### Software Dependencies

1. **Quantum ESPRESSO**: [Download and install QE](https://www.quantum-espresso.org/)
2. **Wannier90**: [Download and install Wannier90](https://github.com/wannier-developers/wannier90)
3. **cif2cell**: [Download and install cif2cell](https://sourceforge.net/projects/cif2cell/)

### Installation

```bash
git clone https://github.com/wannier-utils-dev/cif2qewan.git
cd cif2qewan
pip install .
```

This installs the Python dependencies and the `cif2qewan`, `band_comp` and
`wannier_conv` commands. The pseudopotential tables are installed with the
package; point `pp_list_path` in `cif2qewan.toml` at the one you want, for
example `cif2qewan/pp_psl_rrkj.csv` in the clone.

Without installing, the tools can also be run from a clone as
`python -m cif2qewan.cif2qewan`, `python -m cif2qewan.band_comp` and
`python -m cif2qewan.wannier_conv`.

## Quick Start

1. **Configure the system** by editing `cif2qewan.toml`:

```toml
# Path to cif2cell executable
cif2cell_path = "/path/to/cif2cell"

# Directory containing pseudopotentials
pseudo_dir = "/path/to/pseudopotentials"

# Path to pseudopotential list CSV
pp_list_path = "/path/to/pp_list.csv"

# K-point resolution for SCF (1/Å)
scf_k_resolution = 0.15

# Gaussian smearing (Ry)
degauss = 0.01

# pw2wannier90 configuration
[pw2wan]
write_unk = ".true."
```

2. **Run the complete workflow**:

```bash
# Generate input files and run calculations
./submit_all.sh

# Or run step by step
cif2qewan structure.cif cif2qewan.toml   # or: python -m cif2qewan.cif2qewan structure.cif cif2qewan.toml
# ... run QE and Wannier90 calculations ...
python -m cif2qewan.band_comp -o ./
```

## Configuration

### TOML Configuration File

The `cif2qewan.toml` file contains all necessary configuration parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `cif2cell_path` | Path to cif2cell executable | Required |
| `pseudo_dir` | Directory containing pseudopotentials | Required |
| `pp_list_path` | Path to pseudopotential list CSV | Required |
| `scf_k_resolution` | K-point resolution for SCF (1/Å) | 0.15 |
| `degauss` | Gaussian smearing (Ry) | 0.01 |
| `pw2wan.write_unk` | Write UNK files for Wannier90 | ".true." |

### Pseudopotential List Format

Create a CSV file with the following columns:

```csv
atom,pp_file_name,nexclude,orbitals,ecutwfc,ecutrho
Fe,Fe.pbe-n-rrkjus_psl.1.0.0.UPF,0,spd,40.0,200.0
O,O.pbe-n-rrkjus_psl.1.0.0.UPF,0,sp,40.0,200.0
```

## Usage

### Basic Usage

```bash
# Generate input files from CIF
cif2qewan structure.cif cif2qewan.toml

# With spin-orbit coupling
cif2qewan structure.cif cif2qewan.toml --so

# With magnetic calculations
cif2qewan structure.cif cif2qewan.toml --mag
```

### Command Line Options

| Option | Description |
|--------|-------------|
| `--so` | Include spin-orbit coupling |
| `--mag` | Perform magnetic calculations |

### Generated Files

The script generates the following input files:

- `scf.in` - SCF calculation input
- `nscf.in` - NSCF calculation input
- `pw2wan.in` - pw2wannier90 interface input
- `pwscf.win` - Wannier90 input
- `band/` - Band structure calculation files
- `check_wannier/` - Convergence check files

## Workflow

### Complete Automated Workflow

```bash
# Run the complete workflow
./submit_all.sh
```

This script performs the following steps:

1. **Generate input files** from CIF structure
2. **Run SCF calculation** for ground state
3. **Run NSCF calculation** for Wannier90
4. **Run Wannier90 preprocessing** and interpolation
5. **Check convergence** by comparing energies
6. **Generate band structure** plots
7. **Compare DFT and Wannier90** band structures

### Manual Workflow

```bash
# Step 1: Generate input files
cif2qewan structure.cif cif2qewan.toml

# Step 2: Run SCF calculation
mpirun -n 16 pw.x < scf.in > scf.out

# Step 3: Run NSCF calculation
mpirun -n 16 pw.x < nscf.in > nscf.out

# Step 4: Run Wannier90 preprocessing
wannier90.x -pp pwscf

# Step 5: Run pw2wannier90 interface
pw2wannier90.x < pw2wan.in > pw2wan.out

# Step 6: Set frozen window and run Wannier90
# Edit dis_froz_max in pwscf.win (recommended: EF + 1-3 eV)
wannier90.x pwscf

# Step 7: Check convergence
cd check_wannier
mpirun -n 16 pw.x < nscf.in > nscf.out
cd ..
wannier_conv -e 5.0 -o ./

# Step 8: Generate band structure
cd band
mpirun -n 16 pw.x < nscf.in > nscf.out
mpirun -n 16 bands.x < band.in > band.out
cd ..

# Step 9: Compare band structures
python -m cif2qewan.band_comp -o ./
```

## Band Structure Analysis

### Generate Band Structure Plots

```bash
# Run band structure calculation
cd band
mpirun -n 16 pw.x < nscf.in > nscf.out
mpirun -n 16 bands.x < band.in > band.out
cd ..

# Generate comparison plot
python -m cif2qewan.band_comp -o ./
```

### Output Files

- `band_compare.png` - Band structure comparison plot (PNG format)
- `band_compare.eps` - Band structure comparison plot (EPS format)

The plot shows:
- **Red lines**: DFT band structure
- **Black lines**: Wannier90 interpolated band structure
- **Vertical lines**: High-symmetry points
- **Energy axis**: Relative to Fermi energy

## Convergence Checking

### Check Wannier90 Convergence

```bash
# Run convergence check
wannier_conv -e 5.0 -o ./ -i ./check_wannier/nscf.out
```

### Convergence Metrics

The script calculates two convergence metrics:

1. **Average difference**: $\delta_{avg} = \sqrt{\frac{1}{N} \sum_{n,k} (E_{n,k}^{DFT} - E_{n,k}^{Wannier})^2}$

2. **Maximum difference**: $\delta_{max} = \max_{n,k} |E_{n,k}^{DFT} - E_{n,k}^{Wannier}|$

### Output Files

- `CONV_5.0` - Convergence results for energy window up to 5 eV above Fermi level

### Interpretation

- **Good convergence**: $\delta_{avg} < 0.01$ eV
- **Acceptable convergence**: $\delta_{avg} < 0.1$ eV
- **Poor convergence**: $\delta_{avg} > 0.1$ eV

## Examples

### Example 1: Iron (Fe) Crystal

```bash
# Generate input files (see examples/PSLibrary/Fe/ for the CIF file)
cif2qewan mp-13_Fe.cif cif2qewan.toml --mag

# Run calculations
./submit_all.sh
```

### Example 2: Spin-Orbit Coupling

```bash
# Generate input with SOC
cif2qewan structure.cif cif2qewan.toml --so

# Run calculations
./submit_all.sh
```

### Example 3: Custom Configuration

```toml
# cif2qewan.toml
cif2cell_path = "/usr/local/bin/cif2cell"
pseudo_dir = "/home/user/pseudopotentials"
pp_list_path = "/home/user/pp_list.csv"
scf_k_resolution = 0.20
degauss = 0.02

[pw2wan]
write_unk = ".false."
```

## Troubleshooting

### Common Issues

1. **cif2cell not found**
   ```bash
   # Install cif2cell
   # Add to PATH or update cif2cell_path in config
   ```

2. **Pseudopotentials not found**
   ```bash
   # Check pseudo_dir path
   # Ensure pseudopotential files exist
   ```

3. **Wannier90 convergence issues**
   ```bash
   # Adjust dis_froz_max in pwscf.win
   # Check projection settings
   # Increase k-point mesh density
   ```

4. **Memory issues**
   ```bash
   # Reduce number of MPI processes
   # Use smaller k-point mesh
   # Check available memory
   ```

### Debug Mode

```bash
# Run with verbose output
cif2qewan structure.cif cif2qewan.toml --verbose

# Check intermediate files
ls -la work/
ls -la check_wannier/
ls -la band/
```

### Performance Optimization

1. **Parallel execution**: Use appropriate number of MPI processes
2. **Memory management**: Monitor memory usage during calculations
3. **Disk space**: Ensure sufficient disk space for work directories
4. **Network**: For cluster calculations, use fast interconnect

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone repository
git clone https://github.com/wannier-utils-dev/cif2qewan.git
cd cif2qewan

# Install with the test dependencies
pip install -e '.[test]'

# Run tests
pytest

# Run linting
flake8 cif2qewan tests
```

### Reporting Issues

Please report issues on our [GitHub Issues](https://github.com/wannier-utils-dev/cif2qewan/issues) page.

## References

### Scientific Papers

1. **Iron-based binary ferromagnets for transverse thermoelectric conversion**
   - A. Sakai, S. Minami, T. Koretsune et al.
   - Nature 581, 53-57 (2020)
   - [DOI: 10.1038/s41586-020-2230-z](https://doi.org/10.1038/s41586-020-2230-z)
   - *The database of anomalous Hall conductivity and anomalous Nernst conductivity is generated using cif2qewan.py.*

2. **Systematic first-principles study of the on-site spin-orbit coupling in crystals**
   - Phys. Rev. B 102, 045109 (2020)
   - [DOI: 10.1103/PhysRevB.102.045109](https://doi.org/10.1103/PhysRevB.102.045109)
   - *The spin-orbit couplings are extracted from the tight-binding models generated by cif2qewan.py.*

### Software References

- [Quantum ESPRESSO](https://www.quantum-espresso.org/)
- [Wannier90](https://github.com/wannier-developers/wannier90)
- [cif2cell](https://sourceforge.net/projects/cif2cell/)
- [Materials Project](https://materialsproject.org/)

## License

This project is licensed under the GNU General Public License v2 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Quantum ESPRESSO developers
- Wannier90 developers
- cif2cell developers
- Materials Project team
- pymatgen developers

---

For more information, please visit our [GitHub repository](https://github.com/wannier-utils-dev/cif2qewan) or contact the maintainers.
