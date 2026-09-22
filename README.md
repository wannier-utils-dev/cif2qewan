# cif2qewan

A comprehensive Python toolkit for generating Quantum ESPRESSO and Wannier90 input files from CIF (Crystallographic Information File) structures, with automated workflow management and band structure analysis capabilities.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
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
- Quantum ESPRESSO (QE) and Wannier90 to run the generated inputs
- cif2cell (optional): the default structure reader runs it;
  `--reader pymatgen` reads the CIF without it

### Software Dependencies

1. **Quantum ESPRESSO**: [Download and install QE](https://www.quantum-espresso.org/)
2. **Wannier90**: [Download and install Wannier90](https://github.com/wannier-developers/wannier90)
3. **cif2cell** (optional): `pip install cif2cell`, or see
   [cif2cell](https://sourceforge.net/projects/cif2cell/)

### Installation

```bash
git clone https://github.com/wannier-utils-dev/cif2qewan.git
cd cif2qewan
pip install .              # or: pip install '.[cif2cell]' to install cif2cell too
```

This installs the Python dependencies, the `cif2qewan`, `band_comp` and
`wannier_conv` commands, and the pseudopotential tables (`pp_psl_rrkj.csv`
for PSLibrary, `nc-sr-05_pbe_standard_upf.csv` and
`nc-sr-05_pbe_stringent_upf.csv` for PseudoDojo). The PSLibrary table is
used by default; select another one by its file name with `pp_list_path`.

Without installing, the tools can also be run from a clone as
`python -m cif2qewan.cli`, `python -m cif2qewan.band_comp` and
`python -m cif2qewan.wannier_conv`.

## Quick Start

1. **Configure the system** by editing `cif2qewan.toml`:

```toml
# Directory containing pseudopotentials
pseudo_dir = "/path/to/pseudopotentials"

# Pseudopotential table; default: the bundled PSLibrary table
# pp_list_path = "nc-sr-05_pbe_standard_upf.csv"

# cif2cell executable; default: "cif2cell" on PATH
# cif2cell_path = "/path/to/cif2cell"

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
cif2qewan structure.cif cif2qewan.toml   # or: python -m cif2qewan.cli structure.cif cif2qewan.toml
# ... run QE and Wannier90 calculations ...
python -m cif2qewan.band_comp -o ./
```

## Configuration

### TOML Configuration File

The `cif2qewan.toml` file contains all necessary configuration parameters:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `pseudo_dir` | Directory containing pseudopotentials | Required |
| `scf_k_resolution` | K-point resolution for SCF (1/Å) | Required (0.15 in the sample) |
| `degauss` | Gaussian smearing (Ry) | Required (0.01 in the sample) |
| `pw2wan.write_unk` | Write UNK files for Wannier90 | Required (".true." in the sample) |
| `pp_list_path` | Pseudopotential table: the file name of a bundled table, or a path | `"pp_psl_rrkj.csv"` |
| `cif2cell_path` | cif2cell executable (default reader only) | `"cif2cell"` on PATH |

When `pw2wan.write_unk` is false, the generated `pwscf.win` also sets
`wannier_plot = .false.` because Wannier-function plots require the UNK files.

### Pseudopotential Tables

Three tables are installed with the package and selected by file name:
`pp_psl_rrkj.csv` (PSLibrary, the default), `nc-sr-05_pbe_standard_upf.csv`
and `nc-sr-05_pbe_stringent_upf.csv` (PseudoDojo). A `pp_list_path` with a
directory part (`./my_table.csv`, `/path/to/table.csv`) is used as a path,
so you can also write your own table with the following columns:

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
| `--reader {cif2cell,pymatgen}` | How to read the structure. `cif2cell` (default) runs cif2cell as in earlier versions; `pymatgen` reads the CIF directly and needs no cif2cell |
| `--cif2cell-output FILE` | Reuse an existing cif2cell output (`cif_scf.in`) instead of running cif2cell |
| `--output-dir DIR` | Write the inputs into `DIR` (default: the current directory) |
| `--version` | Print the version |

cif2cell is run in a temporary directory and its output is saved as
`cif_scf.in` next to the generated inputs. Unlike earlier versions, an
existing `cif_scf.in` is not picked up automatically; pass it with
`--cif2cell-output` if you want to reuse it. Errors are reported as
`cif2qewan: error: ...` with exit status 1.

The `pymatgen` reader accepts only fully occupied sites; mixed or partial
occupancies are rejected because the generated QE inputs cannot represent them.

### Magnetic structures (MagCIF)

With `--reader pymatgen` a MagCIF (`.mcif`) file is read including the site
moments. Sites of one element with different moments become separate QE
species (`Mn1`, `Mn2`, ...), and `starting_magnetization(i)`, `angle1(i)` and
`angle2(i)` are set from the moments (relative to the largest one). A
collinear structure is run like `--mag` (collinear SCF, noncollinear NSCF
with `lforcet`), a noncollinear one is noncollinear from the SCF on; add
`--so` for spin-orbit coupling. The Wannier functions are spinors in both
cases. Moments in the file take precedence over `--mag`.

```bash
cif2qewan Mn3Sn.mcif cif2qewan.toml --so --reader pymatgen
```

### Generated Files

The script generates the following input files:

- `cif_scf.in` - the cif2cell output (not written with `--reader pymatgen`)
- `scf.in` - SCF calculation input
- `nscf.in` - NSCF calculation input
- `pw2wan.in` - pw2wannier90 interface input
- `pwscf.win` - Wannier90 input
- `band/` - Band structure calculation files (`nscf.in`, `band.in`, `proj.in`, `pp.in`)
- `check_wannier/` - Convergence check files (`nscf.in`)

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
# Generate input files (see examples/PSLibrary/Fe*/ for the CIF file and
# the reference outputs without options, with --so, and with --so --mag)
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
pseudo_dir = "/home/user/pseudopotentials"
pp_list_path = "/home/user/pp_list.csv"
cif2cell_path = "/usr/local/bin/cif2cell"
scf_k_resolution = 0.20
degauss = 0.02

[pw2wan]
write_unk = ".false."
```

## Troubleshooting

### Common Issues

1. **cif2cell not found**
   ```bash
   pip install cif2cell          # puts the cif2cell command on PATH
   # or set cif2cell_path in cif2qewan.toml, or use --reader pymatgen
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

Bug reports and pull requests are welcome. Work on a branch off `develop`
and open the pull request against `develop`; CI runs the tests on Python 3.9
and 3.12.

```bash
git clone https://github.com/wannier-utils-dev/cif2qewan.git
cd cif2qewan
pip install -e '.[test]'
pytest                                           # no QE, Wannier90 or cif2cell needed
flake8 cif2qewan tests --select=E9,F63,F7,F82    # what CI checks
black cif2qewan tests                            # formatting
```

The code follows PEP 8 with Black formatting, NumPy-style docstrings and type
hints. The package is a pipeline reader -> workflow builder -> renderers ->
output; the scientific choices (cutoffs, band counts, k meshes, spin
settings) live in `cif2qewan/workflow/builder.py` only, and changing one
changes the generated inputs. Such a change must be stated in the pull
request together with the regenerated reference examples under `examples/`,
which the tests compare byte for byte.

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

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Quantum ESPRESSO developers
- Wannier90 developers
- cif2cell developers
- Materials Project team
- pymatgen developers

---

For more information, please visit our [GitHub repository](https://github.com/wannier-utils-dev/cif2qewan) or contact the maintainers.
