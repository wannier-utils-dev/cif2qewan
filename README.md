# cif2qewan

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

cif2qewan generates Quantum ESPRESSO and Wannier90 input files from a CIF (or
MagCIF) structure: the SCF and NSCF runs, the pw2wannier90 and Wannier90
inputs, and the inputs of a band-structure run and of a convergence check.
Pseudopotentials, cutoffs, band counts, Wannier projections, k meshes and the
band path are chosen from bundled tables for PSLibrary and PseudoDojo. Two
scripts, `wannier_conv` and `band_comp`, compare the Wannier interpolation
with the DFT bands.

## Installation

Python 3.9 or later. Quantum ESPRESSO and Wannier90 are needed to run the
generated inputs, not to generate them. The structure is read by cif2cell by
default (optional: `pip install cif2cell`, or `--reader pymatgen` needs no
external program).

```bash
git clone https://github.com/wannier-utils-dev/cif2qewan.git
cd cif2qewan
pip install .              # or: pip install '.[cif2cell]' to install cif2cell too
```

This installs the Python dependencies (numpy, pymatgen, seekpath, spglib,
matplotlib, toml), the commands `cif2qewan`, `wannier_conv` and `band_comp`,
and the pseudopotential tables. Without installing, the commands can be run
from a clone as `python -m cif2qewan.cli`, `python -m cif2qewan.wannier_conv`
and `python -m cif2qewan.band_comp`.

## Quick start

Write a `cif2qewan.toml` (see the commented sample `cif2qewan.toml` in the
repository; only these keys are required):

```toml
pseudo_dir = "/path/to/pseudopotentials"   # PSLibrary files for the default table
scf_k_resolution = 0.15                    # SCF k-point spacing (1/Å)
degauss = 0.01                             # smearing width (Ry)

[pw2wan]
write_unk = true
```

Generate the inputs:

```bash
cif2qewan structure.cif cif2qewan.toml     # options: --so, --mag, --reader pymatgen, -v
```

This writes `scf.in`, `nscf.in`, `pw2wan.in`, `pwscf.win` and the inputs of
the band-structure run (`band/`) and of the convergence check
(`check_wannier/`). `wannier_conv` and `band_comp` then compare the Wannier
interpolation with the DFT bands. `submit_all.sh` runs the whole chain, from
`cif2qewan` to the comparison plots; set the paths at its top first.

See the [user guide](docs/usage.md) for the command-line options, the
configuration keys, the pseudopotential tables (PSLibrary and PseudoDojo),
MagCIF magnetism, `use_ibrav`, `use_symwan`, the step-by-step workflow and
the Python API.

## Contributing

Bug reports and pull requests are welcome. Work on a branch off `develop`
and open the pull request against `develop`; CI runs the tests on Python 3.9
and 3.12.

```bash
pip install -e '.[test]'
pytest                                           # no QE, Wannier90 or cif2cell needed
flake8 cif2qewan tests --select=E9,F63,F7,F82    # what CI checks
black cif2qewan tests
```

The scientific choices (cutoffs, band counts, k meshes, spin settings) live
in `cif2qewan/workflow/builder.py`; changing one changes the generated
inputs, which the tests compare byte for byte with the reference examples
under `examples/`. Such a change must be stated in the pull request together
with the regenerated examples.

## References

1. A. Sakai, S. Minami, T. Koretsune et al., *Iron-based binary ferromagnets
   for transverse thermoelectric conversion*, Nature 581, 53-57 (2020),
   [10.1038/s41586-020-2230-z](https://doi.org/10.1038/s41586-020-2230-z):
   the anomalous Hall and Nernst conductivity database was generated with
   cif2qewan.
2. *Systematic first-principles study of the on-site spin-orbit coupling in
   crystals*, Phys. Rev. B 102, 045109 (2020),
   [10.1103/PhysRevB.102.045109](https://doi.org/10.1103/PhysRevB.102.045109):
   the spin-orbit couplings were extracted from tight-binding models
   generated with cif2qewan.
3. T. Koretsune, *Construction of maximally-localized Wannier functions using
   crystal symmetry*, Comput. Phys. Commun. 285, 108645 (2023),
   [10.1016/j.cpc.2022.108645](https://doi.org/10.1016/j.cpc.2022.108645)
   (symWannier).

Software: [Quantum ESPRESSO](https://www.quantum-espresso.org/),
[Wannier90](https://github.com/wannier-developers/wannier90),
[symWannier](https://github.com/wannier-utils-dev/symWannier),
[cif2cell](https://sourceforge.net/projects/cif2cell/),
[pymatgen](https://pymatgen.org/),
[seekpath](https://github.com/giovannipizzi/seekpath),
[spglib](https://spglib.readthedocs.io/),
[PSLibrary](https://dalcorso.github.io/pslibrary/),
[PseudoDojo](https://www.pseudo-dojo.org/).

## License

MIT, see [LICENSE](LICENSE).
