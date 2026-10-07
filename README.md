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

## Usage

```bash
cif2qewan structure.cif cif2qewan.toml [--so] [--mag] [--reader pymatgen] [-v]
```

| Option | Description |
|--------|-------------|
| `--so` | Include spin-orbit coupling |
| `--mag` | Ferromagnetic calculation (collinear SCF, noncollinear NSCF with `lforcet`) |
| `--reader {cif2cell,pymatgen}` | `cif2cell` (default) runs cif2cell as in earlier versions; `pymatgen` reads the file directly (CIF, MagCIF and the other pymatgen formats) |
| `--cif2cell-output FILE` | Reuse an existing cif2cell output (`cif_scf.in`) instead of running cif2cell |
| `--output-dir DIR`, `-o DIR` | Write the inputs into `DIR` (default: the current directory) |
| `-v`, `--verbose` | Report the decisions taken (pseudopotentials, cutoffs, band counts, k meshes, `ibrav`, the cif2cell command) |
| `--version` | Print the version |

cif2cell runs in a temporary directory and its output is saved as
`cif_scf.in` next to the generated inputs; an existing `cif_scf.in` is not
picked up unless given with `--cif2cell-output`. The pymatgen reader reduces
the structure to its primitive cell (keeping the input cell when that would
fold sites with different moments) and rejects partially occupied sites.
Errors are reported as `cif2qewan: error: ...` with exit status 1. The 0.2.x
entry point `python -m cif2qewan.cif2qewan` still works but is deprecated and
will be removed in 0.4.0.

Generated files: `scf.in`, `nscf.in`, `pw2wan.in`, `pwscf.win`,
`check_wannier/nscf.in` (shifted k mesh for `wannier_conv`), `band/nscf.in`,
`band/band.in`, `band/proj.in`, `band/pp.in`, and `cif_scf.in` on the
cif2cell path. The reference outputs under `examples/` show them for bcc Fe
without options, with `--so`, and with `--so --mag`.

### Magnetic structures (MagCIF)

With `--reader pymatgen` a MagCIF (`.mcif`) file is read including the site
moments. Sites of one element with different moments become separate QE
species (`Mn1`, `Mn2`, ...), and `starting_magnetization(i)` (the moment in
Bohr magneton, as Quantum ESPRESSO 7.3 and later read it; older versions
clamp it to the fully polarized value), `angle1(i)` and `angle2(i)` are set
from the moments. A collinear order is run like `--mag` (collinear SCF,
noncollinear NSCF with `lforcet`, the axis angles written for every species
because `lforcet` uses those of atomic type 1), a noncollinear one is
noncollinear from the SCF on; add `--so` for spin-orbit coupling. The
Wannier functions are spinors in both cases. Moments in the file take
precedence over `--mag`. When every moment is below 1 Bohr magneton, QE
reads the values as polarization per valence electron and cif2qewan warns.

```bash
cif2qewan Mn3Sn.mcif cif2qewan.toml --so --reader pymatgen
```

## Configuration

```toml
pseudo_dir = "/path/to/pseudopotentials"   # written into the QE inputs
scf_k_resolution = 0.15                    # SCF k-point spacing (1/Å)
degauss = 0.01                             # smearing width (Ry)
# pp_list_path = "nc-sr-05_pbe_standard_upf.csv"   # default: pp_psl_rrkj.csv
# cif2cell_path = "/path/to/cif2cell"              # default: cif2cell on PATH
# use_ibrav = true                                 # default: false
# use_symwan = true                                # default: false

[pw2wan]
write_unk = true                           # UNK files for Wannier-function plots
```

| Key | Description | Default |
|-----|-------------|---------|
| `pseudo_dir` | Directory containing the pseudopotential files | required |
| `scf_k_resolution` | k-point spacing of the SCF mesh in 1/Å (`round(|b_i| / resolution)`, cif2cell's rule) | required |
| `degauss` | Gaussian smearing in Ry | required |
| `pw2wan.write_unk` | Write UNK files (`true`/`false`, or `".true."`/`".false."`); when false, `wannier_plot` is switched off in `pwscf.win` | required |
| `pw2wan.wannier_plot_supercell` | Copied into `pw2wan.in` | not written |
| `pp_list_path` | Pseudopotential table: the file name of a bundled table, or a path | `pp_psl_rrkj.csv` |
| `cif2cell_path` | cif2cell executable (cif2cell reader only) | `cif2cell` on PATH |
| `use_ibrav` | Write the cell as QE's `ibrav` with `A`, `B`, `C`, `cosAB`, ... (see below) | `false` |
| `use_symwan` | Prepare the Wannier90 NSCF run for symWannier (see below) | `false` |

The NSCF (Wannier90) mesh is the SCF mesh clamped to 4..8 points per
direction. The fixed settings (`nbnd`, `conv_thr`, `dis_num_iter`,
`dis_froz_max = -200`, ...) are listed in `cif2qewan/workflow/builder.py`;
`dis_froz_max` is meant to be set after the NSCF run (recommended:
E_F + 1 eV to E_F + 3 eV), as `submit_all.sh` does.

### Pseudopotential tables

The tables installed with the package are selected by file name:

| Table | Pseudopotentials |
|-------|------------------|
| `pp_psl_rrkj.csv` (default) | PSLibrary ultrasoft (`rrkjus`), PBE |
| `nc-sr-05_pbe_standard_upf.csv`, `nc-sr-05_pbe_stringent_upf.csv` | PseudoDojo NC v0.5, PBE |
| `nc-sr-04_pbe_standard_upf.csv`, `nc-sr-04_pbe_stringent_upf.csv` | PseudoDojo NC v0.4, PBE |
| `nc-sr-04_pbesol_standard_upf.csv`, `nc-sr-04_pbesol_stringent_upf.csv` | PseudoDojo NC v0.4, PBEsol |

With `--so` the fully relativistic file is used: `.pbe` -> `.rel-pbe`
(PSLibrary) and `_sr` -> `_fr` (PseudoDojo). The PseudoDojo tables are named
after the archives on [pseudo-dojo.org](https://www.pseudo-dojo.org/) and
list the files as `<element>_sr`, so rename `X.upf` of the scalar-relativistic
set to `X_sr.UPF` and, for `--so`, `X.upf` of the matching `nc-fr-04` set to
`X_fr.UPF` (v0.5 has no fully relativistic set; La and Lu have none in the
v0.4 standard sets). Their `ecutwfc` is twice the "high" cutoff hint of the
set (Ha -> Ry), `ecutrho` four times `ecutwfc`, and `nexclude` counts the
semicore bands below the projected shells; `tools/pseudodojo_table.py`
regenerates a table from the `upf` and `djrepo` archives.

A `pp_list_path` with a directory part is used as a path, so you can write
your own table. Columns: `atom`, `pp_file_name` (without `.UPF`),
`nexclude` (low-lying bands excluded from the Wannier fit), `orbitals`
(projections, letters of `s`, `p`, `d`, `f`), `ecutwfc`, `ecutrho` (Ry); an
empty `pp_file_name` marks an unsupported element.

```csv
atom,pp_file_name,nexclude,orbitals,ecutwfc,ecutrho
Fe,Fe.pbe-spn-rrkjus_psl.0.2.1,4,spd,64.0,782.0
O,O.pbe-n-rrkjus_psl.0.1,1,p,47.0,323.0
```

### Cell representation (`use_ibrav`)

By default the cell is written as `ibrav = 0` with the lattice vectors in
`CELL_PARAMETERS {alat}`. With `use_ibrav = true` the inputs use QE's
Bravais-lattice index instead, as mcif2qewan and cif2x do: the lattice type
is determined with spglib from the space group, the `&system` namelist gets
`ibrav` and `A`, `B`, `C`, `cosAB`, `cosAC`, `cosBC`, and no
`CELL_PARAMETERS` card is written. The structure is re-expressed in exactly
the vectors `pw.x` builds from these parameters (new fractional coordinates,
the crystal and its moments rigidly rotated), `pwscf.win` uses the same
vectors, and the SCF mesh is derived from them. The values written are `1`,
`2`, `-3`, `4`, `5`, `6`, `7`, `8`, `9`, `91`, `10`, `11`, `-12`, `-13` and
`14`, as defined in QE 7.x (`ibrav = -13` changed in QE 6.4.1, so QE 6.4.1 or
later is required). Moments that only break point symmetry leave the lattice
type unchanged; a magnetic supercell (antiferromagnet) keeps its own lower
lattice type.

### symWannier (`use_symwan`)

With `use_symwan = true` the Wannier90 NSCF run is prepared for
[symWannier](https://github.com/wannier-utils-dev/symWannier): `nscf.in` has
no `nosym` and uses `K_POINTS {automatic}` on the Wannier90 mesh, so `pw.x`
computes the irreducible k points only, and `pw2wan.in` gets
`irr_bz = .true.` (QE 7.3 or later). Run `symwannier expand pwscf` after
`pw2wannier90.x` to produce the full-mesh `pwscf.mmn`, `pwscf.amn` and
`pwscf.eig`, then `wannier90.x pwscf` as usual; `pwscf.win` is unchanged.
UNK files are written for the irreducible points only and are not expanded,
so `wannier_plot` is switched off (with a warning when `write_unk` is true).

## Workflow

`submit_all.sh` runs the whole chain (set `ESPRESSO_DIR`, `WANNIER90_DIR`,
`TOML_FILE`, `MPI_PREFIX` and, for symWannier, `USE_SYMWAN=1` at its top).
Step by step:

```bash
cif2qewan structure.cif cif2qewan.toml
mpirun -n 16 pw.x < scf.in > scf.out
cp -r work check_wannier/ && cp -r work band/
mpirun -n 16 pw.x < nscf.in > nscf.out
wannier90.x -pp pwscf
mpirun -n 16 pw2wannier90.x < pw2wan.in > pw2wan.out
# symwannier expand pwscf                     # only with use_symwan = true
# set dis_froz_max in pwscf.win (E_F + 1 eV; E_F from nscf.out)
wannier90.x pwscf

# convergence check: DFT vs Wannier energies on a shifted k mesh
(cd check_wannier && mpirun -n 16 pw.x < nscf.in > nscf.out)
wannier_conv -e 5.0 -o ./ -i ./check_wannier/nscf.out   # writes CONV_5.0

# band structure: DFT bands and the Wannier interpolation in one plot
(cd band && mpirun -n 16 pw.x < nscf.in > nscf.out && mpirun -n 16 bands.x < band.in > band.out)
band_comp -o ./                                          # band_compare.png / .eps
```

`wannier_conv` prints the average and maximum difference between the DFT and
the Wannier-interpolated energies up to `-e` eV above the Fermi level.
`band_comp` reads `scf.out`, `pwscf.win`, `pwscf_band.dat` (Wannier90 with
`bands_plot`) and `band/bands.out.gnu` (`bands.x`) from the current
directory.

## Python API

The command is a thin layer over a pipeline of typed models: a reader returns
a `NormalizedStructure` (lattice in angstrom, fractional coordinates,
Cartesian moments), `build_plan` turns it and a `Config` into a
`CalculationPlan` holding every input as a model, and `render_plan` /
`write_rendered` turn the plan into files.

```python
from cif2qewan.config import Config
from cif2qewan.structure.readers import Cif2cellReader, PymatgenReader
from cif2qewan.workflow.builder import build_plan, plan_from_cif2cell_output
from cif2qewan.workflow.render import render_plan
from cif2qewan.workflow.output import write_rendered

config = Config.from_toml("cif2qewan.toml", so=True, mag=True)

output = Cif2cellReader(config.cif2cell_path, config.scf_k_resolution).run("Fe.cif")
plan = plan_from_cif2cell_output(output, config)   # the cif2cell path
plan = build_plan(PymatgenReader().read("Fe.cif"), config)   # the pymatgen path

plan.nscf.system.entries["nbnd"], plan.wannier90.num_wann
write_rendered(render_plan(plan), "run")           # {relative path: text} -> files
```

`plan.files()` maps the file names to their models (`PwInput`,
`NamelistInput`, `Wannier90Input`). Errors raised on purpose derive from
`cif2qewan.exceptions.Cif2qewanError`.

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
