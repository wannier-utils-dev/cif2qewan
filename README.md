# cif2qewan
cif2qewan.py is a simple python script to create quantum-ESPRESSO and wannier90 inputs from cif files.

## Usage ######################################
  1. Prepare cif2cell. (See for details, https://sourceforge.net/projects/cif2cell/)
  
  2. Prepare pseudopotentials in PSLibrary.
  
  3. Download or clone the github repository, e.g.
      % git clone https://github.com/wannier-utils-dev/cif2qewan
  
  
  4. Edit cif2cell_path and pseudo_dir in cif2qewan.py.
  
  5. Run.
      % pw.x < scf.in > scf.out
      % pw.x < nscf.in > nscf.out
      % wannier90.x -pp pwscf
      % pw2wannier90.x < pw2wan.in

  6. Edit dis_froz_max in pwscf.win. Recommended value is around EF+1eV ~ EF+3eV.
  
  7. Wannierize.
  > % wannier90.x pwscf
