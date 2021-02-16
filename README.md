# cif2qewan
cif2qewan.py is a simple python script to create quantum-ESPRESSO and wannier90 inputs from cif files.

## Usage ######################################
  1. Prepare cif2cell. (See for details, https://sourceforge.net/projects/cif2cell/)
  
  2. Prepare pseudopotentials in PSLibrary.
  
  3. Download or clone the github repository, e.g.
  
      % git clone https://github.com/wannier-utils-dev/cif2qewan
  
  
  4. Edit cif2cell_path and pseudo_dir in cif2qewan.py.
  
  5. Run.
      % python cif2qewan.py **.cif
  
      % pw.x < scf.in > scf.out
      
      % pw.x < nscf.in > nscf.out
      
      % wannier90.x -pp pwscf
      
      % pw2wannier90.x < pw2wan.in

  6. Edit dis_froz_max in pwscf.win. Recommended value is around EF+1eV ~ EF+3eV.

  7. Wannierize.
  
      % wannier90.x pwscf


## Compare band structures of DFT and wannier90 #####
cif2qewan.py prepares band calculation input files in directory "band".

     % cd band

     % pw.x < ../scf.in > scf.out

	 % pw.x < nscf.in > nscf.out

	 % bands.x < band.in > band.out

	 % cd ..

	 % python band_comp.py

Then, you can get the band structure plot of DFT and wannier90.

## Compare band energy of DFT and wannier90 #####
cif2qewan.py prepares nscf input file for energy diffierence.
Wannier90 Hamiltonian should reproduce the band energy on the kmesh for wannierzation. (For example, 8x8x8 mesh including gamma point (8 8 8 0 0 0 in QE expression).)
Here, the code checks the energy difference of DFT and wannier90 on the shifted kmesh. (c.f., 8 8 8 1 1 1 in QE expression))

	% cd check_wannier

	% pw.x < ../scf.in > scf.out

	% pw.x < nscf.in > nscf.out

	% cd ..

	% python wannier_conv.py

	% cat check_wannier/CONV

wannier_conv.py calculates the energy differences and outputs the result in check_wannier/CONV.
