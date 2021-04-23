#!/bin/bash

MPI_PREFIX="mpirun -n 16"
ESPRESSO_DIR=/path/to/espresso_dir
WANNIER90_DIR=/path/to/wannier90_dir
CIF2QEWAN_DIR=/path/to/cif2qewan_dir
TOML_FILE=/path/to/cif2qewan.toml

python $CIF2QEWAN_DIR/cif2qewan.py *.cif $TOML_FILE

$MPI_PREFIX $ESPRESSO_DIR/bin/pw.x < scf.in > scf.out
cp -r work check_wannier/
cp -r work band/

$MPI_PREFIX $ESPRESSO_DIR/bin/pw.x < nscf.in > nscf.out &&
$MPI_PREFIX $WANNIER90_DIR/wannier90.x -pp pwscf &&
$MPI_PREFIX $ESPRESSO_DIR/bin/pw2wannier90.x < pw2wan.in > pw2wan.out &&
rm -r work

# set dis_froz_max = EF+1eV
ef=$(grep Fermi nscf.out | cut -c27-35)
ef1=$(bc -l <<< "$ef + 1")
sed -i "s/dis_froz_max .*/dis_froz_max = $ef1/g" pwscf.win

$MPI_PREFIX $WANNIER90_DIR/wannier90.x pwscf

cd check_wannier &&
$MPI_PREFIX $ESPRESSO_DIR/bin/pw.x < nscf.in > nscf.out &&
rm -r work &&
cd ../
python $CIF2QEWAN_DIR/wannier_conv.py

cd band &&
$MPI_PREFIX $ESPRESSO_DIR/bin/pw.x < nscf.in > nscf.out &&
$MPI_PREFIX $ESPRESSO_DIR/bin/bands.x < band.in > band.out &&
rm -r work &&
cd ../

python $CIF2QEWAN_DIR/band_comp.py
