#!/bin/bash
#
# Complete Quantum ESPRESSO + Wannier90 workflow from a CIF file: input
# generation, SCF, NSCF, (symWannier,) Wannier90, the convergence check and
# the band-structure comparison. The script stops at the first failing
# step (set -e); every command writes its output next to its input.
#
# Set the paths below, put the CIF file into an empty directory and run the
# script there. Requires cif2qewan installed in the active Python
# environment (the cif2qewan, wannier_conv and band_comp commands).

set -euo pipefail

# =============================================================================
# Configuration
# =============================================================================

MPI_PREFIX="mpirun -n 16"               # MPI launcher and number of processes
ESPRESSO_DIR=/path/to/espresso_dir      # Quantum ESPRESSO installation directory
WANNIER90_DIR=/path/to/wannier90_dir    # Wannier90 installation directory
TOML_FILE=/path/to/cif2qewan.toml       # Configuration file
USE_SYMWAN=0                            # 1 if use_symwan = true in the TOML file (symwannier installed)

PW="$MPI_PREFIX $ESPRESSO_DIR/bin/pw.x"

step() { echo "== $*"; }
trap 'echo "submit_all.sh: step failed (line $LINENO); see the .out file of the last command" >&2' ERR

# =============================================================================
# Step 1: input files from the CIF file
# =============================================================================

step "Step 1: generating the input files"
cif2qewan ./*.cif "$TOML_FILE"

# =============================================================================
# Step 2: SCF
# =============================================================================

step "Step 2: SCF"
$PW < scf.in > scf.out
cp -r work check_wannier/
cp -r work band/

# =============================================================================
# Step 3: NSCF for Wannier90, Wannier90 preprocessing, pw2wannier90
# =============================================================================

step "Step 3: NSCF for Wannier90"
$PW < nscf.in > nscf.out
step "Wannier90 preprocessing"
$MPI_PREFIX "$WANNIER90_DIR/wannier90.x" -pp pwscf
step "pw2wannier90"
$MPI_PREFIX "$ESPRESSO_DIR/bin/pw2wannier90.x" < pw2wan.in > pw2wan.out
if [ "$USE_SYMWAN" = "1" ]; then
  # irr_bz = .true.: expand Mmn, Amn and Eig from the irreducible k points
  step "symWannier: expanding to the full k mesh"
  symwannier expand pwscf
fi
rm -r work

# =============================================================================
# Step 4: frozen window and Wannier90
# =============================================================================

step "Step 4: Wannier90"
# dis_froz_max = E_F + 1 eV (recommended: E_F + 1 eV to E_F + 3 eV)
ef=$(grep Fermi nscf.out | tail -1 | cut -c27-35)
ef1=$(bc -l <<< "$ef + 1")
echo "Fermi energy $ef eV, dis_froz_max = $ef1 eV"
sed -i "s/dis_froz_max .*/dis_froz_max = $ef1/" pwscf.win
$MPI_PREFIX "$WANNIER90_DIR/wannier90.x" pwscf

# =============================================================================
# Step 5: convergence check on a shifted k mesh
# =============================================================================

step "Step 5: convergence check"
(cd check_wannier && $PW < nscf.in > nscf.out && rm -r work)
wannier_conv -e 5.0 -o ./ -i ./check_wannier/nscf.out

# =============================================================================
# Step 6: band structure
# =============================================================================

step "Step 6: band structure"
(
  cd band
  $PW < nscf.in > nscf.out
  $MPI_PREFIX "$ESPRESSO_DIR/bin/bands.x" < band.in > band.out
  rm -r work
)

# =============================================================================
# Step 7: DFT vs Wannier90 band structures
# =============================================================================

step "Step 7: band-structure comparison"
band_comp -o ./

echo "Workflow finished:"
echo "  CONV_5.0               Wannier90 convergence results"
echo "  band_compare.png/.eps  band-structure comparison plots"
