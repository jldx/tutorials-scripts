#!/bin/bash

################################################################################
# LMOD VIBRATIONAL ANALYSIS USING AMBER SANDER
################################################################################
#
# Purpose:
#   This script performs Low Mode (LMOD) vibrational analysis using AMBER's
#   sander module. LMOD explores protein conformational space by generating
#   molecular dynamics trajectories along low-frequency normal modes, useful
#   for studying collective protein motions and essential dynamics.
#
# Workflow:
#   1. Generate AMBER topology files using LEaP (make_topology.leap)
#   2. Create sander input file (lmod-vib.in) with LMOD parameters
#   3. Execute sander with MPI parallelization
#   4. Generate output files (trajectories, conformer library)
#
# Output Files:
#   - lmod-vib_trajectory: Binary trajectory file of LMOD sampling
#   - lmod-vib_conflib: Conformer library (generated structures)
#   - lmod-vib.out: Logfile with LMOD iteration details
#   - <prefix>-lmod-vib.ipcrd: Final coordinates after LMOD
#
# Prerequisites:
#   - AMBER software suite installed (sander compiled with MPI)
#   - LEaP topology files (<prefix>.parm7, <prefix>.inpcrd)
#   - Run make_topology.leap first to generate topology
#
# LMOD Theory:
#   - Samples conformational space along low-frequency normal modes
#   - Explores energetically favorable conformations
#   - Generates ensemble of structures for analysis
#   - Useful for proteins with large-scale motions
#
# Usage:
#   ./sander_lmod.sh <prefix>
#   
#   Example:
#   ./sander_lmod.sh protein_structure
#
# Author: Essential Dynamics Tutorial
# Last Modified: 2026-04-27
#
################################################################################

# ===== COMMAND-LINE ARGUMENT PARSING =====
# Extract input/output prefix from command-line argument
# Usage: sander_lmod.sh <prefix>
# The prefix should be the base filename without .parm7 or .inpcrd extension
# Example: if your files are protein.parm7 and protein.inpcrd, use: protein
prefix=${1}

# Validate that prefix was provided
if [ -z "$prefix" ]; then
    echo "Usage: $0 <prefix>"
    echo ""
    echo "Arguments:"
    echo "  <prefix>: Base filename without extension (e.g., protein)"
    echo ""
    echo "Prerequisites:"
    echo "  - Run make_topology.leap first to generate <prefix>.parm7 and <prefix>.inpcrd"
    echo ""
    exit 1
fi

# Verify that topology files exist
if [ ! -f "${prefix}.parm7" ] || [ ! -f "${prefix}.inpcrd" ]; then
    echo "Error: Topology files not found!"
    echo "  Expected: ${prefix}.parm7 and ${prefix}.inpcrd"
    echo ""
    echo "Please run LEaP first:"
    echo "  tleap -f make_topology.leap"
    exit 1
fi


# ===== STEP 1: VERIFY TOPOLOGY FILES =====
# Note: Before running this script, ensure make_topology.leap has been executed
# with the correct <prefix> variable substituted:
#   Replace <prefix> in make_topology.leap with your structure name
#   Then execute: tleap -f make_topology.leap
# This generates:
#   - <prefix>.parm7: AMBER topology file
#   - <prefix>.inpcrd: AMBER coordinate/restart file


# ===== STEP 2: CREATE SANDER INPUT FILE FOR LMOD =====
# Generate the sander control file (lmod-vib.in) with LMOD-specific parameters
# This file controls the vibrational sampling via LMOD

cat > lmod-vib.in <<EOF
 Input for sander LMOD vibrational analysis
 &cntrl
   cut     = 999,          ! Non-bonded cutoff (999 = no cutoff, implicit solvent)
   rgbmax  = 999,          ! GB radius cutoff
   ntx     = 1,            ! Read input coordinates only (no velocities)
   irest   = 0,            ! No restart (new MD run)
   ipol    = 0,            ! No polarizable force field
   ntb     = 0,            ! No periodic boundary conditions (implicit solvent)
   igb     = 8,            ! GB implicit solvent model (GBNeck2)
   imin    = 1,            ! Minimization (required for LMOD)
   maxcyc  = 500,          ! Maximum minimization cycles
   ntmin   = 4,            ! Minimization method (LBFGS with line search)
   drms    = 1e-11,        ! Convergence criterion (very tight for LMOD)
   ntpr    = 1,            ! Print frequency (every step)
 /
 &lmod
   ! ===== LMOD ITERATION PARAMETERS =====
   number_lmod_iterations = 1,      ! Number of LMOD sampling iterations
   conflib_size = 20,               ! Maximum conformer library size
   total_low_modes = 20,            ! Total low modes to sample
   explored_low_modes = 20,         ! Low modes to explicitly explore
   number_free_rotrans_modes = 6,   ! Rigid body modes (6 = translation + rotation)
   lmod_minimize_grms = 0.5,        ! Minimization gradient threshold
   
   ! ===== EIGENVALUE/EIGENVECTOR CALCULATION =====
   arnoldi_dimension = 50,          ! Krylov subspace dimension (higher = more accurate)
   matrix_vector_product_method = 'forward',  ! Method for Hessian-vector products
   frequency_eigenvector_recalc = 1,          ! Recalculate eigenvectors per iteration
   
   ! ===== LMOD STEP PARAMETERS =====
   number_lmod_moves = 10,          ! Number of moves along low modes per iteration
   lmod_step_size_max = 0.5,        ! Maximum step size along normal modes (Å)
   lmod_step_size_min = 0.5,        ! Minimum step size along normal modes
   lmod_relax_grms = 1.0,           ! Relaxation gradient threshold (looser than minimize)
   
   ! ===== LIGAND PARAMETERS =====
   number_ligands = 0,              ! Number of ligands (0 for protein only)
   
   ! ===== OUTPUT FILE PARAMETERS =====
   conflib_filename = 'lmod-vib_conflib',        ! Conformer library output
   lmod_trajectory_filename = 'lmod-vib_trajectory',  ! Trajectory output
   lmod_job_title = 'lmod-vib',                       ! Job identifier
   
   ! ===== VERBOSITY AND OPTIMIZATION =====
   lmod_verbosity = 5,              ! Verbosity level (0-5, higher = more detail)
   lbfgs_memory_depth = 5,          ! L-BFGS memory depth for minimization
   xmin_method = 'TNCG',            ! Minimization method (Truncated Newton CG)
   xmin_verbosity = 1,              ! Minimization output verbosity
 /
EOF

cat "LMOD sander input file created: lmod-vib.in\n"


# ===== STEP 3: EXECUTE SANDER WITH LMOD =====
# Run AMBER sander with LMOD vibrational analysis
# Uses MPI for parallel execution if available

cat "Executing sander with LMOD vibrational analysis...\n"
cat "  Input topology: ${prefix}.parm7\n"
cat "  Input coordinates: ${prefix}.inpcrd\n"
cat "  Output trajectory: lmod-vib_trajectory\n"
cat "  Output conformer library: lmod-vib_conflib\n\n"

# Remove old conformer library if present
rm -f conflib.dat

# Execute sander with parameters:
#   -O: Overwrite output files
#   -i: Input control file
#   -o: Output/logfile
#   -c: Coordinate input (starting coordinates)
#   -ref: Reference structure (for distance restraints/analysis)
#   -p: Topology file
#   -r: Output restart/coordinate file
mpirun sander -O -i lmod-vib.in -o lmod-vib.out \
   -c ${prefix}.inpcrd -ref ${prefix}.inpcrd \
   -p ${prefix}.parm7 -r ${prefix}-lmod-vib.ipcrd

# Capture exit status
EXIT_STATUS=$?


# ===== RESULTS SUMMARY =====
# Report completion status and provide next steps

if [ $EXIT_STATUS -eq 0 ]; then
    echo ""
    echo "==============================================="
    echo "LMOD Vibrational Analysis Complete - Success"
    echo "==============================================="
    echo ""
    echo "Output files generated:"
    echo "  - lmod-vib.out: Logfile with analysis details"
    echo "  - lmod-vib_trajectory: Binary trajectory of LMOD sampling"
    echo "  - lmod-vib_conflib: Conformer library (generated structures)"
    echo "  - ${prefix}-lmod-vib.ipcrd: Final coordinates"
    echo ""
    echo "Next steps:"
    echo "  1. Review lmod-vib.out for LMOD iteration details"
    echo "  2. Convert trajectory to readable format for visualization"
    echo "  3. Analyze conformational ensemble from lmod-vib_conflib"
    echo "  4. Perform clustering or principal component analysis on structures"
    echo ""
else
    echo ""
    echo "==============================================="
    echo "Error: LMOD Analysis Failed"
    echo "==============================================="
    echo "Exit status: $EXIT_STATUS"
    echo "Check lmod-vib.out for error details"
    exit $EXIT_STATUS
fi
   
   
