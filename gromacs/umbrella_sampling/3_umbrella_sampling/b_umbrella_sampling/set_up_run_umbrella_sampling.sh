#!/bin/bash
# Umbrella Sampling Simulation Submission Script
#
# Purpose:
# This script automates the submission of umbrella sampling simulations
# for each window (conformation) extracted from a previous equilibration step.
# It creates directories for each window, copies the necessary input files,
# customizes the submission script, and submits the job to a queue system (e.g., PBS).
#
# Usage:
# 1. Ensure all equilibrated conformations (e.g., conf*_eq.gro) are in the ../a_equilibration/conf*/ directory.
# 2. Ensure all necessary simulation files (e.g., .mdp, .top, .ndx, forcefield, submission scripts) are in the raw_files/ directory.
# 3. Run this script to submit umbrella sampling jobs for each window.

# List of conformations (windows) to process.
# These are the frame numbers or identifiers for each window.
# Each value corresponds to an equilibrated conformation file (e.g., conf47_eq.gro, conf241_eq.gro, etc.).
CONF=(0   47  241  331  390  502  683  717  823 1027 1131 1239 1301 1347 1386 1473 1546)

# Loop over each conformation
for c in ${CONF[@]}; do
    # Create a dedicated directory for the current window
    # This keeps files organized and avoids conflicts between simulations
    echo "Setting up umbrella sampling simulation for window $c..."
    mkdir -p conf$c

    # Copy the equilibrated conformation file (e.g., conf${c}_eq.gro) to the window directory
    # This file contains the starting structure for the umbrella sampling simulation
    cp -r ../a_equilibration/conf$c/conf${c}_eq.gro conf$c/

    # Copy all necessary files for umbrella sampling simulation
    # This includes:
    # - .mdp files (simulation parameters)
    # - Forcefield files
    # - index.ndx (index file for GROMACS)
    # - topol.top (topology file)
    # - Submission script (e.g., run_umbrella.pbs)
    cp -r raw_files/* conf$c/

    # Navigate into the window directory
    cd conf$c || exit 1

    # Customize the submission script for the current window
    # Replace 'conf0' with 'conf$c' in the submission script (e.g., run_umbrella.pbs)
    # This ensures:
    # - Unique job names for each window
    # - Correct input files for GROMACS commands (grompp, mdrun)
    echo "Customizing submission script for window $c..."
    sed -i -e "s/conf0/conf$c/g" run_umbrella.pbs

    # Submit the job to the queue system (e.g., PBS)
    # This sends the umbrella sampling job to the cluster for execution
    echo "Submitting umbrella sampling job for window $c..."
    qsub run_umbrella.pbs

    # Return to the parent directory
    cd ..
done

echo "Umbrella sampling simulation submission complete for all windows."