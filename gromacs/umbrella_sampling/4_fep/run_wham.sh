#!/bin/bash
# Weighted Histogram Analysis Method (WHAM) Execution Script
# Author: Julie LEDOUX
# Organization: CNRS
# Date: April 23, 2026
#
# Purpose:
# This script prepares and runs the Weighted Histogram Analysis Method (WHAM)
# to calculate the potential of mean force (PMF) from umbrella sampling simulations.
# It collects the necessary input files (TPR, pullx.xvg, and pullf.xvg) for each window,
# organizes them, and executes the GROMACS WHAM tool.
#
# Usage:
# 1. Ensure all umbrella sampling output files (conf*_umb.tpr, conf*_pullx.xvg, conf*_pullf.xvg)
#    are available in the specified directory.
# 2. Run this script to prepare the input files and execute WHAM.

# List of conformations (windows) to process.
# These are the frame numbers or identifiers for each window.
# Each value corresponds to an umbrella sampling simulation output.
CONF=(0   47  241  331  390  502  683  717  823 1027 1131 1239 1301 1347 1386 1473 1546)

# Create a directory to store all WHAM input files
mkdir -p data
cd data || exit 1

# Loop over each conformation
for c in ${CONF[@]}; do
    # Copy the TPR file for the current window
    # TPR files contain the window definitions and simulation parameters
    cp ../3_umbrella_sampling/b_umbrella_sampling/conf$c/conf$c_umb.tpr .

    # Copy the pullx.xvg file for the current window
    # pullx.xvg files contain the collective variable (CV) trajectories
    cp ../3_umbrella_sampling/b_umbrella_sampling/conf$c/conf$c_umb_pullx.xvg .

    # Copy the pullf.xvg file for the current window
    # pullf.xvg files contain the force data for the CV
    cp ../3_umbrella_sampling/b_umbrella_sampling/conf$c/conf$c_umb_pullf.xvg .
done

# Create a list of TPR files for WHAM input
# This file is required by GROMACS WHAM to identify the window definitions
echo "Creating list of TPR files..."
ls conf*.tpr > tpr.dat

# Create a list of pullx.xvg files for WHAM input
# This file is required by GROMACS WHAM to identify the CV trajectories
echo "Creating list of pullx.xvg files..."
ls conf*_pullx.xvg > pullx.dat

# Optional: Create a list of pullf.xvg files for WHAM input
# This file is required if force data is needed for WHAM analysis
echo "Creating list of pullf.xvg files..."
ls conf*_pullf.xvg > pullf.dat

# Verify that the files are found and listed correctly
echo "Verifying TPR files..."
cat tpr.dat

echo "Verifying pullx.xvg files..."
cat pullx.dat

# Execute GROMACS WHAM
# -it: Input TPR file list
# -ix: Input pullx.xvg file list
# -xvg none: Disable XMGrace output
# -nBootstrap 200: Number of bootstrap samples for error estimation
# -v: Verbose output
echo "Running WHAM..."
gmx wham -it tpr.dat -ix pullx.dat -xvg none -nBootstrap 200 -v

