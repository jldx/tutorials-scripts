#!/bin/bash
# ============================================================
# GROMACS Trajectory Coordinate Extraction
# ============================================================
# Usage:
#   bash extract_coords.sh <topology> <trajectory> <index> <output>
#
# Arguments:
#   $1  topology    – .tpr run file (or .gro/.pdb); provides atom types,
#                     bonds, and box dimensions
#   $2  trajectory  – .xtc or .trr file; contains per-frame coordinates
#   $3  index       – .ndx index file; defines the atom group to extract
#   $4  output      – output name (without .xvg); time-series coordinates
#                     will be written to <output>.xvg
#
# Example:
#   bash extract_coords.sh md.tpr md.xtc index.ndx protein_coords
# ============================================================

# --- Argument validation ------------------------------------
if [ "$#" -ne 4 ]; then
    echo "Error: expected 4 arguments, got $#."
    echo "Usage: $0 <topology> <trajectory> <index> <output>"
    exit 1
fi

TOPOLOGY="$1"
TRAJECTORY="$2"
INDEX="$3"
OUTPUT="$4"

# Check that input files exist before calling GROMACS
for FILE in "$TOPOLOGY" "$TRAJECTORY" "$INDEX"; do
    if [ ! -f "$FILE" ]; then
        echo "Error: file not found: $FILE"
        exit 1
    fi
done

# --- Run ---------------------------------------------------
gmx_mpi traj \
    -s  "$TOPOLOGY"   \  # topology / run input
    -f  "$TRAJECTORY" \  # trajectory
    -n  "$INDEX"      \  # index file (group selected interactively at runtime)
    -ox "$OUTPUT".xvg    # coordinate output (time [ps], x y z [nm] per atom)