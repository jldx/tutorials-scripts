#!/bin/bash

################################################################################
# Post-Process Molecular Dynamics Trajectories using GROMACS
################################################################################
#
# DESCRIPTION:
#   This script performs comprehensive post-processing of GROMACS molecular 
#   dynamics trajectories that were simulated in chunks. It handles:
#   - Periodic boundary condition (PBC) corrections (no-jump, whole, cluster)
#   - Molecular centering and translation
#   - Trajectory concatenation
#   - Fitting to reference structure
#
# USAGE:
#   ./post_process.sh
#
# DEPENDENCIES:
#   - GROMACS (gmx_mpi or gmx command)
#   - Index file (.ndx) with group definitions
#   - GROMACS structure (.tpr) and trajectory (.xtc) files
#
# INPUT FILES REQUIRED:
#   - ${system_name}_${ts}.tpr   (topology file for each chunk)
#   - ${system_name}_${ts}.xtc   (trajectory file for each chunk)
#   - index.ndx                  (GROMACS index file with group definitions)
#
# OUTPUT FILES GENERATED:
#   - ${system_name}_${ts}_processed.xtc  (processed individual chunks)
#
################################################################################
#
# /!\ Modify the groups numbers according to the index file.
#
################################################################################

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directory containing raw simulation files
raw_dir='.'

# Directory where processed trajectories will be saved
chunks_dir='.'

# Overall trajectory length in nanoseconds (format: start_end, used for naming)
traj_length='0_100' # In ns

# Array of simulation time chunks to process (format: start_end in ns)
traj_steps=(
    0_50
    50_100
)

# System name (used as prefix for all input/output files)
system_name='CH4_H2O_rest'

# String of paths to the processed subtrajectories (built during loop)
files=''

ts0 = ${traj_steps[0]}

# ============================================================================
# MAIN PROCESSING LOOP: Process each trajectory chunk
# ============================================================================

for ts in ${traj_steps[@]}; do
    echo "Processing chunk: ${ts} ns"

    # ======================================================================
    # Step 1: Extract starting time from chunk notation (e.g., 0 from "0_50")
    # ======================================================================
    # This is needed to reassign correct timestamps in the output trajectory
    IFS='_' read -r -a array <<< ${ts}
    t0="${array[0]}000"  # Convert ns to ps (multiply by 1000)
    
    # ======================================================================
    # Step 2: Remove jumps due to periodic boundary conditions
    # ======================================================================
    # nojump: removes PBC jumps to ensure molecules don't teleport
    #   Input:  ${system_name}_${ts}.tpr, ${system_name}_${ts}.xtc
    #   Output: tmp_nojump.xtc
    #   Group:  system
    echo "0" | gmx_mpi trjconv \
        -s ${raw_dir}/${system_name}_${ts0}.tpr \
        -f ${raw_dir}/${system_name}_${ts0}.xtc \
        -o tmp_nojump.xtc \
        -pbc nojump \
        -n index.ndx

    # ======================================================================
    # Step 3: Make molecules whole (reassemble split molecules)
    # ======================================================================
    # whole: ensures molecules broken across PBC are reassembled
    #   Input:  tmp_nojump.xtc
    #   Output: tmp_whole.xtc
    #   Group:  system
    echo "0" | gmx_mpi trjconv \
        -s ${raw_dir}/${system_name}_${ts0}.tpr \
        -f tmp_nojump.xtc \
        -o tmp_whole.xtc \
        -pbc whole \
        -n index.ndx

    # ======================================================================
    # Step 4: Cluster molecules around reference structure
    # ======================================================================
    # cluster: minimizes distance by clustering molecules around reference
    #   Input:  tmp_whole.xtc
    #   Output: tmp_cluster.xtc
    #   Group 1: protein/reference structure - cluster around this)
    #   Group 2: system
    echo "1 0" | gmx_mpi trjconv \
        -s ${raw_dir}/${system_name}_${ts0}.tpr \
        -f tmp_whole.xtc \
        -o tmp_cluster.xtc \
        -pbc cluster \
        -n index.ndx

    # ======================================================================
    # Step 5: Center membrane and consolidate box
    # ======================================================================
    # -pbc mol:      keep molecules intact
    # -center:       center selected group in the box
    # -ur compact:   use non-rectangular unit cell
    # -t0:           set initial time (corrects time after concatenation)
    #   Input:  tmp_cluster.xtc
    #   Output: tmp_center.xtc
    #   Group 1: membrane
    #   Group 2: system
    echo "4 0" | gmx_mpi trjconv \
        -s ${raw_dir}/${system_name}_${ts0}.tpr \
        -f tmp_cluster.xtc \
        -o ${chunks_dir}/tmp_center.xtc \
        -pbc mol \
        -center \
        -ur compact \
        -n index.ndx \
        -t0 ${t0}

    # ============================================================================
    # Step 6: Align to reference structure
    # ============================================================================
    # fit rot+trans: performs rotation and translation fit to reference
    #   Input:  ${system_name}_${traj_length}_tmp.xtc
    #   Output: ${system_name}_${traj_length}_processed.xtc
    #   Group 1: what to fit to
    #   Group 2: system
    echo "3 4" | gmx_mpi trjconv \
        -s ${raw_dir}/${system_name}_${ts0}.tpr \
        -f tmp_center.xtc \
        -o ${chunks_dir}/${system_name}_${ts}_processed.xtc \
        -fit rot+trans \
        -n ${raw_dir}/index.ndx

    # ======================================================================
    # Cleanup: Remove temporary intermediate files
    # ======================================================================
    rm tmp_nojump.xtc tmp_whole.xtc tmp_cluster.xtc tmp_mol.xtc tmp_center.xtc
    
done

echo "Post-processing complete!"