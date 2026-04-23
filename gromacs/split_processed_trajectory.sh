#!/bin/bash

################################################################################
# Split Post-Processed Trajectories into Smaller Sub-Chunks
################################################################################
#
# DESCRIPTION:
#   This script further subdivides already post-processed GROMACS molecular 
#   dynamics trajectories into smaller time windows. It is typically used after
#   post_process_trajectory.sh and serves multiple purposes:
#
#   - Analysis convenience: Smaller files are easier to handle and analyze
#   - Parallel processing: Each chunk can be analyzed independently
#   - Storage optimization: Compress or archive individual chunks
#   - Targeted analysis: Focus on specific time windows of interest
#   - Memory efficiency: Load/process small chunks in memory
#
# WORKFLOW CONTEXT:
#   1. Raw simulation output (.xtc from mdrun)
#           |
#           v
#   2. post_process_trajectory.sh
#      (PBC correction, centering, fitting → processed chunks)
#           |
#           v
#   3. split_trajectory.sh (THIS SCRIPT)
#      (Further subdivide into smaller time windows)
#           |
#           v
#   4. Analysis tools (GROMACS, VMD, custom Python scripts)
#
# USAGE:
#   ./split_trajectory.sh
#
# DEPENDENCIES:
#   - GROMACS (gmx_mpi or gmx command)
#   - Index file (.ndx) with group definitions
#   - GROMACS structure (.gro or .tpr) file
#   - Post-processed trajectory files (output from post_process_trajectory.sh)
#
# INPUT FILES REQUIRED:
#   - ${system}_${TRAJ_STEP}_processed.xtc  (from post_process_trajectory.sh)
#   - ${top}.gro or .tpr                     (topology/structure file)
#   - index.ndx                              (group definitions)
#
# OUTPUT FILES GENERATED:
#   - ${subchunks_dir}/${system}_${chunk}_${sub_start}_${sub_end}.xtc  (format: system_chunk_start_end.xtc)
#
# KEY PARAMETERS:
#   - TRAJ_STEP: Array of coarse trajectory chunks (e.g., "0_25", "25_50")
#   - chunk_size: Size of final output sub-chunks in nanoseconds
#
################################################################################

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directory where processed trajectories are saved
chunks_dir='.'

# Directory where the splited trajectories will be saved
subchunks_dir='.'

# Reference structure file (used for topology, box info, etc.)
# Can be .gro (structure), .tpr (binary topology), or .pdb
top="crystal_membrane_solv_ions.gro"

# System name (prefix of processed trajectory files)
# Expects files named: ${system}_${TRAJ_STEP}_processed.xtc
system='H5cm'

# Size of output sub-chunks in nanoseconds
# Smaller chunks = more files but easier to handle
# Typical values: 5 ns (frequent output), 10 ns (moderate), 50 ns (coarse)
chunk_size=5

# Array of coarse trajectory chunks from post_process_trajectory.sh
# Format: "start_end" in nanoseconds
# These are the input processed trajectory chunks
TRAJ_STEP=(
0_25
25_50
50_75
75_100
)

# ============================================================================
# MAIN PROCESSING LOOP
# ============================================================================
# For each coarse chunk, further subdivide into smaller sub-chunks

for chunk in "${TRAJ_STEP[@]}"; do
    
    # ==================================================================
    # Parse coarse chunk notation (e.g., "0_25" → start=0, end=25)
    # ==================================================================
    # This tells us the time span of the current processed trajectory
    IFS='_' read -r start end <<< "$chunk"

    # ==================================================================
    # Loop through coarse chunk in steps of chunk_size
    # ==================================================================
    # For each coarse chunk, create multiple sub-chunks
    # Output naming: ${system}_${chunk}_${sub_start}_${sub_end}
    # Example: For chunk "0_25" with chunk_size=5 produces:
    #   H5cm_0_25_0_5.xtc   (frames 0-5 ns within chunk)
    #   H5cm_0_25_5_10.xtc  (frames 5-10 ns within chunk)
    #   H5cm_0_25_10_15.xtc (frames 10-15 ns within chunk)
    #   H5cm_0_25_15_20.xtc (frames 15-20 ns within chunk)
    #   H5cm_0_25_20_25.xtc (frames 20-25 ns within chunk)
    for ((i = 0; i < end - start; i += chunk_size)); do
    
        # ==============================================================
        # Define input trajectory file (from post_process_trajectory.sh)
        # ==============================================================
        # File format: ${system}_${chunk}_processed.xtc
        # Example: H5cm_0_25_processed.xtc
        traj=${chunks_dir}/${system}_${chunk}_processed.xtc
        
        # ==============================================================
        # Define extraction time window (in nanoseconds, relative to chunk)
        # ==============================================================
        # begin_traj: where to start extracting within the coarse chunk (ns)
        # end_traj:   where to stop extracting within the coarse chunk (ns)
        # These are offsets relative to the start of the current chunk
        begin_traj=$i
        end_traj=$((i + chunk_size))
    
        # ==============================================================
        # Calculate sub-chunk start and end for output file naming
        # ==============================================================
        # Convert relative offsets to absolute time within the coarse chunk
        # Example: For chunk "0_25", i=10 gives:
        #   sub_start = 10 ns
        #   sub_end = 15 ns
        #   Output: H5cm_0_25_10_15.xtc
        sub_start=$i
        sub_end=$((i + chunk_size))
    
        # ==============================================================
        # Safety check: Ensure end boundary doesn't exceed coarse chunk
        # ==============================================================
        # If the calculated end exceeds the coarse chunk boundary,
        # truncate to the actual end of the coarse chunk
        if [ "$sub_end" -gt $((end - start)) ]; then
            sub_end=$((end - start))
        fi
        
        # ==============================================================
        # Define output file path
        # ==============================================================
        # Create sub-chunk trajectory in "chunks/" directory
        # Naming convention: ${system}_${chunk}_${sub_start}_${sub_end}.xtc
        # Example output: chunks/H5cm_0_25_0_5.xtc
        chunk_traj=${subchunks_dir}/${system}_${chunk}_${sub_start}_${sub_end}.xtc

        # ==============================================================
        # Extract sub-chunk using gmx trjconv
        # ==============================================================
        # Parameters:
        #   -b: begin time (ns relative to coarse chunk start)
        #   -e: end time (ns relative to coarse chunk start)
        #   -tu ns: time unit in nanoseconds
        #   -s: structure/topology file (for PBC information)
        #   -f: input trajectory file
        #   -o: output trajectory file
        #   -n: index file (group definitions)
        #
        # Output: Group 0 (typically "System" or "All" from index.ndx)
        # Note: Pipe "0" to select group (non-interactive mode)
        echo "0" | gmx_mpi trjconv \
            -b ${begin_traj} \
            -e ${end_traj} \
            -tu ns \
            -s ${top} \
            -f ${traj} \
            -o ${chunk_traj} \
            -n index.ndx
        
    done
done

echo "Splitting complete!"