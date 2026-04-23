#!/usr/bin/env python3
"""
GROMACS MD Production Chunk Job Generator

This script generates SLURM job submission scripts for running GROMACS
molecular dynamics (MD) production simulations in time chunks on HPC clusters.

The script divides a long simulation into smaller, manageable chunks and creates
a chain of SLURM jobs where each job automatically submits the next one upon
successful completion. This approach manages walltime limitations on HPC
clusters while maintaining simulation continuity.

Typical workflow:
  - Chunk 1: 0–50 ns   (starting from final equilibration structure)
  - Chunk 2: 50–100 ns (starting from chunk 1 output)
  - Chunk 3: 100–150 ns (starting from chunk 2 output)
  - ... and so on

Each chunk automatically submits the next upon successful completion.

Usage:
    1. Edit configuration below (system name, paths, time windows)
    2. Run: python3 generate_gromacs_production_chunks.py
    3. Submit first job: sbatch run_prod_0_50.slurm
"""

# ==============================================================================
# CONFIGURATION PARAMETERS
# ==============================================================================

# System name (used as prefix for output files)
system = "<system_name>"

# Path to topology/reference structure
topology = "<path_to_topology_or_reference>"

# Simulation time parameters (nanoseconds)
current_step = 0          # Starting time of production (ns)
end_step = 1000           # Total production time to run (ns)
steps_per_job = 50        # Duration of each chunk/job (ns)

# Path to production MDP file
mdp_file = "<path_to_production.mdp>"

# Path to final equilibration output structure (.gro or .pdb)
last_equi_file = "<path_to_final_equilibration_output.gro>"

# ==============================================================================
# GENERATE PRODUCTION CHUNK SCRIPTS
# ==============================================================================

# Iterate through time chunks
for current_step in range(current_step, end_step, steps_per_job):
    next_step = current_step + steps_per_job
    previous_step = current_step - steps_per_job
    next_next_step = next_step + steps_per_job
    
    # Output prefix for this chunk
    chunk_name = f'{system}_{current_step}_{next_step}'
    
    # Determine input structure
    if current_step == 0:
        # First chunk: use equilibration output
        input_structure = last_equi_file
    else:
        # Subsequent chunks: use previous chunk's output
        input_structure = f'{system}_{previous_step}_{current_step}.gro'
    
    # Determine job chaining
    total_chunks = (end_step - current_step) // steps_per_job + 1
    chunk_number = (current_step // steps_per_job) + 1
    
    script_filename = f'run_prod_{current_step}_{next_step}.slurm'
    
    with open(script_filename, 'w') as f:
        f.write(f"""#!/bin/bash

[[HEADER, MODULES AND EXPORTED VARIABLES FOR SUBMISSION - (see stand-alone slurm submission scripts)]]

# ============================================================================

# GROMACS Production Run: {current_step}–{next_step} ns
# Chunk {chunk_number} of {total_chunks}

# Step 1: Pre-processing with grompp
#   -f: input MDP file (simulation parameters)
#   -o: output TPR file (portable binary run input)
#   -c: input structure file (starting coordinates)
#   -r: reference structure (for restraints, if any)
#   -p: topology file
#   -n: index file (atom groups)
#   -maxwarn: warning tolerance

gmx grompp -f {mdp_file} -o {chunk_name}.tpr -c {input_structure} -r {topology} -p topol.top -n index.ndx -maxwarn 1

# Step 2: Run the simulation (mdrun)
#   Add partition-specific flags (GPU/CPU) as needed

if gmx mdrun -v -deffnm {chunk_name} [[ SPECIFICATION FOR SUBMISSION ]]; then
    # Current chunk completed successfully
    # Submit next chunk
    sbatch run_prod_{next_step}_{next_next_step}.slurm
fi
""")
    
    print(f"✓ Generated: {script_filename}")

print(f"\n✅ Generated {(end_step - current_step) // steps_per_job + 1} production chunk scripts.")
print(f"Submit the first chunk with: sbatch run_prod_0_{steps_per_job}.slurm")
    