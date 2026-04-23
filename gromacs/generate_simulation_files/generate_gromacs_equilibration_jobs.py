#!/usr/bin/env python3
"""
GROMACS MD Equilibration Stage Job Generator

This script generates SLURM job submission scripts for a series of GROMACS
equilibration stages, each building upon the previous with decreasing
position restraints.

Typical equilibration order:
  - mini         : steepest descent energy minimization
  - npt_0        : NPT with strong restraints (solute & membrane)
  - npt_1..npt_4 : intermediate NPT stages with progressively weaker restraints
  - npt_canonical: final NPT without restraints (true canonical ensemble)

Each stage automatically submits the next stage upon successful completion.

Usage:
    1. Edit configuration below (system name, paths)
    2. Run: python3 generate_gromacs_equilibration_jobs.py
    3. Submit first job: sbatch equilibrate_01_mini.slurm
"""

# ==============================================================================
# CONFIGURATION PARAMETERS
# ==============================================================================

# System name (used as prefix for output files)
system = "<system_name>"

# Path to initial structure file (input for minimization)
initial_structure = "<path_to_initial_structure.gro>"

# Path to topology/reference structure
topology = "<path_to_topology_or_reference>"

# Directory containing MDP files (mini.mdp, npt_0.mdp, ..., npt_canonical.mdp)
mdp_base_path = "./"

# ==============================================================================
# EQUILIBRATION STAGES
# ==============================================================================

equilibration_stages = [
    ("mini", "mini.mdp"),
    ("npt_0", "npt_0.mdp"),
    ("npt_1", "npt_1.mdp"),
    ("npt_2", "npt_2.mdp"),
    ("npt_3", "npt_3.mdp"),
    ("npt_4", "npt_4.mdp"),
    ("npt_canonical", "npt_canonical.mdp"),
]

# ==============================================================================
# GENERATE SLURM JOB SCRIPTS
# ==============================================================================

total_stages = len(equilibration_stages)

for stage_idx, (stage_name, mdp_file) in enumerate(equilibration_stages):
    stage_number = stage_idx + 1
    
    # Input structure for this stage
    if stage_idx == 0:
        input_structure_file = initial_structure
    else:
        prev_stage_name = equilibration_stages[stage_idx - 1][0]
        input_structure_file = f"{system}_{prev_stage_name}"
    
    # Next stage (for job chaining)
    if stage_idx < total_stages - 1:
        next_stage_name = equilibration_stages[stage_idx + 1][0]
        next_script = f"equilibrate_{stage_idx + 2:02d}_{next_stage_name}.slurm"
        chain_command = f"sbatch {next_script}"
    else:
        chain_command = "echo 'All equilibration stages completed!'"
    
    output_prefix = f"{system}_{stage_name}"
    script_filename = f"equilibrate_{stage_number:02d}_{stage_name}.slurm"
    
    with open(script_filename, 'w') as f:
        f.write(f"""#!/bin/bash

[[HEADER, MODULES AND EXPORTED VARIABLES FOR SUBMISSION - (see stand-alone slurm submission scripts)]]

# ============================================================================

# GROMACS Equilibration Stage {stage_number}/{total_stages}: {stage_name}

# Step 1: Prepare the simulation (grompp)
gmx grompp -f {mdp_base_path}{mdp_file} -o {output_prefix}.tpr -c {input_structure_file} -r {topology} -p topol.top -n index.ndx -maxwarn 1

# Step 2: Run the simulation (mdrun)
if gmx mdrun -v -deffnm {output_prefix} [[ SPECIFICATION FOR SUBMISSION ]]; then
    # Current stage completed successfully
    # Chain to the next stage
    {chain_command}
fi
""")
    
    print(f"✓ Generated: {script_filename}")

print(f"\n✅ Generated {total_stages} equilibration job scripts.")
print(f"Submit the first job with: sbatch equilibrate_01_mini.slurm")
