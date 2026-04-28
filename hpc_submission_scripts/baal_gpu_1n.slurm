#!/bin/bash
################################################################################
# GROMACS MD Production Run on Baal - GPU Partition (Local Lab Cluster)
################################################################################
#
# BAAL Cluster Information (Local Lab Cluster - GPU):
# - Compute Nodes: 40-core Intel/AMD CPU + GPU accelerators (1–2 GPUs per node)
# - Submission System: PBS/Torque (not SLURM)
# - Network: InfiniBand (optional, for multi-node jobs)
# - Storage: Local working directory (/PBS_O_WORKDIR)
# - Queue: gpu_40c_1n (single node, GPU-enabled)
#
# Baal GPU Partition Optimization Tips:
# 1. Single node with GPU acceleration: typically 1–2 GPUs per compute node
# 2. For single-node: use PSM local/shared memory (no IB overhead)
# 3. GPU acceleration: offload non-bonded, PME, and bonded forces
# 4. MPI context sharing enabled for efficient CPU-GPU communication
# 5. Job runs in original directory (no temp/scratch cleanup needed)
# 6. Useful for accelerated MD simulations on local infrastructure
#
################################################################################

#PBS -S /bin/bash                               # Shell for script execution
#PBS -N [XXX]                                   # Job name (required - replace XXX)
#PBS -o [XXX].out                               # Standard output file
#PBS -e [XXX].err                               # Standard error file
#PBS -q gpu_40c_1n                              # Queue: GPU-enabled, 1 node, 40 cores
#PBS -l nodes=node075.cluster.lbt:ppn=40:ppn=40 # Request 1 node with 40 cores (change 075 for different node)
#PBS -l walltime=12:00:00                       # Max execution time (12 hours)
#PBS -A simlab_project                          # Project account code
#PBS -m abe                                     # Email notifications: abort, begin, end
#PBS -M ledoux@ibpc.fr                          # Email address for notifications
#PBS -l epilogue=/shared/scripts/ADMIN__epilogue-qsub.example  # Cluster epilogue script

################################################################################
# System Environment Setup (Baal GPU-specific)
################################################################################

# Get working directory and node information
WORKDIR=$PBS_O_WORKDIR
NUM_NODES=$(cat $PBS_NODEFILE|uniq|wc -l)

# Validate PBS environment
if [ ! -n "$PBS_O_HOME" ] || [ ! -n "$PBS_JOBID" ]; then
    echo "ERROR: Required PBS variables not defined. Contact cluster administrator."
    exit 1
else
    echo "✓ Using original directory for job execution"
fi

# Display job environment information
echo "Current directory:  $PBS_O_WORKDIR"
echo "Working directory:  $WORKDIR"
echo "Number of nodes:    $NUM_NODES"
echo "Number of cores:    $(cat $PBS_NODEFILE|wc -l)"
echo "Execution nodes:    $(cat $PBS_NODEFILE|uniq|while read line; do printf "%s" "$line "; done)"

# Change to working directory
cd $WORKDIR

################################################################################
# MPI and InfiniBand Configuration
################################################################################

# For single-node jobs: use shared memory and local communication
# This avoids unnecessary InfiniBand overhead for local message passing
if [ $NUM_NODES -eq 1 ]; then
    # Single-node configuration: efficient local communication
    export PSM_DEVICES=self,shm              # Use shared memory (PSM: Performance Scaled Messaging)
    export OMPI_MCA_mtl=^psm                 # Disable PSM at MTL level (use shared memory instead)
    export OMPI_MCA_btl=shm,self             # Use shared memory and self for communication
    echo "✓ Single-node mode: PSM optimized for shared memory communication"
else
    # Multi-node configuration: manage InfiniBand PSM contexts
    # Each node has a single IB card with max 16 PSM contexts
    # With >16 processes per node, contexts must be shared
    
    CIN=$(cat /proc/cpuinfo | grep -i processor | wc -l)  # Core count per node
    
    if [ $(($CIN/16)) -ge 2 ]; then
        # More cores than PSM contexts available: must share contexts
        PPN=$(grep $HOSTNAME $PBS_NODEFILE|wc -l)  # Processes per node
        
        if [ $CIN -eq 40 ]; then
            # 40-core node: share contexts among 4 processes
            export PSM_SHAREDCONTEXTS_MAX=$(($PPN/4))
        elif [ $CIN -eq 32 ]; then
            # 32-core node: share contexts among 2 processes
            export PSM_SHAREDCONTEXTS_MAX=$(($PPN/2))
        else
            echo "ERROR: This computing node (${CIN} cores) is not supported by this script"
            exit 1
        fi
        
        echo "✓ Multi-node mode: PSM_SHAREDCONTEXTS_MAX = $PSM_SHAREDCONTEXTS_MAX"
    else
        echo "✓ Multi-node mode: PSM contexts sufficient (no sharing needed)"
    fi
fi

################################################################################
# GPU Detection and Configuration
################################################################################

# Detect which GPUs are allocated to this job
# Returns GPU ID(s) for CUDA_VISIBLE_DEVICES or similar GPU configuration
function get_gpu-ids() {
    # If all node cores allocated to job, use all GPUs (typically 0,1)
    if [ $PBS_NUM_PPN -eq $(cat /proc/cpuinfo | grep -cE "^processor.*:") ]; then
        echo "0,1" && return
    fi

    # Try to find GPU assignment from cpuset (newer PBS systems)
    if [ -e /dev/cpuset/torque/$PBS_JOBID/cpus ]; then
        FILE="/dev/cpuset/torque/$PBS_JOBID/cpus"
    elif [ -e /dev/cpuset/torque/$PBS_JOBID/cpuset.cpus ]; then
        FILE="/dev/cpuset/torque/$PBS_JOBID/cpuset.cpus"
    else
        FILE=""
    fi
    
    # Determine GPU ID based on first CPU core assigned
    if [ -e $FILE ]; then
        if [ $(cat $FILE | sed -r 's/^([0-9]).*$/\1/') -eq 0 ]; then
            echo "0" && return  # Cores 0–19 → GPU 0
        else
            echo "1" && return  # Cores 20–39 → GPU 1
        fi
    else
        # Default fallback: use both GPUs if allocation unknown
        echo "0,1" && return
    fi
}

gpus=$(get_gpu-ids)
echo "✓ GPU(s) assigned to this job: $gpus"

# Set CUDA GPU visibility (if GROMACS is CUDA-enabled)
export CUDA_VISIBLE_DEVICES=$gpus

################################################################################
# User Simulation Commands
################################################################################

module load gromacs/gnu/2024.4
module load gcc/12.2.0

# Add your GROMACS commands below (grompp, mdrun, etc.)
# Example with GPU acceleration:
#   gmx grompp -f simulation.mdp -o simulation.tpr -c input.gro -p topol.top
#   gmx mdrun -v -deffnm simulation -nb gpu -bonded gpu -pme gpu -pmefft gpu

[[YOUR GROMACS SCRIPTS]]

################################################################################
# Job Cleanup (Baal GPU-specific)
################################################################################
# Note: Jobs run in original directory, no temporary files or sync needed
# All output remains in the original working directory

if [ $NUM_NODES -le 1 ]; then
    echo "✓ Job completed successfully in original directory"
    echo "GPU output files ready in: $WORKDIR"
fi













