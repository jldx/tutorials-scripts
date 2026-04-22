#!/bin/bash
################################################################################
# GROMACS MD Production Run on TGCC Jean-Zay - CPU Partition (AMD EPYC)
################################################################################
#
# JEAN-ZAY Cluster Information (CPU Partition):
# - Compute Nodes: AMD EPYC 7742 (64 cores/node, 2-socket, 512GB RAM per node)
# - Network: InfiniBand HDR100 (200 Gbps)
# - Storage: High-performance Lustre filesystem
# - Architecture: Modern NUMA-aware CPU design
#
# CPU Partition Optimization Tips:
# 1. Each node has 64 cores across 2 sockets (32 cores per socket)
# 2. NUMA topology: Use -pin to bind threads to cores for optimal cache locality
# 3. MPI + OpenMP hybrid parallelism is recommended (e.g., 2 MPI × 32 OpenMP)
# 4. EPYC 7742 has 3.4 GHz boost and excellent memory bandwidth
#
################################################################################

#SBATCH --job-name=[XXX]                     # Job name (required)
#SBATCH --output=[XXX].out                   # Standard output file
#SBATCH --error=[XXX].err                    # Standard error file
#SBATCH --account=ade@cpu                    # CPU project/account (update as needed)
#SBATCH --nodes=1                            # Use 1 node (AMD EPYC, 64 cores)
#SBATCH --ntasks-per-node=1                  # 1 MPI task per node
#SBATCH --cpus-per-task=64                   # All 64 CPU cores for OpenMP
#SBATCH --hint=multithread                   # Optimize for multithreading
#SBATCH --qos=qos_cpu-t3                     # QoS for CPU jobs (max 20h; use qos_cpu-t4 for ≤100h)
#SBATCH --time=20:00:00                      # Max execution time (adjust as needed)

################################################################################
# Module Setup
################################################################################

# Clean environment
module purge

# Load required modules for AMD EPYC + GROMACS
module load arch/amd
module load gromacs/2023.4

################################################################################
# Environment Variables for CPU Optimization
################################################################################

# OpenMP threading optimization for NUMA architecture
export OMP_PROC_BIND=close                   # Bind threads to cores for better cache locality
export OMP_PLACES=cores                      # Use physical cores (not hardware threads)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK  # Match OpenMP threads to SLURM allocation
export OMP_STACKSIZE=1G                      # Increase stack size for large systems

# GROMACS-specific optimizations
export GMX_MAXBACKUP=-1                      # Disable backup files (saves I/O)
export GMX_FORCE_UPDATE_DEFAULT_GPU=0        # Disable GPU (CPU-only run)

set -x

################################################################################
# GROMACS Commands
################################################################################

# Note: Using 2 MPI tasks × 32 OpenMP threads for optimal EPYC partition usage
# (each MPI task gets one socket of the dual-socket node)

# 1. Pre-processing: Generate topology and run parameters (single-threaded)
gmx grompp \
  -f [MDP FILE].mdp \
  -o [STEP].tpr \
  -c [PREV STEP].gro \
  -r [TOPOLOGY].gro \
  -p topol.top \
  -n index.ndx \
  -maxwarn 1

# 2. MD Production Run: Hybrid MPI + OpenMP for CPU partition
#    - Uses 2 MPI tasks (one per socket) × 32 OpenMP threads each
#    - Pin cores to threads for optimal NUMA performance
#    - Force CPU computation (no GPU)
srun gmx_mpi mdrun \
  -v \
  -deffnm [STEP] \
  -ntmpi 2 \
  -ntomp 32 \
  -pin on \
  -pinoffset 0 \
  -nb cpu \
  -pme cpu \
  -bonded cpu \
  -dlb auto
