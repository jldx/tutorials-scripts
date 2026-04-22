#!/bin/bash
################################################################################
# GROMACS MD Production Run on TGCC Jean-Zay - A100 GPU Partition
################################################################################
#
# JEAN-ZAY Cluster Information (A100 Partition):
# - GPU: NVIDIA A100 (40 GB, NVLink/PCIe, 312 TFLOPS peak, Tensor Cores)
# - CPU: AMD EPYC 7742 (64 cores/node, 2-socket, 512GB RAM per node)
# - Network: InfiniBand HDR100 (200 Gbps)
# - Storage: High-performance Lustre filesystem
# - Typical: 1 GPU per compute node
#
# A100 GPU Optimization Tips:
# 1. One A100 GPU per node with 64-core EPYC CPU
# 2. Use 1 MPI task × 64 OpenMP threads (all cores available for CPU work)
# 3. A100 has NVLink 3.0 for fast host-device communication
# 4. GPU acceleration: Non-bonded, PME, and bonded forces on GPU
# 5. A100 is excellent for mixed-precision workloads and tensor operations
#
################################################################################

#SBATCH --job-name=[XXX]                     # Job name (required)
#SBATCH --output=[XXX].out                   # Standard output file
#SBATCH --error=[XXX].err                    # Standard error file
#SBATCH --account=ade@a100                   # A100 project/account (update as needed)
#SBATCH --nodes=1                            # Use 1 node
#SBATCH --ntasks-per-node=1                  # 1 MPI task for 1 GPU
#SBATCH --gres=gpu:1                         # 1 A100 GPU
#SBATCH --cpus-per-task=64                   # All 64 CPU cores for OpenMP threads
#SBATCH -C a100                              # A100 GPU partition constraint
#SBATCH --hint=multithread                   # Optimize for multithreading
#SBATCH --qos=qos_gpu_a100-t3                # Quality of service for A100
#SBATCH --time=20:00:00                      # Max execution time (≤ 20 hours)

################################################################################
# Module Setup
################################################################################

set -x

module purge

# Load A100-specific environment and GROMACS with CUDA support
module load arch/a100
module load gromacs/2023.4-cuda              # GROMACS compiled with CUDA support

################################################################################
# Environment Variables for GPU Optimization
################################################################################

# GPU-specific GROMACS optimizations for A100
export GMX_GPU_PME_PP_COMMS=true             # Enable PME-PP GPU communication (especially effective on A100)
export GMX_FORCE_UPDATE_DEFAULT_GPU=1        # Enable GPU acceleration
export GMX_GPU_DD_COMMS=true                 # Enable domain decomposition GPU communication
export GMX_MAXBACKUP=-1                      # Disable backup files (saves I/O)

# OpenMP threading optimization for EPYC CPU
export OMP_PROC_BIND=close                   # Bind threads to cores for better cache locality
export OMP_PLACES=cores                      # Use physical cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK  # Match OpenMP threads to SLURM allocation
export OMP_STACKSIZE=1G                      # Increase stack size for large systems

################################################################################
# GROMACS Commands
################################################################################

# Note: Using 1 MPI task × 64 OpenMP threads optimized for 1 A100 GPU
# GPU handles: non-bonded, PME, and bonded forces
# CPU assists with communication and data preparation

# 1. Pre-processing: Generate topology and run parameters
gmx grompp \
  -f [MDP FILE].mdp \
  -o [STEP].tpr \
  -c [PREV STEP].gro \
  -r [TOPOLOGY].gro \
  -p topol.top \
  -n index.ndx \
  -maxwarn 1

# 2. MD Production Run: GPU-accelerated execution on A100
#    - 1 MPI task pinned to 1 A100 GPU
#    - 64 OpenMP threads leverage full EPYC CPU for CPU-bound work
#    - All computational kernels offloaded to GPU
#    - PME-PP communication enabled for better GPU utilization
gmx mdrun \
  -v \
  -deffnm [STEP] \
  -ntmpi 1 \
  -ntomp 64 \
  -npme 0 \
  -pin on \
  -nb gpu \
  -bonded gpu \
  -pme gpu \
  -pmefft gpu


