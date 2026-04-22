#!/bin/bash
################################################################################
# GROMACS MD Production Run on TGCC Jean-Zay - H100 GPU Partition
################################################################################
#
# JEAN-ZAY Cluster Information (H100 Partition):
# - GPU: NVIDIA H100 (80 GB, NVLink 4.0, 1,457 TFLOPS peak, Transformer Engine)
# - CPU: AMD EPYC 9654 (96 cores/node, 2-socket, 768GB RAM per node)
# - Network: InfiniBand NDR200 (400 Gbps)
# - Storage: High-performance Lustre filesystem
# - Next-gen architecture with exceptional performance
#
# H100 GPU Optimization Tips:
# 1. One H100 GPU per node with 96-core EPYC 9654 CPU (next-gen EPYC)
# 2. Use 1 MPI task × 96 OpenMP threads (all cores available for CPU work)
# 3. H100 has NVLink 4.0 for ultra-fast host-device communication
# 4. GPU acceleration: Non-bonded, PME, and bonded forces on GPU
# 5. H100 offers highest performance: ideal for production runs on large systems
# 6. 80 GB HBM3 memory allows larger simulation systems directly on GPU
#
################################################################################

#SBATCH --job-name=[XXX]                     # Job name (required)
#SBATCH --output=[XXX].out                   # Standard output file
#SBATCH --error=[XXX].err                    # Standard error file
#SBATCH --account=ade@h100                   # H100 project/account (update as needed)
#SBATCH --nodes=1                            # Use 1 node
#SBATCH --ntasks-per-node=1                  # 1 MPI task for 1 GPU
#SBATCH --gres=gpu:1                         # 1 H100 GPU
#SBATCH --cpus-per-task=96                   # All 96 CPU cores for OpenMP threads
#SBATCH -C h100                              # H100 GPU partition constraint
#SBATCH --hint=multithread                   # Optimize for multithreading
#SBATCH --qos=qos_gpu_h100-t3                # Quality of service for H100
#SBATCH --time=20:00:00                      # Max execution time (≤ 20 hours)

################################################################################
# Module Setup
################################################################################

set -x

module purge

# Load H100-specific environment and GROMACS with CUDA support
module load arch/h100
module load gromacs/2023.4-cuda              # GROMACS compiled with CUDA support

################################################################################
# Environment Variables for GPU Optimization
################################################################################

# GPU-specific GROMACS optimizations for H100
export GMX_GPU_PME_PP_COMMS=true             # Enable PME-PP GPU communication (critical for H100 performance)
export GMX_FORCE_UPDATE_DEFAULT_GPU=1        # Enable GPU acceleration
export GMX_GPU_DD_COMMS=true                 # Enable domain decomposition GPU communication
export GMX_MAXBACKUP=-1                      # Disable backup files (saves I/O)
export GMX_MAXCONSTRWARN=-1                  # Allow constraint warnings without stopping

# OpenMP threading optimization for EPYC 9654 CPU
export OMP_PROC_BIND=close                   # Bind threads to cores for better cache locality
export OMP_PLACES=cores                      # Use physical cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK  # Match OpenMP threads to SLURM allocation
export OMP_STACKSIZE=1G                      # Increase stack size for large systems

################################################################################
# GROMACS Commands
################################################################################

# Note: Using 1 MPI task × 96 OpenMP threads optimized for 1 H100 GPU
# GPU handles: non-bonded, PME, and bonded forces
# CPU assists with communication and data preparation
# H100's NVLink 4.0 ensures minimal GPU-CPU communication bottlenecks

# 1. Pre-processing: Generate topology and run parameters
gmx grompp \
  -f [MDP FILE].mdp \
  -o [STEP].tpr \
  -c [PREV STEP].gro \
  -r [TOPOLOGY].gro \
  -p topol.top \
  -n index.ndx \
  -maxwarn 1

# 2. MD Production Run: GPU-accelerated execution on H100
#    - 1 MPI task pinned to 1 H100 GPU
#    - 96 OpenMP threads leverage full EPYC 9654 CPU for CPU-bound work
#    - All computational kernels offloaded to GPU
#    - PME-PP communication enabled for optimal GPU utilization
#    - Suitable for large biomolecular systems (>200k atoms)
gmx mdrun \
  -v \
  -deffnm [STEP] \
  -ntmpi 1 \
  -ntomp 96 \
  -npme 0 \
  -pin on \
  -nb gpu \
  -bonded gpu \
  -pme gpu \
  -pmefft gpu
