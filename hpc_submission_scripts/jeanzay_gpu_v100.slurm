#!/bin/bash
################################################################################
# GROMACS MD Production Run on TGCC Jean-Zay - V100 GPU Partition
################################################################################
#
# JEAN-ZAY Cluster Information (V100 Partition):
# - GPU: NVIDIA V100 (32 GB or 16 GB options, PCIe, 125 TFLOPS peak)
# - CPU: Intel Xeon Gold 6148 (24 cores/node, 2-socket, 384GB RAM per node)
# - Network: InfiniBand EDR (100 Gbps)
# - Storage: High-performance Lustre filesystem
# - Typical: 1 GPU per compute node
#
# V100 GPU Optimization Tips:
# 1. One V100 GPU per node (PCIe bandwidth: ~32 GB/s)
# 2. Use 1 MPI task × 24 OpenMP threads (matches 1 GPU per node setup)
# 3. GPU acceleration: Non-bonded, PME, and bonded forces on GPU
# 4. NVIDIA V100: Strong for mixed-precision and tensor operations
# 5. Check available memory: -C v100-32g for 32GB, -C v100-16g for 16GB
#
################################################################################

#SBATCH --job-name=[XXX]                     # Job name (required)
#SBATCH --output=[XXX].out                   # Standard output file
#SBATCH --error=[XXX].err                    # Standard error file
#SBATCH --account=ade@v100                   # V100 project/account (update as needed)
#SBATCH --nodes=1                            # Use 1 node
#SBATCH --ntasks-per-node=1                  # 1 MPI task for 1 GPU
#SBATCH --gres=gpu:1                         # 1 V100 GPU
#SBATCH --cpus-per-task=24                   # All 24 CPU cores for OpenMP threads
#SBATCH -C v100-32g                          # Request 32 GB V100 GPUs (use v100-16g for 16 GB)
#SBATCH --hint=multithread                   # Optimize for multithreading
#SBATCH --qos=qos_gpu-t3                     # Quality of service for V100
#SBATCH --time=20:00:00                      # Max execution time (≤ 100 hours)

################################################################################
# Module Setup
################################################################################

set -x

module purge

# Load V100-specific environment and GROMACS with CUDA support
module load arch/v100
module load gromacs/2023.4-cuda             # GROMACS compiled with CUDA support

################################################################################
# Environment Variables for GPU Optimization
################################################################################

# GPU-specific GROMACS optimizations
export GMX_GPU_PME_PP_COMMS=true             # Enable PME-PP GPU communication
export GMX_FORCE_UPDATE_DEFAULT_GPU=1        # Enable GPU acceleration
export GMX_GPU_DD_COMMS=true                 # Enable domain decomposition GPU communication
export GMX_MAXBACKUP=-1                      # Disable backup files (saves I/O)

# OpenMP threading optimization (for CPU-bound work)
export OMP_PROC_BIND=close                   # Bind threads to cores for better cache locality
export OMP_PLACES=cores                      # Use physical cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK  # Match OpenMP threads to SLURM allocation

################################################################################
# GROMACS Commands
################################################################################

# Note: Using 1 MPI task × 24 OpenMP threads optimized for 1 V100 GPU
# GPU handles: non-bonded, PME, and bonded forces

# 1. Pre-processing: Generate topology and run parameters
gmx grompp \
  -f [MDP FILE].mdp \
  -o [STEP].tpr \
  -c [PREV STEP].gro \
  -r [TOPOLOGY].gro \
  -p topol.top \
  -n index.ndx \
  -maxwarn 1

# 2. MD Production Run: GPU-accelerated execution
#    - 1 MPI task pinned to 1 GPU
#    - 24 OpenMP threads for CPU work (non-GPU tasks)
#    - All computational kernels offloaded to GPU
gmx mdrun \
  -v \
  -deffnm [STEP] \
  -ntmpi 1 \
  -ntomp 24 \
  -npme 0 \
  -pin on \
  -nb gpu \
  -bonded gpu \
  -pme gpu \
  -pmefft gpu


