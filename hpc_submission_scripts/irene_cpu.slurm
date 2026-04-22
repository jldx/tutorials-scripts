#!/bin/bash
################################################################################
# GROMACS MD Production Run on TGCC Irene - Rome Partition (AMD EPYC 7H12)
################################################################################
#
# IRENE Cluster Information:
# - Rome partition: AMD EPYC 7H12 (128 cores/node, 2-socket, 256GB RAM per node)
# - Network: InfiniBand HDR100 (200 Gbps)
# - Interconnect: Torus 3D topology for excellent scaling
# - Storage: Scratch & Working directories on high-performance storage
#
# Rome Partition Optimization Tips:
# 1. Each node has 128 cores across 2 sockets (64 cores per socket)
# 2. NUMA topology: Use -pin to bind threads to cores for optimal cache locality
# 3. MPI + OpenMP hybrid parallelism is recommended (e.g., 2 MPI × 64 OpenMP)
# 4. Rome has 3MB L3 cache per core and fast memory bandwidth
#
################################################################################

#MSUB -r [XXX]              # Job name (required)
#MSUB -o [XXX].o                 # Standard output file
#MSUB -e [XXX].e                 # Standard error file
#MSUB -T 3600                    # Time limit in seconds (1h = 3600s; adjust as needed)
#MSUB -m scratch                 # Use scratch storage (temporary working space)
#MSUB -x                         # Exclusive node access (no job sharing)
#MSUB -N 1                       # Number of nodes (1 node = 128 cores)
#MSUB -n 128                     # Total number of MPI tasks (128 = all cores on 1 node)
#MSUB -A gen6818                 # Project account code (update as needed)
#MSUB -q rome                    # Queue/partition: Rome (AMD EPYC)
#MSUB -@ [MAIL ADRESS]:begin,end # Email notifications at job begin & end

################################################################################
# Module Setup
################################################################################

module purge
module load gromacs

################################################################################
# Environment Variables for Rome Partition Optimization
################################################################################

# OpenMP threading optimization for NUMA architecture
export OMP_PROC_BIND=close              # Bind threads to cores (better cache locality)
export OMP_PLACES=cores                 # Use physical cores (not hardware threads)
export OMP_NUM_THREADS=64               # 64 threads per MPI task (half node per MPI task)
export OMP_STACKSIZE=1G                 # Increase stack size for large molecular systems

# GROMACS-specific optimizations
export GMX_MAXBACKUP=-1                 # Disable backup files (saves I/O)
export GMX_FORCE_UPDATE_DEFAULT_GPU=0   # Ensure CPU-only computation (no GPU)

# MPI optimization for Rome's InfiniBand HDR100 network
export OMPI_MCA_btl=^openib             # Use efficient fabric (auto-detected)

################################################################################
# GROMACS Commands
################################################################################

# Note: Using 2 MPI tasks × 64 OpenMP threads for optimal Rome partition usage
# (each MPI task gets one socket of the dual-socket node)

# 1. Pre-processing: Generate topology and run parameters (MPI not needed)
ccc_mprun -n 1 gmx_mpi grompp \
  -f [MDP FILE].mdp \
  -o [STEP].tpr \
  -c [PREV STEP].gro \
  -r [TOPOLOGY].gro \
  -p topol.top \
  -n index.ndx \
  -maxwarn 1

# 2. MD Production Run: Hybrid MPI + OpenMP for Rome partition
#    - Uses 2 MPI tasks (one per socket) × 64 OpenMP threads each
#    - Pin cores to threads for optimal NUMA performance
#    - Force CPU computation (no GPU)
ccc_mprun -n 2 gmx_mpi mdrun \
  -v \
  -deffnm [STEP] \
  -ntomp 64 \
  -pin on \
  -pinoffset 0 \
  -nb cpu \
  -pme cpu \
  -bonded cpu \
  -dlb auto