#! /bin/bash
#SBATCH --job-name=hbonds_0_250                  # Job name (required)
#SBATCH --output=hbonds_0_250.out                   # Standard output file
#SBATCH --error=hbonds_0_250.err                    # Standard error file
#SBATCH --account=ade@cpu                    # CPU project/account (update as needed)
#SBATCH --nodes=1                            # Use 1 node (AMD EPYC, 64 cores)
#SBATCH --ntasks-per-node=1                  # 1 MPI task per node
#SBATCH --cpus-per-task=64                   # All 64 CPU cores for OpenMP
#SBATCH --hint=multithread                   # Optimize for multithreading
#SBATCH --qos=qos_cpu-dev                     # QoS for CPU jobs (max 20h; use qos_cpu-t4 for ≤100h)
#SBATCH --time=00:10:00                      # Max execution time (adjust as needed)


module purge
module load vmd/1.9.4a43-text
module load gcc/9.1.0 
module load netcdf-c/4.7.2

export VMDFORCECPUCOUNT=64
export WKFFORCECPUCOUNT=64
export RTFORCECPUCOUNT=64


TRAJ_STEP=(
0_25
25_50
50_75
75_100
100_125
125_150
150_175
175_200
200_225
225_250
)


for chunk in "${TRAJ_STEP[@]}"; do 

    vmd -dispdev none -e run_hbonds_${chunk}.tcl

done


