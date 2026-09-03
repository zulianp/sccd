#!/usr/bin/env bash

#SBATCH --job-name=sccd-bench
#SBATCH --account=c40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --time=04:00:00
#SBATCH --output=slurm-sccd-bench-%j.out
#SBATCH --error=slurm-sccd-bench-%j.err
#SBATCH --exclusive
#SBATCH --partition=normal

set -euo pipefail

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true

export PATH=/capstor/scratch/cscs/zulianp/installations/smesh/bin:$PATH

echo "#---------------#"
date
echo "#---------------#"


set -x

 SCCD_BENCH_EXECUTION_SPACE=device srun ./bench.sh 

set +x

echo "#---------------#"
date
echo "#---------------#"
