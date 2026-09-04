#!/usr/bin/env bash
#
# Run the full benchmark sweep as a batch job on Alps.
#
#     sbatch benchmark/scripts/launch.sbatch.sh
#
# bench.sh does the work; this only supplies the allocation and the environment.
# Point SCCD_SMESH_BIN at an smesh installation's bin directory if the harness
# needs its tools on PATH; nothing here assumes one particular account's scratch.

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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-72}"
export OMP_PROC_BIND=true

if [ -n "${SCCD_SMESH_BIN:-}" ]; then
    export PATH="${SCCD_SMESH_BIN}:${PATH}"
fi

echo "#---------------#"
date
echo "#---------------#"

set -x
SCCD_BENCH_EXECUTION_SPACE=device srun "${SCRIPT_DIR}/bench.sh"
set +x

echo "#---------------#"
date
echo "#---------------#"
