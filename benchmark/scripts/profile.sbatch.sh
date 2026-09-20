#!/bin/bash
# Hardware counters for the two narrow phases, Grace and Hopper, in one job.
#
# The paper reports arithmetic intensity, occupancy and warp utilisation for the
# device kernel (all of them derivable from the kernel's own instrumentation) and
# nothing at all about memory throughput or about the host. This collects both.
# Everything it needs is a counter read; it changes no result and gates nothing.
#
#   sbatch benchmark/scripts/profile.sbatch.sh
#
# Writes benchmark/results/profile/ as one text file per measurement, which is
# the directory wip/paper/scripts/make_inputs.py would read.
#
#SBATCH --account=c40
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --uenv=prgenv-gnu/24.11:v2
#SBATCH --view=default
#SBATCH --job-name=sccd-profile

set -euo pipefail
cd "$(dirname "$0")/../.."
out=benchmark/results/profile
mkdir -p "$out"

BUILD=${BUILD:-build-hopper}
SCENE=${SCENE:-cloth-funnel}
DATA=${DATA:-data/$SCENE}

# ---------------------------------------------------------------- Hopper ----
# Three metric groups, kept small on purpose: each replays the kernel, and the
# narrow phase is long enough that asking for everything at once is an hour.
#
#   dram__bytes             DRAM traffic, both directions
#   lts__t_sectors          L2 sectors, to separate L2 hits from DRAM
#   l1tex__t_sectors        L1/texture sectors
#   gpu__time_duration      so a throughput can be formed
#
# The achieved-occupancy and warp-execution-efficiency metrics are what the
# paper's 12.5% and 20.4% should be checked against: the first is measured here
# rather than derived from the register count, the second is the hardware's own
# view of the divergence the kernel instrumentation counts.
NCU=${NCU:-ncu}

echo "== hopper: memory throughput =="
srun "$NCU" --target-processes all \
  --kernel-name-base demangled \
  --kernel-name 'regex:narrow_phase' \
  --launch-count 20 \
  --metrics \
dram__bytes_read.sum,dram__bytes_write.sum,\
dram__throughput.avg.pct_of_peak_sustained_elapsed,\
lts__t_sectors.sum,l1tex__t_sectors.sum,\
gpu__time_duration.sum \
  "$BUILD/sccd_bench" --dataset "$DATA" --space device --mode 2 \
  > "$out/hopper-memory.txt" 2>&1 || echo "ncu memory pass failed, see the file"

echo "== hopper: occupancy and divergence =="
srun "$NCU" --target-processes all \
  --kernel-name-base demangled \
  --kernel-name 'regex:narrow_phase' \
  --launch-count 20 \
  --metrics \
sm__warps_active.avg.pct_of_peak_sustained_active,\
smsp__thread_inst_executed_per_inst_executed.ratio,\
smsp__inst_executed.avg.per_cycle_active,\
launch__occupancy_limit_registers,launch__registers_per_thread \
  "$BUILD/sccd_bench" --dataset "$DATA" --space device --mode 2 \
  > "$out/hopper-occupancy.txt" 2>&1 || echo "ncu occupancy pass failed, see the file"

echo "== hopper: stall reasons =="
srun "$NCU" --target-processes all \
  --kernel-name-base demangled \
  --kernel-name 'regex:narrow_phase' \
  --launch-count 20 \
  --metrics \
smsp__average_warps_issue_stalled_barrier_per_issue_active.ratio,\
smsp__average_warps_issue_stalled_long_scoreboard_per_issue_active.ratio,\
smsp__average_warps_issue_stalled_wait_per_issue_active.ratio \
  "$BUILD/sccd_bench" --dataset "$DATA" --space device --mode 2 \
  > "$out/hopper-stalls.txt" 2>&1 || echo "ncu stall pass failed, see the file"

# ----------------------------------------------------------------- Grace ----
# perf on the host narrow phase. The lane-packed kernel should be issue-bound
# rather than memory-bound if the per-lane geometry replication is doing its job,
# so the numbers to look at are IPC and the LLC miss rate: a high IPC with few
# LLC misses is the intended behaviour, and a low IPC with many is the gather the
# replication exists to avoid.
#
# Counter names vary by kernel and by Neoverse revision; the fallback list is the
# architectural set, which is present everywhere.
echo "== grace: host narrow phase =="
EVENTS=${EVENTS:-cycles,instructions,cache-references,cache-misses,branch-misses,stalled-cycles-frontend,stalled-cycles-backend}
srun perf stat -e "$EVENTS" -- \
  "$BUILD/sccd_bench" --dataset "$DATA" --space host --mode 2 \
  > "$out/grace-narrow.txt" 2>&1 || echo "perf failed, see the file"

echo "== grace: memory bandwidth, host broad phase =="
srun perf stat -e "$EVENTS" -- \
  "$BUILD/sccd_bench" --dataset "$DATA" --space host --mode 2 --phase broad \
  > "$out/grace-broad.txt" 2>&1 || echo "perf failed, see the file"

echo "wrote $out"
ls -la "$out"
