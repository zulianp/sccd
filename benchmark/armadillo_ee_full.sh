#!/bin/bash
# Every armadillo-rollers case, mode 0 against mode 2, per case.
#
# This settles "mode 2 is about 100x slower on armadillo edge-edge", which the
# variant sweep in assess.sbatch.sh could not: that sweep samples 16 cases evenly
# across the trajectory, which is the right shape for "is A faster than B" and the
# wrong shape for "is there one case where A explodes". One case in four hundred
# is invisible to a sum. So this runs the whole scene and emits per-case rows.
#
# Output feeds benchmark/assessment/armadillo-ee-full.csv, cited by the
# "Withdrawn" section of wip/ASSESSMENT.md.
#
# Usage, inside an allocation with a GH200 node (about 12 minutes for two modes):
#   srun --account=... --partition=debug --nodes=1 --ntasks=1 --gpus-per-task=1 \
#        --cpus-per-task=72 --time=00:29:00 --uenv=prgenv-gnu/24.11:v2 \
#        --view=default benchmark/armadillo_ee_full.sh
set -u
: "${SCRATCH:=/capstor/scratch/cscs/zulianp}"
: "${EE_ROOT:=$SCRATCH/sccd-npgap}"
: "${BENCH:=$EE_ROOT/build-hopper/sccd_bench}"
: "${DATA:=$EE_ROOT/data}"
: "${OUT:=$EE_ROOT/ee}"
: "${REPEAT:=1}"

mkdir -p "$OUT"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-72}" OMP_PROC_BIND=close OMP_PLACES=cores
export SCCD_BENCH_EXECUTION_SPACE=host
# Unset, every case runs. That is the entire point of this script.
unset SCCD_BENCH_MAX_CASES

echo "node $(hostname), repeat $REPEAT, all cases, host, $(date)"
for mode in 0 2; do
    echo "=== mode $mode ==="
    # Line-buffered: a case that really did blow up would run into the wall clock,
    # and block-buffered output would lose every row up to that point.
    SCCD_NARROWPHASE_MODE=$mode stdbuf -oL timeout 700 \
        "$BENCH" "$DATA" armadillo-rollers > "$OUT/r${REPEAT}_mode$mode.csv" 2>"$OUT/r${REPEAT}_mode$mode.err"
    echo "  rc=$? rows=$(wc -l < "$OUT/r${REPEAT}_mode$mode.csv")"
done
