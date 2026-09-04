#!/bin/bash
# One batch job, one CSV: the measurements every keep-or-demote call in
# wip/ASSESSMENT.md has to cite.
#
# Grace and Hopper are the same node on Alps (GH200), which is convenient rather
# than a compromise: the whole sweep runs inside one allocation, and that is a
# hard requirement here rather than a nicety. Between srun allocations this
# harness varies by about 40%, so a difference measured across two of them means
# nothing unless it is larger than that. Variants are therefore interleaved --
# the outer loop is the repeat, the inner loop the variant -- so any drift over
# the job's lifetime lands on every variant equally instead of on whichever one
# ran last.
#
# Usage:
#   sbatch benchmark/assess.sbatch.sh
#   ASSESS_REPEATS=5 ASSESS_CASES=24 sbatch benchmark/assess.sbatch.sh
#
# Output:
#   $ASSESS_OUT/assessment.csv    one row per (variant, repeat, phase)
#   $ASSESS_OUT/assess.log        stdout of every driver invocation
#
#SBATCH --account=c40
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=72
#SBATCH --time=03:00:00
#SBATCH --job-name=sccd-assess
#SBATCH --output=assess-%j.out
#SBATCH --uenv=prgenv-gnu/24.11:v2
#SBATCH --view=default

set -uo pipefail

if [ -z "${SCRATCH:-}" ]; then
    echo "error: \$SCRATCH is not set. On Alps the environment sets it; elsewhere," >&2
    echo "       set it to a directory holding the build trees and the datasets." >&2
    exit 1
fi
: "${ASSESS_ROOT:=$SCRATCH/sccd-assess}"
: "${ASSESS_OUT:=$ASSESS_ROOT/benchmark/assessment}"
: "${ASSESS_DATA:=$ASSESS_ROOT/data}"
: "${BUILD_GRACE:=$ASSESS_ROOT/build-grace}"
: "${BUILD_HOPPER:=$ASSESS_ROOT/build-hopper}"
: "${UENV_IMAGE:=prgenv-gnu/24.11:v2}"

# Repeats exist to separate a real difference from run-to-run noise. Three is
# the minimum that lets a median mean anything; raise it for a close call.
: "${ASSESS_REPEATS:=3}"
# Cases per scene. A full scene is 372-1996 cases; the sweep multiplies that by
# every variant, so it is subsampled evenly across the trajectory. The scenes
# themselves are real frames -- synthetic geometry misled this project three
# separate times (a motion whose contact time scaled with the step, a box
# distribution with no pairs, and a scattered volume that made the cell list look
# uniformly bad), so nothing here is generated except the quad scaling study,
# which is explicitly labelled as such.
: "${ASSESS_CASES:=16}"
: "${ASSESS_SCENES:=cloth-ball armadillo-rollers cloth-funnel}"
: "${ASSESS_REFINE_LEVELS:=2}"
# No single driver invocation may consume the whole job. A variant that exceeds
# this is recorded as a timeout row rather than silently starving the rest --
# "this configuration does not finish" is itself a result, and it is the one the
# armadillo edge-edge blowup produces.
: "${ASSESS_TIMEOUT:=600}"

mkdir -p "$ASSESS_OUT"
CSV="$ASSESS_OUT/assessment.csv"
LOG="$ASSESS_OUT/assess.log"

echo "hardware,scene,level,component,variant,phase,repeat,ms,pairs,notes" > "$CSV"
: > "$LOG"

log() { echo "=== $* ===" | tee -a "$LOG"; }

# ---------------------------------------------------------------------------
# Driver wrappers
#
# Each one runs a driver and reduces its stdout to the numbers the CSV wants.
# Pair counts are emitted beside every timing on purpose: a broad phase that
# finds a different pair set is not faster, it is wrong, and a table of times
# alone cannot tell those apart.
# ---------------------------------------------------------------------------

# sccd_bench prints
#   dataset,mode,case,type,queries,prep_ms,broad_ms,narrow_ms,query_narrow_ms,fp,fn,broad_fp,broad_fn
# Sum the per-case milliseconds and queries over the subsample; report broad and
# narrow as separate phases because they are separately optimisable and a single
# total hides which one moved.
run_bench() { # hardware scene component variant repeat  [env assignments...]
    local hw="$1" scene="$2" comp="$3" var="$4" rep="$5"; shift 5
    local build="$BUILD_GRACE"
    [ "$hw" = "hopper" ] && build="$BUILD_HOPPER"

    # A subshell rather than env(1): the variable assignments arrive as separate
    # words and export handles them directly, and the driver then inherits the
    # job's uenv view without a nested uenv run.
    local out
    out=$( export "$@" SCCD_BENCH_MAX_CASES="$ASSESS_CASES"
           timeout "$ASSESS_TIMEOUT" "$build/sccd_bench" "$ASSESS_DATA" "$scene" 2>>"$LOG" )
    local rc=$?
    echo "$out" >> "$LOG"
    if [ "$rc" -ne 0 ]; then
        local why="FAILED rc=$rc"
        [ "$rc" -eq 124 ] && why="TIMEOUT after ${ASSESS_TIMEOUT}s"
        echo "$hw,$scene,-,$comp,$var,broad,$rep,,,$why" >> "$CSV"
        echo "$hw,$scene,-,$comp,$var,narrow,$rep,,,$why" >> "$CSV"
        log "$why  $hw/$scene/$comp/$var rep $rep"
        return
    fi

    echo "$out" | awk -F, -v hw="$hw" -v scene="$scene" -v comp="$comp" \
                        -v var="$var" -v rep="$rep" '
        NR > 1 && NF >= 13 {
            q += $5; broad += $7; narrow += $8; fp += $10; fn += $11; n++
        }
        END {
            if (n == 0) { print "MISSING" ; exit }
            printf "%s,%s,-,%s,%s,broad,%s,%.3f,%d,cases=%d\n", hw, scene, comp, var, rep, broad, q, n
            printf "%s,%s,-,%s,%s,narrow,%s,%.3f,%d,fp=%d;fn=%d\n", hw, scene, comp, var, rep, narrow, q, fp, fn
        }' >> "$CSV"
}

# sccd_refine_scaling prints a fixed-width table:
#   level faces vf_pairs ee_pairs prep_ms bp_fv_ms bp_ee_ms broad_ms narrow_ms ns/pair
run_refine() { # hardware topology component variant repeat  [env assignments...]
    local hw="$1" topo="$2" comp="$3" var="$4" rep="$5"; shift 5
    local build="$BUILD_GRACE"
    [ "$hw" = "hopper" ] && build="$BUILD_HOPPER"

    local out
    out=$( export "$@"
           timeout "$ASSESS_TIMEOUT" "$build/sccd_refine_scaling" "$ASSESS_REFINE_LEVELS" 2>>"$LOG" )
    local rc=$?
    echo "$out" >> "$LOG"
    # A failure here is a result, not an accident to be swallowed: an explicit
    # row puts the gap in the data instead of leaving it as an unexplained
    # absence. Any hopper quad rows already in assessment.csv reading
    # "FAILED rc=134" predate the device vertex-quad narrow phase and should be
    # re-run rather than trusted -- quads now dispatch on the device for both
    # phases.
    if [ "$rc" -ne 0 ]; then
        local why="FAILED rc=$rc"
        [ "$rc" -eq 124 ] && why="TIMEOUT after ${ASSESS_TIMEOUT}s"
        echo "$hw,refine-$topo,-,$comp,$var,narrow,$rep,,,$why" >> "$CSV"
        log "$why  $hw/refine-$topo/$comp/$var rep $rep"
        return
    fi

    echo "$out" | awk -v hw="$hw" -v topo="$topo" -v comp="$comp" \
                      -v var="$var" -v rep="$rep" '
        /^[ ]*[0-9]+[ ]/ && NF == 10 {
            lvl = $1; faces = $2; pairs = $3 + $4; broad = $8; narrow = $9
            printf "%s,refine-%s,%s,%s,%s,broad,%s,%.3f,%d,faces=%d\n", hw, topo, lvl, comp, var, rep, broad, pairs, faces
            printf "%s,refine-%s,%s,%s,%s,narrow,%s,%.3f,%d,faces=%d\n", hw, topo, lvl, comp, var, rep, narrow, pairs, faces
        }' >> "$CSV"
}

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-72}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

log "node $(hostname)  threads=$OMP_NUM_THREADS  repeats=$ASSESS_REPEATS  cases=$ASSESS_CASES"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader >> "$LOG" 2>&1

for build in "$BUILD_GRACE" "$BUILD_HOPPER"; do
    [ -x "$build/sccd_bench" ] || echo "warning: $build/sccd_bench missing" | tee -a "$LOG"
done

# ---------------------------------------------------------------------------
# Warmup. Creating the CUDA context costs hundreds of milliseconds and lands on
# whichever device run happens to go first: a smoke run measured 584 ms for a
# broad phase that takes 0.95 ms once warm. Burning it here keeps it out of
# repeat 1 rather than leaving a 600x outlier for the median to hide.
# ---------------------------------------------------------------------------
log "warmup"
( export SCCD_BENCH_MAX_CASES=1 SCCD_BENCH_EXECUTION_SPACE=device
  timeout "$ASSESS_TIMEOUT" "$BUILD_HOPPER/sccd_bench" "$ASSESS_DATA" \
      "${ASSESS_SCENES%% *}" ) >> "$LOG" 2>&1
( export SCCD_BENCH_MAX_CASES=1
  timeout "$ASSESS_TIMEOUT" "$BUILD_GRACE/sccd_bench" "$ASSESS_DATA" \
      "${ASSESS_SCENES%% *}" ) >> "$LOG" 2>&1

# ---------------------------------------------------------------------------
# The sweep. Repeat is the OUTER loop -- see the note at the top on why.
# ---------------------------------------------------------------------------
for rep in $(seq 1 "$ASSESS_REPEATS"); do
    log "repeat $rep of $ASSESS_REPEATS"

    for scene in $ASSESS_SCENES; do

        # Narrow-phase mode. Mode 0 is a duplicate of mode 2's job but is not
        # dominated by it, so both ship and both are measured here.
        #
        # Only 0 and 2. SCCD_NARROWPHASE_MODE=1 warns and runs Relaxed, so
        # sweeping it would measure mode 0 twice under two names.
        for mode in 0 2; do
            run_bench grace "$scene" narrowphase "mode$mode" "$rep" \
                SCCD_NARROWPHASE_MODE="$mode"
        done

        # Broad phase. Both ship today, selected by broadphase_strategy.hpp; this
        # is the evidence for keeping both, or for dropping one.
        for bp in sweep cell2d; do
            run_bench grace "$scene" broadphase "$bp" "$rep" \
                SCCD_BROADPHASE="$bp" SCCD_NARROWPHASE_MODE=2
        done

        # Host versus device, same scene and same cases.
        run_bench hopper "$scene" device "cuda" "$rep" \
            SCCD_BENCH_EXECUTION_SPACE=device SCCD_NARROWPHASE_MODE=2
        run_bench grace "$scene" device "host" "$rep" \
            SCCD_BENCH_EXECUTION_SPACE=host SCCD_NARROWPHASE_MODE=2
    done

    # Refinement scaling, triangles and quads. Quads are a supported topology and
    # need their own path to optimality in the root finder, so they are measured
    # here rather than assumed to track the triangle numbers. These rows are the
    # baseline that quad root-finder work is judged against.
    # Only mode 0 here, for two independent reasons.
    #
    # The narrow-phase mode comparison belongs on the scene rows above, which use
    # real frames. refine_scaling synthesizes its motion by reflecting the mesh
    # through its centre, which is violent enough that every swept box overlaps
    # every other -- its own usage text says so -- and mode 2 on that all-pairs
    # geometry does not finish in any useful time. A smoke run stalled there for
    # minutes on a two-level cube.
    #
    # And for the quad rows it would be measuring the same code twice: the quad
    # path has exactly one root-finder variant and never consults the mode enum,
    # so SCCD_NARROWPHASE_MODE is silently ignored for quads. Recording that as
    # two variants would manufacture a comparison that does not exist.
    #
    # What these rows are for is scaling with element count, and the quad
    # baseline that quad root-finder work gets judged against.
    for topo in tri quad; do
        topo_env=""
        [ "$topo" = "quad" ] && topo_env="quad"
        run_refine grace "$topo" refine "host" "$rep" \
            SCCD_TOPOLOGY="$topo_env" SCCD_NARROWPHASE_MODE=0
        run_refine hopper "$topo" refine "cuda" "$rep" \
            SCCD_TOPOLOGY="$topo_env" SCCD_BENCH_EXECUTION_SPACE=device SCCD_NARROWPHASE_MODE=0
    done
done

log "done: $(wc -l < "$CSV") rows in $CSV"
