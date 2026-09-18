#!/bin/bash
# One comparison chunk: every library over one case range of one scene, once.
#
#   run_comparison.sh <scene> <begin> <end> <out.csv>
#
# `sweep_comparison.sh` calls this inside a Slurm job and owns the chunking; it
# can also be run by hand inside an allocation. Every binary reads the same
# case-range variables and prints the same schema, so the chunk is one CSV. Each
# library runs in its own process, so one allocator's and cache's state cannot
# colour the next one's timings.
#
# What each run contributes:
#
#   sccd_bench           both processors, Relaxed (0) and Tight (2), sweep broad
#                        phase. `narrow_ms` is the earliest-time-of-impact path,
#                        `narrow_ms_s1` the per-pair path, and the accuracy
#                        columns score the curated queries.
#   scalable_ccd_bench   device: its earliest-time-of-impact pipeline per case,
#                        as shipped. host: its broad phase alone.
#   accd_bench           per pair over SCCD's host sweep candidates, and scored on
#                        the curated queries.
#
# The CSV is written through a temporary file and renamed, so an existing
# <out.csv> is a chunk that finished. Any library exiting non-zero leaves the
# chunk unfinished: a partial chunk merges silently and thins one library.
set -u
set -o pipefail

if [[ $# -ne 4 ]]; then
    echo "usage: $0 <scene> <begin> <end> <out.csv>" >&2
    exit 2
fi
scene="$1" begin="$2" end="$3" out="$4"

SCRATCH=${SCRATCH:-/ritom/scratch/cscs/zulianp}
# The smesh whose geom_t the build and the prepared data agree on. These
# benchmarks want the double one: the dataset's exact roots belong to the
# rational query geometry, and a float32 mesh answers about a different geometry,
# which shows up as times of impact that look late and are not.
SMESH_PREFIX=${SCCD_SMESH_PREFIX:-$SCRATCH/installations/smesh-f64}
export PATH=$SMESH_PREFIX/bin:$PATH
D=${SCCD_DATA_DIR:-$SCRATCH/sccd-data}
B=${SCCD_COMPETITOR_BUILD:-$SCRATCH/sccd/build-comp-f64}

# One Grace, 72 threads. OMP_NUM_THREADS covers SCCD; Scalable CCD's host broad
# phase runs on oneTBB, which follows the CPU affinity mask instead, so the job
# has to be bound to one Grace (--cpus-per-task=72). The line records what the
# chunk actually had.
export OMP_NUM_THREADS=72
echo "$(hostname) cpus=$(nproc) OMP_NUM_THREADS=$OMP_NUM_THREADS GPU=${CUDA_VISIBLE_DEVICES:-unset} $scene [$begin,$end)"

tmp="${out}.partial"
mkdir -p "$(dirname "$out")"
"$B/sccd_bench" --header > "$tmp" || { echo "error: sccd_bench printed no header" >&2; exit 1; }

range=(SCCD_BENCH_CASE_BEGIN="$begin" SCCD_BENCH_CASE_END="$end")

run() {  # $1 label  $2.. command
    local label="$1"; shift
    local t0 rc
    t0=$(date +%s)
    "$@" 2>"${tmp%.partial}.$label.err" | tail -n +2 >> "$tmp"
    rc=$?
    printf '  %-26s exit=%-3d %5ds\n' "$label" "$rc" "$(( $(date +%s) - t0 ))"
    if [[ $rc -ne 0 ]]; then
        head -5 "${tmp%.partial}.$label.err" >&2
        echo "error: $label failed; leaving $out unfinished" >&2
        rm -f "$tmp"
        exit 1
    fi
}

for space in device host; do
    for mode in 0 2; do
        run "sccd-$space-m$mode" env "${range[@]}" \
            SCCD_BENCH_EXECUTION_SPACE=$space SCCD_NARROWPHASE_MODE=$mode SCCD_BROADPHASE=sweep \
            "$B/sccd_bench" "$D" "$scene"
    done
done

for space in device host; do
    run "scalable-$space" env "${range[@]}" SCCD_BENCH_EXECUTION_SPACE=$space \
        "$B/benchmark/competitors/scalable_ccd_bench" "$D" "$scene"
done

# The same sweep broad phase the SCCD rows use, so both per-pair narrow phases
# are handed the same candidates.
run "accd" env "${range[@]}" SCCD_BROADPHASE=sweep \
    "$B/benchmark/competitors/accd_bench" "$D" "$scene"

if [[ "$(wc -l < "$tmp")" -le 1 ]]; then
    echo "error: $out produced no rows; leaving it unfinished" >&2
    rm -f "$tmp"
    exit 1
fi
mv "$tmp" "$out"
echo "done $out ($(( $(wc -l < "$out") - 1 )) rows)"
