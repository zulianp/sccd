#!/bin/bash
# The competitor comparison over whole scenes, as resumable Slurm chunks on Alps.
#
# A chunk is (scene, repeat, case range) and runs every library over that range
# once, through run_comparison.sh, in one job on one GH200 module: one Hopper and
# one Grace bound to 72 CPUs. The debug partition caps a job at 30 minutes and
# the account runs one job at a time, so a whole-scene comparison is a sequence
# of short jobs, any of which can be interrupted. A chunk whose CSV exists is
# skipped, so re-invoking resumes; nothing completed is recomputed.
#
# Repeats are separate chunks: run-to-run spread is only visible across
# processes, and it matters here, because an earliest-time-of-impact search that
# prunes against a running bound can answer differently between runs.
#
# Usage, on the Alps login node:
#   sweep_comparison.sh               submit every outstanding chunk, one at a time
#   sweep_comparison.sh --dry-run     list the chunks and what is outstanding
#   sweep_comparison.sh --merge       concatenate finished chunks into compare.csv
#
# Options:
#   --out DIR       chunk tree       (default $SCRATCH/sccd-compare)
#   --scenes "..."  scenes           (default the three with verified ground truth)
#   --repeats N     independent runs (default 3)
#   --chunk N       cases per chunk  (default 400)
#   --time T        per-job limit    (default 00:29:00)
#   --partition P   (default debug)   --account A (default c40)
#
# Submitting blocks until the last job ends, so run it detached (`setsid nohup`):
# a plain background job dies with the ssh session and leaves its chunk jobs
# orphaned. Then
#   compare_table.py $SCRATCH/sccd-compare/compare.csv
set -u

SCRATCH=${SCRATCH:-/ritom/scratch/cscs/zulianp}
OUT_DIR="$SCRATCH/sccd-compare"
SCENES="armadillo-rollers cloth-ball cloth-funnel"
REPEATS=3
CHUNK=400
TIME_LIMIT="00:29:00"
PARTITION="debug"
ACCOUNT="c40"
UENV="prgenv-gnu/24.11:v2"
DRY_RUN=0
MERGE_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out) OUT_DIR="$2"; shift 2 ;;
        --scenes) SCENES="$2"; shift 2 ;;
        --repeats) REPEATS="$2"; shift 2 ;;
        --chunk) CHUNK="$2"; shift 2 ;;
        --time) TIME_LIMIT="$2"; shift 2 ;;
        --partition) PARTITION="$2"; shift 2 ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --merge) MERGE_ONLY=1; shift ;;
        -h|--help) sed -n '2,31p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "error: unknown argument $1" >&2; exit 2 ;;
    esac
done

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
D=${SCCD_DATA_DIR:-$SCRATCH/sccd-data}
B=${SCCD_COMPETITOR_BUILD:-$SCRATCH/sccd/build-comp-f64}

# Runnable cases, counted the way the harnesses count them: a boxes/<key>/ with
# its pair arrays and a matching queries/<key>.csv.
case_count() {
    local n=0 dir key
    for dir in "$D/$1"/boxes/*/; do
        [[ -d "$dir" ]] || continue
        key="$(basename "$dir")"
        [[ -f "$dir/c0.int32" && -f "$dir/c1.int32" && -f "$D/$1/queries/$key.csv" ]] || continue
        case "$key" in *vf|*ee) n=$((n + 1)) ;; esac
    done
    echo "$n"
}

merge() {
    local merged="$OUT_DIR/compare.csv" n=0 f
    mkdir -p "$OUT_DIR"
    first="$(find "$OUT_DIR" -mindepth 3 -maxdepth 3 -path "$OUT_DIR/*/r*/*.csv" | sort | head -1)"
    [[ -n "$first" ]] || { echo "no finished chunks under $OUT_DIR"; return 1; }
    head -1 "$first" > "$merged"
    while IFS= read -r f; do
        tail -n +2 "$f" >> "$merged"
        n=$((n + 1))
    done < <(find "$OUT_DIR" -mindepth 3 -maxdepth 3 -path "$OUT_DIR/*/r*/*.csv" | sort)
    echo "merged $n chunks -> $merged ($(( $(wc -l < "$merged") - 1 )) rows)"
}

if [[ $MERGE_ONLY -eq 1 ]]; then
    merge
    exit $?
fi

keys=()
for scene in $SCENES; do
    total="$(case_count "$scene")"
    if [[ "$total" -eq 0 ]]; then
        echo "note: $scene has no runnable case under $D; skipping" >&2
        continue
    fi
    for ((r = 1; r <= REPEATS; ++r)); do
        for ((b = 0; b < total; b += CHUNK)); do
            e=$((b + CHUNK)); [[ $e -gt $total ]] && e=$total
            keys+=("$scene|$r|$b|$e")
        done
    done
done

chunk_path() {  # scene r begin end
    printf '%s/%s/r%s/%06d-%06d.csv' "$OUT_DIR" "$1" "$2" "$3" "$4"
}

todo=()
for k in "${keys[@]}"; do
    IFS='|' read -r scene r b e <<< "$k"
    [[ -s "$(chunk_path "$scene" "$r" "$b" "$e")" ]] || todo+=("$k")
done
echo "${#keys[@]} chunks: $(( ${#keys[@]} - ${#todo[@]} )) done, ${#todo[@]} outstanding"

if [[ $DRY_RUN -eq 1 ]]; then
    for k in "${todo[@]:-}"; do
        [[ -n "$k" ]] || continue
        IFS='|' read -r scene r b e <<< "$k"
        printf '  TODO %-20s r%-2s cases [%s, %s)\n' "$scene" "$r" "$b" "$e"
    done
    exit 0
fi

for bin in "$B/sccd_bench" "$B/benchmark/competitors/scalable_ccd_bench" "$B/benchmark/competitors/accd_bench"; do
    [[ -x "$bin" ]] || { echo "error: $bin missing -- build with -DSCCD_ENABLE_COMPETITORS=ON" >&2; exit 1; }
done

# The QOS caps how many jobs one user may have submitted at once, and this
# account is shared with the user's own work. A rejection on that ground means
# "the queue is full right now", not "this chunk failed", so wait for a slot
# rather than burning the chunk -- `benchmark/scripts/sweep.sh` learned the same
# lesson. Anything else is a real error.
submit_with_retry() {
    local attempt=0 err
    while :; do
        err="$(sbatch --wait "$@" 2>&1 >/dev/null)" && return 0
        if [[ "$err" != *QOSMaxSubmitJobPerUserLimit* ]]; then
            printf '%s\n' "$err" >&2
            return 1
        fi
        attempt=$((attempt + 1))
        if [[ "$attempt" -gt "${SWEEP_SUBMIT_RETRIES:-240}" ]]; then
            echo "giving up: the submit limit stayed full for $attempt attempts" >&2
            return 1
        fi
        [[ "$attempt" -eq 1 ]] && echo "    queue full (submit limit); waiting for a slot" >&2
        sleep "${SWEEP_SUBMIT_WAIT:-60}"
    done
}

mkdir -p "$OUT_DIR/logs"
failures=0
for k in "${todo[@]:-}"; do
    [[ -n "$k" ]] || continue
    IFS='|' read -r scene r b e <<< "$k"
    out="$(chunk_path "$scene" "$r" "$b" "$e")"
    name="$scene-r$r-$b"
    echo "==> $scene r$r [$b, $e)  $(date +%T)"
    # --wait serialises the chunks, which is what one running job per user
    # amounts to anyway, and lets a failure be reported against its chunk.
    if ! submit_with_retry --account="$ACCOUNT" --partition="$PARTITION" \
            --nodes=1 --ntasks=1 --cpus-per-task=72 --gpus-per-task=1 \
            --time="$TIME_LIMIT" --uenv="$UENV" --view=default \
            --job-name="cmp-$name" \
            --output="$OUT_DIR/logs/$name.out" --error="$OUT_DIR/logs/$name.err" \
            --wrap="bash '$HERE/run_comparison.sh' '$scene' '$b' '$e' '$out'"; then
        echo "FAILED $name (see $OUT_DIR/logs/$name.err)" >&2
        failures=$((failures + 1))
    fi
done

echo "submitted ${#todo[@]} chunks, $failures failed"
[[ $failures -eq 0 ]] && merge
exit $failures
