#!/bin/bash
# SCCD against Additive CCD on whatever machine this is, host only.
#
# The companion to run_comparison.sh, which needs a GPU and a Slurm allocation.
# Scalable CCD's *broad* phase is host code and runs here; only its narrow phase
# is CUDA, so its rows carry `broad_ms` and leave `narrow_ms` empty. The narrow
# phase comparison against it needs the GPU driver, run_comparison.sh.
#
# Both binaries read the same case-range variables and print the same schema, and
# repeats are separate processes so one run's allocator and cache state cannot
# colour the next one's timings.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${SCCD_BUILD_DIR:-${ROOT}/build_competitors}"
DATA="${SCCD_DATA_DIR:-${ROOT}/data}"
OUT="${SCCD_OUT:-${ROOT}/benchmark/competitors/results/host-$(date +%Y-%m-%d)}"
mkdir -p "$OUT"

SCENES="${SCENES:-armadillo-rollers cloth-ball}"
BEGIN=${BEGIN:-0}
END=${END:-60}
REPEATS=${REPEATS:-3}

SCCD_BIN="$BUILD/sccd_bench"
ACCD_BIN="$BUILD/benchmark/competitors/accd_bench"
SCALABLE_BIN="$BUILD/benchmark/competitors/scalable_ccd_bench"
for b in "$SCCD_BIN" "$ACCD_BIN" "$SCALABLE_BIN"; do
    [ -x "$b" ] || { echo "missing $b -- build with -DSCCD_ENABLE_COMPETITORS=ON"; exit 1; }
done

"$SCCD_BIN" --header > "$OUT/all.csv"

run() {
    local label="$1"; shift
    local t0 t1
    t0=$(date +%s)
    "$@" 2>"$OUT/$label.err" | tail -n +2 >> "$OUT/all.csv"
    local rc=${PIPESTATUS[0]}
    t1=$(date +%s)
    printf '  %-26s exit=%-3d %4ds\n' "$label" "$rc" "$((t1 - t0))"
    [ "$rc" -ne 0 ] && head -3 "$OUT/$label.err"
    return 0
}

echo "host comparison: $(uname -sm), cases [$BEGIN,$END) x$REPEATS"
for scene in $SCENES; do
    echo "=== $scene ==="
    for r in $(seq 1 "$REPEATS"); do
        for mode in 0 2; do
            run "sccd-m$mode-r$r" env \
                SCCD_BENCH_EXECUTION_SPACE=host SCCD_NARROWPHASE_MODE=$mode \
                SCCD_BROADPHASE=sweep SCCD_BENCH_CASE_BEGIN=$BEGIN SCCD_BENCH_CASE_END=$END \
                "$SCCD_BIN" "$DATA" "$scene"
        done
        run "scalable-host-r$r" env \
            SCCD_BENCH_CASE_BEGIN=$BEGIN SCCD_BENCH_CASE_END=$END \
            "$SCALABLE_BIN" "$DATA" "$scene"

        # The same sweep broad phase the SCCD rows use, so both per-pair narrow
        # phases are handed the identical candidate list.
        run "accd-r$r" env \
            SCCD_BROADPHASE=sweep SCCD_BENCH_CASE_BEGIN=$BEGIN SCCD_BENCH_CASE_END=$END \
            "$ACCD_BIN" "$DATA" "$scene"
    done
done

echo
echo "rows: $(( $(wc -l < "$OUT/all.csv") - 1 ))"
echo "csv:  $OUT/all.csv"
