#!/usr/bin/env bash
set -euo pipefail

# The goal of this file is to perform a wide benchmark of the SCCD library using the NYU datasets
# We benchmarki the broad-phase and narrow-phase collision detection:
# 1) for performance in timing
# 2) for accuracy compared to the ground truth
# 3) Number of false positives (should be as low as possible for the narrow-phase) and negatives (should be 0)
# For cuda on Alps: OMP_NUM_THREADS=72 SCCD_BENCH_EXECUTION_SPACE=device srun ./bench.sh 

# Narrow-phase modes to sweep (see src/narrowphase/sccd_narrowphase_mode.hpp):
#   0 relaxed   2 tight
# Those are the two that exist. The binary is run once per mode and the rows are
# concatenated; every row carries its mode, so the postprocessor keeps the series
# apart. 1 and 3 were retired: setting either warns and runs Relaxed, so a sweep
# that included them measured mode 0 twice under two different names.
export SCCD_BENCH_MODES="${SCCD_BENCH_MODES:-0 2}"

# Datasets to run. Overridable from the environment so a quick check can target
# a single dataset without editing this file.
export SCCD_ENABLE_ARMADILLO_ROLLERS="${SCCD_ENABLE_ARMADILLO_ROLLERS:-1}"
export SCCD_ENABLE_CLOTH_BALL="${SCCD_ENABLE_CLOTH_BALL:-1}"
export SCCD_ENABLE_CLOTH_FUNNEL="${SCCD_ENABLE_CLOTH_FUNNEL:-1}"
export SCCD_ENABLE_N_BODY_SIMULATION="${SCCD_ENABLE_N_BODY_SIMULATION:-0}"
export SCCD_ENABLE_PUFFER_BALL="${SCCD_ENABLE_PUFFER_BALL:-0}"
export SCCD_ENABLE_ROD_TWIST="${SCCD_ENABLE_ROD_TWIST:-0}"

# Skip the dataset download when the data is already in place.
export SCCD_SKIP_DOWNLOAD="${SCCD_SKIP_DOWNLOAD:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

exec 3>&1
exec 1>&2

DATA_DIR="${SCCD_DATA_DIR:-"${BENCHMARK_DIR}/../data"}"
SCCD_BUILD_DIR="${SCCD_BUILD_DIR:-"${BENCHMARK_DIR}/../build_benchmark"}"
PYTHON="${PYTHON:-python3}"
ROOT_DIR="${BENCHMARK_DIR}/.."

is_enabled() {
    case "${1:-0}" in
        1|ON|on|On|TRUE|true|True|YES|yes|Yes) return 0 ;;
        *) return 1 ;;
    esac
}

parallel_jobs() {
    if command -v nproc >/dev/null 2>&1; then
        nproc
    elif command -v sysctl >/dev/null 2>&1; then
        sysctl -n hw.ncpu
    else
        echo 1
    fi
}

# Downloading, converting and verifying the datasets is prepare_data.sh's job.
# Every step there skips work that is already current, so this costs a few
# seconds on a repeat run -- and it exits non-zero if the ground truth is
# incomplete, which stops a benchmark from quietly scoring unconverted roots as
# "no collision" the way cloth-funnel's half-converted oracle did.
"${SCRIPT_DIR}/prepare_data.sh"

datasets=()
is_enabled "${SCCD_ENABLE_ARMADILLO_ROLLERS}" && datasets+=("armadillo-rollers")
is_enabled "${SCCD_ENABLE_CLOTH_BALL}" && datasets+=("cloth-ball")
is_enabled "${SCCD_ENABLE_CLOTH_FUNNEL}" && datasets+=("cloth-funnel")
is_enabled "${SCCD_ENABLE_N_BODY_SIMULATION}" && datasets+=("n-body-simulation")
is_enabled "${SCCD_ENABLE_PUFFER_BALL}" && datasets+=("puffer-ball")
is_enabled "${SCCD_ENABLE_ROD_TWIST}" && datasets+=("rod-twist")

cmake_bench_args=(-DCMAKE_BUILD_TYPE=Release -DSCCD_ENABLE_SMESH=ON -DSCCD_ENABLE_OPENMP=ON -DSCCD_ENABLE_TIGHT_INCLUSION=ON)
if [[ -n "${SCCD_SMESH_DIR:-}" ]]; then
    cmake_bench_args+=("-Dsmesh_DIR=${SCCD_SMESH_DIR}")
elif [[ -n "${smesh_DIR:-}" ]]; then
    cmake_bench_args+=("-Dsmesh_DIR=${smesh_DIR}")
elif [[ -f "${ROOT_DIR}/build_release/CMakeCache.txt" ]]; then
    cached_smesh_dir="$(sed -n 's/^smesh_DIR[^=]*=//p' "${ROOT_DIR}/build_release/CMakeCache.txt" | tail -n 1)"
    if [[ -n "${cached_smesh_dir}" && -d "${cached_smesh_dir}" ]]; then
        cmake_bench_args+=("-Dsmesh_DIR=${cached_smesh_dir}")
    fi
fi

cmake -S "${ROOT_DIR}" -B "${SCCD_BUILD_DIR}" "${cmake_bench_args[@]}"
cmake --build "${SCCD_BUILD_DIR}" --config Release --target sccd_bench --parallel "$(parallel_jobs)"
# The accuracy oracle is gated on TightInclusion, which this build enables.
cmake --build "${SCCD_BUILD_DIR}" --config Release --target ti_oracle --parallel "$(parallel_jobs)" || true

if [[ -z "${SCCD_DB_TO_RAW:-}" ]]; then
    if command -v db_to_raw >/dev/null 2>&1; then
        export SCCD_DB_TO_RAW="$(command -v db_to_raw)"
    elif [[ -n "${SCCD_SMESH_DIR:-}" && -x "${SCCD_SMESH_DIR%/}/../../../bin/db_to_raw" ]]; then
        export SCCD_DB_TO_RAW="$(cd "${SCCD_SMESH_DIR%/}/../../../bin" && pwd)/db_to_raw"
    elif [[ -n "${smesh_DIR:-}" && -x "${smesh_DIR%/}/../../../bin/db_to_raw" ]]; then
        export SCCD_DB_TO_RAW="$(cd "${smesh_DIR%/}/../../../bin" && pwd)/db_to_raw"
    elif [[ -n "${cached_smesh_dir:-}" && -x "${cached_smesh_dir%/}/../../../bin/db_to_raw" ]]; then
        export SCCD_DB_TO_RAW="$(cd "${cached_smesh_dir%/}/../../../bin" && pwd)/db_to_raw"
    fi
fi

SCCD_BENCH="${SCCD_BUILD_DIR}/sccd_bench"
if [[ ! -x "${SCCD_BENCH}" && -x "${SCCD_BUILD_DIR}/Release/sccd_bench" ]]; then
    SCCD_BENCH="${SCCD_BUILD_DIR}/Release/sccd_bench"
fi

if [[ -z "${OMP_NUM_THREADS:-}" ]]; then
    export OMP_NUM_THREADS="$(parallel_jobs)"
fi
if [[ -z "${OMP_PROC_BIND:-}" ]]; then
    export OMP_PROC_BIND=true
fi

exec 1>&3

# Everything this run produces lands in one directory, so it can be inspected,
# archived or deleted as a unit and does not accumulate in the source tree
# alongside the harness that wrote it. Every path below is still individually
# overridable; SCCD_BENCH_OUT_DIR just moves the default set.
BENCH_OUT_DIR="${SCCD_BENCH_OUT_DIR:-"${BENCHMARK_DIR}/out"}"

BENCH_CSV="${SCCD_BENCH_CSV:-"${BENCH_OUT_DIR}/bench.csv"}"
BENCH_MISSING_PAIRS_CSV="${SCCD_MISSING_PAIRS_CSV:-"${BENCH_OUT_DIR}/bench_missing_pairs.csv"}"
mkdir -p "$(dirname "${BENCH_CSV}")" "$(dirname "${BENCH_MISSING_PAIRS_CSV}")"

# The header comes from the driver, never from here.
#
# This used to be a hardcoded 13-column string while sccd_bench printed 25, and
# the run appended the driver's rows with `tail -n +2`. Columns 14-25 --
# narrow_ms_s1 and every toi_* accuracy column -- therefore landed unnamed, and
# csv.DictReader filed them under the None key where no report ever saw them.
# The driver prints its header before it touches a dataset, so invoking it with
# no dataset arguments yields exactly the current header and nothing else.
bench_header() {
    "${SCCD_BENCH}" --header 2>/dev/null
}

mode_label() {
    case "$1" in
        0) echo "relaxed" ;;
        2) echo "tight" ;;
        *) echo "mode$1" ;;
    esac
}

BENCH_HEADER="$(bench_header)"
if [[ -z "${BENCH_HEADER}" ]]; then
    printf 'error: %s printed no header; refusing to write a CSV with no schema\n' "${SCCD_BENCH}" >&2
    exit 1
fi

if [[ "${#datasets[@]}" -gt 0 ]]; then
    : > "${BENCH_CSV}"
    printf '%s\n' "${BENCH_HEADER}" >> "${BENCH_CSV}"
    for mode in ${SCCD_BENCH_MODES}; do
        printf '==> narrow-phase mode %s (%s)\n' "${mode}" "$(mode_label "${mode}")" >&2
        # Each mode gets its own process so one mode's allocator and cache state
        # cannot colour the next one's timings.
        SCCD_NARROWPHASE_MODE="${mode}" \
        SCCD_MISSING_PAIRS_CSV="${BENCH_MISSING_PAIRS_CSV}" \
            "${SCCD_BENCH}" "${DATA_DIR}" ${datasets[@]+"${datasets[@]}"} | tail -n +2 >> "${BENCH_CSV}"
    done
    cat "${BENCH_CSV}"
else
    printf '%s\n' "${BENCH_HEADER}" | tee "${BENCH_CSV}"
    printf 'dataset,case,type,phase,query_id,c0,c1\n' > "${BENCH_MISSING_PAIRS_CSV}"
fi

exec 1>&2

# --- accuracy against TightInclusion, per mode -----------------------------
# Timing alone cannot tell you whether a mode is safe to use. The oracle checks
# every query against TightInclusion and fails when a mode misses a collision or
# reports a time of impact after the true one.
TI_ORACLE="${SCCD_BUILD_DIR}/ti_oracle"
if [[ ! -x "${TI_ORACLE}" && -x "${SCCD_BUILD_DIR}/Release/ti_oracle" ]]; then
    TI_ORACLE="${SCCD_BUILD_DIR}/Release/ti_oracle"
fi

ORACLE_DIR="${SCCD_BENCH_ORACLE_DIR:-"${BENCH_OUT_DIR}/oracle"}"
ORACLE_CSV="${ORACLE_DIR}/oracle.csv"
if [[ -x "${TI_ORACLE}" ]]; then
    mkdir -p "${ORACLE_DIR}"
    : > "${ORACLE_CSV}"
    oracle_header_written=0
    for dataset in ${datasets[@]+"${datasets[@]}"}; do
        [[ -d "${DATA_DIR}/${dataset}/queries" ]] || continue
        per_dataset="${ORACLE_DIR}/${dataset}.csv"
        # --no-strict: collect every dataset before deciding pass or fail.
        "${TI_ORACLE}" "${DATA_DIR}/${dataset}" \
            --no-strict \
            ${SCCD_ORACLE_MAX_FILES:+--max-files "${SCCD_ORACLE_MAX_FILES}"} \
            --csv "${per_dataset}" \
            --violations-csv "${ORACLE_DIR}/${dataset}-violations.csv" || true
        if [[ -f "${per_dataset}" ]]; then
            if [[ "${oracle_header_written}" -eq 0 ]]; then
                cat "${per_dataset}" >> "${ORACLE_CSV}"
                oracle_header_written=1
            else
                tail -n +2 "${per_dataset}" >> "${ORACLE_CSV}"
            fi
        fi
    done
else
    printf 'note: ti_oracle was not built; skipping the accuracy comparison\n' >&2
fi

# --- report -----------------------------------------------------------------
# One module over the two CSVs this run produced: figures as PDF and PNG,
# booktabs tables, and the Markdown summary. It exits non-zero if any mode
# reported a late time of impact, so a bad run fails here rather than being
# read off a table by someone who did not think to look.
( cd "${BENCHMARK_DIR}" && "${PYTHON}" -m report \
    "${BENCH_CSV}" "${BENCH_OUT_DIR}/report" "${ORACLE_CSV}" ) \
  && printf 'report: %s\n' "${BENCH_OUT_DIR}/report/summary.md" >&2

# Everything needed to reproduce and to read the run: the two CSVs it was
# computed from, and the report generated off them.
BENCH_ARCHIVE="${BENCH_OUT_DIR}/sccd-benchmark-$(date +%Y-%m-%d).tar.gz"
tar -czf "${BENCH_ARCHIVE}" \
    -C "$(dirname "${BENCH_CSV}")" "$(basename "${BENCH_CSV}")" \
    -C "$(dirname "${BENCH_MISSING_PAIRS_CSV}")" "$(basename "${BENCH_MISSING_PAIRS_CSV}")" \
    -C "${BENCH_OUT_DIR}" "$(basename "${ORACLE_DIR}")" \
    -C "${BENCH_OUT_DIR}" report
printf 'archive: %s\n' "${BENCH_ARCHIVE}" >&2
