#!/usr/bin/env bash
set -euo pipefail

# The goal of this file is to perform a wide benchmark of the SCCD library using the NYU datasets
# We benchmarki the broad-phase and narrow-phase collision detection:
# 1) for performance in timing
# 2) for accuracy compared to the ground truth
# 3) Number of false positives (should be as low as possible for the narrow-phase) and negatives (should be 0)
# For cuda on Alps: OMP_NUM_THREADS=72 SCCD_BENCH_EXECUTION_SPACE=device srun ./bench.sh 

# Narrow-phase modes to sweep (see src/narrowphase/sccd_narrowphase_mode.hpp):
#   0 fast   2 tight
# Those are the two that exist. The binary is run once per mode and the rows are
# concatenated; every row carries its mode, so the postprocessor keeps the series
# apart. 1 and 3 were retired: setting either warns and runs Fast, so a sweep
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

BENCHMARK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec 3>&1
exec 1>&2

DATA_DIR="${SCCD_DATA_DIR:-"${BENCHMARK_DIR}/../data"}"
JSON_PROJECT_DIR="${BENCHMARK_DIR}/json"
JSON_BUILD_DIR="${SCCD_JSON_BUILD_DIR:-"${BENCHMARK_DIR}/../build_json"}"
SCCD_BUILD_DIR="${SCCD_BUILD_DIR:-"${BENCHMARK_DIR}/../build_benchmark"}"
PYTHON_DIR="${BENCHMARK_DIR}/../python"
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

if is_enabled "${SCCD_SKIP_DOWNLOAD}"; then
    printf 'note: SCCD_SKIP_DOWNLOAD is set; assuming %s is already populated\n' "${DATA_DIR}" >&2
else
    "${BENCHMARK_DIR}/download_datasets.sh"
fi

cmake -S "${JSON_PROJECT_DIR}" -B "${JSON_BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${JSON_BUILD_DIR}" --config Release --target boxes_json_to_raw mma_bool_json_to_raw --parallel "$(parallel_jobs)"

BOXES_JSON_TO_RAW="${JSON_BUILD_DIR}/boxes_json_to_raw"
if [[ ! -x "${BOXES_JSON_TO_RAW}" && -x "${JSON_BUILD_DIR}/Release/boxes_json_to_raw" ]]; then
    BOXES_JSON_TO_RAW="${JSON_BUILD_DIR}/Release/boxes_json_to_raw"
fi
MMA_BOOL_JSON_TO_RAW="${JSON_BUILD_DIR}/mma_bool_json_to_raw"
if [[ ! -x "${MMA_BOOL_JSON_TO_RAW}" && -x "${JSON_BUILD_DIR}/Release/mma_bool_json_to_raw" ]]; then
    MMA_BOOL_JSON_TO_RAW="${JSON_BUILD_DIR}/Release/mma_bool_json_to_raw"
fi

datasets=()
is_enabled "${SCCD_ENABLE_ARMADILLO_ROLLERS}" && datasets+=("armadillo-rollers")
is_enabled "${SCCD_ENABLE_CLOTH_BALL}" && datasets+=("cloth-ball")
is_enabled "${SCCD_ENABLE_CLOTH_FUNNEL}" && datasets+=("cloth-funnel")
is_enabled "${SCCD_ENABLE_N_BODY_SIMULATION}" && datasets+=("n-body-simulation")
is_enabled "${SCCD_ENABLE_PUFFER_BALL}" && datasets+=("puffer-ball")
is_enabled "${SCCD_ENABLE_ROD_TWIST}" && datasets+=("rod-twist")

for dataset in "${datasets[@]}"; do
    boxes_dir="${DATA_DIR}/${dataset}/boxes"
    [[ -d "${boxes_dir}" ]] || continue
    find "${boxes_dir}" -maxdepth 1 -name '*.json' -print0 | xargs -0 sh -c '
        if [ "$#" -gt 0 ]; then
            "$0" "$@"
        fi
    ' "${BOXES_JSON_TO_RAW}"
done

"${PYTHON}" "${BENCHMARK_DIR}/roots_to_raw.py" "${DATA_DIR}" "${PYTHON_DIR}" "${datasets[@]}"

is_enabled "${SCCD_ENABLE_CLOTH_FUNNEL}" && "${PYTHON}" "${PYTHON_DIR}"/sccd_strip_nonascii.py "${DATA_DIR}"/cloth-funnel/frames/*.ply

for dataset in "${datasets[@]}"; do
    mma_bool_dir="${DATA_DIR}/${dataset}/mma_bool"
    [[ -d "${mma_bool_dir}" ]] || continue
    find "${mma_bool_dir}" -maxdepth 1 -name '*_mma_bool.json' -print0 | xargs -0 sh -c '
        if [ "$#" -gt 0 ]; then
            "$0" "$@"
        fi
    ' "${MMA_BOOL_JSON_TO_RAW}"
done

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

BENCH_CSV="${SCCD_BENCH_CSV:-"${BENCHMARK_DIR}/bench.csv"}"
BENCH_AGG_CSV="${SCCD_BENCH_AGG_CSV:-"${BENCHMARK_DIR}/bench_aggregate.csv"}"
BENCH_PAIRED_CSV="${BENCH_AGG_CSV/_aggregate/_paired}"
if [[ "${BENCH_PAIRED_CSV}" == "${BENCH_AGG_CSV}" ]]; then
    BENCH_PAIRED_CSV="$(dirname "${BENCH_AGG_CSV}")/$(basename "${BENCH_AGG_CSV%.*}")_paired.${BENCH_AGG_CSV##*.}"
fi
BENCH_NP_QUERY_TIMING_CSV="${BENCH_AGG_CSV/_aggregate/_np_query_timing}"
if [[ "${BENCH_NP_QUERY_TIMING_CSV}" == "${BENCH_AGG_CSV}" ]]; then
    BENCH_NP_QUERY_TIMING_CSV="$(dirname "${BENCH_AGG_CSV}")/$(basename "${BENCH_AGG_CSV%.*}")_np_query_timing.${BENCH_AGG_CSV##*.}"
fi
BENCH_TOI_ERROR_CSV="${BENCH_AGG_CSV/_aggregate/_toi_error}"
if [[ "${BENCH_TOI_ERROR_CSV}" == "${BENCH_AGG_CSV}" ]]; then
    BENCH_TOI_ERROR_CSV="$(dirname "${BENCH_AGG_CSV}")/$(basename "${BENCH_AGG_CSV%.*}")_toi_error.${BENCH_AGG_CSV##*.}"
fi
BENCH_MISSING_PAIRS_CSV="${SCCD_MISSING_PAIRS_CSV:-"${BENCH_AGG_CSV/_aggregate/_missing_pairs}"}"
if [[ "${BENCH_MISSING_PAIRS_CSV}" == "${BENCH_AGG_CSV}" ]]; then
    BENCH_MISSING_PAIRS_CSV="$(dirname "${BENCH_AGG_CSV}")/$(basename "${BENCH_AGG_CSV%.*}")_missing_pairs.${BENCH_AGG_CSV##*.}"
fi
BENCH_FIGURE_DIR="${SCCD_BENCH_FIGURE_DIR:-"${BENCHMARK_DIR}/figures"}"
BENCH_REPORT_TEX="${SCCD_BENCH_REPORT_TEX:-"${BENCHMARK_DIR}/bench_report.tex"}"
mkdir -p "$(dirname "${BENCH_CSV}")" "$(dirname "${BENCH_AGG_CSV}")" "$(dirname "${BENCH_MISSING_PAIRS_CSV}")" "${BENCH_FIGURE_DIR}" "$(dirname "${BENCH_REPORT_TEX}")"

BENCH_HEADER='dataset,mode,case,type,queries,prep_ms,broad_ms,narrow_ms,query_narrow_ms,fp,fn,broad_fp,broad_fn'

mode_label() {
    case "$1" in
        0) echo "fast" ;;
        2) echo "tight" ;;
        *) echo "mode$1" ;;
    esac
}

if [[ "${#datasets[@]}" -gt 0 ]]; then
    : > "${BENCH_CSV}"
    printf '%s\n' "${BENCH_HEADER}" >> "${BENCH_CSV}"
    for mode in ${SCCD_BENCH_MODES}; do
        printf '==> narrow-phase mode %s (%s)\n' "${mode}" "$(mode_label "${mode}")" >&2
        # Each mode gets its own process so one mode's allocator and cache state
        # cannot colour the next one's timings.
        SCCD_NARROWPHASE_MODE="${mode}" \
        SCCD_MISSING_PAIRS_CSV="${BENCH_MISSING_PAIRS_CSV}" \
            "${SCCD_BENCH}" "${DATA_DIR}" "${datasets[@]}" | tail -n +2 >> "${BENCH_CSV}"
    done
    cat "${BENCH_CSV}"
else
    printf '%s\n' "${BENCH_HEADER}" | tee "${BENCH_CSV}"
    printf 'dataset,case,type,phase,query_id,c0,c1\n' > "${BENCH_MISSING_PAIRS_CSV}"
fi

exec 1>&2

"${PYTHON}" "${BENCHMARK_DIR}/bench_postprocess.py" \
    "${BENCH_CSV}" "${BENCH_AGG_CSV}" "${BENCH_FIGURE_DIR}" "${BENCH_REPORT_TEX}" "${DATA_DIR}"

# --- accuracy against TightInclusion, per mode -----------------------------
# Timing alone cannot tell you whether a mode is safe to use. The oracle checks
# every query against TightInclusion and fails when a mode misses a collision or
# reports a time of impact after the true one.
TI_ORACLE="${SCCD_BUILD_DIR}/ti_oracle"
if [[ ! -x "${TI_ORACLE}" && -x "${SCCD_BUILD_DIR}/Release/ti_oracle" ]]; then
    TI_ORACLE="${SCCD_BUILD_DIR}/Release/ti_oracle"
fi

ORACLE_DIR="${SCCD_BENCH_ORACLE_DIR:-"${BENCHMARK_DIR}/oracle/run"}"
ORACLE_CSV="${ORACLE_DIR}/oracle.csv"
if [[ -x "${TI_ORACLE}" ]]; then
    mkdir -p "${ORACLE_DIR}"
    : > "${ORACLE_CSV}"
    oracle_header_written=0
    for dataset in "${datasets[@]}"; do
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

# --- HTML report ------------------------------------------------------------
BENCH_REPORT_HTML="${SCCD_BENCH_REPORT_HTML:-"${BENCHMARK_DIR}/bench_report.html"}"
"${PYTHON}" "${BENCHMARK_DIR}/bench_report_html.py" \
    --aggregate "${BENCH_AGG_CSV}" \
    --toi-error "${BENCH_TOI_ERROR_CSV}" \
    --oracle "${ORACLE_CSV}" \
    --out "${BENCH_REPORT_HTML}" \
  && printf 'report: %s\n' "${BENCH_REPORT_HTML}" >&2

BENCH_ARCHIVE="${BENCHMARK_DIR}/sccd-benchmark-$(date +%Y-%m-%d).tar.gz"
tar -czf "${BENCH_ARCHIVE}" \
    -C "$(dirname "${BENCH_CSV}")" "$(basename "${BENCH_CSV}")" \
    -C "$(dirname "${BENCH_AGG_CSV}")" "$(basename "${BENCH_AGG_CSV}")" \
    -C "$(dirname "${BENCH_PAIRED_CSV}")" "$(basename "${BENCH_PAIRED_CSV}")" \
    -C "$(dirname "${BENCH_NP_QUERY_TIMING_CSV}")" "$(basename "${BENCH_NP_QUERY_TIMING_CSV}")" \
    -C "$(dirname "${BENCH_TOI_ERROR_CSV}")" "$(basename "${BENCH_TOI_ERROR_CSV}")" \
    -C "$(dirname "${BENCH_MISSING_PAIRS_CSV}")" "$(basename "${BENCH_MISSING_PAIRS_CSV}")" \
    -C "$(dirname "${BENCH_REPORT_TEX}")" "$(basename "${BENCH_REPORT_TEX}")" \
    -C "$(dirname "${BENCH_FIGURE_DIR}")" "$(basename "${BENCH_FIGURE_DIR}")" \
    -C "$(dirname "${BENCH_REPORT_HTML}")" "$(basename "${BENCH_REPORT_HTML}")"
