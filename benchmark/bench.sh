#!/usr/bin/env bash
set -euo pipefail

# The goal of this file is to perform a wide benchmark of the SCCD library using the NYU datasets
# We benchmarki the broad-phase and narrow-phase collision detection:
# 1) for performance in timing
# 2) for accuracy compared to the ground truth
# 3) Number of false positives (should be as low as possible for the narrow-phase) and negatives (should be 0)

export  SCCD_ENABLE_ARMADILLO_ROLLERS=1
export  SCCD_ENABLE_CLOTH_BALL=0
export  SCCD_ENABLE_CLOTH_FUNNEL=0
export  SCCD_ENABLE_N_BODY_SIMULATION=0
export  SCCD_ENABLE_PUFFER_BALL=0
export  SCCD_ENABLE_ROD_TWIST=0

BENCHMARK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${BENCHMARK_DIR}/download_datasets.sh"

DATA_DIR="${SCCD_DATA_DIR:-"${BENCHMARK_DIR}/../data"}"
JSON_PROJECT_DIR="${BENCHMARK_DIR}/../external/json"
JSON_BUILD_DIR="${SCCD_JSON_BUILD_DIR:-"${BENCHMARK_DIR}/../build_json"}"

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

cmake -S "${JSON_PROJECT_DIR}" -B "${JSON_BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${JSON_BUILD_DIR}" --config Release --target boxes_json_to_raw --parallel "$(parallel_jobs)"

BOXES_JSON_TO_RAW="${JSON_BUILD_DIR}/boxes_json_to_raw"
if [[ ! -x "${BOXES_JSON_TO_RAW}" && -x "${JSON_BUILD_DIR}/Release/boxes_json_to_raw" ]]; then
    BOXES_JSON_TO_RAW="${JSON_BUILD_DIR}/Release/boxes_json_to_raw"
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

# TODO:
# 2) Use read_wxf to generate the impact exact times for each query (e.g., data/armadillo-rollers/roots/0ee_roots.tar) and put them in the dedicated folder (e.g., roots/0ee) as toi.float64


# 3) Do the same for the mma_bool files and put them in the dedicated folder (e.g., mma_bool/0ee) as mma_bool.uint8

# TODO:
# 4) Create a bench.exe.cpp that reads the meshes and scans the folder boxes and reads the raw files, times the CCD for each file collision files pair (names the trace file after the case and folder e.g., SMESH_TRACE_FILE=armadillo-rollers/0ee)
# the timings are collected as milliseconds and stored in a unique raw binary file for the whole case and collision type e.g., armadillo-rollers-fv.float64 and armadillo-rollers-ee.float64 
# The benchmark also collects the accuracy metrics: number of false positives and negatives for the narrow-phase, and number of false positives and negatives for the broad-phase. writes sccd_toi.float64, sccd_fp.uint8, sccd_fn.uint8, sccd_fp_broad.uint8, sccd_fn_broad.uint8, it also makes sure that the data is aligned/ordeded with the query file and the roots file.
