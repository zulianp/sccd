#!/usr/bin/env bash
set -euo pipefail

# Prepare the NYU benchmark datasets: download, extract, convert, verify.
#
# This used to live inside bench.sh, which meant the only way to prepare data was
# to start a benchmark, and the only way to find out whether the data was any
# good was to read the numbers that came out. Preparation is slow, incremental
# and worth doing once; benchmarking is fast and worth doing often. They are
# separate scripts for that reason.
#
# Every conversion step skips work that is already current, so an interrupted run
# is resumed by running it again -- no flag, and nothing is thrown away.
#
# Usage:
#   prepare_data.sh                 download, convert and verify the enabled scenes
#   prepare_data.sh --verify        verify what is already on disk; convert nothing
#   prepare_data.sh --all           enable all six scenes (large: puffer-ball
#                                   alone is over a gigabyte of queries)
#
# Scenes are selected with the same SCCD_ENABLE_* variables bench.sh uses, so
#   SCCD_ENABLE_CLOTH_BALL=1 SCCD_ENABLE_ARMADILLO_ROLLERS=0 prepare_data.sh
# prepares one scene.

VERIFY_ONLY=0
ENABLE_ALL=0
CHECK_CSV=0
for arg in "$@"; do
    case "${arg}" in
        --verify) VERIFY_ONLY=1 ;;
        --all) ENABLE_ALL=1 ;;
        --check-csv) CHECK_CSV=1 ;;
        -h|--help) sed -n '3,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) printf 'error: unknown argument %s\n' "${arg}" >&2; exit 2 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${BENCHMARK_DIR}/.." && pwd)"
DATA_DIR="${SCCD_DATA_DIR:-"${ROOT_DIR}/data"}"
PYTHON_DIR="${ROOT_DIR}/python"
PYTHON="${PYTHON:-python3}"
JSON_PROJECT_DIR="${BENCHMARK_DIR}/json"
JSON_BUILD_DIR="${SCCD_JSON_BUILD_DIR:-"${ROOT_DIR}/build_json"}"

# The three scenes with complete, verified ground truth are on by default. The
# other three are opt-in: rod-twist and puffer-ball are large downloads, and
# puffer-ball additionally has no runnable case until its frames are extracted.
default_on=1
if [[ "${ENABLE_ALL}" -eq 1 ]]; then default_on=1; fi
export SCCD_ENABLE_ARMADILLO_ROLLERS="${SCCD_ENABLE_ARMADILLO_ROLLERS:-${default_on}}"
export SCCD_ENABLE_CLOTH_BALL="${SCCD_ENABLE_CLOTH_BALL:-${default_on}}"
export SCCD_ENABLE_CLOTH_FUNNEL="${SCCD_ENABLE_CLOTH_FUNNEL:-${default_on}}"
if [[ "${ENABLE_ALL}" -eq 1 ]]; then
    export SCCD_ENABLE_N_BODY_SIMULATION="${SCCD_ENABLE_N_BODY_SIMULATION:-1}"
    export SCCD_ENABLE_PUFFER_BALL="${SCCD_ENABLE_PUFFER_BALL:-1}"
    export SCCD_ENABLE_ROD_TWIST="${SCCD_ENABLE_ROD_TWIST:-1}"
else
    export SCCD_ENABLE_N_BODY_SIMULATION="${SCCD_ENABLE_N_BODY_SIMULATION:-0}"
    export SCCD_ENABLE_PUFFER_BALL="${SCCD_ENABLE_PUFFER_BALL:-0}"
    export SCCD_ENABLE_ROD_TWIST="${SCCD_ENABLE_ROD_TWIST:-0}"
fi
export SCCD_SKIP_DOWNLOAD="${SCCD_SKIP_DOWNLOAD:-0}"

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

scenes=()
is_enabled "${SCCD_ENABLE_ARMADILLO_ROLLERS}" && scenes+=("armadillo-rollers")
is_enabled "${SCCD_ENABLE_CLOTH_BALL}" && scenes+=("cloth-ball")
is_enabled "${SCCD_ENABLE_CLOTH_FUNNEL}" && scenes+=("cloth-funnel")
is_enabled "${SCCD_ENABLE_N_BODY_SIMULATION}" && scenes+=("n-body-simulation")
is_enabled "${SCCD_ENABLE_PUFFER_BALL}" && scenes+=("puffer-ball")
is_enabled "${SCCD_ENABLE_ROD_TWIST}" && scenes+=("rod-twist")

if [[ "${#scenes[@]}" -eq 0 ]]; then
    echo "error: every scene is disabled; nothing to prepare" >&2
    exit 2
fi

verify() {
    local args=("${DATA_DIR}")
    [[ "${CHECK_CSV}" -eq 1 ]] && args=(--check-csv "${args[@]}")
    "${PYTHON}" "${BENCHMARK_DIR}/verify_oracle.py" "${args[@]}" "${scenes[@]}"
}

if [[ "${VERIFY_ONLY}" -eq 1 ]]; then
    verify
    exit $?
fi

printf 'preparing: %s\n' "${scenes[*]}"

# --- download -------------------------------------------------------------
if is_enabled "${SCCD_SKIP_DOWNLOAD}"; then
    printf 'note: SCCD_SKIP_DOWNLOAD is set; assuming %s is already populated\n' \
        "${DATA_DIR}" >&2
else
    "${SCRIPT_DIR}/download_datasets.sh"
fi

# --- the two JSON converters ----------------------------------------------
# A CMake cache remembers the source directory it was generated from and refuses
# to be reused with another -- a hard error, not a reconfigure. Drop a build tree
# that points elsewhere rather than making the caller work that out from CMake's
# "does not match the source used to generate cache".
if [[ -f "${JSON_BUILD_DIR}/CMakeCache.txt" ]]; then
    cached_home="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' \
        "${JSON_BUILD_DIR}/CMakeCache.txt" | tail -n 1)"
    if [[ -n "${cached_home}" && "${cached_home}" != "${JSON_PROJECT_DIR}" ]]; then
        printf 'note: %s was configured from %s; reconfiguring\n' \
            "${JSON_BUILD_DIR}" "${cached_home}" >&2
        rm -rf "${JSON_BUILD_DIR}"
    fi
fi

cmake -S "${JSON_PROJECT_DIR}" -B "${JSON_BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${JSON_BUILD_DIR}" --config Release \
    --target boxes_json_to_raw mma_bool_json_to_raw --parallel "$(parallel_jobs)"

pick_tool() {
    local name="$1"
    if [[ -x "${JSON_BUILD_DIR}/${name}" ]]; then
        printf '%s\n' "${JSON_BUILD_DIR}/${name}"
    elif [[ -x "${JSON_BUILD_DIR}/Release/${name}" ]]; then
        printf '%s\n' "${JSON_BUILD_DIR}/Release/${name}"
    else
        printf 'error: %s was not built\n' "${name}" >&2
        return 1
    fi
}

BOXES_JSON_TO_RAW="$(pick_tool boxes_json_to_raw)"
MMA_BOOL_JSON_TO_RAW="$(pick_tool mma_bool_json_to_raw)"

convert_json() {
    local dir="$1" pattern="$2" tool="$3"
    [[ -d "${dir}" ]] || return 0
    find "${dir}" -maxdepth 1 -name "${pattern}" -print0 | xargs -0 sh -c '
        if [ "$#" -gt 0 ]; then
            "$0" "$@"
        fi
    ' "${tool}"
}

for scene in "${scenes[@]}"; do
    printf '  %s: boxes\n' "${scene}"
    convert_json "${DATA_DIR}/${scene}/boxes" '*.json' "${BOXES_JSON_TO_RAW}"
    printf '  %s: mma_bool\n' "${scene}"
    convert_json "${DATA_DIR}/${scene}/mma_bool" '*_mma_bool.json' "${MMA_BOOL_JSON_TO_RAW}"
done

# --- roots ----------------------------------------------------------------
# The slow one: every stale archive is deserialized and its roots evaluated
# through sympy. Converted archives are skipped on mtime, so this resumes.
printf '  roots (skipping what is current)\n'
"${PYTHON}" "${BENCHMARK_DIR}/roots_to_raw.py" "${DATA_DIR}" "${PYTHON_DIR}" "${scenes[@]}"

# cloth-funnel's PLY headers carry non-ASCII bytes that the frame reader rejects.
if is_enabled "${SCCD_ENABLE_CLOTH_FUNNEL}" \
        && compgen -G "${DATA_DIR}/cloth-funnel/frames/*.ply" >/dev/null; then
    "${PYTHON}" "${PYTHON_DIR}/sccd_strip_nonascii.py" "${DATA_DIR}"/cloth-funnel/frames/*.ply
fi

# --- verify ---------------------------------------------------------------
echo
verify
