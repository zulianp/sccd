#!/usr/bin/env bash
set -euo pipefail

# Run the benchmark as a set of resumable chunks.
#
# bench.sh truncates its CSV and re-runs everything, so an interrupted sweep
# loses all of it. That is the wrong shape for this cluster: the debug partition
# caps a job at 30 minutes, `normal` has no idle nodes, and the account's QOS
# runs one job at a time -- so a full sweep is necessarily many short jobs, and
# any of them can be interrupted.
#
# So the sweep is cut into chunks keyed by (scene, space, repeat, case range).
# Each chunk writes one CSV, and a chunk whose CSV already exists is skipped.
# Re-invoking after an interruption resumes; nothing completed is recomputed and
# nothing is truncated. A chunk writes through a temporary and renames, so a CSV
# that exists is a CSV that finished.
#
# Repeats are separate invocations on purpose: the driver measures each case once
# per run, so run-to-run spread is only visible across processes. It is not
# cosmetic -- the Earliest-output time of impact genuinely varies between runs of
# the same binary, because the parallel search prunes against a running best.
#
# Usage:
#   sweep.sh                      run every chunk that is not already done
#   sweep.sh --dry-run            list the chunks and what is outstanding
#   sweep.sh --merge              concatenate the chunk CSVs and exit
#   sweep.sh --local              run chunks here instead of submitting them
#
# Options:
#   --out DIR        chunk tree              (default benchmark/sweep)
#   --repeats N      independent runs        (default 5)
#   --chunk N        cases per chunk, 0 = whole scene   (default 0)
#   --scenes "..."   scenes to sweep         (default the three verified ones)
#   --spaces "..."   host and/or device      (default host)
#   --modes "..."    narrow-phase modes      (default "0 2")
#   --time HH:MM:SS  per-job limit           (default 00:29:00)
#   --partition P    Slurm partition         (default debug)
#   --account A      Slurm account           (default c40)

OUT_DIR=""
REPEATS=5
CHUNK=0
SCENES="armadillo-rollers cloth-ball cloth-funnel"
SPACES="host"
MODES="0 2"
TIME_LIMIT="00:29:00"
PARTITION="debug"
ACCOUNT="c40"
UENV="prgenv-gnu/24.11:v2"
DRY_RUN=0
MERGE_ONLY=0
LOCAL=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out) OUT_DIR="$2"; shift 2 ;;
        --repeats) REPEATS="$2"; shift 2 ;;
        --chunk) CHUNK="$2"; shift 2 ;;
        --scenes) SCENES="$2"; shift 2 ;;
        --spaces) SPACES="$2"; shift 2 ;;
        --modes) MODES="$2"; shift 2 ;;
        --time) TIME_LIMIT="$2"; shift 2 ;;
        --partition) PARTITION="$2"; shift 2 ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --uenv) UENV="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --merge) MERGE_ONLY=1; shift ;;
        --local) LOCAL=1; shift ;;
        -h|--help) sed -n '3,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) printf 'error: unknown argument %s\n' "$1" >&2; exit 2 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARK_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
ROOT_DIR="$(cd "${BENCHMARK_DIR}/.." && pwd)"
DATA_DIR="${SCCD_DATA_DIR:-"${ROOT_DIR}/data"}"
BUILD_DIR="${SCCD_BUILD_DIR:-"${ROOT_DIR}/build_benchmark"}"
OUT_DIR="${OUT_DIR:-"${BENCHMARK_DIR}/sweep"}"
SCCD_BENCH="${SCCD_BENCH:-"${BUILD_DIR}/sccd_bench"}"

if [[ ! -x "${SCCD_BENCH}" ]]; then
    printf 'error: %s is not executable; build the sccd_bench target first\n' \
        "${SCCD_BENCH}" >&2
    exit 1
fi

header() { "${SCCD_BENCH}" --header; }

# The number of runnable cases in a scene: a boxes/<key>/ with its converted
# arrays and a matching queries/<key>.csv, which is what the driver counts.
case_count() {
    local scene="$1" n=0 key
    local boxes="${DATA_DIR}/${scene}/boxes"
    [[ -d "${boxes}" ]] || { echo 0; return; }
    for dir in "${boxes}"/*/; do
        [[ -d "${dir}" ]] || continue
        key="$(basename "${dir}")"
        [[ -f "${dir}/c0.int32" && -f "${dir}/c1.int32" ]] || continue
        [[ -f "${DATA_DIR}/${scene}/queries/${key}.csv" ]] || continue
        case "${key}" in *vf|*ee) n=$((n + 1)) ;; esac
    done
    echo "${n}"
}

# --- enumerate the chunks --------------------------------------------------
chunk_keys=()
chunk_paths=()
for scene in ${SCENES}; do
    total="$(case_count "${scene}")"
    if [[ "${total}" -eq 0 ]]; then
        printf 'note: %s has no runnable case; skipping\n' "${scene}" >&2
        continue
    fi
    span="${CHUNK}"
    [[ "${span}" -le 0 ]] && span="${total}"
    for space in ${SPACES}; do
        for ((r = 1; r <= REPEATS; ++r)); do
            for ((begin = 0; begin < total; begin += span)); do
                end=$((begin + span))
                [[ "${end}" -gt "${total}" ]] && end="${total}"
                chunk_keys+=("${scene}|${space}|${r}|${begin}|${end}")
                chunk_paths+=("${OUT_DIR}/${scene}/${space}/r${r}/$(printf '%06d-%06d' "${begin}" "${end}").csv")
            done
        done
    done
done

merge() {
    local merged="${OUT_DIR}/bench.csv"
    mkdir -p "${OUT_DIR}"
    header > "${merged}"
    local n=0
    while IFS= read -r -d '' f; do
        tail -n +2 "${f}" >> "${merged}"
        n=$((n + 1))
    done < <(find "${OUT_DIR}" -name '*.csv' -not -name 'bench.csv' -print0 | sort -z)
    printf 'merged %d chunk files -> %s (%d rows)\n' \
        "${n}" "${merged}" "$(( $(wc -l < "${merged}") - 1 ))"
}

if [[ "${MERGE_ONLY}" -eq 1 ]]; then
    merge
    exit 0
fi

todo=()
done_count=0
for i in "${!chunk_keys[@]}"; do
    if [[ -s "${chunk_paths[$i]}" ]]; then
        done_count=$((done_count + 1))
    else
        todo+=("$i")
    fi
done

printf '%d chunks: %d done, %d outstanding\n' \
    "${#chunk_keys[@]}" "${done_count}" "${#todo[@]}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
    for i in "${todo[@]:-}"; do
        [[ -n "${i}" ]] || continue
        IFS='|' read -r scene space r begin end <<< "${chunk_keys[$i]}"
        printf '  TODO %-20s %-6s r%-2s cases [%s, %s)\n' \
            "${scene}" "${space}" "${r}" "${begin}" "${end}"
    done
    exit 0
fi

if [[ "${#todo[@]}" -eq 0 ]]; then
    printf 'nothing to do; use --merge to assemble the CSV\n'
    exit 0
fi

# --- run one chunk ---------------------------------------------------------
# Every mode goes into the same chunk file, each in its own process so one
# mode's allocator and cache state cannot colour the next one's timings.
run_chunk_body() {
    local scene="$1" space="$2" begin="$3" end="$4" out="$5"
    local tmp="${out}.partial"
    mkdir -p "$(dirname "${out}")"
    header > "${tmp}"
    for mode in ${MODES}; do
        SCCD_NARROWPHASE_MODE="${mode}" \
        SCCD_BENCH_EXECUTION_SPACE="${space}" \
        SCCD_BENCH_CASE_BEGIN="${begin}" \
        SCCD_BENCH_CASE_END="${end}" \
            "${SCCD_BENCH}" "${DATA_DIR}" "${scene}" | tail -n +2 >> "${tmp}"
    done
    mv "${tmp}" "${out}"
}

failures=0
for i in "${todo[@]}"; do
    IFS='|' read -r scene space r begin end <<< "${chunk_keys[$i]}"
    out="${chunk_paths[$i]}"
    label="${scene}/${space}/r${r}/${begin}-${end}"

    if [[ "${LOCAL}" -eq 1 ]]; then
        printf '==> %s\n' "${label}"
        if ! run_chunk_body "${scene}" "${space}" "${begin}" "${end}" "${out}"; then
            printf 'FAILED %s\n' "${label}" >&2
            failures=$((failures + 1))
        fi
        continue
    fi

    # One job at a time: the account's QOS will not run two, so submitting the
    # whole sweep at once only fills the queue. --wait blocks until this chunk
    # finishes, and because a finished chunk is skipped on the next invocation,
    # killing the sweep here costs at most the chunk in flight.
    printf '==> submitting %s\n' "${label}"
    gpu_args=()
    [[ "${space}" == "device" ]] && gpu_args=(--gpus-per-task=1)
    if ! sbatch --wait \
            --account="${ACCOUNT}" --partition="${PARTITION}" \
            --nodes=1 --ntasks=1 "${gpu_args[@]}" \
            --time="${TIME_LIMIT}" \
            --uenv="${UENV}" --view=default \
            --job-name="sccd-${scene}-r${r}" \
            --output="${OUT_DIR}/logs/${scene}-${space}-r${r}-${begin}.out" \
            --error="${OUT_DIR}/logs/${scene}-${space}-r${r}-${begin}.err" \
            --wrap="$(declare -f run_chunk_body); \
                     SCCD_BENCH='${SCCD_BENCH}' DATA_DIR='${DATA_DIR}' MODES='${MODES}'; \
                     header() { '${SCCD_BENCH}' --header; }; \
                     run_chunk_body '${scene}' '${space}' '${begin}' '${end}' '${out}'"; then
        printf 'FAILED %s (see %s/logs)\n' "${label}" "${OUT_DIR}" >&2
        failures=$((failures + 1))
    fi
done

printf '\n%d chunks attempted, %d failed\n' "${#todo[@]}" "${failures}"
merge
[[ "${failures}" -eq 0 ]]
