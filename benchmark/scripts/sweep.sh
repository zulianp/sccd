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
#   --no-oracle      skip the accuracy stage (timings only)

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
ORACLE=1
TI_ORACLE=""

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
        --no-oracle) ORACLE=0; shift ;;
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
TI_ORACLE="${SCCD_TI_ORACLE:-"${BUILD_DIR}/ti_oracle"}"

if [[ ! -x "${SCCD_BENCH}" ]]; then
    printf 'error: %s is not executable; build the sccd_bench target first\n' \
        "${SCCD_BENCH}" >&2
    exit 1
fi

# The driver is the only authority on its own schema, so the header is asked of
# it rather than duplicated here. On Alps the binary links CUDA and cannot even
# print a string outside its uenv, so fall back to asking it through one before
# giving up -- a wrong header silently misfiles every accuracy column.
header() {
    local out
    if out="$("${SCCD_BENCH}" --header 2>/dev/null)" && [[ -n "${out}" ]]; then
        printf '%s\n' "${out}"
        return 0
    fi
    if command -v uenv >/dev/null 2>&1; then
        if out="$(uenv run --view=default "${UENV}" -- "${SCCD_BENCH}" --header 2>/dev/null)" \
                && [[ -n "${out}" ]]; then
            printf '%s\n' "${out}"
            return 0
        fi
    fi
    printf 'error: %s printed no header; refusing to write a CSV with no schema\n' \
        "${SCCD_BENCH}" >&2
    return 1
}

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

# The accuracy stage. ti_oracle runs every mode over every query of a scene and
# checks each answer against the dataset's exact roots, so it is one job per
# (scene, repeat) rather than per case range -- it has no case selector, and its
# own timing column is what makes TightInclusion a timed reference rather than
# only a correctness oracle.
oracle_keys=()
oracle_paths=()
if [[ "${ORACLE}" -eq 1 ]]; then
    for scene in ${SCENES}; do
        [[ -d "${DATA_DIR}/${scene}/queries" ]] || continue
        for ((r = 1; r <= REPEATS; ++r)); do
            oracle_keys+=("${scene}|${r}")
            oracle_paths+=("${OUT_DIR}/oracle/${scene}/r${r}.csv")
        done
    done
fi

merge() {
    local merged="${OUT_DIR}/bench.csv"
    mkdir -p "${OUT_DIR}"
    header > "${merged}"
    local n=0
    while IFS= read -r -d '' f; do
        tail -n +2 "${f}" >> "${merged}"
        n=$((n + 1))
    done < <(find "${OUT_DIR}" -path "${OUT_DIR}/oracle" -prune -o \
                  -name '*.csv' -not -name 'bench.csv' -not -name 'oracle.csv' \
                  -print0 | sort -z)
    printf 'merged %d chunk files -> %s (%d rows)\n' \
        "${n}" "${merged}" "$(( $(wc -l < "${merged}") - 1 ))"

    # The oracle CSV has its own schema and is merged into its own file rather
    # than forced into the timing one.
    local oracle_merged="${OUT_DIR}/oracle.csv"
    local first=1 m=0
    for f in "${OUT_DIR}"/oracle/*/r*.csv; do
        [[ -s "${f}" ]] || continue
        if [[ "${first}" -eq 1 ]]; then
            head -1 "${f}" > "${oracle_merged}"
            first=0
        fi
        tail -n +2 "${f}" >> "${oracle_merged}"
        m=$((m + 1))
    done
    if [[ "${m}" -gt 0 ]]; then
        printf 'merged %d oracle files -> %s (%d rows)\n' \
            "${m}" "${oracle_merged}" "$(( $(wc -l < "${oracle_merged}") - 1 ))"
    fi
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

oracle_todo=()
oracle_done=0
for i in "${!oracle_keys[@]}"; do
    if [[ -s "${oracle_paths[$i]}" ]]; then
        oracle_done=$((oracle_done + 1))
    else
        oracle_todo+=("$i")
    fi
done

printf '%d timing chunks: %d done, %d outstanding\n' \
    "${#chunk_keys[@]}" "${done_count}" "${#todo[@]}"
if [[ "${#oracle_keys[@]}" -gt 0 ]]; then
    printf '%d accuracy runs: %d done, %d outstanding\n' \
        "${#oracle_keys[@]}" "${oracle_done}" "${#oracle_todo[@]}"
fi

if [[ "${DRY_RUN}" -eq 1 ]]; then
    for i in "${todo[@]:-}"; do
        [[ -n "${i}" ]] || continue
        IFS='|' read -r scene space r begin end <<< "${chunk_keys[$i]}"
        printf '  TODO timing   %-20s %-6s r%-2s cases [%s, %s)\n' \
            "${scene}" "${space}" "${r}" "${begin}" "${end}"
    done
    for i in "${oracle_todo[@]:-}"; do
        [[ -n "${i}" ]] || continue
        IFS='|' read -r scene r <<< "${oracle_keys[$i]}"
        printf '  TODO accuracy %-20s        r%-2s\n' "${scene}" "${r}"
    done
    exit 0
fi

if [[ "${#todo[@]}" -eq 0 && "${#oracle_todo[@]}" -eq 0 ]]; then
    printf 'nothing to do; use --merge to assemble the CSV\n'
    exit 0
fi

# --- run one chunk ---------------------------------------------------------
# Every mode goes into the same chunk file, each in its own process so one
# mode's allocator and cache state cannot colour the next one's timings.
run_oracle_body() {
    local scene="$1" out="$2"
    local tmp="${out}.partial"
    mkdir -p "$(dirname "${out}")"
    # --no-strict so one scene's violation does not abort the sweep before the
    # rest is measured. The violations are in the CSV either way, and the report
    # fails on them; losing the remaining scenes as well helps nobody.
    "${TI_ORACLE}" "${DATA_DIR}/${scene}" --csv "${tmp}" --no-strict >/dev/null 2>&1 || true
    [[ -s "${tmp}" ]] || return 1
    mv "${tmp}" "${out}"
}

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

# The QOS caps how many jobs one user may have submitted at once, and this
# account is shared with the user's own work. A rejection on that ground means
# "the queue is full right now", not "this chunk failed", so wait for a slot
# instead of burning the chunk. Anything else is a real error.
submit_with_retry() {
    local attempt=0 err
    while :; do
        err="$(sbatch --wait "$@" 2>&1 >/dev/null)" && return 0
        if [[ "${err}" != *QOSMaxSubmitJobPerUserLimit* ]]; then
            printf '%s\n' "${err}" >&2
            return 1
        fi
        attempt=$((attempt + 1))
        if [[ "${attempt}" -gt "${SWEEP_SUBMIT_RETRIES:-240}" ]]; then
            printf 'giving up: the submit limit stayed full for %s attempts\n' \
                "${attempt}" >&2
            return 1
        fi
        if [[ "${attempt}" -eq 1 ]]; then
            printf '    queue full (submit limit); waiting for a slot\n' >&2
        fi
        sleep "${SWEEP_SUBMIT_WAIT:-60}"
    done
}

failures=0
for i in "${todo[@]:-}"; do
    [[ -n "${i}" ]] || continue
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

    if ! submit_with_retry \
            --account="${ACCOUNT}" --partition="${PARTITION}" \
            --nodes=1 --ntasks=1 "${gpu_args[@]}" \
            --time="${TIME_LIMIT}" \
            --uenv="${UENV}" --view=default \
            --job-name="sccd-${scene}-r${r}" \
            --output="${OUT_DIR}/logs/${scene}-${space}-r${r}-${begin}.out" \
            --error="${OUT_DIR}/logs/${scene}-${space}-r${r}-${begin}.err" \
            --wrap="$(declare -f run_chunk_body); \
                     export SCCD_BENCH='${SCCD_BENCH}' DATA_DIR='${DATA_DIR}' MODES='${MODES}'; \
                     export OMP_NUM_THREADS=\"\$(nproc)\"; \
                     header() { '${SCCD_BENCH}' --header; }; \
                     run_chunk_body '${scene}' '${space}' '${begin}' '${end}' '${out}'"; then
        printf 'FAILED %s (see %s/logs)\n' "${label}" "${OUT_DIR}" >&2
        failures=$((failures + 1))
    fi
done

for i in "${oracle_todo[@]:-}"; do
    [[ -n "${i}" ]] || continue
    IFS='|' read -r scene r <<< "${oracle_keys[$i]}"
    out="${oracle_paths[$i]}"
    label="accuracy ${scene}/r${r}"

    if [[ "${LOCAL}" -eq 1 ]]; then
        printf '==> %s\n' "${label}"
        if ! run_oracle_body "${scene}" "${out}"; then
            printf 'FAILED %s\n' "${label}" >&2
            failures=$((failures + 1))
        fi
        continue
    fi

    printf '==> submitting %s\n' "${label}"
    if ! submit_with_retry \
            --account="${ACCOUNT}" --partition="${PARTITION}" \
            --nodes=1 --ntasks=1 --gpus-per-task=1 \
            --time="${TIME_LIMIT}" \
            --uenv="${UENV}" --view=default \
            --job-name="sccd-oracle-${scene}-r${r}" \
            --output="${OUT_DIR}/logs/oracle-${scene}-r${r}.out" \
            --error="${OUT_DIR}/logs/oracle-${scene}-r${r}.err" \
            --wrap="$(declare -f run_oracle_body); \
                     export TI_ORACLE='${TI_ORACLE}' DATA_DIR='${DATA_DIR}'; \
                     export OMP_NUM_THREADS=\"\$(nproc)\"; \
                     run_oracle_body '${scene}' '${out}'"; then
        printf 'FAILED %s (see %s/logs)\n' "${label}" "${OUT_DIR}" >&2
        failures=$((failures + 1))
    fi
done

printf '\n%d runs attempted, %d failed\n' \
    "$(( ${#todo[@]} + ${#oracle_todo[@]} ))" "${failures}"
merge
[[ "${failures}" -eq 0 ]]
