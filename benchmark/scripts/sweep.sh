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
#   --broadphases ".."  sweep, cell2d and/or auto (default cell2d)
#   --modes "..."    narrow-phase modes      (default "0 2")
#   --time HH:MM:SS  per-job limit           (default 00:29:00)
#   --partition P    Slurm partition         (default debug)
#   --account A      Slurm account           (default c40)
#   --no-oracle      skip the accuracy stage (timings only)
#   --oracle-chunk N query files per accuracy job, 0 = whole scene
#   --pack N         chunks per Slurm job
#   --jobs N         Slurm jobs in flight at once (default 1)    (default 1)
#   --max-cases N    sweep an evenly spread subsample of N cases (default: all)
#
# --pack matters more than it looks. The account runs one job at a time and
# shares the queue with other work, so wall-clock time is dominated by waiting,
# not by computing: thirty one-chunk jobs each wait their turn, while six
# five-chunk jobs wait six times. Pack as much as fits inside --time.

OUT_DIR=""
REPEATS=5
CHUNK=0
SCENES="armadillo-rollers cloth-ball cloth-funnel"
SPACES="host"
MODES="0 2"
# Broad-phase strategies to sample. "auto" races them per scene, which is what
# a caller gets by default; naming them explicitly is what makes the two
# comparable, because a raced run reports whichever won and not which ran.
BROADPHASES="cell2d"
TIME_LIMIT="00:29:00"
PARTITION="debug"
ACCOUNT="c40"
UENV="prgenv-gnu/24.11:v2"
DRY_RUN=0
MERGE_ONLY=0
LOCAL=0
ORACLE=1
ORACLE_CHUNK=0
PACK=1
# How many submitted jobs may be in flight at once. One is right when the
# account's QOS admits a single running job, which is what the debug queue
# amounted to; on a partition with idle nodes it leaves the whole sweep
# serialised behind one worker for no reason.
JOBS=1
MAX_CASES=0
TI_ORACLE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out) OUT_DIR="$2"; shift 2 ;;
        --repeats) REPEATS="$2"; shift 2 ;;
        --chunk) CHUNK="$2"; shift 2 ;;
        --scenes) SCENES="$2"; shift 2 ;;
        --spaces) SPACES="$2"; shift 2 ;;
        --broadphases) BROADPHASES="$2"; shift 2 ;;
        --modes) MODES="$2"; shift 2 ;;
        --time) TIME_LIMIT="$2"; shift 2 ;;
        --partition) PARTITION="$2"; shift 2 ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --uenv) UENV="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --merge) MERGE_ONLY=1; shift ;;
        --local) LOCAL=1; shift ;;
        --no-oracle) ORACLE=0; shift ;;
        --pack) PACK="$2"; shift 2 ;;
        --jobs) JOBS="$2"; shift 2 ;;
        --oracle-chunk) ORACLE_CHUNK="$2"; shift 2 ;;
        --max-cases) MAX_CASES="$2"; shift 2 ;;
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

# The driver converts a scene's PLY frames to raw arrays the first time it needs
# them, by shelling out to smesh's db_to_raw. That tool is not on PATH inside a
# Slurm job, and a scene whose frames_raw/ does not already exist then produces a
# chunk with a header and no rows -- which merges silently into an otherwise
# healthy CSV. Find it the way bench.sh does and pass it into the job.
if [[ -z "${SCCD_DB_TO_RAW:-}" ]]; then
    if command -v db_to_raw >/dev/null 2>&1; then
        SCCD_DB_TO_RAW="$(command -v db_to_raw)"
    else
        for cmake_dir in "${SCCD_SMESH_DIR:-}" "${smesh_DIR:-}" \
                "$(sed -n 's/^smesh_DIR[^=]*=//p' "${BUILD_DIR}/CMakeCache.txt" \
                   2>/dev/null | tail -n 1)"; do
            [[ -n "${cmake_dir}" ]] || continue
            candidate="${cmake_dir%/}/../../../bin/db_to_raw"
            if [[ -x "${candidate}" ]]; then
                SCCD_DB_TO_RAW="$(cd "$(dirname "${candidate}")" && pwd)/db_to_raw"
                break
            fi
        done
    fi
fi
export SCCD_DB_TO_RAW="${SCCD_DB_TO_RAW:-}"

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
    # SCCD_BENCH_MAX_CASES narrows the driver's list before the range is applied,
    # so the chunking has to narrow with it. Without this, every range past N
    # addresses cases that are no longer in the list and yields an empty chunk.
    if [[ "${MAX_CASES}" -gt 0 && "${total}" -gt "${MAX_CASES}" ]]; then
        total="${MAX_CASES}"
    fi
    span="${CHUNK}"
    [[ "${span}" -le 0 ]] && span="${total}"
    for space in ${SPACES}; do
        for bp in ${BROADPHASES}; do
            for ((r = 1; r <= REPEATS; ++r)); do
                for ((begin = 0; begin < total; begin += span)); do
                    end=$((begin + span))
                    [[ "${end}" -gt "${total}" ]] && end="${total}"
                    chunk_keys+=("${scene}|${space}|${bp}|${r}|${begin}|${end}")
                    chunk_paths+=("${OUT_DIR}/${scene}/${space}/${bp}/r${r}/$(printf '%06d-%06d' "${begin}" "${end}").csv")
                done
            done
        done
    done
done

# The accuracy stage. ti_oracle runs every mode over every query of a scene and
# checks each answer against the dataset's exact roots, so it is one job per
# (scene, repeat) rather than per case range -- it has no case selector, and its
# own timing column is what makes TightInclusion a timed reference rather than
# only a correctness oracle.
# ti_oracle applies a file range to each phase's list separately, so the number
# of chunks follows the larger of the two phases.
# wc -l pads its output on some systems, and a command substitution inside an
# arithmetic expansion swallows the result when it does; strip it here instead.
phase_files() {
    local scene="$1" suffix="$2"
    # -L because a prepared data tree may reach its query set through a symlink
    # -- the float64 tree shares one copy of the queries with the float32 tree --
    # and plain find does not descend into one, so the count silently comes back
    # zero and the scene looks like it has no work.
    find -L "${DATA_DIR}/${scene}/queries" -maxdepth 1 -name "*${suffix}.csv" \
        2>/dev/null | wc -l | tr -d ' \n'
}

oracle_keys=()
oracle_paths=()
if [[ "${ORACLE}" -eq 1 ]]; then
    for scene in ${SCENES}; do
        [[ -d "${DATA_DIR}/${scene}/queries" ]] || continue
        n_vf="$(phase_files "${scene}" vf)"
        n_ee="$(phase_files "${scene}" ee)"
        widest="${n_vf:-0}"
        [[ "${n_ee:-0}" -gt "${widest}" ]] && widest="${n_ee}"
        span="${ORACLE_CHUNK}"
        [[ "${span}" -le 0 || "${span}" -gt "${widest}" ]] && span="${widest}"
        [[ "${span}" -le 0 ]] && span=1
        for ((r = 1; r <= REPEATS; ++r)); do
            for ((fb = 0; fb < widest; fb += span)); do
                oracle_keys+=("${scene}|${r}|${fb}|${span}")
                # The span is in the name, not just the start: a file named
                # only by where it begins would be silently reused after a
                # change of --oracle-chunk, as covering a range it does not.
                oracle_paths+=("${OUT_DIR}/oracle/${scene}/r${r}-$(printf '%06d-%06d' "${fb}" "$((fb + span))").csv")
            done
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
    # Only files matching the chunk layout,
    # <scene>/<space>/<broadphase>/r<N>/<range>.csv.
    # Excluding a directory by name instead ("everything but oracle/") let a
    # sibling holding CSVs of another schema be swept in: 26 accuracy files
    # parked under the sweep tree merged into the timing CSV and added 250 rows
    # to it. Naming what belongs cannot go wrong that way.
    done < <(find "${OUT_DIR}" -mindepth 5 -maxdepth 5 \
                  -path "${OUT_DIR}/*/*/*/r*/*.csv" -print0 | sort -z)
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
        IFS='|' read -r scene space bp r begin end <<< "${chunk_keys[$i]}"
        printf '  TODO timing   %-20s %-6s %-6s r%-2s cases [%s, %s)\n' \
            "${scene}" "${space}" "${bp}" "${r}" "${begin}" "${end}"
    done
    for i in "${oracle_todo[@]:-}"; do
        [[ -n "${i}" ]] || continue
        IFS='|' read -r scene r fb span <<< "${oracle_keys[$i]}"
        printf '  TODO accuracy %-20s r%-2s files [%s, %s)\n' \
            "${scene}" "${r}" "${fb}" "$((fb + span))"
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
    local scene="$1" out="$2" file_begin="${3:-0}" span="${4:-0}"
    local tmp="${out}.partial"
    local range=()
    [[ "${file_begin}" -gt 0 ]] && range+=(--file-begin "${file_begin}")
    [[ "${span}" -gt 0 ]] && range+=(--max-files "${span}")
    mkdir -p "$(dirname "${out}")"
    # --no-strict so one scene's violation does not abort the sweep before the
    # rest is measured. The violations are in the CSV either way, and the report
    # fails on them; losing the remaining scenes as well helps nobody.
    "${TI_ORACLE}" "${DATA_DIR}/${scene}" --csv "${tmp}" --no-strict \
        "${range[@]}" >/dev/null 2>&1 || true
    [[ -s "${tmp}" ]] || return 1
    mv "${tmp}" "${out}"
}

# Empty means "every case", which is what the driver reads an unset variable as.
MAX_CASES_ENV=""
[[ "${MAX_CASES}" -gt 0 ]] && MAX_CASES_ENV="${MAX_CASES}"

run_chunk_body() {
    local scene="$1" space="$2" bp="$3" begin="$4" end="$5" out="$6"
    local tmp="${out}.partial"
    mkdir -p "$(dirname "${out}")"
    header > "${tmp}"
    for mode in ${MODES}; do
        SCCD_NARROWPHASE_MODE="${mode}" \
        SCCD_BENCH_EXECUTION_SPACE="${space}" \
        SCCD_BROADPHASE="${bp}" \
        SCCD_BENCH_CASE_BEGIN="${begin}" \
        SCCD_BENCH_CASE_END="${end}" \
        SCCD_BENCH_MAX_CASES="${MAX_CASES_ENV}" \
            "${SCCD_BENCH}" "${DATA_DIR}" "${scene}" | tail -n +2 >> "${tmp}"
    done

    # A chunk that produced no rows is a failure, not an empty result. The driver
    # writes its header before it touches a case, so a run where every case
    # failed still leaves a plausible-looking file -- and a header-only file is
    # non-empty, so it would be counted as done, never retried, and merged
    # silently into an otherwise healthy CSV. Refuse to publish it.
    if [[ "$(wc -l < "${tmp}")" -le 1 ]]; then
        printf 'error: %s produced no rows; leaving it unfinished\n' "${out}" >&2
        rm -f "${tmp}"
        return 1
    fi
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
pack_calls=()
pack_labels=()

flush_pack() {
    [[ "${#pack_calls[@]}" -gt 0 ]] || return 0
    local label="${pack_labels[0]}"
    [[ "${#pack_labels[@]}" -gt 1 ]] && \
        label="${pack_labels[0]} .. ${pack_labels[${#pack_labels[@]} - 1]} (${#pack_labels[@]})"
    printf '==> submitting %s\n' "${label}"
    local body
    body="$(printf '%s; ' "${pack_calls[@]}")"

    # With --jobs > 1 the blocking submit runs in a subshell so several packs
    # are in flight at once. Each records its own exit status in a file, because
    # a background job's status cannot be attributed to its label by `wait -n`
    # alone, and a silently dropped failure is the one outcome a sweep must not
    # have.
    if [[ "${JOBS}" -gt 1 ]]; then
        while [[ "$(jobs -rp | wc -l)" -ge "${JOBS}" ]]; do wait -n 2>/dev/null || true; done
        mkdir -p "${OUT_DIR}/logs"
        local status_file="${OUT_DIR}/logs/${pack_name}.status"
        (
            if submit_with_retry \
                    --account="${ACCOUNT}" --partition="${PARTITION}" \
                    --nodes=1 --ntasks=1 "${pack_gpu_args[@]}" \
                    --time="${TIME_LIMIT}" \
                    --uenv="${UENV}" --view=default \
                    --job-name="sccd-${pack_name}" \
                    --output="${OUT_DIR}/logs/${pack_name}.out" \
                    --error="${OUT_DIR}/logs/${pack_name}.err" \
                    --wrap="$(declare -f run_chunk_body); \
                             export SCCD_BENCH='${SCCD_BENCH}' DATA_DIR='${DATA_DIR}' MODES='${MODES}'; \
                             export SCCD_DB_TO_RAW='${SCCD_DB_TO_RAW}'; \
                             export MAX_CASES_ENV='${MAX_CASES_ENV}'; \
                             export OMP_NUM_THREADS=\"\$(nproc)\"; \
                             ${body}"; then
                printf 'ok %s\n' "${label}" > "${status_file}"
            else
                printf 'FAILED %s\n' "${label}" > "${status_file}"
            fi
        ) &
        pack_calls=(); pack_labels=()
        return 0
    fi

    if ! submit_with_retry \
            --account="${ACCOUNT}" --partition="${PARTITION}" \
            --nodes=1 --ntasks=1 "${pack_gpu_args[@]}" \
            --time="${TIME_LIMIT}" \
            --uenv="${UENV}" --view=default \
            --job-name="sccd-${pack_name}" \
            --output="${OUT_DIR}/logs/${pack_name}.out" \
            --error="${OUT_DIR}/logs/${pack_name}.err" \
            --wrap="$(declare -f run_chunk_body); \
                     export SCCD_BENCH='${SCCD_BENCH}' DATA_DIR='${DATA_DIR}' MODES='${MODES}'; \
                     export SCCD_DB_TO_RAW='${SCCD_DB_TO_RAW}'; \
                     export MAX_CASES_ENV='${MAX_CASES_ENV}'; \
                     export OMP_NUM_THREADS=\"\$(nproc)\"; \
                     header() { '${SCCD_BENCH}' --header; }; \
                     ${body}"; then
        printf 'FAILED %s (see %s/logs)\n' "${label}" "${OUT_DIR}" >&2
        failures=$((failures + ${#pack_calls[@]}))
    fi
    pack_calls=()
    pack_labels=()
}

for i in "${todo[@]:-}"; do
    [[ -n "${i}" ]] || continue
    IFS='|' read -r scene space bp r begin end <<< "${chunk_keys[$i]}"
    out="${chunk_paths[$i]}"
    label="${scene}/${space}/${bp}/r${r}/${begin}-${end}"

    if [[ "${LOCAL}" -eq 1 ]]; then
        printf '==> %s\n' "${label}"
        if ! run_chunk_body "${scene}" "${space}" "${bp}" "${begin}" "${end}" "${out}"; then
            printf 'FAILED %s\n' "${label}" >&2
            failures=$((failures + 1))
        fi
        continue
    fi

    pack_gpu_args=()
    [[ "${space}" == "device" ]] && pack_gpu_args=(--gpus-per-task=1)
    pack_name="timing-${space}-${bp}-r${r}-${begin}"
    pack_calls+=("run_chunk_body '${scene}' '${space}' '${bp}' '${begin}' '${end}' '${out}'")
    pack_labels+=("${label}")
    if [[ "${#pack_calls[@]}" -ge "${PACK}" ]]; then
        flush_pack
    fi
    continue

done
# Whatever did not fill a pack still has to be submitted.
flush_pack

oracle_calls=()
oracle_labels=()
oracle_name="oracle"

flush_oracle_pack() {
    [[ "${#oracle_calls[@]}" -gt 0 ]] || return 0
    local label="${oracle_labels[0]}"
    [[ "${#oracle_labels[@]}" -gt 1 ]] && \
        label="${oracle_labels[0]} .. ${oracle_labels[${#oracle_labels[@]} - 1]} (${#oracle_labels[@]})"
    printf '==> submitting %s\n' "${label}"
    local body
    body="$(printf '%s; ' "${oracle_calls[@]}")"
    if ! submit_with_retry \
            --account="${ACCOUNT}" --partition="${PARTITION}" \
            --nodes=1 --ntasks=1 --gpus-per-task=1 \
            --time="${TIME_LIMIT}" \
            --uenv="${UENV}" --view=default \
            --job-name="sccd-${oracle_name}" \
            --output="${OUT_DIR}/logs/${oracle_name}.out" \
            --error="${OUT_DIR}/logs/${oracle_name}.err" \
            --wrap="$(declare -f run_oracle_body); \
                     export TI_ORACLE='${TI_ORACLE}' DATA_DIR='${DATA_DIR}'; \
                     export OMP_NUM_THREADS=\"\$(nproc)\"; \
                     ${body}"; then
        printf 'FAILED %s (see %s/logs)\n' "${label}" "${OUT_DIR}" >&2
        failures=$((failures + ${#oracle_calls[@]}))
    fi
    oracle_calls=()
    oracle_labels=()
}

for i in "${oracle_todo[@]:-}"; do
    [[ -n "${i}" ]] || continue
    IFS='|' read -r scene r fb span <<< "${oracle_keys[$i]}"
    out="${oracle_paths[$i]}"
    label="accuracy ${scene}/r${r}/${fb}"

    if [[ "${LOCAL}" -eq 1 ]]; then
        printf '==> %s\n' "${label}"
        if ! run_oracle_body "${scene}" "${out}" "${fb}" "${span}"; then
            printf 'FAILED %s\n' "${label}" >&2
            failures=$((failures + 1))
        fi
        continue
    fi

    oracle_calls+=("run_oracle_body '${scene}' '${out}' '${fb}' '${span}'")
    oracle_labels+=("${label}")
    oracle_name="oracle-${scene}-r${r}-${fb}"
    if [[ "${#oracle_calls[@]}" -ge "${PACK}" ]]; then
        flush_oracle_pack
    fi
done
flush_oracle_pack

# Concurrent submits are still running at this point; the merge must not see a
# half-written tree, and a failure recorded in a background subshell has to be
# counted here or it is lost.
if [[ "${JOBS}" -gt 1 ]]; then
    wait
    while IFS= read -r sf; do
        [[ -n "${sf}" ]] || continue
        if [[ "$(head -c 6 "${sf}")" == FAILED ]]; then
            cat "${sf}" >&2
            failures=$((failures + 1))
        fi
        rm -f "${sf}"
    done < <(find "${OUT_DIR}/logs" -name '*.status' 2>/dev/null)
fi

printf '\n%d runs attempted, %d failed\n' \
    "$(( ${#todo[@]} + ${#oracle_todo[@]} ))" "${failures}"
merge
[[ "${failures}" -eq 0 ]]
