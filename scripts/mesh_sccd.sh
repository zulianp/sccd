#!/usr/bin/env bash
#
# Run the mesh demo on two frames of a benchmark scene, on the CPU or the GPU.
#
#     scripts/mesh_sccd.sh <scene> <frame0> <frame1> [options]
#
#     scripts/mesh_sccd.sh armadillo-rollers 326 327
#     scripts/mesh_sccd.sh cloth-ball cloth_ball92 cloth_ball93 --cuda --repeats 4
#
# Options:
#   --cuda           run mesh_sccd_cuda instead of mesh_sccd
#   --binary NAME    run some other demo taking <mesh_t0> <mesh_t1>, such as the
#                    spike mesh_lower_bound_cuda
#   --repeats N      run N times, keeping each run's smesh.trace.csv (default 1)
#   --build DIR      build directory (default: build, or $BUILD_DIR)
#   --data DIR       dataset directory (default: data, or $SCCD_DATA_DIR)
#   --out DIR        where to put the traces (default: traces-<cpu|cuda>)
#
# The demo needs a build with -DSCCD_ENABLE_SMESH=ON, and --cuda additionally
# needs -DSCCD_ENABLE_CUDA=ON.
#
# $LAUNCH is prefixed to the binary, for running under a job launcher:
#
#     LAUNCH="srun --account=c40 --partition=debug --nodes=1 --ntasks=1 \
#             --gpus-per-task=1 --time=00:10:00" \
#         scripts/mesh_sccd.sh cloth-ball cloth_ball92 cloth_ball93 --cuda
#
# Frames come from this repository's own data/, which
# benchmark/scripts/download_datasets.sh populates.

set -euo pipefail

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR="$( cd -- "${SCRIPTPATH}/.." >/dev/null 2>&1 ; pwd -P )"

# The comment block at the top of this file is the usage message; print it
# rather than keeping a second copy that can drift from it.
usage() {
    awk 'NR > 2 && /^#/ { sub(/^# ?/, ""); print; next } NR > 2 { exit }' "$0"
    exit "${1:-1}"
}

case "${1:-}" in -h|--help|"") usage 0 ;; esac
[ $# -ge 3 ] || { echo "error: need <scene> <frame0> <frame1>" >&2; usage 1; }

scene="$1"; frame0="$2"; frame1="$3"; shift 3

use_cuda=0
binary_override=""
repeats=1
build_dir="${BUILD_DIR:-${PROJECT_DIR}/build}"
data_dir="${SCCD_DATA_DIR:-${PROJECT_DIR}/data}"
out_dir=""

while [ $# -gt 0 ]; do
    case "$1" in
        --cuda)    use_cuda=1; shift ;;
        --binary)  binary_override="$2"; shift 2 ;;
        --repeats) repeats="$2"; shift 2 ;;
        --build)   build_dir="$2"; shift 2 ;;
        --data)    data_dir="$2"; shift 2 ;;
        --out)     out_dir="$2"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) echo "error: unknown option $1" >&2; usage 1 ;;
    esac
done

if [ -n "${binary_override}" ]; then
    binary_name="${binary_override}"
    # An explicit binary carries its own answer to "is this the GPU one", which
    # only decides the default thread count and the trace directory's name.
    case "${binary_name}" in *cuda*) use_cuda=1 ;; esac
else
    binary_name=$([ "${use_cuda}" -eq 1 ] && echo mesh_sccd_cuda || echo mesh_sccd)
fi
binary="${build_dir}/${binary_name}"
: "${out_dir:=${PWD}/traces-$([ "${use_cuda}" -eq 1 ] && echo cuda || echo cpu)}"

if [ ! -x "${binary}" ]; then
    echo "error: ${binary} is not built." >&2
    echo "       cmake -S ${PROJECT_DIR} -B ${build_dir} -DSCCD_ENABLE_SMESH=ON$([ "${use_cuda}" -eq 1 ] && echo ' -DSCCD_ENABLE_CUDA=ON')" >&2
    echo "       cmake --build ${build_dir} -j --target ${binary_name}" >&2
    exit 1
fi

frames="${data_dir}/${scene}/frames"
if [ ! -d "${frames}" ]; then
    echo "error: no frames at ${frames}." >&2
    echo "       ${PROJECT_DIR}/benchmark/scripts/download_datasets.sh fetches them." >&2
    exit 1
fi

# smesh reads a directory of raw arrays, not a PLY file. benchmark/ply_to_smesh.py
# does that conversion with nothing but the standard library; smesh's own
# db_to_raw is used instead when it is on PATH, since it is the reference.
convert() {
    local ply="$1" out="$2"
    [ -f "${ply}" ] || { echo "error: no such frame: ${ply}" >&2; exit 1; }
    if command -v db_to_raw >/dev/null 2>&1; then
        db_to_raw "${ply}" "${out}"
    else
        python3 "${PROJECT_DIR}/benchmark/ply_to_smesh.py" "${ply}" "${out}"
    fi
}

work_dir="$(mktemp -d "${TMPDIR:-/tmp}/mesh_sccd.XXXXXX")"
trap 'rm -rf "${work_dir}"' EXIT

convert "${frames}/${frame0}.ply" "${work_dir}/${frame0}"
convert "${frames}/${frame1}.ply" "${work_dir}/${frame1}"

export OMP_PROC_BIND="${OMP_PROC_BIND:-true}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$( [ "${use_cuda}" -eq 1 ] && echo 72 || echo 10 )}"

mkdir -p "${out_dir}"
for ((i = 1; i <= repeats; i++)); do
    echo "==> ${binary_name} ${scene} ${frame0} ${frame1}  (run ${i}/${repeats})"
    ( set -x; ${LAUNCH:-} "${binary}" "${work_dir}/${frame0}" "${work_dir}/${frame1}" )
    if [ -f smesh.trace.csv ]; then
        mv smesh.trace.csv "${out_dir}/${binary_name}_${i}.csv"
    fi
done

echo "traces in ${out_dir}"
