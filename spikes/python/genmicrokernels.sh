#!/usr/bin/env bash
#
# Run the SymPy generators and write their output to a directory.
#
#     ./genmicrokernels.sh <output-dir>
#
# This used to write straight into src/. It does not any more, and it refuses
# to: the headers it once produced have been corrected by hand since, and
# regenerating them would revert the min(max_abs, 1) clamp in the numerical
# error bound and resurrect eight objective functions that were demoted to
# spikes/src/dead.hpp. See README.md in this directory for the full divergence.
#
# What comes out is a derivation to read against the shipped header, not a
# replacement for it.

set -euo pipefail

SCRIPTPATH="$(cd -- "$(dirname "$0")" >/dev/null 2>&1 && pwd -P)"
PROJECT_DIR="$(cd -- "${SCRIPTPATH}/../.." >/dev/null 2>&1 && pwd -P)"

if [ $# -ne 1 ]; then
    echo "usage: $0 <output-dir>" >&2
    exit 1
fi

OUT_DIR="$1"
mkdir -p "${OUT_DIR}"
OUT_DIR="$(cd "${OUT_DIR}" && pwd -P)"

case "${OUT_DIR}/" in
    "${PROJECT_DIR}/src/"*)
        echo "error: refusing to write into src/. The headers there are hand-corrected;" >&2
        echo "       see $(dirname "$0")/README.md." >&2
        exit 1
        ;;
esac

PYTHON="${PYTHON:-python3}"

"${PYTHON}" "${SCRIPTPATH}/codegen_numerical_error.py" > "${OUT_DIR}/numerical_error.inc"
"${PYTHON}" "${SCRIPTPATH}/codegen_tolerance.py"       > "${OUT_DIR}/tolerance.inc"
"${PYTHON}" "${SCRIPTPATH}/codegen_objective.py"       > "${OUT_DIR}/objective.inc"
"${PYTHON}" "${SCRIPTPATH}/codegen_interval.py"        > "${OUT_DIR}/interval.inc"
"${PYTHON}" "${SCRIPTPATH}/codegen_bisect.py"          > "${OUT_DIR}/bisect.inc"

echo "wrote 5 files to ${OUT_DIR}"
