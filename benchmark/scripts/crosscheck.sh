#!/usr/bin/env bash
#
# Run python/sccd_crosscheck.py over the datasets: every query through the C ABI,
# compared against the datasets' exact roots. Fails on a missed collision or a
# time of impact after the exact root; false positives are reported and allowed.
#
# This is not the test suite -- that is `ctest --test-dir build`, which needs no
# data. This needs the benchmark datasets under data/, which are gigabytes.

set -e

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR="$( cd -- "${SCRIPTPATH}/../.." >/dev/null 2>&1 ; pwd -P )"

# shellcheck source=/dev/null
source "${SCRIPTPATH}/venv.sh"

BUILD_DIR="${BUILD_DIR:-${PROJECT_DIR}/build}"
cmake --build "${BUILD_DIR}" -j

DB_DIR="${PROJECT_DIR}/data"
CSV_DIR="${PROJECT_DIR}/research/csv"
mkdir -p "${CSV_DIR}"

datasets=(armadillo-rollers n-body-simulation puffer-ball)

export SCCD_TOL=1e-16
export SCCD_LIB_PATH="${BUILD_DIR}"

for ds in "${datasets[@]}"; do
	echo "---------------------------------------------------------"
	echo "Cross-checking $ds"
	python3 "${PROJECT_DIR}/python/sccd_crosscheck.py" "${DB_DIR}/${ds}" ee
	python3 "${PROJECT_DIR}/python/sccd_crosscheck.py" "${DB_DIR}/${ds}" vf
	mv ./*_table.csv "${CSV_DIR}/"
	echo "---------------------------------------------------------"
done
