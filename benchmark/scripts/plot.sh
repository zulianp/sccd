#!/usr/bin/env bash
#
# Render the figures from the cross-check tables produced by scripts/crosscheck.sh.
# Reads the tables crosscheck.sh wrote and renders the figures beside them.

set -e

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR="$( cd -- "${SCRIPTPATH}/../.." >/dev/null 2>&1 ; pwd -P )"

# shellcheck source=/dev/null
source "${SCRIPTPATH}/venv.sh"

OUT_DIR="${SCCD_BENCH_OUT_DIR:-${PROJECT_DIR}/benchmark/out}/crosscheck"
CSV_DIR="${OUT_DIR}"
FIG_DIR="${OUT_DIR}/figures"
mkdir -p "${FIG_DIR}"

set -x

python3 "${PROJECT_DIR}/python/sccd_plot_toi_error.py" "${FIG_DIR}/armadillo-rollers_ee"  "${CSV_DIR}/armadillo-rollers_ee_table.csv"
python3 "${PROJECT_DIR}/python/sccd_plot_toi_error.py" "${FIG_DIR}/armadillo-rollers_fv"  "${CSV_DIR}/armadillo-rollers_vf_table.csv"
python3 "${PROJECT_DIR}/python/sccd_plot_toi_error.py" "${FIG_DIR}/n-body-simulation_vf"  "${CSV_DIR}/n-body-simulation_vf_table.csv"
python3 "${PROJECT_DIR}/python/sccd_plot_toi_error.py" "${FIG_DIR}/puffer-ball_vf"        "${CSV_DIR}/puffer-ball_vf_table.csv"

# Bar charts: total false positives and false negatives per scene.
python3 "${PROJECT_DIR}/python/sccd_plot_fp_fn.py" "${FIG_DIR}/vf_fp_fn_per_scene" \
	armadillo-rollers "${CSV_DIR}/armadillo-rollers_vf_table.csv" \
	n-body-simulation "${CSV_DIR}/n-body-simulation_vf_table.csv" \
	puffer-ball       "${CSV_DIR}/puffer-ball_vf_table.csv"

python3 "${PROJECT_DIR}/python/sccd_plot_fp_fn.py" "${FIG_DIR}/ee_fp_fn_per_scene" \
	armadillo-rollers "${CSV_DIR}/armadillo-rollers_ee_table.csv" \
	n-body-simulation "${CSV_DIR}/n-body-simulation_ee_table.csv" \
	puffer-ball       "${CSV_DIR}/puffer-ball_ee_table.csv"
