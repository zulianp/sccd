#!/usr/bin/env bash

set -e 

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR=$SCRIPTPATH/..

source $PROJECT_DIR/data/venv/bin/activate

set -x



python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/armadillo-rollers_ee  $PROJECT_DIR/csv/armadillo-rollers_ee_table.csv
python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/armadillo-rollers_fv  $PROJECT_DIR/csv/armadillo-rollers_vf_table.csv

# python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/n-body-simulation_ee  $PROJECT_DIR/csv/n-body-simulation_ee_table.csv
python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/n-body-simulation_vf  $PROJECT_DIR/csv/n-body-simulation_vf_table.csv

# python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/puffer-ball_ee 	   	 $PROJECT_DIR/csv/puffer-ball_ee_table.csv
python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/puffer-ball_vf 	     $PROJECT_DIR/csv/puffer-ball_vf_table.csv

# Bar charts: total false positives / false negatives per simulation scene (VF tables)
python3 $PROJECT_DIR/python/plot_ccd_fp_fn_bars.py $PROJECT_DIR/figures/vf_fp_fn_per_scene \
	armadillo-rollers $PROJECT_DIR/csv/armadillo-rollers_vf_table.csv \
	n-body-simulation $PROJECT_DIR/csv/n-body-simulation_vf_table.csv \
	puffer-ball $PROJECT_DIR/csv/puffer-ball_vf_table.csv

python3 $PROJECT_DIR/python/plot_ccd_fp_fn_bars.py $PROJECT_DIR/figures/ee_fp_fn_per_scene \
	armadillo-rollers $PROJECT_DIR/csv/armadillo-rollers_ee_table.csv \
	n-body-simulation $PROJECT_DIR/csv/n-body-simulation_ee_table.csv \
	puffer-ball $PROJECT_DIR/csv/puffer-ball_ee_table.csv
