#!/usr/bin/env bash

source $CODE_DIR/merge_git_repos/sfem/venv/bin/activate

set -e 

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR=$SCRIPTPATH/..

python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/armadillo-rollers_ee  $PROJECT_DIR/csv/armadillo-rollers_ee_table.csv
python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/armadillo-rollers_fv  $PROJECT_DIR/csv/armadillo-rollers_vf_table.csv

python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/n-body-simulation_ee  $PROJECT_DIR/csv/n-body-simulation_ee_table.csv
python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/n-body-simulation_vf  $PROJECT_DIR/csv/n-body-simulation_vf_table.csv
python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/puffer-ball_ee 	   	 $PROJECT_DIR/csv/puffer-ball_ee_table.csv
python3 $PROJECT_DIR/python/plot_ccd.py $PROJECT_DIR/figures/puffer-ball_vf 	     $PROJECT_DIR/csv/puffer-ball_vf_table.csv
