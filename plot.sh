#!/usr/bin/env bash

set -e 

python3 plot_ccd.py armadillo-rollers_ee  armadillo-rollers_ee_table.csv
python3 plot_ccd.py armadillo-rollers_fv  armadillo-rollers_vf_table.csv

python3 plot_ccd.py  n-body-simulation_ee n-body-simulation_ee_table.csv
python3 plot_ccd.py  n-body-simulation_vf n-body-simulation_vf_table.csv
python3 plot_ccd.py  puffer-ball_ee puffer-ball_ee_table.csv
python3 plot_ccd.py  puffer-ball_vf puffer-ball_vf_table.csv
