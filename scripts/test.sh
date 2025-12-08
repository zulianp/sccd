#!/usr/bin/env bash

set -e

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR=$SCRIPTPATH/..

source $PROJECT_DIR/data/venv/bin/activate

cd $PROJECT_DIR/build
make -j8

db_dir="$PROJECT_DIR/data"
datasets=(armadillo-rollers n-body-simulation puffer-ball)

cd -
export SCCD_LIB_PATH=$PROJECT_DIR/build/libsccd.dylib
for ds in ${datasets[@]}
do
	echo "---------------------------------------------------------"
	echo "Testing on $ds"
	python3 $PROJECT_DIR/python/ccd_test.py $db_dir/$ds ee
	python3 $PROJECT_DIR/python/ccd_test.py $db_dir/$ds vf
	mv *.csv $PROJECT_DIR/csv
	echo "---------------------------------------------------------"
done
