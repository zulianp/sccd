#!/usr/bin/env bash

source $CODE_DIR/merge_git_repos/sfem/venv/bin/activate

set -e

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
PROJECT_DIR=$SCRIPTPATH/..

cd $PROJECT_DIR/build
make -j8

db_dir="$PROJECT_DIR/data"
datasets=(armadillo-rollers n-body-simulation puffer-ball)

cd -

for ds in ${datasets[@]}
do
	echo "---------------------------------------------------------"
	echo "Testing on $ds"
	SCCD_LIB_PATH=$PROJECT_DIR/build/libsccd.dylib python3 ccd_test.py $db_dir/$ds
	python3 ccd_test.py $db_dir/$ds ee
	python3 ccd_test.py $db_dir/$ds vf
	mv *.csv $PROJECT_DIR/csv
	echo "---------------------------------------------------------"
done
