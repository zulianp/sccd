#!/usr/bin/env bash

set -e

source data/venv/bin/activate 

# cd build
cd /Users/patrickzulian/Desktop/code/installations/sources/Scalable-CCD/build
make -j8

# db_dir=/Users/patrickzulian/Desktop/code/installations/sources/Scalable-CCD/tests/data
# datasets=(armadillo-rollers cloth-ball cloth-funnel puffer-ball rod-twist n-body-simulation)

db_dir="data"
# datasets=(armadillo-rollers)
datasets=(armadillo-rollers n-body-simulation puffer-ball)

cd -

for ds in ${datasets[@]}
do
	echo "---------------------------------------------------------"
	echo "Testing on $ds"
	# SCCD_LIB_PATH=build/libsccd.dylib python3 ccd_test.py $db_dir/$ds
	export SCCD_LIB_PATH=/Users/patrickzulian/Desktop/code/installations/sources/Scalable-CCD/build/src/sccd/libsccd.dylib
	# python3 ccd_test.py $db_dir/$ds ee
	python3 ccd_test.py $db_dir/$ds vf
	echo "---------------------------------------------------------"
done