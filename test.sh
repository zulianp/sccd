#!/usr/bin/env bash

set -e

source data/venv/bin/activate 

cd build
make -j8

db_dir=/Users/patrickzulian/Desktop/code/installations/sources/Scalable-CCD/tests/data
datasets=(armadillo-rollers  cloth-ball cloth-funnel n-body-simulation puffer-ball rod-twist)

cd -

for ds in ${datasets[@]}
do
	echo "Testing on $ds"
	SCCD_LIB_PATH=build/libsccd.dylib python3 ccd_test.py $db_dir/$ds
done