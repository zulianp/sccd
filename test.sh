#!/usr/bin/env bash

set -e

source data/venv/bin/activate 

cd build
make -j8

cd -
SCCD_LIB_PATH=build/libsccd.dylib python3 ccd_test.py