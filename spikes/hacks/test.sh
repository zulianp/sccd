#!/usr/bin/env bash

set -e

make -j8

echo "CASE: cloth-ball"
CASE=0 ./src/sccd/hacks/sccd_tests

echo "CASE: rod-twist"
CASE=1 ./src/sccd/hacks/sccd_tests

echo "CASE: armadillo"
CASE=2 ./src/sccd/hacks/sccd_tests

echo "CASE: cloth-funnel"
CASE=3 ./src/sccd/hacks/sccd_tests

echo "CASE: n-body-sim"
CASE=4 ./src/sccd/hacks/sccd_tests

echo "CASE: puffer-ball"
CASE=5 ./src/sccd/hacks/sccd_tests
