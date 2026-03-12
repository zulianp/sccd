#!/usr/bin/env bash


set -e

# case=../data/n-body-simulation/frames/
# mesh_t0=balls16_1
# mesh_t1=balls16_2

case=../data/armadillo-rollers/frames 
mesh_t0="296"
mesh_t1="297"

# ls -1 $case

set -x

$CODE_DIR/smesh/python/smesh/db_to_raw.py $case/"$mesh_t0".ply "$mesh_t0"
$CODE_DIR/smesh/python/smesh/db_to_raw.py $case/"$mesh_t1".ply "$mesh_t1"

$CODE_DIR/smesh/python/smesh/raw_to_db.py "$mesh_t0" mesh.00.vtk 
$CODE_DIR/smesh/python/smesh/raw_to_db.py "$mesh_t1" mesh.01.vtk 

./mesh_sccd $mesh_t0 $mesh_t1
