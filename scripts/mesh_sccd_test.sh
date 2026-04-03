#!/usr/bin/env bash


set -e

export PATH=$INSTALL_DIR/smesh/bin:$PATH


case=../data/n-body-simulation/frames/
mesh_t0=balls16_1
mesh_t1=balls16_2

# case=../data/armadillo-rollers/frames 
# mesh_t0="296"
# mesh_t1="297"

set -x

db_to_raw $case/"$mesh_t0".ply "$mesh_t0"
db_to_raw $case/"$mesh_t1".ply "$mesh_t1"

refine $mesh_t0 $mesh_t0
refine $mesh_t1 $mesh_t1

refine $mesh_t0 $mesh_t0 
refine $mesh_t1 $mesh_t1 

# refine $mesh_t0 $mesh_t0 
# refine $mesh_t1 $mesh_t1 

./mesh_sccd $mesh_t0 $mesh_t1
