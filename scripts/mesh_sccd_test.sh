#!/usr/bin/env bash


set -e

export PATH=$INSTALL_DIR/smesh/bin:$PATH
DATA_PATH=/Users/patrickzulian/Desktop/code/installations/sources/Scalable-CCD/tests/data


# case=../data/n-body-simulation/frames/
# mesh_t0=balls16_1
# mesh_t1=balls16_2

# case=../data/armadillo-rollers/frames 
# mesh_t0="296"
# mesh_t1="297"

# case=../data/cloth-funnel/frames/
# mesh_t0=400
# mesh_t1=401

case=$DATA_PATH/armadillo-rollers/frames
mesh_t0=326
mesh_t1=327

# case=$DATA_PATH/cloth-ball/frames
# mesh_t0=cloth_ball92
# mesh_t1=cloth_ball93

# case=$DATA_PATH/cloth-funnel/frames/
# mesh_t0=227
# mesh_t1=228

# case=$DATA_PATH/n-body-simulation/frames/
# mesh_t0=balls16_18
# mesh_t1=balls16_19

# case=$DATA_PATH/rod-twist/frames/
# mesh_t0=3036
# mesh_t1=3037

set -x

db_to_raw $case/"$mesh_t0".ply "$mesh_t0"
db_to_raw $case/"$mesh_t1".ply "$mesh_t1"

# refine $mesh_t0 $mesh_t0
# refine $mesh_t1 $mesh_t1

# refine $mesh_t0 $mesh_t0 
# refine $mesh_t1 $mesh_t1 

# refine $mesh_t0 $mesh_t0 
# refine $mesh_t1 $mesh_t1 

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=10
# export SCCD_TOL=1e-16
# export SCCD_TOL=1e-6
# export SCCD_MAX_ITER=10
# export  SCCD_USE_TI=1

./mesh_sccd $mesh_t0 $mesh_t1
