#!/usr/bin/env bash


set -e

export PATH=$INSTALL_DIR/smesh/bin:$PATH

####################
# DO not use to compare!!!
####################

# case=../data/n-body-simulation/frames/
# mesh_t0=balls16_1
# mesh_t1=balls16_2

# case=../data/armadillo-rollers/frames 
# mesh_t0="296"
# mesh_t1="297"

####################
# Use to compare!!!
####################

# case=/capstor/scratch/cscs/zulianp/Scalable-CCD/tests/data/armadillo-rollers/frames
# mesh_t0=326
# mesh_t1=327

case=/capstor/scratch/cscs/zulianp/Scalable-CCD/tests/data/cloth-ball/frames
mesh_t0=cloth_ball92
mesh_t1=cloth_ball93

# case=/capstor/scratch/cscs/zulianp/Scalable-CCD/tests/data/cloth-funnel/frames/
# mesh_t0=227
# mesh_t1=228

# case=/capstor/scratch/cscs/zulianp/Scalable-CCD/tests/data/n-body-simulation/frames/
# mesh_t0=balls16_18
# mesh_t1=balls16_19

# case=/capstor/scratch/cscs/zulianp/Scalable-CCD/tests/data/rod-twist/frames/
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

# refine $mesh_t0 $mesh_t0 
# refine $mesh_t1 $mesh_t1 

set +x

export OMP_NUM_THREADS=72
export OMP_PROC_BIND=true

ITERATIONS=4

echo "CUDA: "
rm -rf cuda_traces
mkdir -p cuda_traces
for ((i=1; i<=ITERATIONS; i++)); do
    echo "Iteration $i"
    $LAUNCH ./mesh_sccd_cuda $mesh_t0 $mesh_t1
    mv smesh.trace.csv cuda_traces/ccd_GPU_$i.csv
done

echo "CPU: "
rm -rf cpu_traces
mkdir -p cpu_traces
for ((i=1; i<=ITERATIONS; i++)); do
    echo "Iteration $i"
    SCCD_NEWTON_SPLIT_VF=1 SCCD_TOL=1e-12 SCCD_MAX_ITER=35 ./mesh_sccd      $mesh_t0 $mesh_t1
    mv smesh.trace.csv cpu_traces/ccd_CPU_$i.csv
done

tar -czvf cuda_traces.tar.gz cuda_traces
tar -czvf cpu_traces.tar.gz cpu_traces