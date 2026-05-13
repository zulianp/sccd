
# The goal of this file is to perform a wide benchmark of the SCCD library using the NYU datasets
# We start with benchmarking the broad-phase collision detection

# TODO: Create a separate bash script to download the datasets (currently download 2 of them for testing and enable each dataset with an env var. e.g., SCCD_ENABLE_ARMADILLO_ROLLERS=1)
# https://archive.nyu.edu/handle/2451/74508

# TODO: 
# 1) create parser using https://github.com/simdjson/simdjson.git (separate executable and compilation project) in external/json
# The goal is to parse the boxes json files e.g.,  data/armadillo-rollers/boxes/0ee.json)
# and generate to binary (raw) collision files c0.int32 c1.int32 (SoA) put in the correct folder (e.g., boxes/0ee/)=

# TODO:
# 2) Create a bench.exe.cpp that reads the meshes and scans the folder boxes and reads the raw files, times the CCD for each file collision files pair (names the trace file after the case and folder e.g., SMESH_TRACE_FILE=armadillo-rollers/0ee)
# the timings are collected as milliseconds and stored in a unique raw binary file for the whole case and collision type e.g., armadillo-rollers-fv.float64 and armadillo-rollers-ee.float64 



