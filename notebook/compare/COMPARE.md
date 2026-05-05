<link rel="stylesheet" href="style.css">

# Compare Against Scalable-CCD

This page compares runtime measurements against `Scalable-CCD` using the same target scenarios.

## Build Scalable-CCD (Hopper)

```bash
git clone  https://github.com/Continuous-Collision-Detection/Scalable-CCD.git
cd Scalable-CCD
mkdir build 
cd build
cmake .. -DSCALABLE_CCD_BUILD_TESTS=ON  -DSCALABLE_CCD_WITH_PROFILER=ON -DSCALABLE_CCD_WITH_CUDA=ON  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_ARCHITECTURES=90
make -j
```

## Broad-Phase

Run a specific broad-phase scenario:

```bash
TBB_NUM_THREADS=72 ./tests/scalable_ccd_tests  "Test CUDA broad phase" -c "Armadillo-Rollers" --log 0
# Or
TBB_NUM_THREADS=72 ./tests/scalable_ccd_tests  "Test CUDA broad phase" -c "N-Body"            --log 0
# Or
TBB_NUM_THREADS=72 ./tests/scalable_ccd_tests  "Test CUDA broad phase" -c "Rod-Twist"         --log 0
```

The key broad-phase timing in `Scalable-CCD` is `sweep_and_tiniest_queue` (it excludes output buffer allocation timing, included in `sccd`):

```json
// Armadillo-Rollers
float : 32.821696281433105 [ms]
double : 15 [ms]

// N-Body
float: 109.64368057250977
double: 55.731937408447266

// Rod-Twist
float: 14.969567775726318
double: 8.026240110397339

```

In `sccd`:

```c
// Armadillo-Rollers
E2E: 0.00470495 [s]
F2V: 0.00219941 [s]
Total 6.90436 [ms] # 2.1x speed-up (faster)

// N-Body (1)
EE: 0.0443335 [s]
FV: 0.0248334 [s]
Total: 69.17 [ms]   # 0.79x speed-up (slower)

// N-Body (2)
E2E: 0.019861 [s]
F2V: 0.0115311 [s]
Total: 31.3921 [ms]  # 1.77 speed-up (faster)

// N-Body (3)
EE: 0.0209968 [s]
FV: 0.0100455 [s]
Total: 31 [ms]

// Rod-Twist (1)
E2E: 0.00593066
F2V: 0.00706148
Total: 12.99214 [ms] # 0.61x speed-up (slower)

// Rod-Twist (2)
EE: 0.00230765 [s]
FV: 0.00108385 [s]
Total: 3.3915 [ms] # 2.3  speed-up (faster)
```

## Narrow-Phase

```bash
TBB_NUM_THREADS=72 ./tests/scalable_ccd_tests  "Test CUDA narrow phase" --log 0
```

In single precision, the test currently fails due to `TOI=0`, but profiling still reports:

```c
// Float
FV: 142.94204711914063 [ms]
// Double
EE: 20.52681541442871 [ms]
FV: 38.93027114868164 [ms]
create_ccd_data: 35.338175773620605 [ms]
Total: 96.79945755004883 [ms], TOI: 3.814697265625e-06
```

For `cloth-ball`, steps 92-93:

In `sccd`:

```c
EE: 33.4055 [ms], TOI: 7.3699e-06
FV: TBA     [ms]
```

# Single vs Double (Rod-Twist)

Timings are in seconds

Single 

```csv
name,calls,total,avg
Broadphase,1,0.00670576,0.00670576
Broadphase: E2E,1,0.00222325,0.00222325
Broadphase: F2V,1,0.000902653,0.000902653
Mesh::node_to_node_graph_upper_triangular,1,0.01232,0.01232
Mesh::read,2,0.00465608,0.00232804
Narrow phase,1,0.0392153,0.0392153
Narrow phase: E2E,1,0.0392082,0.0392082
Narrow phase: F2V,1,1.43051e-06,1.43051e-06
Sorting AABBs,1,0.00126553,0.00126553
create_crs_graph_upper_triangular_from_element,1,0.0122893,0.0122893
create_n2e,1,0.00225544,0.00225544
mesh_sccd_cuda.exe,1,2.50109,2.50109
```

Double

```csv
name,calls,total,avg
Broadphase,1,0.203465,0.203465
Broadphase: E2E,1,0.00271726,0.00271726
Broadphase: F2V,1,0.00123119,0.00123119
Mesh::node_to_node_graph_upper_triangular,1,0.0125632,0.0125632
Mesh::read,2,0.00413275,0.00206637
Narrow phase,1,0.0265396,0.0265396
Narrow phase: E2E,1,0.0265284,0.0265284
Narrow phase: F2V,1,4.29153e-06,4.29153e-06
Sorting AABBs,1,0.00160575,0.00160575
create_crs_graph_upper_triangular_from_element,1,0.0125022,0.0125022
create_n2e,1,0.00219655,0.00219655
mesh_sccd_cuda.exe,1,2.54673,2.54673
```

Double precision seems to be faster in the narrow phase