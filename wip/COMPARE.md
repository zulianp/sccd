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
# Or
TBB_NUM_THREADS=72 ./tests/scalable_ccd_tests  "Test CUDA broad phase" -c "Cloth-Ball"        --log 0
```

The key broad-phase timing in `Scalable-CCD` is `sweep_and_tiniest_queue` (it excludes output buffer allocation timing, included in `sccd`):

```c
// Armadillo-Rollers
17.772 [ms]

// N-Body
56.204 [ms]

// Rod-Twist 
7.806 [ms]

// Cloth-Ball
38.910 [ms]
```


In `sccd`:

```c
// Armadillo-Rollers
6.904 [ms] 

// N-Body
31 [ms]

// Rod-Twist
3.391 [ms] 

// Cloth-ball
28.212 [ms]
```

## Narrow-Phase

```bash
TBB_NUM_THREADS=72 ./tests/scalable_ccd_tests  "Test CUDA narrow phase" --log 0
```

In single precision, the test currently fails due to `TOI=0`, but profiling still reports:

```c
// Float
FV: 142.942 [ms]
// Double
EE: 20.526 [ms]
FV: 38.930 [ms]
create_ccd_data: 35.338 [ms]
Total: 96.799 [ms], TOI: 3.814697265625e-06
```

For `cloth-ball`, steps 92-93:

In `sccd`:

```c
EE: 17.664     [ms]
FV: 8.029      [ms]
EE+FV: 25.693 [ms], TOI:  7.49206e-06
```











# SCCD

Timings are in seconds


```c
NP: 31.207 [ms]
EE: 21.329 [ms]
FV: 9.8733 [ms]
```

Double precision seems to be faster in the narrow phase


# Uniform vs Adaptive Earliest Time of Impact

## Uniform

```c
./mesh_sccd cloth_ball92 cloth_ball93
#faces 92230 #edges 138825 $nodes 46598, 0.16784 [s], toi 7.49198e-06

10.707 [ms]
11.418 [ms]
10.054 [ms]
 9.499 [ms]
```

## Adaptive

```c
./mesh_sccd cloth_ball92 cloth_ball93
#faces 92230 #edges 138825 $nodes 46598, 0.138413 [s], toi 7.49198e-06
10.313 [ms]
10.556 [ms]
10.122 [ms]
10.386 [ms]
```

# Uniform vs Adaptive All Times of Impact


## Uniform

```c
./mesh_sccd cloth_ball92 cloth_ball93
#faces 92230 #edges 138825 $nodes 46598, #e2e 5197332 #f2v 1655541, 0.737226 [s], 
toi 7.49198e-06, toi_vf 0.000147666, toi_ee 7.49198e-06

202.091 [ms]
122.469 [ms]
120.55 [ms]
200.227 [ms]
```

## Adaptive

```c
./mesh_sccd cloth_ball92 cloth_ball93
#faces 92230 #edges 138825 $nodes 46598, #e2e 5197332 #f2v 1655541, 0.693039 [s], 
toi 7.49198e-06, toi_vf 0.000147666, toi_ee 7.49198e-06

140.732 [ms]
137.916 [ms]
137.783 [ms]
137.233 [ms]
```