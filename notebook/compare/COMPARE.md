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
```

The key broad-phase timing in `Scalable-CCD` is `sweep_and_tiniest_queue` (it excludes output buffer allocation timing, included in `sccd`):

```json
// Armadillo-Rollers
"sweep_and_tiniest_queue": {
    "on_cpu": false,
    "on_gpu": true,
    "time_ms": 32.821696281433105
}
// N-Body
"sweep_and_tiniest_queue": {
    "on_cpu": false,
    "on_gpu": true,
    "time_ms": 109.64368057250977
}

// Rod-Twist
"sweep_and_tiniest_queue": {
    "on_cpu": false,
    "on_gpu": true,
    "time_ms": 14.969567775726318
}

```

In `sccd`:

```c
// Armadillo-Rollers
EE: 0.00227499  [s] 
FV: 0.000754833 [s]
Total 3.029823 [ms]  # ~10x speed-up

// N-Body
EE: 0.0443335 [s]
FV: 0.0248334 [s]
Total: 69.17 [ms]   # 1.5x speed-up

// Rod-Twist
E2E: 0.00593066
F2V: 0.00706148
Total: 12.99214 [ms] # 1.15x speed-up
```

## Narrow-Phase

```bash
TBB_NUM_THREADS=72 ./tests/scalable_ccd_tests  "Test CUDA narrow phase" --log 0
```

The test currently fails due to `TOI=0`, but profiling still reports:

```json
"FV CCD": {
    "on_cpu": false,
    "on_gpu": true,
    "time_ms": 142.94204711914063
}
```

For `cloth-ball`, steps 92-93:

In `sccd`:

```c
EE: 33.4055 [ms], TOI: 7.3699e-06
FV: TBA     [ms]
```