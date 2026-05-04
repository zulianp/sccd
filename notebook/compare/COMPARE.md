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
