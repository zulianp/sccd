# SCCD — Simple Continuous Collision Detection

## Installation

SCCD uses CMake and requires a C++17 compiler.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build -j
cmake --install build
```

Useful optional CMake flags:

```sh
-DSCCD_ENABLE_TBB=ON|OFF
-DSCCD_ENABLE_OPENMP=ON|OFF
-DSCCD_ENABLE_NATIVE_ARCH=ON|OFF
-DSCCD_ENABLE_TIGHT_INCLUSION=ON|OFF
-DSCCD_USE_VNARROW_PHASE_DEFAULT=ON|OFF
-DSCCD_VNARROWPHASE_TI_COMPAT_DEFAULT=ON|OFF
-DSCCD_ENABLE_SMESH=ON|OFF
-DSCCD_ENABLE_CUDA=ON|OFF
```

The vectorized vertex-face narrow phase can also be selected at runtime:

```sh
SCCD_USE_VNARROW_PHASE=1
SCCD_VNARROWPHASE_TI_COMPAT=1
```

`SCCD_ENABLE_NATIVE_ARCH` (on by default) adds `-march=native`. SCCD is mostly
header-only and its AABB overlap kernels are guarded on `__AVX2__`,
`__AVX512F__` and `__ARM_NEON`, so without a host-architecture flag the hottest
broadphase kernel silently falls back to scalar code -- in your translation
units as well as ours. Turn it off only for a portable binary, and then pass an
explicit `-mavx2` (or similar) in `CMAKE_CXX_FLAGS`.

`SCCD_VNARROWPHASE_TI_COMPAT` requires a build with
`SCCD_ENABLE_TIGHT_INCLUSION=ON`. It runs the normal vector kernel and then
corrects its predicate and time-of-impact outputs with TightInclusion's DFS.
This strict comparison mode reproduces TightInclusion results; it is not the
performance path.

## Approaches

SCCD is organized around a two-stage continuous collision detection pipeline for moving triangle meshes:

- **Broadphase candidate generation:** computes swept vertex, edge, and face AABBs, sorts them along the axis with the largest center variance, and scans sorted ranges to collect possible vertex-face and edge-edge overlaps. The CPU path uses parallel loops, prefix sums, cumulative maxima, and vectorized AABB rejection.
- **Narrowphase time-of-impact solving:** evaluates each candidate with vertex-face or edge-edge root finders over `(t, u, v)`, using uniform or adaptive domain splitting and optional refinement. When built with Tight Inclusion support, the narrowphase can dispatch to the Tight Inclusion root finder instead.
- **CUDA acceleration:** `src/cuda` mirrors the AABB setup, sorting, broadphase scans, lower-bound helpers, and narrowphase DFS kernels on the GPU, with CUB-based sorting/scans and warp/block-level reductions.
- **Integration and experiments:** the `smesh` integration wires the pipeline into mesh IO and host/device execution, while the Python scripts and data directories provide numerical root-finding checks, plotting, and benchmark inputs.
