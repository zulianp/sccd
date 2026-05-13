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
-DSCCD_ENABLE_TIGHT_INCLUSION=ON|OFF
-DSCCD_ENABLE_SMESH=ON|OFF
-DSCCD_ENABLE_CUDA=ON|OFF
```

## Approaches

SCCD is organized around a two-stage continuous collision detection pipeline for moving triangle meshes:

- **Broadphase candidate generation:** computes swept vertex, edge, and face AABBs, sorts them along the axis with the largest center variance, and scans sorted ranges to collect possible vertex-face and edge-edge overlaps. The CPU path uses parallel loops, prefix sums, cumulative maxima, and vectorized AABB rejection.
- **Narrowphase time-of-impact solving:** evaluates each candidate with vertex-face or edge-edge root finders over `(t, u, v)`, using uniform or adaptive domain splitting and optional refinement. When built with Tight Inclusion support, the narrowphase can dispatch to the Tight Inclusion root finder instead.
- **CUDA acceleration:** `src/cuda` mirrors the AABB setup, sorting, broadphase scans, lower-bound helpers, and narrowphase DFS kernels on the GPU, with CUB-based sorting/scans and warp/block-level reductions.
- **Integration and experiments:** the `smesh` integration wires the pipeline into mesh IO and host/device execution, while the Python scripts and data directories provide numerical root-finding checks, plotting, and benchmark inputs.
