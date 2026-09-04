# SCCD

Continuous collision detection for moving triangle and quad meshes, on CPU and
GPU.

A reported time of impact is **conservative**: at or before the true one, never
after. "No collision" means no collision exists. Reporting early, or reporting a
contact that does not exist, costs work and never correctness.

Measured across twelve configurations and 101,164 exact roots: zero missed
collisions, zero late times of impact. See [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md).

## Build

```sh
cmake -S . -B build && cmake --build build -j
ctest --test-dir build
```

No external dependencies, no network, no options, no environment variables.
`CMAKE_BUILD_TYPE` defaults to `Release`.

| Option | Default | Effect |
|---|---|---|
| `SCCD_ENABLE_NATIVE_ARCH` | `ON` | `-march=native`. Off, the AVX2/AVX-512/NEON AABB kernels fall back to scalar — in your translation units as well as ours. |
| `SCCD_ENABLE_OPENMP` | `ON` | Parallelise host loops. |
| `SCCD_ENABLE_TBB` | `OFF` | TBB instead of OpenMP. |
| `SCCD_ENABLE_CUDA` | `OFF` | Device kernels. Set `SCCD_CUDA_ARCHITECTURES` for your GPU; it defaults to `90`. |
| `SCCD_ENABLE_SMESH` | `OFF` | The `CCD<T>` mesh interface and demos. Needs [smesh](https://github.com/zulianp/smesh). |
| `SCCD_ENABLE_TIGHT_INCLUSION` | `OFF` | Validation only — builds the accuracy oracle. Not a performance option. |
| `SCCD_ENABLE_SPIKES` | `OFF` | Demoted and dead code under `spikes/`. Not installed. |

## Use

SCCD is kernels over plain arrays — no container type, no mesh type, no required
dependency. Structure-of-arrays geometry in, times of impact out:

```cpp
#include "sccd_narrowphase.hpp"

std::vector<double> x0 = {0.0, 1.0, 0.0, 0.25};   // one array per axis
std::vector<double> y0 = {0.0, 0.0, 1.0, 0.25};
std::vector<double> z0 = {0.0, 0.0, 0.0, 1.00};
double* v0[3] = {x0.data(), y0.data(), z0.data()};
// ... v1 likewise for the end of the step

sccd::narrow_phase_vf<3, double, int>(n_pairs, vertex_of_pair, face_of_pair,
                                      v0, v1, /*face_stride=*/1, faces,
                                      /*max_toi=*/1.0, toi.data(),
                                      /*max_depth=*/69, /*tol=*/3e-8,
                                      /*toi_stride=*/1);
```

`demo/sccd_minimal.exe.cpp` is both stages end to end — two disconnected
triangles, swept AABBs, the sweep-and-prune broad phase, then the narrow phase —
and builds with no options:

```sh
cmake --build build -j --target sccd_minimal && ./build/sccd_minimal
```

It reports `0.499999999` for a contact whose exact time is `0.5`, and fails if
that is ever late. Early by 7e-10 is the safe direction, and the guarantee in one
line.

For meshes, `CCD<T>` wires the broad and narrow phases together over an smesh
mesh; that path is optional and needs `SCCD_ENABLE_SMESH=ON`. Both are in
[`docs/API.md`](docs/API.md).

### Quads

Quads are first class — their own inclusion function, host root finder and
device kernel.

**Nodes are in lexicographic order, not cyclic.** The bilinear form weights node
`k` by `(1-u)(1-v)`, `u(1-v)`, `(1-u)v`, `uv`, so the nodes are the corners
`(0,0)`, `(1,0)`, `(0,1)`, `(1,1)`. Winding them around the quad swaps the last
two and the solver searches a surface you did not mean — silently, since it still
returns a valid time of impact for that surface.

## Layout

```
src/core/          base, math, parallel, aabb
src/broadphase/    sweep-and-prune, 2D cell list, run-time strategy
src/narrowphase/   modes, root finders, tolerances, error bounds
src/cuda/          device kernels
src/api/           C ABI
python/            ctypes binding and the analysis tools
benchmark/         harness, drivers, and committed result tables
spikes/            demoted and dead code; not built, not installed
wip/               open work, decision records, retractions
```

## Documentation

| | |
|---|---|
| [`docs/BENCHMARKS.md`](docs/BENCHMARKS.md) | Full evaluation: timings, accuracy, occupancy, work counts. |
| [`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) | How the guarantee is obtained. |
| [`docs/API.md`](docs/API.md) | The kernels, the broad phase, the C ABI, Python, `CCD<T>`. |
| [`python/README.md`](python/README.md) | The `ctypes` binding and the analysis tools. |
| [`benchmark/README.md`](benchmark/README.md) | How the numbers are produced. |
| [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md) | Every environment variable. None is needed. |
| [`wip/`](wip/) | Open items and the measurement record, including withdrawn claims. |

## License

See [`LICENSE`](LICENSE).
