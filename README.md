# SCCD — Simple Continuous Collision Detection

Continuous collision detection for moving triangle and quad meshes, with a
guarantee: a reported time of impact is **conservative** — at or before the true
one, never after — and "no collision" means no collision exists. Reporting early,
or reporting a contact that does not exist, costs work and never correctness.

## Build

```sh
cmake -S . -B build && cmake --build build -j
```

That is the whole thing. The default build has **no external dependencies** and
needs no network, no options and no environment variables. `CMAKE_BUILD_TYPE`
defaults to `Release` if you do not set it.

To run the tests:

```sh
ctest --test-dir build
```

To install:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build -j && cmake --install build
```

### Options

| Option | Default | What it does |
|---|---|---|
| `SCCD_ENABLE_NATIVE_ARCH` | `ON` | Build for the host CPU (`-march=native`). See the warning below. |
| `SCCD_ENABLE_OPENMP` | `ON` | Parallelise the host loops. |
| `SCCD_ENABLE_TBB` | `OFF` | Use TBB instead of OpenMP. |
| `SCCD_ENABLE_CUDA` | `OFF` | Build the device kernels. Needs a CUDA toolkit. |
| `SCCD_CUDA_ARCHITECTURES` | `90` | Target architecture. 90 is GH200; **set this for your GPU.** |
| `SCCD_ENABLE_SMESH` | `OFF` | Build the `CCD<T>` mesh interface and the demos. Needs [smesh](https://github.com/zulianp/smesh). |
| `SCCD_ENABLE_TIGHT_INCLUSION` | `OFF` | Validation only — builds the oracle that checks the kernels. Not a performance option. |
| `SCCD_ENABLE_SPIKES` | `OFF` | Build demoted code under `spikes/`. Not installed, not tested, deletable. |

**`SCCD_ENABLE_NATIVE_ARCH` matters more than it looks.** SCCD is mostly
header-only and its AABB overlap kernels are guarded on `__AVX2__`, `__AVX512F__`
and `__ARM_NEON`, so without a host-architecture flag the hottest broad-phase
kernel silently falls back to scalar code — in your translation units as well as
ours. Turn it off only for a portable binary, and then pass an explicit `-mavx2`
(or similar) in `CMAKE_CXX_FLAGS`.

## Using it

Three entry points, in decreasing order of how much they do for you.

**`CCD<T>`** — the mesh interface, and what most callers want. Needs
`SCCD_ENABLE_SMESH=ON`.

```cpp
#include "sccd_smesh_CCD.hpp"

sccd::CCD<double> ccd(mesh);
double toi;
ccd.find_earliest_impact_time(points_t0, points_t1, toi);   // one number for the step
ccd.find_impact_times(points_t0, points_t1, vf_tois, ee_tois);  // one per candidate
```

The broad and narrow phases are also callable separately if you want to
interleave your own logic between them. See
`src/integrations/smesh/sccd_smesh_CCD.hpp`.

**The C ABI** — seven functions, one query each, no broad phase and no batching.
Use it from C, or from another language; `python/sccd_py.py` is a `ctypes`
binding on top of it. Documented in [`docs/C_API.md`](docs/C_API.md).

**The kernels** — `sccd::narrow_phase_vf`, `narrow_phase_ee`, `narrow_phase_vq`
and the broad phases, if you have your own data layout. They take
structure-of-arrays geometry; the signatures are in `src/narrowphase.hpp`,
`src/narrowphase_vertex_quad.hpp` and `src/broadphase.hpp`.

### Quads

Quads are first class, not triangles with an extra node: their own inclusion
function, their own host root finder, their own device kernel.

**The four nodes are in lexicographic order, not cyclic order.** The bilinear
form weights node `k` by

```
w1 = (1-u)(1-v)    w2 = u(1-v)    w3 = (1-u)v    w4 = uv
```

so node 3 is the `(0,1)` corner and node 4 is `(1,1)`. Winding the nodes around
the quad swaps the last two and puts the surface somewhere the solver is not
looking, which shows up as contacts that are silently not found rather than as
an error.

## What is in the box

- **Broad phase** — a sweep-and-prune and a 2D cell list, both kept because
  neither wins everywhere; they produce identical pair sets.
  `broadphase_strategy.hpp` picks between them by racing them at run time.
  `SCCD_BROADPHASE=sweep|cell2d` forces one.
- **Narrow phase** — branch and bound over boxes in `(t, u, v)` with Gauss–Newton
  adaptive splitting. Two kernels ship: a scalar reference and a vectorised
  TightInclusion-exact one. Root finding computes in double regardless of the
  storage type, because in single precision the certified error bound and the
  tolerances that terminate the search are too close together for the guarantee
  to survive.
- **CUDA** — broad phase (sweep and cell list), narrow phase (triangles and
  quads), and AABB construction. The device broad phase beats 72 Grace cores by
  1.3–4.8×.

## Documentation

| | |
|---|---|
| [`docs/OVERVIEW.md`](docs/OVERVIEW.md) | What ships, what each part is for, and what was demoted with the reason. **Start here.** |
| [`docs/C_API.md`](docs/C_API.md) | The C ABI. |
| [`docs/ENVIRONMENT.md`](docs/ENVIRONMENT.md) | Every environment variable SCCD reads. None is needed to use the library. |
| [`benchmark/ASSESSMENT.md`](benchmark/ASSESSMENT.md) | The measurements behind the keep-or-demote decisions, including the ones that were retracted. |
| [`spikes/README.md`](spikes/README.md) | Demoted code, and why each piece is there. |

## A note on the numbers

Nothing in this repository claims a speedup that is not measured, and several
claims have been withdrawn after being measured properly — a broad-phase default
that a micro-benchmark supported and real scenes did not, and a host-versus-device
comparison that turned out to be running different kernels on the two sides. The
retractions are kept next to the claims rather than edited out.
