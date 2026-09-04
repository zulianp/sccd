# API

SCCD is a library of kernels over plain arrays. It has no container type, no mesh
type and no required dependency — you hand it structure-of-arrays geometry and it
hands back times of impact. Everything below is that, in increasing order of how
much is done for you.

Every entry point reports a time of impact that is **conservative**: at or before
the true one, never after. A returned "no collision" means none exists.

## The kernels

The narrow phase is the library. One call, three arrays of geometry per time
step, and one time of impact per candidate pair.

```cpp
template <int nxe, typename T, typename I>
int sccd::narrow_phase_vf(size_t n_pairs,
                          const I* vertex_of_pair, const I* face_of_pair,
                          T** v0, T** v1,          // start / end geometry, SoA
                          size_t face_stride, I** faces,
                          T max_toi, T* toi,       // out
                          int max_depth, T tol,
                          int toi_stride = 0);
```

`T** v0` is three row pointers — x, y and z — not a matrix type. `I** faces` is
one row per vertex slot of the element, indexed by face. That is the whole data
contract, and `std::vector` satisfies it:

```cpp
std::vector<double> x0 = {0.0, 1.0, 0.0, 0.25};
std::vector<double> y0 = {0.0, 0.0, 1.0, 0.25};
std::vector<double> z0 = {0.0, 0.0, 0.0, 1.00};
double* v0[3] = {x0.data(), y0.data(), z0.data()};
```

`toi_stride` selects what comes back: `1` writes one time of impact per pair,
`0` writes a single earliest time to `toi[0]` and lets every query prune against
the running minimum, which is markedly cheaper.

| call | header | query |
|---|---|---|
| `sccd::narrow_phase_vf<3>` | `sccd_narrowphase.hpp` | vertex against triangle |
| `sccd::narrow_phase_ee` | `sccd_narrowphase.hpp` | edge against edge |
| `sccd::narrow_phase_vq<4>` | `sccd_narrowphase_quad.hpp` | vertex against quad |
| `sccd::device::narrow_phase_*` | `sccd_narrowphase.cuh`, `sccd_narrowphase_vq.cuh` | the same, on the GPU |

Device entry points take **device** arrays of device pointers, not host arrays.

Which kernel runs is chosen by `SCCD_NARROWPHASE_MODE`; see
[`ENVIRONMENT.md`](ENVIRONMENT.md) for the modes and
[`BENCHMARKS.md`](BENCHMARKS.md) for what each costs and how tight it is.

### A complete working example

[`demo/sccd_minimal.exe.cpp`](../demo/sccd_minimal.exe.cpp) is the above end to
end, with `std::vector` and no other dependency. It builds in the default
configuration and runs as part of `ctest`:

```sh
cmake -S . -B build && cmake --build build -j --target sccd_minimal
./build/sccd_minimal
```

```
query  vertex  face  time of impact
  0       3      0    0.249999998
  1       4      0    no collision

exact time of impact for query 0 is 0.250000000
reported 0.249999998, which is at or before it, as guaranteed (early by 1.966e-09)
```

The scene is built so the exact answer is known: a vertex crossing the triangle's
plane at `t = 1/4`, and a second vertex that never touches it. The reported value
is 2e-9 early — the safe direction, and the whole point of the guarantee.

## The broad phase

To go from a mesh to candidate pairs, compute swept AABBs and run one of the two
broad phases. Both produce identical pair sets; `sccd_broadphase_strategy.hpp`
picks between them at run time, and `SCCD_BROADPHASE=sweep|cell2d` forces one.

| | header |
|---|---|
| swept AABBs | `sccd_aabb.hpp` (`sccd::compute_aabbs`) |
| sweep and prune | `sccd_broadphase_sweep.hpp` |
| 2D cell list | `sccd_broadphase_cell2d.hpp` |
| run-time choice | `sccd_broadphase_strategy.hpp` |

## The C ABI

For callers that are not C++. `src/api/sccd_c_api.cpp` is SCCD's only compiled
translation unit; `python/sccd_py.py` is a `ctypes` binding on top of it.

It is a thin layer. Each entry point takes one query, runs the branch-and-bound
search on it, and writes back a time of impact and the parameter coordinates of
the contact. There is no broad phase here and no batching.

### Entry points

All of them share a shape:

```c
int sccd_find_root_<kind>(int max_iter, <T> tol,
                          const <T> a0[3], const <T> a1[3], const <T> a2[3], const <T> a3[3],
                          const <T> b0[3], const <T> b1[3], const <T> b2[3], const <T> b3[3],
                          <T>* t, <T>* u, <T>* v);
```

The first four points are the configuration at the start of the step and the
next four the configuration at the end. The return value is non-zero when a
collision was found.

| Function | Query | Scalar |
|---|---|---|
| `sccd_find_root_vf_f` | vertex–face | `float` |
| `sccd_find_root_vf_d` | vertex–face | `double` |
| `sccd_find_root_ee_d` | edge–edge | `double` |
| `sccd_find_root_bisection_vf_f` | vertex–face, bisecting splitter | `float` |
| `sccd_find_root_bisection_vf_d` | vertex–face, bisecting splitter | `double` |
| `sccd_find_root_tight_inclusion_vf_d` | vertex–face, TightInclusion | `double` |
| `sccd_find_root_tight_inclusion_ee_d` | edge–edge, TightInclusion | `double` |

The last two exist only when the build has `SCCD_ENABLE_TIGHT_INCLUSION`. They
are for validation, not for use.

For vertex–face the points are the vertex and the face's three corners, and
`(u, v)` are barycentric on the face. For edge–edge they are the two endpoints of
each edge, and `u` and `v` parameterise the first and second edge.

### What the result means

`*t` is a time of impact in `[0, 1]`, and it is **conservative**: it is at or
before the true one, never after. Reporting early costs a shorter step; reporting
late lets a simulation pass through the contact, which is why the whole search is
built to round the other way. A returned "no collision" means no collision
exists — misses are not a tolerance this library trades against speed.

`max_iter` caps the subdivision depth and `tol` sets the codomain tolerance.
Tightening either sharpens `*t` and costs time; neither can make the answer
unsafe, because the rejection test pads by a certified numerical error bound that
these do not touch.


## With smesh

`CCD<T>` wires the phases together over an smesh mesh, so you hand it point
buffers rather than assembling the SoA arrays yourself. It is **optional** —
build with `SCCD_ENABLE_SMESH=ON`.

```cpp
#include "sccd_smesh_ccd.hpp"

sccd::CCD<double> ccd(mesh);

double toi = 1.0;
ccd.find_earliest_impact_time(points_t0, points_t1, toi);       // one value per step
ccd.find_impact_times(points_t0, points_t1, vf_tois, ee_tois);  // one per candidate
```

| call | writes | cost |
|---|---|---|
| `find_earliest_impact_time` | one time of impact for the whole step | every query prunes against the running minimum |
| `find_impact_times` | one per candidate pair | no shared bound; 1.6–4.4× the above on the host, more on the device |

The phases are callable separately, to interleave your own logic:

```cpp
ccd.broad_phase_prep(points_t0, points_t1);
ccd.broad_phase_fv_step(v_overlap, f_overlap);   // or broad_phase_ee_step
ccd.narrow_phase_fv(toi, vf_tois, max_depth, tol, /*toi_stride=*/0);
```

`narrow_phase_fv` and `narrow_phase_ee` take `max_toi` by reference: it bounds
the search on the way in and, for `toi_stride == 0`, comes back holding the
earliest time of impact.

One caveat worth knowing: **smesh stores coordinates as `float`**
(`typedef float geom_t`). The root finders still compute in double, but the
geometry they are handed has been rounded on the way into the mesh. If you are
comparing against an exact reference, compare against roots computed for the
coordinates the mesh actually holds.

## Building against it

The library installs as `libsccd`, and the C declarations are in
`src/api/sccd_c_api.cpp`; there is no separate installed C header yet, so a consumer
declares the entry points it needs. `python/sccd_py.py` shows the shape, and
`python/ccd_test.py` exercises it.
