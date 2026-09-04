# API

Three entry points, in decreasing order of how much they do for you. All of them
report a time of impact that is **conservative**: at or before the true one,
never after.

## `CCD<T>` — the mesh interface

Needs `SCCD_ENABLE_SMESH=ON`. Declared in
`src/integrations/smesh/sccd_smesh_ccd.hpp`.

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

The phases are also callable separately, to interleave your own logic:

```cpp
ccd.broad_phase_prep(points_t0, points_t1);
ccd.broad_phase_fv_step(v_overlap, f_overlap);   // or broad_phase_ee_step
ccd.narrow_phase_fv(toi, vf_tois, max_depth, tol, /*toi_stride=*/0);
```

`narrow_phase_fv` and `narrow_phase_ee` take `max_toi` by reference: it bounds
the search on the way in and, for `toi_stride == 0`, comes back holding the
earliest time of impact.

## The kernels

If you have your own data layout, call the kernels directly. They take
structure-of-arrays geometry — `T** v0` is three row pointers, one per axis.

| function | header |
|---|---|
| `sccd::narrow_phase_vf<nxe>` | `sccd_narrowphase.hpp` |
| `sccd::narrow_phase_ee` | `sccd_narrowphase.hpp` |
| `sccd::narrow_phase_vq<nxe>` | `sccd_narrowphase_quad.hpp` |
| `sccd::device::narrow_phase_*` | `sccd_narrowphase.cuh`, `sccd_narrowphase_vq.cuh` |

Device entry points take **device** arrays of device pointers, not host arrays.

Which kernel runs is chosen by `SCCD_NARROWPHASE_MODE`; see
[`ENVIRONMENT.md`](ENVIRONMENT.md) and the accuracy table in
[`BENCHMARKS.md`](BENCHMARKS.md).

## The C ABI

`src/api/sccd_c_api.cpp` is SCCD's only compiled translation unit. It exposes a
small C ABI so the root finders can be called from languages that are not C++;
`python/sccd_py.py` is the binding built on it, via `ctypes`.

It is a thin layer. Each entry point takes one query, runs the branch-and-bound
search on it, and writes back a time of impact and the parameter coordinates of
the contact. There is no broad phase here and no batching — for that, use the C++
`CCD<T>` interface in `src/integrations/smesh/sccd_smesh_ccd.hpp`.

## Entry points

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

## What the result means

`*t` is a time of impact in `[0, 1]`, and it is **conservative**: it is at or
before the true one, never after. Reporting early costs a shorter step; reporting
late lets a simulation pass through the contact, which is why the whole search is
built to round the other way. A returned "no collision" means no collision
exists — misses are not a tolerance this library trades against speed.

`max_iter` caps the subdivision depth and `tol` sets the codomain tolerance.
Tightening either sharpens `*t` and costs time; neither can make the answer
unsafe, because the rejection test pads by a certified numerical error bound that
these do not touch.

## Building against it

The library installs as `libsccd`, and the C declarations are in
`src/api/sccd_c_api.cpp`; there is no separate installed C header yet, so a consumer
declares the entry points it needs. `python/sccd_py.py` shows the shape, and
`python/ccd_test.py` exercises it.
