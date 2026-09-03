# The C API

`src/sccd.cpp` is SCCD's only compiled translation unit. It exposes a small C ABI
so the root finders can be called from languages that are not C++;
`python/sccd_py.py` is the binding built on it, via `ctypes`.

It is a thin layer. Each entry point takes one query, runs the branch-and-bound
search on it, and writes back a time of impact and the parameter coordinates of
the contact. There is no broad phase here and no batching — for that, use the C++
`CCD<T>` interface in `src/integrations/smesh/sccd_smesh_CCD.hpp`.

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

## Two notes on the history of this file

Both are the kind of thing that only shows up once something is actually tested,
which this ABI was not until `src/tests/c_abi_test.exe.cpp` was written.

**It shipped the losing splitter.** Every C entry point defaulted
`SCCD_ADAPTIVE_SPLIT` to `0`, so the installed ABI ran uniform interval splitting
while the C++ path ran the adaptive splitter. The assessment measured uniform
behind adaptive on every real scene, and it has since been demoted to a spike.
The C entry points now use the adaptive splitter, like everything else.

**The bisecting variants could miss a collision.** At the subdivision cap, the
search discarded the box instead of accepting it — an unsound rejection, and the
one way this algorithm can lose a root. It was reachable: a vertex dropping
straight through a stationary triangle was missed by
`sccd_find_root_bisection_vf_d` while `sccd_find_root_bisection_vf_f` found it,
because in single precision the box met a tolerance condition before the cap and
in double it did not. Exhaustion now accepts at the box's `t` lower bound, as
every other termination path in that loop already did.

## Building against it

The library installs as `libsccd`, and the C declarations are in
`src/sccd.cpp`; there is no separate installed C header yet, so a consumer
declares the entry points it needs. `python/sccd_py.py` shows the shape, and
`python/ccd_test.py` exercises it.
