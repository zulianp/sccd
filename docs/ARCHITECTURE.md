# Architecture

The surface that survived the cleanup, and what each part is for. The keep bar
was: a component ships if it satisfies the conservativeness invariant **and** it
is the only implementation of its job. Where several competed, the assessment in
`wip/ASSESSMENT.md` decided; everything demoted is in `spikes/` with the
reason recorded.

## The guarantee

One property runs through all of it. A reported time of impact is **conservative**:
at or before the true one, never after. A "no collision" answer means no
collision exists — misses are not traded against speed. Reporting early, or
reporting a contact that does not exist, costs work and never correctness, and
the whole search is built to round in that direction.

This is structural rather than a per-kernel promise. The narrow phase is branch
and bound over boxes in `(t, u, v)`, and for the multilinear function these
queries use, the hull of the eight corner values is the exact range over a box.
So accepting a box is always safe, refining is always safe, and the only way to
lose a root is an unsound rejection — discarding a box that contains one. That
is why every kernel pads its rejection test with a certified numerical error
bound, and why exhaustion (a depth cap, a full stack) **accepts** rather than
drops.

## Broad phase

Two implementations, both kept, because neither wins everywhere and they produce
**identical pair sets**.

| | `sccd_broadphase_sweep.hpp` | `sccd_broadphase_cell2d.hpp` |
|---|---|---|
| method | sweep and prune along a sorted axis | uniform 2D cell list, no sorting |
| wins | armadillo-rollers 1.59×, cloth-ball 1.36× | cloth-funnel; dense synthetic input by 4–7× |

Two dimensions, not three: a surface does not fill a volume, so a third axis buys
cells rather than selectivity.

`sccd_broadphase_strategy.hpp` chooses between them by **racing them** — one on a
step, the other on the next, keep the faster, re-probe periodically. Five
attempts to predict the winner from geometry all failed (anisotropy, expected
sweep window, estimated pair density, mesh size, and a fixed default), and they
are documented there so they are not retried. `SCCD_BROADPHASE=sweep|cell2d`
forces one.

Both handle triangles and quads, on host and device.

## Narrow phase

| mode | `NarrowPhaseMode` | kernel |
|---|---|---|
| 0 | `Fast` | scalar; codomain widths against domain tolerances |
| 2 | `Tight` | lane-packed; domain widths against domain tolerances |

Two modes, and both ship because neither dominates: `Fast` wins cloth-funnel on
speed, `Tight` wins cloth-ball, and `Tight` is 69× tighter at the median error.
They differ in accuracy and speed, never in safety — every mode is conservative.

Modes 1 and 3 used to exist and are gone. Mode 1 was the vectorised form of the
`Fast` predicate and was the **slowest** kernel in the library on every scene
measured — 15.9 ms against 6.4 and 8.1 on cloth-funnel, 90.9 against 19.4 and
17.5 on armadillo-rollers. It is in `spikes/`. Mode 3 ran mode 1 and corrected
each answer with TightInclusion; a hybrid of this library's kernel and the
reference is neither, so it is gone too. **To get TightInclusion's answer, call
TightInclusion**: `SCCD_USE_TI=1` dispatches straight to it, and the
`sccd_find_root_tight_inclusion_*` entry points are in the C ABI. Both need
`SCCD_ENABLE_TIGHT_INCLUSION=ON`.

Splitting is Gauss–Newton adaptive. Uniform splitting was a complete second
implementation, roughly 550 lines, that never won on any real scene; it is a
spike now.

**Quads** have their own path throughout, not a triangle path adapted: their own
inclusion function over the quad's own parameter domain, their own optimised host
root finder, and their own device kernel. `SCCD_NARROWPHASE_MODE` is still
ignored for quads — there is one root-finder variant, so the enum has nothing to
select between.

## CUDA

`src/cuda/` ships five kernels: broad phase (sweep), 2D cell list, narrow phase
(triangles), narrow phase (quads), and AABB construction.

The device broad phase is a clear win, 1.3× to 4.8× over 72 Grace cores — it is
the phase whose shape suits a GPU, being count, prefix sum and scatter with no
sequential window walk.

The device **narrow** phase depends entirely on which kernel you ask for, and an
earlier reading of the assessment got this wrong. Mode for mode against 72 Grace
cores, the default kernel is 2.0× *ahead* on cloth-ball and 2.2–3.6× behind on
the two smaller scenes; end to end, the whole pipeline on the GPU beats the whole
pipeline on the host by 3.2× on cloth-ball and is within 35% on the others. The
retracted claim — that it loses on every scene by up to 87× — came from a sweep
that set `SCCD_NARROWPHASE_MODE=2` on both sides, which is not the same kernel on
both sides: mode 2 is the host's *fastest* path and the device's *slowest* one.

What the numbers do show is a device-internal gap. Being conservative costs the
host 1.15× and the device 26×, on the same scene with the same tolerances. Part
of that is the price of the guarantee, which the host pays too; part of it is
that the device re-evaluates all eight corners of both children at every split
where the host inherits four of them from the parent. `wip/ASSESSMENT.md`
has the measurement and the decomposition.

Root finding computes in double regardless of the storage type: in single
precision the certified error bound and the tolerances that terminate the search
are too close together for the guarantee to survive. Float results are narrowed
toward negative infinity.

## Interfaces

- **`CCD<T>`** (`src/integrations/smesh/sccd_smesh_ccd.hpp`) — the main interface.
  `find_earliest_impact_time` and `find_impact_times` for one-shot use, or the
  staged broad/narrow calls to interleave your own logic. Needs smesh.
- **C ABI** (`src/api/sccd_c_api.cpp`, documented in `docs/API.md`) — seven exports for
  single queries, and the basis for `python/sccd.py`.
- **Headers** — sixteen installed, listed explicitly in `CMakeLists.txt`. A
  configure-time check fails on any header under `src/` that nobody classified,
  so the list cannot drift from the tree.

## Correctness gate

`ti_oracle` compares every mode against the datasets' exact roots and exits
non-zero on a missed collision or a late time of impact. It runs as part of
`ctest` when TightInclusion and the datasets are present.

`sccd_narrowphase_cuda_test` is the device narrow phase's gate. It builds queries
whose time of impact is known -- the primitive stays in a plane and the vertex
descends through it, so the contact time is one division, solved in `long double`
and rounded up -- then requires every kernel, host and device, triangle and quad,
to report a contact at or before it. It found that the device quad entry point
dereferenced device pointers on the host, which had never been called.

Twelve test executables. Test counts by configuration: 5 default, 7 with smesh,
9 with TightInclusion, 9 with CUDA. The default build needs no options and no
network.

## Building

```sh
cmake -S . -B build && cmake --build build      # zero dependencies
```

Options: `SCCD_ENABLE_SMESH`, `SCCD_ENABLE_CUDA`, `SCCD_ENABLE_TIGHT_INCLUSION`
(validation only), `SCCD_ENABLE_SPIKES` (demoted code, off).
`docs/ENVIRONMENT.md` lists every environment variable; none is needed to use the
library.

## Where the numbers are

This document describes the mechanism. Every measured figure lives in
[`BENCHMARKS.md`](BENCHMARKS.md); the record of which alternatives were tried,
kept, demoted or withdrawn lives in [`../wip/ASSESSMENT.md`](../wip/ASSESSMENT.md).
