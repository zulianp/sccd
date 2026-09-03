# What SCCD ships

The surface that survived the cleanup, and what each part is for. The keep bar
was: a component ships if it satisfies the conservativeness invariant **and** it
is the only implementation of its job. Where several competed, the assessment in
`benchmark/ASSESSMENT.md` decided; everything demoted is in `spikes/` with the
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

| | `broadphase.hpp` | `cell2d_broadphase.hpp` |
|---|---|---|
| method | sweep and prune along a sorted axis | uniform 2D cell list, no sorting |
| wins | armadillo-rollers 1.59×, cloth-ball 1.36× | cloth-funnel; dense synthetic input by 4–7× |

Two dimensions, not three: a surface does not fill a volume, so a third axis buys
cells rather than selectivity.

`broadphase_strategy.hpp` chooses between them by **racing them** — one on a
step, the other on the next, keep the faster, re-probe periodically. Five
attempts to predict the winner from geometry all failed (anisotropy, expected
sweep window, estimated pair density, mesh size, and a fixed default), and they
are documented there so they are not retried. `SCCD_BROADPHASE=sweep|cell2d`
forces one.

Both handle triangles and quads, on host and device.

## Narrow phase

| mode | kernel | status |
|---|---|---|
| 0 | scalar reference | **ships** — wins cloth-funnel by 1.27× |
| 2 | TightInclusion-exact, vectorised | **ships** — wins cloth-ball by 1.62× |
| 1, 3 | fast-vector, TightInclusion-compat | validation only, require a TightInclusion build |

Modes 0 and 2 are not duplicates in the sense the keep bar means: each wins on
real scenes the other loses. Modes 1 and 3 both enter the same vectorised kernel;
mode 3 corrects its own answers with TightInclusion, which makes it an oracle
rather than a code path, and mode 1 is a duplicate that loses on every scene.
Asking for either without TightInclusion gets mode 0.

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
sequential window walk. The device **narrow** phase is not: it loses to the host
on every scene measured, by 21× and 87× on two of them. That is too large to be
a tuning gap, and the configuration worth supporting is likely device broad phase
feeding host narrow phase. It is measured in `benchmark/ASSESSMENT.md` rather
than presented as the device story.

Root finding computes in double regardless of the storage type: in single
precision the certified error bound and the tolerances that terminate the search
are too close together for the guarantee to survive. Float results are narrowed
toward negative infinity.

## Interfaces

- **`CCD<T>`** (`src/integrations/smesh/sccd_smesh_CCD.hpp`) — the main interface.
  `find_earliest_impact_time` and `find_impact_times` for one-shot use, or the
  staged broad/narrow calls to interleave your own logic. Needs smesh.
- **C ABI** (`src/sccd.cpp`, documented in `docs/C_API.md`) — seven exports for
  single queries, and the basis for `python/sccd_py.py`.
- **Headers** — sixteen installed, listed explicitly in `CMakeLists.txt`. A
  configure-time check fails on any header under `src/` that nobody classified,
  so the list cannot drift from the tree.

## Correctness gate

`ti_oracle` compares every mode against the datasets' exact roots and exits
non-zero on a missed collision or a late time of impact. It runs as part of
`ctest` when TightInclusion and the datasets are present.

Eleven test executables. Test counts by configuration: 5 default, 7 with smesh,
9 with TightInclusion, 8 with CUDA. The default build needs no options and no
network.

## Building

```sh
cmake -S . -B build && cmake --build build      # zero dependencies
```

Options: `SCCD_ENABLE_SMESH`, `SCCD_ENABLE_CUDA`, `SCCD_ENABLE_TIGHT_INCLUSION`
(validation only), `SCCD_ENABLE_SPIKES` (demoted code, off).
`docs/ENVIRONMENT.md` lists every environment variable; none is needed to use the
library.

## What was demoted, and where it went

`spikes/` — not built by default, not installed, not covered by the gate,
deletable without notice. Nothing shipped may include a spike header.

`hacks/` (an out-of-tree Scalable-CCD adapter that cannot configure, plus the
three headers only it used), a 3D cell list superseded by the 2D one, a
lower-bound broad phase with no caller, a warp broad phase called by nothing, an
all-to-all lower-bound kernel reachable only from its own demo, uniform interval
splitting, and an `external/json/` that was never in the build.
