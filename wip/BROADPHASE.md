# Broad phase: measurements and where the time goes

Measured with `sccd_refine_scaling` on cloth-ball frames 40 → 41, host, refined
to 1.48M triangles. See `REFINE_SCALING.md` for the harness.

## The broad phase is the collision step

| level | faces | vf_pairs | ee_pairs | prep_ms | bp_fv_ms | bp_ee_ms | broad_ms | narrow_ms |
|---|---|---|---|---|---|---|---|---|
| 0 | 92,230 | 12,104 | 56,881 | 33 | 9 | 22 | 64 | 2.5 |
| 1 | 368,920 | 122,626 | 505,058 | 59 | 127 | 250 | 436 | 20 |
| 2 | 1,475,680 | 2,080,634 | 6,481,554 | 237 | 1,755 | 2,312 | 4,305 | 182 |

At 1.5M triangles the broad phase is **96%** of the step, and the two sweeps are
**94%** of the broad phase. Prep -- the AABBs and the sort -- scales cleanly with
element count (1.8x then 4.0x per 4x refinement), so it is not what makes the
curve bend. It is still 237 ms of `O(n log n)` sorting that a cell grid would not
need at all, which is a second reason to want one.

Per pair, the sweep costs far more than resolving that pair costs:

| level | ns to *find* a vf pair | ns to *find* an ee pair | ns to *resolve* a pair |
|---|---|---|---|
| 0 | 729 | 383 | 36 |
| 1 | 1,039 | 495 | 32 |
| 2 | 843 | 357 | 21 |

**Finding a candidate costs about 22x what resolving it costs.** Any further
narrow-phase optimisation is worth at most 4% of the step on a scene like this.

## Why: the sweep scans ~2,100 candidates per pair it finds

Instrumenting `count_overlaps` to record the candidate window it walks against
the pairs it emits:

| level | candidates scanned | pairs found | ratio |
|---|---|---|---|
| 0 | 30,823,436 | 12,104 | 2,547 |
| 1 | 347,985,182 | 122,626 | 2,838 |
| 2 | 4,366,188,947 | 2,080,634 | 2,099 |

4.4 **billion** AABB tests to find 2.1M pairs. The sort-axis interval overlap is
almost no filter at all on a draped cloth: everything overlaps everything along
any single axis, and the 3D test then rejects 99.95% of what the sweep walked.

**This is not a bad axis choice.** Forcing each axis at level 1:

| sort axis | ratio | bp_fv_ms |
|---|---|---|
| 0 (what `choose_axis` picks) | 2,882 | 118 |
| 1 | **126,063** | **3,587** |
| 2 | 2,838 | 120 |

The variance heuristic is already picking a near-best axis and is worth keeping --
it avoids a 44x disaster. The limit is that one-dimensional pruning cannot
separate this geometry, whichever dimension is chosen.

## Proposed variants

**1. A two-dimensional cell list, with no sorting anywhere.** This is the one
worth doing, and the reason is as much the sort as the scan ratio. Today's prep
sorts three AABB lists along the chosen axis -- `O(n log n)`, 237 ms at level 2 --
and the sweep gets a 2,100:1 scan ratio for it. Binning is `O(n)`.

**Two dimensions, not three.** The geometry here is surfaces, so with N elements
of size `h` in a box of side `L`, `N ~ (L/h)^2`. A 2D grid at that resolution has
~N cells and is mostly occupied; a 3D grid has `(L/h)^3 = N^1.5`:

| N | 2D cells | 3D cells | 3D occupancy |
|---|---|---|---|
| 92,230 | 9.2e4 | 2.8e7 | 0.33% |
| 368,920 | 3.7e5 | 2.2e8 | 0.16% |
| 1,475,680 | 1.5e6 | 1.8e9 | 0.08% |

1.8 billion cells at 0.08% occupancy buys pruning along an axis a surface barely
occupies. The third axis does not need a grid dimension: bin on two, and let the
ordinary AABB test inside the cell reject on the third. Full 3D rejection at 2D
cost.

`src/cell_broadphase.hpp` is a useful starting point but **does not get there as
written**, and the gap is worth being precise about:

* Its binning is already right in kind, just one axis short. `cell_populate` is a
  counting sort -- `cell_count` builds `cellptr`, then a bookkeeping pass scatters
  into it. `O(n)`, no comparison sort, and itself a count-then-fill.
* Its *query* still needs sorted input. `cell_count_overlaps` takes a
  `cell_list_axis` *and* a `sort_axis`: it bins on one, then inside each cell
  walks a window on another and breaks out on the first non-overlapping entry.
  That break is only correct if the cell's contents are ordered along `sort_axis`,
  so the sort it was meant to replace is still required.

So the work is: bin on the two best axes instead of one, and replace the
break-based window scan with a test over every box in the touched cells. Dropping
the break is what removes the ordering requirement, and it costs nothing when a
cell holds a handful of boxes. The existing SIMD block test and the CRS output
carry over unchanged.

**2. Do not replace count-then-fill.** An earlier revision of this file proposed
collecting into per-thread growable buffers and concatenating, to avoid sweeping
twice. That is wrong and is retracted. Counting then filling is a stable
performance pattern: the allocation is exact and known before anything is
written, there is no reallocation mid-sweep, no per-thread capacity to tune, no
merge whose cost depends on how pairs happened to distribute across threads, and
the CRS order is deterministic. The same structure works unchanged on the GPU,
where count / prefix-sum / fill is the standard idiom. The second sweep is the
price of predictable behaviour. Savings should come from *what the sweep walks*,
which is what (1) attacks.

**3. Morton ordering: demoted.** Sorting by Morton code would prune in three
dimensions while keeping the current machinery, but it is still a sort, so it
keeps the `O(n log n)` prep that (1) removes. Only worth revisiting if a grid
turns out to be unsuitable -- very non-uniform element sizes would be the reason.

**4. The armadillo-rollers edge-edge blowup stays a priority.** It is a separate
problem from scaling -- a pathological input where the conservative kernel runs
~25x slower than the same search on a CPU -- and being 4% of a refined cloth
scene does not make it less worth fixing. See `oracle/README.md` for what has
been ruled out so far.

## Implemented: `src/broadphase/sccd_broadphase_cell2d.hpp`

Both passes now have a cell-list path, selected by `SCCD_BROADPHASE=cell2d`. The
sweep remains the default. No sorting is done on this path at all -- vertices,
faces and edges are all binned, so `sort_along_axis` is not called.

Same scene, same levels, host:

| level | faces | | prep_ms | bp_fv_ms | bp_ee_ms | **broad_ms** |
|---|---|---|---|---|---|---|
| 0 | 92,230 | sweep | 48 | 9.3 | 26.0 | 84 |
| | | cell2d | 16 | 6.4 | 19.4 | **42** |
| 1 | 368,920 | sweep | 81 | 101.5 | 202.8 | 385 |
| | | cell2d | 49 | 41.3 | 96.1 | **186** |
| 2 | 1,475,680 | sweep | 234 | 1,730.1 | 2,503.2 | 4,467 |
| | | cell2d | 191 | **320.2** | **551.9** | **1,063** |

**4.2x on the broad phase at 1.5M triangles** -- face-vertex 5.4x, edge-edge 4.5x --
and the gain grows with size (2.0x, 2.1x, 4.2x), which is what replacing an
`O(n·k)` scan with `O(n)` binning should look like. The whole collision step goes
from 4,691 ms to 1,267 ms, **3.7x**.

Pair counts are **identical** at every level: 2,080,634 vertex-face and 6,481,554
edge-edge either way. That is the property that matters -- this is an
optimisation, so producing a different answer would make it worthless.

### Correctness

`src/tests/cell2d_broadphase_test` compares the two pair sets directly, on cases
picked for where they could disagree: boxes spanning many cells (the case
duplicate suppression exists for), boxes larger than the domain, fully dense
overlap, degenerate extents, and the self-overlap form. All eight agree exactly.
It needs neither smesh nor TightInclusion, so it runs in any build.

Duplicate suppression is worth understanding before touching this code. A box
spans several cells, so a pair can be met in several of them. The pair is
attributed to the cell holding the minimum corner of the two boxes' *overlap*:
that corner is inside both boxes, so both are binned in that cell, and the cell
is unique. Two clamps, no hash and no mark array. The self-overlap form adds
`j > i` on top, which the sweep gets for free by starting its window at `i + 1`
in sorted order.

### Not done

Edge-edge and face-vertex use separate grids, each sized from its own AABBs;
sharing one would save a binning pass. `prep_ms` improved only 1.23x because
binning replaces sorting rather than removing work outright.

The device kernels exist and are tested, but are **not yet wired into**
`sccd_smesh_ccd.hpp`'s device path -- `broad_phase_prep_device_` and friends still
sort and sweep. Connecting them is the same edit already made on the host side.

## Caveats

The scan-ratio and end-to-end tables are host only, one scene, one frame pair.
The device numbers above are synthetic distributions rather than a real scene,
because the device path is not yet wired up. The GPU sweep is a separate
implementation (`src/cuda/sccd_broadphase.cu`) and its scan ratio was not measured
directly, though `benchmark/alps/bench_aggregate.csv` shows the same qualitative
split -- broad phase dominating narrow by 4-75x across scenes. The instrumentation
used for the scan-ratio table was a scratch build, not committed: it counts only
the vertex-face `count_overlaps`, and the edge-edge sweep is assumed similar
because its cost per found pair tracks.
