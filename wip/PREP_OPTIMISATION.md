# Preparation phase: what was changed and what it measured

Section 8.7 of the article reports that host preparation stops improving at
about `2.4x` and stays there from sixteen cores to seventy-two, and names three
serial passes as the cause. This records the work done against that finding, and
what it is worth.

## What the passes actually cost

Profiling `broad_phase_prep_host_` on cloth-funnel (9,450 nodes, 18,484 faces,
27,930 edges) found the three named serial passes to be real but not the whole
story. Two other costs were larger:

- `compute_aabbs` opened **three parallel regions per call**, one per dimension,
  and prep calls it three times. Nine forks and nine joins for a mesh whose
  per-element work is measured in microseconds.
- `sort_along_axis` opened eight regions per call, and the comparison sort
  inside it did not scale: `parallel_sort` reaches only `2.0x` to `2.4x` at ten
  cores, and the sort is `93%` of the sweep path's preparation.

## The changes

| Change | Where |
|---|---|
| One parallel region per `compute_aabbs` call, dimension loop inside it | `src/core/sccd_aabb.hpp` |
| `parallel_for_chunks`, a parallel loop for short ranges | `src/core/sccd_parallel.hpp` |
| `parallel_tiled_reduce`, a deterministic tiled reduction | `src/core/sccd_parallel.hpp` |
| `SCCD_MIN_PARALLEL_N`, a guard so short passes skip the fork | `src/core/sccd_parallel.hpp` |
| `center_variance` as a tiled reduction | `src/broadphase/sccd_broadphase_sweep.hpp` |
| Parallel radix sort of the (key, index) pairs | `src/broadphase/sccd_broadphase_sweep.hpp` |
| Grid sizing as a tiled reduction | `src/broadphase/sccd_broadphase_cell2d.hpp` |
| `Cell2DPartition`: binning partitioned by cell row, so blocks write disjoint slices | `src/broadphase/sccd_broadphase_cell2d.hpp` |
| The unused `choose_axis` removed from the cell-list path; redundant zero-fill of `cellptr` removed; index fills parallelised | `src/integrations/smesh/sccd_smesh_ccd.hpp` |

`parallel_for_br` tiles its range into fixed blocks of 128, so a loop shorter
than that ran on one worker. The first attempt at the partitioned binning was
*slower* for exactly that reason: the block loops collapsed to one thread while
still paying for the extra passes. `parallel_for_chunks` exists for loops whose
range length is the worker count.

## What is guaranteed, and how it is checked

The partitioned binning produces **byte-identical** `cellptr` and `cellidx` to
the serial binning, because blocks cover cell rows in order and boxes keep index
order inside a block. The radix sort produces **exactly** the permutation
`std::sort` produces under the same total order, because it is stable and the
pairs arrive in index order. Both are asserted directly by
`sccd_broadphase_cell2d_test`, at sizes that cross the parallel thresholds, and
the sort case is run on geometry with heavy coordinate ties.

On \textsc{GH200}, every correctness column of the benchmark (false positives,
false negatives, broad-phase false positives and negatives, time-of-impact
count, late count, maximum late, maximum early, median early, ground-truth
earliest, root count) is identical between the two builds on all 64 rows of both
scenes. Zero missed and zero late throughout.

## Measured on one Grace, 72 threads

A GH200 node carries four Grace modules, so `nproc` reports 288 and a harness
that trusts it measures four processors while reporting one. Everything below is
a single Grace at 72 threads, which is the configuration this project reports.

The thread count is not a detail. Preparation behaves completely differently
across that range, in the unmodified code as much as in this one: on
armadillo-rollers it measures `67 ms` at 72 threads and `341 ms` at 288 on the
same build. Work tuned against 288 threads is therefore worse at 72, which is
where this landed once before -- see the grain note below.

Four cases per scene, Tight, tuner on, two repeats, minimum taken; milliseconds
summed over the four cases. `benchmark/results/profile/prep-ab2/`.

| scene | before | after | gain |
|---|---:|---:|---:|
| armadillo-rollers | 67.3 | 38.3 | 1.76x |
| cloth-ball | 232.2 | 60.3 | 3.85x |

The broad-phase traversal is unchanged in both, `79.4` against `79.2` and
`1068.2` against `1055.1`, which is the control that reads this as a preparation
result.

### The grain, and why the elegant form is the wrong one

The thresholds that decide when a pass is worth a parallel region are fixed
element counts. Sizing the team from the work instead -- one worker per few
hundred elements, so a wider machine narrows the team -- looks more principled
and measures four times *slower* at 72 threads on armadillo-rollers, `72.2`
against `300.5 ms`. The preparation passes are tens of thousands of elements, so
a per-worker grain caps the partition at a handful of blocks on exactly the sizes
that matter, and leaves most of the machine idle. Both files record this where
the constants are defined, because the per-worker form is what a reader will
propose.

## Refinement scaling, where the change is largest

`sccd_refine_scaling` refines cloth-ball frames 0 to 1 five times, quadrupling
the element count at each level. Both builds were run on one node with the same
inputs and the same parameters. Preparation, in milliseconds:

| level | elements | before | after | gain |
|---:|---:|---:|---:|---:|
| 0 |     92,230 |    6.12 |   0.92 |  6.7x |
| 1 |    368,920 |   33.70 |   5.28 |  6.4x |
| 2 |  1,475,680 |  141.42 |  14.70 |  9.6x |
| 3 |  5,902,720 |  759.87 |  59.99 | 12.7x |
| 4 | 23,610,880 | 4139.35 | 243.23 | **17.0x** |

The gain grows with the element count, which is the signature of removing serial
passes rather than constant overhead. The whole broad phase at level 4 goes from
`4904.6` to `1020.7 ms`. On the sweep path preparation goes from `2150.3` to
`483.0 ms` at the same level, a factor of `4.5`.

Two controls make that readable as a preparation result and not a traversal one.
The traversal columns do not move: `bp_fv` goes `211.6` to `220.1` and `bp_ee`
`553.6` to `557.4` on the cell list, and on the sweep the two together go
`11{,}133` to `11{,}156 ms`, a fifth of a percent. And the axis diagnostics are
identical at every level, `lambda 0.0/836.8/0.0` with anisotropy `2994247.16` at
level 0 on both builds, so the tiled reduction that replaced the serial
centre-variance pass chooses the same axis it used to.

### A second provenance problem in the same files

The committed `benchmark/results/scaling/*.txt` reports `6,069,281` edge-edge
pairs at level 4. Both the current build and the unmodified baseline report
`6,108,123` on the same frames, so the committed file cannot be reproduced from
this tree by either. The difference is in the refinement, not in the broad
phase, and it predates this work: it is a build-vintage difference that the
files did not record. The regenerated files are produced by one binary from
inputs their header now names.

## The device path

The device builds the same boxes and sorts them with CUB, and at the largest
refinement level its preparation costs `1682 ms` against the host's `243`. The
three sorts account for `739` of that (`129`, `248` and `362` for the vertex,
face and edge lists); the rest is the AABB pass.

`compute_aabbs_face_kernel` accumulated each bound through `aabbs[d][e]` in
global memory, once per vertex per dimension: about thirty memory operations for
a triangle where six writes do, with the row pointers themselves reloaded from
device memory every time they are named. It now reduces over the element's
vertices in registers and writes each bound once, which is the shape the host
kernel already had.

Rounding outward once at the end gives the same box as rounding every vertex,
because `dnextafter_down` and `dnextafter_up` are monotone and the minimum of
the rounded values is the rounded minimum. The box is identical, not merely
close, which is what the conservativeness argument needs.

Both builds in one job, cloth-ball frames 0 to 1 at the largest level:

| | before | after | |
|---|---:|---:|---|
| structure | 3150.78 | 2060.60 | 1.53x |
| traversal, face-vertex | 340.32 | 332.25 | unchanged |
| traversal, edge-edge | 1397.64 | 1398.18 | unchanged |
| broad total | 4888.73 | 3791.03 | 1.29x |

The traversal columns holding still is the control that reads this as an AABB
result. Device timings vary widely between jobs -- an earlier run of the same
unmodified binary reported `1858.5 ms` for the structure where this one reports
`3150.8` -- so only the within-job comparison is worth quoting.

### The answers are unchanged, and the check needed a control

A device benchmark run diffed between the two builds reports differing
`toi_max_early` and `toi_med_early` on several cases. That is not the change:
running the *same* baseline binary twice produces the same differences, on the
same rows. The benchmark's narrow phase uses `Earliest` output, where workers
prune against a shared atomic minimum, so the last digits depend on the order in
which they publish -- the behaviour Section 9 of the article already records.

The safety columns are identical between the builds: false positives, false
negatives, broad-phase false positives and negatives, and late counts all match.
Those are the columns that carry the invariant, and they are deterministic.

## What still limits preparation

On the host, preparation still plateaus, at `3.5x` to `5.5x` at thirty-two
cores. Two costs remain, and neither is one of the three passes Section 8.7
names.

The first is in the harness. `sccd_bench` constructs a fresh `CCD` per case, so
`init()` -- the edge-graph build and every allocation -- lands inside the timed
region for every case after the first. Measured at seventy-two cores on
cloth-funnel, the warmed first case reports `1.32 ms` against `1.46` to `1.49`
for the other three, so `init()` is about `10%` of what the prep column reports.
A solver calling `broad_phase_prep` on a persistent object each step pays it
once, not per step, so the prep column overstates per-step preparation and its
serial fraction.

The second is parallel-region count. At seventy-two cores on cloth-funnel,
preparation is flat at about `1.4 ms` per case whether the case carries `27,918`
queries or `388,206`, while the broad phase over the same cases moves from
`0.79` to `5.42 ms`. Preparation is the same work every case, so a constant is
expected; what the constant says is that the remaining cost is fixed overhead
rather than work over elements.

## Folding this into the article

The article's tables all come from one sweep, and preparation appears in
`tab:processor`, `tab:per-frame`, `tab:broadphase`, `tab:scaling`,
`fig:breakdown` and `fig:strong`. Updating only the strong-scaling figure would
leave the document inconsistent with itself. Absorbing this work means re-running
the whole sweep -- six scenes, two processors, two modes -- and regenerating
every table from it. Nothing in the article has been changed on the strength of
the measurements above.

## Four device defects, and the harness bug that hid them

Regenerating the sweep surfaced two separate problems on the device path: a
harness bug that published truncated chunks, and an intermittent illegal memory
access that killed about one chunk in seventeen. The harness bug is fixed. Four
device defects were found and repaired, one of them confirmed by a device assert
under a diagnostic build; which of them produced the illegal access is being
settled by an exposure run that puts the unfixed build, a build carrying only the
confirmed repair, and the fully repaired build through the same case list on the
same node.

**The harness published truncated chunks.** `run_chunk_body` in
`benchmark/scripts/sweep.sh` wrote each chunk to a temporary file and renamed it
on success, guarding only against a file with nothing but a header. It never
inspected the benchmark's exit status: the binary sits on the left of a pipe
into `tail`, so the status belonged to `tail`, and a run that wrote real rows and
then died was renamed and merged like any other. Eleven chunks were published
that way, all on the device, the worst being
`cloth-ball/device/sweep/r2`, which contributed `12` rows of an expected `158`.

That matters more than a lost chunk. Truncation removes whichever cases follow
the fault, so what survives is a biased subset rather than a thinned one, and
every median and maximum taken over it moves with no sign that anything is
missing. Host combinations were untouched, which is what made the pattern
legible: they match the committed row counts exactly, combination for
combination.

The committed `sweep-gh200-bp.csv` is not affected. Checked the same way, every
one of its (scene, mode, broad phase) combinations carries identical counts
across all four mode-and-space pairs. The defect was latent there and this run
is what triggered it.

The chunk body now checks the exit status of each mode under `pipefail` and
refuses to publish, leaving the chunk unfinished for the next pass.

**The fault had no single detection site, because the detection sites were
meaningless.** `SCCD_CUDA_SYNC_FOR_CHECK()` is a no-op under `NDEBUG`
(`src/cuda/sccd_cuda_base.cuh:145`), so every `SCCD_CHECK_CUDA` in the build the
sweeps use is a launch check rather than a fault check. Kernels run
asynchronously; the first real synchronisation after each stage is a `cudaMemcpy`
read-back, and all three of the read-backs on the device path discarded their
status. An error therefore stayed sticky until whichever `SCCD_*` macro came
next, which is why one fault appeared at `sccd_broadphase.cu:189`,
`sccd_broadphase.cu:495` and `sccd_narrowphase.cu:2699` in different runs. The
line numbers described where the code next asked, not where anything went wrong.

The mechanism that was available all along is in that file's own comment:
*"Define `SCCD_CUDA_SYNC_CHECKS` to get it back when chasing an async fault to
its launch."* It was never wired to anything. It is now a CMake option beside
`SCCD_CUDA_VERBOSE_PTXAS`, and a build that sets it alongside live device asserts
found the defect in forty-five seconds.

**A block reduction read shared memory no warp had written.**

```
sccd_reduce.cuh:54: block_reduce_to_gmem(T, T *, T *) [with T = double]:
block: [156,0,0], thread: [0,0,0] Assertion `acc == acc` failed.
CUDA error: device-side assert triggered in sccd_broadphase.cu:189
```

`acc == acc` is a NaN test, and `sccd_broadphase.cu:189` is the
`SCCD_CUDA_LAST_ERROR()` at the head of `sort_along_axis` -- the next check after
`choose_axis`, exactly as the paragraph above predicts.

`block_reduce_to_gmem` reduces within each warp, has lane 0 of each warp write
its partial sum to `block_accumulator[warp_id]`, and then has warp 0 read every
slot back:

```c++
if (!lid) block_accumulator[warp_id] = acc;
__syncthreads();
if (!warp_id) {
    acc = lid < n_warps ? block_accumulator[lid] : 0;
```

Both `choose_axis` kernels opened with `if (i >= n) return;`. A warp lying wholly
past the end therefore never reached the write, while warp 0 read its slot
regardless -- and `__shared__` is not zero-initialised, so the block's sum
carried whatever that multiprocessor's shared memory happened to hold, normally
the previous launch's partial sums. Blocks are 256 threads and the last block is
partial for almost every mesh, so the variance that chooses the sweep axis was
contaminated on essentially every call. Whether the contamination changed the
answer depended on residue, and residue depends on what ran on that
multiprocessor before: that is where the intermittency came from. It reproduces
within 133 of 200 rod-twist device cases, every time.

The repair is to keep every thread of the block in the reduction and have the
out-of-range ones contribute zero, which leaves the sum identical and writes
every slot. The alternative, naming only the live lanes with `__activemask()`,
fixes the shuffles but not the unwritten slots, so it is the wrong half of the
problem. `src/tests/cuda/sccd_choose_axis_cuda_test.exe.cpp` fills the slots with
values far above anything the probe can produce and then asks for the axis of a
small input whose answer is unambiguous, sixty-four times. On the unfixed code it
fails all sixty-four, returning axis `0` where the answer is `2`; on the repaired
code it passes all sixty-four. A defect that looked intermittent across a sweep
is deterministic once the residue is controlled, which is the point of writing
the test that way rather than repeating a failing case.

**Three further defects, each wrong on its own terms.** None needed the fault to
justify fixing, and two of them are what turned a first error into the reported
signature.

*The global stack committed a capacity it had not allocated.* `grow_stack`
(`src/cuda/sccd_narrowphase.cu`) freed its eight arrays, issued eight unchecked
`cudaMalloc`s and set `gstack_cap = new_cap` whatever happened. The device push's
only guard is that capacity (`sccd_device_dfs_stack.cuh`), so one failed
allocation left null pointers with a capacity saying otherwise and the next
launch scattered through null -- surfacing at `read_cursor`'s copy, which is
detection site `sccd_narrowphase.cu:2699`. Failure was plausible rather than
theoretical because the growth request is `g_request`, which counts every
*dropped* push summed over Pass 1 and every drain round, not the boxes that must
be resident: on the first narrow-phase call of a process the stack is empty, so
every push is dropped and the deficit approximates the whole search tree. At 112
bytes per unit of capacity a deficit of 10^8 asks for 11 GB. It now allocates
into temporaries, commits only once all eight succeed, keeps the working stack
otherwise, and clamps the request against `cudaMemGetInfo`. Declining to grow is
safe: the retry count is bounded and a dropped box can only leave the time of
impact too early.

*The overlap-count read-backs were unchecked.* In
`src/integrations/smesh/sccd_smesh_ccd.hpp` the count is copied into a variable
initialised to `0`. If an earlier kernel had faulted the copy failed, the count
stayed zero, `create_buffer(0)` produced an empty buffer, and the collect kernel
was launched anyway with a null output pointer while `ccdptr` still held real
non-zero offsets -- so every thread with work wrote through null. That is what
made the failures look like a disagreement between the counting and filling
passes, and it is why the two kernels read identically when compared line by
line: they never disagreed. The three read-backs are checked now, and the collect
is skipped when the count is zero.

*The device stack cursors could wrap.* `push_global` and `pop_global` bump their
cursor on every attempt, failed ones included, so a launch attempting more than
2^31 pushes wraps it negative and both `>=` guards pass. The host side already
clamped; the device side does so now.

**What the fix does not change.** Safety is untouched: across all six scenes,
both broad phases, on the device at mode 0, the false positives, false negatives,
broad-phase agreement and late times of impact are identical row for row with and
without the repair, as are the overlap counts. That is expected -- sweep and prune
enumerates the same pairs whichever axis it sorts on, so the axis decides how much
work finding them costs, never what is found.

Cost is within the noise. A single measurement over sixty cases moved
broad-phase time by between `-6.6%` and `+1.9%` depending on the scene, and the
largest movement of the twelve was rod-twist under the **cell list**, which never
calls `choose_axis` at all -- so a swing of that size is what the measurement
does on its own. Repeating armadillo, the one scene whose sweep looked like it
had moved, four times each way settles it:

| broad phase | unfixed | repaired | difference |
| --- | --- | --- | --- |
| sweep and prune | `102.94 +- 3.74` ms | `101.52 +- 3.55` ms | `-1.4%` |
| cell list (control) | `103.33 +- 4.01` ms | `101.46 +- 3.59` ms | `-1.8%` |

The apparent gain under sweep and prune is smaller than the one the cell list
shows, and the cell list cannot have gained anything: it does not call the
repaired code. Both sit well inside one standard deviation of about `3.6%`.

**Settled against the build that produced the published numbers.** The check
above compares the repair against its immediate predecessor; what the article
needs is a comparison against the build the sweep actually ran, including the
CUDA flag changes made alongside the repair. Three binaries over the same cases,
five interleaved repeats each: `pub`, the sweep-time source and flags; `fix`,
`pub` plus the four device repairs; `opt`, `fix` plus the flag changes, which is
what ships. Device, sweep and prune, 200 cases.

| scene | mode | broad `opt-pub` | narrow `opt-pub` |
| --- | --- | --- | --- |
| armadillo-rollers | Relaxed | `-0.4%` (0.1 sd) | `+1.0%` (0.6 sd) |
| armadillo-rollers | Tight | `-0.0%` (0.0 sd) | `-0.0%` (0.0 sd) |
| cloth-ball | Relaxed | `-0.2%` (0.8 sd) | `+1.7%` (1.5 sd) |
| cloth-ball | Tight | `-0.7%` (0.5 sd) | `-1.1%` (0.7 sd) |

Every comparison is inside 1.5 standard deviations, and every minimum is within
`1.7%`. **The published device numbers stand, and no re-run is needed to
correct them.**

Two things that took a wrong answer first are worth recording, because both are
easy to repeat. Arms have to be **interleaved**, not run in blocks: a first
attempt ran all of `pub`, then all of `fix`, then all of `opt`, and reported a
`+3.8%` regression at `3.4` standard deviations on armadillo's broad phase --
which was drift over the job landing on whichever arm went last. Interleaving the
same measurement gives `-0.4%`. And the **minimum** belongs beside the mean: it
is the run least disturbed by whatever else the machine was doing, and in the
blocked attempt the minima already said `0.5%` where the means said `3.8%`.

The narrow-phase stack was the specific worry, since `grow_stack` now clamps its
request against free device memory and allocates before freeing, and it runs on
the first narrow-phase call of every process. It never bound: `0` refused
growths and `0` give-ups across both scenes, both modes and all fifteen runs.

That the cost barely moves is itself informative. Where the contaminated
variance changed the answer it collapsed the argmax to axis `0`, because adding
the same large residue to all three variances makes them indistinguishable in
double and the first comparison then wins. The device was therefore sorting along
axis `0` whatever the geometry, and on these six scenes that happens to cost
about what the correct choice costs. The defect was free here; it was not going
to stay free.

