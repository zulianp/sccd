# TODO

Known open work, with the evidence for each so the next person does not have to
re-derive it. Items leave this file when they are done or when a measurement
refutes them; a refuted item goes to `wip/ASSESSMENT.md` rather than being
deleted, so it is not retried.

---

## Notes that are not work items

* [`wip/ALPS.md`](ALPS.md) — how to run the tests on the CSCS Alps cluster.
  `ctest` under `srun` stalls after the first test; run the binaries directly.
* [`wip/CUDA_NARROWPHASE_PLAN.md`](CUDA_NARROWPHASE_PLAN.md) — the device narrow
  phase's remaining work: the 94× box-count gap on the earliest-impact path, what
  has already been ruled out, and the variants worth trying. Branch `cuda-op-np`.

---

## High

### Nothing tests the install tree

Three defects shipped in the exported CMake package at once, and each was
invisible for the same reason: no test ever consumed an install. The exported
target advertised `include/cuda` in builds that install no `.cuh`, so
`find_package(SCCD)` failed outright on every non-CUDA install; `SCCDConfig.cmake`
never found the CUDA toolkit though the target links `CUDA::cudart`; and Thrust
was linked `PUBLIC`, exporting a target name a consumer cannot define, which CMake
degraded to a literal `-lThrust`.

All three are fixed and were verified by hand -- a C++ consumer through
`find_package` against a non-CUDA install, and a `.cu` consumer against a CUDA
install on GH200. **That verification is not automated.** A ctest that installs
to a temporary prefix and compiles a three-line consumer against it would have
caught all three, and would catch the next one. The header-classification guard
already fails configure on an unclassified header; this is the same idea one step
later.

### Naming, the last of the API pass

Left from the API work, all mechanical, none behavioural:

* The public narrow-phase entry points spell their overlap arrays `voveralp`,
  `foveralp`, `qoveralp`, `e0overalp`, `e1overalp` -- a transposition typo -- while
  the internal kernel they forward to spells the same parameters correctly.
  Both spellings appear in one call chain.
* One concept has five names for stride (`face_stride`, `edge_stride`,
  `quad_stride`, `first_stride`/`second_stride`, `stride`) and the output pair has
  three (`foverlap`/`noverlap`, `first_out`/`second_out`, `out0`/`out1`).
  `noverlap` also reads as a count while colliding with `noverlaps`, which is one.
* `CCD<T>::narrow_phase_fv` against the core's `narrow_phase_vf`.

### Another spatial dimension needs the overlap predicate, not a parameter

`SCCD_DIM` is a constant rather than a parameter because `sccd::disjoint` and
`vaabb_overlap_one_to_many_bits` take x, y and z as separate positional arguments
and their SIMD forms load three min rows and three max rows. Measured at `dim = 2`
with a four-row table, the sweep returns **zero pairs**: it reads row 2, an x
maximum, as a z minimum. Generalising the predicate is the work; it is the broad
phase's hot loop, so it needs its own benchmarking rather than being folded into
an API change.


### ~~Shrink the clone: 4.6 GiB of history~~ — withdrawn, the measurement was wrong

**There was never 4.6 GiB to shrink.** A real `git clone` of this repository
transfers **24 MB**.

The claim came from running `git count-objects -vH` on a local working clone, and
from `git clone --local .`, which hardlinks the *entire local object store*
whether or not it is reachable. Neither measures what a clone from the remote
costs. The check that settles it is a mirror clone from the URL, which nobody ran
until the rewrite was about to happen:

| | |
|---|---:|
| real `git clone` from GitHub | 24 MB |
| `data/` objects in remote history | 2 |
| `data/` objects in the local store | 24,963 |

The 4.56 GiB was entirely local, and almost all of it was pinned by a single
stale `refs/codex/turn-diffs/checkpoints/...` ref left behind by tooling. Dropping
that ref and running `git gc --prune=now` took the local store from 4.56 GiB to
**24.03 MiB** with every branch intact — no rewrite, no force-push, nothing for
anyone else to re-clone.

What remains on the remote is 4 `csv/` blobs in no branch tip, worth perhaps 10 MB
of a 24 MB clone. Not worth a force-push across 12 branches, and not planned.

**The lesson, since it cost a planned history rewrite:** `git count-objects`
measures a local object store, not a repository. Before proposing surgery on
shared history, clone the URL and measure that.

---

## Done

- ~~The vertex-quad path lagged the triangle path~~ — **done.** The certified
  error bound had its clamp inverted, `max(max_coord, 1)^3` where TightInclusion
  uses `min(...)`, which inflated the pad by the cube of the scene size: the
  scale-invariance check drifted 6x at scale 1e3 and 618x at 1e4 and is now flat.
  The tolerance caps from `sccd_tolerance.hpp` are applied, the split scale
  follows the pruned time window (12% at `ToiOutput::Earliest` on a moving quad),
  and six divergences from the triangle path are gone, including a re-created
  seq_cst atomic and a raw OpenMP loop that ran single-threaded in a TBB build.

- ~~The device quad kernel's local DFS stack~~ — **built, measured, rejected.**
  The restructure onto a block-shared pool works and is conservative, and removes
  the spill (8128 B frame and 112/216 B spill become 336 B and none), but wins one
  measured cell of four and halves throughput in another. `wip/DECISIONS.md` has
  the numbers and why. Reverted; the header extraction it needed,
  `sccd_device_dfs_stack.cuh`, stands on its own.

- ~~The assessment discarded every refine row~~ — **fixed.** `run_refine`
  filtered with `awk 'NF == 10'` and `sccd_refine_scaling` prints eleven fields,
  so refine successes were dropped while the `rc != 0` branch kept writing failure
  rows. That is why `FAILED rc=134` survived in the hopper quad rows with nothing
  to contradict it.


- ~~Validation-only narrow-phase modes in the shipped API~~ — **removed.** Mode 1
  was the slowest kernel in the library on every scene (15.9 / 90.9 / 223.4 ms
  against Fast's 6.4 / 19.4 / 190.7) and is now in
  `spikes/src/sccd_narrowphase_fast_vector.hpp`. Mode 3 corrected mode 1 with
  TightInclusion; `SCCD_USE_TI=1` calls the library directly instead. The enum is
  `Relaxed` and `Tight`.


- ~~Device dead code in the shipped CUDA build~~ — **done**, in
  `spikes/src/dead.cuh`. Verified: clean CUDA build and ctest 9/9 on GH200.


- ~~The refactor's CUDA build was unverified~~ — **verified** on GH200,
  2026-09-04, after the certificate was re-signed. Clean compile with zero
  errors; `ctest` 9/9 including `sccd_narrowphase_cuda_test`; and all twelve
  configurations reproduce `benchmark/assessment/mode-stride-matrix.csv`
  exactly — every `fp`, `fn`, root count and `late`, across 101,164 exact roots.
  One timing outlier worth a repeat if anyone cares: `armadillo-rollers` mode 2
  device at `toi_stride=1` measured 3776.8 ms against a 3126.7 ms three-repeat
  median, ~21% high. Its `toi_stride=0` companion is 56.3 ms before and after,
  so this is the block-per-query kernel's variance rather than the refactor.


- ~~Mode 2's earliest-impact answer may be late on armadillo-rollers~~ —
  **withdrawn, the reference was wrong.** smesh stores coordinates as `float`
  (`typedef float geom_t`), so `find_earliest_impact_time` runs on float32-rounded
  geometry, while the exact roots are computed for the dyadic rationals in
  `queries/*.csv`. Only armadillo's coordinates need more than 24 mantissa bits,
  and armadillo is the only scene the check fired on. SCCD's own
  mode-2 answers on the two geometries scatter apart by 5e-8 to 4.6e-5 in *both*
  directions, and on the CSV geometry the kernel is conservative on 16 of 16
  cases. An unsound kernel would err one way and would show on the CSV path too.
  Checking the mesh path needs exact roots for the mesh geometry, which the
  datasets do not ship.

- ~~The device global queue deadlocks under sustained overflow~~ — **fixed**. The
  queue is double-buffered: a launch reads one buffer and writes the other, so
  neither the writer's `atomicCAS` spin nor the reader's commit spin exists any
  more and both operations are wait-free. That let the shared stack shrink from
  1024 entries to 64, which is what actually redistributes heavy queries, and the
  conservative device kernel got 17× faster on cloth-funnel and 5.1× on
  armadillo-rollers. Written up under "Fixed: the queue is double-buffered" in
  `wip/ASSESSMENT.md`, including one unexplained regression at the old
  capacity.

- ~~The armadillo edge-edge blowup~~ — **withdrawn**. All 396 edge-edge cases,
  twice: the worst is 9.4×, none exceeds 10×, and mode 2 is 10% faster over the
  trajectory. Written up under "Withdrawn: mode 2 is about 100× slower on
  armadillo edge-edge" in `wip/ASSESSMENT.md`.

### The benchmark reports `fn=0` on datasets that have no ground truth

`write_case_outputs` sets `expected` from `mma_bool` when it exists and from
`roots/<key>/toi.float64` otherwise. With neither present it is false for every
query, so `fn` is zero by construction and `fp` degenerates into a count of
collisions reported. Nothing in the output says the comparison had nothing to
compare against.

This is not hypothetical: it invalidated every `fn=0` reading in this branch
until the full matrix was run against a tree that does ship roots. A
conservativeness claim rests on that column, so the benchmark should refuse to
print it rather than print a zero — emit `fn=-1`, or a `gt=none` marker, or fail
outright, but not a clean zero that reads as a pass. `ti_oracle` already gets
this right; `sccd_bench` does not.

Cheap to fix and worth doing before the next accuracy claim is made from a
benchmark run.


