# TODO

Known open work, with the evidence for each so the next person does not have to
re-derive it. Items leave this file when they are done or when a measurement
refutes them; a refuted item goes to `wip/ASSESSMENT.md` rather than being
deleted, so it is not retried.

---

## High

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

- ~~Validation-only narrow-phase modes in the shipped API~~ — **removed.** Mode 1
  was the slowest kernel in the library on every scene (15.9 / 90.9 / 223.4 ms
  against Fast's 6.4 / 19.4 / 190.7) and is now in
  `spikes/src/sccd_narrowphase_fast_vector.hpp`. Mode 3 corrected mode 1 with
  TightInclusion; `SCCD_USE_TI=1` calls the library directly instead. The enum is
  `Fast` and `Tight`.


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


