# TODO

Known open work, with the evidence for each so the next person does not have to
re-derive it. Items leave this file when they are done or when a measurement
refutes them; a refuted item goes to `benchmark/ASSESSMENT.md` rather than being
deleted, so it is not retried.

---

## High

### Shrink the clone: 4.6 GiB of history for a 1.9 MB checkout

**The problem.** `git clone` of this repository transfers **4.56 GiB**. The
tracked checkout is **1.9 MB across 134 files**, and no tracked file exceeds
132 KB. All of the weight is history: datasets and benchmark archives that were
committed, then later covered by `.gitignore` in `e67d5105` without being removed
from past commits. Ignoring a path stops new blobs; it does not shrink the ones
already in the object store.

**Where it is.** Two path prefixes account for essentially all of it
(uncompressed blob totals over all reachable commits):

| Path | In history | Blobs |
|---|---:|---:|
| `data/` | 6,229 MB | 20,569 |
| `benchmark/alps/` | 744 MB | 133 |

The single largest objects:

| Size | Object |
|---:|---|
| 190 MB | `data/n-body-simulation_ee_table.csv` |
| 141 MB | `data/puffer-ball/queries/21ee.csv` |
| 130 MB | `data/puffer-ball/queries/22ee.csv` |
| 104 MB | `data/.downloads/cloth-funnel.tar.gz` |
| 98 MB | `data/puffer-ball_ee_table.csv` |
| 83 MB | `data/puffer-ball/roots/21ee_roots.tar.gz` |
| 82 MB | `benchmark/alps/sccd-benchmark-2026-05-26.tar.gz` |

**Why it is worth doing.** This is the first thing a new user experiences, and it
is a multi-gigabyte download before they can read the README. It is also paid
again by every CI job, every fork, and every clone on a cluster. Nothing in those
6.9 GB is needed to build, test or use the library: the default build has no
dependencies and no data, and `ctest` skips the dataset-backed tests when the
data is absent — which on a fresh clone it always is.

**What it takes.**

```sh
git filter-repo --path data --path benchmark/alps --invert-paths
```

Both prefixes are already in `.gitignore`, so nothing that ships is touched, and
the working tree after the rewrite is byte-identical to the one before it.

**Why it is not done yet, and what the decision is.** This rewrites every commit
hash on **12 published branches** (`clean-up`, `debugging`, `main`, `narrow`
through `narrow8`, `perf-optimizations`) and needs a force-push. Everyone
re-clones; any fork, open pull request, or CI pin against an old SHA breaks; the
`Claude-Session` and commit references in existing messages still resolve, but
SHAs quoted anywhere outside the repository will not.

That is a coordination call, not a technical one. Suggested sequence when it is
made:

1. Tag the current tip of every branch (`pre-rewrite/<branch>`) and push the tags,
   so the old history is recoverable by SHA for as long as anyone needs it.
2. Confirm no open pull requests.
3. Run `filter-repo` on a fresh mirror clone, not on a working checkout.
4. Verify: `git count-objects -vH`, then a fresh clone builds and passes `ctest`.
5. Force-push all branches and tell anyone with a clone to re-clone rather than
   pull.

**Where the data should live instead.** `data/README.md` already describes the
datasets; the download step belongs there, so a user who wants the benchmarks
fetches them and a user who wants the library does not.

---

---

## Done

- ~~The device global queue deadlocks under sustained overflow~~ — **fixed**. The
  queue is double-buffered: a launch reads one buffer and writes the other, so
  neither the writer's `atomicCAS` spin nor the reader's commit spin exists any
  more and both operations are wait-free. That let the shared stack shrink from
  1024 entries to 64, which is what actually redistributes heavy queries, and the
  conservative device kernel got 17× faster on cloth-funnel and 5.1× on
  armadillo-rollers. Written up under "Fixed: the queue is double-buffered" in
  `benchmark/ASSESSMENT.md`, including one unexplained regression at the old
  capacity.

- ~~The armadillo edge-edge blowup~~ — **withdrawn**. All 396 edge-edge cases,
  twice: the worst is 9.4×, none exceeds 10×, and mode 2 is 10% faster over the
  trajectory. Written up under "Withdrawn: mode 2 is about 100× slower on
  armadillo edge-edge" in `benchmark/ASSESSMENT.md`.

### Mode 2's earliest-impact answer may be late on armadillo-rollers edge-edge

**Possible conservativeness violation. Unresolved, and it is the one open item
that touches the invariant directly.**

`find_earliest_impact_time` (`toi_stride = 0`) at `SCCD_NARROWPHASE_MODE=2`
reports a time of impact *later* than the smallest exact root in the case's query
set, on 5 of 16 armadillo-rollers cases on the host and 6 of 16 on the device.
All five are edge-edge. Mode 0 is clean on the same cases, as is every scene at
either mode on cloth-funnel and cloth-ball.

| case | mode 2 | true earliest | over by | mode 0 |
|---|---|---|---|---|
| 100ee | 0.00103297 | 0.001031 | +2.0e-6 | 0.00101829 |
| 123ee | 0.0110207 | 0.0109835 | +3.7e-5 | 0.0109409 |
| 25ee | 0.00726114 | 0.00725506 | +6.1e-6 | 0.00724742 |
| 304ee | 0.000270735 | 0.000269097 | +1.6e-6 | 0.000263057 |
| 370ee | 0.00220261 | 0.00220123 | +1.4e-6 | 0.00218988 |

Reproduce with `sccd_bench` against a dataset that ships exact roots and read the
`s0_late` / `s0_toi` / `gt_earliest` columns.

**Two candidate explanations, and they need separating before anything is
concluded.** Either the TightInclusion-exact kernel's zero-stride path really can
report late — which the per-query path does not, since `toi_late` is 0 for the
same mode on the same queries — or the reference is wrong for these cases,
because the curated query set contains a pair that SCCD's broad phase never
examines. `broad_fn` is 0, but that is computed against a set built from the same
query list through `edge_id_map`, so a pair that fails to map could be missing
from both sides and leave the count clean. Checking whether the pair achieving
`gt_earliest` is in SCCD's broad-phase output settles it.

Until it is settled, do not claim the earliest-impact path is verified
conservative at mode 2.

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

## Lower

- **The device conservative search classifies 111 boxes per query where the host
  classifies 1.2.** Note this is a count of *work*, not of time: double-buffering
  the queue cut the time gap to 5.1×, 3.4× and better-than-host across the three
  scenes without changing the search at all. The remaining item is to make the
  search smaller, not better balanced. On the path `find_earliest_impact_time` uses: 97.8M boxes
  against 1.04M on cloth-funnel, **94×**. Same acceptance test, same tolerances,
  same queries, and not the price of the guarantee — the host is conservative too.
  Counted with `SCCD_NP_COUNT_BOXES` (off by default).

  A **separate and larger** defect sits on the `find_impact_times` path, which
  runs the block-per-query kernel: 274 queries there cost 706M boxes, 2.58 million
  each against the host's 1,844, a ratio of **1397×**. Start with its 128-way
  seeding dice, and re-run the measurement that kept the dice (single-root seeding
  was recorded as worse, 803 → 871 ms on armadillo edge-edge, which is hard to
  credit against 2.58M boxes per query).

  Ruled out along the way: the depth cap (neither side ever accepts at it), the
  split rule (both take the largest width/tolerance ratio), the acceptance test
  and the tolerances (transcribed). What is left is the `t` bound each query
  searches under and the order boxes leave a stack that mixes queries — the device
  tree *broadens* to 4.2M boxes at level 16 where the host has 7,000.

  Corner reuse, which this item used to propose, is **not available**: the kernel
  uses 238 of 255 registers with no spills, and carrying eight corners for three
  components needs 48 more. Three other explanations were measured and refuted —
  sequential batching, atomic contention on the shared time of impact, and
  per-box arithmetic — and are recorded so they are not retried.

  Counted per query, the gap **scales with query difficulty**: the two are within
  a few percent on queries needing 0–1 boxes, 2× apart at 8, 32× at 64 and 462× at
  16,384. So the device is not starting wrong, it is failing to converge on the
  queries that need real search.

- **`SCCD_NARROWPHASE_MODE` does not reach the quad path.** There is one quad
  root-finder variant, so the enum has nothing to select between. It now says so
  once on stderr rather than being silently ignored, which closes the trap but
  not the gap.

- **`broadphase_strategy.hpp` has no test.** It picks between the sweep and the
  cell list by racing them at run time, and nothing checks that the race
  converges on the faster one or that switching mid-run keeps the pair sets
  identical.

- **Two CUDA tests are built but never registered.** `mesh_sccd_cuda_test` and
  `sccd_cuda_cpu_gpu_parity` need mesh data that a fresh clone does not have.
  `sccd_add_raw_mesh_test` now reports what is missing instead of returning
  silently, so the gap is visible, but the tests still do not run anywhere.
