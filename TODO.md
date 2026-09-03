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

## Done

- ~~The armadillo edge-edge blowup~~ — **withdrawn**. All 396 edge-edge cases,
  twice: the worst is 9.4×, none exceeds 10×, and mode 2 is 10% faster over the
  trajectory. Written up under "Withdrawn: mode 2 is about 100× slower on
  armadillo edge-edge" in `benchmark/ASSESSMENT.md`.

## Lower

- **The device conservative search classifies 944 boxes per query where the host
  classifies 1.8.** Same acceptance test, same tolerances, same queries: 796.3M
  boxes against 1.49M on cloth-funnel, a **535×** work difference. That is the
  whole of the 26× time gap and then some, and it is not the price of the
  guarantee — the host is conservative too, on the same inputs, in under two boxes
  per query. Counted with `SCCD_NP_COUNT_BOXES` (off by default); the numbers and
  the method are under "Counted" in `benchmark/ASSESSMENT.md`.

  Corner reuse, which this item used to propose, is **not available**: the kernel
  uses 238 of 255 registers with no spills, and carrying eight corners for three
  components needs 48 more. Three other explanations were measured and refuted —
  sequential batching, atomic contention on the shared time of impact, and
  per-box arithmetic — and are recorded so they are not retried. Start from the
  535×, not from the arithmetic.

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
