# Decision record

Claims this project made, measured, and in several cases withdrew. Kept because a
retraction is only useful if it is as findable as the claim was, and because
every one of these was believed at the time on evidence that looked sufficient.

## Withdrawn or corrected

- [Withdrawn: "mode 2 is about 100× slower on armadillo edge-edge"](ASSESSMENT.md)
- [Retracted: "the device narrow phase loses on every scene"](ASSESSMENT.md)
- [Fixed on the way past: an unsound rejection in the device's mode-0 kernel](ASSESSMENT.md)
- [Withdrawn: "mode 2's earliest impact is late on armadillo-rollers"](ASSESSMENT.md)
- [Retracted: every `fn=0` in this document before the matrix above](ASSESSMENT.md)
- [Resolved: the sweep dropped touching pairs](ASSESSMENT.md)
- Reversed: demoting `external/json` to `spikes/` — see below
- Withdrawn: "the quad device kernel's local stack costs it 255 registers and a
  spill" — see below
- Corrected: "Scalable CCD misses 74% of curated contacts" and its timings — see
  section 8
- Withdrawn, then withdrawn again: the per-pair comparison with Additive CCD —
  see section 9
- Demoted: a centroid-binned cell list, at one level and at two — see below
- Corrected twice: "binning self pairs at the minimum corner loses pairs", then
  "it is sound but too costly" — it is sound over a total order on cells, and
  with the cell bounds used to prune it tests half to a fifth of the candidates
  the shipped scheme does — see below

The full argument and the numbers behind each sit in
[`ASSESSMENT.md`](ASSESSMENT.md).

## What kept going wrong

Five distinct failures, all of the same shape: a measurement that answered a
different question than the one being asked.

1. **Comparing across code paths.** "The device narrow phase loses by 87x" raced
   the host's *fastest* kernel against the device's *slowest*, because
   `SCCD_NARROWPHASE_MODE=2` names a different kernel on each side.
2. **Comparing across populations.** "The device does 520x the host's work" was
   an average over two paths differing by 15x; split by `toi_stride` it is 94x on
   the interface anyone calls and 1397x on a path handling 274 queries.
3. **Comparing against a reference that does not apply.** "Mode 2 reports late on
   armadillo-rollers" compared a mesh-path result against roots computed for the
   query geometry. smesh stores coordinates as `float`, so those are different
   numbers, and only armadillo's coordinates need more than 24 mantissa bits.
4. **Comparing against nothing at all.** Every `fn=0` in this branch before the
   full matrix was vacuous: that data tree ships no exact roots, so `expected`
   was false for every query and a false negative was impossible by construction.

Two independent instruments agreeing is not enough — in (2) they agreed because
both were making the same mistake. The check that has actually caught things is
asking what the number would look like if the code were wrong.

## 5. Measuring the wrong object entirely

Added after the fact, because it nearly cost a history rewrite across 12
branches.

"`git clone` of this repository transfers 4.6 GiB" was recorded as the
top-priority item, with a filter-repo plan and a force-push sequence. It came
from `git count-objects -vH` on a local working clone and from
`git clone --local .`, which hardlinks the whole local object store whether or not
it is reachable. A real clone from the URL is **24 MB**; the remote history holds
2 `data/` objects against 24,963 locally.

The 4.56 GiB was local, and nearly all of it was pinned by one stale
`refs/codex/...` ref left by tooling. Dropping it and running `git gc` gave
24.03 MiB with every branch intact.

Same shape as the other four: a number that answered a different question than
the one being asked. It survived longer than the others because it was never
challenged — it was measured once, written down with a table of largest blobs
that made it look thoroughly established, and cited unchallenged for several
turns.

## 6. Checking one caller and calling it none

`external/json` was demoted to `spikes/` in `e313109` on the stated grounds that
it was "never `add_subdirectory`'d, so it was not in the build at all". That
sentence is true and the conclusion drawn from it is not. The main
`CMakeLists.txt` does not reference it, which is what I checked; `benchmark/scripts/bench.sh`
configures it as a **standalone CMake project**, which I did not.

Its two programs turn the datasets' `boxes/*.json` and `mma_bool/*.json` into the
raw arrays `bench.exe.cpp` reads — `boxes/<key>/c0.int32`, `c1.int32` and
`mma_bool/<key>/mma_bool.uint8`. Without them `sccd_bench` has no expected pair
sets and no hit/miss booleans, so the whole accuracy half of the benchmark is
gone. And `bench.sh` did not degrade gracefully: it configures that path
unconditionally, so from `e313109` until now the harness failed on its second
step with "The source directory .../external/json does not exist".

It is now `benchmark/json/`, beside its only caller. It stays a separate CMake
project rather than joining the main one, because it fetches simdjson at
configure time and the shipped library must build with nothing fetched.

The shape is the same as the other five, one level up: the *check* answered a
different question than the one being asked. "Is it in the build?" is not "does
anything use it?", and a grep of `CMakeLists.txt` cannot tell them apart.

## 7. Reading a symptom off the wrong line of the compiler's output

`narrow_phase_vq_kernel` reported, at `-Xptxas -v`, an 8128-byte stack frame,
112 bytes of spill stores, 216 of spill loads and 255 registers — the only
kernel in the build that spilled. Its `Domain stack[140]` was 7840 of those
8128 bytes, so the diagnosis wrote itself: the thread-local stack is the
problem, move it to a block-shared pool with a global overflow queue the way
the triangle kernel does. That was the plan, and it was approved.

It was wrong on every count, and one compile settled it. Rebuilding at stack
caps from 4 to 513 entries:

    cap    frame      registers   spill
      4      448 B    254           0/0
      8      736 B    255       112/216
     16     1184 B    255       112/216
     32     2080 B    255       112/216
    140     8128 B    255       112/216
    513    29024 B    255       112/216

The registers and the spill do not move with the array at all. They are the
inclusion function's working set — thirty coordinates, two `Frame`s, eight
corner triples — and no stack change touches them. The frame is the array, and
the frame is not the spill; they had been read as one number.

Nor does the frame cost time. Only the entries a search actually reaches are
ever written; the rest is address space. On GH200, 400k queries at depth 69:
cap 140 runs 7.589 ms, cap 257 runs 7.441 ms, cap 513 runs 7.355 ms — the
largest is the fastest, which is to say it is noise.

**The near-miss.** A first timing run appeared to confirm the restructure
handsomely: cap 140 at 7.613 ms against cap 32 at 3.681 ms, a 2.07x win for the
smaller stack. It measured the wrong thing. `kMaxDeviceDepth` is *derived* from
the cap, so shrinking the array shortens the search: cap 32 means depth 15, not
depth 69. Holding depth fixed at 15, cap 140 and cap 32 run 3.719 ms and
3.680 ms to a bit-identical answer. This is failure mode 5 again, measuring the
wrong object entirely, and it would have "validated" a rewrite that could not
have delivered what it promised.

What the exercise did produce is the opposite change from the one planned.
Since depth is what costs and headroom for it is free, the stack is now sized
*from* a depth (`SCCD_VQ_MAX_DEPTH`, 128) rather than being an entry count that
a depth is derived from. The old 140 entries capped the device at 69, which
happens to equal the host's default `SCCD_MAX_DEPTH` — so the two agreed by
coincidence, and raising `SCCD_MAX_DEPTH` gave a host that searched deeper and
a device that quietly did not. The clamp is safe, since exhaustion accepts at
the box's `t` lower bound, but it is an accuracy divergence and it is now
reported on stderr instead of being silent.

The rule to carry forward: **when a plan rests on a number, re-derive the
number before executing the plan, and check that the knob you are varying moves
only the thing you think it moves.**

### Then it was built anyway, and measured

The reasoning above is sound about *sizes* and wrong about one conclusion, which
is worth separating. Varying `SCCD_VQ_STACK_CAP` does not change the register
count or the spill — that part holds. But the spill was not caused by the array's
size; it was caused by there being a dynamically-indexed array at all. The
restructured kernel, which keeps one box per thread and puts the rest in a
block-shared pool, reports:

    before   8128 B frame   112/216 spill   255 registers
    after     336 B frame       0/0 spill   216-226 registers, 3.6 KB smem

So the restructure does deliver the static improvement the plan claimed. It was
built in full — a seed kernel, a drain kernel, a shared body, a host grow-and-retry
loop, all mirroring `narrow_phase_dfs_zero_stride_body` — and it is conservative:
`sccd_narrowphase_cuda_test` passes all 20 configurations, and at `tol = 1e-16`
the device vertex-quad row is slightly *better* than before (6.675e-14 against
8.601e-14 of earliness).

It still does not ship, because throughput does not follow. 400k queries on
GH200, best of 5, every configuration returning an identical answer:

    scene / stride          before   order only   restructured
    stationary, stride 1    7.393      4.456          8.181
    stationary, stride 0    0.972      0.608          0.739
    moving,     stride 1   16.367     18.012         34.734
    moving,     stride 0    1.796      1.900          1.501

The restructure wins one cell of four (1.20x) and loses more than half its
throughput in another (0.47x). The loss is not the global queue overflowing:
raising the shared pool from 64 to 512 entries moves moving/stride 1 only from
34.4 ms to 32.1 ms. It is the design itself. In `per_query` mode every query is
uniformly deep, so work sharing has no imbalance to correct, while the block
still pays a `__syncthreads()` per DFS iteration and a thirty-coordinate reload
plus a `prepare()` every time a thread picks up a box belonging to another
query. That is the fourth experiment in this project to land on the same
finding: **the lanes are finished, not waiting.** Best-first ordering, per-query
bounds and 128-way dicing on the triangle kernel all lost for the same reason.

### The ordering lead, and why it also does not ship

The `order only` column above is a separate change and looked like the better
prize: ten lines against the existing kernel — buffer the surviving sub-boxes and
push them in reverse, so the LIFO pop follows the earliest `t` first and tightens
the bound sooner. On the synthetic scenes it is worth 1.60-1.66x on a stationary
quad and costs 5-9% on a moving one.

It was taken to the real workload, which meant getting smesh into the CUDA build
on Alps (it was present but never installed, so its build tree had no
`smeshTargets.cmake`; `cmake --install` to a prefix fixes it) and running
`sccd_refine_scaling` with `SCCD_TOPOLOGY=quad` on the device. Level 4 —
6144 faces, 37.7M vertex-face and 301.8M edge-edge pairs — narrow-phase
milliseconds, interleaved repeats:

    rep      before     order
      1    10521.25  10874.16
      2     9763.41   9420.46
      3     9599.97   9654.04
      4     9651.51   9770.79
    median     9707      9713

Identical within noise, and the reported time of impact is 0.6666653904 in all
eight runs. **The synthetic win does not transfer at all.**

The reason is the query mix. Following the earliest `t` first only pays when
there is a root to find early: it reaches one sooner and the tighter bound then
prunes everything else. The synthetic scene is one crossing vertex per quad, so
every query has a root. The real scene is dominated by candidate pairs with no
root at all, and a box holding no root has to be exhausted whatever order its
children are visited in. Traversal order buys nothing on the boxes that cost the
most.

So the knob is not kept. A tuning parameter measured as worthless on the workload
does not earn a place in the shipped kernel any more than a kernel does.

**Note on determinism.** These runs also make it explicit that the quad narrow
phase is not run-to-run reproducible at `toi_stride == 0`: the shared minimum is
an atomic under a dynamic schedule, so which thread prunes first varies and the
reported time of impact moves in the last few digits between runs of the *same
binary* (0.6666666588 / 0.6666664327 / 0.6666653904 were all observed at
level 2-4 on the host). Every value is at or before the true 2/3, so the
guarantee holds; it is the exact answer that is not stable, and a test that
pins a quad time of impact to more than about 1e-6 will flake.

**A stale record corrected by the same run.** `benchmark/assessment/assessment.csv`
carries `FAILED rc=134` for the hopper quad rows, and `assess.sbatch.sh` used to
explain them as the device narrow phase not existing for QUADSHELL4. It does
exist, and `sccd_refine_scaling` with `SCCD_TOPOLOGY=quad` now completes all five
levels on the device. Those rows predate the kernel and should be re-run.

What survives from the whole exercise is the header extraction
(`sccd_device_dfs_stack.cuh`), which stands on its own, and the stack now being
sized from a depth rather than the other way round.

## 8. Measuring a competitor in a configuration its authors do not ship

The first GH200 comparison (`benchmark/competitors/results/compare-gh200-2026-09-17.csv`)
reported Scalable CCD missing `15,393` of `20,668` curated queries a pass on
armadillo-rollers and a step-level error of `+0.99`, and timed its narrow phase at
`81.8` ns per candidate against SCCD's `21.4`. The miss count is real in the sense
that the library produced it. What the harness said about it was not.

**What the harness got wrong.**

1. `SCALABLE_CCD_TOI_PER_QUERY` was forced on for every build. It is off by
   default, and it changes the algorithm: on, each query keeps its own bound and
   nothing prunes across queries. Every timing and every step-level answer was
   therefore taken from a search doing strictly more work than the shipped one.
2. The harness comments and `results/README.md` described a global prune
   (`min_t >= *toi`) that the per-query build does not compile. The "lifted bound"
   route built on that description changed nothing, and the explanation of why
   the batched call "did no work" was a reading of code that was not running.
3. The per-query accuracy numbers came from one query per kernel launch. That
   sizes the subdivision buffer at two domains, which is the worst case for the
   defect below, so it maximised the miss count.
4. The step loop never flushed its last step, so each scene's final step row was
   dropped, and the cases ran edge-edge before vertex-face, the reverse of
   `cuda::ccd`.

**What is actually wrong in the library.** `CCDBuffer::push` checks `is_full()`
and then increments the tail as two separate operations, from hundreds of device
threads at once. When a level of the search needs more room than the buffer has,
two threads can both pass the check for the last slot; the tail then wraps onto
the head, `is_full()` reads false again, no overflow is flagged, and unprocessed
domains are overwritten. `MemoryHandler::handleNarrowPhase` sizes the buffer at
`2 * n` domains for `n` queries, so short lists of colliding queries overflow
constantly. Exhaustion by overwrite drops boxes, and a dropped box is a missed
root.

Established three ways, on armadillo-rollers:

* **Padding the list fixes it.** The same curated queries followed by a million
  pairs rejected on their first test (`--buffer-probe`) go from `0` of `1299`
  edge-edge contacts at step 107 to all `1299`, none late. The padding changes
  nothing but the buffer size.
* **The answer changes between identical runs.** At a padding of `10,000`, step
  100 edge-edge found `7` of `65` on one run and `65` on the next.
* **A race-free check removes every miss.** A copy of the library that refuses,
  on the host, to launch a level unless capacity exceeds three times the level
  (the most a level can push) finds all `20,668` queries over cases `[0, 120)`,
  identically on repeat. That copy is diagnosis only and is not what the
  comparison reports.

The library as shipped (per-query off) is affected too: over the same cases,
`24` of `63` steps with a contact report a time of impact after the true one, two
of them no contact at all, and `cuda::ccd` itself returns `0.0546875` for step
151 against a root at `0.0140006` -- and `1` on another run.

**Settled.** With the race removed the search is conservative and nothing else is
wrong with it. A first run of the patched copy left `47` late per-query answers,
worst `3.0e-3`, which looked like a second defect; they were the float32 mesh
artefact of section 10. Rebuilt against the double-precision smesh, the same
patched copy finds all `54,664` curated contacts of armadillo's first `400` cases
with zero missed and zero late. The buffer race accounts for the whole of the
library's observed non-conservativeness.

**What the comparison does now.** Each competitor is compared on the question it
answers, in its own table: Scalable CCD on the earliest time of impact per case,
both libraries running their earliest-time-of-impact path from a fresh bound over
identical candidates; ACCD per collision pair, both narrow phases timed on SCCD's
broad-phase candidates and scored on the curated queries at the coordinates their
exact roots were computed for. The per-query build of Scalable CCD is no longer
part of it.

**The lesson.** A competitor's compile-time options are part of what is being
measured. A build flag chosen for the harness's convenience silently turned the
cost comparison into a comparison against a different algorithm, and the
accuracy result -- although it pointed at a real defect -- was reported with a
mechanism nobody had checked against the code that was compiled.

## 9. Timing one thread against seventy-two

The first per-pair comparison reported additive CCD at `332` ns per candidate
against SCCD Tight's `33.1` on armadillo-rollers, and read that as the cost of
conservative advancement. Both numbers were measured, and the conclusion drawn
from them was wrong twice over.

**The harness ran ACCD serially.** `ipc::AdditiveCCD::point_triangle_ccd` and
`edge_edge_ccd` answer one pair at a time and carry no parallelism, so a plain
loop over candidates uses one core while SCCD's narrow phase uses 72. The
toolkit does not use them that way: `Candidates::compute_collision_free_stepsize`
wraps the per-candidate call in a `tbb::parallel_for`, which is the loop the
harness should have had. With it, the same chunk of armadillo gives `9.2` ns per
candidate against SCCD Tight's `45.0` -- so per pair ACCD is about five times
*cheaper*, and the trade against it is tightness, not speed: its median error is
`1e-2` where SCCD Tight's is `1e-6`.

**And it scored a different candidate list.** `bench.exe.cpp` computes in
`double`; the ACCD harness had `using scalar_t = smesh::geom_t`, which is
`float`. Both read the same float32 mesh, but the outward-ULP box rounding is
coarser in float, so SCCD's broad phase produced `225,490` candidates for the
ACCD rows against `225,426` for the SCCD rows of the same case -- a `0.03%`
difference, small enough to look like nothing and large enough to make "the same
list" false.

**The lesson.** A cost comparison has to state, and equalise, the resources both
sides get: threads first of all. A per-pair kernel is not the unit a library is
used in, and timing it as if it were measures the harness's loop rather than the
method. The check that would have caught it is the one that caught the second
fault: assert that the candidate counts agree case by case, and ask what hardware
each number was produced on.

### Withdrawn in turn: "per pair additive CCD is about five times cheaper"

That `9.2` ns per candidate is itself a units error, and the correction above is
wrong in the direction it corrected.

The harness times additive CCD twice. `narrow_ms` is its pass over the
broad-phase candidates, and `query_narrow_ms` its pass over the curated queries,
which for armadillo-rollers are `85,296,282` and `131,441` -- six hundred times
apart. `9.2` is the second divided by the first's denominator. Dividing
`query_narrow_ms` by the candidate count today gives `4.03` ns, the same order,
on a chunk of the scene; dividing `narrow_ms` by it gives `338.18`.

Measured over the whole of armadillo-rollers, both libraries answering one time
of impact per candidate over a candidate list that agrees case by case on all
`6,394` cases of the benchmark, on the same 72 threads:

| | ns per candidate |
|---|---|
| SCCD Tight, `narrow_ms_s1` | **32.07** |
| additive CCD, `narrow_ms` | **338.18** |

So per pair SCCD is about ten times cheaper, and the trade against additive CCD
is not speed against tightness: SCCD is ahead on both, its median earliness
`3.39e-06` against `0.0286`. Additive CCD's cost is also insensitive to the
thread binding -- 131.7, 130.9 and 131.4 ns per candidate unbound, under
`numactl --cpunodebind=0` and under `taskset -c 0-71` -- so the original serial
reading of `332` was a serial loop and nothing else.

Two things follow for the harness. The `queries` column is the candidate count,
so any per-pair figure has to come from a column measured over candidates:
`narrow_ms_s1` for SCCD, `narrow_ms` for additive CCD. And SCCD's
`query_narrow_ms` is not a per-query cost at all -- `8413.8` ms over `131,441`
curated queries is `64` microseconds each, `10.8` ms per case over 781 cases,
which is per-case buffer setup rather than search.

## 10. Explaining an artefact instead of removing it

The competitor comparison reported SCCD Tight as late on `137` of `394`
armadillo-rollers steps, and the write-up explained it: the mesh path reads
float32 coordinates while the dataset's exact roots were computed for the
rational query geometry, so the two describe different geometries and the
comparison is not a conservativeness test. The explanation was correct, and it
was still the wrong thing to publish. Patrick's instruction settled it in one
line: **smesh needs to be compiled with double `geom_t` for these benchmarks.**

Built that way, with the frames converted to match, the same sweep reports zero
late steps on all three scenes in every mode, and the worst signed step error is
negative everywhere. The number that needed a paragraph of explanation simply
went away, and what remains is a measurement that means what it appears to mean.
Scalable CCD's `265` late steps are unchanged by the switch, which is what makes
that finding worth reporting.

Two practical notes for rebuilding it. The float64 explicit instantiations are on
smesh's `sideset` branch; `master` has `f32` only and leaves
`mesh_from_folder<int, double>` undefined inside `libsmesh.a`, so nothing that
reads a mesh links. And on GCC the build needs `SMESH_ENABLE_DEV_MODE=OFF`, plus
`#include <cstddef>` in `src/mesh/geometry/smesh_fff.hpp`, which clang supplies
transitively and GCC does not.

**The lesson.** An artefact that is understood is still an artefact. When the
cheaper move is to remove its cause rather than to document it, document nothing
and remove the cause -- a reader should not have to hold a caveat in mind to read
a table correctly.

## A note on the C ABI

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

## Two levels do not beat one on these scenes

A second binning rule was built, measured, and demoted to
`spikes/src/broadphase_hgrid2.hpp`. It is recorded here because it is the
textbook structure and will be proposed again by anyone who reads the
broad-phase literature.

**The rule.** A box is binned into the one cell holding its centroid, and a cell
is at least as wide as any box in it. Two overlapping boxes then have centroids
at most one cell apart on each axis, so a fixed `3x3` stencil finds every
partner, walking the five cells of index at least a box's own reports each
unordered pair exactly once, and there is no duplicate test at all. That is
strictly more elegant than the shipped minimum-corner rule.

**Why one level cannot have it.** The cell has to hold the *widest* box, so every
box pays the worst case. On armadillo-rollers the widest swept edge spans 11% of
the scene, which caps the grid at `8x9` and puts 36,363 edges into 72 cells --
505 to a cell. Measured host-side against the sweep, the broad phase went from
179 ms to 2,126 ms on armadillo-rollers and from 1,177 ms to 11,483 ms on
cloth-ball, up to 11.9x. The shipped rule's cost per box follows that box's own
size; this one's follows the largest box in the scene.

**Two levels recover most of it and still lose.** Splitting the boxes by size
across a fine grid and a coarse one, with the cut chosen per step by costing
every candidate against a histogram of box sizes, gives a correct broad phase --
identical pair sets to the sweep on every synthetic case, and zero missed
collisions, zero late times of impact and identical false-positive counts on
three scenes. It is 1.3x to 1.5x behind the minimum-corner cell list:

| scene | sweep | cell list | two levels |
|---|---|---|---|
| armadillo-rollers | 586 ms | 372 ms | 480 ms |
| cloth-ball | 1787 ms | 1801 ms | 2395 ms |
| cloth-funnel | 468 ms | 370 ms | 546 ms |

Whole-scene broad-phase totals including preparation, 72 threads on one Grace,
120 cases per scene.

**Why two is the wrong number.** The size ratios on these scenes run from 4.5x to
24x, which is two to five octaves, and one cut has to straddle all of it. The
cost model is left trading the fine level's resolution away to keep the coarse
level's population down, and its choices show it: on armadillo's edges it picks a
`33x45` fine grid holding 31,739 boxes, 21 to a cell, when the cell cap would
have allowed 145,000 cells. It refuses to go finer because every step finer
pushes more boxes onto a coarse level stuck at `8x9`. On cloth-funnel's vertices,
at a 16x ratio, it gives up and makes both levels the same grid -- `163x231`
against `162x230`, half the scene on each.

The structure that would work is a geometric ladder of levels with the level
count taken from the size distribution, which is what the polydisperse-particle
literature does (Ogarko and Luding, *A fast multilevel algorithm for contact
detection of arbitrarily polydisperse objects*). That was not built.

## Binning self pairs at the minimum corner: sound, and cheaper than what ships

A specialisation for the edge-edge broad phase, where one list is queried against
itself: bin each box into the single cell holding its **minimum** corner rather
than into every cell its extent touches. One entry per box, a partner met at most
once, and the shipped duplicate rule -- attribute the pair to the cell holding
the minimum corner of the overlap -- is not needed at all.

This entry has been wrong twice. It first recorded the idea as unsound, then as
sound but too expensive. Both conclusions came from a weaker implementation than
the proposal deserved, and the record is kept in this shape because that is the
failure worth remembering.

**What "forward" means decides soundness.** Reading the box's own footprint
rectangle -- its columns crossed with its rows -- is *not* complete. It loses
exactly the pairs whose minimum-corner cells are incomparable, one box ahead on
the first axis and behind on the second, because then neither footprint holds the
other's bin cell and neither ever reads the other. Componentwise order on cells
is partial and a forward walk needs a total one. Two boxes at a cell width of one:

    A = [0.5, 1.5] x [1.5, 2.5]   cell (0, 1)
    B = [1.2, 2.2] x [0.8, 1.8]   cell (1, 0)

Reading the **row-major linear index** instead is complete, because that order is
total: for overlapping `i` and `j`, `L(cell(j.min))` never passes
`L(cell(i.max))`, so whichever has the smaller `L` holds the other in its band.
Equal `L` means one cell, where the index breaks the tie -- the only index test
still needed anywhere.

**The cost objection was an artefact of not pruning.** Walked naively the band is
whole rows and costs a grid width per row spanned, which is where the second
wrong conclusion came from. It need not be walked naively. Each cell carries the
largest upper bound of the boxes binned in it; a prefix maximum along each row
then identifies the columns holding nothing that reaches back to this box, and a
binary search skips them in one step -- the sweep's cummax trick, applied per row.
A second per-cell bound does the same on the other axis.

Measured by `spikes/src/min_corner_self_probe.probe.cpp` against brute force and
against the shipped full-extent binning on the same boxes:

| boxes | grid | candidates vs shipped | cell entries vs shipped |
|---|---|---|---|
| 4,000, uniform | 34x34 | 0.58x | 0.26x |
| 4,000, uniform, small boxes | 204x200 | 0.51x | 0.25x |
| 20,000, uniform | 67x67 | 0.60x | 0.25x |
| 4,000, 1% of boxes at 20x the mean | 25x24 | 0.66x | 0.21x |
| 4,000, 5% at 20x | 12x12 | 0.60x | 0.24x |
| 4,000, 1% at 60x | 47x46 | **0.21x** | **0.11x** |

Nothing is missed on any shape. It tests a half to a fifth of the candidates and
holds a quarter to a ninth of the cell entries, and a **heavy tail favours it
more, not less**: a wide box costs the shipped binning a place in every cell it
crosses and costs this scheme one entry, which is exactly the size distribution
the real scenes have.

**Still unverified.** The probe is two-dimensional, where the shipped broad phase
runs a 2D grid over 3D boxes and culls the third axis inside the cell, so the
per-candidate cost differs. The inputs are uniformly scattered, where the real
scenes are surfaces, and a prefix maximum is a statistic a clustered distribution
can defeat. Neither the parallel binning nor the device port has been thought
through. What the probe establishes is that the idea is complete and that the
cost objection recorded here twice does not hold; it does not establish a
speed-up on a real scene.

## The per-row search earns its place, and the own-row search never did

Two questions about the edge-edge walk, both settled by measurement rather than
argument: what the per-row binary search buys, and whether the walk can be held
to the querying box's own footprint.

**The footprint bound is incomplete, and the loss is large.** Holding the shipped
kernel to `c0 = c0b` on later rows drops pairs on 10 of 13 self cases in
`sccd_broadphase_cell2d_test`: 3,698 of 23,756 on "many cells", 882 of 6,971 on
a heavy-tailed spread, 7,460 of 113,500 on "dense". Two overlapping boxes can be
incomparable -- one binned a column to the right and a row below the other -- so
each footprint covers one of the other's coordinates and misses the other.
`tikz/incomparable.tex` draws the case to scale.

**The search removes about 25 column visits for every one it leaves.**
Instrumenting the kernel with counters and running six real cloth-ball cases on
an M1:

```
boxes=1,388,250  rows/box=1.99  multirow=80.2%  searches=0.99/box
cols_skipped=103,667,394  cols_visited=4,288,912  ratio=24.2  cells_read=5,568,886
```

A cloth-ball edge box spans almost exactly two rows, so nearly every box runs one
search, and each jumps 75 columns to leave 3.

**In time, removing it costs 3x to 7x.** One Grace at 72 threads, mode 2,
edge-edge queries only, `broad_ms` summed over the cases; `old` is the guard
inside the inner loop, `new` is the own cell read first, `nosearch` is `new` with
the per-row search replaced by a scan from column 0.

| scene | strategy | cases | old | new | nosearch | new/old | new/nosearch |
|---|---|---|---|---|---|---|---|
| cloth-ball | Cell2DMin | 86 | 822.1 ms | 782.8 ms | 2,336.3 ms | 1.050x | 2.98x |
| cloth-ball | Cell2DMinSort | 86 | 601.8 ms | 603.9 ms | 1,166.7 ms | 0.997x | 1.93x |
| rod-twist | Cell2DMin | 218 | 1,474.6 ms | 1,359.0 ms | 9,368.8 ms | 1.085x | 6.89x |
| rod-twist | Cell2DMinSort | 218 | 1,296.3 ms | 1,158.0 ms | 4,832.1 ms | 1.119x | 4.17x |

Per frame, the maxima move the same way: cloth-ball Cell2DMin 16.500 ms, 15.496
ms and 46.432 ms; rod-twist Cell2DMin 11.007 ms, 10.416 ms and 61.180 ms.

**The own-row search, by contrast, could never fire.** The row's prefix maximum
at column `c0b` includes the querying box's own maximum on that axis, so the
search there can never return a later column than the box's own cell. Asserting
that inside the kernel over the whole test suite found no case where it did, so
the own row now starts at `c0b + 1` and pays no search. That, with the own-cell
guard hoisted out of the inner loop, is the `new` column above.

## Face-vertex: binning the faces and querying with the vertex segment

The vertex side of a face-vertex query is the one primitive whose swept volume is
not a box but a segment: vertices move affinely, so a vertex over a step is the
line from its position at $t=0$ to its position at $t=1$. The proposal is to
build the cell list over the faces and query it with that segment, walking only
the cells the line crosses.

**It is complete.** A face-vertex contact needs the vertex to lie inside the
moving triangle at some $t$, so its position at that $t$ lies in the face's swept
box, so the segment meets the face's swept box, so some cell the segment crosses
holds that face -- provided the faces are binned by extent, into every cell their
projected box touches. The filter is also strictly tighter than the shipped one:
segment-meets-box implies box-meets-box and not the reverse.

**The duplicate rule has to change.** A face met along the segment may be met in
several cells, so the overlap-corner rule of `alg:cell` does not apply. The
analogue is: compute the parameter at which the segment enters the face's box and
emit only from the cell holding the segment at that parameter. Same shape, still
$O(1)$, still stateless.

**Measured on the shipped data, the gain is small.** Instrumenting the host cell
list to count both designs on the same frames (six to eight steps per scene, M1):

| scene | vertex box cells | segment cells | box/seg | face box cells |
|---|---|---|---|---|
| cloth-ball | 2.06 | 1.86 | 1.11x | 3.95 |
| armadillo-rollers | 1.46 | 1.41 | 1.03x | 4.25 |
| cloth-funnel | 1.35 | 1.30 | 1.04x | 5.05 |

A vertex crosses barely more than one cell, because the grid is sized to the mean
box extent and a step moves a vertex less than that: only 25% to 65% of vertices
leave their own cell at all. The segment is therefore hardly tighter than the box
it replaces.

Cell *visits* look dramatic and are misleading. Today a face box is walked over a
grid sized by vertex displacement, so it covers 8.4 cells on cloth-ball and 26.4
on armadillo-rollers, against 1.9 and 1.4 for a segment on the face grid -- 8x to
38x fewer visits. Occupancy cancels almost all of it, because a face-sized cell
holds around forty faces where a displacement-sized cell holds one vertex. On
candidates actually examined, which is the work:

| scene | today | faces binned, vertex box query | faces binned, vertex segment query |
|---|---|---|---|
| cloth-ball | 27,108,652 | 26,456,011 (1.02x) | 23,220,596 (1.17x) |
| armadillo-rollers | 2,375,127 | 2,128,810 (1.12x) | 2,059,618 (1.15x) |
| cloth-funnel | 2,008,912 | 2,027,987 (0.99x) | 1,887,937 (1.06x) |

So 6% to 17% fewer candidate tests, against a cell array that grows from 2.1
entries per vertex to about 4.3 entries per face over twice as many faces -- four
times the binning work and memory. On these scenes the build cost plausibly eats
the query gain.

**Where it would win.** The whole argument turns on displacement against cell
size. A step that moves a vertex several cells makes its box quadratically worse
than its segment, $(d_{0}+1)(d_{1}+1)$ cells against $1 + d_{0} + d_{1}$, and the
gain grows without bound. Our scenes do not do this; a solver taking larger steps
would.

**The measurement not yet made** is the one most likely to matter: the pairs
*emitted*, not the candidates examined. Segment-against-box is a strictly tighter
predicate than box-against-box, so it hands the narrow phase fewer pairs while
staying conservative, and the narrow phase is the expensive side. The counters
above measure the broad phase's own work only.

### Measured, after the query was made to pay what a box query pays

The first run of `Cell2DSeg` lost everywhere, and the cause was the inner loop.
The slab test divided by each component of the step per candidate face, and the
cell's boxes were gathered through `cellidx`, which no SIMD kernel can see as
lanes. Hoisting the reciprocals to the vertex, packing each cell's boxes into
cell order so the query runs the sweep's 32-at-a-time kernel, and putting the
segment's own hull in front as an exact pre-pass turned the query around: the
face-vertex query went from `0.56x`-`0.99x` to `1.09x`-`2.19x` against the box
query. Output is unchanged, and no scene misses a collision.

Counting the structure build, which is where this design pays for itself, the
result is scene-dependent. One Grace at 72 threads, mode 2, face-vertex rows
only, paired by case:

| scene | broad (with binning) | narrow | total | pairs |
|---|---|---|---|---|
| armadillo-rollers | 165.0 -> 171.2 ms (0.964x) | 1.061x | **0.984x** | 23.3% fewer |
| cloth-ball | 565.9 -> 361.8 ms (1.564x) | 0.988x | **1.391x** | 16.9% fewer |
| cloth-funnel | 91.9 -> 104.9 ms (0.876x) | 0.838x | **0.869x** | 16.5% fewer |
| n-body-simulation | 2645.7 -> 2157.8 ms (1.226x) | 1.159x | **1.203x** | 25.8% fewer |
| rod-twist | 206.5 -> 298.2 ms (0.692x) | 1.050x | **0.798x** | 1.6% fewer |
| all five | 5331.8 -> 4568.7 ms | | **1.167x** | |

The split is by how much query there is to save. Binning the faces and packing
their boxes costs 15 to 111 ms more per scene than binning the vertices, and the
query returns 2 to 558 ms of it:

| scene | extra preparation | query saved | net |
|---|---|---|---|
| cloth-ball | +59.5 ms | -263.6 ms | **-202.3 ms** |
| n-body-simulation | +70.5 ms | -558.4 ms | **-666.9 ms** |
| armadillo-rollers | +44.6 ms | -38.4 ms | +3.5 ms |
| cloth-funnel | +15.3 ms | -2.3 ms | +17.2 ms |
| rod-twist | +110.7 ms | -19.0 ms | +85.4 ms |

So it wins on the two scenes whose face-vertex query is large and loses on the
three where the preparation is most of the cost. The aggregate favours it, and
the aggregate is carried by n-body-simulation.

The narrow phase does not repay the pair reduction. cloth-ball sheds 16.9% of
pairs for `0.988x`, cloth-funnel 16.5% for `0.838x`. A tighter broad phase
removes pairs that are far apart, which are the ones the narrow phase rejects in
its first box test; the expensive pairs are near-contacts and the segment filter
keeps every one. cloth-funnel being slower suggests the pair order also matters,
since the pairs now arrive grouped by vertex where they used to arrive grouped by
face.

**The obvious next move** is to fold the packing into `cell2d_fill`, which
already scatters one value per cell entry and could scatter the six box
coordinates in the same pass. That removes a whole gather over the cell array
from the preparation, and it is the larger half of the penalty on rod-twist.

**Second, the packing is not specific to this query.** The shipped `cell2d`
face-vertex and edge-edge queries gather through `cellidx` too, so they cannot
vectorise their inner loops either. If the packed layout pays here it should be
tried under them.

## The cell list gets a layout its queries can vectorise

Sweep and prune was the only broad phase here with a vectorised inner loop,
because it sorts its boxes and its window is therefore contiguous. The cell list
stores indices, so a query reached a box by gathering six coordinates through
`cellidx` and no SIMD kernel could see consecutive lanes. `cell2d_pack_boxes`
writes each cell's boxes out in cell order once per step, after which every query
tests $32$ of them at a time with `vaabb_overlap_one_to_many_bits` and visits the
survivors by a bit scan.

**Fusing the packing into the scatter is slower, and by a lot.** Writing the box
inside `cell2d_fill`, where the scatter already holds it and already knows where
the entry lands, costs six scattered stores per entry, one per array. A separate
pass costs six sequential stores and one gather load. Isolated on the face-vertex
query, where only the packing differs (ms, prep + broad): armadillo-rollers
$130.2 \to 165.5$, cloth-ball $391.8 \to 483.4$, cloth-funnel $111.3 \to 151.4$,
n-body-simulation $1921.2 \to 2045.9$, rod-twist $270.7 \to 372.3$. The comment
on `cell2d_pack_boxes` records it.

**It pays according to how many entries the list holds per box.** One Grace at
72 threads, mode 2, whole broad phase, every figure a median over three repeats
of the case, from `benchmark/assessment/bp-packing.csv`:

| scene | Cell2DMin ee | Cell2DMin vf | Cell2D ee | Cell2D vf |
|---|---|---|---|---|
| armadillo-rollers | 1.007x | 1.069x | 1.018x | 0.925x |
| cloth-ball | 1.292x | 1.461x | 1.197x | 1.254x |
| cloth-funnel | 1.011x | 0.949x | 0.960x | 0.879x |
| n-body-simulation | 1.252x | 1.296x | 0.995x | 1.209x |
| puffer-ball | 1.335x | -- | 0.992x | -- |
| rod-twist | 0.837x | 0.836x | 1.122x | 0.792x |
| all scenes | **1.310x** | **1.265x** | **0.998x** | **1.151x** |

It earns its place under the shipped default and nowhere else: `Cell2DMin` gains
$1.31\times$ on edge-edge and $1.27\times$ on face-vertex, while the
extent-binned `cell2d` is level on edge-edge. The minimum-corner list holds one
entry per box where extent binning holds two to four, so the layout costs a
quarter to a half as much for the same query. rod-twist loses on three of the
four columns, and it is the scene with the least query per entry.

**Every earlier figure in this section was a single sample.** The scripts that
produced them indexed rows by case, and with more than one repeat that keeps
whichever repeat came last rather than the median. Two runs of the same binary
disagreed by $1.5\times$ on puffer-ball, which is how it surfaced. The table
above is a median of three; `benchmark/assessment/bp-packing-report.py` is the
corrected analysis and says so in its header.

## Sorting the cells wins only once the structure is counted with it

The six-scene run reported `Cell2DMinSort` ahead of `Cell2DMin` on all six
scenes, $1.098\times$ to $1.350\times$. That was `broad_ms` alone. Sorting each
cell is a pass over the cell array and belongs to the broad phase like any other
part of building its structure, and re-read with `prep_ms` included the same run
says:

| scene | broad only | prep + broad |
|---|---|---|
| armadillo-rollers | 1.243x | 0.991x |
| cloth-ball | 1.292x | 1.121x |
| cloth-funnel | 1.350x | 0.925x |
| n-body-simulation | 1.098x | 0.850x |
| puffer-ball | 1.203x | 1.170x |
| rod-twist | 1.278x | 1.222x |
| all scenes | 1.199x | **1.139x** |

Still ahead in aggregate on that reading, but on three scenes behind, and the
margin is a third of what the query-only figure claimed. Both columns are single
samples per case and neither is to be trusted.

Measured again through the report pipeline, which takes the median over repeats,
the sorted variant loses outright. Over the six scenes `cell2dmin` costs
$16.8\,\mathrm{s}$ against `cell2dminsort`'s $20.0$, and edge-edge alone
$9.4\,\mathrm{s}$ against $11.9$. `tab:bpvariant` reports it and the host runs
`cell2dmin`.

## The minimum-corner walk on the device

`cell2dmin` is the host default, and the device refused it, so the two
processors could not run the same strategy. The port is
`cell2dmin_setup_and_count`, `cell2dmin_fill`, `cell2dmin_bounds` and the two
self-overlap entry points in `src/cuda/sccd_cell2d_broadphase.cu`.

**It is correct.** On GH200 over four scenes it emits exactly the edge-edge pair
counts the extent-binned device kernel does -- 5,692,480 on armadillo-rollers,
259,683,882 on cloth-ball, 3,636,886 on cloth-funnel and 54,083,436 on rod-twist
-- with no missed collisions, and the host kernel agrees with both.

**A thread to a box does not pay there; a warp to a box does.** Edge-edge,
prep + broad, median over repeats, against the extent-binned cell list on the
same processor:

| scene | host | device, thread per box | device, warp per box |
|---|---|---|---|
| armadillo-rollers | 1.436x | 1.523x | **6.243x** |
| cloth-ball | 2.589x | 0.963x | **2.798x** |
| cloth-funnel | 1.702x | 0.673x | **3.464x** |
| rod-twist | 1.557x | 1.710x | **2.706x** |
| all four | **2.205x** | 0.996x | **3.086x** |

The thread-per-box form diverges on everything: the number of rows a box spans,
where its binary search lands, and how many entries each cell holds all differ
between neighbouring boxes, so the lanes of a warp are mostly masked off and
each lane reads a different part of the cell array. A warp to a box inverts it.
Every lane computes the same walk, redundantly and branch-free, and the lanes
split each cell's entries between them -- the only loop left that can diverge,
and the one whose loads are then consecutive. Redundant arithmetic on 32 lanes
buys regular control flow and coalesced reads, which is the right trade on a GPU
and the wrong one on a host: the two kernels differ in shape now rather than
being transliterations of each other.

Two details carry the correctness. The early return is on `fi = tid / 32`, which
is uniform over the warp, so a warp leaves together and the ballots are taken
converged; a per-lane return would break them. And the entry loop strides a
shared trip count rather than starting each lane at `begin + lane`, so every
lane reaches the ballot even when a cell's entry count is not a multiple of 32.
Counting reduces with a ballot popcount and collecting writes at each lane's
rank within the ballot, so neither pass needs an atomic.

It was the walk itself, and the warp-per-box form above is the answer.

## Confirming the device default with one binary

The competitor comparison of 28 September and the one of 21 September were
placed side by side to see what the minimum-corner list costs the device, and
they disagreed with everything else: puffer-ball read 21.9 s against 10.4 for
the extent-binned list and n-body 4.5 against 1.9, while armadillo-rollers,
cloth-funnel and rod-twist went the other way. Two runs a week apart, two
builds, so the comparison was between builds as much as between binnings.

Measuring both with the same binary answers it. Device, mode 2, prep + broad
over the first 100 cases of each scene, median over two repeats
(`benchmark/assessment/device-binning.csv`):

| scene | cell2dmin | cell2d | ratio |
|---|---|---|---|
| armadillo-rollers | 148.7 ms | 400.5 | **0.37x** |
| cloth-ball | 690.2 | 977.7 | **0.71x** |
| cloth-funnel | 295.4 | 525.9 | **0.56x** |
| n-body-simulation | 2670.2 | 4323.3 | **0.62x** |
| puffer-ball | 8922.4 | 10951.7 | **0.81x** |
| rod-twist | 193.8 | 411.3 | **0.47x** |
| all six | **12920.7** | 17590.3 | **0.73x** |

The minimum-corner list leads on every scene of the device, so `Cell2DMin` as
the single default holds on both processors and nothing measured with it needs
revising.

What the cross-run comparison actually saw is a device regression in the
extent-binned kernel between the two builds: n-body cost it 12.7 ms per case in
the 21 September build and 43.2 ms here. `cell2d` is a shipped option rather
than the default, so this does not affect any reported number, and the cause is
not yet identified.

A comparison taken from two runs of different builds measures the difference
between the builds. Only the same binary answers a question about two
strategies.

## Vertex-face: binning the faces at their minimum corner

The edge-edge query was moved to minimum-corner binning and left the vertex-face
one alone, so vertex-face is now the larger half of the broad phase: over six
scenes at `tight`, host 18,807.7 ms against edge-edge's 16,097.5 (53.9%), device
11,323.6 against 17,829.4 (38.8%), reaching 74.6% on n-body host and 82.6% on
armadillo device. It is a third to an eighth of the queries and 29% to 83% of the
time, so it costs two to seven times more per emitted pair.

What ships bins the **vertices** by extent and iterates over the **faces**.
`cell2d_setup` sizes cells from the list it bins, so the cells are vertex-sized
and the larger element walks a rectangle of them, paying the overlap-corner
duplicate rule on every surviving candidate. Binning the vertices at their
minimum corner instead buys little: they already occupy about two cells each, and
a prefix maximum over small boxes rises steadily, so the search would skip
almost nothing. The lever is which list is indexed, not the binning rule.

`Cell2DSeg` already established that the mirror image has the cheaper query, and
lost on the structure: extent-binning faces costs about 4.3 entries each over
twice as many elements. **Minimum-corner binning removes exactly that cost**, one
entry per face, which is why this was measured again.

Completeness differs from the edge-edge walk, which leans on a self-symmetry that
two lists do not have. A face overlaps a vertex only if `f.min <= v.max`
componentwise, so `cell(f.min) <= cell(v.max)` in both axes and nothing right of
or above the vertex's maximum cell is ever read. Leftwards, the per-row prefix
maximum gives the first column holding anything that reaches back, by the same
binary search the edge-edge walk uses. Downwards, a face reaching the vertex on
the second axis has `row(f.max) >= row(v.min)` and spans at most `K1` rows, so it
is binned no lower than `row(v.min) - K1`. Anchoring that bound on `row(v.max)`
instead loses pairs whenever the vertex straddles a boundary -- 81 of 5,543 in
the first run of the probe.

`spikes/src/cellmin_vf_probe.probe.cpp` counts both schemes on the same boxes,
with each grid sized from the list it bins and the `4n` cell cap applied, which
brings the occupancies to 2.44 and 4.49 cells against the 2.06 and 3.95 measured
on the scenes. Proposed over shipped:

| shape | entries | cells read | candidates |
|---|---|---|---|
| base | 0.82x | **0.23x** | 0.57x |
| + 8 faces at 20x the mean | 0.82x | 0.25x | 0.62x |
| + 80 faces at 20x | 0.82x | 0.30x | 0.83x |
| five times the elements | 0.59x | 0.28x | 0.68x |
| five times, + 200 at 20x | 0.59x | 0.34x | 0.96x |
| two vertices per face | 0.19x | 0.33x | 0.73x |

Neither scheme misses a pair on any of eight shapes. Cells read fall three to
four times over, and that is the counter the case rests on: the candidate count
degrades towards parity under a heavy tail (0.96x with 200 oversized faces)
because a wide face reaches into many vertices however it is indexed.

The cost the table does not show is the binary search per row scanned. `K1` is a
maximum, so one face twenty times the mean takes the row range from 3 to 21, and
every vertex then pays 21 searches whether or not a row yields a cell. The
skipped-column counts say those rows are almost entirely empty (778,094 columns
skipped without the tail, 4,423,950 with it), so the searches are doing their
job -- but an implementation should count them, and a per-row upper bound on the
second axis would let a row be rejected before its search rather than after.

This clears the bar to implement and measure. It is not yet a result about time.

### What it measured

`Cell2DMinFV` against `Cell2DMin` from one binary, host, 60 cases per scene,
`tight`, three repeats, on one Grace at 72 threads under
`numactl --cpunodebind=0 --membind=0`, median over repeats per case and summed
over cases. Both keep the same edge-edge query, so the vertex-face columns carry
the whole difference. Taken after the subsample was fixed to spread over each
query type, which is what put puffer-ball's vertex-face cases in the run at all.

| vertex-face, ms | prep | | | query | | | net | |
|---|---|---|---|---|---|---|---|---|
| | min | minfv | delta | min | minfv | delta | | |
| cloth-funnel | 24.3 | 26.5 | +2.2 | 7.7 | 5.7 | -2.0 | +0.2 | 0.99x |
| armadillo-rollers | 32.3 | 25.9 | -6.3 | 35.5 | 9.9 | -25.6 | **-31.9** | 1.89x |
| cloth-ball | 37.4 | 35.0 | -2.4 | 199.7 | 68.3 | -131.4 | **-133.8** | 2.30x |
| rod-twist | 38.9 | 43.7 | +4.7 | 48.0 | 28.3 | -19.7 | **-15.0** | 1.21x |
| n-body-simulation | 43.7 | 41.3 | -2.4 | 923.1 | 147.0 | -776.1 | **-778.5** | 5.13x |
| puffer-ball | 468.9 | 286.1 | -182.8 | 2951.4 | 993.4 | -1958.0 | **-2140.8** | 2.67x |
| all six | 645.5 | 458.5 | **-187.0** | 4165.4 | 1252.6 | **-2912.8** | **-3099.8** | **2.81x** |

The prediction the plan rested on is the prep column. Extent-binning the faces
cost `Cell2DSeg` between +15.3 ms and +110.7 ms of prep per scene, which is what
sank it; one entry per face costs nothing anywhere and pays 182.8 ms back on
puffer-ball, where the shipped scheme's 468.9 ms of vertex binning is a sixth of
the whole broad phase. The query falls 3.33x on top of that.

Cloth-funnel is a wash at +0.2 ms: its vertex-face query is 7.7 ms to begin with,
too small for the structure to pay for itself, and too small to matter.

Attributing only the vertex-face delta, the broad phase over the six scenes goes
from 8219.9 ms to 5120.0 ms, **1.61x**:

| BP full, ms | min | minfv | |
|---|---|---|---|
| cloth-funnel | 88.5 | 88.8 | 1.00x |
| armadillo-rollers | 118.1 | 86.2 | 1.37x |
| cloth-ball | 429.6 | 295.8 | 1.45x |
| rod-twist | 249.2 | 234.3 | 1.06x |
| n-body-simulation | 1338.6 | 560.1 | 2.39x |
| puffer-ball | 5995.8 | 3854.9 | 1.56x |
| all six | **8219.9** | **5120.0** | **1.61x** |

The edge-edge halves are held at the measured `cell2dmin` figure in that table
rather than carried across from the other run, because both strategies execute
identical edge-edge code and its columns move several percent between runs. That
spread is run-to-run noise and is not read as a result in either direction.

### Correctness

Across 2160 rows: no missed pair, no missed collision, no time of impact later
than the dataset's exact root, under either strategy. False-positive counts are
identical per scene (3, 6, 3, 60, 6, 372), and the candidate count agrees for
every scene and query type, so the two emit the same pair set. Brute force checks
the same kernel on every case of
`src/tests/sccd_broadphase_cell2d_test.exe.cpp`, quads included.

The reported earliest time of impact differs in the last digits on 292 of 1080
cases, at most 2.6e-05 absolute, and symmetrically -- 141 earlier and 151 later.
Changing the order in which pairs reach the shared atomic minimum sharpens it
differently; neither run is late against ground truth, which is the property that
has to hold.

Measured with the double-`geom_t` smesh at `$SCRATCH/installations/smesh-f64`. A
float build reports six late answers on armadillo under **both** strategies, so
local runs settle A/B questions only and say nothing about conservativeness.

### 11. Sampling the query type along with the trajectory

The even-spread subsample took `cases[k * cases.size() / MAX_CASES]` from a list
sorted by case key. Keys end in `ee` or `vf`, so a scene holding equal numbers of
each sorts into interleaved pairs with edge-edge on every even index. Puffer-ball
has exactly 240 such cases, so `--max-cases 60` strides by 4, every index it
picks is even, and the subsample is entirely edge-edge. Index 0 is `100ee` and
index 4 is `102ee`, which is what the result CSV's first two rows are.

The other five scenes escaped by accident, their counts being unequal enough for
the stride to drift across both types -- cloth-ball 43 and 36, armadillo 396 and
385.

Nothing in the output said so. Puffer-ball appeared in the table as a row of
zeros, which reads as a scene that issues no vertex-face queries and was written
down as exactly that. The dataset holds 130 converted vertex-face box directories
and 120 vertex-face query CSVs for it.

The subsample now spreads over each query type separately and keeps the two in
the proportion the dataset has them, so a balanced case list yields 30 and 30
where it yielded 60 and 0. A type the dataset holds is always represented, and
neither share exceeds what that type has.

What it cost to miss: puffer-ball is the largest scene in the set, and on the
rerun it is the largest absolute win of the change, -2140.8 ms, with the biggest
prep saving of any scene. The table above was published without it and read
1.67x where it now reads 1.61x on a set twice the size -- the direction held, but
the scene carrying the most weight had been dropped silently.

The generalisation is the one this document already makes five times over: a
measurement that answers a different question than the one asked. The particular
form is worth naming on its own, because the tell was visible in the result file
and went unread -- a count of rows by query type shows every scene with a mix and
one scene without, which is not a shape a scene produces.

## The minimum-corner vertex-face query on the device

A warp to a vertex, the shape `warp_forward_self_partner` already uses for the
edge-edge walk and for the same reason: the rows a box spans, the columns its
search leaves and the entries in each cell all differ between neighbouring
boxes, so every lane computes the same walk and the lanes split each cell's
entries between them. The cross-list form drops the own cell and the index test,
because the lists are distinct and a face has exactly one cell, so a pair is met
exactly once. In exchange every row pays a search, where the self walk skips it
on its own row.

**It is correct.** The device test now runs both vertex-face forms on the same
boxes and they return the same pair set on all nine shapes -- triangles and
quads, both lopsided list ratios, and the degenerate spreads. Over the six
scenes, 2160 rows: no missed pair, no missed collision, no time of impact later
than the dataset's exact roots, false positives identical per scene
(3, 6, 3, 60, 9, 387) and candidate counts equal everywhere.

Six scenes, 60 cases each, `tight`, three repeats, medians, one binary:

| vertex-face, ms | prep min | prep fv | delta | query min | query fv | delta | net | |
|---|---|---|---|---|---|---|---|---|
| cloth-funnel | 6.5 | 7.3 | +0.8 | 91.8 | 7.3 | -84.5 | **-83.7** | 6.75x |
| armadillo-rollers | 6.9 | 6.3 | -0.6 | 69.9 | 8.4 | -61.5 | **-62.1** | 5.22x |
| cloth-ball | 7.9 | 9.3 | +1.3 | 322.2 | 55.2 | -267.0 | **-265.6** | 5.12x |
| rod-twist | 9.3 | 7.2 | -2.1 | 30.3 | 16.9 | -13.3 | **-15.4** | 1.64x |
| n-body-simulation | 16.0 | 14.1 | -1.9 | 694.7 | 349.9 | -344.9 | **-346.8** | 1.95x |
| puffer-ball | 65.3 | 68.1 | +2.8 | 1418.8 | 1176.1 | -242.7 | **-239.9** | 1.19x |
| all six | 111.9 | 112.3 | **+0.4** | 2627.7 | 1613.7 | **-1013.9** | **-1013.6** | **1.59x** |

The broad phase over the six goes from 6869.5 ms to 5855.9 ms, **1.17x**, ahead
on every scene: 2.81x on armadillo, 2.33x on cloth-ball, 2.14x on cloth-funnel,
1.27x on n-body, 1.18x on rod-twist and 1.06x on puffer-ball.

### The prep column went the other way, and it was the grid sizing

The table above is the corrected one. On the first device run prep came out
**dearer** by 87.7 ms over six scenes, 79.7 of it puffer-ball, where the host
measurement of the same structure had it *cheaper* by 187.0. Two measurements of
the same thing disagreeing in sign meant one of them was being charged for
something the other was not.

Two guesses, both wrong, both worth recording as such:

1. **The row span round trip.** The walk took it as a launch argument, so every
   step fetched it back with a `cudaMemcpy` that synchronises the device. Keeping
   it in device memory moved the gap from +97.7 ms to +87.7 -- about a tenth. The
   change is kept, because a host round trip for a device-computed value is wrong
   whatever it costs, but it was not the explanation.
2. **`min_row_prefix_kernel`.** The extent path never runs it, it is one thread
   to a row scanning serially, and the face grid is the larger of the two, so it
   looked like the answer. It is **0.78 ms**, 0.0% of GPU time.

`nsys` on the shipped binary settled in one run what two readings of the source
had not. Puffer-ball, 12 cases, per-kernel, shipped against minimum-corner:

| prep kernel | cell2dmin | cell2dminfv | delta |
|---|---|---|---|
| `grid_stats_kernel` | 152.12 ms | 191.18 ms | **+39.06** |
| `min_cell_bounds_kernel` | 1.94 | 3.51 | +1.57 |
| all the binning, prefix and span kernels together | 2.42 | 2.93 | +0.51 |

**94% of the gap, and 10.1% of all GPU time, in the kernel that sizes the grid.**
Nothing about minimum-corner binning is expensive. `grid_stats_kernel` launched
`<<<1, 256>>>` -- one block reducing the whole array on a device with 132 SMs, so
its cost is linear in the element count with none of the machine working. The
minimum-corner path sizes its grid from the faces, which is the longer list, and
the kernel charges by length.

Split into a grid-stride pass over 512 blocks writing per-block partials and a
single-block pass combining them, it falls from 152.12 ms to **1.33 ms**, and the
difference between the two strategies from +39.06 ms to **+0.17**. The kernel's
own note had said its cost was "not worth optimising", which was true when it was
a fixed small cost and stopped being true once it was a tenth of the device.

It helps both strategies, so the incumbent gains too: `cell2dmin`'s broad phase
over the six scenes falls from 7563.9 ms to 6869.5. The vertex-face prep gap is
now +0.4 ms over six scenes, which is the host's answer as well -- one entry per
face costs nothing, on either processor.
