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
- Withdrawn: "SCCD is 6x to 24x cheaper per collision pair than Additive CCD" —
  see section 9
- Demoted: a centroid-binned cell list, at one level and at two — see below
- Corrected: "binning self pairs at the minimum corner loses pairs" — it is
  sound over a total order on cells; only the footprint walk fails — see below

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

## Binning self pairs at the minimum corner: sound, and it costs a grid width

A specialisation for the edge-edge broad phase, where one list is queried
against itself: bin each box into the single cell holding its **minimum** corner
rather than into every cell its extent touches, and read only cells that come
after its own. Each box is then in one cell, so a partner is met at most once,
the minimum-corner duplicate test disappears, and the cell array shrinks from one
entry per covered cell to one per box.

**It is sound, and everything turns on what "after" means.** The first reading
tried here -- read the cells of the box's own footprint, its columns crossed with
its rows -- is not complete, and that was recorded here as a refutation of the
whole idea. That was too strong, and the correction is the point of this entry.

The footprint walk loses exactly the pairs whose minimum-corner cells are
**incomparable**: one box ahead on the first axis and behind on the second. A
concrete one, at a cell width of one:

    A = [0.5, 1.5] x [1.5, 2.5]   min cell (0, 1)
    B = [1.2, 2.2] x [0.8, 1.8]   min cell (1, 0)

They overlap on both axes, B sits one cell left of A and one cell above it, and
neither footprint contains the other's bin cell, so neither ever reads the other.
No index test is involved -- the inner loop never yields the other index.
Componentwise ordering on cells is only **partial**, and a forward walk needs a
total one.

**The row-major linear index is a total order, and over it the walk is
complete.** Read the cells whose linear index lies between the box's own cell and
the cell of its maximum corner. For overlapping `i` and `j`, `L(cell(j.min))` is
never past `L(cell(i.max))`: if the rows differ then the row term settles it, and
if they are equal then `j.min` not past `i.max` on the first axis settles the
column. So whichever of the two has the smaller `L` holds the other inside its
band. Equal `L` means the same cell, where the index breaks the tie -- which is
the one place an index test is still needed, and the only one.

`spikes/src/min_corner_self_probe.probe.cpp` runs both walks against brute force;
build with `-DSCCD_ENABLE_SPIKES=ON` and run
`min_corner_self_probe [boxes] [box extent]`. The band walk misses nothing on any
shape tried, from a 5x5 grid to 204x200 and box extents from 1 to 40. The
footprint walk misses 4,488 of 26,950 pairs on the default shape.

**What keeps it out of the shipped code is the cost.** "Forward" in a total order
on cells means whole rows, so a box spanning `R` rows reads about `R` times the
grid width rather than its own footprint. Candidates tested, band against the
shipped full-extent binning on the same boxes:

| boxes | grid | shipped | band |
|---|---|---|---|
| 4,000, extent 6 | 34x34 | 107,926 | 488,018 |
| 4,000, extent 1 | 204x200 | 3,081 | 80,049 |
| 1,000, extent 40 | 5x5 | 210,023 | 207,314 |

The band is competitive only where the grid is narrow enough that a row costs
little -- the last row, where it draws level. On a grid fine enough to be worth
having it is four to twenty-six times worse, because each row the box spans drags
in two hundred columns it has no geometric reason to read. The shipped rule buys
its completeness from the geometry rather than from an ordering, and pays two
clamps per surviving pair for it.
