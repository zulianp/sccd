# Assessment

Measurements behind the keep-or-demote calls in the repository cleanup. Every
claim here cites a row of `benchmark/assessment/assessment.csv`; regenerate the
tables with `python3 benchmark/assess_report.py benchmark/assessment/assessment.csv`.

**Run.** Alps job 4592444, one GH200 node, 6m10s. Grace at
`OMP_NUM_THREADS=72`, Hopper (sm_90) on the same node so the whole sweep sits in
one allocation. Three repeats, variants interleaved (repeat is the outer loop),
16 cases per scene sampled evenly across the trajectory. Scenes are real frames:
cloth-ball, armadillo-rollers, cloth-funnel.

Numbers below are medians in milliseconds, summed over the 16 sampled cases,
with the observed run-to-run spread. **A gap inside the spread is not a
result** — several conclusions in this project were retracted for exactly that,
and they are marked here rather than reported as ratios.

`fn=0` on every scene row: no configuration missed a collision.

---

## Narrow-phase mode

| scene | mode 0 scalar | mode 1 fast-vector | mode 2 TI-exact | spread |
|---|---:|---:|---:|---:|
| armadillo-rollers | 19.4 | 90.9 | **17.5** | 2–4% |
| cloth-ball | 190.7 | 223.4 | **117.7** | 1–2% |
| cloth-funnel | **6.4** | 15.9 | 8.1 | 6–10% |

**Mode 1 is demoted.** It loses on all three scenes and by 5.2× on
armadillo-rollers, which settles the "expected to lose, include it so the
demotion cites a number" item: it is a duplicate that loses, and it is
vertex-face only.

**Mode 0 and mode 2 both stay.** Neither dominates. Mode 2 wins cloth-ball by
1.62× and armadillo-rollers narrowly; mode 0 wins cloth-funnel by 1.27×, a 21%
gap against 6–10% spread, so that is a real win and not noise. Two kernels that
each win a third of the real scenes are not duplicates in the sense the keep bar
means, and mode 0 stays as the fallback.

**One prior claim is not reproduced**, and is now withdrawn — see the section
below. The plan carried "mode 2 is about 100× slower on armadillo edge-edge" as
the reason mode 0 had to survive. On this 16-case sample mode 2 is *faster* on
armadillo-rollers (17.5 vs 19.4), which at the time settled nothing: 16 cases
sampled across the trajectory can miss a pathological one, and the earlier
measurement was of the edge-edge phase alone rather than of a mixed sample.

## Withdrawn: "mode 2 is about 100× slower on armadillo edge-edge"

The way to settle this was never a trajectory average, so this is the targeted
run the section above asked for: **every** armadillo-rollers case, not a
subsample — 396 edge-edge and 385 vertex-face — for mode 0 and mode 2, on Grace
at `OMP_NUM_THREADS=72`, at the shipped defaults (`max_depth=69`, `tol=3e-8`).
Two repeats. Per-case rows, because the question is whether *one* case blows up
and a sum cannot see that.

**There is no 100× case. There is no 10× case.**

| | repeat 1 | repeat 2 |
|---|---:|---:|
| total edge-edge narrow, mode 0 | 526.8 ms | 485.0 ms |
| total edge-edge narrow, mode 2 | 473.2 ms | 438.6 ms |
| ratio | 0.90× | 0.90× |
| worst single case | 9.46× | 9.37× |
| cases above 10× | **0 of 396** | **0 of 396** |
| cases above 100× | **0 of 396** | **0 of 396** |

Summed over the whole trajectory mode 2 is 10% *faster* at edge-edge, the same
figure in both repeats. The distribution: median ratio 0.79, p95 2.71, p99 5.78.
Mode 2 is slower at all on 144 of 396 cases, slower by more than 2× on 46, by
more than 5× on 7, and by more than 10× on none.

**The tail is real, it is just an order of magnitude smaller than the claim.**
The worst case reproduces — `389ee` at 9.46× and 9.37× across the two repeats,
0.52 ms against 4.94 ms — as do `153ee` (5.27×, 5.60×) and `181ee` (5.09×,
4.95×). Others in the worst-ten do not: `3ee` reads 8.01× in one repeat and
0.16× in the other, which is what a sub-millisecond case measured once looks
like, and is the reason for running it twice. So there genuinely are edge-edge
configurations where the TightInclusion-exact kernel costs several times the
scalar reference; none of them costs a hundred times it, and all of them are
half-millisecond cases whose absolute cost is negligible against the 28 ms
`326ee` — where mode 2 is *faster* (17.1 ms against 28.1 ms).

`fn=0` for both modes on all 781 cases.

**Effect on the keep decision: none.** Mode 0 stays, on the evidence that already
supported it — cloth-funnel by 1.27×, outside the spread. What changes is that
the reason cited for it in the plan was wrong, and anyone reaching for "mode 2
explodes on edge-edge" as a design constraint should stop.

## Split rule

| scene | uniform | adaptive | gap | spread |
|---|---:|---:|---:|---:|
| armadillo-rollers | 19.3 | **18.7** | 3% | 4% |
| cloth-ball | 117.7 | **115.1** | 2% | 1–2% |
| cloth-funnel | 8.9 | **7.2** | 19% | 6–10% |

**Uniform is demoted.** Adaptive wins all three. Two of the three margins are
inside the spread and prove nothing on their own, but cloth-funnel's 19% is
outside it, and uniform never wins anywhere. A duplicate that never wins, reached
only through an undocumented environment variable, does not meet the keep bar.

## Broad phase

| scene | sweep | cell2d | spread |
|---|---:|---:|---:|
| armadillo-rollers | **22.1** | 35.1 | 4–8% |
| cloth-ball | **218.7** | 296.9 | 1–6% |
| cloth-funnel | 30.8 | **22.8** | 8–24% |

Pair counts agree exactly between the two on every scene, so these are
comparisons of the same work.

**Both stay** — that part of the plan holds. **But the default is wrong.**
`choose_broadphase_strategy` resolves `Auto` to Cell2D unconditionally, and on
two of the three real scenes the sweep is faster: by 1.59× on armadillo-rollers
and 1.36× on cloth-ball, both well outside the spread. Cell2D's win on
cloth-funnel is 26% against a 24% spread, which is marginal at best.

This does not contradict the earlier cell-list measurements, and the difference
is worth stating precisely because it decides what to do. Those were raw
box-list micro-benchmarks — 400k to 1.5M boxes, self-overlap, cell2d ahead by
4–7× — where the whole cost is the overlap search. These are end-to-end scene
steps, where the broad phase also pays for binning or sorting and the pair
density is whatever the frame happens to be. **The micro-benchmark and the
scene disagree, and the scene is what ships.**

Recommended: `Auto` should not be a constant. Either default to the sweep, which
wins the majority of real scenes here, or make `Auto` measure. What it must not
do is keep asserting Cell2D on the strength of a micro-benchmark that the
end-to-end numbers do not support. This is the one place where the plan's
"settled" label was premature; the keep decision survives, the default does not.

## Host versus device

| scene | phase | Grace host | Hopper CUDA | |
|---|---|---:|---:|---|
| armadillo-rollers | broad | 34.5 | **15.5** | 2.2× device |
| cloth-ball | broad | 293.0 | **61.2** | 4.8× device |
| cloth-funnel | broad | 23.1 | **17.4** | 1.3× device |
| armadillo-rollers | narrow | **17.6** | 371.0 | 21× host |
| cloth-ball | narrow | **115.8** | 145.3 | 1.25× host |
| cloth-funnel | narrow | **7.3** | 637.3 | 87× host |

Spreads 0–9%. `sccd_bench` performs its own warmup call, so the scene rows above
exclude CUDA context creation.

**The refinement rows do not**, and that matters when reading them:
`refine_scaling` has no warmup, each driver invocation is a fresh process, and
the job-level warmup only warms the page cache. The level-0 device broad phase
measures 2,039 ms against 0.92 ms at level 1 in the same process — that first
number is context creation, not broad-phase cost, and must not be read as one.

**The device broad phase is a clear win** — 1.3× to 4.8×, and it is the phase
whose shape suits a GPU (count, prefix sum, scatter, no sequential window walk).

**The device narrow phase is not** — but the two rows above do not say what
they were read as saying, and the section below retracts that reading.

## Retracted: "the device narrow phase loses on every scene"

The narrow rows above were run with `SCCD_NARROWPHASE_MODE=2` on both sides,
which looks like the fair comparison and is not one. **The mode enum names a
different kernel on each side.** On the host, mode 2 is the vectorised
TightInclusion-exact kernel in `src/sccd_vnarrowphase_ti.hpp` — the *fast* one,
which wins cloth-ball. On the device, mode 2 selects
`narrow_phase_generic<..., conservative=true>`, a scalar per-thread branch and
bound that is the device's *slowest* path. The sweep therefore raced the host's
best kernel against the device's worst and reported the difference as a hardware
verdict. Mode 0 — the default, and what a user gets without setting anything —
was never measured on the device at all.

Crossing the two factors, three repeats in one allocation on a GH200 node,
`OMP_NUM_THREADS=72`, 16 cases per scene, medians in milliseconds:

| scene | phase | host m0 | device m0 | | host m2 | device m2 | |
|---|---|---:|---:|---|---:|---:|---|
| cloth-funnel | narrow | **6.2** | 22.6 | 3.6× host | **7.2** | 592.2 | 83× host |
| armadillo-rollers | narrow | **17.9** | 40.1 | 2.2× host | **17.4** | 366.5 | 21× host |
| cloth-ball | narrow | 185.5 | **91.1** | 2.0× device | **115.4** | 121.1 | 1.05× host |

Spread 0–9%, and the mode-2 device column reproduces the original rows
(366.5 against 371.0, 121.1 against 145.3, 592.2 against 637.3), which is what
identifies the confound rather than a code change since the sweep.

**What actually holds.** With the default kernel the device narrow phase is
competitive: 2.0× *ahead* on the largest scene and 2.2–3.6× behind on the two
small ones. End to end, broad plus narrow at mode 0:

| scene | Grace host | Hopper CUDA | |
|---|---:|---:|---|
| cloth-funnel | **29.0** | 39.0 | 1.34× host |
| armadillo-rollers | **48.1** | 56.0 | 1.16× host |
| cloth-ball | 478.6 | **150.8** | 3.2× device |

So "a user is better served running the narrow phase on the host even when the
broad phase ran on the GPU" is not supported. On the largest scene the whole
pipeline on the GPU is 3.2× faster than the whole pipeline on 72 Grace cores,
and on the other two the device is within 35%. The recommendation to build and
measure a device-broad/host-narrow configuration stands on its own merits, but
not as a consequence of these numbers.

## The gap that is real: the conservative device kernel

Note what survives the retraction. Mode 0 is fast on both machines and mode 2 is
fast only on the host, so the cost of being conservative is:

| | mode 0 | mode 2 | conservative costs |
|---|---:|---:|---|
| host, cloth-funnel | 6.2 | 7.2 | **1.15×** |
| device, cloth-funnel | 22.6 | 592.2 | **26×** |

Same algorithm, same inputs, same tolerances — `compute_face_vertex_tolerance`
and the `(vf ? 30 : 28) * eps * min(max_coord, 1)^3` bound are transcribed on
both sides. A conservative search that costs 15% on one machine and 26× on the
other is a defect in one of the two implementations, not a property of
conservativeness. This is the device narrow phase's actual problem, and it is
device-internal rather than host-versus-device.

Two structural causes, both visible in the source:

1. **The device recomputes every corner it already has.** The host carries the
   eight corner values with each box (`TiBox::corner[3][8]`), and on a split
   evaluates only the four mid-face corners, inheriting the other four from the
   parent: 12 corner evaluations per split. The device's `Domain` is six bounds
   and nothing else, so `evaluate_cell_3d_policy` re-evaluates all eight corners
   of both children on all three axes: 48 per split, 4× the arithmetic for the
   same tree.
2. **Mode 0 accepts far earlier, so it never builds the tree.** Its acceptance
   conditions compare a *codomain* width against a *domain* tolerance, which
   fires long before TightInclusion's certified conditions do. Accepting a box is
   always safe however loose the test — it reports the box's `t` lower bound —
   so this is an accuracy difference, not a soundness one, and it is the price
   the conservative kernel pays for a tighter answer. The host pays it too.

Cause 1 is worth roughly 4× and is a contained change. It does not close 26× on
its own, so the tree-size difference between the two device kernels should be
counted before more is claimed.

### Counted. It is not the arithmetic, and cause 1 is not implementable

The counting turned both of those claims over, so they are corrected here rather
than left standing.

`SCCD_NP_COUNT_BOXES` compiles a counter into both the device kernel and the host
TightInclusion kernel that ticks once per box classified — the same unit on both
sides, which is what makes them comparable. It is off by default and a global
atomic when on, so an instrumented build's *timings* are meaningless; only the
counts are. On cloth-funnel, 16 cases, 843,140 queries:

| | boxes classified | per query |
|---|---:|---:|
| host conservative (mode 2) | 1,488,233 | **1.8** |
| device conservative (mode 2) | 796,343,700 | **944** |
| device mode 0 | 3,035,988 | 3.6 |

**The device conservative search does 535× the host's work on identical queries.**
Both use the same three acceptance conditions — origin containment padded by the
certified bound, whole range inside the error box, domain widths within the
domain tolerances — and the same `compute_face_vertex_tolerance`. The host decides
under two boxes per query, which is to say it classifies the root and stops. The
device classifies 944.

Two consequences:

- **The 26× time gap is work, not per-box cost.** Mode 2 does 259× more corner
  evaluations than mode 0 and takes 26× longer, so per evaluation it is about ten
  times *cheaper* — mode 0's time is dominated by fixed per-query setup, not by
  search. The conservative kernel is genuinely evaluation-bound.
- **Corner reuse cannot be built.** Measured with `-Xptxas=-v` on CUDA 12.6 for
  sm_90, the zero-stride kernel uses **238 of the 255 available registers** at
  `conservative=true` (224 at `conservative=false`), with no spills. Carrying
  eight corner values for three components is 24 doubles, 48 more registers, into
  a budget with 17 left. The shared stack cannot take them either: at 57,376 bytes
  per block it is already the second binding constraint. So the item is not "worth
  4×, contained" — it is not available at all in this kernel's shape.

Both modes sit at about 12.5% occupancy (238 registers × 128 threads leaves room
for two blocks per SM out of the 65,536 registers an SM has), so occupancy is not
what separates them either.

### Three things that looked like the cause and are not

Recorded so they are not retried.

**Sequential batching, to recover the host's collapsing time of impact.** The host
processes queries with a running global minimum, so later queries start against a
tiny `t` interval; the device launches everything at once against `max_toi`.
`SCCD_BATCH_SIZE` makes the device process candidates in sequential batches that
inherit the previous batch's bound, which should reproduce it. It does not
reproduce it, and it is not usable as a default:

| scene | mode 2, one batch | batch=1024 | |
|---|---:|---:|---|
| cloth-funnel | 611.6 ms | **48.0 ms** | 12.7× better |
| armadillo-rollers | 328.9 ms | 8699.2 ms | 26× worse |
| cloth-ball | 115.1 ms | 1899.8 ms | 16.5× worse |

The batch loop pays two blocking device-to-host copies per batch, and cloth-ball
runs two million queries per case, so batch=1024 is thirty-two thousand launches.
Worse, the 12.7× on cloth-funnel is **not** pruning: counted, batching moves the
work from 720,588,918 boxes to 671,597,526, a 7% reduction against a 12.7× time
change. Whatever earns that speedup, it is not the mechanism the knob was reached
for, and it does not generalise.

*(An earlier commit in this branch recorded "SCCD_BATCH_SIZE is monotonically
worse, 4211 ms at batch=256 on cloth-ball" as a general refutation. That
measurement was taken at the default mode 0, which averages 3.6 boxes per query
and so has no tree to prune. It is right about mode 0 and wrong as a general
statement, and the table above supersedes it.)*

**Contention on the shared time of impact.** Every resident block hits the single
word `toi[0]` with an atomic on every loop iteration, which at a few hundred
blocks and 944 boxes per query looks like an obvious serialisation. Refreshing
less often instead is monotonically worse everywhere:

| refresh every | cloth-funnel | armadillo-rollers | cloth-ball |
|---|---:|---:|---:|
| 1 iteration (shipped) | **638.7 ms** | **361.5 ms** | **134.4 ms** |
| 4 | 669.6 | 357.4 | 144.2 |
| 16 | 938.4 | 434.0 | 161.8 |
| 64 | 967.2 | 687.4 | 166.9 |
| 256 | 3377.5 | 1539.3 | 181.7 |

The pruning a fresh bound buys is worth more than the atomic costs, on every
scene. The knob was removed rather than shipped; `fn=0` throughout, as expected,
since a stale bound only ever prunes less.

### Where the 516× sits: the gap scales with query difficulty

The mean hides the shape, and the shape is the finding. Counting per query rather
than in total, cloth-funnel, 881,694 queries, host and device on the same set:

| | queries | boxes | per query | worst single query |
|---|---:|---:|---:|---:|
| host conservative | 881,694 | 1,545,887 | 1.8 | 48,561 |
| device conservative | 881,694 | 797,257,768 | 904.2 | **19,737,992** |

It is not a handful of queries exploding, and it is not a uniform slowdown.
**The two agree on the easy queries and diverge further the harder the query
gets.**

| boxes needed | host queries | device queries | ratio |
|---|---:|---:|---:|
| 0–1 | 617,471 | 587,953 | 1.0× |
| ≥ 8 | 42,406 | 103,502 | 2× |
| ≥ 64 | 1,073 | 33,902 | 32× |
| ≥ 1,024 | 164 | 11,581 | 71× |
| ≥ 16,384 | 3 | 1,386 | 462× |
| ≥ 1,048,576 | **0** | **208** | — |

The host's hardest query in the scene costs 48,561 boxes. The device has 208
queries costing over a million, and its worst costs 19.7 million — 406× the
host's worst. Below eight boxes the two distributions are the same to within a
few percent, so the device is not starting wrong: it is failing to *converge* on
the queries that need real search, while the host converges on all of them.

That rules out several things at once. It is not the seeding, which would show up
on the easy queries too. It is not the acceptance test or the tolerances, which
are transcribed and would shift the whole distribution rather than stretch its
tail. And it is not the split rule in the obvious sense: `bisect_ti_axis` picks
the axis with the largest width/tolerance ratio, the same rule as the host's
`ti_split_axis`.

What is left is what happens to a box once the search is deep — the interaction
between the depth cap, the `t` bound each query searches under, and the order in
which boxes come off the stack. The host pops in an order that lets `toi_q[q]`
tighten as it goes; the device's DFS pops from a shared stack that mixes queries,
and its bound is a single global minimum. Whether that is enough to explain a
distribution with a millionfold tail is the next thing to establish, and it wants
a single hard query traced box by box, not another aggregate.

*(A note on the instruments, since one of them nearly produced a wrong number.
The per-query counter initially covered only the zero-stride body, while the
global counter also saw the block-per-query kernel the benchmark's separate
per-query path uses. The two disagreed by 7×, which is what caught it. With every
evaluation site instrumented they agree — 820,651,600 against a bucket range of
600M–1,200M on the same run — and the 516× here is consistent with the 535×
measured before, so that figure stands.)*

## Fixed on the way past: an unsound rejection in the device's mode-0 kernel

Looking at mode 0 closely enough to explain the gap turned up a real defect in
it, and this is the third time this project has found the same one.

Rejecting a box is the only operation in the search that can lose a root, and it
is sound only when the origin-containment test is padded by at least the
certified numerical error bound, `(vf ? 30 : 28) * eps * min(max_coord, 1)^3`.
The device's `evaluate_cell_3d` padded it with the caller's distance tolerance
instead. The host's mode 0 does not — `srootfinder.hpp:525` pads by
`numerical_error[d]` — and neither does the device's quad kernel, whose comment
names this exact substitution as "the unsound rejection this project has already
found once". The triangle device kernel was the one still carrying it.

**Severity: latent, and no case was found that trips it.** The kernel computes in
double, so the bound is at most `30 * 2.22e-16 = 6.7e-15`, and every measurement
here passes `tol = 3e-8`, four hundred thousand times wider. The pad was
therefore too wide rather than too narrow on every scene in this document.

It was not reachable at `tol = 1e-16` either. `narrowphase_cuda_test` was written
partly to gate this, running 1,400 constructed queries -- full-mantissa
coordinates, a quarter of them grazing contacts where the root is
ill-conditioned -- at a tolerance below the bound, and the pre-fix kernel passes
it. The reason looks structural rather than a matter of not trying hard enough:
losing a root needs the *computed* fmin to land above the pad while the true one
is at or below zero, which needs the box to have shrunk until F is within
rounding error of zero across it, and by then the acceptance conditions have
long since fired. So the honest statement is that the pad was narrower than the
bound the guarantee is stated in terms of, not that a collision was being
dropped. It did fire historically in single precision, where the bound is
`3.6e-6` and a 3e-8 tolerance is *narrower* than the error in the corner values.

**Fix anyway.** The pad is now `max(tol, aerr[d])`, with `aerr` computed for both
device paths instead of only the conservative one. The acceptance conditions keep
using `tol`, since accepting is safe at any looseness. Soundness here should rest
on the certified bound rather than on the caller happening to pass a tolerance
wider than it, and the change costs nothing: at any realistic tolerance the max
selects `tol` and the kernel behaves exactly as before.

**Verified.** Rebuilt on GH200 and re-run: false positives and false negatives
are identical on all three scenes (274/0, 5642/0, 95424/0) and the timings are
unchanged within the 0–9% spread, which is what a pad that was already wider than
the new floor should do.

**One consequence for how the mode-0 rows are read.** Even with the rejection
fixed, mode 0's *acceptance* is looser than TightInclusion's, so the device's
mode-0 win on cloth-ball is a win by the less accurate kernel. Every
configuration here reports `fn=0` against the datasets' ground truth, but the
invariant this project holds is signed — a time of impact later than the true one
is the same failure as a miss — and `fn` cannot see that. The
conservative-to-conservative comparison is the mode-2 pair, which the device
loses everywhere.

## Quads

**Quads had no device narrow phase when this ran.** The three `hopper/refine-quad`
rows are `FAILED rc=134` — `sccd_smesh_CCD.hpp` raised `SMESH_ERROR` and aborted.
Recorded as explicit rows so the gap is in the data rather than an absence. It has
since been implemented (`src/cuda/sccd_narrowphase_vq.cu`) and these rows are kept
as the before picture; they have not been re-measured.

Host scaling, `SCCD_TOPOLOGY=quad` (hex cube skinned to QUADSHELL4) against the
tetrahedral cube skinned to TRISHELL3:

| level | tri faces | tri pairs | tri narrow ms | quad faces | quad pairs | quad narrow ms |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 48 | 3,324 | 1.26 | 24 | 4,512 | 1.78 |
| 1 | 192 | 58,152 | 0.79 | 96 | 80,112 | 0.36 |
| 2 | 768 | 951,384 | 3.20 | 384 | 1,315,632 | 4.28 |

These are the baseline for quad root-finder work, not a verdict. The geometry is
synthetic — the driver reflects the mesh through its centre, which makes every
swept box overlap every other — and the two columns are not a fair race: a hex
cube skins to half as many quads as a tet cube does triangles, yet produces
*more* pairs (1.32M against 951k at level 2), because each quad is larger and
its swept box overlaps more of its neighbours.

Per pair the two are level: 3.36 ns for triangles against 3.25 ns for quads at
level 2. That is better than "quads are untested" and weaker than an optimality
claim — it says the quad path is not pathologically slow on this geometry, not
that it is well optimised. Establishing the latter needs the oracle coverage and
the dedicated root finder, both of which are scheduled work.

What these rows do establish is that the quad path runs, scales, and can now be
measured, which it could not be before.

The quad gaps stay scheduled work, not demotion candidates: no device narrow
phase, `SCCD_NARROWPHASE_MODE` silently ignored (one root-finder variant, so the
conservative kernel is unreachable for quads), no oracle coverage, and no
optimised root finder of their own.

## Resolved: the sweep dropped touching pairs

The pair-count disagreement this run surfaced was **not** host versus device. It
was cell2d versus sweep, and the sweep was wrong.

The device runs the sweep (the device cell list is not wired into the CCD path),
while the host defaults to Cell2D, so the comparison labelled the two
implementations as two machines. Forcing `SCCD_BROADPHASE` on the host
reproduces it in one command, and isolates it to the edge-edge step:

| level | cell2d ee | sweep ee | vf (both) |
|---:|---:|---:|---:|
| 0 | 2,220 | 2,200 | 1,104 |
| 1 | 39,912 | 39,808 | 18,240 |

Exactly the 20 and 104 seen between Grace and Hopper.

**The defect.** The AABB predicate treats touching as overlapping — it rejects
only on a strict `amin > bmax` — but the sweep's candidate window advanced past
any box with `xmax <= fimin`:

```c
for (; begin < second_count; ++begin) {
    if (fimin < second_xmax[begin]) break;   // strict: skips xmax == fimin
}
```

So the window discarded pairs the predicate it feeds would have accepted. Boxes
are sorted by `xmin`, so for `j > fi` we always have
`xmax[j] >= xmin[j] >= fimin`, and the strict test can only fail on equality —
`xmin[j] == xmax[j] == fimin`, a **zero-extent box sitting exactly at another
box's lower bound on the sort axis**. That is what an axis-aligned face produces
when its swept AABB is flat, which is why a refined cube triggers it and the
three real scenes did not.

**Severity.** Missed pairs, not extra ones: a collision the narrow phase never
gets the chance to see. That is a false negative, which the conservativeness
invariant does not permit at any cost. It is also worse than the 0.6% headline
suggests — when the *sort axis* is the degenerate one, the loss is total. A
regression case of 2,000 flat coincident boxes finds 64,922 pairs with the cell
list and **zero** with the unfixed sweep. On the refined cube it showed up as
0.6% only because `choose_axis` picked a non-degenerate axis; nothing guarantees
that.

**Fix.** The window comparison is now inclusive (`fimin <= second_xmax[begin]`),
in `broadphase.hpp` and in both CUDA sweeps. cell2d and sweep now agree exactly
at every refinement level.

**Regression test.** `cell2d_broadphase_test` gained flat-coincident-box cases
that run the sweep on each of the three axes explicitly. Forcing the axis is
what makes the test work: `choose_axis` avoids the degenerate axis, so the
original random cases passed against the buggy sweep and would have kept
passing. Verified to fail before the fix and pass after.

**Verified on the hardware that showed it.** Rebuilt for Grace and Hopper and
re-run on a GH200 node, the two now agree at every refinement level:

| level | Grace host (vf / ee) | Hopper device (vf / ee) |
|---:|---|---|
| 0 | 1,104 / 2,220 | 1,104 / 2,220 |
| 1 | 18,240 / 39,912 | 18,240 / 39,912 |
| 2 | 294,144 / 657,240 | 294,144 / 657,240 |

`cell2d_broadphase_cuda_test` also passes, so the device cell list and the device
sweep agree as well.

**Effect on the numbers above.** None. The three real scenes had identical pair
counts for both broad phases before the fix, so the bug never fired there and
the timings stand.

## Summary of decisions

| item | decision | evidence |
|---|---|---|
| narrow-phase mode 1 | demote | loses all three scenes, 5.2× on armadillo |
| narrow-phase mode 0 | **keep** | wins cloth-funnel by 1.27×, outside spread |
| "mode 2 is 100× slower on armadillo edge-edge" | **withdrawn** | all 396 ee cases, twice: worst is 9.4×, none above 10×, mode 2 10% faster overall |
| narrow-phase mode 2 | keep | wins cloth-ball by 1.62× |
| uniform split | demote | never wins; adaptive takes all three |
| sweep broad phase | keep | wins armadillo 1.59×, cloth-ball 1.36× |
| cell2d broad phase | keep | wins cloth-funnel; wins micro-benchmarks 4–7× |
| `Auto` → Cell2D default | **revisit** | sweep wins 2 of 3 real scenes |
| device broad phase | keep, and wire into the CCD path | 1.3–4.8× over Grace |
| device narrow phase | keep; the "loses by up to 87×" reading is retracted | mode for mode: 2.0× ahead on cloth-ball, 2.2–3.6× behind on the others |
| device conservative kernel | **optimise** | being conservative costs the host 1.15× and the device 26× |
| quads | ship; gaps are scheduled work | device path aborted at the time of this run |
