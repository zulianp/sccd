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

**One prior claim is not reproduced.** The plan carried "mode 2 is about 100×
slower on armadillo edge-edge" as the reason mode 0 had to survive. On this
sample mode 2 is *faster* on armadillo-rollers (17.5 vs 19.4). That does not
refute the blowup: this is 16 cases sampled across 1600, so the pathological
case may simply not be in the sample, and the earlier measurement was of the
edge-edge phase alone rather than of a mixed sample. **The armadillo edge-edge
blowup needs a targeted run against the specific case, not a trajectory
average.** Until then, treat the 100× figure as unconfirmed rather than as
either established or withdrawn — and note that the conclusion it supported
(keep mode 0) is independently supported by cloth-funnel above.

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

**Severity: latent, not firing.** The kernel computes in double, so the bound is
at most `30 * 2.22e-16 = 6.7e-15`, and every measurement here passes
`tol = 3e-8`, four hundred thousand times wider. The pad was therefore too wide
rather than too narrow on every scene in this document, which is why nothing
showed it. It fires for any caller who asks for a tolerance below the bound, and
it fired historically when the kernel computed in single precision, where the
bound is `3.6e-6` and a 3e-8 tolerance is *narrower* than the error in the corner
values.

**Fix.** The pad is now `max(tol, aerr[d])`, with `aerr` computed for both
device paths instead of only the conservative one. The acceptance conditions
keep using `tol`, since accepting is safe at any looseness.

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
| narrow-phase mode 2 | keep | wins cloth-ball by 1.62× |
| uniform split | demote | never wins; adaptive takes all three |
| sweep broad phase | keep | wins armadillo 1.59×, cloth-ball 1.36× |
| cell2d broad phase | keep | wins cloth-funnel; wins micro-benchmarks 4–7× |
| `Auto` → Cell2D default | **revisit** | sweep wins 2 of 3 real scenes |
| device broad phase | keep, and wire into the CCD path | 1.3–4.8× over Grace |
| device narrow phase | keep; the "loses by up to 87×" reading is retracted | mode for mode: 2.0× ahead on cloth-ball, 2.2–3.6× behind on the others |
| device conservative kernel | **optimise** | being conservative costs the host 1.15× and the device 26× |
| quads | ship; gaps are scheduled work | device path aborted at the time of this run |
