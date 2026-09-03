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

**The device narrow phase is not.** It loses to 72 Grace cores on every scene,
by 21× and 87× on two of them. That is far too large to be a tuning gap and
points at the branch-and-bound structure rather than at constants: divergent
per-query stack depth, the global-stack overflow path, or occupancy. The plan
ships it as "the only conservative device path", which is still true, but these
numbers say a user is better served running the narrow phase on the host even
when the broad phase ran on the GPU. **The obvious configuration to support is
device broad phase feeding host narrow phase**, and it should be measured before
the CUDA narrow phase is presented as the device story.

## Quads

**Quads have no device narrow phase at all.** The three `hopper/refine-quad`
rows are `FAILED rc=134` — `sccd_smesh_CCD.hpp` raises `SMESH_ERROR` and aborts.
Recorded as explicit rows so the gap is in the data rather than an absence.

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

## Open problem: host and device disagree on the pair set

On the synthetic refinement geometry, and only there, the two find different
numbers of pairs:

| level | Grace | Hopper | missing |
|---:|---:|---:|---:|
| 0 | 3,324 | 3,304 | 20 |
| 1 | 58,152 | 58,048 | 104 |
| 2 | 951,384 | 950,920 | 464 |

About 0.05–0.6%, consistently with the **device finding fewer**. The three real
scenes agree exactly, and every scene row reports `fn=0`, so nothing here shows a
missed collision yet.

It still has to be explained before these rows are used, and it is the dangerous
direction: a pair the device does not report is a collision its narrow phase
never gets the chance to find, and the conservativeness invariant makes that a
defect rather than a tolerance. The leading candidate is precision — smesh's
`geom_t` is `float`, so device and host AABBs can straddle a boundary
differently — which would be benign only if the device rounds outward. Worth
noting that the geometry showing it is the degenerate all-pairs case, where
nearly every box touches nearly every other and borderline overlaps are
therefore everywhere.

Next step: compare the pair sets themselves rather than their counts, on the
smallest level, and check whether the missing pairs are boxes that touch within
a float ulp.

---

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
| device narrow phase | keep, but do not present as the device story | loses by up to 87× |
| quads | ship; gaps are scheduled work | device path aborts today |
