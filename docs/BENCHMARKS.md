# Benchmarks

All figures below are measured, not estimated. Reproduce with `sccd_bench`; the
raw rows are committed under `benchmark/assessment/`.

**Platform.** One GH200 node: Grace (72 threads, `OMP_NUM_THREADS=72`) and Hopper
sm_90 in the same allocation. `tol = 3e-8`, `max_depth = 69`, 16 cases per scene,
three repeats, medians reported. Run-to-run spread is 0–9%; a difference inside
that band is not a result.

**Scenes.** Real frames, not synthetic motion.

| scene | broad-phase pairs | exact roots |
|---|---:|---:|
| cloth-funnel | 843,140 | 104 |
| armadillo-rollers | 3,205,151 | 5,636 |
| cloth-ball | 33,329,729 | 95,424 |

**Modes.** `Relaxed` (0) is the scalar search with the looser acceptance test.
`Tight` (2) compares domain widths against domain tolerances. On the CPU
`Relaxed` is the cheaper of the two; **on the GPU it is not** — it costs 1.5×
`Tight`'s boxes on armadillo-rollers at `toi_stride=1` and about the same on
cloth-ball. A looser acceptance ends a box sooner but reports a time of impact
further before the true one, and a looser bound prunes less, so the queries that
follow do more work. Both ship; modes 1
and 3 are validation-only and need `SCCD_ENABLE_TIGHT_INCLUSION=ON`.

**Strides.** `toi_stride=0` is `find_earliest_impact_time` — one time of impact
for the step, so every query prunes against the running minimum. `toi_stride=1`
is `find_impact_times` — one result per candidate, no shared bound.

---

## Narrow phase, milliseconds

| scene | mode | CPU s0 | GPU s0 | CPU s1 | GPU s1 |
|---|---|---:|---:|---:|---:|
| cloth-funnel | Relaxed | **6.2** | 24.3 | 9.7 | **22.2** |
| cloth-funnel | Tight | **7.2** | 27.2 | **12.2** | 33.4 |
| armadillo-rollers | Relaxed | **17.9** | 32.4 | 62.9 | **66.8** |
| armadillo-rollers | Tight | **17.5** | 37.5 | 77.1 | **65.2** |
| cloth-ball | Relaxed | 185.9 | **105.6** | 309.6 | **154.8** |
| cloth-ball | Tight | 115.2 | **103.0** | 241.6 | **158.4** |

Each cell is the sum over the scene's 16 cases within one repeat, then the median
over three repeats; vertex-face and edge-edge cases are summed together. Timings
come from a build **without** `-DSCCD_NP_COUNT_BOXES`, which puts a global atomic
on the hot path and makes an instrumented build 20× slower on the host.

Stride 1 costs the CPU 1.6–4.4× over stride 0 — the value of the shared running
minimum. It costs the GPU 2.6–6.6×, because stride 1 runs a different kernel: one
block per query, seeded with a 128-way dice.

The GPU stride-1 column improved by an order of magnitude against the previous
measurement of this table, from two changes:

| `Tight`, GPU s1 | was | now | |
|---|---:|---:|---:|
| armadillo-rollers | 3126.7 | **65.2** | **48×** |
| cloth-funnel | 789.9 | **33.4** | **24×** |
| cloth-ball | 902.8 | **158.4** | **5.7×** |

**The per-query path now runs one thread per query**, as the earliest-impact path
always has. It used to give each query a whole block, which is bound by scheduling
the blocks rather than by the search — 843,414 of them on cloth-funnel, each
classifying 63 boxes across 128 threads. On its own that change is 4.5×, 2.8× and
4.4×. `SCCD_NP_S1_BLOCK_PER_QUERY=1` restores the old kernel.

The device's edge-edge domain tolerances had been assigned to the wrong axes: `u`
got the `t` tolerance, up to 286× too large, and the `v` denominator was never
computed. The split rule picks the axis with the largest width/tolerance ratio, so
`u` was starved of splits and the search refined `t` and `v` against a `u`
interval pinned at full width. Nothing unsafe followed — a larger domain tolerance
only accepts sooner, which is the conservative direction, and every gate stayed
green — but on one query it cost 53 million boxes where the host needed 55.

The CPU columns reproduce the previous measurement to within 1–7% and are
unaffected.

## Broad phase, milliseconds

| scene | CPU | GPU | |
|---|---:|---:|---|
| cloth-funnel | 22.4 | **17.4** | 1.3× GPU |
| armadillo-rollers | 30.1 | **15.7** | 1.9× GPU |
| cloth-ball | 291.3 | **59.5** | 4.9× GPU |

Two implementations ship — sweep-and-prune and a 2D cell list — and they produce
identical pair sets. `sccd_broadphase_strategy.hpp` picks between them by racing
them at run time, because which one wins is not predictable from the geometry.

## End to end, each side at its best mode

| scene | CPU | GPU | |
|---|---:|---:|---|
| cloth-funnel | **29.2** | 38.4 | 1.3× CPU |
| armadillo-rollers | 47.8 | 47.3 | parity |
| cloth-ball | 408.8 | **159.4** | 2.6× GPU |

> These predate the narrow-phase fix above. The GPU column is a sum of prep,
> broad phase and narrow phase, and its narrow-phase term has fallen; refreshing
> it needs the broad-phase timings measured in the same run, which the run behind
> the table above did not capture.

## Accuracy against exact roots

Signed error over queries with a real collision. **Zero late times of impact and
zero missed collisions in all twelve configurations**, across 101,164 roots. Late
would be a correctness failure; early only costs a solver step size.

| scene | roots | Relaxed CPU | Relaxed GPU | Tight CPU | Tight GPU |
|---|---:|---:|---:|---:|---:|
| cloth-funnel | 104 | 1.18e-01 | 1.60e-01 | **1.72e-03** | 9.40e-02 |
| armadillo-rollers | 5,636 | 4.04e-04 | 5.56e-04 | **1.83e-05** | **1.83e-05** |
| cloth-ball | 95,424 | 4.08e-07 | 1.69e-06 | 3.85e-07 | **1.72e-07** |

Median earliness; times of impact are in `[0, 1]` over the step. Worst cases for
`Relaxed` are 9.41e-01, 4.16e-02 and 1.29e-04.

`Tight` is **69× tighter** than `Relaxed` at the median on cloth-funnel and **22×**
on armadillo-rollers. That is what the two modes trade: `Relaxed` accepts boxes
sooner, and on grazing contacts that costs most of the step.

Accuracy is measured on the datasets' curated query sets, which carry exact
roots and whose coordinates are exact dyadic rationals. It cannot be measured
through the mesh path: smesh stores coordinates as `float`, so the mesh geometry
is a rounded copy of the geometry those roots belong to.

## Hit and miss, against TightInclusion and against exact roots

`ti_oracle`, every query of every curated query set, edge-edge, one GH200
allocation. `FP`/`FN` are against TightInclusion's own answer; `gtMISS`/`gtLATE`
are against the datasets' exact roots and are the gate.

| scene | queries | TI hits | mode | hits | FP | FN | gtMISS | gtLATE |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| cloth-funnel | 6,751 | 6,259 | `tight` CPU | 6,259 | **0** | **0** | 0 | 0 |
| | | | `tight` GPU | 6,259 | **0** | **0** | 0 | 0 |
| | | | `relaxed` CPU | 6,700 | 441 | 0 | 0 | 0 |
| | | | `relaxed` GPU | 6,734 | 475 | 0 | 0 | 0 |
| armadillo-rollers | 99,104 | 98,761 | `tight` CPU | 98,761 | **0** | **0** | 0 | 0 |
| | | | `tight` GPU | 98,761 | **0** | **0** | 0 | 0 |
| | | | `relaxed` CPU | 98,895 | 134 | 0 | 0 | 0 |
| | | | `relaxed` GPU | 98,930 | 169 | 0 | 0 | 0 |
| cloth-ball | 557,683 | 557,668 | `tight` CPU | 557,668 | **0** | **0** | 0 | 0 |
| | | | `tight` GPU | 557,668 | **0** | **0** | 0 | 0 |
| | | | `relaxed` CPU | 557,669 | 1 | 0 | 0 | 0 |
| | | | `relaxed` GPU | 557,669 | 1 | 0 | 0 | 0 |

**`Tight` reproduces TightInclusion's hit set exactly on both machines**, over
663,538 queries: no false positive, no false negative, no missed collision and no
late time of impact. The gate exits zero on all three scenes.

`Relaxed` accepts sooner by design, so it reports more hits — 441 and 475 extra on
cloth-funnel out of 6,751. Those are false positives, which cost work and never
safety, and the two implementations of the looser predicate do not agree query for
query, which is expected: `Relaxed` is a different search on each machine, unlike
`Tight`.

One difference worth naming rather than hiding. The host's `Tight` is bit-identical
to TightInclusion (`abserr_max` 0.000e+00 on every scene); the device's reproduces
the hit set exactly but its times of impact sit up to 5.25e-3, 4.43e-4 and 1.76e-6
*after* TightInclusion's. TightInclusion's answer is itself a lower bound on the
truth, so being later than it is not a violation — `gtLATE` is zero — and the
`lateTI` column that counts it is documented as over-reporting for that reason.

## Device narrow phase, occupancy

Measured with `-Xptxas=-v` on CUDA 12.6 for sm_90, zero-stride kernel:

| | registers | shared memory | blocks/SM |
|---|---:|---:|---:|
| `Relaxed` | 224 / 255 | 3,584 B | 2 |
| `Tight` | 238 / 255 | 3,584 B | 2 |

No spills. Registers bind occupancy at roughly 12.5%, not shared memory, so
shrinking the block-local stack buys load balance rather than residency.

## Work, in boxes classified

Counted with `-DSCCD_NP_COUNT_BOXES` (off by default), cloth-funnel, `Tight`:

| | boxes | per query |
|---|---:|---:|
| host | 1,038,395 | 1.2 |
| device, `toi_stride=0` | 97,757,836 | 111 |
| device, `toi_stride=1` | 705,986,614 | 2,576,594 |

These counts are from the same pre-fix build as the timing table above. After the
edge-edge tolerance fix the same cloth-funnel run classifies 1.35 M boxes on the
earliest-impact path and 129 M on the per-query path, against 4.08 M and 1.78 G
before. The investigation is in
[`../wip/CUDA_NARROWPHASE_PLAN.md`](../wip/CUDA_NARROWPHASE_PLAN.md).
