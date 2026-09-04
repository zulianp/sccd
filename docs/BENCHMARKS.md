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
`Tight` (2) compares domain widths against domain tolerances. Both ship; modes 1
and 3 are validation-only and need `SCCD_ENABLE_TIGHT_INCLUSION=ON`.

**Strides.** `toi_stride=0` is `find_earliest_impact_time` — one time of impact
for the step, so every query prunes against the running minimum. `toi_stride=1`
is `find_impact_times` — one result per candidate, no shared bound.

---

## Narrow phase, milliseconds

| scene | mode | CPU s0 | GPU s0 | CPU s1 | GPU s1 |
|---|---|---:|---:|---:|---:|
| cloth-funnel | Fast | **6.2** | 23.1 | **9.7** | 56.3 |
| cloth-funnel | Tight | **7.1** | 36.5 | **12.3** | 789.9 |
| armadillo-rollers | Fast | **18.0** | 31.7 | **62.8** | 218.8 |
| armadillo-rollers | Tight | **17.4** | 56.3 | **77.3** | 3126.7 |
| cloth-ball | Fast | 184.1 | **102.1** | **308.4** | 612.2 |
| cloth-ball | Tight | 114.9 | **107.9** | **241.2** | 902.8 |

Stride 1 costs the CPU 1.6–4.4× over stride 0 — the value of the shared running
minimum. It costs the GPU 2.4–55×, because stride 1 runs a different kernel: one
block per query, seeded with a 128-way dice, which classifies 2.58 million boxes
per query against the host's 1,844.

## Broad phase, milliseconds

| scene | CPU | GPU | |
|---|---:|---:|---|
| cloth-funnel | 22.4 | **17.4** | 1.3× GPU |
| armadillo-rollers | 30.1 | **15.7** | 1.9× GPU |
| cloth-ball | 291.3 | **59.5** | 4.9× GPU |

Two implementations ship — sweep-and-prune and a 2D cell list — and they produce
identical pair sets. `sccd_broadphase_strategy.hpp` picks between them by racing
them at run time; five attempts to predict the winner from geometry all failed.

## End to end, each side at its best mode

| scene | CPU | GPU | |
|---|---:|---:|---|
| cloth-funnel | **29.2** | 38.4 | 1.3× CPU |
| armadillo-rollers | 47.8 | 47.3 | parity |
| cloth-ball | 408.8 | **159.4** | 2.6× GPU |

## Accuracy against exact roots

Signed error over queries with a real collision. **Zero late times of impact and
zero missed collisions in all twelve configurations**, across 101,164 roots. Late
would be a correctness failure; early only costs a solver step size.

| scene | roots | Fast CPU | Fast GPU | Tight CPU | Tight GPU |
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

The device does 94× the host's work on the earliest-impact path. The per-query
path is worse by an order of magnitude and is the largest open item; see
`wip/TODO.md`.
