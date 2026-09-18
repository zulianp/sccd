# Comparison: SCCD, Scalable CCD and Additive CCD

`compare-gh200-full-2026-09-18.csv.gz`, produced by `../sweep_comparison.sh` and
tabulated by `../compare_table.py` (`gunzip -k` it first; the raw CSV is 4.4 MB).
One GH200 module per job -- one Hopper, and one Grace bound to 72 CPUs -- over
whole scenes: every prepared case of armadillo-rollers (781), cloth-ball (79)
and cloth-funnel (577), three repeats each, as 15 resumable Slurm chunks.

**Double precision throughout.** smesh is built with `SMESH_GEOM_TYPE=float64`
and the frames are converted to match, so the mesh path and the dataset's exact
roots describe the same geometry. That is what makes the per-step comparison a
conservativeness test rather than a comparison of two geometries: under a float32
mesh the coordinates disagree at about `1e-7`, which near a grazing contact moves
the true time of impact by up to `1e-4`, and a TightInclusion-identical kernel
then reports a time of impact after the dataset's root on roughly a third of
armadillo's steps. In double that count is zero.

Every narrow phase runs on the same hardware: SCCD's host path and ACCD on 72 CPU
threads (OpenMP and oneTBB respectively), SCCD's device path and Scalable CCD on
the Hopper. All of them score identical broad-phase candidate lists, case by case.

Each competitor is compared on the question it answers, so there are two tables
per scene.

## Earliest time of impact: SCCD against Scalable CCD

Both libraries run their earliest-time-of-impact path once per case -- one query
kind of one simulation step, from a bound of 1. `narrow ms` is SCCD's
`ToiOutput::Earliest` call and Scalable CCD's default build, the library as its
authors ship it.

`err` is the earliest time of impact reported for a step minus the step's
earliest exact root, scored for every repeat; `late` counts steps where it is
positive, per pass.

## Per collision pair: SCCD against Additive CCD

`narrow ms` is the per-pair path over the broad-phase candidates: SCCD's
`ToiOutput::PerPair` call and ACCD over SCCD's host sweep output. Additive CCD
answers one pair at a time and has no parallelism of its own, so the harness
loops over the candidates with the `tbb::parallel_for` the toolkit's own
`Candidates::compute_collision_free_stepsize` uses; a serial loop would time one
thread against SCCD's 72.

`fp`, `fn` and the error columns score every curated query at the coordinates its
exact root was computed for, which is where `sccd_bench` scores SCCD.

Timings are per case in milliseconds, median / maximum. `ns/q` is narrow-phase
time over the candidates handed to it.

```
### armadillo-rollers -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |    4.607 |   1.485 /    7.70 |   0.600 /   27.28 |    6.661 |    7.948 |   -1.43e-05 |      -1.4e-07 |    0 |
| SCCD host Tight       |    3.758 |   1.461 /    7.87 |   0.720 /   17.91 |    6.086 |    8.368 |   -1.24e-06 |        -2e-08 |    0 |
| SCCD device Relaxed   |    1.630 |   1.942 /   14.85 |   1.672 /    3.55 |    5.177 |    15.75 |    -3.6e-05 |     -1.04e-06 |    0 |
| SCCD device Tight     |    1.645 |   1.954 /   14.50 |   1.950 /    9.75 |    5.503 |    18.52 |   -1.32e-06 |     -1.92e-08 |    0 |
| Scalable CCD device   |    1.293 |  12.533 /   27.45 |   4.223 /  108.61 |   18.559 |    57.27 |     0.00892 |         0.999 |  265 |
| Scalable CCD host     |    0.646 |   5.290 /   17.40 |   0.000 /    0.00 |    5.940 |        - |           - |             - |    - |

  steps with a known contact: 394

### armadillo-rollers -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    4.607 |   1.485 /    7.70 |   1.788 /   50.65 |    8.101 |    26.78 |   232 |     0 |    0 |  -4.66e-05 |          0 |      -0.941 |
| SCCD host Tight       |    3.758 |   1.461 /    7.87 |   2.084 /   49.43 |    7.638 |    32.56 |    24 |     0 |    0 |  -3.39e-06 |          0 |      -0.016 |
| SCCD device Relaxed   |    1.630 |   1.942 /   14.85 |   2.673 /   29.49 |    6.350 |    27.52 |   322 |     0 |    0 |  -9.22e-05 |          0 |      -0.942 |
| SCCD device Tight     |    1.645 |   1.954 /   14.50 |   2.979 /   22.95 |    6.670 |     31.4 |    24 |     0 |    0 |  -3.41e-06 |          0 |      -0.016 |
| ACCD host             |    4.648 |   1.432 /    9.89 |   0.600 /    7.49 |    6.534 |    6.612 |   444 |     0 |    0 |    -0.0286 |          0 |      -0.996 |

  curated queries with a contact per pass: 130859

### cloth-ball -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |    7.917 |  13.009 /   30.70 |   8.505 /  159.34 |   29.659 |    6.661 |    -3.8e-07 |             0 |    0 |
| SCCD host Tight       |    8.016 |  12.973 /   30.51 |   5.476 /   23.00 |   26.573 |    3.656 |   -1.66e-07 |             0 |    0 |
| SCCD device Relaxed   |    4.556 |   6.129 /   23.30 |   2.220 /    4.13 |   13.209 |   0.9954 |   -1.45e-06 |             0 |    0 |
| SCCD device Tight     |    4.574 |   6.282 /   22.61 |   2.320 /    4.50 |   13.470 |     1.04 |   -1.69e-07 |             0 |    0 |
| Scalable CCD device   |    1.900 |  19.343 /   28.32 |  28.040 /  220.87 |   50.105 |    21.83 |   -9.97e-06 |     -1.66e-07 |    0 |
| Scalable CCD host     |    1.194 |  21.727 /   28.54 |   0.000 /    0.00 |   22.893 |        - |           - |             - |    - |

  steps with a known contact: 43

### cloth-ball -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    7.917 |  13.009 /   30.70 |  15.860 /  118.12 |   36.004 |    10.81 |     2 |     0 |    0 |  -3.06e-07 |          0 |    -0.00054 |
| SCCD host Tight       |    8.016 |  12.973 /   30.51 |  17.903 /   87.92 |   40.790 |    9.169 |     1 |     0 |    0 |  -2.72e-07 |          0 |   -0.000193 |
| SCCD device Relaxed   |    4.556 |   6.129 /   23.30 |   7.866 /   72.41 |   18.047 |    5.205 |     2 |     0 |    0 |   -1.1e-06 |          0 |   -0.000898 |
| SCCD device Tight     |    4.574 |   6.282 /   22.61 |   8.827 /   38.97 |   18.549 |    4.906 |     1 |     0 |    0 |   -2.7e-07 |          0 |   -0.000193 |
| ACCD host             |    8.067 |  13.044 /   31.01 |   3.003 /    8.85 |   24.451 |     1.44 |    21 |     0 |    0 |    -0.0656 |          0 |       -0.94 |

  curated queries with a contact per pass: 664919

### cloth-funnel -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |    2.842 |   1.025 /    6.42 |   0.246 /    3.49 |    4.188 |    8.141 |           0 |             0 |    0 |
| SCCD host Tight       |    2.976 |   1.043 /    6.39 |   0.272 /   80.62 |    4.453 |    13.13 |           0 |             0 |    0 |
| SCCD device Relaxed   |    1.385 |   1.194 /    7.29 |   1.138 /    5.95 |    3.794 |    27.08 |           0 |             0 |    0 |
| SCCD device Tight     |    1.377 |   1.151 /    7.26 |   1.389 /    3.65 |    4.039 |    32.99 |           0 |             0 |    0 |
| Scalable CCD device   |    1.111 |  13.960 /   22.04 |   2.183 /   12.55 |   17.486 |    56.31 |           0 |             1 |    8 |
| Scalable CCD host     |    0.607 |   4.090 /   13.30 |   0.000 /    0.00 |    4.707 |        - |           - |             - |    - |

  steps with a known contact: 363

### cloth-funnel -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    2.842 |   1.025 /    6.42 |   0.482 /    4.38 |    4.487 |    14.28 |   687 |     0 |    0 |  -0.000357 |          0 |          -1 |
| SCCD host Tight       |    2.976 |   1.043 /    6.39 |   0.579 /  222.51 |    4.796 |    28.92 |    15 |     0 |    0 |  -6.87e-06 |          0 |          -1 |
| SCCD device Relaxed   |    1.385 |   1.194 /    7.29 |   1.176 /    4.00 |    3.875 |    28.76 |   742 |     0 |    0 |  -0.000522 |          0 |          -1 |
| SCCD device Tight     |    1.377 |   1.151 /    7.26 |   1.950 /   16.34 |    4.582 |    42.21 |    15 |     0 |    0 |  -6.83e-06 |          0 |          -1 |
| ACCD host             |    3.080 |   1.004 /    6.72 |   0.370 /   55.80 |    4.598 |    14.37 |   422 |     0 |    0 |   -0.00857 |          0 |      -0.956 |

  curated queries with a contact per pass: 6773

  Time-of-impact error is signed: reported minus exact root, so negative is
  conservative and a positive value is a contact reported after the true one.
  Counts are per pass over the case list. The mesh path and the exact roots
  describe the same geometry here, because smesh is built with double geom_t;
  a float32 mesh answers about different coordinates and makes correct kernels
  look late.
```

## What the numbers say

**Conservativeness.** SCCD reports no missed contact and no time of impact after
the true root anywhere: zero late steps on all three scenes in every mode, and
zero late answers over `130,859`, `664,919` and `6,773` curated contacts per
pass. Its worst signed step error is negative on every scene. ACCD is likewise
never late and never misses.

Scalable CCD's earliest time of impact lands after the true one on `265` of `394`
armadillo steps per pass, by up to `0.999` -- reporting no contact at all for
steps that have one -- and on `8` of `363` cloth-funnel steps. On cloth-ball it is
conservative throughout. The cause is a race in its subdivision buffer's overflow
check, described in `../README.md` and diagnosed in `wip/DECISIONS.md` section 8;
it bites when a candidate list is short relative to the search it starts, which is
armadillo's situation and not cloth-ball's, and its answers vary between runs.
Double precision does not change this: the same run under a float32 mesh gave
`271`.

**Against Scalable CCD, on cost.** At the median, end to end, SCCD on the device
is `3.4x` faster on armadillo (`5.50` against `18.56` ms per case), `3.7x` on
cloth-ball (`13.47` against `50.11`) and `4.3x` on cloth-funnel (`4.04` against
`17.49`). Per candidate its narrow phase costs `18.5` against `57.3` ns on
armadillo, `1.04` against `21.8` on cloth-ball and `33.0` against `56.3` on
cloth-funnel. On the host, where Scalable CCD has only a broad phase, SCCD's prep
plus query is `5.22` against `5.94` ms on armadillo, `21.0` against `22.9` on
cloth-ball and `4.02` against `4.70` on cloth-funnel.

**Against ACCD, per collision pair, the cost goes the other way.** Over the same
candidates on the same 72 threads, ACCD is the cheaper narrow phase: `6.6`
against SCCD Tight's `32.6` ns per candidate on armadillo, `1.44` against `9.17`
on cloth-ball, `14.4` against `28.9` on cloth-funnel -- two to six times.
Conservative advancement is cheap per pair precisely because it stops at a bound
instead of isolating a root.

**Tightness separates them by three to four orders of magnitude.** Per pair, SCCD
Tight's median error is `-3.4e-6` on armadillo, `-2.7e-7` on cloth-ball and
`-6.9e-6` on cloth-funnel; ACCD's is `-2.9e-2`, `-6.6e-2` and `-8.6e-3`. That is
`conservative_rescaling = 0.9` doing what it is for, and it shows in the
false-positive counts too: `444` against `24` on armadillo. Within SCCD, Relaxed
trades the same way against Tight -- `232` false positives against `24` on
armadillo, `687` against `15` on cloth-funnel.

**cloth-funnel starts in contact.** `3,570` of the scene's `7,552` curated
queries have an exact root of `0`: the primitives already touch when the step
begins. ACCD reports `0` for them and says so ("Initial distance 0 ≤ d_min=0"),
which is conservative. Counting its broad-phase candidates too, that warning
fires about `7,000` times per pass and fills the scene's `.accd.err` files; the
other two scenes produce none. It is also why `worst early` reaches `-1` there:
a query whose root is late in the step but whose primitives touch at the start.

## What this run does not establish

Three scenes -- the ones with verified ground truth -- on one machine, three
repeats. n-body-simulation, rod-twist and puffer-ball are not prepared here;
puffer-ball has no runnable case until its frames are extracted.

The two host narrow phases are threaded by different runtimes, SCCD's by OpenMP
and ACCD's by oneTBB, both over 72 threads of one Grace. That is the closest
available comparison, not an identical one.

With Scalable CCD's buffer race removed in a patched copy, `47` per-query times
of impact on armadillo were still late by up to `3.0e-3`, where SCCD and ACCD on
the same inputs have none. That was measured under the float32 mesh and has no
established cause; the patched copy is not what this table measures.
