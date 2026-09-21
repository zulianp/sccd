# Comparison: SCCD, Scalable CCD and Additive CCD

`compare-gh200-full-2026-09-21.csv.gz`, produced by `../sweep_comparison.sh` and
tabulated by `../compare_table.py` (`gunzip -k` it first; the raw CSV is 26 MB).
One GH200 module per job -- one Hopper, and one Grace bound to 72 CPUs -- over
whole scenes: every prepared case of armadillo-rollers (781), cloth-ball (79),
cloth-funnel (577), n-body-simulation (146), puffer-ball (240) and rod-twist
(4,571), `6,394` in all, three repeats each, as 57 resumable Slurm chunks.

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

**Scalable CCD is compared on the device only.** Its narrow phase is CUDA, so its
host side is a broad phase with no pipeline behind it and not something a caller
runs end to end. SCCD's device runs twice, over the cell list and over the sweep,
so both of our broad-phase strategies and Scalable CCD's sit in one allocation
and can be put on one axis; between allocations this harness varies by about
40%. The host stays on the sweep, because ACCD is handed the candidates it
produces.

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
| SCCD host Relaxed     |    4.663 |   1.464 /    7.86 |   0.597 /   27.74 |    6.732 |    7.923 |   -1.44e-05 |      -1.4e-07 |    0 |
| SCCD host Tight       |    5.028 |   1.520 /    8.21 |   0.737 /   17.72 |    7.188 |    8.515 |   -1.24e-06 |        -2e-08 |    0 |
| SCCD device Relaxed   |    1.322 |   1.547 /   14.69 |   1.841 /   25.85 |    4.900 |    21.29 |   -3.59e-05 |     -1.04e-06 |    0 |
| SCCD device Tight     |    1.289 |   1.555 /   14.48 |   2.113 /   27.43 |    5.200 |    23.81 |   -1.32e-06 |     -1.92e-08 |    0 |
| Scalable CCD device   |    1.329 |  12.713 /   21.62 |   4.178 /  109.31 |   18.514 |    55.36 |     0.00976 |             1 |  272 |

  steps with a known contact: 394

### armadillo-rollers -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    4.663 |   1.464 /    7.86 |   1.812 /   51.07 |    8.360 |    26.76 |   232 |     0 |    0 |  -4.66e-05 |          0 |      -0.941 |
| SCCD host Tight       |    5.028 |   1.520 /    8.21 |   2.124 /   49.31 |    8.908 |    32.86 |    24 |     0 |    0 |  -3.39e-06 |          0 |      -0.016 |
| SCCD device Relaxed   |    1.322 |   1.547 /   14.69 |   2.761 /   30.85 |    5.588 |    28.32 |   322 |     0 |    0 |  -9.29e-05 |          0 |      -0.942 |
| SCCD device Tight     |    1.289 |   1.555 /   14.48 |   3.018 /   21.98 |    5.907 |    31.94 |    24 |     0 |    0 |   -3.4e-06 |          0 |      -0.016 |
| ACCD host             |    3.836 |   1.409 /    9.79 |   0.610 /    6.59 |    5.961 |    6.622 |   444 |     0 |    0 |    -0.0286 |          0 |      -0.996 |

  curated queries with a contact per pass: 130859

### cloth-ball -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |    7.307 |  13.005 /   31.55 |   8.413 /   37.90 |   28.843 |     5.81 |    -3.8e-07 |             0 |    0 |
| SCCD host Tight       |    7.479 |  13.101 /   31.84 |   5.455 /   22.95 |   26.122 |    3.629 |   -1.66e-07 |             0 |    0 |
| SCCD device Relaxed   |    2.830 |   6.669 /   22.68 |   2.702 /   14.07 |   14.210 |    2.013 |   -1.47e-06 |             0 |    0 |
| SCCD device Tight     |    2.788 |   6.721 /   22.65 |   2.871 /   14.21 |   14.232 |    2.058 |    -1.8e-07 |             0 |    0 |
| Scalable CCD device   |    2.066 |  19.351 /   28.54 |  27.877 /  224.90 |   50.246 |    22.13 |   -9.97e-06 |     -1.66e-07 |    0 |

  steps with a known contact: 43

### cloth-ball -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    7.307 |  13.005 /   31.55 |  16.266 /  108.45 |   34.955 |    10.29 |     2 |     0 |    0 |  -3.06e-07 |          0 |    -0.00054 |
| SCCD host Tight       |    7.479 |  13.101 /   31.84 |  12.278 /   87.62 |   31.561 |      7.9 |     1 |     0 |    0 |  -2.72e-07 |          0 |   -0.000193 |
| SCCD device Relaxed   |    2.830 |   6.669 /   22.68 |   7.695 /   73.15 |   16.064 |    5.089 |     2 |     0 |    0 |   -1.1e-06 |          0 |   -0.000923 |
| SCCD device Tight     |    2.788 |   6.721 /   22.65 |   8.762 /   39.73 |   16.821 |    4.896 |     1 |     0 |    0 |   -2.7e-07 |          0 |   -0.000193 |
| ACCD host             |    4.460 |  12.770 /   31.42 |   3.285 /    8.14 |   21.456 |    1.526 |    21 |     0 |    0 |    -0.0656 |          0 |       -0.94 |

  curated queries with a contact per pass: 664919

### cloth-funnel -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |    2.880 |   1.012 /    6.21 |   0.247 /    3.34 |    4.264 |        8 |           0 |             0 |    0 |
| SCCD host Tight       |    2.626 |   1.013 /    6.09 |   0.262 /   81.52 |    4.350 |    12.92 |           0 |             0 |    0 |
| SCCD device Relaxed   |    0.963 |   1.805 /    7.36 |   1.354 /    5.25 |    4.274 |    42.83 |           0 |             0 |    0 |
| SCCD device Tight     |    0.951 |   1.807 /    7.36 |   1.609 /    5.97 |    4.375 |    47.87 |           0 |             0 |    0 |
| Scalable CCD device   |    1.068 |  14.210 /   24.88 |   2.126 /   18.11 |   17.627 |     54.9 |           0 |             1 |   10 |

  steps with a known contact: 363

### cloth-funnel -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    2.880 |   1.012 /    6.21 |   0.471 /    4.37 |    4.523 |    14.14 |   687 |     0 |    0 |  -0.000357 |          0 |          -1 |
| SCCD host Tight       |    2.626 |   1.013 /    6.09 |   0.542 /  227.01 |    4.631 |    28.57 |    15 |     0 |    0 |  -6.87e-06 |          0 |          -1 |
| SCCD device Relaxed   |    0.963 |   1.805 /    7.36 |   1.140 /    3.56 |    3.894 |    28.18 |   742 |     0 |    0 |  -0.000539 |          0 |          -1 |
| SCCD device Tight     |    0.951 |   1.807 /    7.36 |   1.852 /   16.41 |    4.424 |    40.46 |    15 |     0 |    0 |   -6.9e-06 |          0 |          -1 |
| ACCD host             |    3.893 |   1.066 /    7.84 |   0.380 /   55.13 |    5.482 |    13.98 |   422 |     0 |    0 |   -0.00857 |          0 |      -0.956 |

  curated queries with a contact per pass: 6773

### n-body-simulation -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |    7.778 |  72.604 /  147.72 |  54.891 /   97.48 |  131.660 |    2.707 |           0 |             0 |    0 |
| SCCD host Tight       |    7.650 |  72.178 /  161.35 |  44.276 /   81.23 |  118.919 |    2.264 |           0 |             0 |    0 |
| SCCD device Relaxed   |    4.093 |  18.241 /  106.37 |   4.315 /  103.63 |   53.906 |     1.26 |           0 |             0 |    0 |
| SCCD device Tight     |    4.078 |  18.308 /  105.31 |   4.106 /  108.64 |   53.575 |    1.275 |           0 |             0 |    0 |
| Scalable CCD device   |    2.945 |  39.883 /   53.79 | 215.337 /  470.45 |  254.285 |    10.65 |           0 |             0 |    0 |

  steps with a known contact: 74

### n-body-simulation -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    7.778 |  72.604 /  147.72 | 149.907 /  305.57 |  227.337 |    6.925 |     9 |     0 |    0 |  -1.82e-08 |          0 |    -0.00127 |
| SCCD host Tight       |    7.650 |  72.178 /  161.35 |  89.396 /  264.41 |  179.368 |    7.018 |    11 |     0 |    0 |  -3.02e-08 |          0 |    -0.00024 |
| SCCD device Relaxed   |    4.093 |  18.241 /  106.37 |  34.625 /  128.94 |   50.121 |    2.009 |    28 |     0 |    0 |  -1.04e-07 |          0 |    -0.00263 |
| SCCD device Tight     |    4.078 |  18.308 /  105.31 |  33.312 /   88.43 |   48.720 |    1.859 |    12 |     0 |    0 |  -3.02e-08 |          0 |    -0.00024 |
| ACCD host             |    9.223 |  72.090 /  161.29 |  21.500 /   46.05 |  105.921 |    1.474 |   107 |     0 |    0 |    -0.0643 |          0 |      -0.897 |

  curated queries with a contact per pass: 2947611

### puffer-ball -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |   36.954 | 1496.290 / 4847.71 | 111.661 / 1043.15 | 1650.859 |    4.719 |     -0.0325 |       -0.0005 |    0 |
| SCCD host Tight       |   38.277 | 1495.055 / 4584.35 |  69.825 /  759.53 | 1605.304 |    2.993 |   -1.75e-05 |     -1.37e-06 |    0 |
| SCCD device Relaxed   |   27.141 |  54.663 /  998.64 |   8.838 /  763.56 |  160.951 |    1.216 |     -0.0584 |      -0.00251 |    0 |
| SCCD device Tight     |   27.085 |  54.978 /  991.15 |  11.841 /  775.02 |  161.786 |    1.272 |   -1.91e-05 |        -7e-07 |    0 |
| Scalable CCD device   |   10.517 | 858.211 / 2093.04 | 157.623 / 1586.40 | 1013.981 |    6.979 |     -0.0584 |      -0.00251 |    0 |

  steps with a known contact: 120

### puffer-ball -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |   36.954 | 1496.290 / 4847.71 | 117.065 / 1177.17 | 1660.480 |    5.028 | 27380 |     0 |    0 |    -0.0398 |          0 |      -0.975 |
| SCCD host Tight       |   38.277 | 1495.055 / 4584.35 | 128.516 / 1798.88 | 1776.293 |    5.998 |   536 |     0 |    0 |     -4e-05 |          0 |      -0.261 |
| SCCD device Relaxed   |   27.141 |  54.663 /  998.64 |  30.193 /  553.82 |  136.452 |    1.357 | 27380 |     0 |    0 |    -0.0416 |          0 |      -0.975 |
| SCCD device Tight     |   27.085 |  54.978 /  991.15 |  31.554 /  846.93 |  146.688 |    1.462 |   558 |     0 |    0 |  -4.02e-05 |          0 |      -0.261 |
| ACCD host             |   52.913 | 1475.890 / 5436.89 |  26.604 /  222.13 | 1564.583 |     1.19 | 27374 |     0 |    0 |     -0.105 |          0 |      -0.976 |

  curated queries with a contact per pass: 1486790

### rod-twist -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |    6.466 |   3.122 /   11.96 |   5.507 /   29.06 |   15.707 |    9.475 |     -0.0119 |     -5.57e-06 |    0 |
| SCCD host Tight       |    6.362 |   3.112 /   14.02 |   5.183 /  127.47 |   14.357 |    8.423 |   -7.32e-05 |      -3.3e-07 |    0 |
| SCCD device Relaxed   |    2.465 |   1.584 /    9.46 |   1.920 /   11.14 |    6.602 |    2.903 |     -0.0147 |     -9.92e-06 |    0 |
| SCCD device Tight     |    2.503 |   1.567 /    9.03 |   3.035 /   22.90 |    7.488 |    4.864 |   -7.26e-05 |     -2.15e-07 |    0 |
| Scalable CCD device   |    1.679 |  10.566 /   24.61 |  13.909 /  121.49 |   26.481 |    22.04 |     -0.0202 |         0.985 |   36 |

  steps with a known contact: 2481

### rod-twist -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    6.466 |   3.122 /   11.96 |   8.082 /  490.18 |   18.842 |    16.98 | 228608 |     0 |    0 |    -0.0268 |          0 |      -0.996 |
| SCCD host Tight       |    6.362 |   3.112 /   14.02 |   7.714 / 2002.34 |   17.760 |    21.28 |  1243 |     0 |    0 |  -0.000275 |          0 |      -0.993 |
| SCCD device Relaxed   |    2.465 |   1.584 /    9.46 |   3.164 /   34.62 |    7.643 |    5.079 | 244973 |     0 |    0 |    -0.0294 |          0 |      -0.996 |
| SCCD device Tight     |    2.503 |   1.567 /    9.03 |   4.145 /   81.66 |    8.484 |    7.091 |  1245 |     0 |    0 |  -0.000275 |          0 |      -0.993 |
| ACCD host             |    5.787 |   3.017 /   10.14 |   1.293 /   28.46 |   10.355 |    1.818 | 16708 |     0 |    0 |    -0.0232 |          0 |      -0.994 |

  curated queries with a contact per pass: 285431

  Time-of-impact error is signed: reported minus exact root, so negative is
  conservative and a positive value is a contact reported after the true one.
  Counts are per pass over the case list. The mesh path and the exact roots
  describe the same geometry here, because smesh is built with double geom_t;
  a float32 mesh answers about different coordinates and makes correct kernels
  look late.
```

## What the numbers say

**Conservativeness.** SCCD reports no missed contact and no time of impact after
the true root anywhere: zero late steps on all six scenes in every mode, and zero
late answers over the `5,522,383` curated contacts scored per pass. Its worst
signed step error is negative on every scene. ACCD is likewise never late and
never misses.

Scalable CCD's earliest time of impact lands after the true one on `272` of `394`
armadillo steps per pass, by up to `1.000` -- reporting no contact at all for
steps that have one -- on `36` of rod-twist's `2,481` and on `10` of
cloth-funnel's `363`. On cloth-ball, n-body-simulation and puffer-ball it is
conservative throughout. The cause is a race in its subdivision buffer's overflow
check, described in `../README.md` and diagnosed in `wip/DECISIONS.md` section 8.
The buffer is sized from the number of queries it is handed, so it bites when a
candidate list is short relative to the search it starts: the three scenes it
fails are the three with the fewest candidate pairs per step, and the three it
answers correctly are those with millions. Its answers vary between runs. Double
precision does not change this: the same run under a float32 mesh gave `271` on
armadillo.

**Against Scalable CCD, on cost.** The figures below are per-case medians. The
paper's tables sum whole scenes instead, which weights the heavy cases and puts
the same comparison at `3.0x` to `10.3x`.

At the median, end to end, SCCD on the device is `3.6x` faster on armadillo
(`5.20` against `18.51` ms per case), `3.5x` on cloth-ball (`14.23` against
`50.25`), `4.0x` on cloth-funnel (`4.38` against `17.63`), `4.7x` on n-body
(`53.58` against `254.29`), `6.3x` on puffer-ball (`161.79` against `1013.98`)
and `3.5x` on rod-twist (`7.49` against `26.48`). Per candidate its narrow phase
costs `23.8` against `55.4` ns on armadillo, `2.06` against `22.1` on cloth-ball,
`47.9` against `54.9` on cloth-funnel, `1.28` against `10.7` on n-body, `1.27`
against `6.98` on puffer-ball and `4.86` against `22.0` on rod-twist.

The broad phase is where the device gap is widest, and the run measures both of
our strategies so the algorithm and the implementation can be told apart. Summed
over the benchmark our device broad phase costs `30.3 s` against Scalable CCD's
`313.8 s`, and on puffer-ball `10.4 s` against `227.2 s`. Our own sweep sits
between the two at `111.5 s`, and on n-body it is actually the slower of the
two implementations -- `7.4 s` against Scalable CCD's `6.3` -- where our cell
list takes the scene at `1.9`.

**Against ACCD, per collision pair, the cost goes the other way.** Over the same
candidates on the same 72 threads, ACCD is the cheaper narrow phase on every
scene: `6.6` against SCCD Tight's `32.6` ns per candidate on armadillo, `1.44`
against `9.17` on cloth-ball, `14.4` against `28.9` on cloth-funnel, `1.47`
against `9.08` on n-body, `1.22` against `6.07` on puffer-ball and `1.84` against
`21.3` on rod-twist -- two to twelve times. Conservative advancement is cheap per
pair precisely because it stops at a bound instead of isolating a root.

**Tightness separates them by two to six orders of magnitude.** Per pair, SCCD
Tight's median error is `-3.4e-6` on armadillo, `-2.7e-7` on cloth-ball, `-6.9e-6`
on cloth-funnel, `-3.0e-8` on n-body, `-4.0e-5` on puffer-ball and `-2.8e-4` on
rod-twist; ACCD's is `-2.9e-2`, `-6.6e-2`, `-8.6e-3`, `-6.4e-2`, `-1.1e-1` and
`-2.3e-2`. That is `conservative_rescaling = 0.9` doing what it is for, and it
shows in the false-positive counts too: `444` against `24` on armadillo and
`27,374` against `536` on puffer-ball. Within SCCD, Relaxed trades the same way
against Tight -- `232` false positives against `24` on armadillo, `687` against
`15` on cloth-funnel, `228,608` against `1,243` on rod-twist, and on puffer-ball
`27,380`, which puts ACCD at its shipped rescaling about where our looser
acceptance test sits.

**cloth-funnel starts in contact.** `3,570` of the scene's `7,552` curated
queries have an exact root of `0`: the primitives already touch when the step
begins. ACCD reports `0` for them and says so ("Initial distance 0 ≤ d_min=0"),
which is conservative. Counting its broad-phase candidates too, that warning
fires about `7,000` times per pass and fills the scene's `.accd.err` files. It is
also why `worst early` reaches `-1` there: a query whose root is late in the step
but whose primitives touch at the start.

## What this run does not establish

One machine, three repeats. Every scene of the benchmark is covered, at every
prepared case, so the gaps left are in hardware and in repetition.

The two host narrow phases are threaded by different runtimes, SCCD's by OpenMP
and ACCD's by oneTBB, both over 72 threads of one Grace. That is the closest
available comparison, not an identical one.

The Scalable CCD result rests on a diagnosis rather than only on symptoms: a copy
of the library whose overflow check is made race-free on the host, run over the
same cases in double precision, finds all `54,664` curated contacts of
armadillo's first `400` cases with nothing missed and nothing late. So the buffer
race accounts for every miss and every late answer, and the library's search is
otherwise conservative. That copy is diagnosis only; the tables above measure the
library as published.
