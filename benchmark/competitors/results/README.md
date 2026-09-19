# Comparison: SCCD, Scalable CCD and Additive CCD

`compare-gh200-full-2026-09-19.csv.gz`, produced by `../sweep_comparison.sh` and
tabulated by `../compare_table.py` (`gunzip -k` it first; the raw CSV is 20 MB).
One GH200 module per job -- one Hopper, and one Grace bound to 72 CPUs -- over
whole scenes: every prepared case of armadillo-rollers (781), cloth-ball (79),
cloth-funnel (577), n-body-simulation (146), puffer-ball (240) and rod-twist
(4,571), `6,394` in all, three repeats each, as 66 resumable Slurm chunks.

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

### n-body-simulation -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |   10.021 |  71.595 /  160.36 |  56.538 /   99.30 |  135.141 |    2.798 |           0 |             0 |    0 |
| SCCD host Tight       |   10.582 |  71.449 /  170.30 |  44.001 /   80.43 |  121.829 |    2.266 |           0 |             0 |    0 |
| SCCD device Relaxed   |    6.942 |  52.256 /  104.34 |   2.919 /    5.29 |   62.476 |   0.1864 |           0 |             0 |    0 |
| SCCD device Tight     |    6.874 |  53.527 /  104.63 |   2.913 /    4.42 |   63.528 |   0.1854 |           0 |             0 |    0 |
| Scalable CCD device   |    3.007 |  39.981 /   54.35 | 222.156 /  470.49 |  260.500 |     10.7 |           0 |             0 |    0 |
| Scalable CCD host     |    1.916 |  53.291 /   89.36 |   0.000 /    0.00 |   55.181 |        - |           - |             - |    - |

  steps with a known contact: 74

### n-body-simulation -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |   10.021 |  71.595 /  160.36 | 150.401 /  306.75 |  230.759 |    6.993 |     9 |     0 |    0 |  -1.82e-08 |          0 |    -0.00127 |
| SCCD host Tight       |   10.582 |  71.449 /  170.30 | 140.621 /  320.61 |  241.023 |    9.076 |    11 |     0 |    0 |  -3.02e-08 |          0 |    -0.00024 |
| SCCD device Relaxed   |    6.942 |  52.256 /  104.34 |  36.667 /  130.37 |   96.829 |    2.068 |    28 |     0 |    0 |  -1.04e-07 |          0 |    -0.00263 |
| SCCD device Tight     |    6.874 |  53.527 /  104.63 |  34.395 /   87.75 |   95.467 |    1.936 |    12 |     0 |    0 |  -3.02e-08 |          0 |    -0.00024 |
| ACCD host             |   10.303 |  71.557 /  151.35 |  20.933 /   44.78 |  105.386 |    1.468 |   107 |     0 |    0 |    -0.0643 |          0 |      -0.897 |

  curated queries with a contact per pass: 2947611

### puffer-ball -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |   34.696 | 1448.995 / 4348.33 | 111.825 / 3734.16 | 1593.714 |    5.399 |     -0.0325 |       -0.0005 |    0 |
| SCCD host Tight       |   35.607 | 1472.220 / 4411.97 |  69.182 /  645.97 | 1581.523 |    2.993 |   -1.76e-05 |     -1.37e-06 |    0 |
| SCCD device Relaxed   |   49.488 | 191.803 /  938.98 |   2.900 /  188.09 |  243.261 |   0.1204 |     -0.0584 |      -0.00251 |    0 |
| SCCD device Tight     |   49.706 | 193.499 /  951.18 |   5.100 /   38.37 |  246.353 |   0.1995 |   -1.85e-05 |      -7.6e-07 |    0 |
| Scalable CCD device   |   10.764 | 856.337 / 2090.40 | 157.686 / 1561.35 | 1014.892 |    6.986 |     -0.0584 |      -0.00251 |    0 |
| Scalable CCD host     |   11.091 | 4632.390 /14896.80 |   0.000 /    0.00 | 4644.581 |        - |           - |             - |    - |

  steps with a known contact: 120

### puffer-ball -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |   34.696 | 1448.995 / 4348.33 | 115.887 / 3336.12 | 1601.138 |    5.612 | 27380 |     0 |    0 |    -0.0398 |          0 |      -0.975 |
| SCCD host Tight       |   35.607 | 1472.220 / 4411.97 | 131.416 / 2021.25 | 1765.613 |    6.072 |   536 |     0 |    0 |     -4e-05 |          0 |      -0.261 |
| SCCD device Relaxed   |   49.488 | 191.803 /  938.98 |  33.160 /  475.57 |  271.319 |    1.447 | 27380 |     0 |    0 |     -0.042 |          0 |      -0.975 |
| SCCD device Tight     |   49.706 | 193.499 /  951.18 |  32.599 /  656.20 |  273.630 |    1.521 |   558 |     0 |    0 |  -4.03e-05 |          0 |      -0.261 |
| ACCD host             |   33.612 | 1481.095 / 4732.60 |  27.346 /  228.13 | 1567.887 |    1.219 | 27374 |     0 |    0 |     -0.105 |          0 |      -0.976 |

  curated queries with a contact per pass: 1486790

### rod-twist -- earliest time of impact

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |     err med |     err worst | late |
|-----------------------|----------|------------------|------------------|----------|----------|-------------|---------------|------|
| SCCD host Relaxed     |    6.244 |   3.133 /   12.46 |   5.469 /   28.62 |   15.587 |    9.427 |      -0.012 |     -5.11e-06 |    0 |
| SCCD host Tight       |    6.482 |   3.079 /   10.22 |   5.204 /  126.97 |   14.724 |    8.425 |   -7.32e-05 |      -3.3e-07 |    0 |
| SCCD device Relaxed   |    4.009 |   1.783 /    9.97 |   1.845 /   58.52 |    7.392 |    2.316 |     -0.0146 |     -9.92e-06 |    0 |
| SCCD device Tight     |    4.020 |   1.851 /   24.82 |   2.938 /  142.51 |    8.588 |    4.373 |   -7.23e-05 |     -2.15e-07 |    0 |
| Scalable CCD device   |    1.691 |  10.438 /   30.86 |  13.858 /  137.57 |   26.363 |    22.07 |     -0.0203 |         0.969 |   33 |
| Scalable CCD host     |    1.166 |  14.001 /   30.37 |   0.000 /    0.00 |   15.198 |        - |           - |             - |    - |

  steps with a known contact: 2481

### rod-twist -- per collision pair

| library               |  prep ms |         broad ms |        narrow ms |    total |     ns/q |    fp |    fn | late |    err med | worst late | worst early |
|-----------------------|----------|------------------|------------------|----------|----------|-------|-------|------|------------|------------|-------------|
| SCCD host Relaxed     |    6.244 |   3.133 /   12.46 |   8.061 /  492.85 |   18.624 |    16.93 | 228608 |     0 |    0 |    -0.0268 |          0 |      -0.996 |
| SCCD host Tight       |    6.482 |   3.079 /   10.22 |   7.733 / 1978.78 |   18.372 |    21.29 |  1243 |     0 |    0 |  -0.000275 |          0 |      -0.993 |
| SCCD device Relaxed   |    4.009 |   1.783 /    9.97 |   3.153 /   34.72 |    8.758 |    5.055 | 244973 |     0 |    0 |    -0.0294 |          0 |      -0.996 |
| SCCD device Tight     |    4.020 |   1.851 /   24.82 |   4.142 /   80.14 |    9.713 |    7.083 |  1245 |     0 |    0 |  -0.000275 |          0 |      -0.993 |
| ACCD host             |    7.183 |   3.038 /   10.06 |   1.304 /   26.72 |   11.429 |    1.841 | 16708 |     0 |    0 |    -0.0232 |          0 |      -0.994 |

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

Scalable CCD's earliest time of impact lands after the true one on `265` of `394`
armadillo steps per pass, by up to `0.999` -- reporting no contact at all for
steps that have one -- on `33` of rod-twist's `2,481` and on `8` of
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
the same comparison at `3.2x` to `4.5x`.

At the median, end to end, SCCD on the device
is `3.4x` faster on armadillo (`5.50` against `18.56` ms per case), `3.7x` on
cloth-ball (`13.47` against `50.11`), `4.3x` on cloth-funnel (`4.04` against
`17.49`), `4.1x` on n-body (`63.53` against `260.50`), `4.1x` on puffer-ball
(`246.35` against `1014.89`) and `3.1x` on rod-twist (`8.59` against `26.36`).
Per candidate its narrow phase costs `18.5` against `57.3` ns on armadillo, `1.04`
against `21.8` on cloth-ball, `33.0` against `56.3` on cloth-funnel, `0.19`
against `10.7` on n-body, `0.20` against `6.99` on puffer-ball and `4.37` against
`22.1` on rod-twist. On the host, where Scalable CCD has only a broad phase,
SCCD's prep plus query is `5.22` against `5.94` ms on armadillo, `21.0` against
`22.9` on cloth-ball, `4.02` against `4.71` on cloth-funnel and `9.56` against
`15.20` on rod-twist; the sweep leads on the other two, `55.2` against `82.0` on
n-body and `4644.6` against `1507.8` on puffer-ball, which is the same split the
broad-phase study reports between our own two strategies.

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
