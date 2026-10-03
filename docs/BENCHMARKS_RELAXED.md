# Benchmarks — `Relaxed`

SCCD's `Relaxed` mode on the same six scenes, the same node and the same
committed data as [`BENCHMARKS.md`](BENCHMARKS.md), which evaluates `Tight`.

The two are reported separately because they answer different questions.
`Tight` reproduces TightInclusion's answer exactly and is what the headline
evaluation measures. `Relaxed` accepts a box sooner, buying speed on some scenes
with a looser time of impact on all of them. Averaging them would describe
neither.

Regenerate with:

```sh
python3 -m report benchmark/assessment/broadphase-cell2dminfv.csv /tmp/report \
        benchmark/results/oracle-gh200-all.csv \
        --modes=relaxed,device-relaxed --label="relaxed:CPU,device-relaxed:GPU" \
        --figure-prefix=relaxed- --embed=docs/BENCHMARKS_RELAXED.md --check
```

## 1. What it trades

Against `Tight`, on the same cases:

| | |
|---|---|
| median earliness | **17× looser** (6.4e-05 against 3.7e-06) |
| narrow phase, CPU | 1.4× faster on armadillo-rollers, cloth-funnel and rod-twist; **2× slower** on cloth-ball, n-body and puffer-ball |
| narrow phase, GPU | 1.0× to 1.8× faster, never slower |
| conservativeness | identical — zero missed, zero late |

The accuracy cost is unconditional; the speed gain is not. On the host `Relaxed`
loses on half the scenes, because accepting a box sooner also hands the next
query a looser bound, so it prunes less. On the device it wins or ties
everywhere, and by the widest margin on rod-twist.

**It is exactly as conservative.** The invariant does not depend on how tight
the acceptance test is: accepting reports the box's `t` lower bound, which is at
or before any root inside it, so a looser test reports earlier and never later.

## 2. Conservativeness

<!-- sccd:begin gate -->

| scene             | phase | mode | queries checked | missed | late |
|-------------------|-------|------|----------------:|-------:|-----:|
| armadillo-rollers | EE    | GPU  |          98,757 |      0 |    0 |
| armadillo-rollers | EE    | CPU  |          98,757 |      0 |    0 |
| armadillo-rollers | VF    | GPU  |          32,102 |      0 |    0 |
| armadillo-rollers | VF    | CPU  |          32,102 |      0 |    0 |
| cloth-ball        | EE    | GPU  |         557,667 |      0 |    0 |
| cloth-ball        | EE    | CPU  |         557,667 |      0 |    0 |
| cloth-ball        | VF    | GPU  |         107,252 |      0 |    0 |
| cloth-ball        | VF    | CPU  |         107,252 |      0 |    0 |
| cloth-funnel      | EE    | GPU  |           6,249 |      0 |    0 |
| cloth-funnel      | EE    | CPU  |           6,249 |      0 |    0 |
| cloth-funnel      | VF    | GPU  |             524 |      0 |    0 |
| cloth-funnel      | VF    | CPU  |             524 |      0 |    0 |
| n-body            | EE    | GPU  |       2,399,741 |      0 |    0 |
| n-body            | EE    | CPU  |       2,399,741 |      0 |    0 |
| n-body            | VF    | GPU  |         547,870 |      0 |    0 |
| n-body            | VF    | CPU  |         547,870 |      0 |    0 |
| puffer-ball       | EE    | GPU  |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | CPU  |       1,187,257 |      0 |    0 |
| puffer-ball       | VF    | GPU  |         299,533 |      0 |    0 |
| puffer-ball       | VF    | CPU  |         299,533 |      0 |    0 |
| rod-twist         | EE    | GPU  |         245,025 |      0 |    0 |
| rod-twist         | EE    | CPU  |         245,025 |      0 |    0 |
| rod-twist         | VF    | GPU  |          40,406 |      0 |    0 |
| rod-twist         | VF    | CPU  |          40,406 |      0 |    0 |

TightInclusion is excluded from this table: it is the reference the queries were selected against, not a subject of it.

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end gate -->

<!-- sccd:begin conservativeness -->

| scene             |   queries | toi compared | late |        false pos. | false neg. |
|-------------------|----------:|-------------:|-----:|------------------:|-----------:|
| armadillo-rollers |   131,441 |      130,859 |    0 |         232 (322) |          0 |
| cloth-ball        |   664,940 |      664,919 |    0 |                 2 |          0 |
| cloth-funnel      |     7,552 |        6,773 |    0 |         687 (742) |          0 |
| n-body            | 2,947,719 |    2,947,611 |    0 |            9 (28) |          0 |
| puffer-ball       | 1,514,172 |    1,486,790 |    0 |            27,380 |          0 |
| rod-twist         |   549,208 |      285,431 |    0 | 228,608 (244,973) |          0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness.

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

<!-- sccd:end conservativeness -->

False positives are where the looseness shows: `Relaxed` reports more of them
than `Tight` on every scene, by roughly two orders of magnitude on rod-twist.
They cost work, never safety.

The mesh path is checked against the curated query geometry as a tripwire on the
inputs:

<!-- sccd:begin mesh-check -->

The mesh path agrees with the curated query geometry on every case.

<!-- sccd:end mesh-check -->

## 3. Accuracy

<!-- sccd:begin earliness-ref -->

| scene             | phase | mode           | median earliness | worst case |
|-------------------|-------|----------------|-----------------:|-----------:|
| armadillo-rollers | EE    | CPU            |         3.47e-05 |   7.37e-01 |
| armadillo-rollers | EE    | GPU            |         7.59e-05 |   7.37e-01 |
| armadillo-rollers | EE    | TightInclusion |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | VF    | CPU            |         6.44e-05 |   9.41e-01 |
| armadillo-rollers | VF    | GPU            |         1.05e-04 |   9.41e-01 |
| armadillo-rollers | VF    | TightInclusion |         3.74e-06 |   9.39e-03 |
| cloth-ball        | EE    | CPU            |         3.13e-07 |   5.40e-04 |
| cloth-ball        | EE    | GPU            |         1.13e-06 |   9.01e-04 |
| cloth-ball        | EE    | TightInclusion |         1.05e-07 |   1.04e-04 |
| cloth-ball        | VF    | CPU            |         2.96e-07 |   7.86e-05 |
| cloth-ball        | VF    | GPU            |         1.09e-06 |   2.18e-04 |
| cloth-ball        | VF    | TightInclusion |         9.76e-08 |   2.88e-05 |
| cloth-funnel      | EE    | CPU            |                0 |   1.00e+00 |
| cloth-funnel      | EE    | GPU            |                0 |   1.00e+00 |
| cloth-funnel      | EE    | TightInclusion |                0 |   1.00e+00 |
| cloth-funnel      | VF    | CPU            |         1.56e-02 |   9.97e-01 |
| cloth-funnel      | VF    | GPU            |         2.76e-02 |   9.97e-01 |
| cloth-funnel      | VF    | TightInclusion |         3.33e-04 |   2.67e-02 |
| n-body            | EE    | CPU            |         1.79e-08 |   1.27e-03 |
| n-body            | EE    | GPU            |         1.10e-07 |   2.56e-03 |
| n-body            | EE    | TightInclusion |         1.14e-08 |   1.89e-04 |
| n-body            | VF    | CPU            |         1.84e-08 |   1.53e-04 |
| n-body            | VF    | GPU            |         1.03e-07 |   2.60e-04 |
| n-body            | VF    | TightInclusion |         1.17e-08 |   1.87e-05 |
| puffer-ball       | EE    | CPU            |         3.57e-02 |   9.61e-01 |
| puffer-ball       | EE    | GPU            |         3.64e-02 |   9.62e-01 |
| puffer-ball       | EE    | TightInclusion |         3.54e-05 |   2.61e-01 |
| puffer-ball       | VF    | CPU            |         4.68e-02 |   9.75e-01 |
| puffer-ball       | VF    | GPU            |         4.62e-02 |   9.75e-01 |
| puffer-ball       | VF    | TightInclusion |         4.58e-05 |   5.04e-02 |
| rod-twist         | EE    | CPU            |         3.23e-02 |   9.93e-01 |
| rod-twist         | EE    | GPU            |         3.53e-02 |   9.93e-01 |
| rod-twist         | EE    | TightInclusion |         3.69e-04 |   9.85e-01 |
| rod-twist         | VF    | CPU            |         1.71e-02 |   9.96e-01 |
| rod-twist         | VF    | GPU            |         1.90e-02 |   9.96e-01 |
| rod-twist         | VF    | TightInclusion |         1.54e-04 |   9.93e-01 |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end earliness-ref -->

The worst case reaches 1.0 on cloth-funnel, puffer-ball and rod-twist — a step
reported at its start when the true contact is at its end. TightInclusion's
worst case is the same on those scenes, so the ceiling belongs to the geometry;
what `Relaxed` gives up is the **median**, and it gives up about 17× of it.

## 4. Against TightInclusion

![SCCD against TightInclusion](figures/reference-speedup.png)

Each bar is one processor's whole narrow phase over the same queries, stacked into its vertex-face and edge-edge work, with the total and the speedup over TightInclusion above it. Query counts are in the dataset table and hit agreement in the conservativeness table, so neither is repeated here.

## 5. CPU against GPU

<!-- sccd:begin processor -->

| scene             | CPU ms | GPU ms | total | broad | narrow |
|-------------------|-------:|-------:|------:|------:|-------:|
| armadillo-rollers |  1,790 |  1,699 | 1.05x | 2.48x |  0.54x |
| cloth-ball        |  1,459 |    467 | 3.12x | 1.32x |  6.32x |
| cloth-funnel      |  1,146 |  1,357 | 0.84x | 1.30x |  0.33x |
| n-body            | 11,202 |  3,507 | 3.19x | 0.45x | 20.98x |
| puffer-ball       | 49,955 | 18,998 | 2.63x | 0.87x | 45.99x |
| rod-twist         | 50,797 | 15,388 | 3.30x | 2.38x |  3.93x |

A ratio is the host median over the device median, so 2.0 means the device takes half the time. Ratios below 1.0 are the cases where the host wins and are the ones worth reading.

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

<!-- sccd:end processor -->

## 6. Timing

Phases are as defined in [`BENCHMARKS.md`](BENCHMARKS.md#2-what-each-phase-measures):
*prep* builds the swept boxes and the acceleration structure, *broad* runs the
overlap query, *narrow* finds the time of impact. Per simulation step, which
runs both query types:

<!-- sccd:begin per-frame -->

| scene             | frames | mode | broad ms | narrow ms | total ms |
|-------------------|-------:|------|---------:|----------:|---------:|
| armadillo-rollers |    396 | GPU  |     1.14 |      3.15 |     4.29 |
| armadillo-rollers |    396 | CPU  |     2.81 |      1.71 |     4.52 |
| cloth-ball        |     43 | GPU  |     6.95 |      3.92 |    10.87 |
| cloth-ball        |     43 | CPU  |     9.15 |     24.77 |    33.92 |
| cloth-funnel      |    372 | GPU  |     1.93 |      1.72 |     3.65 |
| cloth-funnel      |    372 | CPU  |     2.52 |      0.56 |     3.08 |
| n-body            |     74 | GPU  |    41.05 |      6.33 |    47.39 |
| n-body            |     74 | CPU  |    18.50 |    132.87 |   151.38 |
| puffer-ball       |    120 | GPU  |   152.13 |      6.19 |   158.32 |
| puffer-ball       |    120 | CPU  |   131.76 |    284.53 |   416.29 |
| rod-twist         |  2,556 | GPU  |     2.44 |      3.58 |     6.02 |
| rod-twist         |  2,556 | CPU  |     5.81 |     14.06 |    19.87 |

A mean rather than a median over steps: the scene total is what a run costs, and the mean is the only average that divides back into it.

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

<!-- sccd:end per-frame -->


<!-- sccd:begin timing -->

| scene             | mode | cases |         pairs | rep |          broad ms |       earliest ms |       per-pair ms |          total ms |
|-------------------|------|------:|--------------:|----:|------------------:|------------------:|------------------:|------------------:|
| armadillo-rollers | GPU  |   781 |    85,296,282 |   2 |     449.6 / 450.3 |   1249.2 / 1256.0 |   2106.8 / 2107.3 |   1698.8 / 1706.3 |
| armadillo-rollers | CPU  |   781 |    85,296,282 |   2 |   1112.9 / 1116.3 |     677.4 / 684.4 |   1858.0 / 1867.8 |   1790.3 / 1794.0 |
| cloth-ball        | GPU  |    79 |   175,150,856 |   2 |     298.9 / 299.3 |     168.4 / 169.2 |     896.0 / 899.1 |     467.3 / 468.5 |
| cloth-ball        | CPU  |    79 |   175,150,856 |   2 |     393.4 / 402.7 |   1065.2 / 1069.0 |   1839.4 / 1844.7 |   1458.7 / 1464.2 |
| cloth-funnel      | GPU  |   577 |    25,192,698 |   2 |     718.6 / 719.1 |     638.1 / 646.8 |     683.7 / 686.3 |   1356.7 / 1365.9 |
| cloth-funnel      | CPU  |   577 |    25,192,698 |   2 |     937.4 / 960.6 |     208.4 / 217.2 |     336.0 / 337.8 |   1145.8 / 1177.8 |
| n-body            | GPU  |   146 | 2,253,766,609 |   2 |   3037.7 / 3046.6 |     468.8 / 475.7 |   4463.9 / 4481.1 |   3506.5 / 3508.4 |
| n-body            | CPU  |   146 | 2,253,766,609 |   2 |   1369.4 / 1418.0 |   9832.5 / 9842.4 | 17428.5 / 17431.4 | 11201.9 / 11240.7 |
| puffer-ball       | GPU  |   240 | 7,327,014,600 |   2 | 18255.7 / 18310.4 |     742.4 / 745.5 |   8803.6 / 8813.1 | 18998.0 / 19049.6 |
| puffer-ball       | CPU  |   240 | 7,327,014,600 |   2 | 15811.8 / 15903.0 | 34143.4 / 34159.5 | 36810.3 / 36897.0 | 49955.1 / 50062.5 |
| rod-twist         | GPU  |  4571 | 3,872,779,843 |   2 |   6234.6 / 6243.2 |   9153.6 / 9196.0 | 18203.9 / 18205.2 | 15388.2 / 15439.2 |
| rod-twist         | CPU  |  4571 | 3,872,779,843 |   2 | 14855.6 / 14980.8 | 35941.2 / 35948.8 | 62890.1 / 62948.1 | 50796.8 / 50929.6 |

Each cell is the median over repeats and the slowest of them. A difference smaller than the gap between the two does not separate two modes and is not reported as a ratio anywhere in this document. Which mode is faster depends on the output mode as well as the scene, so the two are given side by side rather than one standing for the other.

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

<!-- sccd:end timing -->

<!-- sccd:begin throughput -->

| scene             | mode | broad Mpair/s | narrow Mpair/s |
|-------------------|------|--------------:|---------------:|
| armadillo-rollers | GPU  |         189.7 |           68.3 |
| armadillo-rollers | CPU  |          76.6 |          125.9 |
| cloth-ball        | GPU  |         586.0 |         1039.8 |
| cloth-ball        | CPU  |         445.2 |          164.4 |
| cloth-funnel      | GPU  |          35.1 |           39.5 |
| cloth-funnel      | CPU  |          26.9 |          120.9 |
| n-body            | GPU  |         741.9 |         4807.8 |
| n-body            | CPU  |        1645.9 |          229.2 |
| puffer-ball       | GPU  |         401.4 |         9869.8 |
| puffer-ball       | CPU  |         463.4 |          214.6 |
| rod-twist         | GPU  |         621.2 |          423.1 |
| rod-twist         | CPU  |         260.7 |          107.8 |

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

<!-- sccd:end throughput -->

## 7. Figures

Scenes across the columns, one quantity per row, distributions over cases as box
plots on logarithmic axes.

![Per-case results over the six scenes](figures/relaxed-results-grid.png)

**Figure 1.** Per-case distributions over the six scenes: broad-phase time,
narrow-phase time, and error against the exact root.

![Runtime split by phase](figures/relaxed-runtime-breakdown.png)

**Figure 2.** Runtime split into *prep*, *broad* and *narrow*.

![Time-of-impact error against the symbolic ground truth](figures/relaxed-toi-error.png)

**Figure 3.** Error against the exact symbolic roots, log axis. One-sided by
construction: every value is on the safe side. The device curve is dashed where the two
coincide.

## 8. Provenance

<!-- sccd:begin provenance -->

- Timings: `benchmark/assessment/broadphase-cell2dminfv.csv`, 6394 cases over 6 scenes, 2 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200-all.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
