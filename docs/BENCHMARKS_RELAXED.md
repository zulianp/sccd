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
python3 -m report benchmark/results/sweep-gh200-bp.csv /tmp/report \
        benchmark/results/oracle-gh200-all.csv \
        --modes=relaxed,device-relaxed --embed=docs/BENCHMARKS_RELAXED.md --check
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

| scene             | phase | mode          | queries checked | missed | late |
|-------------------|-------|---------------|----------------:|-------:|-----:|
| armadillo-rollers | EE    | Relaxed (GPU) |          98,757 |      0 |    0 |
| armadillo-rollers | EE    | Relaxed       |          98,757 |      0 |    0 |
| armadillo-rollers | VF    | Relaxed (GPU) |          32,102 |      0 |    0 |
| armadillo-rollers | VF    | Relaxed       |          32,102 |      0 |    0 |
| cloth-ball        | EE    | Relaxed (GPU) |         557,667 |      0 |    0 |
| cloth-ball        | EE    | Relaxed       |         557,667 |      0 |    0 |
| cloth-ball        | VF    | Relaxed (GPU) |         107,252 |      0 |    0 |
| cloth-ball        | VF    | Relaxed       |         107,252 |      0 |    0 |
| cloth-funnel      | EE    | Relaxed (GPU) |           6,249 |      0 |    0 |
| cloth-funnel      | EE    | Relaxed       |           6,249 |      0 |    0 |
| cloth-funnel      | VF    | Relaxed (GPU) |             524 |      0 |    0 |
| cloth-funnel      | VF    | Relaxed       |             524 |      0 |    0 |
| n-body            | EE    | Relaxed (GPU) |       2,399,741 |      0 |    0 |
| n-body            | EE    | Relaxed       |       2,399,741 |      0 |    0 |
| n-body            | VF    | Relaxed (GPU) |         547,870 |      0 |    0 |
| n-body            | VF    | Relaxed       |         547,870 |      0 |    0 |
| puffer-ball       | EE    | Relaxed (GPU) |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | Relaxed       |       1,187,257 |      0 |    0 |
| puffer-ball       | VF    | Relaxed (GPU) |         299,533 |      0 |    0 |
| puffer-ball       | VF    | Relaxed       |         299,533 |      0 |    0 |
| rod-twist         | EE    | Relaxed (GPU) |         245,025 |      0 |    0 |
| rod-twist         | EE    | Relaxed       |         245,025 |      0 |    0 |
| rod-twist         | VF    | Relaxed (GPU) |          40,406 |      0 |    0 |
| rod-twist         | VF    | Relaxed       |          40,406 |      0 |    0 |

TightInclusion is excluded from this table: it is the reference the queries were selected against, not a subject of it.

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end gate -->

<!-- sccd:begin conservativeness -->

| scene             | mode          |   queries | toi compared | late | false pos. | false neg. | mesh-path divergence |
|-------------------|---------------|----------:|-------------:|-----:|-----------:|-----------:|---------------------:|
| armadillo-rollers | Relaxed (GPU) |   131,261 |      130,681 |    0 |        322 |          0 |                    0 |
| armadillo-rollers | Relaxed       |   131,261 |      130,681 |    0 |        232 |          0 |                    0 |
| cloth-ball        | Relaxed (GPU) |   664,939 |      664,918 |    0 |          2 |          0 |                    0 |
| cloth-ball        | Relaxed       |   664,939 |      664,918 |    0 |          2 |          0 |                    0 |
| cloth-funnel      | Relaxed (GPU) |     7,544 |        6,766 |    0 |        741 |          0 |                    0 |
| cloth-funnel      | Relaxed       |     7,544 |        6,766 |    0 |        686 |          0 |                    0 |
| n-body            | Relaxed (GPU) | 2,921,071 |    2,920,965 |    0 |         28 |          0 |                    0 |
| n-body            | Relaxed       | 2,921,071 |    2,920,965 |    0 |          9 |          0 |                    0 |
| puffer-ball       | Relaxed (GPU) | 1,513,963 |    1,486,587 |    0 |     27,374 |          0 |                    0 |
| puffer-ball       | Relaxed       | 1,513,963 |    1,486,587 |    0 |     27,374 |          0 |                    0 |
| rod-twist         | Relaxed (GPU) |   546,081 |      283,655 |    0 |    243,717 |          0 |                    0 |
| rod-twist         | Relaxed       |   546,081 |      283,655 |    0 |    227,441 |          0 |                    0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness. The last column is not part of the gate. It counts cases where the earliest-impact answer computed over the *mesh* is later than the earliest exact root of the *curated queries*, which are two separately stored geometries: the mesh is read from PLY, the queries are exact dyadic rationals. It is a measure of the agreement between those two inputs rather than of the kernel, and with the mesh stored in double it is zero everywhere.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end conservativeness -->

False positives are where the looseness shows: `Relaxed` reports more of them
than `Tight` on every scene, by roughly two orders of magnitude on rod-twist.
They cost work, never safety.

## 3. Accuracy

<!-- sccd:begin earliness-ref -->

| scene             | phase | mode           | median earliness | worst case |
|-------------------|-------|----------------|-----------------:|-----------:|
| armadillo-rollers | EE    | Relaxed        |         3.47e-05 |   7.37e-01 |
| armadillo-rollers | EE    | Relaxed (GPU)  |         7.47e-05 |   7.37e-01 |
| armadillo-rollers | EE    | TightInclusion |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | VF    | Relaxed        |         6.44e-05 |   9.41e-01 |
| armadillo-rollers | VF    | Relaxed (GPU)  |         1.06e-04 |   9.42e-01 |
| armadillo-rollers | VF    | TightInclusion |         3.74e-06 |   9.39e-03 |
| cloth-ball        | EE    | Relaxed        |         3.13e-07 |   5.40e-04 |
| cloth-ball        | EE    | Relaxed (GPU)  |         1.13e-06 |   9.22e-04 |
| cloth-ball        | EE    | TightInclusion |         1.05e-07 |   1.04e-04 |
| cloth-ball        | VF    | Relaxed        |         2.96e-07 |   7.86e-05 |
| cloth-ball        | VF    | Relaxed (GPU)  |         1.09e-06 |   2.15e-04 |
| cloth-ball        | VF    | TightInclusion |         9.76e-08 |   2.88e-05 |
| cloth-funnel      | EE    | Relaxed        |                0 |   1.00e+00 |
| cloth-funnel      | EE    | Relaxed (GPU)  |                0 |   1.00e+00 |
| cloth-funnel      | EE    | TightInclusion |                0 |   1.00e+00 |
| cloth-funnel      | VF    | Relaxed        |         1.56e-02 |   9.97e-01 |
| cloth-funnel      | VF    | Relaxed (GPU)  |         2.82e-02 |   9.97e-01 |
| cloth-funnel      | VF    | TightInclusion |         3.33e-04 |   2.67e-02 |
| n-body            | EE    | Relaxed        |         1.79e-08 |   1.27e-03 |
| n-body            | EE    | Relaxed (GPU)  |         1.10e-07 |   2.51e-03 |
| n-body            | EE    | TightInclusion |         1.14e-08 |   1.89e-04 |
| n-body            | VF    | Relaxed        |         1.84e-08 |   1.53e-04 |
| n-body            | VF    | Relaxed (GPU)  |         1.03e-07 |   2.64e-04 |
| n-body            | VF    | TightInclusion |         1.17e-08 |   1.87e-05 |
| puffer-ball       | EE    | Relaxed        |         3.57e-02 |   9.61e-01 |
| puffer-ball       | EE    | Relaxed (GPU)  |         3.64e-02 |   9.62e-01 |
| puffer-ball       | EE    | TightInclusion |         3.54e-05 |   2.61e-01 |
| puffer-ball       | VF    | Relaxed        |         4.68e-02 |   9.75e-01 |
| puffer-ball       | VF    | Relaxed (GPU)  |         4.60e-02 |   9.75e-01 |
| puffer-ball       | VF    | TightInclusion |         4.58e-05 |   5.04e-02 |
| rod-twist         | EE    | Relaxed        |         2.75e-02 |   9.93e-01 |
| rod-twist         | EE    | Relaxed (GPU)  |         2.94e-02 |   9.93e-01 |
| rod-twist         | EE    | TightInclusion |         2.47e-04 |   9.85e-01 |
| rod-twist         | VF    | Relaxed        |         1.41e-02 |   9.96e-01 |
| rod-twist         | VF    | Relaxed (GPU)  |         1.64e-02 |   9.96e-01 |
| rod-twist         | VF    | TightInclusion |         1.28e-04 |   9.93e-01 |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end earliness-ref -->

The worst case reaches 1.0 on cloth-funnel, puffer-ball and rod-twist — a step
reported at its start when the true contact is at its end. TightInclusion's
worst case is the same on those scenes, so the ceiling belongs to the geometry;
what `Relaxed` gives up is the **median**, and it gives up about 17× of it.

## 4. Against TightInclusion

<!-- sccd:begin reference -->

| scene             | phase | mode           |   queries |      hits |        time ms |     vs. TI |
|-------------------|-------|----------------|----------:|----------:|---------------:|-----------:|
| armadillo-rollers | EE    | Relaxed (GPU)  |    99,104 |    98,933 |    2335 / 2349 |       4.5× |
| armadillo-rollers | EE    | Relaxed        |    99,104 |    98,895 |    3377 / 3384 |       3.1× |
| armadillo-rollers | EE    | TightInclusion |    99,104 |    98,761 |  10499 / 10604 | 1.0× (ref) |
| armadillo-rollers | VF    | Relaxed (GPU)  |    32,337 |    32,248 |    3469 / 3476 |       2.2× |
| armadillo-rollers | VF    | Relaxed        |    32,337 |    32,196 |    2219 / 2242 |       3.4× |
| armadillo-rollers | VF    | TightInclusion |    32,337 |    32,122 |    7609 / 7686 | 1.0× (ref) |
| cloth-ball        | EE    | Relaxed (GPU)  |   557,683 |   557,669 |      645 / 647 |       2.4× |
| cloth-ball        | EE    | Relaxed        |   557,683 |   557,669 |      402 / 404 |       3.8× |
| cloth-ball        | EE    | TightInclusion |   557,683 |   557,668 |    1523 / 1529 | 1.0× (ref) |
| cloth-ball        | VF    | Relaxed (GPU)  |   107,257 |   107,252 |    2212 / 2367 |       0.6× |
| cloth-ball        | VF    | Relaxed        |   107,257 |   107,252 |      302 / 307 |       4.2× |
| cloth-ball        | VF    | TightInclusion |   107,257 |   107,252 |    1255 / 1267 | 1.0× (ref) |
| cloth-funnel      | EE    | Relaxed (GPU)  |     6,751 |     6,734 |      768 / 793 |       4.1× |
| cloth-funnel      | EE    | Relaxed        |     6,751 |     6,700 |    1111 / 1149 |       2.8× |
| cloth-funnel      | EE    | TightInclusion |     6,751 |     6,259 |    3130 / 3137 | 1.0× (ref) |
| cloth-funnel      | VF    | Relaxed (GPU)  |       801 |       781 |    2124 / 2599 |       0.6× |
| cloth-funnel      | VF    | Relaxed        |       801 |       760 |      402 / 410 |       3.3× |
| cloth-funnel      | VF    | TightInclusion |       801 |       529 |    1321 / 1321 | 1.0× (ref) |
| n-body            | EE    | Relaxed (GPU)  | 2,399,812 | 2,399,762 |    1106 / 1351 |       2.4× |
| n-body            | EE    | Relaxed        | 2,399,812 | 2,399,747 |    6576 / 6626 |       0.4× |
| n-body            | EE    | TightInclusion | 2,399,812 | 2,399,746 |    2633 / 2657 | 1.0× (ref) |
| n-body            | VF    | Relaxed (GPU)  |   547,907 |   547,877 |    2883 / 2900 |       0.5× |
| n-body            | VF    | Relaxed        |   547,907 |   547,873 |      423 / 424 |       3.7× |
| n-body            | VF    | TightInclusion |   547,907 |   547,873 |    1577 / 1587 | 1.0× (ref) |
| puffer-ball       | EE    | Relaxed (GPU)  | 1,206,952 | 1,206,951 |      346 / 346 |       7.6× |
| puffer-ball       | EE    | Relaxed        | 1,206,952 | 1,206,951 |      252 / 254 |      10.5× |
| puffer-ball       | EE    | TightInclusion | 1,206,952 | 1,187,650 |    2642 / 2662 | 1.0× (ref) |
| puffer-ball       | VF    | Relaxed (GPU)  |   307,220 |   307,219 |    2518 / 2549 |       0.6× |
| puffer-ball       | VF    | Relaxed        |   307,220 |   307,219 |      209 / 209 |       7.4× |
| puffer-ball       | VF    | TightInclusion |   307,220 |   299,676 |    1544 / 1556 | 1.0× (ref) |
| rod-twist         | EE    | Relaxed (GPU)  |   492,120 |   474,114 |    2898 / 4493 |      28.6× |
| rod-twist         | EE    | Relaxed        |   492,120 |   458,641 |   8100 / 22112 |      10.2× |
| rod-twist         | EE    | TightInclusion |   492,120 |   246,132 | 82987 / 249504 | 1.0× (ref) |
| rod-twist         | VF    | Relaxed (GPU)  |    57,088 |    56,290 |    3618 / 4208 |       6.5× |
| rod-twist         | VF    | Relaxed        |    57,088 |    55,398 |    2648 / 2974 |       8.9× |
| rod-twist         | VF    | TightInclusion |    57,088 |    40,542 |  23528 / 43830 | 1.0× (ref) |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end reference -->

## 5. CPU against GPU

<!-- sccd:begin processor -->

| scene             |  CPU ms | GPU ms | total | broad | narrow |
|-------------------|--------:|-------:|------:|------:|-------:|
| armadillo-rollers |  12,188 |  4,442 | 2.74x | 3.96x |  0.50x |
| cloth-ball        |   2,531 |  1,318 | 1.92x | 2.54x |  0.91x |
| cloth-funnel      |   8,422 |  2,718 | 3.10x | 4.03x |  0.68x |
| n-body            |  13,164 |  7,004 | 1.88x | 5.42x |  0.74x |
| puffer-ball       |  58,281 | 82,220 | 0.71x | 1.00x |  0.32x |
| rod-twist         | 109,057 | 42,687 | 2.55x | 4.26x |  1.37x |

A ratio is the host median over the device median, so 2.0 means the device takes half the time. Ratios below 1.0 are the cases where the host wins and are the ones worth reading.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end processor -->

## 6. Timing

<!-- sccd:begin timing -->

| scene             | mode          | cases |         pairs | rep |           prep ms |          broad ms |       earliest ms |       per-pair ms |            total ms |
|-------------------|---------------|------:|--------------:|----:|------------------:|------------------:|------------------:|------------------:|--------------------:|
| armadillo-rollers | Relaxed (GPU) |   779 |    85,015,700 |   3 |   2069.6 / 2071.3 |     838.6 / 841.4 |   1533.7 / 1536.7 |   2324.8 / 2333.5 |     4393.4 / 4446.4 |
| armadillo-rollers | Relaxed       |   779 |    85,015,700 |   3 |   8106.5 / 8135.7 |   3322.1 / 3354.3 |     759.8 / 771.9 |   2942.3 / 2959.0 |   12188.4 / 12262.0 |
| cloth-ball        | Relaxed (GPU) |    78 |   175,057,694 |   3 |     525.1 / 526.1 |     367.1 / 369.8 |     425.6 / 428.9 |     920.2 / 935.3 |     1316.7 / 1324.8 |
| cloth-ball        | Relaxed       |    78 |   175,057,694 |   3 |   1212.5 / 1218.4 |     932.9 / 937.5 |     385.5 / 387.8 |     789.9 / 802.9 |     2530.8 / 2543.6 |
| cloth-funnel      | Relaxed (GPU) |   575 |    25,141,046 |   3 |   1358.4 / 1376.9 |     540.9 / 545.7 |     818.3 / 837.9 |     727.4 / 728.4 |     2694.1 / 2760.5 |
| cloth-funnel      | Relaxed       |   575 |    25,141,046 |   3 |   5687.6 / 5764.2 |   2181.2 / 2191.1 |     553.7 / 568.0 |   1026.8 / 1038.2 |     8422.4 / 8523.4 |
| n-body            | Relaxed (GPU) |   145 | 2,233,533,498 |   3 |     967.8 / 978.7 |   1326.4 / 1330.1 |   4709.9 / 4717.8 |   4335.3 / 4455.1 |     6983.6 / 7015.0 |
| n-body            | Relaxed       |   145 | 2,233,533,498 |   3 |   2482.9 / 2492.2 |   7183.8 / 7189.1 |   3497.6 / 3554.2 | 10002.6 / 10064.4 |   13164.3 / 13228.7 |
| puffer-ball       | Relaxed (GPU) |   239 | 7,298,095,145 |   3 | 16853.8 / 16990.3 | 27239.6 / 27256.6 | 38127.0 / 38131.4 | 10667.8 / 10669.8 |   81789.3 / 82356.9 |
| puffer-ball       | Relaxed       |   239 | 7,298,095,145 |   3 | 18928.3 / 19590.4 | 27289.2 / 27356.9 | 12063.1 / 12188.8 | 14353.1 / 14739.7 |   58172.0 / 58854.7 |
| rod-twist         | Relaxed (GPU) |  4559 | 3,861,367,133 |   3 | 26191.6 / 26628.1 |   6780.7 / 7013.1 |   9715.2 / 9723.6 | 19749.6 / 19844.2 |   42687.5 / 43286.8 |
| rod-twist         | Relaxed       |  4559 | 3,861,367,133 |   3 | 66897.1 / 66993.6 | 28867.0 / 28918.5 | 13293.1 / 13373.1 | 39888.5 / 39892.8 | 109044.9 / 109285.3 |

Each cell is the median over repeats and the slowest of them. A difference smaller than the gap between the two does not separate two modes and is not reported as a ratio anywhere in this document. Which mode is faster depends on the output mode as well as the scene, so the two are given side by side rather than one standing for the other.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end timing -->

<!-- sccd:begin throughput -->

| scene             | mode          | broad Mpair/s | narrow Mpair/s |
|-------------------|---------------|--------------:|---------------:|
| armadillo-rollers | Relaxed (GPU) |         101.4 |           55.4 |
| armadillo-rollers | Relaxed       |          25.6 |          111.9 |
| cloth-ball        | Relaxed (GPU) |         476.9 |          411.3 |
| cloth-ball        | Relaxed       |         187.7 |          454.1 |
| cloth-funnel      | Relaxed (GPU) |          46.5 |           30.7 |
| cloth-funnel      | Relaxed       |          11.5 |           45.4 |
| n-body            | Relaxed (GPU) |        1683.9 |          474.2 |
| n-body            | Relaxed       |         310.9 |          638.6 |
| puffer-ball       | Relaxed (GPU) |         267.9 |          191.4 |
| puffer-ball       | Relaxed       |         267.4 |          605.0 |
| rod-twist         | Relaxed (GPU) |         569.5 |          397.5 |
| rod-twist         | Relaxed       |         133.8 |          290.5 |

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end throughput -->

## 7. Figures

Scenes across the columns, one quantity per row, distributions over cases as box
plots on logarithmic axes.

![Per-case results over the six scenes](figures/relaxed-results-grid.png)

**Figure 1.** Per-case distributions over the six scenes: broad-phase time,
narrow-phase time, and error against the exact root. `*` marks a parallel CPU
mode, `†` a GPU one.

![Runtime split by phase](figures/relaxed-runtime-breakdown.png)

**Figure 2.** Runtime split into *prep*, *broad* and *narrow*, one bar per
configuration: `*` is a parallel CPU mode, `†` a GPU one.

![Time-of-impact error against the symbolic ground truth](figures/relaxed-toi-error.png)

**Figure 3.** Error against the exact symbolic roots, log axis. One-sided by
construction: every value is on the safe side. `*` marks a parallel CPU mode,
`†` a GPU one; the device curve is dashed where the two coincide.

## 8. Provenance

<!-- sccd:begin provenance -->

- Timings: `benchmark/results/sweep-gh200-bp.csv`, 6375 cases over 6 scenes, 3 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200-all.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
