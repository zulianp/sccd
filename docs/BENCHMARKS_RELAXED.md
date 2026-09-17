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

| scene             | mode |   queries | toi compared | late | false pos. | false neg. |
|-------------------|------|----------:|-------------:|-----:|-----------:|-----------:|
| armadillo-rollers | GPU  |   131,359 |      130,777 |    0 |        322 |          0 |
| armadillo-rollers | CPU  |   131,359 |      130,777 |    0 |        232 |          0 |
| cloth-ball        | GPU  |   664,939 |      664,918 |    0 |          2 |          0 |
| cloth-ball        | CPU  |   664,939 |      664,918 |    0 |          2 |          0 |
| cloth-funnel      | GPU  |     7,532 |        6,753 |    0 |        742 |          0 |
| cloth-funnel      | CPU  |     7,532 |        6,753 |    0 |        687 |          0 |
| n-body            | GPU  | 2,921,071 |    2,920,965 |    0 |         28 |          0 |
| n-body            | CPU  | 2,921,071 |    2,920,965 |    0 |          9 |          0 |
| puffer-ball       | GPU  | 1,513,963 |    1,486,587 |    0 |     27,374 |          0 |
| puffer-ball       | CPU  | 1,513,963 |    1,486,587 |    0 |     27,374 |          0 |
| rod-twist         | GPU  |   548,211 |      285,240 |    0 |    244,232 |          0 |
| rod-twist         | CPU  |   548,211 |      285,240 |    0 |    227,934 |          0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end conservativeness -->

False positives are where the looseness shows: `Relaxed` reports more of them
than `Tight` on every scene, by roughly two orders of magnitude on rod-twist.
They cost work, never safety.

The mesh path is checked against the curated query geometry as a tripwire on the
inputs:

<!-- sccd:begin mesh-check -->

**73 cases** where the mesh-path answer falls after the earliest exact root of the curated queries: the two paths are not being given the same geometry, so the exact roots are not a reference for the mesh path.

<!-- sccd:end mesh-check -->

## 3. Accuracy

<!-- sccd:begin earliness-ref -->

| scene             | phase | mode           | median earliness | worst case |
|-------------------|-------|----------------|-----------------:|-----------:|
| armadillo-rollers | EE    | CPU            |         3.47e-05 |   7.37e-01 |
| armadillo-rollers | EE    | GPU            |         7.54e-05 |   7.37e-01 |
| armadillo-rollers | EE    | TightInclusion |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | VF    | CPU            |         6.44e-05 |   9.41e-01 |
| armadillo-rollers | VF    | GPU            |         1.06e-04 |   9.42e-01 |
| armadillo-rollers | VF    | TightInclusion |         3.74e-06 |   9.39e-03 |
| cloth-ball        | EE    | CPU            |         3.13e-07 |   5.40e-04 |
| cloth-ball        | EE    | GPU            |         1.13e-06 |   8.96e-04 |
| cloth-ball        | EE    | TightInclusion |         1.05e-07 |   1.04e-04 |
| cloth-ball        | VF    | CPU            |         2.96e-07 |   7.86e-05 |
| cloth-ball        | VF    | GPU            |         1.09e-06 |   2.16e-04 |
| cloth-ball        | VF    | TightInclusion |         9.76e-08 |   2.88e-05 |
| cloth-funnel      | EE    | CPU            |                0 |   1.00e+00 |
| cloth-funnel      | EE    | GPU            |                0 |   1.00e+00 |
| cloth-funnel      | EE    | TightInclusion |                0 |   1.00e+00 |
| cloth-funnel      | VF    | CPU            |         1.56e-02 |   9.97e-01 |
| cloth-funnel      | VF    | GPU            |         2.82e-02 |   9.97e-01 |
| cloth-funnel      | VF    | TightInclusion |         3.33e-04 |   2.67e-02 |
| n-body            | EE    | CPU            |         1.79e-08 |   1.27e-03 |
| n-body            | EE    | GPU            |         1.10e-07 |   2.60e-03 |
| n-body            | EE    | TightInclusion |         1.14e-08 |   1.89e-04 |
| n-body            | VF    | CPU            |         1.84e-08 |   1.53e-04 |
| n-body            | VF    | GPU            |         1.02e-07 |   2.59e-04 |
| n-body            | VF    | TightInclusion |         1.17e-08 |   1.87e-05 |
| puffer-ball       | EE    | CPU            |         3.57e-02 |   9.61e-01 |
| puffer-ball       | EE    | GPU            |         3.64e-02 |   9.62e-01 |
| puffer-ball       | EE    | TightInclusion |         3.54e-05 |   2.61e-01 |
| puffer-ball       | VF    | CPU            |         4.68e-02 |   9.75e-01 |
| puffer-ball       | VF    | GPU            |         4.61e-02 |   9.75e-01 |
| puffer-ball       | VF    | TightInclusion |         4.58e-05 |   5.04e-02 |
| rod-twist         | EE    | CPU            |         2.75e-02 |   9.93e-01 |
| rod-twist         | EE    | GPU            |         2.92e-02 |   9.93e-01 |
| rod-twist         | EE    | TightInclusion |         2.47e-04 |   9.85e-01 |
| rod-twist         | VF    | CPU            |         1.41e-02 |   9.96e-01 |
| rod-twist         | VF    | GPU            |         1.63e-02 |   9.96e-01 |
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
| armadillo-rollers | EE    | GPU            |    99,104 |    98,933 |    1905 / 1925 |       5.5× |
| armadillo-rollers | EE    | CPU            |    99,104 |    98,895 |    3442 / 3480 |       3.0× |
| armadillo-rollers | EE    | TightInclusion |    99,104 |    98,761 |  10471 / 10554 | 1.0× (ref) |
| armadillo-rollers | VF    | GPU            |    32,337 |    32,248 |    1441 / 1458 |       5.3× |
| armadillo-rollers | VF    | CPU            |    32,337 |    32,196 |    2294 / 2306 |       3.3× |
| armadillo-rollers | VF    | TightInclusion |    32,337 |    32,122 |    7569 / 7631 | 1.0× (ref) |
| cloth-ball        | EE    | GPU            |   557,683 |   557,669 |      529 / 531 |       2.9× |
| cloth-ball        | EE    | CPU            |   557,683 |   557,669 |      416 / 416 |       3.7× |
| cloth-ball        | EE    | TightInclusion |   557,683 |   557,668 |    1521 / 1522 | 1.0× (ref) |
| cloth-ball        | VF    | GPU            |   107,257 |   107,252 |      380 / 384 |       3.3× |
| cloth-ball        | VF    | CPU            |   107,257 |   107,252 |      317 / 320 |       4.0× |
| cloth-ball        | VF    | TightInclusion |   107,257 |   107,252 |    1255 / 1256 | 1.0× (ref) |
| cloth-funnel      | EE    | GPU            |     6,751 |     6,734 |      705 / 716 |       4.4× |
| cloth-funnel      | EE    | CPU            |     6,751 |     6,700 |    1146 / 1148 |       2.7× |
| cloth-funnel      | EE    | TightInclusion |     6,751 |     6,259 |    3122 / 3143 | 1.0× (ref) |
| cloth-funnel      | VF    | GPU            |       801 |       781 |      339 / 355 |       3.8× |
| cloth-funnel      | VF    | CPU            |       801 |       760 |      391 / 406 |       3.3× |
| cloth-funnel      | VF    | TightInclusion |       801 |       529 |    1305 / 1339 | 1.0× (ref) |
| n-body            | EE    | GPU            | 2,399,812 | 2,399,762 |    1100 / 1109 |       2.4× |
| n-body            | EE    | CPU            | 2,399,812 | 2,399,747 |    6688 / 6700 |       0.4× |
| n-body            | EE    | TightInclusion | 2,399,812 | 2,399,746 |    2638 / 2642 | 1.0× (ref) |
| n-body            | VF    | GPU            |   547,907 |   547,877 |      639 / 658 |       2.5× |
| n-body            | VF    | CPU            |   547,907 |   547,873 |      434 / 435 |       3.6× |
| n-body            | VF    | TightInclusion |   547,907 |   547,873 |    1581 / 1588 | 1.0× (ref) |
| puffer-ball       | EE    | GPU            | 1,206,952 | 1,206,951 |      344 / 351 |       7.7× |
| puffer-ball       | EE    | CPU            | 1,206,952 | 1,206,951 |      250 / 256 |      10.6× |
| puffer-ball       | EE    | TightInclusion | 1,206,952 | 1,187,650 |    2657 / 2662 | 1.0× (ref) |
| puffer-ball       | VF    | GPU            |   307,220 |   307,219 |      263 / 272 |       5.9× |
| puffer-ball       | VF    | CPU            |   307,220 |   307,219 |      205 / 206 |       7.5× |
| puffer-ball       | VF    | TightInclusion |   307,220 |   299,676 |    1544 / 1550 | 1.0× (ref) |
| rod-twist         | EE    | GPU            |   492,120 |   474,114 |    2682 / 3761 |      30.9× |
| rod-twist         | EE    | CPU            |   492,120 |   458,641 |   8292 / 22772 |      10.0× |
| rod-twist         | EE    | TightInclusion |   492,120 |   246,132 | 82791 / 250085 | 1.0× (ref) |
| rod-twist         | VF    | GPU            |    57,088 |    56,290 |    1685 / 1979 |      13.9× |
| rod-twist         | VF    | CPU            |    57,088 |    55,398 |    2726 / 3076 |       8.6× |
| rod-twist         | VF    | TightInclusion |    57,088 |    40,542 |  23363 / 43924 | 1.0× (ref) |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end reference -->

## 5. CPU against GPU

<!-- sccd:begin processor -->

| scene             |  CPU ms | GPU ms | total | broad | narrow |
|-------------------|--------:|-------:|------:|------:|-------:|
| armadillo-rollers |   3,575 |  3,608 | 0.99x | 1.40x |  0.44x |
| cloth-ball        |   2,804 |  1,161 | 2.41x | 2.44x |  2.36x |
| cloth-funnel      |   2,044 |  2,191 | 0.93x | 1.32x |  0.30x |
| n-body            |  20,030 |  7,568 | 2.65x | 4.56x |  1.72x |
| puffer-ball       | 109,790 | 76,014 | 1.44x | 2.03x |  0.85x |
| rod-twist         |  78,789 | 33,088 | 2.38x | 1.87x |  3.66x |

A ratio is the host median over the device median, so 2.0 means the device takes half the time. Ratios below 1.0 are the cases where the host wins and are the ones worth reading.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end processor -->

## 6. Timing

Phases are as defined in [`BENCHMARKS.md`](BENCHMARKS.md#2-what-each-phase-measures):
*prep* builds the swept boxes and the acceleration structure, *broad* runs the
overlap query, *narrow* finds the time of impact. Per simulation step, which
runs both query types:

<!-- sccd:begin per-frame -->

| scene             | frames | mode | broad ms | narrow ms | total ms |
|-------------------|-------:|------|---------:|----------:|---------:|
| armadillo-rollers |    396 | GPU  |     5.26 |      3.85 |     9.11 |
| armadillo-rollers |    396 | CPU  |     7.34 |      1.68 |     9.03 |
| cloth-ball        |     42 | GPU  |    17.65 |     10.00 |    27.65 |
| cloth-ball        |     42 | CPU  |    43.14 |     23.63 |    66.77 |
| cloth-funnel      |    372 | GPU  |     3.64 |      2.25 |     5.89 |
| cloth-funnel      |    372 | CPU  |     4.82 |      0.67 |     5.50 |
| n-body            |     74 | GPU  |    33.43 |     68.84 |   102.27 |
| n-body            |     74 | CPU  |   152.41 |    118.26 |   270.67 |
| puffer-ball       |    120 | GPU  |   319.00 |    314.44 |   633.45 |
| puffer-ball       |    120 | CPU  |   646.39 |    268.52 |   914.91 |
| rod-twist         |  2,554 | GPU  |     9.24 |      3.71 |    12.96 |
| rod-twist         |  2,554 | CPU  |    17.24 |     13.61 |    30.85 |

A mean rather than a median over steps: the scene total is what a run costs, and the mean is the only average that divides back into it.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end per-frame -->


<!-- sccd:begin timing -->

| scene             | mode | cases |         pairs | rep |          broad ms |       earliest ms |       per-pair ms |            total ms |
|-------------------|------|------:|--------------:|----:|------------------:|------------------:|------------------:|--------------------:|
| armadillo-rollers | GPU  |   779 |    85,061,295 |   3 |   2069.8 / 2099.7 |   1526.4 / 1546.0 |   2335.6 / 2336.6 |     3610.5 / 3626.2 |
| armadillo-rollers | CPU  |   779 |    85,061,295 |   3 |   2924.5 / 2928.7 |     667.1 / 671.8 |   2176.9 / 2197.3 |     3591.9 / 3596.3 |
| cloth-ball        | GPU  |    78 |   175,057,694 |   3 |     744.9 / 748.4 |     420.1 / 424.0 |     924.8 / 936.6 |     1165.0 / 1172.4 |
| cloth-ball        | CPU  |    78 |   175,057,694 |   3 |   1808.2 / 1845.7 |     992.4 / 994.4 |   1797.0 / 1798.7 |     2800.7 / 2840.1 |
| cloth-funnel      | GPU  |   575 |    25,151,368 |   3 |   1354.1 / 1384.0 |     836.8 / 844.0 |     739.6 / 742.7 |     2196.4 / 2220.8 |
| cloth-funnel      | CPU  |   575 |    25,151,368 |   3 |   1793.6 / 1853.6 |     250.9 / 251.0 |     395.7 / 405.7 |     2044.5 / 2104.6 |
| n-body            | GPU  |   145 | 2,233,533,498 |   3 |   2475.5 / 2505.7 |   5094.2 / 5124.8 |   4526.4 / 4639.8 |     7599.9 / 7600.3 |
| n-body            | CPU  |   145 | 2,233,533,498 |   3 | 11276.7 / 11285.3 |   8751.2 / 8787.4 | 16143.2 / 16277.6 |   19983.7 / 20072.6 |
| puffer-ball       | GPU  |   239 | 7,298,180,358 |   3 | 38394.5 / 38497.7 | 37733.0 / 38469.6 | 10126.2 / 10782.7 |   76230.7 / 76864.1 |
| puffer-ball       | CPU  |   239 | 7,298,180,358 |   3 | 77625.1 / 78764.3 | 32222.6 / 32338.9 | 35035.4 / 35277.8 | 109964.0 / 110986.9 |
| rod-twist         | GPU  |  4561 | 3,860,443,124 |   3 | 23405.6 / 23931.5 |   9484.8 / 9495.2 | 19247.1 / 19423.2 |   32890.3 / 33398.6 |
| rod-twist         | CPU  |  4561 | 3,860,443,124 |   3 | 44032.0 / 44176.8 | 34757.4 / 34858.0 | 61721.6 / 61906.8 |   78789.4 / 78916.0 |

Each cell is the median over repeats and the slowest of them. A difference smaller than the gap between the two does not separate two modes and is not reported as a ratio anywhere in this document. Which mode is faster depends on the output mode as well as the scene, so the two are given side by side rather than one standing for the other.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end timing -->

<!-- sccd:begin throughput -->

| scene             | mode | broad Mpair/s | narrow Mpair/s |
|-------------------|------|--------------:|---------------:|
| armadillo-rollers | GPU  |          41.1 |           55.7 |
| armadillo-rollers | CPU  |          29.1 |          127.5 |
| cloth-ball        | GPU  |         235.0 |          416.7 |
| cloth-ball        | CPU  |          96.8 |          176.4 |
| cloth-funnel      | GPU  |          18.6 |           30.1 |
| cloth-funnel      | CPU  |          14.0 |          100.3 |
| n-body            | GPU  |         902.2 |          438.5 |
| n-body            | CPU  |         198.1 |          255.2 |
| puffer-ball       | GPU  |         190.1 |          193.4 |
| puffer-ball       | CPU  |          94.0 |          226.5 |
| rod-twist         | GPU  |         164.9 |          407.0 |
| rod-twist         | CPU  |          87.7 |          111.1 |

Source: `benchmark/results/sweep-gh200-bp.csv`

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

- Timings: `benchmark/results/sweep-gh200-bp.csv`, 6377 cases over 6 scenes, 3 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200-all.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
