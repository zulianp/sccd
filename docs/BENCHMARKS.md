# Benchmarks

SCCD's `Tight` mode, evaluated against TightInclusion on the NYU CCD dataset:
six scenes, 6,394 simulation steps, on one GH200 node.

`Relaxed`, the mode that trades accuracy for speed, is evaluated separately in
[`BENCHMARKS_RELAXED.md`](BENCHMARKS_RELAXED.md). The two answer different
questions and are not averaged together.

Every number here is generated from a committed CSV by a committed script:

```sh
python3 -m report benchmark/results/sweep-gh200-bp.csv /tmp/report \
        benchmark/results/oracle-gh200-all.csv \
        --modes=tight,device-tight --embed=docs/BENCHMARKS.md --check
```

`--check` returns non-zero if this document is not what the data reproduces.
Drop it to refresh. To measure from scratch, `benchmark/scripts/prepare_data.sh`
then `benchmark/scripts/sweep.sh`; the first refuses to finish on an incomplete
oracle, the second is resumable at chunk granularity.

## 1. What is measured

SCCD is a **conservative** detector, so one property is an invariant and the
rest are trade-offs:

- A reported time of impact is **at or before** the true one. Later would let a
  simulation step through the contact.
- A collision that exists is reported.
- Reporting a collision that does not exist, or a time of impact earlier than
  the true one, is acceptable — it costs work and step size, not safety.

The conservativeness tables are checked against the dataset's exact symbolic
roots, not against another implementation. Repeated measurements are reported as
`median / slowest`, and distributions over cases carry their worst case: a
median describes the typical case, and this search is paid for in the worst one.

## 2. Platform

One GH200 node on CSCS Alps — Grace (288 threads) and Hopper sm_90 in the same
allocation — built with `prgenv-gnu/24.11:v2`, `-O3`,
`CMAKE_CUDA_ARCHITECTURES=90`, at the shipped default search parameters. Mesh
geometry is stored in double (`SMESH_GEOM_TYPE=float64`); the root finder
computes in double regardless.

## 3. Dataset

The NYU CCD benchmark (archive 2451/74508): real simulation frames, each case a
step with a curated query set whose coordinates are exact dyadic rationals and
whose roots were computed symbolically.

<!-- sccd:begin dataset -->

| scene             | cases | candidate pairs/step |   queries | with a root |
|-------------------|------:|---------------------:|----------:|------------:|
| armadillo-rollers |   781 |              109,214 |   131,441 |     130,859 |
| cloth-ball        |    79 |            2,217,099 |   664,940 |     664,919 |
| cloth-funnel      |   577 |               43,661 |     7,552 |       6,773 |
| n-body            |   146 |           15,436,757 | 2,947,719 |   2,947,611 |
| puffer-ball       |   240 |           30,529,227 | 1,514,172 |   1,486,790 |
| rod-twist         | 4,571 |              847,250 |   549,208 |     285,431 |

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end dataset -->

Coverage is complete on all six scenes: a query carries a root exactly when
`mma_bool` marks it as colliding, and `benchmark/verify_oracle.py` gates on that
agreement query for query. The proportion varies by scene because it is the
fraction of curated queries that collide — rod-twist's 52% and cloth-ball's
99.997% are both complete.

## 4. Conservativeness

The invariant, checked against the exact symbolic roots for every scene on both
processors. **Zero missed collisions and zero late times of impact.**

<!-- sccd:begin gate -->

| scene             | phase | mode        | queries checked | missed | late |
|-------------------|-------|-------------|----------------:|-------:|-----:|
| armadillo-rollers | EE    | Tight (GPU) |          98,757 |      0 |    0 |
| armadillo-rollers | EE    | Tight       |          98,757 |      0 |    0 |
| armadillo-rollers | VF    | Tight (GPU) |          32,102 |      0 |    0 |
| armadillo-rollers | VF    | Tight       |          32,102 |      0 |    0 |
| cloth-ball        | EE    | Tight (GPU) |         557,667 |      0 |    0 |
| cloth-ball        | EE    | Tight       |         557,667 |      0 |    0 |
| cloth-ball        | VF    | Tight (GPU) |         107,252 |      0 |    0 |
| cloth-ball        | VF    | Tight       |         107,252 |      0 |    0 |
| cloth-funnel      | EE    | Tight (GPU) |           6,249 |      0 |    0 |
| cloth-funnel      | EE    | Tight       |           6,249 |      0 |    0 |
| cloth-funnel      | VF    | Tight (GPU) |             524 |      0 |    0 |
| cloth-funnel      | VF    | Tight       |             524 |      0 |    0 |
| n-body            | EE    | Tight (GPU) |       2,399,741 |      0 |    0 |
| n-body            | EE    | Tight       |       2,399,741 |      0 |    0 |
| n-body            | VF    | Tight (GPU) |         547,870 |      0 |    0 |
| n-body            | VF    | Tight       |         547,870 |      0 |    0 |
| puffer-ball       | EE    | Tight (GPU) |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | Tight       |       1,187,257 |      0 |    0 |
| puffer-ball       | VF    | Tight (GPU) |         299,533 |      0 |    0 |
| puffer-ball       | VF    | Tight       |         299,533 |      0 |    0 |
| rod-twist         | EE    | Tight (GPU) |         245,025 |      0 |    0 |
| rod-twist         | EE    | Tight       |         245,025 |      0 |    0 |
| rod-twist         | VF    | Tight (GPU) |          40,406 |      0 |    0 |
| rod-twist         | VF    | Tight       |          40,406 |      0 |    0 |

TightInclusion is excluded from this table: it is the reference the queries were selected against, not a subject of it.

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end gate -->

<!-- sccd:begin conservativeness -->

| scene             | mode        |   queries | toi compared | late | false pos. | false neg. | mesh-path divergence |
|-------------------|-------------|----------:|-------------:|-----:|-----------:|-----------:|---------------------:|
| armadillo-rollers | Tight (GPU) |   131,441 |      130,859 |    0 |         24 |          0 |                    0 |
| armadillo-rollers | Tight       |   131,441 |      130,859 |    0 |         24 |          0 |                    0 |
| cloth-ball        | Tight (GPU) |   664,940 |      664,919 |    0 |          1 |          0 |                    0 |
| cloth-ball        | Tight       |   664,940 |      664,919 |    0 |          1 |          0 |                    0 |
| cloth-funnel      | Tight (GPU) |     7,552 |        6,773 |    0 |         15 |          0 |                    0 |
| cloth-funnel      | Tight       |     7,552 |        6,773 |    0 |         15 |          0 |                    0 |
| n-body            | Tight (GPU) | 2,947,719 |    2,947,611 |    0 |         12 |          0 |                    0 |
| n-body            | Tight       | 2,947,719 |    2,947,611 |    0 |         11 |          0 |                    0 |
| puffer-ball       | Tight (GPU) | 1,514,172 |    1,486,790 |    0 |        558 |          0 |                    0 |
| puffer-ball       | Tight       | 1,514,172 |    1,486,790 |    0 |        536 |          0 |                    0 |
| rod-twist         | Tight (GPU) |   549,208 |      285,431 |    0 |      1,245 |          0 |                    0 |
| rod-twist         | Tight       |   549,208 |      285,431 |    0 |      1,243 |          0 |                    0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness. The last column is not part of the gate. It counts cases where the earliest-impact answer computed over the *mesh* is later than the earliest exact root of the *curated queries*, which are two separately stored geometries: the mesh is read from PLY, the queries are exact dyadic rationals. It is a measure of the agreement between those two inputs rather than of the kernel, and with the mesh stored in double it is zero everywhere.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end conservativeness -->

The last column is not part of the gate. It counts cases where the earliest
impact computed over the **mesh** lands after the earliest exact root of the
**curated queries** — two separately stored copies of the same frames. With the
mesh stored in double it is zero everywhere.

## 5. Against TightInclusion

TightInclusion is the reference implementation of a certified conservative
narrow phase. Both sides are given the same task and the same machine: SCCD at
`ToiOutput::PerPair` and TightInclusion unbounded over the whole step, both
scheduled through `sccd::parallel_for_br_dynamic`, with TightInclusion used
exactly as released and run depth-first, which is its faster configuration when
there is no bound to exploit.

<!-- sccd:begin reference -->

| scene             | phase | mode           |   queries |      hits |        time ms |     vs. TI |
|-------------------|-------|----------------|----------:|----------:|---------------:|-----------:|
| armadillo-rollers | EE    | Tight (GPU)    |    99,104 |    98,761 |    2908 / 2910 |       3.6× |
| armadillo-rollers | EE    | Tight          |    99,104 |    98,761 |    6474 / 6502 |       1.6× |
| armadillo-rollers | EE    | TightInclusion |    99,104 |    98,761 |  10499 / 10604 | 1.0× (ref) |
| armadillo-rollers | VF    | Tight (GPU)    |    32,337 |    32,122 |    1987 / 2011 |       3.8× |
| armadillo-rollers | VF    | Tight          |    32,337 |    32,122 |    2717 / 2743 |       2.8× |
| armadillo-rollers | VF    | TightInclusion |    32,337 |    32,122 |    7609 / 7686 | 1.0× (ref) |
| cloth-ball        | EE    | Tight (GPU)    |   557,683 |   557,668 |      609 / 612 |       2.5× |
| cloth-ball        | EE    | Tight          |   557,683 |   557,668 |      955 / 961 |       1.6× |
| cloth-ball        | EE    | TightInclusion |   557,683 |   557,668 |    1523 / 1529 | 1.0× (ref) |
| cloth-ball        | VF    | Tight (GPU)    |   107,257 |   107,252 |      428 / 429 |       2.9× |
| cloth-ball        | VF    | Tight          |   107,257 |   107,252 |      408 / 414 |       3.1× |
| cloth-ball        | VF    | TightInclusion |   107,257 |   107,252 |    1255 / 1267 | 1.0× (ref) |
| cloth-funnel      | EE    | Tight (GPU)    |     6,751 |     6,259 |    1358 / 1370 |       2.3× |
| cloth-funnel      | EE    | Tight          |     6,751 |     6,259 |    1743 / 1746 |       1.8× |
| cloth-funnel      | EE    | TightInclusion |     6,751 |     6,259 |    3130 / 3137 | 1.0× (ref) |
| cloth-funnel      | VF    | Tight (GPU)    |       801 |       529 |      451 / 456 |       2.9× |
| cloth-funnel      | VF    | Tight          |       801 |       529 |      394 / 395 |       3.4× |
| cloth-funnel      | VF    | TightInclusion |       801 |       529 |    1321 / 1321 | 1.0× (ref) |
| n-body            | EE    | Tight (GPU)    | 2,399,812 | 2,399,746 |     890 / 1077 |       3.0× |
| n-body            | EE    | Tight          | 2,399,812 | 2,399,746 |    1176 / 1184 |       2.2× |
| n-body            | EE    | TightInclusion | 2,399,812 | 2,399,746 |    2633 / 2657 | 1.0× (ref) |
| n-body            | VF    | Tight (GPU)    |   547,907 |   547,874 |      628 / 790 |       2.5× |
| n-body            | VF    | Tight          |   547,907 |   547,873 |      507 / 511 |       3.1× |
| n-body            | VF    | TightInclusion |   547,907 |   547,873 |    1577 / 1587 | 1.0× (ref) |
| puffer-ball       | EE    | Tight (GPU)    | 1,206,952 | 1,187,650 |      941 / 943 |       2.8× |
| puffer-ball       | EE    | Tight          | 1,206,952 | 1,187,650 |    1614 / 1629 |       1.6× |
| puffer-ball       | EE    | TightInclusion | 1,206,952 | 1,187,650 |    2642 / 2662 | 1.0× (ref) |
| puffer-ball       | VF    | Tight (GPU)    |   307,220 |   299,698 |      634 / 636 |       2.4× |
| puffer-ball       | VF    | Tight          |   307,220 |   299,676 |      586 / 595 |       2.6× |
| puffer-ball       | VF    | TightInclusion |   307,220 |   299,676 |    1544 / 1556 | 1.0× (ref) |
| rod-twist         | EE    | Tight (GPU)    |   492,120 |   246,132 |    4913 / 8176 |      16.9× |
| rod-twist         | EE    | Tight          |   492,120 |   246,132 | 37118 / 121624 |       2.2× |
| rod-twist         | EE    | TightInclusion |   492,120 |   246,132 | 82987 / 249504 | 1.0× (ref) |
| rod-twist         | VF    | Tight (GPU)    |    57,088 |    40,544 |    3014 / 3236 |       7.8× |
| rod-twist         | VF    | Tight          |    57,088 |    40,542 |   6382 / 10301 |       3.7× |
| rod-twist         | VF    | TightInclusion |    57,088 |    40,542 |  23528 / 43830 | 1.0× (ref) |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end reference -->

### 5.1. Accuracy

TightInclusion is a subject here, not the reference: its own answer is a
conservative lower bound, so the question is how close all three get to the
exact roots.

<!-- sccd:begin earliness-ref -->

| scene             | phase | mode           | median earliness | worst case |
|-------------------|-------|----------------|-----------------:|-----------:|
| armadillo-rollers | EE    | Tight          |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | EE    | Tight (GPU)    |         3.08e-06 |   1.60e-02 |
| armadillo-rollers | EE    | TightInclusion |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | VF    | Tight          |         3.74e-06 |   9.39e-03 |
| armadillo-rollers | VF    | Tight (GPU)    |         3.73e-06 |   9.39e-03 |
| armadillo-rollers | VF    | TightInclusion |         3.74e-06 |   9.39e-03 |
| cloth-ball        | EE    | Tight          |         1.05e-07 |   1.04e-04 |
| cloth-ball        | EE    | Tight (GPU)    |         1.05e-07 |   1.04e-04 |
| cloth-ball        | EE    | TightInclusion |         1.05e-07 |   1.04e-04 |
| cloth-ball        | VF    | Tight          |         9.76e-08 |   2.88e-05 |
| cloth-ball        | VF    | Tight (GPU)    |         9.75e-08 |   2.88e-05 |
| cloth-ball        | VF    | TightInclusion |         9.76e-08 |   2.88e-05 |
| cloth-funnel      | EE    | Tight          |                0 |   1.00e+00 |
| cloth-funnel      | EE    | Tight (GPU)    |                0 |   1.00e+00 |
| cloth-funnel      | EE    | TightInclusion |                0 |   1.00e+00 |
| cloth-funnel      | VF    | Tight          |         3.33e-04 |   2.67e-02 |
| cloth-funnel      | VF    | Tight (GPU)    |         3.33e-04 |   2.67e-02 |
| cloth-funnel      | VF    | TightInclusion |         3.33e-04 |   2.67e-02 |
| n-body            | EE    | Tight          |         1.14e-08 |   1.89e-04 |
| n-body            | EE    | Tight (GPU)    |         1.14e-08 |   1.89e-04 |
| n-body            | EE    | TightInclusion |         1.14e-08 |   1.89e-04 |
| n-body            | VF    | Tight          |         1.17e-08 |   1.87e-05 |
| n-body            | VF    | Tight (GPU)    |         1.17e-08 |   1.87e-05 |
| n-body            | VF    | TightInclusion |         1.17e-08 |   1.87e-05 |
| puffer-ball       | EE    | Tight          |         3.54e-05 |   2.61e-01 |
| puffer-ball       | EE    | Tight (GPU)    |         3.54e-05 |   2.61e-01 |
| puffer-ball       | EE    | TightInclusion |         3.54e-05 |   2.61e-01 |
| puffer-ball       | VF    | Tight          |         4.58e-05 |   5.04e-02 |
| puffer-ball       | VF    | Tight (GPU)    |         4.56e-05 |   5.04e-02 |
| puffer-ball       | VF    | TightInclusion |         4.58e-05 |   5.04e-02 |
| rod-twist         | EE    | Tight          |         2.47e-04 |   9.85e-01 |
| rod-twist         | EE    | Tight (GPU)    |         2.47e-04 |   9.85e-01 |
| rod-twist         | EE    | TightInclusion |         2.47e-04 |   9.85e-01 |
| rod-twist         | VF    | Tight          |         1.28e-04 |   9.93e-01 |
| rod-twist         | VF    | Tight (GPU)    |         1.28e-04 |   9.93e-01 |
| rod-twist         | VF    | TightInclusion |         1.28e-04 |   9.93e-01 |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end earliness-ref -->

**`Tight` does not approximate the reference, it reproduces it.** On all twelve
scene-phases its median earliness and its worst case are identical to
TightInclusion's, to every digit. `Tight (GPU)` matches the worst case exactly
and the median to within a part in a thousand.

The worst cases belong to the problem rather than to SCCD: on cloth-funnel and
rod-twist the largest earliness approaches 1.0 for *every* implementation
including the reference, because a conservative search has no tighter answer
available on a grazing or already-touching configuration.

## 6. CPU against GPU

<!-- sccd:begin processor -->

| scene             |  CPU ms | GPU ms | total | broad | narrow |
|-------------------|--------:|-------:|------:|------:|-------:|
| armadillo-rollers |  12,457 |  4,398 | 2.83x | 4.10x |  0.60x |
| cloth-ball        |   2,390 |  1,341 | 1.78x | 2.56x |  0.49x |
| cloth-funnel      |   8,506 |  2,838 | 3.00x | 4.04x |  0.70x |
| n-body            |  11,570 |  7,037 | 1.64x | 5.51x |  0.37x |
| puffer-ball       |  52,283 | 83,155 | 0.63x | 0.99x |  0.16x |
| rod-twist         | 113,926 | 49,251 | 2.31x | 4.25x |  1.06x |

A ratio is the host median over the device median, so 2.0 means the device takes half the time. Ratios below 1.0 are the cases where the host wins and are the ones worth reading.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end processor -->

The split is consistent and mechanical. The **broad phase** suits a GPU — count,
prefix sum, scatter, no sequential window walk — and runs 2.6× to 5.5× faster on
five of six scenes. The **narrow phase** is a depth-first interval search with a
divergent stack, and runs 0.16× to 1.06× as fast: slower on five, faster only on
rod-twist. Because the broad phase is the larger share on the host, the GPU
still wins end to end on five scenes by 1.6× to 3.0×.

**puffer-ball inverts it**, at 0.63× overall. Its GPU narrow phase is 6× slower
than the host's and its broad phase gives nothing back (0.99×). It is the
largest scene in the set by candidate pairs, 30.5 million per step, so problem
size is not a rule for choosing a processor — measure the scene.

<!-- sccd:begin timing -->

| scene             | mode        | cases |         pairs | rep |           prep ms |          broad ms |       earliest ms |       per-pair ms |            total ms |
|-------------------|-------------|------:|--------------:|----:|------------------:|------------------:|------------------:|------------------:|--------------------:|
| armadillo-rollers | Tight (GPU) |   781 |    85,296,282 |   3 |   1848.8 / 1960.1 |     813.1 / 838.0 |   1735.9 / 1763.2 |   2650.0 / 2663.7 |     4417.8 / 4509.0 |
| armadillo-rollers | Tight       |   781 |    85,296,282 |   3 |   8088.8 / 8212.7 |   3331.1 / 3376.8 |   1036.6 / 1040.4 |   3565.9 / 3587.7 |   12456.6 / 12629.9 |
| cloth-ball        | Tight (GPU) |    79 |   175,150,856 |   3 |     530.4 / 533.7 |     369.1 / 371.3 |     441.2 / 442.2 |     883.1 / 889.4 |     1343.4 / 1343.8 |
| cloth-ball        | Tight       |    79 |   175,150,856 |   3 |   1230.0 / 1231.1 |     945.3 / 946.7 |     214.7 / 216.6 |   1206.8 / 1220.8 |     2388.0 / 2393.3 |
| cloth-funnel      | Tight (GPU) |   577 |    25,192,698 |   3 |   1353.5 / 1358.1 |     540.6 / 541.4 |     943.9 / 969.6 |   1071.1 / 1074.4 |     2836.1 / 2836.5 |
| cloth-funnel      | Tight       |   577 |    25,192,698 |   3 |   5661.3 / 5701.0 |   2185.8 / 2198.8 |     658.6 / 664.2 |   1298.2 / 1298.4 |     8524.3 / 8545.5 |
| n-body            | Tight (GPU) |   146 | 2,253,766,609 |   3 |   1004.8 / 1010.8 |   1311.7 / 1337.1 |   4720.7 / 4722.7 |   3976.9 / 4065.0 |     7043.2 / 7062.6 |
| n-body            | Tight       |   146 | 2,253,766,609 |   3 |   2607.7 / 2618.8 |   7223.9 / 7286.5 |   1738.0 / 1756.1 | 19258.3 / 19443.3 |   11555.4 / 11650.2 |
| puffer-ball       | Tight (GPU) |   240 | 7,327,014,600 |   3 | 16975.7 / 17124.0 | 27521.8 / 27584.2 | 38657.9 / 38794.1 | 10752.3 / 11037.2 |   83171.7 / 83354.0 |
| puffer-ball       | Tight       |   240 | 7,327,014,600 |   3 | 18898.1 / 19172.8 | 27335.1 / 27640.5 |   6049.4 / 6154.7 | 56908.1 / 56970.6 |   52374.1 / 52503.0 |
| rod-twist         | Tight (GPU) |  4571 | 3,872,779,843 |   3 | 25140.3 / 25800.2 |   6742.8 / 6908.7 | 17368.0 / 17514.2 | 27523.3 / 27625.2 |   49397.4 / 50076.9 |
| rod-twist         | Tight       |  4571 | 3,872,779,843 |   3 | 66765.3 / 67399.4 | 28680.1 / 28853.0 | 18480.3 / 18549.3 | 65174.0 / 65384.2 | 113887.7 / 114732.7 |

Each cell is the median over repeats and the slowest of them. A difference smaller than the gap between the two does not separate two modes and is not reported as a ratio anywhere in this document. Which mode is faster depends on the output mode as well as the scene, so the two are given side by side rather than one standing for the other.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end timing -->

<!-- sccd:begin throughput -->

| scene             | mode        | broad Mpair/s | narrow Mpair/s |
|-------------------|-------------|--------------:|---------------:|
| armadillo-rollers | Tight (GPU) |         104.9 |           49.1 |
| armadillo-rollers | Tight       |          25.6 |           82.3 |
| cloth-ball        | Tight (GPU) |         474.5 |          397.0 |
| cloth-ball        | Tight       |         185.3 |          815.7 |
| cloth-funnel      | Tight (GPU) |          46.6 |           26.7 |
| cloth-funnel      | Tight       |          11.5 |           38.3 |
| n-body            | Tight (GPU) |        1718.2 |          477.4 |
| n-body            | Tight       |         312.0 |         1296.8 |
| puffer-ball       | Tight (GPU) |         266.2 |          189.5 |
| puffer-ball       | Tight       |         268.0 |         1211.2 |
| rod-twist         | Tight (GPU) |         574.4 |          223.0 |
| rod-twist         | Tight       |         135.0 |          209.6 |

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end throughput -->

## 7. Broad phase

Two strategies produce the candidate pairs — a sweep over sorted intervals and a
cell list over a uniform grid. They return identical pair sets, so the choice is
purely about cost, and the shipped default races them per scene rather than
fixing a winner. Across 25,538 case-mode combinations the two report identical
candidate pair counts.

`SCCD_BROADPHASE` is read in the three host broad-phase functions and nowhere
else; the device steps have no branch on it, so only host rows are listed.

<!-- sccd:begin broadphase -->

| scene             | mode  | cell2d prep ms | sweep prep ms | cell2d broad ms | sweep broad ms | faster       |
|-------------------|-------|---------------:|--------------:|----------------:|---------------:|--------------|
| armadillo-rollers | Tight |           8089 |         26637 |            3331 |           3231 | sweep 1.03x  |
| cloth-ball        | Tight |           1230 |          3064 |             945 |            697 | sweep 1.36x  |
| cloth-funnel      | Tight |           5661 |         19460 |            2186 |           2340 | cell2d 1.07x |
| n-body            | Tight |           2608 |          5902 |            7224 |           3833 | sweep 1.88x  |
| puffer-ball       | Tight |          18898 |         27213 |           27335 |         110588 | cell2d 4.05x |
| rod-twist         | Tight |          66765 |        169005 |           28680 |          21551 | sweep 1.33x  |

`faster` names the winning strategy and by how much on the whole broad phase. A margin inside the run-to-run spread is reported as a tie rather than a winner.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end broadphase -->

**Neither wins.** The sweep takes armadillo-rollers, cloth-ball, n-body and
rod-twist by 1.03× to 1.88×; the cell list takes cloth-funnel by 1.06× and
puffer-ball by **4.0×** — 27 s against 109 s on the largest scene. A fixed choice
is wrong somewhere, and wrong by a factor of four at the worst point, which is
what makes racing them worth its two probe steps.

## 8. Scaling with element count

`sccd_refine_scaling` refines one surface repeatedly, quadrupling the element
count at each level, over two consecutive cloth-ball frames from 92,230 elements
to 23.6 million.

<!-- sccd:begin scaling -->

| mode                          | level |   elements | candidate pairs | prep ms | step ms | broad ms | narrow ms |    p |
|-------------------------------|------:|-----------:|----------------:|--------:|--------:|---------:|----------:|-----:|
| relaxed / host / tri / cell2d |     0 |     92,230 |          17,982 |    27.9 |     8.8 |     36.7 |       8.9 | 0.85 |
|                               |     1 |    368,920 |          87,771 |    50.8 |    14.6 |     65.4 |       6.2 |      |
|                               |     2 |  1,475,680 |         379,818 |   200.1 |    45.5 |    245.5 |       2.0 |      |
|                               |     3 |  5,902,720 |       1,573,741 |   852.6 |   142.6 |    995.2 |       5.4 |      |
|                               |     4 | 23,610,880 |       6,402,081 |  3884.8 |   515.0 |   4399.8 |      20.2 |      |
| relaxed / host / tri / sweep  |     0 |     92,230 |          17,982 |    53.9 |     7.4 |     61.3 |       7.5 | 0.86 |
|                               |     1 |    368,920 |          87,771 |    68.8 |    13.3 |     82.1 |       5.5 |      |
|                               |     2 |  1,475,680 |         379,818 |   195.8 |    69.0 |    264.8 |       2.1 |      |
|                               |     3 |  5,902,720 |       1,573,741 |   720.5 |   411.2 |   1131.7 |       6.5 |      |
|                               |     4 | 23,610,880 |       6,402,081 |  2671.7 |  4957.1 |   7628.8 |      22.4 |      |
| tight / host / tri / cell2d   |     0 |     92,230 |          17,982 |    27.2 |     8.3 |     35.5 |       5.3 | 0.87 |
|                               |     1 |    368,920 |          87,771 |    49.2 |    14.7 |     63.9 |       9.8 |      |
|                               |     2 |  1,475,680 |         379,818 |   173.8 |    42.7 |    216.5 |       5.1 |      |
|                               |     3 |  5,902,720 |       1,573,741 |   832.5 |   151.4 |    983.9 |       4.2 |      |
|                               |     4 | 23,610,880 |       6,402,081 |  3915.5 |   640.2 |   4555.7 |      15.0 |      |
| tight / host / tri / sweep    |     0 |     92,230 |          17,982 |    52.5 |     8.7 |     61.2 |       5.9 | 0.87 |
|                               |     1 |    368,920 |          87,771 |    72.3 |    13.5 |     85.9 |       8.8 |      |
|                               |     2 |  1,475,680 |         379,818 |   204.5 |    58.0 |    262.5 |       2.6 |      |
|                               |     3 |  5,902,720 |       1,573,741 |   739.2 |   431.6 |   1170.8 |       7.0 |      |
|                               |     4 | 23,610,880 |       6,402,081 |  2726.8 |  5024.1 |   7750.9 |      18.2 |      |

The two frames used here do not come into contact, so the narrow phase has almost no work to do and its column is dominated by noise rather than by element count; what this measures is the broad phase and the preparation that feeds it. Narrow-phase cost against problem size is in the per-case figure, over cases that do collide. Where the exponent is below 1 it is because the fixed cost visible at the smallest size is amortised as the mesh grows. `prep` builds the acceleration structure -- the cell list's grid or the sweep's sorted intervals -- and `step` is the traversal that reports pairs; the two strategies divide the work between those columns quite differently.

Source: `benchmark/results/scaling/host-cell2d-mode0.txt, benchmark/results/scaling/host-sweep-mode0.txt, benchmark/results/scaling/host-cell2d-mode2.txt, benchmark/results/scaling/host-sweep-mode2.txt`

<!-- sccd:end scaling -->

The two strategies divide their cost oppositely: the sweep builds more cheaply
and scales better doing it (exponent 0.73 against 0.92), the cell list traverses
far better (0.75 against 1.19). Traversal is what grows.

![Cost against element count](figures/refine-scaling.png)

**Figure 1.** Broad- and narrow-phase cost against element count, log-log. The
fitted exponent is on the broad phase; the narrow-phase series is flat and noisy
because these two frames do not come into contact.

## 9. Figures

Presentation follows the dataset paper (Belgrod et al., *TOI dataset for CCD and
a scalable conservative algorithm*) over the same six scenes: scenes across the
columns, one quantity per row, distributions over cases as box plots on
logarithmic axes. A **star** marks a parallel CPU mode, a **dagger** a GPU one.

![Per-case results over the six scenes](figures/results-grid.png)

**Figure 2.** Per-case distributions over the six scenes (columns): broad-phase
time, narrow-phase time, false positives, and earliness against the exact root.
Each box spans the first to the third quartile with the median inside, whiskers
reach the furthest case within 1.5 interquartile ranges, and cases beyond are
drawn individually. Counts are on a symmetric-log axis so a case with no false
positive — the common outcome — stays on the axis.

![Runtime split by phase](figures/runtime-breakdown.png)

**Figure 3.** Runtime split into *prep*, *broad* and *narrow*. Preparation
dominates on the host; on the GPU the narrow phase does.

![Time-of-impact error against the symbolic ground truth](figures/toi-error.png)

**Figure 4.** Error against the exact symbolic roots, log axis. The error is
one-sided by construction — how far *before* the true root — so every value is on
the safe side, and nothing falls off the axis on the late side because no such
case exists.

![Narrow-phase time per case against candidate pairs](figures/narrow-per-case.png)

**Figure 5.** Narrow-phase time per case against the candidate pairs the broad
phase handed it, log-log.

## 10. Provenance

<!-- sccd:begin provenance -->

- Timings: `benchmark/results/sweep-gh200-bp.csv`, 6394 cases over 6 scenes, 3 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200-all.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
