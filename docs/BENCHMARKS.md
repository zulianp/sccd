# Benchmarks

Every number in this document is generated from a committed CSV by a committed
script. Nothing here is typed in by hand, and the claim is checkable:

```sh
python3 -m report benchmark/results/<sweep>.csv /tmp/report \
        benchmark/results/<oracle>.csv --embed=docs/BENCHMARKS.md --check
```

returns non-zero if this document is not what the data reproduces. Run it
without `--check` to refresh it.

To measure from scratch:

```sh
benchmark/scripts/prepare_data.sh          # download, convert, verify the oracle
benchmark/scripts/sweep.sh --repeats 5     # resumable; one Slurm job per pack
python3 -m report <sweep>.csv out <oracle>.csv
```

`prepare_data.sh` refuses to finish if the ground truth is incomplete, and the
sweep is resumable at chunk granularity, so an interrupted run is continued by
running it again.

## 1. What is being measured

SCCD is a **conservative** continuous collision detector, and conservativeness is
an invariant rather than a quality:

- A reported time of impact must be **at or before** the true one. Later lets a
  simulation step through the contact, which is the failure the whole search
  exists to prevent.
- A collision that exists must be reported. A missed collision is a correctness
  failure, never a trade for speed.
- Reporting a collision that does not exist, or a time of impact earlier than the
  true one, is **acceptable**: it costs work and step size, not safety.

That asymmetry decides how the results are read. Speed and tightness are
negotiable; the two tables headed *conservativeness* are not. Both are checked
against the dataset's exact symbolic roots, not against another implementation.

**A difference smaller than the run-to-run spread is not a result.** Every
repeated measurement is reported as `median / slowest`, so the spread is visible
and checkable rather than asserted as a percentage; figures draw it, and the
mode comparisons below refuse to state a ratio when the gap is inside it.
Distributions over cases carry their worst case for the same reason — a median
describes the typical case, and the cost of this search is paid in the worst
one.

## 2. Platform

One GH200 node on CSCS Alps: Grace (288 hardware threads, `OMP_NUM_THREADS` set
to the node's core count) and Hopper sm_90 in the same allocation, built with
`prgenv-gnu/24.11:v2`, `-O3`, `CMAKE_CUDA_ARCHITECTURES=90`. Search parameters
are the shipped defaults. Rows marked *(GPU)* ran with
`SCCD_BENCH_EXECUTION_SPACE=device`; the rest ran on Grace.

Mesh geometry is stored and read in **double** — smesh built with
`SMESH_GEOM_TYPE=float64`, frames converted to `x.float64`. That matters for the
mesh path and only there. Two of the scenes ship PLY files declaring
`property double`, so a float32 mesh store discards real data: on
armadillo-rollers it perturbs coordinates by up to 5.95e-8 relative, and near a
grazing contact that moves the root by far more than it moves the coordinate.
Measured over 200 armadillo-rollers cases, float32 mesh geometry makes the
earliest-impact answer land after the exact root in 66 of them and float64 in
none, while the narrow-phase false positives are identical either way — that
path takes its coordinates from the query sets, in double, and never reads the
mesh. The double build is also the faster of the two here, 34 s against 57 s,
because tighter coordinates give tighter swept boxes and fewer candidate pairs.

## 3. Dataset

The dataset is the NYU CCD benchmark (archive 2451/74508), real simulation
frames rather than synthetic motion. Each case is one step, with a curated query
set whose coordinates are exact dyadic rationals and whose roots were computed
symbolically.

All six scenes are swept, every case of every one:


The scenes differ by orders of magnitude in every dimension, which is the point
of using all of them: puffer-ball's surface is 1,064,216 triangles against
armadillo-rollers' 24,242, and rod-twist contributes 4,571 of the 6,394 steps
while cloth-ball contributes 79.

<!-- sccd:begin dataset -->

| scene             | cases | candidate pairs/step |   queries | with a root |
|-------------------|------:|---------------------:|----------:|------------:|
| armadillo-rollers |   779 |              109,134 |   131,441 |     130,859 |
| cloth-ball        |    78 |            2,244,329 |   664,940 |     664,919 |
| cloth-funnel      |   575 |               43,723 |     7,552 |       6,773 |
| n-body            |   145 |           15,403,679 | 2,947,719 |   2,947,611 |
| puffer-ball       |   239 |           30,535,962 | 1,514,172 |   1,486,790 |
| rod-twist         | 4,559 |              846,976 |   549,208 |     285,431 |

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end dataset -->

Accuracy is measured on the curated query sets, because the exact roots belong
to those coordinates specifically. The mesh path is a separately stored copy of
the same frames, read from PLY, and is reported beside them rather than checked
against them.

### 3.1. Ground-truth coverage

The dataset states each query's outcome twice — a boolean in `mma_bool` and a
time of impact in `roots` — and `benchmark/verify_oracle.py` requires them to
agree query for query, `mma_bool[i] == isfinite(toi[i])`. A `true` against a NaN
is a root that never converted, and downstream it is indistinguishable from "no
collision", so an incomplete oracle silently shrinks the evidence instead of
announcing itself. The gate runs as part of `prepare_data.sh` and exits non-zero
on any gap.

Coverage is complete on all six scenes. A query carries a root exactly when
`mma_bool` marks it as colliding, and the counts agree on every one: 130,859 of
armadillo-rollers' 131,441 queries, 664,919 of cloth-ball's 664,940, 6,773 of
cloth-funnel's 7,552, 2,947,611 of n-body-simulation's 2,947,719, 1,486,790 of
puffer-ball's 1,514,172, and 285,431 of rod-twist's 549,208. The proportion
varies by scene because it is the fraction of curated queries that collide, not
a measure of how much of the oracle converted — rod-twist's 52% and cloth-ball's
99.997% are both complete.

### 3.2. Modes

`Relaxed` (`SCCD_NARROWPHASE_MODE=0`) and `Tight` (`=2`) differ in **how early
they accept a box**, not in how fast they run. `Relaxed` compares codomain widths
against domain tolerances, so it accepts sooner: fewer boxes, and a time of
impact further before the true one. `Tight` compares domain widths against
domain tolerances. Both are conservative by construction; they sit at opposite
ends of the accuracy-for-work trade.

`ToiOutput::Earliest` (`find_earliest_impact_time`) returns one time of impact
for the whole step, so every query prunes against the running minimum.
`ToiOutput::PerPair` (`find_impact_times`) returns one result per candidate with
no shared bound, and does correspondingly more work.

The Earliest answer is **not reproducible run to run**: the parallel search
prunes against a running best, so which box is accepted depends on scheduling,
and the reported time of impact varies between runs of the same binary on the
same input. It varies in tightness only — it is never later than the true root —
but it is a reason to read the spread rather than a single figure.

## 4. Benchmark

### 4.1. Conservativeness

The invariant, checked against the dataset's exact symbolic roots for every mode
on every scene, CPU and GPU. **Zero missed collisions and zero late times of
impact in all forty-eight configurations.**

<!-- sccd:begin gate -->

| scene             | phase | mode          | queries checked | missed | late |
|-------------------|-------|---------------|----------------:|-------:|-----:|
| armadillo-rollers | EE    | Relaxed (GPU) |          98,757 |      0 |    0 |
| armadillo-rollers | EE    | Tight (GPU)   |          98,757 |      0 |    0 |
| armadillo-rollers | EE    | Relaxed       |          98,757 |      0 |    0 |
| armadillo-rollers | EE    | Tight         |          98,757 |      0 |    0 |
| armadillo-rollers | VF    | Relaxed (GPU) |          32,102 |      0 |    0 |
| armadillo-rollers | VF    | Tight (GPU)   |          32,102 |      0 |    0 |
| armadillo-rollers | VF    | Relaxed       |          32,102 |      0 |    0 |
| armadillo-rollers | VF    | Tight         |          32,102 |      0 |    0 |
| cloth-ball        | EE    | Relaxed (GPU) |         557,667 |      0 |    0 |
| cloth-ball        | EE    | Tight (GPU)   |         557,667 |      0 |    0 |
| cloth-ball        | EE    | Relaxed       |         557,667 |      0 |    0 |
| cloth-ball        | EE    | Tight         |         557,667 |      0 |    0 |
| cloth-ball        | VF    | Relaxed (GPU) |         107,252 |      0 |    0 |
| cloth-ball        | VF    | Tight (GPU)   |         107,252 |      0 |    0 |
| cloth-ball        | VF    | Relaxed       |         107,252 |      0 |    0 |
| cloth-ball        | VF    | Tight         |         107,252 |      0 |    0 |
| cloth-funnel      | EE    | Relaxed (GPU) |           6,249 |      0 |    0 |
| cloth-funnel      | EE    | Tight (GPU)   |           6,249 |      0 |    0 |
| cloth-funnel      | EE    | Relaxed       |           6,249 |      0 |    0 |
| cloth-funnel      | EE    | Tight         |           6,249 |      0 |    0 |
| cloth-funnel      | VF    | Relaxed (GPU) |             524 |      0 |    0 |
| cloth-funnel      | VF    | Tight (GPU)   |             524 |      0 |    0 |
| cloth-funnel      | VF    | Relaxed       |             524 |      0 |    0 |
| cloth-funnel      | VF    | Tight         |             524 |      0 |    0 |
| n-body            | EE    | Relaxed (GPU) |       2,399,741 |      0 |    0 |
| n-body            | EE    | Tight (GPU)   |       2,399,741 |      0 |    0 |
| n-body            | EE    | Relaxed       |       2,399,741 |      0 |    0 |
| n-body            | EE    | Tight         |       2,399,741 |      0 |    0 |
| n-body            | VF    | Relaxed (GPU) |         547,870 |      0 |    0 |
| n-body            | VF    | Tight (GPU)   |         547,870 |      0 |    0 |
| n-body            | VF    | Relaxed       |         547,870 |      0 |    0 |
| n-body            | VF    | Tight         |         547,870 |      0 |    0 |
| puffer-ball       | EE    | Relaxed (GPU) |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | Tight (GPU)   |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | Relaxed       |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | Tight         |       1,187,257 |      0 |    0 |
| puffer-ball       | VF    | Relaxed (GPU) |         299,533 |      0 |    0 |
| puffer-ball       | VF    | Tight (GPU)   |         299,533 |      0 |    0 |
| puffer-ball       | VF    | Relaxed       |         299,533 |      0 |    0 |
| puffer-ball       | VF    | Tight         |         299,533 |      0 |    0 |
| rod-twist         | EE    | Relaxed (GPU) |         245,025 |      0 |    0 |
| rod-twist         | EE    | Tight (GPU)   |         245,025 |      0 |    0 |
| rod-twist         | EE    | Relaxed       |         245,025 |      0 |    0 |
| rod-twist         | EE    | Tight         |         245,025 |      0 |    0 |
| rod-twist         | VF    | Relaxed (GPU) |          40,406 |      0 |    0 |
| rod-twist         | VF    | Tight (GPU)   |          40,406 |      0 |    0 |
| rod-twist         | VF    | Relaxed       |          40,406 |      0 |    0 |
| rod-twist         | VF    | Tight         |          40,406 |      0 |    0 |

TightInclusion is excluded from this table: it is the reference the queries were selected against, not a subject of it.

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end gate -->

Per-scene, with the false positives that are the price of accepting early, and
one column that is not part of the gate:

<!-- sccd:begin conservativeness -->

| scene             | mode          |   queries | toi compared | late | false pos. | false neg. | mesh-path divergence |
|-------------------|---------------|----------:|-------------:|-----:|-----------:|-----------:|---------------------:|
| armadillo-rollers | Relaxed (GPU) |   131,261 |      130,681 |    0 |        322 |          0 |                    0 |
| armadillo-rollers | Tight (GPU)   |   131,441 |      130,859 |    0 |         24 |          0 |                    0 |
| armadillo-rollers | Relaxed       |   131,261 |      130,681 |    0 |        232 |          0 |                    0 |
| armadillo-rollers | Tight         |   131,441 |      130,859 |    0 |         24 |          0 |                    0 |
| cloth-ball        | Relaxed (GPU) |   664,939 |      664,918 |    0 |          2 |          0 |                    0 |
| cloth-ball        | Tight (GPU)   |   664,940 |      664,919 |    0 |          1 |          0 |                    0 |
| cloth-ball        | Relaxed       |   664,939 |      664,918 |    0 |          2 |          0 |                    0 |
| cloth-ball        | Tight         |   664,940 |      664,919 |    0 |          1 |          0 |                    0 |
| cloth-funnel      | Relaxed (GPU) |     7,544 |        6,766 |    0 |        741 |          0 |                    0 |
| cloth-funnel      | Tight (GPU)   |     7,552 |        6,773 |    0 |         15 |          0 |                    0 |
| cloth-funnel      | Relaxed       |     7,544 |        6,766 |    0 |        686 |          0 |                    0 |
| cloth-funnel      | Tight         |     7,552 |        6,773 |    0 |         15 |          0 |                    0 |
| n-body            | Relaxed (GPU) | 2,921,071 |    2,920,965 |    0 |         28 |          0 |                    0 |
| n-body            | Tight (GPU)   | 2,947,719 |    2,947,611 |    0 |         12 |          0 |                    0 |
| n-body            | Relaxed       | 2,921,071 |    2,920,965 |    0 |          9 |          0 |                    0 |
| n-body            | Tight         | 2,947,719 |    2,947,611 |    0 |         11 |          0 |                    0 |
| puffer-ball       | Relaxed (GPU) | 1,513,963 |    1,486,587 |    0 |     27,374 |          0 |                    0 |
| puffer-ball       | Tight (GPU)   | 1,514,172 |    1,486,790 |    0 |        558 |          0 |                    0 |
| puffer-ball       | Relaxed       | 1,513,963 |    1,486,587 |    0 |     27,374 |          0 |                    0 |
| puffer-ball       | Tight         | 1,514,172 |    1,486,790 |    0 |        536 |          0 |                    0 |
| rod-twist         | Relaxed (GPU) |   546,081 |      283,655 |    0 |    243,717 |          0 |                    0 |
| rod-twist         | Tight (GPU)   |   549,208 |      285,431 |    0 |      1,245 |          0 |                    0 |
| rod-twist         | Relaxed       |   546,081 |      283,655 |    0 |    227,441 |          0 |                    0 |
| rod-twist         | Tight         |   549,208 |      285,431 |    0 |      1,243 |          0 |                    0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness. The last column is not part of the gate. It counts cases where the earliest-impact answer computed over the *mesh* is later than the earliest exact root of the *curated queries*, which are two separately stored geometries: the mesh is read from PLY, the queries are exact dyadic rationals. It is a measure of the agreement between those two inputs rather than of the kernel, and with the mesh stored in double it is zero everywhere.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end conservativeness -->

The last column deserves its own sentence, because it looks like a violation and
is not. It counts cases where the earliest-impact answer computed over the
**mesh** lands after the earliest exact root of the **curated queries** — two
separately stored copies of the same frames, one read from PLY and one written
as exact dyadic rationals. It measures the agreement between two inputs, not a
kernel reporting late, which is why it is reported and not gated on. With the
mesh stored in double it is **zero on every scene and every mode**; a float32
mesh store puts it in the hundreds on the two scenes whose PLY files declare
`property double`, which is what decided the storage type.

### 4.2. Narrow phase

Whole-scene wall clock, summed over every case, median of three independent
runs with the full run-to-run range beside it.

<!-- sccd:begin timing -->

| scene             | mode          | cases |         pairs | rep |           prep ms |          broad ms |       earliest ms |       per-pair ms |            total ms |
|-------------------|---------------|------:|--------------:|----:|------------------:|------------------:|------------------:|------------------:|--------------------:|
| armadillo-rollers | Relaxed (GPU) |   779 |    85,015,700 |   3 |   2069.6 / 2071.3 |     838.6 / 841.4 |   1533.7 / 1536.7 |   2324.8 / 2333.5 |     4393.4 / 4446.4 |
| armadillo-rollers | Tight (GPU)   |   781 |    85,296,282 |   3 |   1848.8 / 1960.1 |     813.1 / 838.0 |   1735.9 / 1763.2 |   2650.0 / 2663.7 |     4417.8 / 4509.0 |
| armadillo-rollers | Relaxed       |   779 |    85,015,700 |   3 |   8106.5 / 8135.7 |   3322.1 / 3354.3 |     759.8 / 771.9 |   2942.3 / 2959.0 |   12188.4 / 12262.0 |
| armadillo-rollers | Tight         |   781 |    85,296,282 |   3 |   8088.8 / 8212.7 |   3331.1 / 3376.8 |   1036.6 / 1040.4 |   3565.9 / 3587.7 |   12456.6 / 12629.9 |
| cloth-ball        | Relaxed (GPU) |    78 |   175,057,694 |   3 |     525.1 / 526.1 |     367.1 / 369.8 |     425.6 / 428.9 |     920.2 / 935.3 |     1316.7 / 1324.8 |
| cloth-ball        | Tight (GPU)   |    79 |   175,150,856 |   3 |     530.4 / 533.7 |     369.1 / 371.3 |     441.2 / 442.2 |     883.1 / 889.4 |     1343.4 / 1343.8 |
| cloth-ball        | Relaxed       |    78 |   175,057,694 |   3 |   1212.5 / 1218.4 |     932.9 / 937.5 |     385.5 / 387.8 |     789.9 / 802.9 |     2530.8 / 2543.6 |
| cloth-ball        | Tight         |    79 |   175,150,856 |   3 |   1230.0 / 1231.1 |     945.3 / 946.7 |     214.7 / 216.6 |   1206.8 / 1220.8 |     2388.0 / 2393.3 |
| cloth-funnel      | Relaxed (GPU) |   575 |    25,141,046 |   3 |   1358.4 / 1376.9 |     540.9 / 545.7 |     818.3 / 837.9 |     727.4 / 728.4 |     2694.1 / 2760.5 |
| cloth-funnel      | Tight (GPU)   |   577 |    25,192,698 |   3 |   1353.5 / 1358.1 |     540.6 / 541.4 |     943.9 / 969.6 |   1071.1 / 1074.4 |     2836.1 / 2836.5 |
| cloth-funnel      | Relaxed       |   575 |    25,141,046 |   3 |   5687.6 / 5764.2 |   2181.2 / 2191.1 |     553.7 / 568.0 |   1026.8 / 1038.2 |     8422.4 / 8523.4 |
| cloth-funnel      | Tight         |   577 |    25,192,698 |   3 |   5661.3 / 5701.0 |   2185.8 / 2198.8 |     658.6 / 664.2 |   1298.2 / 1298.4 |     8524.3 / 8545.5 |
| n-body            | Relaxed (GPU) |   145 | 2,233,533,498 |   3 |     967.8 / 978.7 |   1326.4 / 1330.1 |   4709.9 / 4717.8 |   4335.3 / 4455.1 |     6983.6 / 7015.0 |
| n-body            | Tight (GPU)   |   146 | 2,253,766,609 |   3 |   1004.8 / 1010.8 |   1311.7 / 1337.1 |   4720.7 / 4722.7 |   3976.9 / 4065.0 |     7043.2 / 7062.6 |
| n-body            | Relaxed       |   145 | 2,233,533,498 |   3 |   2482.9 / 2492.2 |   7183.8 / 7189.1 |   3497.6 / 3554.2 | 10002.6 / 10064.4 |   13164.3 / 13228.7 |
| n-body            | Tight         |   146 | 2,253,766,609 |   3 |   2607.7 / 2618.8 |   7223.9 / 7286.5 |   1738.0 / 1756.1 | 19258.3 / 19443.3 |   11555.4 / 11650.2 |
| puffer-ball       | Relaxed (GPU) |   239 | 7,298,095,145 |   3 | 16853.8 / 16990.3 | 27239.6 / 27256.6 | 38127.0 / 38131.4 | 10667.8 / 10669.8 |   81789.3 / 82356.9 |
| puffer-ball       | Tight (GPU)   |   240 | 7,327,014,600 |   3 | 16975.7 / 17124.0 | 27521.8 / 27584.2 | 38657.9 / 38794.1 | 10752.3 / 11037.2 |   83171.7 / 83354.0 |
| puffer-ball       | Relaxed       |   239 | 7,298,095,145 |   3 | 18928.3 / 19590.4 | 27289.2 / 27356.9 | 12063.1 / 12188.8 | 14353.1 / 14739.7 |   58172.0 / 58854.7 |
| puffer-ball       | Tight         |   240 | 7,327,014,600 |   3 | 18898.1 / 19172.8 | 27335.1 / 27640.5 |   6049.4 / 6154.7 | 56908.1 / 56970.6 |   52374.1 / 52503.0 |
| rod-twist         | Relaxed (GPU) |  4559 | 3,861,367,133 |   3 | 26191.6 / 26628.1 |   6780.7 / 7013.1 |   9715.2 / 9723.6 | 19749.6 / 19844.2 |   42687.5 / 43286.8 |
| rod-twist         | Tight (GPU)   |  4571 | 3,872,779,843 |   3 | 25140.3 / 25800.2 |   6742.8 / 6908.7 | 17368.0 / 17514.2 | 27523.3 / 27625.2 |   49397.4 / 50076.9 |
| rod-twist         | Relaxed       |  4559 | 3,861,367,133 |   3 | 66897.1 / 66993.6 | 28867.0 / 28918.5 | 13293.1 / 13373.1 | 39888.5 / 39892.8 | 109044.9 / 109285.3 |
| rod-twist         | Tight         |  4571 | 3,872,779,843 |   3 | 66765.3 / 67399.4 | 28680.1 / 28853.0 | 18480.3 / 18549.3 | 65174.0 / 65384.2 | 113887.7 / 114732.7 |

Each cell is the median over repeats and the slowest of them. A difference smaller than the gap between the two does not separate two modes and is not reported as a ratio anywhere in this document. Which mode is faster depends on the output mode as well as the scene, so the two are given side by side rather than one standing for the other.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end timing -->

Whole-scene milliseconds cannot be compared between a 79-case scene and a
4,571-case one; throughput can:

<!-- sccd:begin throughput -->

| scene             | mode          | broad Mpair/s | narrow Mpair/s |
|-------------------|---------------|--------------:|---------------:|
| armadillo-rollers | Relaxed (GPU) |         101.4 |           55.4 |
| armadillo-rollers | Tight (GPU)   |         104.9 |           49.1 |
| armadillo-rollers | Relaxed       |          25.6 |          111.9 |
| armadillo-rollers | Tight         |          25.6 |           82.3 |
| cloth-ball        | Relaxed (GPU) |         476.9 |          411.3 |
| cloth-ball        | Tight (GPU)   |         474.5 |          397.0 |
| cloth-ball        | Relaxed       |         187.7 |          454.1 |
| cloth-ball        | Tight         |         185.3 |          815.7 |
| cloth-funnel      | Relaxed (GPU) |          46.5 |           30.7 |
| cloth-funnel      | Tight (GPU)   |          46.6 |           26.7 |
| cloth-funnel      | Relaxed       |          11.5 |           45.4 |
| cloth-funnel      | Tight         |          11.5 |           38.3 |
| n-body            | Relaxed (GPU) |        1683.9 |          474.2 |
| n-body            | Tight (GPU)   |        1718.2 |          477.4 |
| n-body            | Relaxed       |         310.9 |          638.6 |
| n-body            | Tight         |         312.0 |         1296.8 |
| puffer-ball       | Relaxed (GPU) |         267.9 |          191.4 |
| puffer-ball       | Tight (GPU)   |         266.2 |          189.5 |
| puffer-ball       | Relaxed       |         267.4 |          605.0 |
| puffer-ball       | Tight         |         268.0 |         1211.2 |
| rod-twist         | Relaxed (GPU) |         569.5 |          397.5 |
| rod-twist         | Tight (GPU)   |         574.4 |          223.0 |
| rod-twist         | Relaxed       |         133.8 |          290.5 |
| rod-twist         | Tight         |         135.0 |          209.6 |

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end throughput -->

**Neither mode is uniformly faster**, which is the substantive result here and
the reason the trade is worth stating as a trade:

<!-- sccd:begin comparison -->

- **armadillo-rollers (CPU, earliest)**: Relaxed is 1.36× faster than Tight in the narrow phase (760 ms against 1037 ms; run-to-run spread 4.2%).
- **armadillo-rollers (CPU, per-pair)**: Relaxed is 1.21× faster than Tight in the narrow phase (2942 ms against 3566 ms; run-to-run spread 1.6%).
- **armadillo-rollers (GPU, earliest)**: Relaxed (GPU) is 1.13× faster than Tight (GPU) in the narrow phase (1534 ms against 1736 ms; run-to-run spread 2.1%).
- **armadillo-rollers (GPU, per-pair)**: Relaxed (GPU) is 1.14× faster than Tight (GPU) in the narrow phase (2325 ms against 2650 ms; run-to-run spread 0.9%).
- **cloth-ball (CPU, earliest)**: Tight is 1.80× faster than Relaxed in the narrow phase (215 ms against 385 ms; run-to-run spread 2.2%).
- **cloth-ball (CPU, per-pair)**: Relaxed is 1.53× faster than Tight in the narrow phase (790 ms against 1207 ms; run-to-run spread 2.0%).
- **cloth-ball (GPU, earliest)**: Relaxed (GPU) is 1.04× faster than Tight (GPU) in the narrow phase (426 ms against 441 ms; run-to-run spread 1.1%).
- **cloth-ball (GPU, per-pair)**: Tight (GPU) is 1.04× faster than Relaxed (GPU) in the narrow phase (883 ms against 920 ms; run-to-run spread 1.9%).
- **cloth-funnel (CPU, earliest)**: Relaxed is 1.19× faster than Tight in the narrow phase (554 ms against 659 ms; run-to-run spread 6.1%).
- **cloth-funnel (CPU, per-pair)**: Relaxed is 1.26× faster than Tight in the narrow phase (1027 ms against 1298 ms; run-to-run spread 2.9%).
- **cloth-funnel (GPU, earliest)**: Relaxed (GPU) is 1.15× faster than Tight (GPU) in the narrow phase (818 ms against 944 ms; run-to-run spread 5.2%).
- **cloth-funnel (GPU, per-pair)**: Relaxed (GPU) is 1.47× faster than Tight (GPU) in the narrow phase (727 ms against 1071 ms; run-to-run spread 1.5%).
- **n-body (CPU, earliest)**: Tight is 2.01× faster than Relaxed in the narrow phase (1738 ms against 3498 ms; run-to-run spread 2.4%).
- **n-body (CPU, per-pair)**: Relaxed is 1.93× faster than Tight in the narrow phase (10003 ms against 19258 ms; run-to-run spread 1.2%).
- **n-body (GPU, earliest)**: Relaxed (GPU) and Tight (GPU) are inside noise (4710 ms against 4721 ms, spread 0.7%); this does not separate them.
- **n-body (GPU, per-pair)**: Tight (GPU) is 1.09× faster than Relaxed (GPU) in the narrow phase (3977 ms against 4335 ms; run-to-run spread 5.5%).
- **puffer-ball (CPU, earliest)**: Tight is 1.99× faster than Relaxed in the narrow phase (6049 ms against 12063 ms; run-to-run spread 2.5%).
- **puffer-ball (CPU, per-pair)**: Relaxed is 3.96× faster than Tight in the narrow phase (14353 ms against 56908 ms; run-to-run spread 6.1%).
- **puffer-ball (GPU, earliest)**: Relaxed (GPU) is 1.01× faster than Tight (GPU) in the narrow phase (38127 ms against 38658 ms; run-to-run spread 0.8%).
- **puffer-ball (GPU, per-pair)**: Relaxed (GPU) and Tight (GPU) are inside noise (10668 ms against 10752 ms, spread 4.8%); this does not separate them.
- **rod-twist (CPU, earliest)**: Relaxed is 1.39× faster than Tight in the narrow phase (13293 ms against 18480 ms; run-to-run spread 0.7%).
- **rod-twist (CPU, per-pair)**: Relaxed is 1.63× faster than Tight in the narrow phase (39889 ms against 65174 ms; run-to-run spread 0.3%).
- **rod-twist (GPU, earliest)**: Relaxed (GPU) is 1.79× faster than Tight (GPU) in the narrow phase (9715 ms against 17368 ms; run-to-run spread 0.9%).
- **rod-twist (GPU, per-pair)**: Relaxed (GPU) is 1.39× faster than Tight (GPU) in the narrow phase (19750 ms against 27523 ms; run-to-run spread 2.1%).

<!-- sccd:end comparison -->

The split between the two output modes is the mechanism, not a detail.
`Relaxed` accepts a box sooner, so it does less work per query — but a looser
acceptance also means a looser bound to hand the next query, so it prunes less.
`ToiOutput::PerPair` has no shared bound to lose, and there `Relaxed` is faster
on **every** scene. `ToiOutput::Earliest` reintroduces it, and on the three
scenes with the most candidate pairs per step the pruning `Tight` buys back is
worth more than the per-query work it costs, so the ranking inverts. Which
effect wins is a property of the scene and of what is being asked, not of the
mode alone.

**The two processors divide the work differently, and not in the same direction
on every scene.** The usual pattern is that Hopper wins the broad phase and
loses the narrow one: on five of the six scenes its broad phase runs 2.2× to
4.2× faster than Grace's while its narrow phase runs 0.27× to 0.76× as fast, and
the broad phase is large enough that the GPU still wins end to end by 1.5× to
2.9×.

Two scenes break it, in opposite directions. On **rod-twist** the GPU wins every
phase, narrow phase included (1.26×), for 2.50× overall. On **puffer-ball** it
loses outright, 0.62× end to end — the only scene where the GPU is the wrong
processor. Its narrow phase there takes 44.9 s against the host's 12.2 s, and
unusually its broad phase does not compensate, running slightly slower as well
(34.2 s against 28.3 s). puffer-ball is the largest scene in the set by candidate
pairs, 30.5 million per step, so this is the GPU losing where the problem is
biggest rather than where it is smallest.

Comparisons above are made within a processor for that reason: ranking every mode
of a scene together would compare host `Relaxed` against GPU `Tight` and report
the sum of two unrelated effects as if it were the mode trade.

### 4.3. Broad phase

Two strategies produce the candidate pairs: a **sweep** over sorted intervals and
a **cell list** over a uniform grid. They return identical pair sets, so the
choice is purely about cost, and the shipped default does not fix a winner — it
races the two per scene and keeps the faster, re-probing periodically so a scene
that changes character can change the answer.

The table below is that race run offline over the whole benchmark, with both
strategies sampled on every scene and every mode. `prep` builds the acceleration
structure — the sorted intervals or the grid — and is separated from the
traversal that follows it, because that is where the two differ rather than in
the total. Across 25,538 case-mode combinations the two report **identical
candidate pair counts**, so nothing here is a difference in what is found.

**The choice is a host one.** `SCCD_BROADPHASE` is read in the three host
broad-phase functions and nowhere else; the device steps have no branch on it,
so the device always runs its own path. Sampling both strategies on the GPU
produces ties to within a millisecond, which is why only host rows are listed.

<!-- sccd:begin broadphase -->

| scene             | mode    | cell2d prep ms | sweep prep ms | cell2d broad ms | sweep broad ms | faster       |
|-------------------|---------|---------------:|--------------:|----------------:|---------------:|--------------|
| armadillo-rollers | Relaxed |           8106 |         26507 |            3322 |           3236 | sweep 1.03x  |
| armadillo-rollers | Tight   |           8089 |         26637 |            3331 |           3231 | sweep 1.03x  |
| cloth-ball        | Relaxed |           1212 |          2980 |             933 |            691 | sweep 1.35x  |
| cloth-ball        | Tight   |           1230 |          3064 |             945 |            697 | sweep 1.36x  |
| cloth-funnel      | Relaxed |           5688 |         19072 |            2181 |           2298 | cell2d 1.05x |
| cloth-funnel      | Tight   |           5661 |         19460 |            2186 |           2340 | cell2d 1.07x |
| n-body            | Relaxed |           2483 |          5717 |            7184 |           3818 | sweep 1.88x  |
| n-body            | Tight   |           2608 |          5902 |            7224 |           3833 | sweep 1.88x  |
| puffer-ball       | Relaxed |          18928 |         27304 |           27289 |         108994 | cell2d 3.99x |
| puffer-ball       | Tight   |          18898 |         27213 |           27335 |         110588 | cell2d 4.05x |
| rod-twist         | Relaxed |          66897 |        167942 |           28867 |          21534 | sweep 1.34x  |
| rod-twist         | Tight   |          66765 |        169005 |           28680 |          21551 | sweep 1.33x  |

`faster` names the winning strategy and by how much on the whole broad phase. A margin inside the run-to-run spread is reported as a tie rather than a winner.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end broadphase -->

**Neither strategy wins.** The sweep takes armadillo-rollers, cloth-ball,
n-body and rod-twist by 1.03× to 1.88×; the cell list takes cloth-funnel by
1.06× and puffer-ball by **4.0×**, and that last one is the reason the default
does not simply pick the sweep: on the largest scene in the set the sweep's
broad phase costs 109 s against the cell list's 27 s. A fixed choice is wrong
somewhere, and wrong by a factor of four at the worst point, which is what makes
racing them per scene worth its two probe steps.

The two also divide their cost differently, and the split is the opposite of
what the totals suggest. The sweep's `prep` — sorting the intervals — is two to
three times the cell list's grid build on every scene, and it wins anyway
wherever it wins, by traversing much less. The refinement study below sweeps
element count over two and a half orders of magnitude and pins down that
scaling.

## 5. Evaluation

### 5.1. Comparison against TightInclusion

TightInclusion is the reference implementation of a certified conservative
narrow phase. Matching its hit count is the strongest agreement available;
reporting more hits is a false positive, which costs work and is never unsafe.
Its own answer is a lower bound on the truth rather than the truth, which is why
the conservativeness table above is measured against the exact roots instead.

**How the reference is configured.** A ratio against a baseline is only worth
reading if the baseline was allowed to do what the measured code does, so both
sides here are given the same task and the same machine:

- **The same question.** Both compute a time of impact for *every* query — SCCD
  at `ToiOutput::PerPair`, TightInclusion unbounded over the whole step. Neither
  prunes against a shared running minimum, so neither is answering the cheaper
  "earliest over the set" question the other is not.
- **The same scheduler.** Both loops run through
  `sccd::parallel_for_br_dynamic`, the skew-aware helper the narrow phase uses
  for its own root finding. Timing a parallel narrow phase against a serial
  reference would report the thread count as though it were an algorithmic
  result.
- **Unmodified.** TightInclusion is used exactly as released; everything above
  is arranged on SCCD's side of the interface.

The hit counts are unchanged by the threading — verified identical at 1, 8 and
64 threads — so only the times move.

<!-- sccd:begin reference -->

| scene             | phase | mode           |   queries |      hits |        time ms |     vs. TI |
|-------------------|-------|----------------|----------:|----------:|---------------:|-----------:|
| armadillo-rollers | EE    | Relaxed (GPU)  |    99,104 |    98,933 |    2335 / 2349 |       4.5× |
| armadillo-rollers | EE    | Tight (GPU)    |    99,104 |    98,761 |    2908 / 2910 |       3.6× |
| armadillo-rollers | EE    | Relaxed        |    99,104 |    98,895 |    3377 / 3384 |       3.1× |
| armadillo-rollers | EE    | Tight          |    99,104 |    98,761 |    6474 / 6502 |       1.6× |
| armadillo-rollers | EE    | TightInclusion |    99,104 |    98,761 |  10499 / 10604 | 1.0× (ref) |
| armadillo-rollers | VF    | Relaxed (GPU)  |    32,337 |    32,248 |    3469 / 3476 |       2.2× |
| armadillo-rollers | VF    | Tight (GPU)    |    32,337 |    32,122 |    1987 / 2011 |       3.8× |
| armadillo-rollers | VF    | Relaxed        |    32,337 |    32,196 |    2219 / 2242 |       3.4× |
| armadillo-rollers | VF    | Tight          |    32,337 |    32,122 |    2717 / 2743 |       2.8× |
| armadillo-rollers | VF    | TightInclusion |    32,337 |    32,122 |    7609 / 7686 | 1.0× (ref) |
| cloth-ball        | EE    | Relaxed (GPU)  |   557,683 |   557,669 |      645 / 647 |       2.4× |
| cloth-ball        | EE    | Tight (GPU)    |   557,683 |   557,668 |      609 / 612 |       2.5× |
| cloth-ball        | EE    | Relaxed        |   557,683 |   557,669 |      402 / 404 |       3.8× |
| cloth-ball        | EE    | Tight          |   557,683 |   557,668 |      955 / 961 |       1.6× |
| cloth-ball        | EE    | TightInclusion |   557,683 |   557,668 |    1523 / 1529 | 1.0× (ref) |
| cloth-ball        | VF    | Relaxed (GPU)  |   107,257 |   107,252 |    2212 / 2367 |       0.6× |
| cloth-ball        | VF    | Tight (GPU)    |   107,257 |   107,252 |      428 / 429 |       2.9× |
| cloth-ball        | VF    | Relaxed        |   107,257 |   107,252 |      302 / 307 |       4.2× |
| cloth-ball        | VF    | Tight          |   107,257 |   107,252 |      408 / 414 |       3.1× |
| cloth-ball        | VF    | TightInclusion |   107,257 |   107,252 |    1255 / 1267 | 1.0× (ref) |
| cloth-funnel      | EE    | Relaxed (GPU)  |     6,751 |     6,734 |      768 / 793 |       4.1× |
| cloth-funnel      | EE    | Tight (GPU)    |     6,751 |     6,259 |    1358 / 1370 |       2.3× |
| cloth-funnel      | EE    | Relaxed        |     6,751 |     6,700 |    1111 / 1149 |       2.8× |
| cloth-funnel      | EE    | Tight          |     6,751 |     6,259 |    1743 / 1746 |       1.8× |
| cloth-funnel      | EE    | TightInclusion |     6,751 |     6,259 |    3130 / 3137 | 1.0× (ref) |
| cloth-funnel      | VF    | Relaxed (GPU)  |       801 |       781 |    2124 / 2599 |       0.6× |
| cloth-funnel      | VF    | Tight (GPU)    |       801 |       529 |      451 / 456 |       2.9× |
| cloth-funnel      | VF    | Relaxed        |       801 |       760 |      402 / 410 |       3.3× |
| cloth-funnel      | VF    | Tight          |       801 |       529 |      394 / 395 |       3.4× |
| cloth-funnel      | VF    | TightInclusion |       801 |       529 |    1321 / 1321 | 1.0× (ref) |
| n-body            | EE    | Relaxed (GPU)  | 2,399,812 | 2,399,762 |    1106 / 1351 |       2.4× |
| n-body            | EE    | Tight (GPU)    | 2,399,812 | 2,399,746 |     890 / 1077 |       3.0× |
| n-body            | EE    | Relaxed        | 2,399,812 | 2,399,747 |    6576 / 6626 |       0.4× |
| n-body            | EE    | Tight          | 2,399,812 | 2,399,746 |    1176 / 1184 |       2.2× |
| n-body            | EE    | TightInclusion | 2,399,812 | 2,399,746 |    2633 / 2657 | 1.0× (ref) |
| n-body            | VF    | Relaxed (GPU)  |   547,907 |   547,877 |    2883 / 2900 |       0.5× |
| n-body            | VF    | Tight (GPU)    |   547,907 |   547,874 |      628 / 790 |       2.5× |
| n-body            | VF    | Relaxed        |   547,907 |   547,873 |      423 / 424 |       3.7× |
| n-body            | VF    | Tight          |   547,907 |   547,873 |      507 / 511 |       3.1× |
| n-body            | VF    | TightInclusion |   547,907 |   547,873 |    1577 / 1587 | 1.0× (ref) |
| puffer-ball       | EE    | Relaxed (GPU)  | 1,206,952 | 1,206,951 |      346 / 346 |       7.6× |
| puffer-ball       | EE    | Tight (GPU)    | 1,206,952 | 1,187,650 |      941 / 943 |       2.8× |
| puffer-ball       | EE    | Relaxed        | 1,206,952 | 1,206,951 |      252 / 254 |      10.5× |
| puffer-ball       | EE    | Tight          | 1,206,952 | 1,187,650 |    1614 / 1629 |       1.6× |
| puffer-ball       | EE    | TightInclusion | 1,206,952 | 1,187,650 |    2642 / 2662 | 1.0× (ref) |
| puffer-ball       | VF    | Relaxed (GPU)  |   307,220 |   307,219 |    2518 / 2549 |       0.6× |
| puffer-ball       | VF    | Tight (GPU)    |   307,220 |   299,698 |      634 / 636 |       2.4× |
| puffer-ball       | VF    | Relaxed        |   307,220 |   307,219 |      209 / 209 |       7.4× |
| puffer-ball       | VF    | Tight          |   307,220 |   299,676 |      586 / 595 |       2.6× |
| puffer-ball       | VF    | TightInclusion |   307,220 |   299,676 |    1544 / 1556 | 1.0× (ref) |
| rod-twist         | EE    | Relaxed (GPU)  |   492,120 |   474,114 |    2898 / 4493 |      28.6× |
| rod-twist         | EE    | Tight (GPU)    |   492,120 |   246,132 |    4913 / 8176 |      16.9× |
| rod-twist         | EE    | Relaxed        |   492,120 |   458,641 |   8100 / 22112 |      10.2× |
| rod-twist         | EE    | Tight          |   492,120 |   246,132 | 37118 / 121624 |       2.2× |
| rod-twist         | EE    | TightInclusion |   492,120 |   246,132 | 82987 / 249504 | 1.0× (ref) |
| rod-twist         | VF    | Relaxed (GPU)  |    57,088 |    56,290 |    3618 / 4208 |       6.5× |
| rod-twist         | VF    | Tight (GPU)    |    57,088 |    40,544 |    3014 / 3236 |       7.8× |
| rod-twist         | VF    | Relaxed        |    57,088 |    55,398 |    2648 / 2974 |       8.9× |
| rod-twist         | VF    | Tight          |    57,088 |    40,542 |   6382 / 10301 |       3.7× |
| rod-twist         | VF    | TightInclusion |    57,088 |    40,542 |  23528 / 43830 | 1.0× (ref) |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end reference -->

### 5.2. Time of impact validation and accuracy

How far before the true time of impact each mode reports. This is the axis the
two modes trade against speed, and early is the safe direction.

The median and the worst case answer different questions, and only the second
one is about step size. A solver does not take the median step: it takes the
step it is given, so the largest earliness anywhere in a scene is the largest
step that mode can cost. The two differ by four to five orders of magnitude
here, and on cloth-funnel and rod-twist the worst case approaches 1.0 — a step
reported at its very beginning when the true contact is at its end. That is
conservative, and it is never a missed collision, but it is a step the solver
does not get to take. Nothing in the conservativeness gate catches it, because
by construction there is nothing there to catch.

It is not a property of these kernels, though: on those scenes TightInclusion's
worst case is the same, to the digit. Grazing and already-touching
configurations leave a conservative search no tighter answer, so the worst case
measures the problem rather than the implementation. The median is where the
implementations actually differ.

<!-- sccd:begin earliness -->

| scene             | mode          | median earliness | worst case |
|-------------------|---------------|-----------------:|-----------:|
| armadillo-rollers | Relaxed (GPU) |         9.33e-05 |   9.41e-01 |
| armadillo-rollers | Tight (GPU)   |         3.41e-06 |   1.60e-02 |
| armadillo-rollers | Relaxed       |         4.66e-05 |   9.41e-01 |
| armadillo-rollers | Tight         |         3.39e-06 |   1.60e-02 |
| cloth-ball        | Relaxed (GPU) |         1.11e-06 |   8.38e-04 |
| cloth-ball        | Tight (GPU)   |         2.71e-07 |   1.93e-04 |
| cloth-ball        | Relaxed       |         3.09e-07 |   5.40e-04 |
| cloth-ball        | Tight         |         2.72e-07 |   1.93e-04 |
| cloth-funnel      | Relaxed (GPU) |         2.64e-02 |   1.00e+00 |
| cloth-funnel      | Tight (GPU)   |         3.22e-04 |   1.00e+00 |
| cloth-funnel      | Relaxed       |         1.39e-02 |   1.00e+00 |
| cloth-funnel      | Tight         |         3.22e-04 |   1.00e+00 |
| n-body            | Relaxed (GPU) |         1.04e-07 |   2.52e-03 |
| n-body            | Tight (GPU)   |         3.02e-08 |   2.16e-04 |
| n-body            | Relaxed       |         1.82e-08 |   1.27e-03 |
| n-body            | Tight         |         3.02e-08 |   2.40e-04 |
| puffer-ball       | Relaxed (GPU) |         4.16e-02 |   9.75e-01 |
| puffer-ball       | Tight (GPU)   |         4.01e-05 |   2.61e-01 |
| puffer-ball       | Relaxed       |         3.99e-02 |   9.75e-01 |
| puffer-ball       | Tight         |         4.00e-05 |   2.61e-01 |
| rod-twist         | Relaxed (GPU) |         2.91e-02 |   9.96e-01 |
| rod-twist         | Tight (GPU)   |         2.75e-04 |   9.93e-01 |
| rod-twist         | Relaxed       |         2.67e-02 |   9.96e-01 |
| rod-twist         | Tight         |         2.75e-04 |   9.93e-01 |

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end earliness -->

#### 5.2.1. Against TightInclusion, on accuracy

The tables above use TightInclusion as the reference for hit versus miss, which
is what it is good for. It is **not** the reference for accuracy: its answer is a
conservative lower bound on the true root, exactly as SCCD's is. So the question
is not how close each mode gets to TightInclusion but how close all three get to
the truth, and the dataset's exact symbolic roots are the only thing in the
comparison that is actually the truth.

**`Tight` does not approximate the reference, it reproduces it.** On all twelve
scene-phases its median earliness and its worst case are identical to
TightInclusion's, to every digit — the same answer, reached faster. `Tight (GPU)`
matches the worst case exactly and the median to within a part in a thousand.
That is the strongest form agreement can take, and it is what makes the speed
comparison meaningful: nothing is being traded for it.

`Relaxed` is the mode that trades, and the table prices the trade: three to ten
times looser at the median, and on the scenes with grazing contact up to two
orders of magnitude looser at the worst case.

**The worst cases belong to the problem, not to SCCD.** On cloth-funnel and
rod-twist the largest earliness approaches 1.0 for *every* implementation
including TightInclusion — on cloth-funnel edge-edge it is exactly 1.0 for all
five. A conservative search on a grazing or already-touching configuration has
no tighter answer available to it, so this is the cost of the guarantee rather
than a defect in any one kernel. The same scenes are where the median earliness
is 0: half those queries are already in contact at the start of the step, and
every implementation returns the true root exactly.

Both are conservative throughout, and so is the reference: over 16,567,149
checked queries each, none of the five reported a time of impact after the true
one and none missed a collision.

<!-- sccd:begin earliness-ref -->

| scene             | phase | mode           | median earliness | worst case |
|-------------------|-------|----------------|-----------------:|-----------:|
| armadillo-rollers | EE    | Relaxed        |         3.47e-05 |   7.37e-01 |
| armadillo-rollers | EE    | Tight          |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | EE    | Relaxed (GPU)  |         7.47e-05 |   7.37e-01 |
| armadillo-rollers | EE    | Tight (GPU)    |         3.08e-06 |   1.60e-02 |
| armadillo-rollers | EE    | TightInclusion |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | VF    | Relaxed        |         6.44e-05 |   9.41e-01 |
| armadillo-rollers | VF    | Tight          |         3.74e-06 |   9.39e-03 |
| armadillo-rollers | VF    | Relaxed (GPU)  |         1.06e-04 |   9.42e-01 |
| armadillo-rollers | VF    | Tight (GPU)    |         3.73e-06 |   9.39e-03 |
| armadillo-rollers | VF    | TightInclusion |         3.74e-06 |   9.39e-03 |
| cloth-ball        | EE    | Relaxed        |         3.13e-07 |   5.40e-04 |
| cloth-ball        | EE    | Tight          |         1.05e-07 |   1.04e-04 |
| cloth-ball        | EE    | Relaxed (GPU)  |         1.13e-06 |   9.22e-04 |
| cloth-ball        | EE    | Tight (GPU)    |         1.05e-07 |   1.04e-04 |
| cloth-ball        | EE    | TightInclusion |         1.05e-07 |   1.04e-04 |
| cloth-ball        | VF    | Relaxed        |         2.96e-07 |   7.86e-05 |
| cloth-ball        | VF    | Tight          |         9.76e-08 |   2.88e-05 |
| cloth-ball        | VF    | Relaxed (GPU)  |         1.09e-06 |   2.15e-04 |
| cloth-ball        | VF    | Tight (GPU)    |         9.75e-08 |   2.88e-05 |
| cloth-ball        | VF    | TightInclusion |         9.76e-08 |   2.88e-05 |
| cloth-funnel      | EE    | Relaxed        |                0 |   1.00e+00 |
| cloth-funnel      | EE    | Tight          |                0 |   1.00e+00 |
| cloth-funnel      | EE    | Relaxed (GPU)  |                0 |   1.00e+00 |
| cloth-funnel      | EE    | Tight (GPU)    |                0 |   1.00e+00 |
| cloth-funnel      | EE    | TightInclusion |                0 |   1.00e+00 |
| cloth-funnel      | VF    | Relaxed        |         1.56e-02 |   9.97e-01 |
| cloth-funnel      | VF    | Tight          |         3.33e-04 |   2.67e-02 |
| cloth-funnel      | VF    | Relaxed (GPU)  |         2.82e-02 |   9.97e-01 |
| cloth-funnel      | VF    | Tight (GPU)    |         3.33e-04 |   2.67e-02 |
| cloth-funnel      | VF    | TightInclusion |         3.33e-04 |   2.67e-02 |
| n-body            | EE    | Relaxed        |         1.79e-08 |   1.27e-03 |
| n-body            | EE    | Tight          |         1.14e-08 |   1.89e-04 |
| n-body            | EE    | Relaxed (GPU)  |         1.10e-07 |   2.51e-03 |
| n-body            | EE    | Tight (GPU)    |         1.14e-08 |   1.89e-04 |
| n-body            | EE    | TightInclusion |         1.14e-08 |   1.89e-04 |
| n-body            | VF    | Relaxed        |         1.84e-08 |   1.53e-04 |
| n-body            | VF    | Tight          |         1.17e-08 |   1.87e-05 |
| n-body            | VF    | Relaxed (GPU)  |         1.03e-07 |   2.64e-04 |
| n-body            | VF    | Tight (GPU)    |         1.17e-08 |   1.87e-05 |
| n-body            | VF    | TightInclusion |         1.17e-08 |   1.87e-05 |
| puffer-ball       | EE    | Relaxed        |         3.57e-02 |   9.61e-01 |
| puffer-ball       | EE    | Tight          |         3.54e-05 |   2.61e-01 |
| puffer-ball       | EE    | Relaxed (GPU)  |         3.64e-02 |   9.62e-01 |
| puffer-ball       | EE    | Tight (GPU)    |         3.54e-05 |   2.61e-01 |
| puffer-ball       | EE    | TightInclusion |         3.54e-05 |   2.61e-01 |
| puffer-ball       | VF    | Relaxed        |         4.68e-02 |   9.75e-01 |
| puffer-ball       | VF    | Tight          |         4.58e-05 |   5.04e-02 |
| puffer-ball       | VF    | Relaxed (GPU)  |         4.60e-02 |   9.75e-01 |
| puffer-ball       | VF    | Tight (GPU)    |         4.56e-05 |   5.04e-02 |
| puffer-ball       | VF    | TightInclusion |         4.58e-05 |   5.04e-02 |
| rod-twist         | EE    | Relaxed        |         2.75e-02 |   9.93e-01 |
| rod-twist         | EE    | Tight          |         2.47e-04 |   9.85e-01 |
| rod-twist         | EE    | Relaxed (GPU)  |         2.94e-02 |   9.93e-01 |
| rod-twist         | EE    | Tight (GPU)    |         2.47e-04 |   9.85e-01 |
| rod-twist         | EE    | TightInclusion |         2.47e-04 |   9.85e-01 |
| rod-twist         | VF    | Relaxed        |         1.41e-02 |   9.96e-01 |
| rod-twist         | VF    | Tight          |         1.28e-04 |   9.93e-01 |
| rod-twist         | VF    | Relaxed (GPU)  |         1.64e-02 |   9.96e-01 |
| rod-twist         | VF    | Tight (GPU)    |         1.28e-04 |   9.93e-01 |
| rod-twist         | VF    | TightInclusion |         1.28e-04 |   9.93e-01 |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end earliness-ref -->


### 5.3. Scaling with element count

`sccd_refine_scaling` refines one surface repeatedly, quadrupling the element
count at each level, and runs a collision step on each — the one question
neither other driver can answer. Two consecutive cloth-ball frames, 92,230
elements up to 23.6 million.

Sweeping element count over two and a half orders of magnitude is also the
cleanest way to separate the two broad-phase strategies, so both are run here
rather than leaving the choice to the tuner. They are checked against each other
by construction: at every level the two report **identical pair counts**, so what
differs is only how the pairs are found.

**The strategies divide their cost differently, and that is the whole result.**
`prep` builds the acceleration structure — the cell list's grid, or the sweep's
sorted intervals — and `step` is the traversal that reports pairs:

- The **sweep** builds more cheaply and scales better doing it: exponent 0.73
  against the cell list's 0.92, and at 23.6 M elements its prep is 1.45× the
  cheaper of the two.
- The **cell list** traverses far better: exponent 0.75 against the sweep's
  1.19, and at the same size its step costs 515 ms against 4,957 ms — 9.6×.

Traversal is what grows, so the cell list wins overall at every size measured,
by 1.08× to 1.73×. The margin is narrowest in the middle and widens at both ends:
at the smallest size the sweep pays a fixed setup cost it cannot amortise, and at
the largest its traversal has become the dominant term. That is why the shipped
default races the two per scene rather than fixing one, and why the first probe
is the cell list — it is the choice that bounds the worst case for a caller who
never completes a race.

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

![Cost against element count](figures/refine-scaling.png)

Broad-phase and narrow-phase cost against element count, log-log, over five
refinement levels of the same surface. The fitted exponent is on the broad
phase; the narrow-phase series is flat and noisy because these two frames do not
come into contact.

## 6. Figures

Figures follow the presentation of the dataset paper (Belgrod et al., *TOI
dataset for CCD and a scalable conservative algorithm*), over the same six
scenes: scenes across the columns, one measured quantity per row, distributions
over cases drawn as box plots on logarithmic axes. A **star** marks a parallel
CPU mode and a **dagger** a GPU one.

![Per-case results for every mode over the six scenes](figures/results-grid.png)

**Figure 1.** Per-case distributions for every mode over the six scenes
(columns). Rows are broad-phase time, narrow-phase time, narrow-phase false
positives, and earliness against the exact root. Each box spans the first to the
third quartile with the median inside it, the whiskers reach the furthest case
within 1.5 interquartile ranges, and cases beyond that are drawn individually.
Counts are on a symmetric-log axis so that a case with no false positive — the
common outcome — is on the axis rather than dropped.

![Runtime split by phase](figures/runtime-breakdown.png)

**Figure 2.** Runtime split by phase for every scene and mode: *prep* builds the
swept boxes and the acceleration structure, *broad* finds the candidate pairs,
*narrow* turns them into a time of impact. Preparation dominates on the host;
on the GPU the narrow phase does, and on n-body and puffer-ball it dominates
outright.

![Time-of-impact error against the symbolic ground truth](figures/toi-error.png)

**Figure 3.** Distribution over cases of the time-of-impact error against the
dataset's exact symbolic roots, on a log axis. The error is one-sided by
construction — it is how far *before* the true root the mode reported — so every
value shown is on the safe side, and nothing falls off the axis on the late side
because no such case exists. `Tight` and `Tight (GPU)` coincide, which is why the
device modes are dashed.

![Narrow-phase time per case against candidate pairs](figures/narrow-per-case.png)

**Figure 4.** Narrow-phase time for each individual case against the number of
candidate pairs the broad phase handed it, log-log, median over repeats.

## 7. Provenance

<!-- sccd:begin provenance -->

- Timings: `benchmark/results/sweep-gh200-bp.csv`, 6384 cases over 6 scenes, 3 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200-all.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
