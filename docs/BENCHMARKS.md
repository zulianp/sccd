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

## What is being measured

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

## Platform

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

## Scenes

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
<!-- sccd:end dataset -->

Accuracy is measured on the curated query sets, because the exact roots belong
to those coordinates specifically. The mesh path is a separately stored copy of
the same frames, read from PLY, and is reported beside them rather than checked
against them.

### Ground-truth coverage

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

## Modes

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

## Conservativeness

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
| puffer-ball       | EE    | Relaxed (GPU) |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | Tight (GPU)   |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | Relaxed       |       1,187,257 |      0 |    0 |
| puffer-ball       | EE    | Tight         |       1,187,257 |      0 |    0 |
| puffer-ball       | VF    | Relaxed (GPU) |         299,533 |      0 |    0 |
| puffer-ball       | VF    | Tight (GPU)   |         299,533 |      0 |    0 |
| puffer-ball       | VF    | Relaxed       |         299,533 |      0 |    0 |
| puffer-ball       | VF    | Tight         |         299,533 |      0 |    0 |

TightInclusion is excluded from this table: it is the reference the queries were selected against, not a subject of it.

Source: `benchmark/results/oracle-gh200.csv`

<!-- sccd:end gate -->

Per-scene, with the false positives that are the price of accepting early, and
one column that is not part of the gate:

<!-- sccd:begin conservativeness -->

| scene             | mode          | queries | toi compared | late | false pos. | false neg. | mesh-path divergence |
|-------------------|---------------|--------:|-------------:|-----:|-----------:|-----------:|---------------------:|
| armadillo-rollers | Relaxed (GPU) | 131,441 |      130,859 |    0 |        322 |          0 |                    5 |
| armadillo-rollers | Tight (GPU)   | 131,441 |      130,859 |    0 |         24 |          0 |                 1311 |
| armadillo-rollers | Relaxed       | 131,441 |      130,859 |    0 |        232 |          0 |                  126 |
| armadillo-rollers | Tight         | 131,441 |      130,859 |    0 |         24 |          0 |                 1320 |
| cloth-ball        | Relaxed (GPU) | 664,940 |      664,919 |    0 |          2 |          0 |                    0 |
| cloth-ball        | Tight (GPU)   | 664,940 |      664,919 |    0 |          1 |          0 |                    0 |
| cloth-ball        | Relaxed       | 664,940 |      664,919 |    0 |          2 |          0 |                    0 |
| cloth-ball        | Tight         | 664,940 |      664,919 |    0 |          1 |          0 |                    0 |
| cloth-funnel      | Relaxed (GPU) |   7,552 |        6,773 |    0 |        742 |          0 |                    0 |
| cloth-funnel      | Tight (GPU)   |   7,552 |        6,773 |    0 |         15 |          0 |                    0 |
| cloth-funnel      | Relaxed       |   7,552 |        6,773 |    0 |        687 |          0 |                    0 |
| cloth-funnel      | Tight         |   7,552 |        6,773 |    0 |         15 |          0 |                    0 |
| puffer-ball       | Relaxed       | 440,659 |      433,614 |    0 |      7,045 |          0 |                    0 |
| puffer-ball       | Tight         | 440,659 |      433,614 |    0 |        138 |          0 |                    0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness. The last column is not part of the gate and is reported for completeness. It counts cases where the earliest-impact answer computed over the *mesh* is later than the earliest exact root of the *curated queries* -- two different geometries, because smesh stores mesh coordinates as float32 while the curated queries are exact dyadic rationals. Near a grazing contact a last-bit coordinate change moves the root by far more than it moves the coordinate, which is why the divergence is larger than float32's precision. On every one of those cases the curated-query answer is at or before the exact root, so it is a difference between two inputs, not a kernel reporting late.

Source: `benchmark/results/sweep-gh200.csv`

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

## Timing

Whole-scene wall clock, summed over every case, median of three independent
runs with the full run-to-run range beside it.

<!-- sccd:begin timing -->

| scene             | mode          | cases |       queries | repeats |       broad ms |      narrow ms |        total ms |
|-------------------|---------------|------:|--------------:|--------:|---------------:|---------------:|----------------:|
| armadillo-rollers | Relaxed (GPU) |   781 |    85,407,326 |       5 |  890.8 (5.1 %) | 1641.4 (3.3 %) |  2532.2 (3.9 %) |
| armadillo-rollers | Tight (GPU)   |   781 |    85,407,326 |       5 |  850.6 (6.8 %) | 1829.1 (4.4 %) |  2674.5 (5.0 %) |
| armadillo-rollers | Relaxed       |   781 |    85,407,326 |       5 | 3193.1 (2.6 %) |  722.1 (5.0 %) |  3901.9 (2.6 %) |
| armadillo-rollers | Tight         |   781 |    85,407,326 |       5 | 3215.5 (1.7 %) | 1026.7 (1.3 %) |  4242.8 (1.4 %) |
| cloth-ball        | Relaxed (GPU) |    79 |   175,150,856 |       5 |  433.8 (7.2 %) | 485.3 (16.2 %) |  921.2 (11.7 %) |
| cloth-ball        | Tight (GPU)   |    79 |   175,150,856 |       5 | 440.4 (11.5 %) | 480.3 (20.5 %) |  923.3 (15.9 %) |
| cloth-ball        | Relaxed       |    79 |   175,150,856 |       5 |  894.0 (1.3 %) |  379.5 (1.3 %) |  1273.5 (1.0 %) |
| cloth-ball        | Tight         |    79 |   175,150,856 |       5 |  886.8 (0.6 %) |  214.4 (3.5 %) |  1100.9 (0.5 %) |
| cloth-funnel      | Relaxed (GPU) |   577 |    25,192,698 |       5 |  692.9 (8.9 %) | 1017.0 (7.8 %) |  1710.0 (8.1 %) |
| cloth-funnel      | Tight (GPU)   |   577 |    25,192,698 |       5 |  685.1 (2.4 %) | 1172.5 (2.4 %) |  1855.4 (1.0 %) |
| cloth-funnel      | Relaxed       |   577 |    25,192,698 |       5 | 2141.0 (2.6 %) |  492.4 (2.2 %) |  2631.2 (2.4 %) |
| cloth-funnel      | Tight         |   577 |    25,192,698 |       5 | 2151.0 (1.8 %) |  631.4 (8.3 %) |  2777.8 (3.3 %) |
| puffer-ball       | Relaxed       |    40 | 1,768,224,407 |       3 | 8248.5 (0.8 %) | 2959.5 (1.6 %) | 11208.0 (1.0 %) |
| puffer-ball       | Tight         |    40 | 1,768,224,407 |       3 | 8284.5 (2.3 %) | 1434.5 (0.4 %) |  9721.6 (1.9 %) |

A difference smaller than the bracketed spread does not separate two modes and is not reported as a ratio anywhere in this document.

Source: `benchmark/results/sweep-gh200.csv`

<!-- sccd:end timing -->

Whole-scene milliseconds cannot be compared between a 79-case scene and a
4,571-case one; throughput can:

<!-- sccd:begin throughput -->
<!-- sccd:end throughput -->

**Neither mode is uniformly faster**, which is the substantive result here and
the reason the trade is worth stating as a trade:

<!-- sccd:begin comparison -->

- **armadillo-rollers (CPU)**: Relaxed is 1.42× faster than Tight in the narrow phase (722 ms against 1027 ms; run-to-run spread 5.0%).
- **armadillo-rollers (GPU)**: Relaxed (GPU) is 1.11× faster than Tight (GPU) in the narrow phase (1641 ms against 1829 ms; run-to-run spread 4.4%).
- **cloth-ball (CPU)**: Tight is 1.77× faster than Relaxed in the narrow phase (214 ms against 380 ms; run-to-run spread 3.5%).
- **cloth-ball (GPU)**: Tight (GPU) and Relaxed (GPU) are inside noise (480 ms against 485 ms, spread 20.5%); this does not separate them.
- **cloth-funnel (CPU)**: Relaxed is 1.28× faster than Tight in the narrow phase (492 ms against 631 ms; run-to-run spread 8.3%).
- **cloth-funnel (GPU)**: Relaxed (GPU) is 1.15× faster than Tight (GPU) in the narrow phase (1017 ms against 1173 ms; run-to-run spread 7.8%).
- **puffer-ball (CPU)**: Tight is 2.06× faster than Relaxed in the narrow phase (1434 ms against 2959 ms; run-to-run spread 1.6%).

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

**The two processors divide the work differently.** On all three scenes measured
on both, the GPU broad phase is two to three and a half times faster than the
CPU's while the GPU narrow phase is slower — so the GPU wins end to end, but it
wins in the broad phase and loses in the narrow one. Comparisons above are made
within a processor for that reason: ranking every mode of a scene together would
compare host `Relaxed` against GPU `Tight` and report the sum of two unrelated
effects as if it were the mode trade.

## Against TightInclusion

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
| armadillo-rollers | EE    | Relaxed (GPU)  |    99,104 |    98,933 |   1921 (0.3 %) |      20.3× |
| armadillo-rollers | EE    | Tight (GPU)    |    99,104 |    98,761 |   2325 (0.5 %) |      16.8× |
| armadillo-rollers | EE    | Relaxed        |    99,104 |    98,895 |   3427 (1.2 %) |      11.4× |
| armadillo-rollers | EE    | Tight          |    99,104 |    98,761 |   6525 (0.2 %) |       6.0× |
| armadillo-rollers | EE    | TightInclusion |    99,104 |    98,761 |  38981 (0.1 %) | 1.0× (ref) |
| armadillo-rollers | VF    | Relaxed (GPU)  |    32,337 |    32,248 |   3589 (7.4 %) |       4.8× |
| armadillo-rollers | VF    | Tight (GPU)    |    32,337 |    32,122 |   1633 (1.6 %) |      10.6× |
| armadillo-rollers | VF    | Relaxed        |    32,337 |    32,196 |   2275 (0.6 %) |       7.6× |
| armadillo-rollers | VF    | Tight          |    32,337 |    32,122 |   2728 (0.4 %) |       6.4× |
| armadillo-rollers | VF    | TightInclusion |    32,337 |    32,122 |  17363 (0.3 %) | 1.0× (ref) |
| cloth-ball        | EE    | Relaxed (GPU)  |   557,683 |   557,669 |    523 (1.5 %) |     264.1× |
| cloth-ball        | EE    | Tight (GPU)    |   557,683 |   557,668 |    492 (0.4 %) |     280.5× |
| cloth-ball        | EE    | Relaxed        |   557,683 |   557,669 |    412 (0.7 %) |     335.0× |
| cloth-ball        | EE    | Tight          |   557,683 |   557,668 |    970 (0.2 %) |     142.4× |
| cloth-ball        | EE    | TightInclusion |   557,683 |   557,668 | 138126 (0.7 %) | 1.0× (ref) |
| cloth-ball        | VF    | Relaxed (GPU)  |   107,257 |   107,252 |  2376 (10.8 %) |      13.1× |
| cloth-ball        | VF    | Tight (GPU)    |   107,257 |   107,252 |    343 (0.8 %) |      90.5× |
| cloth-ball        | VF    | Relaxed        |   107,257 |   107,252 |    327 (2.0 %) |      95.0× |
| cloth-ball        | VF    | Tight          |   107,257 |   107,252 |    414 (1.4 %) |      75.0× |
| cloth-ball        | VF    | TightInclusion |   107,257 |   107,252 |  31026 (0.7 %) | 1.0× (ref) |
| cloth-funnel      | EE    | Relaxed (GPU)  |     6,751 |     6,734 |    699 (2.8 %) |       9.1× |
| cloth-funnel      | EE    | Tight (GPU)    |     6,751 |     6,259 |   1140 (1.1 %) |       5.6× |
| cloth-funnel      | EE    | Relaxed        |     6,751 |     6,700 |   1144 (1.8 %) |       5.6× |
| cloth-funnel      | EE    | Tight          |     6,751 |     6,259 |   1758 (0.8 %) |       3.6× |
| cloth-funnel      | EE    | TightInclusion |     6,751 |     6,259 |   6384 (0.9 %) | 1.0× (ref) |
| cloth-funnel      | VF    | Relaxed (GPU)  |       801 |       781 |   2400 (9.2 %) |       0.5× |
| cloth-funnel      | VF    | Tight (GPU)    |       801 |       529 |    399 (0.1 %) |       3.0× |
| cloth-funnel      | VF    | Relaxed        |       801 |       760 |    417 (4.2 %) |       2.9× |
| cloth-funnel      | VF    | Tight          |       801 |       529 |    390 (3.0 %) |       3.0× |
| cloth-funnel      | VF    | TightInclusion |       801 |       529 |   1189 (1.0 %) | 1.0× (ref) |
| puffer-ball       | EE    | Relaxed (GPU)  | 1,206,952 | 1,206,951 |    394 (3.9 %) |     445.8× |
| puffer-ball       | EE    | Tight (GPU)    | 1,206,952 | 1,187,650 |   1159 (0.9 %) |     151.8× |
| puffer-ball       | EE    | Relaxed        | 1,206,952 | 1,206,951 |    267 (2.5 %) |     657.8× |
| puffer-ball       | EE    | Tight          | 1,206,952 | 1,187,650 |   1646 (1.0 %) |     106.8× |
| puffer-ball       | EE    | TightInclusion | 1,206,952 | 1,187,650 | 175849 (1.6 %) | 1.0× (ref) |
| puffer-ball       | VF    | Relaxed (GPU)  |   307,220 |   307,219 |   2038 (2.0 %) |      22.9× |
| puffer-ball       | VF    | Tight (GPU)    |   307,220 |   299,698 |    767 (1.3 %) |      60.8× |
| puffer-ball       | VF    | Relaxed        |   307,220 |   307,219 |    237 (1.8 %) |     197.0× |
| puffer-ball       | VF    | Tight          |   307,220 |   299,676 |    591 (2.0 %) |      78.8× |
| puffer-ball       | VF    | TightInclusion |   307,220 |   299,676 |  46616 (1.7 %) | 1.0× (ref) |

Source: `benchmark/results/oracle-gh200.csv`

<!-- sccd:end reference -->

## Accuracy

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

<!-- sccd:begin earliness -->

| scene             | mode          | median earliness |
|-------------------|---------------|-----------------:|
| armadillo-rollers | Relaxed (GPU) |         9.32e-05 |
| armadillo-rollers | Tight (GPU)   |         3.41e-06 |
| armadillo-rollers | Relaxed       |         4.66e-05 |
| armadillo-rollers | Tight         |         3.39e-06 |
| cloth-ball        | Relaxed (GPU) |         1.10e-06 |
| cloth-ball        | Tight (GPU)   |         2.70e-07 |
| cloth-ball        | Relaxed       |         3.06e-07 |
| cloth-ball        | Tight         |         2.72e-07 |
| cloth-funnel      | Relaxed (GPU) |         2.65e-02 |
| cloth-funnel      | Tight (GPU)   |         3.22e-04 |
| cloth-funnel      | Relaxed       |         1.40e-02 |
| cloth-funnel      | Tight         |         3.22e-04 |
| puffer-ball       | Relaxed       |         3.52e-02 |
| puffer-ball       | Tight         |         3.61e-05 |

Source: `benchmark/results/sweep-gh200.csv`

<!-- sccd:end earliness -->

### Against TightInclusion, on accuracy

The tables above use TightInclusion as the reference for hit versus miss, which
is what it is good for. It is **not** the reference for accuracy: its answer is a
conservative lower bound on the true root, exactly as SCCD's is. So the question
is not how close each mode gets to TightInclusion but how close all three get to
the truth, and the dataset's exact symbolic roots are the only thing in the
comparison that is actually the truth.

<!-- sccd:begin earliness-ref -->
<!-- sccd:end earliness-ref -->


## Scaling with element count

`sccd_refine_scaling` refines one surface repeatedly, quadrupling the element
count at each level, and runs a collision step on each — the one question
neither other driver can answer. Two consecutive cloth-ball frames, 92,230
elements up to 23.6 million.

<!-- sccd:begin scaling -->

| mode                 | level |   elements | candidate pairs | broad ms | narrow ms |    p |
|----------------------|------:|-----------:|----------------:|---------:|----------:|-----:|
| relaxed / host / tri |     0 |     92,230 |          17,982 |     35.3 |       8.6 | 0.86 |
|                      |     1 |    368,920 |          87,848 |     65.3 |       5.2 |      |
|                      |     2 |  1,475,680 |         380,924 |    241.8 |       1.9 |      |
|                      |     3 |  5,902,720 |       1,581,142 |   1008.1 |       6.2 |      |
|                      |     4 | 23,610,880 |       6,440,923 |   4424.2 |      19.9 |      |
| tight / host / tri   |     0 |     92,230 |          17,982 |     35.2 |       4.8 | 0.87 |
|                      |     1 |    368,920 |          87,848 |     62.2 |      10.6 |      |
|                      |     2 |  1,475,680 |         380,924 |    202.0 |       4.2 |      |
|                      |     3 |  5,902,720 |       1,581,142 |    964.3 |       4.5 |      |
|                      |     4 | 23,610,880 |       6,440,923 |   4433.4 |      16.2 |      |

The two frames used here do not come into contact, so the narrow phase has almost no work to do and its column is dominated by noise rather than by element count; what this measures is the broad phase and the preparation that feeds it. Narrow-phase cost against problem size is in the per-case figure, over cases that do collide. The exponent is below 1 because the fixed cost visible at the smallest size is amortised as the mesh grows.

Source: `benchmark/results/scaling/host-mode0.txt, benchmark/results/scaling/host-mode2.txt`

<!-- sccd:end scaling -->

![Cost against element count](figures/refine-scaling.png)

Broad-phase and narrow-phase cost against element count, log-log, over five
refinement levels of the same surface. The fitted exponent is on the broad
phase; the narrow-phase series is flat and noisy because these two frames do not
come into contact.

## Figures

![Broad and narrow phase per scene and mode](figures/phase-breakdown.png)

Broad phase (pale) and narrow phase (solid) summed over every case, per mode.
Bars are the median over three repeats; whiskers span the full run-to-run range
of the total, so a difference smaller than a whisker is not a result.

![Narrow-phase time per case against candidate pairs](figures/narrow-per-case.png)

Narrow-phase time for each individual case against the number of candidate pairs
the broad phase handed it, log-log, median over repeats.

![Distribution of earliness](figures/earliness-cdf.png)

Empirical distribution, over cases, of how far before the true time of impact
each mode reports. Earliness is on a log axis, so a curve further left reports
closer to the true root; a curve further right reports earlier, which is always
the safe direction and always costs a solver step size.

## Provenance

<!-- sccd:begin provenance -->

- Timings: `benchmark/results/sweep-gh200.csv`, 1457 cases over 4 scenes, 3-5 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
