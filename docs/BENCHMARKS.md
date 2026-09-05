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

**A difference smaller than the run-to-run spread is not a result.** Timing
tables carry that spread beside every median, figures draw it, and the mode
comparisons below refuse to state a ratio when the gap is inside it.

## Platform

One GH200 node on CSCS Alps: Grace (72 threads, `OMP_NUM_THREADS` set to the
node's core count) and Hopper sm_90 in the same allocation, built with
`prgenv-gnu/24.11:v2`, `-O3`, `CMAKE_CUDA_ARCHITECTURES=90`. Search parameters
are the shipped defaults.

## Scenes

The dataset is the NYU CCD benchmark (archive 2451/74508), real simulation
frames rather than synthetic motion. Each case is one step, with a curated query
set whose coordinates are exact dyadic rationals and whose roots were computed
symbolically.

Three of the six scenes are swept. The other three are not, and the reason is
recorded here rather than left to be inferred:

| scene | swept | why not |
|---|---|---|
| armadillo-rollers | yes | — |
| cloth-ball | yes | — |
| cloth-funnel | yes | — |
| puffer-ball | no | ships boxes, queries and roots but no extracted frames, so it has zero runnable cases; its 240 root archives are also unconverted |
| n-body-simulation | no | not downloaded |
| rod-twist | no | not downloaded |

Accuracy is measured on the curated query sets and cannot be measured through
the mesh path: smesh stores coordinates as `float`, so mesh geometry is a
rounded copy of the geometry those exact roots belong to.

### Ground-truth coverage

The dataset states each query's outcome twice — a boolean in `mma_bool` and a
time of impact in `roots` — and `benchmark/verify_oracle.py` requires them to
agree query for query, `mma_bool[i] == isfinite(toi[i])`. A `true` against a NaN
is a root that never converted, and downstream it is indistinguishable from "no
collision", so an incomplete oracle silently shrinks the evidence instead of
announcing itself. The gate runs as part of `prepare_data.sh` and exits non-zero
on any gap.

Coverage on the swept scenes is complete: 130,859 of armadillo-rollers' 131,441
queries carry a root, 664,919 of cloth-ball's 664,940, and 6,773 of
cloth-funnel's 7,552 — in each case exactly the set `mma_bool` marks as
colliding.

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
impact in all twenty-four configurations.**

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

TightInclusion is excluded from this table: it is the reference the queries were selected against, not a subject of it.

Source: `benchmark/results/oracle-gh200.csv`

<!-- sccd:end gate -->

Per-scene, with the false positives that are the price of accepting early, and
one column that is not part of the gate:

<!-- sccd:begin conservativeness -->

| scene             | mode    | queries | toi compared | late | false pos. | false neg. | mesh-path divergence |
|-------------------|---------|--------:|-------------:|-----:|-----------:|-----------:|---------------------:|
| armadillo-rollers | Relaxed | 131,441 |      130,859 |    0 |        232 |          0 |                  126 |
| armadillo-rollers | Tight   | 131,441 |      130,859 |    0 |         24 |          0 |                 1320 |
| cloth-ball        | Relaxed | 664,940 |      664,919 |    0 |          2 |          0 |                    0 |
| cloth-ball        | Tight   | 664,940 |      664,919 |    0 |          1 |          0 |                    0 |
| cloth-funnel      | Relaxed |   7,552 |        6,773 |    0 |        687 |          0 |                    0 |
| cloth-funnel      | Tight   |   7,552 |        6,773 |    0 |         15 |          0 |                    0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness. The last column is not part of the gate and is reported for completeness. It counts cases where the earliest-impact answer computed over the *mesh* is later than the earliest exact root of the *curated queries* -- two different geometries, because smesh stores mesh coordinates as float32 while the curated queries are exact dyadic rationals. Near a grazing contact a last-bit coordinate change moves the root by far more than it moves the coordinate, which is why the divergence is larger than float32's precision. On every one of those cases the curated-query answer is at or before the exact root, so it is a difference between two inputs, not a kernel reporting late.

Source: `benchmark/results/sweep-gh200-host.csv`

<!-- sccd:end conservativeness -->

The last column deserves its own sentence, because it looks like a violation and
is not. It counts cases where the earliest-impact answer computed over the
**mesh** lands after the earliest exact root of the **curated queries**. Those are
two different geometries: smesh stores mesh coordinates as `float32`, while the
curated queries are exact dyadic rationals. On every one of those cases the
curated-query answer is at or before the exact root. Near a grazing contact a
last-bit change in a coordinate moves the root by far more than it moves the
coordinate, which is why the divergence is larger than `float32`'s own precision.
It measures the distance between two inputs, not a kernel reporting late, so it
is reported and not gated on.

## Timing

Whole-scene wall clock, summed over every case, median of five independent runs
with the full run-to-run range beside it.

<!-- sccd:begin timing -->

| scene             | mode    | cases |     queries | repeats |       broad ms |      narrow ms |       total ms |
|-------------------|---------|------:|------------:|--------:|---------------:|---------------:|---------------:|
| armadillo-rollers | Relaxed |   781 |  85,407,326 |       5 | 3193.1 (2.6 %) |  722.1 (5.0 %) | 3901.9 (2.6 %) |
| armadillo-rollers | Tight   |   781 |  85,407,326 |       5 | 3215.5 (1.7 %) | 1026.7 (1.3 %) | 4242.8 (1.4 %) |
| cloth-ball        | Relaxed |    79 | 175,150,856 |       5 |  894.0 (1.3 %) |  379.5 (1.3 %) | 1273.5 (1.0 %) |
| cloth-ball        | Tight   |    79 | 175,150,856 |       5 |  886.8 (0.6 %) |  214.4 (3.5 %) | 1100.9 (0.5 %) |
| cloth-funnel      | Relaxed |   577 |  25,192,698 |       5 | 2141.0 (2.6 %) |  492.4 (2.2 %) | 2631.2 (2.4 %) |
| cloth-funnel      | Tight   |   577 |  25,192,698 |       5 | 2151.0 (1.8 %) |  631.4 (8.3 %) | 2777.8 (3.3 %) |

A difference smaller than the bracketed spread does not separate two modes and is not reported as a ratio anywhere in this document.

Source: `benchmark/results/sweep-gh200-host.csv`

<!-- sccd:end timing -->

**Neither mode is uniformly faster**, which is the substantive result here and
the reason the trade is worth stating as a trade:

<!-- sccd:begin comparison -->

- **armadillo-rollers**: Relaxed is 1.42× faster than Tight in the narrow phase (722 ms against 1027 ms; run-to-run spread 5.0%).
- **cloth-ball**: Tight is 1.77× faster than Relaxed in the narrow phase (214 ms against 380 ms; run-to-run spread 3.5%).
- **cloth-funnel**: Relaxed is 1.28× faster than Tight in the narrow phase (492 ms against 631 ms; run-to-run spread 8.3%).

<!-- sccd:end comparison -->

`Relaxed` accepts a box sooner, so it does less work per query — but a looser
bound also prunes less, so the queries after it do more. Which effect wins is a
property of the scene, not of the mode.

## Against TightInclusion

TightInclusion is the reference implementation of a certified conservative
narrow phase. Matching its hit count is the strongest agreement available;
reporting more hits is a false positive, which costs work and is never unsafe.
Its own answer is a lower bound on the truth rather than the truth, which is why
the conservativeness table above is measured against the exact roots instead.

<!-- sccd:begin reference -->

| scene             | phase | mode           | queries |    hits |        time ms |     vs. TI |
|-------------------|-------|----------------|--------:|--------:|---------------:|-----------:|
| armadillo-rollers | EE    | Relaxed (GPU)  |  99,104 |  98,933 |   1921 (0.3 %) |      20.3× |
| armadillo-rollers | EE    | Tight (GPU)    |  99,104 |  98,761 |   2325 (0.5 %) |      16.8× |
| armadillo-rollers | EE    | Relaxed        |  99,104 |  98,895 |   3427 (1.2 %) |      11.4× |
| armadillo-rollers | EE    | Tight          |  99,104 |  98,761 |   6525 (0.2 %) |       6.0× |
| armadillo-rollers | EE    | TightInclusion |  99,104 |  98,761 |  38981 (0.1 %) | 1.0× (ref) |
| armadillo-rollers | VF    | Relaxed (GPU)  |  32,337 |  32,248 |   3589 (7.4 %) |       4.8× |
| armadillo-rollers | VF    | Tight (GPU)    |  32,337 |  32,122 |   1633 (1.6 %) |      10.6× |
| armadillo-rollers | VF    | Relaxed        |  32,337 |  32,196 |   2275 (0.6 %) |       7.6× |
| armadillo-rollers | VF    | Tight          |  32,337 |  32,122 |   2728 (0.4 %) |       6.4× |
| armadillo-rollers | VF    | TightInclusion |  32,337 |  32,122 |  17363 (0.3 %) | 1.0× (ref) |
| cloth-ball        | EE    | Relaxed (GPU)  | 557,683 | 557,669 |    523 (1.5 %) |     264.1× |
| cloth-ball        | EE    | Tight (GPU)    | 557,683 | 557,668 |    492 (0.4 %) |     280.5× |
| cloth-ball        | EE    | Relaxed        | 557,683 | 557,669 |    412 (0.7 %) |     335.0× |
| cloth-ball        | EE    | Tight          | 557,683 | 557,668 |    970 (0.2 %) |     142.4× |
| cloth-ball        | EE    | TightInclusion | 557,683 | 557,668 | 138126 (0.7 %) | 1.0× (ref) |
| cloth-ball        | VF    | Relaxed (GPU)  | 107,257 | 107,252 |  2376 (10.8 %) |      13.1× |
| cloth-ball        | VF    | Tight (GPU)    | 107,257 | 107,252 |    343 (0.8 %) |      90.5× |
| cloth-ball        | VF    | Relaxed        | 107,257 | 107,252 |    327 (2.0 %) |      95.0× |
| cloth-ball        | VF    | Tight          | 107,257 | 107,252 |    414 (1.4 %) |      75.0× |
| cloth-ball        | VF    | TightInclusion | 107,257 | 107,252 |  31026 (0.7 %) | 1.0× (ref) |
| cloth-funnel      | EE    | Relaxed (GPU)  |   6,751 |   6,734 |    699 (2.8 %) |       9.1× |
| cloth-funnel      | EE    | Tight (GPU)    |   6,751 |   6,259 |   1140 (1.1 %) |       5.6× |
| cloth-funnel      | EE    | Relaxed        |   6,751 |   6,700 |   1144 (1.8 %) |       5.6× |
| cloth-funnel      | EE    | Tight          |   6,751 |   6,259 |   1758 (0.8 %) |       3.6× |
| cloth-funnel      | EE    | TightInclusion |   6,751 |   6,259 |   6384 (0.9 %) | 1.0× (ref) |
| cloth-funnel      | VF    | Relaxed (GPU)  |     801 |     781 |   2400 (9.2 %) |       0.5× |
| cloth-funnel      | VF    | Tight (GPU)    |     801 |     529 |    399 (0.1 %) |       3.0× |
| cloth-funnel      | VF    | Relaxed        |     801 |     760 |    417 (4.2 %) |       2.9× |
| cloth-funnel      | VF    | Tight          |     801 |     529 |    390 (3.0 %) |       3.0× |
| cloth-funnel      | VF    | TightInclusion |     801 |     529 |   1189 (1.0 %) | 1.0× (ref) |

Source: `benchmark/results/oracle-gh200.csv`

<!-- sccd:end reference -->

## Accuracy

How far before the true time of impact each mode reports. This is the axis the
two modes trade against speed, and early is the safe direction.

<!-- sccd:begin earliness -->

| scene             | mode    | median earliness |
|-------------------|---------|-----------------:|
| armadillo-rollers | Relaxed |         4.66e-05 |
| armadillo-rollers | Tight   |         3.39e-06 |
| cloth-ball        | Relaxed |         3.06e-07 |
| cloth-ball        | Tight   |         2.72e-07 |
| cloth-funnel      | Relaxed |         1.40e-02 |
| cloth-funnel      | Tight   |         3.22e-04 |

Source: `benchmark/results/sweep-gh200-host.csv`

<!-- sccd:end earliness -->

## Figures

![Broad and narrow phase per scene and mode](figures/phase-breakdown.png)

Broad phase (pale) and narrow phase (solid) summed over every case, per mode.
Bars are the median over five repeats; whiskers span the full run-to-run range of
the total, so a difference smaller than a whisker is not a result.

![Narrow-phase time per case against candidate pairs](figures/narrow-per-case.png)

Narrow-phase time for each individual case against the number of candidate pairs
the broad phase handed it, log-log, median over repeats.

![Distribution of earliness](figures/earliness-cdf.png)

How far before the true time of impact each mode reports, over cases. A curve
further left is less accurate and never unsafe.

## Provenance

<!-- sccd:begin provenance -->

- Timings: `benchmark/results/sweep-gh200-host.csv`, 1437 cases over 3 scenes, 5 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
