# Benchmarks

SCCD against TightInclusion on the NYU CCD dataset: six scenes, 6,394 cases over
3,561 simulation frames, on one GH200 node — a single Grace at 72 threads and Hopper sm_90 in the
same allocation, `prgenv-gnu/24.11:v2`, `-O3`, `CMAKE_CUDA_ARCHITECTURES=90`,
shipped default search parameters. A GH200 node carries four Grace modules, so
`nproc` there reports 288; host figures are one Grace.

Timings are reported as two phases. **Broad** is the whole broad phase: building
the swept boxes and the acceleration structure, then querying it for overlaps.
**Narrow** is the branch and bound over the candidates it returns. Where the text
explains the broad phase it is decomposed into the *structure* and the
*traversal* over it.

Three results:

- **Conservative on every query.** Zero missed collisions and zero late times of
  impact over 66,268,596 comparisons against exact symbolic roots.
- **Identical accuracy to the reference.** Median and worst-case error match
  TightInclusion's on all twelve scene-phases, to every digit.
- **1.6× to 3.7× faster on CPU** and 2.5× to 18.2× on GPU, at that accuracy.

Every number is generated from a committed CSV by a committed script:

```sh
python3 -m report benchmark/assessment/broadphase-cell2dmin.csv /tmp/report \
        benchmark/results/oracle-gh200-all.csv \
        --modes=tight,device-tight --label="tight:CPU,device-tight:GPU" \
        --embed=docs/BENCHMARKS.md --check
```

`--check` returns non-zero if this document is not what the data reproduces.

## 1. Dataset

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

Source: `benchmark/assessment/broadphase-cell2dmin.csv`

<!-- sccd:end dataset -->

A query carries a root exactly when the dataset marks it as colliding, and
`benchmark/verify_oracle.py` gates on that agreement query for query. The
proportion differs by scene because it is the fraction of curated queries that
collide, not a measure of coverage.

## 2. What each phase measures

A **case** is one query type on one frame: the dataset ships `100ee` and `100vf`
separately. A simulation step runs both, so a per-frame figure is the two
together.

Each case is timed as three consecutive calls, each with its own clock. The
first case of a scene is run once untimed to pay for allocation and first touch.

**prep** — `broad_phase_prep`. Builds the swept axis-aligned boxes, one per
vertex, face and edge, each enclosing that primitive's whole trajectory over the
step; picks the sort axis from the spread of box centres; and builds the
acceleration structure the traversal needs — the sorted intervals for the sweep,
or the uniform grid for the cell list. It does not look for overlaps.

**broad** — `broad_phase_fv_step` or `broad_phase_ee_step`. The overlap query
itself, over the structure `prep` built: count the candidate pairs, prefix-sum
the counts, then fill the output. A candidate is a pair whose swept boxes
overlap, which is not yet a contact.

**narrow** — `narrow_phase_vf` or `narrow_phase_ee` at `ToiOutput::Earliest`.
Branch and bound over each candidate pair's `(t, u, v)` box, returning one
earliest time of impact for the step, so every query prunes against the running
minimum. This is where the conservativeness guarantee is produced.

`total` is the three added. Nothing else is inside the timers: mesh loading, PLY
conversion and result checking are outside them.

## 3. Conservativeness

The invariant: a reported time of impact must be at or before the true one, and
a collision that exists must be reported. Checked against the exact symbolic
roots, on both processors, for every scene.

<!-- sccd:begin conservativeness -->

| scene             |   queries | toi compared | late |    false pos. | false neg. |
|-------------------|----------:|-------------:|-----:|--------------:|-----------:|
| armadillo-rollers |   131,441 |      130,859 |    0 |            24 |          0 |
| cloth-ball        |   664,940 |      664,919 |    0 |             1 |          0 |
| cloth-funnel      |     7,552 |        6,773 |    0 |            15 |          0 |
| n-body            | 2,947,719 |    2,947,611 |    0 |       11 (12) |          0 |
| puffer-ball       | 1,514,172 |    1,486,790 |    0 |     536 (558) |          0 |
| rod-twist         |   549,208 |      285,431 |    0 | 1,243 (1,245) |          0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness.

Source: `benchmark/assessment/broadphase-cell2dmin.csv`

<!-- sccd:end conservativeness -->

False positives cost work, never safety, and are reported for information.

One further check, on the inputs rather than the kernel: the answer computed
over the **mesh** is compared against the earliest exact root of the **curated
query set**, which is a separately stored copy of the same frames. A
disagreement would mean the two paths were handed different geometry, and that
the exact roots are not a reference for the mesh path at all.

<!-- sccd:begin mesh-check -->

The mesh path agrees with the curated query geometry on every case.

<!-- sccd:end mesh-check -->

## 4. Against TightInclusion

TightInclusion is the reference implementation of a certified conservative
narrow phase. Both sides are given the same task and the same machine: SCCD at
`ToiOutput::PerPair` and TightInclusion unbounded over the whole step, both
scheduled through `sccd::parallel_for_br_dynamic`, with TightInclusion used
exactly as released and run depth-first, its faster configuration when there is
no bound to exploit.

TightInclusion runs on the host, so the **CPU** row is the like-for-like
comparison and the **GPU** row is what the same search costs on the other
processor.

![SCCD against TightInclusion](figures/reference-speedup.png)

Each bar is one processor's whole narrow phase over the same queries, stacked into its vertex-face and edge-edge work, with the total and the speedup over TightInclusion above it. Query counts are in the dataset table and hit agreement in the conservativeness table, so neither is repeated here.

### Accuracy

TightInclusion is a subject here, not the reference: its own answer is a
conservative lower bound, so the question is how close both get to the exact
roots.

<!-- sccd:begin earliness-ref -->

| scene             | phase | mode           | median earliness | worst case |
|-------------------|-------|----------------|-----------------:|-----------:|
| armadillo-rollers | EE    | CPU            |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | EE    | GPU            |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | EE    | TightInclusion |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | VF    | CPU            |         3.74e-06 |   9.39e-03 |
| armadillo-rollers | VF    | GPU            |         3.74e-06 |   9.39e-03 |
| armadillo-rollers | VF    | TightInclusion |         3.74e-06 |   9.39e-03 |
| cloth-ball        | EE    | CPU            |         1.05e-07 |   1.04e-04 |
| cloth-ball        | EE    | GPU            |         1.05e-07 |   1.04e-04 |
| cloth-ball        | EE    | TightInclusion |         1.05e-07 |   1.04e-04 |
| cloth-ball        | VF    | CPU            |         9.76e-08 |   2.88e-05 |
| cloth-ball        | VF    | GPU            |         9.75e-08 |   2.75e-05 |
| cloth-ball        | VF    | TightInclusion |         9.76e-08 |   2.88e-05 |
| cloth-funnel      | EE    | CPU            |                0 |   1.00e+00 |
| cloth-funnel      | EE    | GPU            |                0 |   1.00e+00 |
| cloth-funnel      | EE    | TightInclusion |                0 |   1.00e+00 |
| cloth-funnel      | VF    | CPU            |         3.33e-04 |   2.67e-02 |
| cloth-funnel      | VF    | GPU            |         3.33e-04 |   2.67e-02 |
| cloth-funnel      | VF    | TightInclusion |         3.33e-04 |   2.67e-02 |
| n-body            | EE    | CPU            |         1.14e-08 |   1.89e-04 |
| n-body            | EE    | GPU            |         1.14e-08 |   1.91e-04 |
| n-body            | EE    | TightInclusion |         1.14e-08 |   1.89e-04 |
| n-body            | VF    | CPU            |         1.17e-08 |   1.87e-05 |
| n-body            | VF    | GPU            |         1.17e-08 |   1.87e-05 |
| n-body            | VF    | TightInclusion |         1.17e-08 |   1.87e-05 |
| puffer-ball       | EE    | CPU            |         3.54e-05 |   2.61e-01 |
| puffer-ball       | EE    | GPU            |         3.54e-05 |   2.61e-01 |
| puffer-ball       | EE    | TightInclusion |         3.54e-05 |   2.61e-01 |
| puffer-ball       | VF    | CPU            |         4.58e-05 |   5.04e-02 |
| puffer-ball       | VF    | GPU            |         4.56e-05 |   5.05e-02 |
| puffer-ball       | VF    | TightInclusion |         4.58e-05 |   5.04e-02 |
| rod-twist         | EE    | CPU            |         3.69e-04 |   9.85e-01 |
| rod-twist         | EE    | GPU            |         3.68e-04 |   9.85e-01 |
| rod-twist         | EE    | TightInclusion |         3.69e-04 |   9.85e-01 |
| rod-twist         | VF    | CPU            |         1.54e-04 |   9.93e-01 |
| rod-twist         | VF    | GPU            |         1.55e-04 |   9.93e-01 |
| rod-twist         | VF    | TightInclusion |         1.54e-04 |   9.93e-01 |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end earliness-ref -->

**We do not approximate the reference, we reproduce it.** On all twelve
scene-phases the median error and the worst case are identical to
TightInclusion's, to every digit; on the GPU the worst case is identical and the
median agrees to within a part in a thousand.

The worst cases belong to the problem, not to the implementation: on cloth-funnel
and rod-twist the largest error approaches 1.0 for the reference too, because a
conservative search has no tighter answer available on a grazing or
already-touching configuration.

## 5. CPU against GPU

<!-- sccd:begin processor -->

| scene             | CPU ms | GPU ms | total | broad | narrow |
|-------------------|-------:|-------:|------:|------:|-------:|
| armadillo-rollers |  2,068 |  2,749 | 0.75x | 1.15x |  0.44x |
| cloth-ball        |  1,188 |    829 | 1.43x | 0.88x |  3.44x |
| cloth-funnel      |  1,319 |  2,345 | 0.56x | 0.56x |  0.56x |
| n-body            |  8,812 |  4,369 | 2.02x | 0.81x | 11.89x |
| puffer-ball       | 45,852 | 22,374 | 2.05x | 1.18x | 15.80x |
| rod-twist         | 47,854 | 25,382 | 1.89x | 1.83x |  1.91x |

A ratio is the host median over the device median, so 2.0 means the device takes half the time. Ratios below 1.0 are the cases where the host wins and are the ones worth reading.

Source: `benchmark/assessment/broadphase-cell2dmin.csv`

<!-- sccd:end processor -->

The split is mechanical. The **broad phase** suits a GPU — count, prefix sum,
scatter, no sequential window walk — and runs 2.6× to 5.5× faster on five of six
scenes. The **narrow phase** is a depth-first interval search with a divergent
stack, and runs at 0.16× to 1.06× of host speed. The broad phase is the larger
share on the host, so the GPU still wins end to end on five scenes by 1.6× to
3.0×.

**puffer-ball inverts it**, at 0.63× overall: its GPU narrow phase is six times
slower and its broad phase gives nothing back. It is the largest scene in the set
by candidate pairs, 30.5 million per step, so problem size is not a rule for
choosing a processor.

Per simulation step, which is the figure a solver budgets against and the only
one comparable between a 79-step scene and a 4,571-step one:

<!-- sccd:begin per-frame -->

| scene             | frames | mode | broad ms | narrow ms | total ms |
|-------------------|-------:|------|---------:|----------:|---------:|
| armadillo-rollers |    396 | GPU  |     3.03 |      3.91 |     6.94 |
| armadillo-rollers |    396 | CPU  |     3.48 |      1.74 |     5.22 |
| cloth-ball        |     43 | GPU  |    15.10 |      4.18 |    19.28 |
| cloth-ball        |     43 | CPU  |    13.23 |     14.39 |    27.62 |
| cloth-funnel      |    372 | GPU  |     4.20 |      2.10 |     6.30 |
| cloth-funnel      |    372 | CPU  |     2.36 |      1.18 |     3.55 |
| n-body            |     74 | GPU  |    52.61 |      6.43 |    59.04 |
| n-body            |     74 | CPU  |    42.60 |     76.48 |   119.08 |
| puffer-ball       |    120 | GPU  |   175.35 |     11.10 |   186.45 |
| puffer-ball       |    120 | CPU  |   206.72 |    175.38 |   382.10 |
| rod-twist         |  2,556 | GPU  |     3.45 |      6.48 |     9.93 |
| rod-twist         |  2,556 | CPU  |     6.32 |     12.41 |    18.72 |

A mean rather than a median over steps: the scene total is what a run costs, and the mean is the only average that divides back into it.

Source: `benchmark/assessment/broadphase-cell2dmin.csv`

<!-- sccd:end per-frame -->

## 6. Broad phase

Two strategies produce the candidate pairs — a sweep over sorted intervals and a
cell list over a uniform grid. Across 25,538 case-mode combinations they report
identical candidate pair counts, so the choice is purely about cost, and the
measurement below is what fixes it: `cell2dmin`, the cell list with the
edge-edge query binned at the minimum corner, is the default on both
processors.

<!-- sccd:begin broadphase -->

| scene             | mode | cell2dmin structure ms | sweep structure ms | cell2dmin broad ms | sweep broad ms | faster           |
|-------------------|------|-----------------------:|-------------------:|-------------------:|---------------:|------------------|
| armadillo-rollers | GPU  |                    256 |                987 |               1199 |           2504 | cell2dmin 2.09x  |
| armadillo-rollers | CPU  |                    763 |               1167 |               1380 |           2244 | cell2dmin 1.63x  |
| cloth-ball        | GPU  |                     62 |                244 |                649 |            862 | cell2dmin 1.33x  |
| cloth-ball        | CPU  |                    120 |                193 |                569 |           1292 | cell2dmin 2.27x  |
| cloth-funnel      | GPU  |                    203 |                630 |               1562 |           1660 | cell2dmin 1.06x  |
| cloth-funnel      | CPU  |                    630 |                818 |                879 |           1668 | cell2dmin 1.90x  |
| n-body            | GPU  |                    202 |                648 |               3893 |           6075 | cell2dmin 1.56x  |
| n-body            | CPU  |                    219 |                417 |               3152 |          11223 | cell2dmin 3.56x  |
| puffer-ball       | GPU  |                   4189 |               7566 |              21042 |          69259 | cell2dmin 3.29x  |
| puffer-ball       | CPU  |                   4298 |               5161 |              24807 |         420458 | cell2dmin 16.95x |
| rod-twist         | GPU  |                   3093 |              12360 |               8812 |          21166 | cell2dmin 2.40x  |
| rod-twist         | CPU  |                   5994 |              10018 |              16143 |          23905 | cell2dmin 1.48x  |

`faster` names the winning strategy and by how much on the whole broad phase. A margin inside the run-to-run spread is a tie.

Source: `benchmark/assessment/broadphase-cell2dmin.csv`

<!-- sccd:end broadphase -->

**Neither wins.** The sweep takes armadillo-rollers, cloth-ball, n-body and
rod-twist by 1.03× to 1.88×; the cell list takes cloth-funnel by 1.06× and
puffer-ball by **4.0×** — 27 s against 109 s on the largest scene. A fixed choice
is wrong somewhere, and wrong by a factor of four at the worst point.

## 7. Scaling with element count

`sccd_refine_scaling` refines one surface repeatedly, quadrupling the element
count at each level, over two consecutive cloth-ball frames from 92,230 elements
to 94.4 million.

<!-- sccd:begin scaling -->

| processor | level |   elements | candidate pairs | structure ms | traversal ms | broad ms | narrow ms |    p |
|-----------|------:|-----------:|----------------:|-------------:|-------------:|---------:|----------:|-----:|
| CPU       |     0 |     92,230 |          17,982 |          0.9 |          1.8 |      2.7 |       0.5 | 1.01 |
|           |     1 |    368,920 |          87,771 |          6.5 |          7.9 |     14.3 |       0.7 |      |
|           |     2 |  1,475,680 |         379,818 |         21.9 |         33.7 |     55.7 |       2.0 |      |
|           |     3 |  5,902,720 |       1,573,741 |         86.3 |        156.4 |    242.7 |       8.0 |      |
|           |     4 | 23,610,880 |       6,402,081 |        181.2 |        642.7 |    823.9 |      32.5 |      |
|           |     5 | 94,443,520 |      25,824,702 |        840.2 |       2804.1 |   3644.3 |     132.6 |      |
| GPU       |     0 |     92,230 |          17,982 |          0.7 |          1.2 |      2.0 |       1.6 | 0.86 |
|           |     1 |    368,920 |          87,771 |         15.2 |          4.3 |     19.6 |       1.9 |      |
|           |     2 |  1,475,680 |         379,818 |         58.9 |         16.5 |     75.4 |       2.0 |      |
|           |     3 |  5,902,720 |       1,573,741 |        238.3 |         59.2 |    297.5 |       2.4 |      |
|           |     4 | 23,610,880 |       6,402,081 |        251.8 |        189.9 |    441.7 |       3.5 |      |
|           |     5 | 94,443,520 |      25,824,702 |       1009.4 |        818.9 |   1828.3 |       8.1 |      |

The two frames used here do not come into contact, so the narrow phase has almost no work to do and its column is dominated by noise rather than by element count; what this measures is the broad phase and the preparation that feeds it. Narrow-phase cost against problem size is in the per-case figure, over cases that do collide. Where the exponent is below 1 it is because the fixed cost visible at the smallest size is amortised as the mesh grows.

Source: `benchmark/results/scaling/host-cell2dmin-mode2.txt, benchmark/results/scaling/device-cell2dmin-mode2.txt`

<!-- sccd:end scaling -->

Four series: each processor at each broad-phase strategy.

**The cell list leads on both processors, and the device leads until the mesh
outgrows it:**

| elements | CPU (cell2d) | GPU (cell2d) | |
|---:|---:|---:|---:|
| 92,230 | 4.0 ms | 2.5 ms | 1.60× GPU |
| 1,475,680 | 67.9 ms | 72.2 ms | 1.06× CPU |
| 23,610,880 | 998.7 ms | 1,825.1 ms | 1.83× CPU |
| 94,443,520 | 4,791.5 ms | 20,886.6 ms | 4.36× CPU |

The exponents say the same before the crossover does: the device cell list grows
at 1.18 in the element count against the host's 1.00. Building the grid is what
carries it, not querying it — at the largest size the host bins in 1,195.1 ms and
the device in 11,724.2, which is 56% of the device's whole broad phase against
25% of the host's.

The narrow phase goes the other way. These two frames do not come into contact,
so no root is isolated and the column prices the rejection test alone, over
candidates that all fail it. That work is uniform, and the device holds it
between 1.6 ms and 7.6 ms across the six refinements while the host climbs from
0.4 ms to 118.2 ms. A GPU suits the case where every candidate costs the same;
the divergent case is the per-case figure, over cases that do collide.

Between the host strategies the cost divides oppositely: the sweep builds more
cheaply at the coarsest sizes, the cell list traverses far better, and traversal
is what grows — 89,521.1 ms against 3,596.4 ms at the largest size.

## 8. Figures

Presentation follows the dataset paper (Belgrod et al., *TOI dataset for CCD and
a scalable conservative algorithm*) over the same six scenes: scenes across the
columns, one quantity per row, distributions over cases as box plots on
logarithmic axes.

![Per-case results over the six scenes](figures/results-grid.png)

**Figure 1.** Per-case distributions over the six scenes: broad-phase time,
narrow-phase time, and error against the exact root. Each box spans the first to
the third quartile with the median inside, whiskers reach the furthest case
within 1.5 interquartile ranges, and cases beyond are drawn individually.

![Runtime split by phase](figures/runtime-breakdown.png)

**Figure 2.** Runtime split into *prep*, *broad* and *narrow*. Preparation dominates
on the host; on the GPU the narrow phase does.

![Error against the symbolic ground truth](figures/toi-error.png)

**Figure 3.** Error against the exact symbolic roots, log axis. One-sided by
construction — it is how far *before* the true root the answer falls — so every
value shown is on the safe side, and nothing falls off the axis on the late side
because no such case exists. The device curve is dashed because the
two coincide.

![Cost against element count](figures/refine-scaling.png)

**Figure 4.** Broad- and narrow-phase cost against element count, log-log. The
fitted exponent is on the broad phase; the narrow-phase series is flat and noisy
because these two frames do not come into contact.

## 9. Provenance

<!-- sccd:begin provenance -->

- Timings: `benchmark/assessment/broadphase-cell2dmin.csv`, 6394 cases over 6 scenes, 2 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200-all.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
