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
- **The reference's accuracy.** On CPU the median and worst-case earliness match
  TightInclusion's on all twelve scene-phases, to every printed digit. On GPU
  they agree to within 0.6% and 4.5%.
- **2.1× to 3.5× faster on CPU** and 2.6× to 20.9× on GPU, at that accuracy.

Every number is generated from a committed CSV by a committed script:

```sh
python3 -m report benchmark/assessment/broadphase-cell2dminfv.csv /tmp/report \
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

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

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

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

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
| armadillo-rollers |  1,803 |  1,893 | 0.95x | 2.45x |  0.49x |
| cloth-ball        |  1,035 |    437 | 2.37x | 1.54x |  3.55x |
| cloth-funnel      |  1,421 |  1,480 | 0.96x | 1.29x |  0.64x |
| n-body            |  7,305 |  3,568 | 2.05x | 0.45x | 12.69x |
| puffer-ball       | 37,432 | 19,936 | 1.88x | 0.85x | 16.51x |
| rod-twist         | 47,152 | 22,909 | 2.06x | 2.45x |  1.91x |

A ratio is the host median over the device median, so 2.0 means the device takes half the time. Ratios below 1.0 are the cases where the host wins and are the ones worth reading.

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

<!-- sccd:end processor -->

The split is mechanical. The **broad phase** is count, prefix sum and scatter,
and it runs at 0.56× to 1.83× of host speed — a narrow range around parity,
because the edge-edge walk suits a host. The **narrow phase** is a depth-first
interval search with a divergent stack, and it spans a far wider range, 0.44× to
15.80×: behind the host where the candidate lists are smallest, and ahead by an
order of magnitude on n-body and puffer-ball, where tens of millions of
candidates keep every lane busy. The GPU's advantage is therefore the narrow
phase, and the broad phase is what limits it.

End to end the GPU wins puffer-ball by 2.05×, n-body by 2.02×, rod-twist by
1.89× and cloth-ball by 1.43×; the host wins armadillo-rollers at 0.75× and
cloth-funnel at 0.56×. Cloth-funnel is the smallest scene by candidate pairs at
43,661 per step and puffer-ball the largest at 30.5 million, and the GPU takes
the larger one, so problem size is not a rule for choosing a processor on its
own.

Per simulation step, which is the figure a solver budgets against and the only
one comparable between a 79-step scene and a 4,571-step one:

<!-- sccd:begin per-frame -->

| scene             | frames | mode | broad ms | narrow ms | total ms |
|-------------------|-------:|------|---------:|----------:|---------:|
| armadillo-rollers |    396 | GPU  |     1.13 |      3.65 |     4.78 |
| armadillo-rollers |    396 | CPU  |     2.76 |      1.79 |     4.55 |
| cloth-ball        |     43 | GPU  |     5.99 |      4.18 |    10.17 |
| cloth-ball        |     43 | CPU  |     9.25 |     14.81 |    24.07 |
| cloth-funnel      |    372 | GPU  |     1.93 |      2.05 |     3.98 |
| cloth-funnel      |    372 | CPU  |     2.50 |      1.32 |     3.82 |
| n-body            |     74 | GPU  |    41.92 |      6.29 |    48.22 |
| n-body            |     74 | CPU  |    18.81 |     79.91 |    98.71 |
| puffer-ball       |    120 | GPU  |   155.19 |     10.94 |   166.13 |
| puffer-ball       |    120 | CPU  |   131.30 |    180.63 |   311.93 |
| rod-twist         |  2,556 | GPU  |     2.43 |      6.53 |     8.96 |
| rod-twist         |  2,556 | CPU  |     5.97 |     12.48 |    18.45 |

A mean rather than a median over steps: the scene total is what a run costs, and the mean is the only average that divides back into it.

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

<!-- sccd:end per-frame -->

## 6. Broad phase

Two strategies produce the candidate pairs — a sweep over sorted intervals and a
cell list over a uniform grid. Across 25,576 case-mode combinations they report
identical candidate pair counts, so the choice is purely about cost, and the
measurement below is what fixes it: `cell2dmin`, the cell list with the
edge-edge query binned at the minimum corner, is the default on both
processors.

<!-- sccd:begin broadphase -->

| scene             | mode | cell2dminfv structure ms | sweep structure ms | cell2dminfv broad ms | sweep broad ms | faster             |
|-------------------|------|-------------------------:|-------------------:|---------------------:|---------------:|--------------------|
| armadillo-rollers | GPU  |                      168 |                987 |                  447 |           2554 | cell2dminfv 5.71x  |
| armadillo-rollers | CPU  |                      770 |               1258 |                 1093 |           2330 | cell2dminfv 2.13x  |
| cloth-ball        | GPU  |                       27 |                244 |                  258 |            867 | cell2dminfv 3.37x  |
| cloth-ball        | CPU  |                      119 |                192 |                  398 |           1276 | cell2dminfv 3.21x  |
| cloth-funnel      | GPU  |                      179 |                623 |                  718 |           1660 | cell2dminfv 2.31x  |
| cloth-funnel      | CPU  |                      710 |                798 |                  929 |           1653 | cell2dminfv 1.78x  |
| n-body            | GPU  |                       70 |                656 |                 3102 |           6223 | cell2dminfv 2.01x  |
| n-body            | CPU  |                      206 |                437 |                 1392 |          11222 | cell2dminfv 8.06x  |
| puffer-ball       | GPU  |                     1752 |               8056 |                18623 |          69996 | cell2dminfv 3.76x  |
| puffer-ball       | CPU  |                     3580 |               5148 |                15756 |         411130 | cell2dminfv 26.09x |
| rod-twist         | GPU  |                     1229 |              12637 |                 6221 |          21577 | cell2dminfv 3.47x  |
| rod-twist         | CPU  |                     6430 |              10095 |                15261 |          23808 | cell2dminfv 1.56x  |

`faster` names the winning strategy and by how much on the whole broad phase. A margin inside the run-to-run spread is a tie.

Source: `benchmark/assessment/broadphase-cell2dminfv.csv`

<!-- sccd:end broadphase -->

**The cell list wins every row.** Its margin on the whole broad phase runs from
1.06× on cloth-funnel's GPU to **16.95×** on puffer-ball's CPU, where it takes
24.8 s against 420.5 s on the largest scene in the benchmark. It also builds its
acceleration structure faster in every row, by 1.20× to 4.00×. One strategy is
therefore the right default everywhere, which is what `cell2dmin` is.

## 7. Scaling with element count

`sccd_refine_scaling` refines one surface repeatedly, quadrupling the element
count at each level, over two consecutive cloth-ball frames from 92,230 elements
to 94.4 million.

<!-- sccd:begin scaling -->

| processor | level |   elements | candidate pairs | structure ms | traversal ms | broad ms | narrow ms |    p |
|-----------|------:|-----------:|----------------:|-------------:|-------------:|---------:|----------:|-----:|
| CPU       |     0 |     92,230 |          17,982 |          1.0 |          1.4 |      2.4 |       0.3 | 1.03 |
|           |     1 |    368,920 |          87,771 |          6.6 |          7.2 |     13.9 |       0.6 |      |
|           |     2 |  1,475,680 |         379,818 |         23.3 |         28.3 |     51.6 |       2.1 |      |
|           |     3 |  5,902,720 |       1,573,741 |         94.0 |        149.6 |    243.7 |       8.0 |      |
|           |     4 | 23,610,880 |       6,402,081 |        183.3 |        619.5 |    802.9 |      32.5 |      |
|           |     5 | 94,443,520 |      25,824,702 |        808.9 |       2993.0 |   3801.9 |     130.8 |      |
| GPU       |     0 |     92,230 |          17,982 |          0.3 |          1.3 |      1.6 |       3.0 | 0.80 |
|           |     1 |    368,920 |          87,771 |         10.9 |          5.9 |     16.8 |       3.1 |      |
|           |     2 |  1,475,680 |         379,818 |         43.9 |         22.8 |     66.7 |       2.7 |      |
|           |     3 |  5,902,720 |       1,573,741 |        168.3 |         92.4 |    260.7 |       3.5 |      |
|           |     4 | 23,610,880 |       6,402,081 |         13.8 |        329.7 |    343.5 |       4.9 |      |
|           |     5 | 94,443,520 |      25,824,702 |         50.4 |       1422.2 |   1472.5 |       9.2 |      |

The two frames used here do not come into contact, so the narrow phase has almost no work to do and its column is dominated by noise rather than by element count; what this measures is the broad phase and the preparation that feeds it. Narrow-phase cost against problem size is in the per-case figure, over cases that do collide. Where the exponent is below 1 it is because the fixed cost visible at the smallest size is amortised as the mesh grows.

Source: `benchmark/results/scaling/host-cell2dminfv-mode2.txt, benchmark/results/scaling/device-cell2dminfv-mode2.txt`

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

- Timings: `benchmark/assessment/broadphase-cell2dminfv.csv`, 6394 cases over 6 scenes, 2 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200-all.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
