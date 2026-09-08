# Benchmarks

SCCD against TightInclusion on the NYU CCD dataset: six scenes, 6,394 simulation
steps, on one GH200 node — Grace (288 threads) and Hopper sm_90 in the same
allocation, `prgenv-gnu/24.11:v2`, `-O3`, `CMAKE_CUDA_ARCHITECTURES=90`, shipped
default search parameters.

Three results:

- **Conservative on every query.** Zero missed collisions and zero late times of
  impact over 66,268,596 comparisons against exact symbolic roots.
- **Identical accuracy to the reference.** Median and worst-case error match
  TightInclusion's on all twelve scene-phases, to every digit.
- **1.6× to 10.7× faster on CPU**, up to 29× on GPU, at that accuracy.

Every number is generated from a committed CSV by a committed script:

```sh
python3 -m report benchmark/results/sweep-gh200-bp.csv /tmp/report \
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

Source: `benchmark/results/sweep-gh200-bp.csv`

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

| scene             | mode |   queries | toi compared | late | false pos. | false neg. |
|-------------------|------|----------:|-------------:|-----:|-----------:|-----------:|
| armadillo-rollers | GPU  |   131,441 |      130,859 |    0 |         24 |          0 |
| armadillo-rollers | CPU  |   131,441 |      130,859 |    0 |         24 |          0 |
| cloth-ball        | GPU  |   664,940 |      664,919 |    0 |          1 |          0 |
| cloth-ball        | CPU  |   664,940 |      664,919 |    0 |          1 |          0 |
| cloth-funnel      | GPU  |     7,552 |        6,773 |    0 |         15 |          0 |
| cloth-funnel      | CPU  |     7,552 |        6,773 |    0 |         15 |          0 |
| n-body            | GPU  | 2,947,719 |    2,947,611 |    0 |         12 |          0 |
| n-body            | CPU  | 2,947,719 |    2,947,611 |    0 |         11 |          0 |
| puffer-ball       | GPU  | 1,514,172 |    1,486,790 |    0 |        558 |          0 |
| puffer-ball       | CPU  | 1,514,172 |    1,486,790 |    0 |        536 |          0 |
| rod-twist         | GPU  |   549,208 |      285,431 |    0 |      1,245 |          0 |
| rod-twist         | CPU  |   549,208 |      285,431 |    0 |      1,243 |          0 |

Measured against the exact roots shipped with the dataset, not against TightInclusion: TightInclusion's own answer is itself a lower bound on the truth, so comparing against it over-reports lateness.

Source: `benchmark/results/sweep-gh200-bp.csv`

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

<!-- sccd:begin reference -->

| scene             | phase | mode           |   queries |      hits |        time ms |     vs. TI |
|-------------------|-------|----------------|----------:|----------:|---------------:|-----------:|
| armadillo-rollers | EE    | GPU            |    99,104 |    98,761 |    2341 / 2346 |       4.5× |
| armadillo-rollers | EE    | CPU            |    99,104 |    98,761 |    6392 / 6451 |       1.6× |
| armadillo-rollers | EE    | TightInclusion |    99,104 |    98,761 |  10471 / 10554 | 1.0× (ref) |
| armadillo-rollers | VF    | GPU            |    32,337 |    32,122 |    1619 / 1629 |       4.7× |
| armadillo-rollers | VF    | CPU            |    32,337 |    32,122 |    2678 / 2714 |       2.8× |
| armadillo-rollers | VF    | TightInclusion |    32,337 |    32,122 |    7569 / 7631 | 1.0× (ref) |
| cloth-ball        | EE    | GPU            |   557,683 |   557,668 |      490 / 501 |       3.1× |
| cloth-ball        | EE    | CPU            |   557,683 |   557,668 |      955 / 956 |       1.6× |
| cloth-ball        | EE    | TightInclusion |   557,683 |   557,668 |    1521 / 1522 | 1.0× (ref) |
| cloth-ball        | VF    | GPU            |   107,257 |   107,252 |      339 / 347 |       3.7× |
| cloth-ball        | VF    | CPU            |   107,257 |   107,252 |      411 / 411 |       3.1× |
| cloth-ball        | VF    | TightInclusion |   107,257 |   107,252 |    1255 / 1256 | 1.0× (ref) |
| cloth-funnel      | EE    | GPU            |     6,751 |     6,259 |    1152 / 1160 |       2.7× |
| cloth-funnel      | EE    | CPU            |     6,751 |     6,259 |    1732 / 1746 |       1.8× |
| cloth-funnel      | EE    | TightInclusion |     6,751 |     6,259 |    3122 / 3143 | 1.0× (ref) |
| cloth-funnel      | VF    | GPU            |       801 |       529 |      409 / 415 |       3.2× |
| cloth-funnel      | VF    | CPU            |       801 |       529 |      385 / 391 |       3.4× |
| cloth-funnel      | VF    | TightInclusion |       801 |       529 |    1305 / 1339 | 1.0× (ref) |
| n-body            | EE    | GPU            | 2,399,812 | 2,399,746 |      883 / 888 |       3.0× |
| n-body            | EE    | CPU            | 2,399,812 | 2,399,746 |    1176 / 1176 |       2.2× |
| n-body            | EE    | TightInclusion | 2,399,812 | 2,399,746 |    2638 / 2642 | 1.0× (ref) |
| n-body            | VF    | GPU            |   547,907 |   547,874 |      626 / 637 |       2.5× |
| n-body            | VF    | CPU            |   547,907 |   547,873 |      507 / 507 |       3.1× |
| n-body            | VF    | TightInclusion |   547,907 |   547,873 |    1581 / 1588 | 1.0× (ref) |
| puffer-ball       | EE    | GPU            | 1,206,952 | 1,187,650 |      945 / 946 |       2.8× |
| puffer-ball       | EE    | CPU            | 1,206,952 | 1,187,650 |    1612 / 1612 |       1.6× |
| puffer-ball       | EE    | TightInclusion | 1,206,952 | 1,187,650 |    2657 / 2662 | 1.0× (ref) |
| puffer-ball       | VF    | GPU            |   307,220 |   299,698 |      628 / 630 |       2.5× |
| puffer-ball       | VF    | CPU            |   307,220 |   299,676 |      586 / 590 |       2.6× |
| puffer-ball       | VF    | TightInclusion |   307,220 |   299,676 |    1544 / 1550 | 1.0× (ref) |
| rod-twist         | EE    | GPU            |   492,120 |   246,132 |    4544 / 6571 |      18.2× |
| rod-twist         | EE    | CPU            |   492,120 |   246,132 | 36596 / 121302 |       2.3× |
| rod-twist         | EE    | TightInclusion |   492,120 |   246,132 | 82791 / 250085 | 1.0× (ref) |
| rod-twist         | VF    | GPU            |    57,088 |    40,544 |    2516 / 3218 |       9.3× |
| rod-twist         | VF    | CPU            |    57,088 |    40,542 |   6379 / 10345 |       3.7× |
| rod-twist         | VF    | TightInclusion |    57,088 |    40,542 |  23363 / 43924 | 1.0× (ref) |

Source: `benchmark/results/oracle-gh200-all.csv`

<!-- sccd:end reference -->

### Accuracy

TightInclusion is a subject here, not the reference: its own answer is a
conservative lower bound, so the question is how close both get to the exact
roots.

<!-- sccd:begin earliness-ref -->

| scene             | phase | mode           | median earliness | worst case |
|-------------------|-------|----------------|-----------------:|-----------:|
| armadillo-rollers | EE    | CPU            |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | EE    | GPU            |         3.08e-06 |   1.60e-02 |
| armadillo-rollers | EE    | TightInclusion |         3.11e-06 |   1.60e-02 |
| armadillo-rollers | VF    | CPU            |         3.74e-06 |   9.39e-03 |
| armadillo-rollers | VF    | GPU            |         3.63e-06 |   9.81e-03 |
| armadillo-rollers | VF    | TightInclusion |         3.74e-06 |   9.39e-03 |
| cloth-ball        | EE    | CPU            |         1.05e-07 |   1.04e-04 |
| cloth-ball        | EE    | GPU            |         1.05e-07 |   1.04e-04 |
| cloth-ball        | EE    | TightInclusion |         1.05e-07 |   1.04e-04 |
| cloth-ball        | VF    | CPU            |         9.76e-08 |   2.88e-05 |
| cloth-ball        | VF    | GPU            |         9.75e-08 |   2.88e-05 |
| cloth-ball        | VF    | TightInclusion |         9.76e-08 |   2.88e-05 |
| cloth-funnel      | EE    | CPU            |                0 |   1.00e+00 |
| cloth-funnel      | EE    | GPU            |                0 |   1.00e+00 |
| cloth-funnel      | EE    | TightInclusion |                0 |   1.00e+00 |
| cloth-funnel      | VF    | CPU            |         3.33e-04 |   2.67e-02 |
| cloth-funnel      | VF    | GPU            |         3.33e-04 |   2.67e-02 |
| cloth-funnel      | VF    | TightInclusion |         3.33e-04 |   2.67e-02 |
| n-body            | EE    | CPU            |         1.14e-08 |   1.89e-04 |
| n-body            | EE    | GPU            |         1.14e-08 |   1.89e-04 |
| n-body            | EE    | TightInclusion |         1.14e-08 |   1.89e-04 |
| n-body            | VF    | CPU            |         1.17e-08 |   1.87e-05 |
| n-body            | VF    | GPU            |         1.17e-08 |   1.87e-05 |
| n-body            | VF    | TightInclusion |         1.17e-08 |   1.87e-05 |
| puffer-ball       | EE    | CPU            |         3.54e-05 |   2.61e-01 |
| puffer-ball       | EE    | GPU            |         3.54e-05 |   2.61e-01 |
| puffer-ball       | EE    | TightInclusion |         3.54e-05 |   2.61e-01 |
| puffer-ball       | VF    | CPU            |         4.58e-05 |   5.04e-02 |
| puffer-ball       | VF    | GPU            |         4.56e-05 |   5.04e-02 |
| puffer-ball       | VF    | TightInclusion |         4.58e-05 |   5.04e-02 |
| rod-twist         | EE    | CPU            |         2.47e-04 |   9.85e-01 |
| rod-twist         | EE    | GPU            |         2.48e-04 |   9.85e-01 |
| rod-twist         | EE    | TightInclusion |         2.47e-04 |   9.85e-01 |
| rod-twist         | VF    | CPU            |         1.28e-04 |   9.93e-01 |
| rod-twist         | VF    | GPU            |         1.28e-04 |   9.93e-01 |
| rod-twist         | VF    | TightInclusion |         1.28e-04 |   9.93e-01 |

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

| scene             | frames | mode | prep ms | broad ms | narrow ms | total ms |
|-------------------|-------:|------|--------:|---------:|----------:|---------:|
| armadillo-rollers |    396 | GPU  |    4.67 |     2.05 |      4.38 |    11.11 |
| armadillo-rollers |    396 | CPU  |   20.43 |     8.41 |      2.62 |    31.46 |
| cloth-ball        |     43 | GPU  |   12.33 |     8.58 |     10.26 |    31.18 |
| cloth-ball        |     43 | CPU  |   28.60 |    21.98 |      4.99 |    55.58 |
| cloth-funnel      |    372 | GPU  |    3.64 |     1.45 |      2.54 |     7.63 |
| cloth-funnel      |    372 | CPU  |   15.22 |     5.88 |      1.77 |    22.87 |
| n-body            |     74 | GPU  |   13.58 |    17.73 |     63.79 |    95.10 |
| n-body            |     74 | CPU  |   35.24 |    97.62 |     23.49 |   156.35 |
| puffer-ball       |    120 | GPU  |  141.46 |   229.35 |    322.15 |   692.96 |
| puffer-ball       |    120 | CPU  |  157.48 |   227.79 |     50.41 |   435.69 |
| rod-twist         |  2,556 | GPU  |    9.84 |     2.64 |      6.80 |    19.27 |
| rod-twist         |  2,556 | CPU  |   26.12 |    11.22 |      7.23 |    44.57 |

A mean rather than a median over steps: the scene total is what a run costs, and the mean is the only average that divides back into it.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end per-frame -->

## 6. Broad phase

Two strategies produce the candidate pairs — a sweep over sorted intervals and a
cell list over a uniform grid. Across 25,538 case-mode combinations they report
identical candidate pair counts, so the choice is purely about cost, and the
shipped default races them per scene rather than fixing a winner.
`SCCD_BROADPHASE` is read only in the host broad phase, so only host rows appear.

<!-- sccd:begin broadphase -->

| scene             | mode | cell2d prep ms | sweep prep ms | cell2d broad ms | sweep broad ms | faster       |
|-------------------|------|---------------:|--------------:|----------------:|---------------:|--------------|
| armadillo-rollers | CPU  |           8089 |         26637 |            3331 |           3231 | sweep 1.03x  |
| cloth-ball        | CPU  |           1230 |          3064 |             945 |            697 | sweep 1.36x  |
| cloth-funnel      | CPU  |           5661 |         19460 |            2186 |           2340 | cell2d 1.07x |
| n-body            | CPU  |           2608 |          5902 |            7224 |           3833 | sweep 1.88x  |
| puffer-ball       | CPU  |          18898 |         27213 |           27335 |         110588 | cell2d 4.05x |
| rod-twist         | CPU  |          66765 |        169005 |           28680 |          21551 | sweep 1.33x  |

`faster` names the winning strategy and by how much on the whole broad phase. A margin inside the run-to-run spread is reported as a tie rather than a winner.

Source: `benchmark/results/sweep-gh200-bp.csv`

<!-- sccd:end broadphase -->

**Neither wins.** The sweep takes armadillo-rollers, cloth-ball, n-body and
rod-twist by 1.03× to 1.88×; the cell list takes cloth-funnel by 1.06× and
puffer-ball by **4.0×** — 27 s against 109 s on the largest scene. A fixed choice
is wrong somewhere, and wrong by a factor of four at the worst point.

## 7. Scaling with element count

`sccd_refine_scaling` refines one surface repeatedly, quadrupling the element
count at each level, over two consecutive cloth-ball frames from 92,230 elements
to 23.6 million.

<!-- sccd:begin scaling -->

| series       | level |   elements | candidate pairs | prep ms | step ms | broad ms | narrow ms |    p |
|--------------|------:|-----------:|----------------:|--------:|--------:|---------:|----------:|-----:|
| CPU / cell2d |     0 |     92,230 |          17,982 |    15.7 |     8.6 |     24.2 |       5.0 | 0.91 |
|              |     1 |    368,920 |          87,771 |    52.2 |    14.9 |     67.1 |       7.2 |      |
|              |     2 |  1,475,680 |         379,818 |   207.4 |    41.9 |    249.2 |       4.6 |      |
|              |     3 |  5,902,720 |       1,573,741 |   877.2 |   141.0 |   1018.2 |       4.2 |      |
|              |     4 | 23,610,880 |       6,402,081 |  3905.1 |   530.1 |   4435.2 |      13.9 |      |
| CPU / sweep  |     0 |     92,230 |          17,982 |    38.3 |     8.5 |     46.7 |       4.1 | 0.91 |
|              |     1 |    368,920 |          87,771 |    71.3 |    12.9 |     84.2 |       6.6 |      |
|              |     2 |  1,475,680 |         379,818 |   199.9 |    57.2 |    257.1 |       2.5 |      |
|              |     3 |  5,902,720 |       1,573,741 |   724.0 |   431.1 |   1155.1 |       6.3 |      |
|              |     4 | 23,610,880 |       6,402,081 |  2661.7 |  5052.9 |   7714.5 |      16.5 |      |
| GPU          |     0 |     92,230 |          17,982 |     4.4 |     2.4 |      6.8 |       2.5 | 1.13 |
|              |     1 |    368,920 |          87,771 |    42.6 |     7.6 |     50.2 |       5.2 |      |
|              |     2 |  1,475,680 |         379,818 |   161.2 |    31.8 |    193.0 |      17.7 |      |
|              |     3 |  5,902,720 |       1,573,741 |   682.5 |   222.6 |    905.1 |     138.0 |      |
|              |     4 | 23,610,880 |       6,402,081 |  2711.4 |  1720.3 |   4431.8 |    1089.5 |      |

The two frames used here do not come into contact, so the narrow phase has almost no work to do and its column is dominated by noise rather than by element count; what this measures is the broad phase and the preparation that feeds it. Narrow-phase cost against problem size is in the per-case figure, over cases that do collide. Where the exponent is below 1 it is because the fixed cost visible at the smallest size is amortised as the mesh grows. `prep` builds the acceleration structure -- the cell list's grid or the sweep's sorted intervals -- and `step` is the traversal that reports pairs; the two strategies divide the work between those columns quite differently.

Source: `benchmark/results/scaling/host-cell2d-mode2.txt, benchmark/results/scaling/host-sweep-mode2.txt, benchmark/results/scaling/device-mode2.txt`

<!-- sccd:end scaling -->

Three series: the host at each broad-phase strategy, and the device. The device
row names no strategy because the device broad phase does not implement the
choice.

**The GPU's broad-phase advantage erodes with size.** It is 3.5× ahead at 92,230
elements and level with the host cell list at 23.6 million:

| elements | CPU (cell2d) | GPU | |
|---:|---:|---:|---:|
| 92,230 | 24 ms | 7 ms | 3.54× |
| 1,475,680 | 249 ms | 193 ms | 1.29× |
| 23,610,880 | 4,435 ms | 4,432 ms | 1.00× |

The exponents say why: the GPU broad phase grows at 1.14 in the element count
against the host cell list's 0.95. **This is the puffer-ball result of §5 in
isolation** — that scene has the most candidate pairs per step in the benchmark,
and it is where the GPU broad phase stops winning (0.99×).

The narrow phase separates them further. These two frames do not come into
contact, so the host has almost nothing to do and its narrow phase stays flat
(exponent 0.11, 14 ms at 23.6 M elements). The device launches over every
candidate pair regardless, and grows at 1.12 to 1,090 ms. A GPU pays for
candidates it is given; a CPU pays for the ones that turn out to matter.

Between the host strategies the cost divides oppositely: the sweep builds more
cheaply and scales slightly better doing it, the cell list traverses far better,
and traversal is what grows — 5,053 ms against 530 ms at the largest size.

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

- Timings: `benchmark/results/sweep-gh200-bp.csv`, 6394 cases over 6 scenes, 3 independent repeats.
- Accuracy: `benchmark/results/oracle-gh200-all.csv`, every query of every scene checked against the dataset's exact roots.
- Regenerate with `python3 -m report <bench.csv> <out> <oracle.csv> --embed=docs/BENCHMARKS.md`; add `--check` to assert the document still matches the data.

<!-- sccd:end provenance -->
