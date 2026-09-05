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

Against the dataset's exact roots. Both columns must be zero.

<!-- sccd:begin gate -->
<!-- sccd:end gate -->

## Timing

<!-- sccd:begin timing -->
<!-- sccd:end timing -->

<!-- sccd:begin comparison -->
<!-- sccd:end comparison -->

## Against TightInclusion

TightInclusion is the reference implementation of a certified conservative
narrow phase. Matching its hit count is the strongest agreement available;
reporting more hits is a false positive, which costs work and is never unsafe.
Its own answer is a lower bound on the truth rather than the truth, which is why
the conservativeness table above is measured against the exact roots instead.

<!-- sccd:begin reference -->
<!-- sccd:end reference -->

## Accuracy

How far before the true time of impact each mode reports. This is the axis the
two modes trade against speed, and early is the safe direction.

<!-- sccd:begin earliness -->
<!-- sccd:end earliness -->

## Provenance

<!-- sccd:begin provenance -->
<!-- sccd:end provenance -->
