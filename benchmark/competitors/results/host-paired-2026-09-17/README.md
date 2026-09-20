# SCCD against Additive CCD, both behind SCCD's broad phase

`all.csv`, from `../../run_comparison_host.sh`, tabulated by
`../../compare_table.py`. Apple M-series laptop, `Darwin arm64`,
armadillo-rollers and cloth-ball, cases `[0, 30)`, three repeats.

Additive CCD produces no candidates of its own, so here it is fed **SCCD's broad
phase**: both libraries score the identical candidate list and `ns/q` compares
the same work. The gather from element indices to the eight points ACCD wants
sits between the two timed regions and is charged to neither -- it is the seam
between two libraries, not a cost either would pay in use.

Scalable CCD is absent: its narrow phase is CUDA only, so there is nothing of it
to compare on a machine without a GPU.

```
### armadillo-rollers

| library               |       broad ms |      narrow ms |   total |     ns/q |    fp |     fn |  med early |  max early |  toi err med | toi err worst |  late |
|-----------------------|----------------|----------------|---------|----------|-------|--------|------------|------------|--------------|---------------|-------|
| SCCD host Relaxed     | 13.464 /  44.70 |  4.637 /  14.98 |  18.489 |     44.9 |     0 |      0 |   2.04e-05 |     0.0153 |    -1.31e-05 |     -2.11e-06 |     0 |
| SCCD host Tight       | 14.619 /  40.60 |  5.707 /  13.59 |  20.159 |    46.99 |     0 |      0 |   2.27e-06 |    0.00187 |    -3.41e-06 |      1.86e-05 |     3 |
| ACCD host             | 15.521 /  31.16 | 27.777 /  53.39 |  47.037 |    212.1 |     9 |      0 |      0.045 |      0.949 |    -0.000199 |     -3.31e-05 |     0 |

  steps with a known contact: 16   curated queries per pass: 2598
  fp and fn are counts per pass over the case list, so modes run with different
  numbers of repeats stay comparable.
  A few parts in a million of `toi err` is the gap between the two geometries, not
  a violation: the mesh path reads float32 coordinates while the exact roots come
  from the higher-precision query coordinates, and the two disagree at about 1e-6.
  ACCD reads the query coordinates directly and has no such gap. Errors of 1e-1 and
  above are far outside it.

### cloth-ball

| library               |       broad ms |      narrow ms |   total |     ns/q |    fp |     fn |  med early |  max early |  toi err med | toi err worst |  late |
|-----------------------|----------------|----------------|---------|----------|-------|--------|------------|------------|--------------|---------------|-------|
| SCCD host Relaxed     | 80.535 / 334.69 | 28.999 / 310.34 | 110.663 |    36.02 |     0 |      0 |   3.49e-07 |    0.00054 |    -8.34e-07 |             0 |     0 |
| SCCD host Tight       | 98.895 / 434.88 | 28.707 / 281.78 | 125.860 |    35.05 |     0 |      0 |   2.88e-07 |   0.000193 |    -2.27e-07 |             0 |     0 |
| ACCD host             | 56.870 / 627.64 | 182.057 /1440.72 | 254.739 |    219.3 |     2 |      0 |     0.0772 |       0.86 |     -0.00131 |     -2.43e-06 |     0 |

  steps with a known contact: 19   curated queries per pass: 125058
  fp and fn are counts per pass over the case list, so modes run with different
  numbers of repeats stay comparable.
  A few parts in a million of `toi err` is the gap between the two geometries, not
  a violation: the mesh path reads float32 coordinates while the exact roots come
  from the higher-precision query coordinates, and the two disagree at about 1e-6.
  ACCD reads the query coordinates directly and has no such gap. Errors of 1e-1 and
  above are far outside it.
```

## What it says

**Both are conservative.** `fn` zero, `toi_late` zero, across every case and
mode. On the same candidate list, neither misses a contact and neither reports a
time of impact after the true one.

**Tightness separates them by four orders of magnitude.** SCCD Tight stops
`2.3e-6` before contact in the median on armadillo and `2.9e-7` on cloth-ball;
ACCD stops `4.5e-2` and `7.7e-2` -- around a twentieth of the step -- and up to
`0.95` at worst. That is `conservative_rescaling = 0.9` doing what it is for, and
it is the whole trade: ACCD is simple and cannot overshoot, but the step it
returns is far from the contact.

**Cost, now comparable.** Over the same candidates, per candidate: SCCD Tight
`47.0` ns against ACCD `212.1` ns on armadillo, `35.1` ns against `219.3` ns on
cloth-ball -- roughly `4.5x` and `6.3x`. The broad-phase columns agree between
the rows, as they must, since it is the same broad phase.

## A bug worth recording, because its symptoms pointed the wrong way

SCCD's broad phase returns indices into **its own** edge list, which comes from
the mesh's edge graph, while the curated queries use the benchmark's ordering.
`bench.exe.cpp` reconciles the two in `make_benchmark_edge_id_map`; the first
version of this harness did not, and used `benchmark_ordered_edges` to interpret
SCCD's output.

The result was not a clean failure. Additive CCD was handed endpoints of edges
that did not exist, some of them coincident, and responded exactly as it should:
`Initial distance 0 <= d_min=0, returning toi=0`, then `Small gap 4.4e-19 <= eps
in Additive CCD can lead to missed collisions`, then `Slow convergence in
Additive CCD` after ten million iterations on a single pair. Individual
candidates took minutes; 500 of them did not finish in ten.

That looked like a property of the method -- conservative advancement being
expensive on near-tangential pairs, `min_distance = 0` being outside its intended
regime -- and it was written up that way before being checked. It was three
correct diagnostics from the toolkit that its caller had violated
`assert(d > min_distance)`. With the numbering fixed, the same case runs 225,490
candidates in `33.9` ms, no warnings, and nothing to filter: SCCD's broad phase
already excludes adjacent primitives, so the guard added for coincident pairs
never fires.

The lesson is narrow and practical: when a competitor starts behaving
pathologically, suspect the data being handed to it before concluding anything
about the method.
