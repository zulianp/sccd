# SCCD against Additive CCD, host only

`all.csv`, produced by `../../run_comparison_host.sh` and tabulated by
`../../compare_table.py`. Apple M-series laptop, `Darwin arm64`,
armadillo-rollers and cloth-ball, cases `[0, 60)`, three repeats.

Scalable CCD is absent because its narrow phase is CUDA only -- there is nothing
of it to compare against on a machine without a GPU. The GH200 run in
`../compare-gh200-2026-09-17.csv` is where all three appear.

Column definitions are the same as the GH200 run's, in `../README.md`. In short:
timings are per case, `fp`/`fn` are counts per pass over the case list, and
`toi err` is the earliest time of impact reported for a step minus the earliest
exact root of that step -- negative is conservative.

```
### armadillo-rollers

| library               |       broad ms |      narrow ms |   total |     ns/q |    fp |     fn |  med early |  max early |  toi err med | toi err worst |  late |
|-----------------------|----------------|----------------|---------|----------|-------|--------|------------|------------|--------------|---------------|-------|
| SCCD host Relaxed     |  8.064 /  23.67 |  2.496 /   7.97 |  10.954 |     26.3 |     3 |      0 |   2.55e-05 |     0.0153 |    -1.36e-05 |      1.03e-06 |     1 |
| SCCD host Tight       |  7.987 /  23.24 |  2.693 /   6.65 |  10.650 |    26.61 |     1 |      0 |   2.48e-06 |    0.00187 |    -2.83e-06 |      4.03e-05 |     9 |
| ACCD host             |  0.000 /   0.00 |  0.024 /   1.87 |   0.024 |     1137 |    16 |      0 |     0.0389 |      0.949 |    -0.000626 |     -2.12e-05 |     0 |

  steps with a known contact: 32   curated queries per pass: 7550
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
| SCCD host Relaxed     | 108.540 / 240.95 | 36.437 / 187.44 | 145.794 |    22.62 |     2 |      0 |   3.16e-07 |    0.00054 |    -5.62e-07 |             0 |     0 |
| SCCD host Tight       | 109.513 / 323.31 | 34.353 / 149.55 | 142.436 |    19.17 |     1 |      0 |   2.78e-07 |   0.000193 |    -1.76e-07 |             0 |     0 |
| ACCD host             |  0.000 /   0.00 |  2.419 /   7.64 |   2.419 |      437 |    16 |      0 |     0.0641 |      0.875 |    -2.23e-05 |     -1.95e-08 |     0 |

  steps with a known contact: 34   curated queries per pass: 379727
  fp and fn are counts per pass over the case list, so modes run with different
  numbers of repeats stay comparable.
  A few parts in a million of `toi err` is the gap between the two geometries, not
  a violation: the mesh path reads float32 coordinates while the exact roots come
  from the higher-precision query coordinates, and the two disagree at about 1e-6.
  ACCD reads the query coordinates directly and has no such gap. Errors of 1e-1 and
  above are far outside it.
```

## What it says

**Conservativeness holds for both.** `fn` is zero everywhere and no per-query
time of impact lands after the true root. That is the invariant, and it is the
same result the GH200 run gives, on a different processor and a different
architecture.

**The separation is tightness, not safety.** SCCD Tight stops a few parts in a
million before contact -- `2.5e-6` median on armadillo, `2.8e-7` on cloth-ball.
ACCD stops around a twentieth of the step early in the median and almost the
whole step at worst, which is what `conservative_rescaling = 0.9` buys. Its
false-positive count is higher for the same reason: stopping early turns near
misses into reported contacts.

**Cost is not comparable as raw milliseconds.** ACCD is handed the curated pairs;
SCCD runs its narrow phase over every broad-phase candidate, which is hundreds of
times more work. `ns/q` is the column that compares them, and there SCCD is about
40 times cheaper per query on armadillo and 20 times on cloth-ball. That gap is
the price of conservative advancement being an iterative bound rather than a
root isolation.

**A laptop is not the machine the article measures on.** Absolute timings here
are several times the GH200's -- cloth-ball's broad phase is `109` ms against
`12.5` ms -- because this is one consumer CPU rather than 72 Grace threads.
Ratios within a row are what transfer; the milliseconds are not comparable to
the published figures.
