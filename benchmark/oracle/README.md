# Narrow-phase accuracy oracle

`ti_oracle` compares each narrow-phase mode against TightInclusion, query by
query, on the CCD benchmark sets under `data/`. It is the acceptance gate for
the TI-accurate vectorized kernel: no kernel change lands without a run here.

## The invariant it enforces

SCCD must be **conservative**:

* a collision must never be missed (a false negative is a correctness failure);
* a reported time of impact must never be **later** than the true one, because a
  solver would then step straight through the contact;
* a false positive, or a time of impact that is too **early**, is acceptable --
  it only costs work.

The benchmark datasets ship exact roots (`data/<scene>/roots/<key>/toi.float64`,
NaN where there is no collision), so the invariant is checked directly:
`gt_missed` counts queries with a real root that a mode reported as no
collision, and `gt_late` counts times of impact after the true one. Those two
are the gate.

**Do not gate on TightInclusion.** Its answer is a *conservative lower bound* on
the true time of impact, not the truth, so a result landing between its answer
and the real root is inside the safe band and is not a violation. The `late`
column measures exactly that and over-reports; it is kept for context only. The
same caveat applies to `false_negative`: TightInclusion reports hits on queries
the exact roots say have no collision at all.

`ti_oracle` exits non-zero on `gt_missed` or `gt_late`, falling back to the
TightInclusion comparison only for datasets that ship no roots. Pass
`--no-strict` to report without failing.

## Build and run

```sh
cmake -S . -B build -DSCCD_ENABLE_TIGHT_INCLUSION=ON
cmake --build build -j --target ti_oracle
./build/ti_oracle data/cloth-funnel --csv benchmark/oracle/cloth-funnel.csv
```

Options: `--phase vf|ee|both`, `--max-files N`, `--tol T`, `--max-depth N`,
`--csv PATH`, `--violations-csv PATH`, `--no-strict`, `--gate MODE`, and with
CUDA `--device-float` and `--bench N` (see "Throughput"). `--float-geometry` narrows the input for every mode; see
"Why there is no float row to gate on".

`--gate ti-vec` restricts the exit code to the conservative kernel, which is
what a CI check should use: modes 0 and 1 are known to violate the invariant, so
the default `--gate all` always fails while they are still present.

It reads `data/<scene>/queries/*.csv` (exact rationals, 8 rows per query, the
same layout `benchmark/bench.exe.cpp` uses) and, when present,
`data/<scene>/mma_bool/<key>/mma_bool.uint8` for the Mathematica ground truth.
It does not need smesh.

## The device rows

Built with `-DSCCD_ENABLE_CUDA=ON`, the oracle adds two rows for the CUDA narrow
phase, run over the same queries and scored by the same gate:

* `device` -- the original kernel, whose acceptance test is the host's mode 0.
* `device-ti` -- the conservative kernel: TightInclusion's predicate (a *domain*
  width against a domain tolerance), its split rule (the axis furthest past its
  tolerance, halved at the midpoint), and its numerical error bound.

They are selected through `SCCD_NARROWPHASE_MODE` like the host modes, so the
device is no longer mode-blind. `--gate device-ti` scopes the exit code.

Exit codes are distinct on purpose: **2** means a mode broke conservativeness,
**3** means a CUDA call failed. CI needs to tell "the GPU is wrong" apart from
"the GPU is missing".

### Measured, on GH200 / CUDA 12.6

| dataset | row | hits vs TI | FP | gtMISS / gtLATE |
|---|---|---|---|---|
| armadillo-rollers | `device` | 4211 / 9588 (TI: 4200 / 9581) | 11 / 7 | 0 / 0 |
| armadillo-rollers | `device-ti` | **4200 / 9581** | **0 / 0** | 0 / 0 |
| cloth-ball (665k queries) | `device` | 107252 / 557669 | 0 / 1 | 0 / 0 |
| cloth-ball (665k queries) | `device-ti` | **107252 / 557668** | **0 / 0** | 0 / 0 |

`device-ti` reproduces TightInclusion's hit set exactly -- no false positives, no
false negatives, and a median relative error of `0.000e+00` -- on all 665k
cloth-ball queries as well as on armadillo-rollers, where the mode-0 kernel
carries 18 spurious hits and an absolute error two orders of magnitude larger.

It is also cheaper: no register spills, and 167 registers against the mode-0
kernel's 220, because it drops the Gauss-Newton `adaptive_split_longest_axis`.

### Throughput

`--bench N` merges every query file into one batch and times the device kernels
with CUDA events, uploads excluded. Use it rather than the per-file `ms` column
below, which is dominated by launch overhead: it made `device-ti` look 4.8x
*faster* than the mode-0 kernel on vertex-face, where batched it is 8% slower.

Best of 5 on GH200, double storage (float storage measured within 1%, as it
must be -- the search is double either way):

| dataset / phase | queries | `device` | `device-ti` |
|---|---|---|---|
| cloth-ball VF | 107k | 65.0 ms | 70.1 ms |
| cloth-ball EE | 558k | 294.9 ms | 473.7 ms |
| armadillo-rollers VF | 15k | 14.8 ms | 18.4 ms |
| **armadillo-rollers EE** | 39k | **26.1 ms** | **11509 ms** |

(armadillo-rollers EE measures 750-1050 ms on the 14.8k-query batch, with the
spread being run-to-run variance rather than a change -- see below.)

The first three are an ordinary 8-60% cost for exactness. The fourth is not: it
is 440x, and it is why **the device default has not been switched to the
conservative kernel**.

### CPU vs GPU, and `toi_stride` 0 vs 1

One Alps node: 72 Grace cores against one GH200. Best of 3, one merged batch,
CUDA events for the device rows and a host clock for the CPU rows, double
storage. `--bench` prints this and works with or without CUDA.

Milliseconds, stride 1 / stride 0:

| mode | | armadillo VF (5.8k) | armadillo EE (14.8k) | cloth-ball VF (107k) | cloth-ball EE (558k) |
|---|---|---|---|---|---|
| `scalar` | CPU | 17.1 / 0.15 | 27.9 / 0.24 | 96.7 / 1.46 | 461 / 7.8 |
| `vector` | CPU | **619 / 3.09** | 27.9 / 0.24 | **924 / 4.46** | 461 / 7.8 |
| `ti-vec` | CPU | 17.6 / 0.46 | 37.8 / 2.74 | 90.5 / 1.18 | 488 / 3.76 |
| `device` | GPU | 5.5 / 1.55 | 10.6 / 1.49 | 64.6 / 1.77 | 299 / 2.41 |
| `device-ti` | GPU | 5.9 / 2.37 | **1293 / 3.04** | 70.6 / 2.41 | 471 / 3.75 |

**Like-for-like, the GPU's advantage is small and not uniform.** Comparing the
same conservative search on both — `ti-vec` against `device-ti`, stride 1 — the
GPU is 3.0x faster on armadillo vertex-face, 1.3x on cloth-ball vertex-face,
level on cloth-ball edge-edge, and 34x *slower* on armadillo edge-edge. One
GH200 against 72 Grace cores is roughly a wash outside the pathological case.
Nobody should assume the device path is a win here without measuring the scene.

**At stride 0 the CPU is competitive or ahead.** `ti-vec` beats `device-ti` on
three of the four (0.46 vs 2.37, 1.18 vs 2.41, 3.76 vs 3.75) and ties on the
fourth. Stride 0 collapses the work so far that the GPU never gets to amortize a
kernel launch.

**`vector` (mode 1) is the outlier, and it is a host-side one.** Investigated:

* One real defect, fixed. `compute_codomain_neon` nested the eight-corner loop
  outside the lane loop, so the eight geometry values a corner needs -- which
  depend on the lane's query, not on the corner -- were gathered *inside* the
  corner loop: 192 two-element gathers per box for values that never changed. The
  running min/max was likewise round-tripped through memory once per corner.
  Both are hoisted; results are bit-identical (same hits, false positives and
  errors, verified against a pre-change build) and it is worth about 20%.
* Lane width is **not** the problem here, which is worth knowing because it is
  next door: 2, 4 and 8 lanes measure 301, 327 and 297 ms on armadillo
  vertex-face, all the same, where the same sweep moves the conservative kernel
  by 17%. `SCCD_VNARROWPHASE_VSIZE` now exists so this can be re-checked; 8 was
  hardcoded before.
* The remaining ~5x against `scalar` is not a micro-optimization. After the hoist
  the profile is flat across the whole search loop, and mode 1 finds more hits and
  more false positives than `scalar` on the same input -- it is exploring a
  different, larger tree. That is an algorithmic difference in its acceptance and
  split rules, not something to tune.

Given mode 1 is documented as "kept for comparison and for the benchmark, not
because it is safe to ship", the open question is whether it is worth keeping at
all: it is slower than the scalar reference it is meant to accelerate *and* than
the conservative kernel, which is also the only correct one.

Original observation: 619 ms against
`scalar`'s 17 ms on armadillo vertex-face, 924 against 97 on cloth-ball
vertex-face -- 36x and 10x *slower* than the scalar reference it is supposed to
accelerate. Note it matches `scalar` exactly on both edge-edge rows, which is the
tell: mode 1 has no edge-edge kernel and falls back to scalar there, so the
anomaly is confined to its vertex-face path. That kernel is not benchmarked
anywhere else and this is the first time it has been timed on a real batch.

### `toi_stride` 0 vs 1

The two strides are different questions answered by different kernels. Stride 1
writes a time of impact per query and gives each query a **block**; stride 0
writes a single global minimum and gives each query a **thread**, with every
query pruning against every other query's progress through one shared `toi`.

Best of 3, one batch, CUDA events, double storage, GH200:

| dataset / phase | queries | mode | stride 1 | stride 0 | ratio |
|---|---|---|---|---|---|
| cloth-ball VF | 107k | `device` | 64.7 ms | 1.77 ms | 37x |
| | | `device-ti` | 71.2 ms | 2.44 ms | 29x |
| cloth-ball EE | 558k | `device` | 294.3 ms | 2.50 ms | 118x |
| | | `device-ti` | 470.0 ms | 3.98 ms | 118x |
| armadillo VF | 5.8k | `device` | 5.5 ms | 1.56 ms | 3.5x |
| | | `device-ti` | 6.0 ms | 1.94 ms | 3.1x |
| armadillo EE | 14.8k | `device` | 9.2 ms | 1.49 ms | 6.2x |
| | | **`device-ti`** | **919.0 ms** | **3.01 ms** | **306x** |

Two things to take from it.

**The conservative kernel's pathological case is a stride-1 case.** The 100x-plus
gap between `device` and `device-ti` on armadillo-rollers edge-edge — the one that
has driven this whole investigation — is 919 ms against 9 ms at stride 1, and
3.0 ms against 1.5 ms at stride 0. At stride 0 the conservative kernel costs
roughly 2x the mode-0 one, which is an ordinary price for exactness. The blowup
belongs to the block-per-query kernel with per-query output, not to
TightInclusion's predicate.

**Do not read the ratio as a speedup.** Stride 0 computes strictly less: one
number instead of `n`. Most of its advantage is the shared upper bound — once any
query finds an early time of impact, every box in every other query whose `t`
lower bound sits above it is pruned on sight. That makes stride 0 fastest exactly
when an early collision exists, which is true of all these datasets. A frame with
no collision has nothing to prune against and would not show this. Choose the
stride by what the caller needs; this table is not an argument for switching.

### What the blowup actually is

`--bench` also times the host `ti-vec` row, which runs the *identical*
TightInclusion search. On the same 14.8k-query armadillo-rollers edge-edge batch:

| | ms |
|---|---|
| `ti-vec` (72 CPU cores, same search) | 37 |
| `device` (GPU, mode 0) | 9 |
| `device-ti` (GPU, same search as `ti-vec`) | 750 - 1050 |

A CPU beating the GPU by roughly 25x at its own algorithm is the finding: this is
not an expensive search, it is an inefficient implementation of it. Always time
the same algorithm on the other processor before concluding a workload is simply
hard.

**Read the spread in that last row before trusting any small comparison here.**
The same binary on the same input varies by about 40% between `srun`
allocations -- 749 ms and 1054 ms on consecutive runs. `--bench` takes a best-of-N
*within* one process, which does nothing about that. Any A/B below the ~40% level
has to be run inside a single allocation to mean anything.

Things that were tried, and what is actually established:

* **Seeding by one root box instead of the uniform dice** (removing the shape
  mismatch with TightInclusion's split rule): 871 ms against 803 ms. Inside the
  noise; nothing established either way.
* **A breadth-first bisection ramp** -- genuine nodes of TightInclusion's tree,
  and enough to fill the block: **much worse, and outside the noise.** It does not
  finish 60 files in eight minutes, and needs the shared stack raised from 1024 to
  1792 entries merely to complete 16. Its seeds sit deep and are the boxes that
  survived pruning, so each is an expensive subtree; they overflow the shared
  stack into the global one and the batch is rerun.
* **Growing the global stack geometrically instead of by the observed deficit:**
  ~950 ms against ~800 ms, so also inside the noise, and reverted for a reason
  that is not: every call memsets the entire stack capacity to restore the empty
  slot marker, so an oversized stack is paid for on every subsequent call rather
  than once.

What did ship, on the grounds that they are correct rather than measured wins:
the from-stack kernel drains many boxes per launch instead of one;
`SCCD_CHECK_CUDA`'s device synchronize is debug-only; and the grow-and-retry loop
is bounded by `SCCD_NP_MAX_RETRY_ROUNDS`. That last one matters most -- it was
unbounded, which is why a kernel producing too many boxes presented as a hang
rather than as a slow answer. Giving up there stays conservative: dropped boxes
can only make the time of impact too large, never too small.

**The cliff is the global-stack overflow path.** Block size, swept properly --
all variants in one `srun`, interleaved, so the 40% allocation variance cancels:

| threads/block | 32 | 64 | 128 | 256 |
|---|---|---|---|---|
| armadillo EE, 6 files | | | 187, 192 | **56, 65** |
| armadillo EE, 60 files | 1997 | 1118 | 875, 942, 1368 | **does not finish in 300 s** |

Read those two rows together, because the reversal is the finding. On a small
input 256 threads per block is **3x faster** than 128, and the trend from 32
upward is monotonic in the same direction. On a large input the same binary falls
off a cliff.

What changes between them is pressure on the global stack. `DfsSplit` dices into
`N` cells, so a 256-thread block starts 256 boxes per query against a shared
stack of 1024 entries -- four per thread. Small inputs never exhaust it and get
the full benefit of the extra parallelism; large ones overflow into the global
stack, and every overflow makes the host grow it and **rerun the entire batch
from its seeds**. That path is expensive enough to erase a 3x win and then some.

So the ordering of the earlier results was misleading, and so was my reading of
them. More parallelism per query genuinely helps -- until it tips the search into
the overflow path, which is where all the time goes. That also explains the
bisection ramp: its seeds are deep boxes that survived pruning, each an expensive
subtree, so it tipped straight into the same path and needed the shared stack
raised from 1024 to 1792 entries merely to complete 16 files.

**The next fix is the overflow path itself, not the search.** Three things make
it costly and all are addressable: a retry reruns the whole batch rather than
resuming, the stack grows by the last observed deficit rather than to a size that
would stop the retries, and every call memsets the entire stack capacity to
restore the empty-slot marker -- which is why growing geometrically to avoid the
retries made things worse rather than better. Fix those and the block size
becomes a free 3x.

**A note on the block-barrier hypothesis, now refuted.** The natural next suspect was
the `__syncthreads()` the DFS loop carries every iteration: 128 threads run at
the pace of the deepest subtree any of them holds, and at 220-240 registers only
2 blocks fit per SM. If that were the cost, a smaller block would help, since the
barrier degenerates toward a warp and more blocks fit. Measured properly -- all
four variants in one `srun`, interleaved, so the 40% allocation variance cancels:

| threads/block | 32 | 64 | 128 |
|---|---|---|---|
| ms | 1997 | 1118 | 942 |

Smaller blocks are monotonically worse, so the block-wide barrier and the low
block count are not what dominates -- the opposite of the prediction. Testing it
cost one line (a `DfsSplit<32>` specialization) instead of the ~150 lines of
warp-intrinsic code the hypothesis called for.

### Depth sensitivity

Holding the dataset fixed and varying `--max-depth`
on armadillo-rollers edge-edge:

| max-depth | 16 | 32 | 64 | 96 |
|---|---|---|---|---|
| ms | 4.7 | 52.8 | 989 | 940 |

TightInclusion's split subdivides until a box is inside the error box or within
its domain tolerance, and edge-edge tolerances on near-parallel edges are small
enough that boxes reach great depth. Deep trees overflow the per-block shared
stack into the global one, and every global-stack round costs a kernel launch and
a blocking readback. The fixes are already planned -- seeding by a bisection ramp
instead of a uniform dice, a thread-per-query triage pass, and warp-aggregated
stack operations -- and this is the case that should drive them.

Note also that `max_depth` here is a *per-box tree depth*, where TightInclusion's
own `max_itr` is a total pop budget. Depth 96 is far more permissive than
anything TightInclusion would run, which is safe but not free.

## Why there is no float row to gate on

The search computes in `double` whatever the interface scalar type is (see
`TC` in `src/cuda/sccd_narrowphase.cu` and `T_HP` in
`src/sccd_vnarrowphase_ti.hpp`). A caller may still *store* geometry as float;
`--device-float` uploads it that way, and `--float-geometry` narrows it for every
mode including the TightInclusion reference.

**Neither is a valid conservativeness test, and both are therefore non-gating.**
The exact roots shipped with the datasets belong to the original rational
geometry. Rounding the input to float moves the true root by around `1e-5`, so a
kernel that correctly brackets the root of the *rounded* geometry is reported
late against roots computed for the *unrounded* geometry. The runs say
`NOT GATED` and exit zero.

The control that settles it: the host `ti-vec` kernel is bit-identical to
TightInclusion and already computes in double. On armadillo-rollers it reports
**zero** late times of impact on the shipped geometry and **4543** under
`--float-geometry` (1438 vertex-face, 3105 edge-edge). The device's float numbers
match it almost exactly (1438 / 3054). That is the input rounding, not either
kernel.

An earlier revision of this file read those numbers as a defect in the CUDA
kernel and claimed single precision broke the invariant. That was the same
measurement error this file already warns about in another form: comparing
against a reference that does not correspond to the inputs actually used. To
judge a kernel, run it on the geometry the roots were computed for.

## Reading the columns

| column | meaning |
|---|---|
| `gtMISS!` | an exact root exists and the mode reported no collision. **Unsafe.** This is the gate. |
| `gtLATE!` | the mode's time of impact is after the true one. **Unsafe.** This is the gate. |
| `FN!` | TightInclusion found a collision, the mode did not. Context only -- TightInclusion over-reports. |
| `FP` | the mode found a collision, TightInclusion did not. Conservative: wasted work, not a failure. |
| `lateTI` | later than TightInclusion's answer but not necessarily later than the truth. Context only. |
| `relerr_*` | relative TOI error, over queries where TI's toi >= 1e-9. |
| `abserr_*` | absolute TOI error over all agreed hits. Carries the queries TI puts at toi == 0, where relative error is undefined. |
| `near_zero_ref` | how many agreed hits had TI toi < 1e-9. |

Relative error alone is misleading here: TI reports `toi == 0` often enough
(geometry already touching at the start of the step) that a naive mean is
dominated by division by zero. Both metrics are reported for that reason.

## Caveat on the `ms` column

Timings are per query *file*, and the benchmark files are small (a handful of
queries each). The vectorized mode batches 8 queries per SIMD block, so on these
files most lanes are idle and its timing is not representative of production
use. Treat `ms` as a coarse work proxy for the accuracy comparison only, and use
`benchmark/bench.exe.cpp` for real throughput numbers.

## Status

**Every host mode is conservative**, by construction and by measurement -- see
"Modes" below for why that is a property of the algorithm rather than of any one
kernel. Measured against the exact roots, modes 0,
1 and 2 report `gt_missed == 0` and `gt_late == 0` on every dataset. Mode 2
additionally reproduces TightInclusion bit for bit -- same hit set, times of
impact equal to the last bit -- for both vertex-face and edge-edge, which is what
makes it useful as a reference implementation rather than what makes it safe.

**On the device, the conservative kernel is the one to use.** See "The device
rows" above: `device-ti` reproduces TightInclusion's hit set exactly on every
dataset with exact roots, where the mode-0 kernel carries spurious hits and a
much larger time-of-impact error. Both pass the gate; `device-ti` is also the
cheaper of the two in registers.

An earlier revision of this file claimed modes 0 and 1 were not conservative.
That was a measurement error: it compared against TightInclusion rather than
against the true roots, and counted results inside the safe band as failures.

## Modes

Selected by `SCCD_NARROWPHASE_MODE` (or the legacy `SCCD_USE_VNARROW_PHASE`):

| value | kernel | differs by |
|---|---|---|
| 0 | scalar reference | accepts on several conditions beyond TightInclusion's, so it stops earlier and reports a less accurate (earlier) time of impact |
| 1 | fast vectorized (VF only; EE falls back to scalar) | same acceptance as 0, lane-packed |
| 2 | TightInclusion-equivalent vectorized, VF and EE | reproduces TightInclusion bit for bit; the most accurate, and the reference the others are measured against |

**They differ in accuracy and speed, not in safety.** An earlier revision of this
table marked 0 and 1 "not conservative", which was wrong twice over: it
contradicted the measured status directly below it, and it misread the algorithm.
Branch and bound with a sound inclusion function is conservative by
construction -- for the multilinear `F` here the hull of the eight corner values
is the *exact* range over a box, so accepting is safe however loose the test
(accept reports the box's `t` lower bound, so a looser test only reports earlier),
refinement is safe, and exhaustion accepts rather than drops.

The one thing that can lose a root is an **unsound rejection**: discarding a box
that contains one, which happens when the origin-containment test is padded by
less than the certified bound in `snumerical_error.hpp`. All three modes pad by
that bound. So a kernel measuring as non-conservative is a defect with three
places to look --

1. does the rejection test pad by the certified error bound, or by something that
   merely resembles it;
2. does accept report the box's `t` lower bound;
3. does depth, format-resolution or stack exhaustion accept rather than drop --

and not a property to document and live with. Two real instances have been found
by asking exactly that:

* **The CUDA kernel** substituted the user's distance tolerance for the bound.
  Larger than it in double, smaller in float.
* **The vertex-quad search** (`codomain_acceptance_vq`) padded with
  `std::numeric_limits<T>::epsilon()`, roughly 30x too small for unit-scale
  geometry -- while `sccd_get_numerical_error_vq_soa`, which computes the right
  bound, sat unused in the same file. Now wired in.

## Notes on the two guards

1. **The "no zero toi" guard is correct and is in place.** `srootfinder.hpp`
   (`accepted = accepted && (tt_min > 0)`, four sites), the vertex-quad finder,
   and three sites in `sccd_vnarrowphase.hpp` refuse to accept a box whose `t`
   lower bound is zero -- TightInclusion's `no_zero_toi` policy, which keeps a
   solver from stalling on a contact reported at exactly `t == 0`. It does not
   cost a collision: with it in place every mode still reports `gt_missed == 0`.
   TightInclusion does report hits at `t == 0` on some of these queries, but the
   exact roots agree with rejecting them.

2. **The edge-edge tolerance was genuinely wrong and is fixed.**
   `compute_edge_edge_tolerance_soa` returned the same value for `tol0` and
   `tol1`, where TightInclusion computes three distinct groupings of the corner
   differences; the old `tol1` and `tol2` were off by factors of up to 2.5x and
   6.9x. It is now written directly from TightInclusion's definition and agrees
   bitwise, pinned by `tolerance_ti_parity_test`.

## Baselines

The CSVs in this directory are the pre-change reference, recorded on the current
`fast` vectorized kernel. Regenerate them after any narrow-phase change and diff.
