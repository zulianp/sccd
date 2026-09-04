# CUDA narrow phase: where it stands, and what to try next

Branch `cuda-op-np`. This is a plan for exploring variants, not a commitment to
any of them. Every number below is measured and sourced; the point of writing
them down is that several plausible-sounding fixes have already been tried and
refuted, and the refutations are as valuable as the open questions.

## Where it stands

After the queue rewrite, mode 2 on GH200 against 72 Grace cores, narrow-phase
milliseconds, 16 cases:

| scene | host `Tight` | device `Tight` | |
|---|---:|---:|---|
| cloth-funnel | 7.4 | 37.4 | 5.1× host |
| armadillo-rollers | 19.1 | 64.3 | 3.4× host |
| cloth-ball | 115.8 | **106.7** | 1.09× device |

The time gap is a **work** gap, not a per-box cost gap. Counted with
`-DSCCD_NP_COUNT_BOXES` on cloth-funnel, split by code path:

| path | queries | host boxes | device boxes | ratio |
|---|---:|---:|---:|---:|
| `find_earliest_impact_time` (`toi_stride=0`) | 881,420 | 1,038,395 | 97,757,836 | **94×** |
| `find_impact_times` (`toi_stride=1`) | **274** | 505,243 | 705,986,614 | **1397×** |

111 boxes per query against the host's 1.2 on the interface everyone calls; 2.58
million against 1,844 on a path that handled 274 queries. Two different problems.

The shape of the 94× is the finding. By DFS level, zero-stride path:

| level | host | device |
|---:|---:|---:|
| 0 | 531,777 | 288,962 |
| 10 | 14,334 | 3,697,397 |
| 16 | 6,992 | **4,195,073** |
| 25 | 11,727 | 1,228,312 |
| 40 | 23,552 | 38,286 |

The host decays from level 0 and carries a thin tail. The device *broadens*,
peaking three orders of magnitude above the host at level 16. In branching terms
the host keeps ~0.7 children per box through the middle depths and the device
~1.8: the host discards most of what it produces, the device keeps nearly all of
it. Same predicate, same split rule, same tolerances — all transcribed.

Per query, the two agree on easy work and diverge with difficulty. Below 8 boxes
the distributions match within a few percent; at ≥1,024 boxes the device has 71×
as many queries, and it has 208 queries costing over a million boxes where the
host's worst in the whole scene costs 48,561.

## Ruled out — do not retry these

Each was measured, and the measurement is in `wip/ASSESSMENT.md`.

| Hypothesis | Result |
|---|---|
| Depth cap truncating or forcing accepts | Zero accepts at the cap on either side, out of 2.3M host and 50M device boxes. Both reach depth ~60. |
| The split rule | `bisect_ti_axis` picks the largest width/tolerance ratio, the same rule as the host's `tight_split_axis`. |
| The seeding dice | Replacing the 128-way dice with a single root seed moved armadillo edge-edge 803 → 871 ms. A breadth-first ramp was worse still. |
| Contention on the shared `toi[0]` | Refreshing less often is monotonically worse on every scene: 638 → 3,377 ms on cloth-funnel at every-256. |
| Sequential batching to recover the host's collapsing bound | 12.7× better on cloth-funnel, 26× and 16.5× *worse* on the other two. And the cloth-funnel gain is not pruning: box count fell 7% against a 12.7× time change, so whatever earns it is not the mechanism. |
| Corner inheritance, as a drop-in | The host carries `corner[3][8]` and evaluates 12 corners per split against the device's 48. It cannot be built in this kernel's shape: 24 doubles is 48 more registers into a budget with 17 left (238 of 255 used, no spills). |
| Occupancy as the differentiator | Both modes sit at ~12.5%, two blocks per SM, register-bound. Shrinking shared memory 57,376 → 3,584 bytes added no blocks. |

**The cause of the 94× is not established.** Everything above says what it is
not. That is why step 0 exists.

## Step 0, first result: the host's 1.2 boxes per query is not its search

Measured with `sccd_np_trace`, cloth-funnel `227vf`, 92 queries, the **same host
`Tight` kernel** on the **same data**, differing only in how the queries are
presented to it:

| | host boxes | per query |
|---|---:|---:|
| one query per call, `max_toi = 1` | 51,509 | **559.88** |
| all 92 in one `toi_stride=0` call | 4,098 | **44.54** |

**12.6× from the shared collapsing bound alone**, with nothing else changed. The
batched call reaches an earliest impact of 0.045, and once the bound is there
almost every remaining query is rejected at its root box.

That is the same mechanism, at 92 queries, that produces 1.2 boxes per query over
881,420. The host is not a better searcher than the device — it is the same
search running behind a bound that has already collapsed. Its advantage is
sequencing: `parallel_for_br_dynamic` hands out chunks of `VSIZE` queries, each
chunk seeds `toi_q` from the global minimum as it stands when that chunk starts,
and a few early chunks establish a bound the remaining hundreds of thousands are
killed by. The device has no such sequence: one launch starts every query at once
against `max_toi`, and its bound can only tighten as the launch runs.

**So the 94× is very likely a bound-collapse-rate difference, not a search-quality
difference.** That is a different problem with a different set of fixes, and it
promotes B1 and C1 below for a reason they did not previously have.

Two supporting measurements, and one caveat:

- **Removing the host's per-chunk seeding** (`SCCD_NP_NO_GLOBAL_SEED=1`, an
  instrumented-build-only knob) over four cases of cloth-funnel moves the
  aggregate only 0.72 → 0.80 boxes per query. But per call it ranges from 0.78×
  to **26.35×**: the aggregate is diluted because `max_toi` already arrives
  collapsed from the previous call, so in most calls the knob has nothing to
  remove. The knob isolates within-call seeding only.
- **`tight_toi >= relaxed_toi` held on all 92 queries**, with zero violations.
  That ordering is what B1's culling argument needs, and it is now measured rather
  than assumed.
## Step 0, second result: in isolation the device is 2.3x the host, not 94x

`sccd_np_trace --device` on a GH200, the same 92 queries, one query per call,
`max_toi = 1`, `toi_stride = 1` — nothing another query found available to prune
either side:

| | boxes | per query |
|---|---:|---:|
| host `Tight` | 51,509 | 559.88 |
| device `Tight` | 118,808 | **1,291.39** |
| | | **2.31×** |

Worst single query 5.0×; the ratio is 3–5× across the costliest few and lower
elsewhere. Not 94×, and not 1397×.

**Step 0 is answered.** The device's per-query search is about twice the host's.
Everything above that in the scene-level numbers is the bound: the host's
sequencing collapses it to the scene minimum within a few chunks and the device
never gets that head start.

That settles the priorities. The per-box work — corner inheritance at 12
evaluations per split against 48, family A — is a 4× arithmetic factor sitting
behind a 2.3× box factor, on a problem whose remaining 40× is somewhere else
entirely. **The lever is the bound collapse rate, and B1 and C1 are the two ways
to pull it.**

The 2.3× residual is still worth one careful look, and it is a small, contained
question rather than an architectural one. The first thing to check is the
counting convention: the device ticks once per `evaluate_cell_3d_policy` call,
which happens for *both* children of every split including ones discarded
immediately, while the host ticks once per box taken off its work list. Some of
2.3× may be that and not real work.

## What is still worth tracing box by box

Only if the device's isolated per-query cost turns out **not** to match the
host's. In that case: *trace a single hard query box by box, not another
aggregate.*

Pick one query from the ≥1,048,576-box set on cloth-funnel. Run it through the
host `Tight` kernel and the device kernel with a per-box trace — `(level, t, u, v)
bounds, the three acceptance flags, the rejection decision, the value of the
bound it was tested against` — and diff the two trees.

Three outcomes, each pointing somewhere different:

- **The trees diverge at a specific box.** Then the predicates are not as
  transcribed as they look, and this is a correctness-adjacent bug in one of
  them, not an architecture problem. Cheapest possible outcome.
- **The trees are identical in shape but the device tests boxes against a looser
  `t` bound.** Then the problem is when the bound arrives, and the fix is
  scheduling — families C1 and C3 below.
- **The trees genuinely differ because the host's DFS reaches a leaf sooner.**
  Then it is search order, and the fix is C1.

This is a day's work and it decides which of the families below is worth
building. Nothing else in this document should start before it finishes.

## Family A — reduce the cost of a box

These do not shrink the search. They are worth pursuing only alongside a C.

**A1. One query per warp instead of one thread per box.** Today every thread
holds its own query's geometry in registers: six `Vec4`s is 48 registers before
`atol`, `aerr` and the box. A warp cooperating on one query's 24 corner
evaluations holds the geometry once, which is what makes A2 possible and lifts
occupancy off its 12.5% floor. This is the largest single change in this
document and it changes the kernel's shape, not a parameter.

**A2. Corner inheritance.** The host's 4× arithmetic advantage per split, which
A1's register budget makes reachable. Not worth doing alone: it is 4× on a 94×
gap.

**A3. Geometry in shared memory rather than registers.** A cheaper approximation
of A1 — the block holds the geometry for the queries it is working on. Bounded
by how many distinct queries a block's stack mixes, which is currently
unbounded and would have to be capped.

## Family B — reduce the number of queries that need the tight search

**B1. `Relaxed` prepass to cull queries.** *Highest-value item, and the estimate
below is why — but note carefully which machine it has to be estimated on.* `Relaxed` averages 3.6 boxes per query
against `Tight`'s 944, and its answer is *at or before* the true time of impact.
So for the earliest-impact query: if a query's `Relaxed` time of impact is later
than the best `Tight` answer found so far, its true time of impact is later too,
and the query cannot hold the earliest impact — it can be discarded outright.
That is sound, exploits the invariant's asymmetry rather than fighting it, and
costs one cheap pass over the candidate set.

### The cull rate, measured

How many queries have a `Relaxed` time of impact below the scene's `Tight`
earliest impact — i.e. how many survive the cull. `sccd_np_trace`, vertex-face,
the largest query file in each scene:

| scene | file | queries | survivors | cull |
|---|---|---:|---:|---:|
| cloth-funnel | 227vf | 92 | 1 | 92× |
| cloth-funnel | 317vf | 24 | 1 | 24× |
| armadillo-rollers | 326vf | 4,652 | 1 | **4,652×** |
| armadillo-rollers | 229vf | 1,597 | 1 | 1,597× |
| cloth-ball | 92vf | 19,034 | 1 | **19,034×** |
| cloth-ball | 91vf | 11,278 | 0 | 11,278× |

One query, or none, survives. The earliest impact in a call is extremely
concentrated, which is exactly the property B1 needs.

### The cost side, and the trap in measuring it on the host

The prepass has to pay for itself. Batched — the way the library runs — the two
modes cost almost the same on the host:

| scene | file | queries | `Relaxed` b/q | `Tight` b/q | |
|---|---|---:|---:|---:|---:|
| cloth-funnel | 227vf | 91 | 6.53 | 44.54 | 6.8× |
| cloth-funnel | 317vf | 23 | 77.46 | 283.12 | 3.7× |
| armadillo-rollers | 326vf | 4,651 | 13.14 | 16.04 | 1.22× |
| armadillo-rollers | 229vf | 1,596 | 28.59 | 38.57 | 1.35× |
| cloth-ball | 92vf | 19,033 | 7.75 | 7.85 | **1.01×** |
| cloth-ball | 91vf | 11,277 | 8.03 | 8.47 | 1.05× |

On the host a `Relaxed` prepass would roughly **double** the work: it costs what
`Tight` costs, and `Tight` still has to run afterwards. The reason is the same
mechanism as everything else here — once the shared bound has collapsed both
modes reject at the root, and a root rejection costs the same in either mode.
`Tight`'s extra cost only appears on the queries that survive the bound, and
there is one of those.

**That is not a refutation of B1, and it would be easy to read it as one.** The
host does not have the problem B1 solves. On the device, where the bound does not
collapse, the same two modes are 3.6 and 111 boxes per query on cloth-funnel's
broad-phase set — **31× apart**. A prepass at 3.6 plus `Tight` on one surviving
query is ~3.6 per query against 111, if the cull rate above carries over from the
curated query sets to the broad-phase candidate sets.

### Priced on the device. B1 buys 1.6x here, not 31x

`sccd_np_trace --device --batch` on a GH200, the same files, one device
`toi_stride=0` call each:

| scene | file | queries | host R | host T | **dev R** | **dev T** | dev T / dev R |
|---|---|---:|---:|---:|---:|---:|---:|
| cloth-funnel | 227vf | 92 | 6.53 | 44.54 | 42.3 | 144.5 | 3.4× |
| cloth-funnel | 317vf | 24 | 77.46 | 283.12 | 231.2 | 1,643.8 | 7.1× |
| armadillo-rollers | 326vf | 4,652 | 13.14 | 16.04 | 85.0 | 165.1 | 1.94× |
| armadillo-rollers | 229vf | 1,597 | 28.59 | 38.57 | 126.0 | 200.4 | 1.59× |
| cloth-ball | 92vf | 19,034 | 7.75 | 7.85 | 94.7 | 152.1 | 1.61× |
| cloth-ball | 91vf | 11,278 | 8.03 | 8.47 | 97.7 | 155.0 | 1.59× |

On the device `Relaxed` is **1.6× to 3.4×** cheaper than `Tight` on the large
files, not 31×. A prepass at 94.7 boxes per query plus `Tight` on the one
surviving query is ~95 against `Tight` alone at 152: **1.6×, for an extra launch
and an extra pass over the candidate set.** Real, but not the order of magnitude
the cull rate suggested.

**Where the 31× came from, and why it should not have been quoted.** It paired
`3.6` boxes per query for device mode 0 with `111` for device mode 2. Those come
from different runs over different populations: 3.6 was measured on the broad-phase
candidate set before the queue rewrite and summed across both `toi_stride` paths,
111 is the post-rewrite zero-stride figure. That is the same error the decision
record already lists twice — comparing across populations — and it produced an
estimate 20× too favourable.

**What is still genuinely open.** These are curated query files: every query in
them is an interesting case. The broad-phase candidate set the library actually
sees is dominated by pairs that are trivially far apart, and `Relaxed` should
dispose of those far more cheaply than `Tight` does. The ratio on *that*
population is what decides B1, it is not measurable with `sccd_np_trace` as it
stands (which reads curated query CSVs), and it wants a clean single-run
measurement of both modes on the broad-phase set through `sccd_bench` — both
strides separated, one allocation, post-queue-rewrite.

Until that exists, **B1 is a 1.6× optimisation with a measured cull rate and an
unmeasured prize on the population that matters.** It should not be built yet.

### A device number worth carrying forward

Batched, on the same call, the device costs 12-19× the host: 94.7 against 7.75 on
cloth-ball 92vf, 152.1 against 7.85 for `Tight`. That is the bound-collapse gap
measured on a single call rather than inferred from a scene aggregate, and it is
the number C1 has to move.

## Family C — reduce the size of the search

**C1. Best-first instead of depth-first.** The goal on the headline path is the
*earliest* time of impact, and DFS finds it last. A frontier ordered by `tlower`
expands the box that could contain the earliest contact first, so the bound
arrives early and prunes the rest — which is precisely what the host appears to
be getting from its ordering and the device is not. On a GPU this wants a binned
or bucketed priority structure rather than a heap: bucket by `tlower`, drain the
lowest non-empty bucket. The double-buffered queue already gives the launch
structure this needs; the buckets replace the single FIFO.

This is the strongest candidate in the document.

### First cut, measured: 2.2x where the search is deep, and it removes the variance

The cheapest form of C1 is in `sccd_narrowphase.cu` behind `SCCD_NP_BEST_FIRST`,
off by default. The block already descends into the child with the smaller
`tlower` and pushes the sibling, so the *descent* is best-first already; what is
not is the **refill**. An idle thread claims the top of the block-shared stack,
which is the most recently pushed box — often a deep box of some other query
whose `t` interval starts far above the answer. `promote_min_tlower` runs one
block-wide argmin over the stack and swaps the minimum to the top before the
claims run, so the next claim takes the most promising box.

GH200, device batched `toi_stride=0`, mode 2, three repeats, medians:

| scene | file | queries | DFS | spread | best-first | spread | boxes |
|---|---|---:|---:|---:|---:|---:|---:|
| cloth-funnel | 317vf | 24 | 1,411.4 | 33% | **646.4** | 1% | **−54%** |
| armadillo-rollers | 326vf | 4,652 | 171.3 | 11% | 162.8 | 4% | −5.0% |
| cloth-funnel | 227vf | 92 | 219.3 | 34% | 211.8 | 0% | −3.4% |
| cloth-ball | 92vf | 19,034 | 152.1 | 1% | 151.8 | 1% | −0.2% |

**Two results, and the second was not the one being looked for.**

The search shrinks by 2.2× where it is deep — 317vf averages 1,411 boxes per
query — and by nothing where the bound already collapses fast. That is the
expected shape: reordering the frontier only matters when there is a frontier.

And **it removes the run-to-run variance**. DFS spreads 33–34% across three
identical runs on the two cloth-funnel files; best-first spreads 0–1% on every
file. Which box a block picks up next currently depends on the order blocks
happened to finish in, and that decides how fast the bound arrives. Ordering the
frontier makes it deterministic. For a project whose retractions are mostly
differences read off inside the noise, that may be worth more than the 2.2×.

`sccd_narrowphase_cuda_test` passes on the best-first build: all 20
configurations conservative, `missed=0 late=0`.

### Time: costs 4–9% where it saves nothing, and the rest is not resolvable

`np_trace --repeat 11` times the device call itself with CUDA events, so process
startup is out. Two independent jobs, medians of 11:

| scene | file | queries | DFS job 1 | DFS job 2 | best-first |
|---|---|---:|---:|---:|---:|
| cloth-ball | 92vf | 19,034 | 2.635 | 2.644 | 2.750 / 2.895 (**+4 to +9%**) |
| cloth-ball | 91vf | 11,278 | 2.498 | 2.380 | 2.430 / 2.496 |
| armadillo-rollers | 326vf | 4,652 | 3.562 | 3.131 | 3.381 / 3.006 |
| cloth-funnel | 317vf | 24 | 2.369 | **3.499** | 1.760 / 1.765 |
| cloth-funnel | 227vf | 92 | 1.550 | 1.643 | 1.653 / 1.671 |

**Only the first row is resolvable.** cloth-ball 92vf is the one file where the
DFS baseline repeats across jobs — 2.635 against 2.644, 0.3% apart. There
best-first costs **4–9% and saves no boxes**, which is the argmin's price with
nothing to buy.

Everywhere else the baseline moves more than the effect. The DFS time on 317vf
went 2.369 → 3.499 ms between two jobs on the same binary and the same 24
queries: **48%**, against a claimed −50%. These kernels are 1.5–3.5 ms on files of
24 to 19,034 queries, and at that size the launch and the drain rounds dominate.
Note that best-first itself repeats to ~0.3% across jobs, consistent with the box
counts — it is the DFS baseline that is unstable, which is the same finding as the
33% box-count spread, seen in time.

A threshold that skips the reorder on shallow stacks
(`SCCD_NP_BEST_FIRST_MIN_TOP`) was tried at 8 and 32 and read *slower* on
cloth-ball — +31% and +16% against +9% for no threshold. That is the wrong sign
for a knob that removes work, it is one sample against an unresolvable baseline,
and it is recorded as unexplained rather than as a result. The default is 1, which
is the configuration actually measured.

### On the real workload it does nothing

`sccd_bench` now builds on Alps (smesh at
`/capstor/scratch/cscs/zulianp/installations/smesh-rel/lib64/cmake/smesh`), so the
broad-phase candidate set can be measured. cloth-funnel, 4 cases, device, mode 2,
two runs of each build:

| run | s0 boxes | s0 ms | s1 boxes | s1 ms | fp | fn |
|---|---:|---:|---:|---:|---:|---:|
| DFS | 4,082,926 | 12.39 | 1,775,276,978 | 160.51 | 75 | 0 |
| best-first | 3,984,686 | 12.03 | 2,274,348,484 | 164.37 | 75 | 0 |
| DFS | 4,049,386 | 11.76 | 2,001,051,794 | 158.07 | 75 | 0 |
| best-first | 4,079,490 | 12.43 | 1,772,120,888 | 165.19 | 75 | 0 |

**No difference.** Stride-0 boxes agree within 2%, time within 5%, and the
stride-1 box count swings 28% between two runs of the *same* binary. False
positives are 75 and false negatives 0 in all four, so the results are identical.

The 54% on cloth-funnel 317vf does not carry over, and in hindsight the curated
files already said so: the reduction appeared only where the search is deep, and
on the broad-phase set stride 0 averages **7.6 boxes per query**. There is no
frontier to reorder. C1's first cut is a no-op on the workload that matters.

**C1 in this form is finished.** A global frontier bucketed by `tlower` would go
further than a block-local swap, but the evidence now says the block-local
frontier is not where the work is, and there is no reason to think a global one
would be either — not on a path averaging 7.6 boxes.

## What the real workload says instead: stride 1 is the whole problem

The same four runs, split by code path:

| path | queries | boxes | boxes/query | ms |
|---|---:|---:|---:|---:|
| `toi_stride=0` | 537,940 | 4.08 M | **7.6** | **12.4** |
| `toi_stride=1` | 510,098 | 1.78 G | **3,480** | **160.5** |

Stride 1 costs **13× the time and 435× the boxes**, and individual calls are far
worse than the average: one call ran 3 queries at 2,060,321 boxes each, another 15
queries at 8,482,418 each.

This is the same imbalance the earlier counting found — 274 queries accounting for
88% of the device's boxes — but measured after the queue rewrite, on the real
candidate sets, with the two paths separated and the time attributed. It is not a
tail on the main path. It is a different kernel, `narrow_phase_dfs_*` with one
block per query and the 128-way seeding dice, and it dominates the narrow phase.

**Everything else in this document is optimising the 12 ms.** The next work is on
the 160 ms:

### ~~C4, the seeding dice~~ — closed. It is not a work multiplier

`SCCD_NP_SINGLE_SEED` replaces the 128-way dice with one root seed on thread 0.
cloth-funnel, 4 cases, device, mode 2, two runs of each:

| build | s0 boxes | s0 ms | s1 boxes | s1 ms | fp | fn |
|---|---:|---:|---:|---:|---:|---:|
| dice | 4,108,982 | 11.75 | 2,030,476,582 | 166.43 | 75 | 0 |
| single seed | 4,071,618 | 11.84 | 2,103,899,876 | **324.50** | 75 | 0 |
| dice | 4,065,004 | 11.42 | 1,748,064,282 | 159.25 | 75 | 0 |
| single seed | 4,083,746 | 10.95 | 2,090,055,042 | **244.97** | 75 | 0 |

**The box counts do not move** — 2.03 G and 1.75 G with the dice against 2.10 G
and 2.09 G without, inside the 28% run-to-run swing. The time is 1.5–2× *worse*,
which is what leaving 127 of 128 threads idle costs.

So the dice does no extra work; it fills the block. The original time-based
verdict was right, and it is now right for a reason that can be stated: removing
the dice changes nothing about the search and only wastes the machine.

It also refutes the hypothesis that motivated re-opening it — that 128 subtrees
starting before any accept exists is what stops the bound arriving. If that were
so, one seed would have cut the box count. It does not move it at all.

Conservativeness holds on the single-seed build: all 20 configurations, `fn=0`.

### The open question, sharpened

Same workload, same 510,098 stride-1 queries, no shared bound on either side:

| | boxes per query |
|---|---:|
| host `Tight`, `toi_stride=1` | **1.22** |
| device `Tight`, `toi_stride=1` | **3,480** |

Neither side has a cross-query bound to help it, so this is not the bound
schedule. It is also not the dice, not the depth cap, not the split rule, not the
acceptance test, and not the seeding — all measured.

And it contradicts the isolated measurement from step 0, where the device was
**2.31×** the host on one query at `max_toi = 1`. Two possibilities, and they are
distinguishable by one look:

1. `sccd_bench` passes the stride-1 call a `max_toi` already collapsed by the
   stride-0 call. The host's root is then `[0, tiny]` and dies immediately, 1.22
   boxes; the device would have to be failing to use the same bound. That is a
   contained bug, and the first thing to check.
2. The two populations differ — the curated query files step 0 used are not the
   broad-phase candidates. That would make 2.31× and 2,852× both correct and
   about different things.

**Checked. Possibility 1 is dead.** `time_narrowphase_per_query` and
`time_narrowphase_zero_stride` in `bench.exe.cpp` both open with
`scalar_t toi = scalar_t(1)`. Both sides receive `max_toi = 1`; nothing is carried
over from the zero-stride call.

So it is a population difference, and the shape of it resolves the paradox. Put
the four numbers together:

| | host | device |
|---|---:|---:|
| curated query files (hard) | 559.88 | 1,291.39 |
| broad-phase candidates (easy) | **1.22** | **3,480** |

The device appears to do *more* work on the easy population than on the hard one,
which cannot be true of a per-query average unless the average is not describing
the population. It is not. Most broad-phase candidates are pairs whose swept boxes
overlap and which never come close; the host rejects them at the root, 1.22 boxes,
and the device spends its 128-box dice on them. If 99.9% of the 510,098 queries
cost 128, that is 65 M of the 1.78 G, and **the remaining ~510 queries carry 3.4
million boxes each**.

That matches the per-call figures directly: one call ran 3 queries at 2,060,321
boxes each, another 15 at 8,482,418. The mean of 3,480 describes nothing.

## Found it: the device's edge-edge domain tolerances were on the wrong axes

The host computes three denominators — the maximum L-infinity displacement of the
eight difference-function corners along `t`, along `u`, and along `v` — and caps
each quotient. The device, still carrying the SymPy generator's output while the
host had been rewritten by hand from TightInclusion, computed **two** and assigned
them `(dl, dl, edge0)`. So `tol1` was the *t* tolerance and `tol2` the *u*
tolerance, and the `v` denominator was never computed at all.

On `worst-query-w1ee`, `tol1` came out **286× too large**.

Nothing unsafe follows. A larger domain tolerance only makes acceptance easier,
which reports an earlier time of impact — the conservative direction — and every
gate stayed green throughout. But **the split rule picks the axis with the largest
width/tolerance ratio**, so a `u` tolerance 286× too large starves `u` of splits.
The search then refines `t` and `v` against a `u` interval pinned at full width
and never converges: 53 million boxes is a binary tree to depth 26, which is what
running until the geometry stops you rather than until the tolerance does looks
like.

### The fix, measured

The device's `compute_edge_edge_tolerance` is now a transcription of the host's,
caps included. One query per call, `max_toi = 1`, stride 1:

| | host | device before | device after | |
|---|---:|---:|---:|---:|
| `w1ee` | 55 | 53,015,756 | **274** | 193,000× fewer |
| `w2ee` | 543 | 47,578,924 | **1,786** | 26,600× fewer |
| `w3ee` | 645 | 2,466,808 | **1,056** | 2,336× fewer |

The device goes from 964,000×, 87,600× and 3,800× the host to **5.0×, 3.3× and
1.6×**.

cloth-funnel, 4 cases, device, mode 2, the whole benchmark:

| | before | after | |
|---|---:|---:|---:|
| stride-0 boxes | 4,082,926 | **1,354,372** | 3.0× fewer |
| stride-0 ms | 12.39 | 12.07 | — |
| stride-1 boxes | 1,775,276,978 | **128,546,876** | **13.8× fewer** |
| stride-1 ms | 160.51 | **37.62** | **4.3× faster** |
| narrow phase, total | 172.90 | **49.69** | **3.5× faster** |
| false positives | 75 | 74 | one fewer |
| false negatives | 0 | 0 | |

`sccd_narrowphase_cuda_test` passes: all 20 configurations conservative. The one
fewer false positive and `w3ee`'s time of impact moving 0.163146973 → 0.174438477
are both the expected sign — a correct, smaller tolerance converges where the old
one accepted early, so the answer is *later* and still at or before the truth.

### After the fix the device and host agree exactly on the classification

Same binary, same cases, same candidate sets, mode 2, only
`SCCD_BENCH_EXECUTION_SPACE` differing:

| scene | type | queries | host fp | device fp | fn |
|---|---|---:|---:|---:|---:|
| armadillo-rollers | vf | 62,420 | **444** | **444** | 0 |
| armadillo-rollers | ee | 409,527 | **170** | **170** | 0 |
| cloth-funnel | ee | 510,022 | **74** | **74** | 0 |

Exact agreement, and zero false negatives on both sides.

This is new. Before the tolerance fix cloth-funnel edge-edge read 75 on the device
against the host's 74 — the extra false positive was the wrong tolerance
accepting a box it should have refined. So the fix is not only 3.5× faster, it
removes a divergence in the answers.

What this does *not* establish: these are counts, not a per-query comparison, and
equal totals could in principle hide two compensating disagreements. `ti_oracle`
compares query by query and gates on it; running it host-against-device on all
three scenes is the stronger check and has not been done since the fix.

### The vertex-face caps: parity, not a win

Both sides use the generated expression there, with three distinct denominators,
so the axes are right. What the device lacked is the **cap**: the host applies
`clamp_domain_tol` to all three and the device applied none.

Added, and measured on armadillo-rollers, 4 cases, device, mode 2:

| | no cap | with caps |
|---|---:|---:|
| vertex-face, stride 0 | 5,749,760 boxes (46.1/q) | **6,113,234 (49.0/q)** |
| vertex-face, stride 1 | 54,131,944 | 53,207,216 |
| narrow ms | 4.03 | 4.06 |
| per-query ms | 36.08 | 33.44 |
| false positives | 444 | **444** |
| false negatives | 0 | 0 |

**It costs 6.3% more stride-0 boxes and buys nothing measurable here.** The false
positive count does not move, and I could not construct a query where the missing
cap changes an answer: a static vertex-face query held clear of the triangle is
rejected by origin containment before the tolerance is consulted at all — 128
boxes and `toi = 1` on both builds.

Kept anyway, and the reason is parity rather than a number. The two sides are
meant to compute the same tolerance, and a divergence between them is exactly
what produced the edge-edge defect above — which sat in every measurement this
project recorded, undetected, because nothing compared the two functions. The
failure mode the cap closes is also real even if these datasets do not reach it:
the host's own note records a near-static query yielding `tol[0] = 1.0e+06`, which
accepts the root box on sight and reports a contact at the root's `t` lower
bound. 0.7% of the cheap path is a low price for removing that.

### The old target, now closed

A few hundred queries per scene cost the device millions of boxes each. The host's
per-query path costs 1,844 boxes on its own hardest queries. Nothing measured so
far explains a query that costs 3.4 million on one machine and thousands on the
other, and the curated query files do not contain such a query — their worst
device query is 4,462 boxes.

**The next step is to get one of them.** `sccd_bench` already has the per-query
counter (`g_np_perq`) and prints a log2 histogram of it; what it does not do is
say *which* query, or let its geometry be dumped. Adding "write the worst N
queries' eight points to a CSV" makes them loadable by `sccd_np_trace`, which can
then run one of them on both machines with a per-box trace. That is the box-by-box
diff the original step 0 called for, now aimed at a query that will actually show
something.

**C2. Bound the search by `t` before refining `u` and `v`.** A 1D pass over `t`
with `u` and `v` at full width is cheap and its rejections are sound; it can
carve the `t` interval down before the 3D search starts. Whether the multilinear
hull is tight enough over full-width `u`,`v` for this to reject anything is an
open question and answerable on paper before any code.

**~~C3. Per-query bound in addition to the global one.~~** Dropped. The device
already prunes against a single global minimum refreshed every loop iteration,
which is *tighter* than the host's per-chunk snapshot, not looser. Reading the two
kernels side by side settles it without a measurement: the device's bound
discipline is the better of the two, and adding a per-query bound on top can only
be weaker than the global one it already has.

**C4. Re-run the seeding measurement on the `toi_stride=1` path.** The single-root
seed was measured at 803 → 871 ms, which is why the dice stayed. That measurement
predates knowing the path costs 2.58 million boxes per query. It should be
repeated with box counts rather than only milliseconds.

## Family D — scheduling

**D1. Make the global queue the primary distributor.** It is now wait-free and
double-buffered, and the shared stack is 64 entries. The remaining question is
whether the shared stack should be a cache in front of the queue rather than the
first-class store it still is.

**~~D2. Understand the capacity-1024 anomaly.~~ It no longer exists.** It was
recorded as 4,059 ms on cloth-funnel against the old queue's 634 at the same
capacity — a 6.4× cliff. Re-measured after the tolerance fix and the
thread-per-query change, on a clean build, mode 2, device, three repeats,
medians, `SCCD_NP_SHARED_STACK_CAP` being a compile-time constant so each is its
own build:

| scene | | 4 | 8 | 16 | 32 | 64 | 256 | 1024 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| cloth-funnel | s0 | 26.5 | **25.3** | 26.8 | 27.4 | 31.1 | 30.3 | 32.0 |
| cloth-funnel | s1 | 33.8 | 32.2 | **28.3** | 30.4 | 33.5 | 38.5 | 38.8 |
| armadillo-rollers | s0 | 34.5 | 32.9 | **31.8** | 34.3 | 37.2 | 47.2 | 46.7 |
| armadillo-rollers | s1 | 79.1 | 70.3 | **60.6** | 64.6 | 68.4 | 108.4 | 134.6 |
| cloth-ball | s0 | 114.7 | 106.3 | 103.8 | **102.7** | 104.5 | — | — |
| cloth-ball | s1 | 170.8 | 182.0 | 172.9 | **158.5** | 159.9 | — | — |

No cliff anywhere. 1024 is 1.2× and 2.0× the shipped capacity on the two smaller
scenes, which is a gentle penalty for hoarding work in a block instead of
spilling it, not a 6.4× regression. False positives are 260 / 5,637 / 95,424 at
every capacity, `fn=0` throughout.

**Why it went away.** The anomaly was almost certainly a symptom of the
edge-edge tolerance defect. That made single queries generate subtrees of
millions of boxes, and a large shared stack keeps such a subtree inside the one
block that started it instead of spilling it to the global queue for every other
block to help with — so the bigger the stack, the more completely one block
starved while the rest idled. With the subtrees three to fourteen times smaller
there is nothing left to hoard.

**The shipped capacity stays at 64.** 16 and 32 look 2–3% better on the sum of
all six cells, and the repeat-to-repeat spread on that same sum is 1–3%: the
gap is the size of the noise. Changing a shipped default on a difference that
cannot be resolved is the exact mistake this project's decision record is full
of. If someone wants it, the measurement to make is many more repeats at 16, 32
and 64 in one allocation — not this one again.

## C: the end-to-end table, and what it shows

The reduction was recovered from `mode-stride-matrix.csv` the same way as the
narrow-phase one: **`broad_ms + narrow_ms_s0`, summed over the 16 cases within a
repeat, median over repeats, minimum over the two modes.** No `prep_ms` — that
was the first thing tried and it overshoots by 25 ms on cloth-funnel.

Freshly measured, and split by stride because the two now say different things:

| | scene | CPU | GPU | |
|---|---|---:|---:|---|
| `find_earliest_impact_time` | cloth-funnel | **29.1** | 41.2 | 1.4× CPU |
| | armadillo-rollers | **47.9** | 48.7 | parity |
| | cloth-ball | 421.6 | **162.8** | 2.6× GPU |
| `find_impact_times` | cloth-funnel | **32.6** | 39.2 | 1.2× CPU |
| | armadillo-rollers | 95.6 | **79.3** | 1.2× GPU |
| | cloth-ball | 547.7 | **215.7** | 2.5× GPU |

**The earliest-impact row has not moved** — 29.1 against a published 29.2, 47.9
against 47.8, 162.8 against 159.4. Everything this branch fixed was on the
per-query path, and a pipeline built on `find_earliest_impact_time` sees none of
it. That is worth saying plainly rather than letting a reader infer a
speedup that is not there.

The per-query row is the one that moved: the GPU went from losing on all three
scenes to winning on two. The broad-phase table is unchanged and reproduces to
within 2%, which is the control that says this run is comparable to the old one.

## D: the residual, and the counter that was supposed to explain it

Two things were suspected: that the two counters count different units, and that
something in the per-query search still differed.

**The counting asymmetry is real and negligible.** The host drops a box in its
refill loop when the bound has already passed it, *before* classifying it, so it
never ticks; the device evaluates both children of a split and ticks them, then
applies the equivalent filter. `g_np_bound_killed` counts exactly those, so
`corner_evals - bound_killed` is the figure that means the same thing on both
machines. Measured on cloth-funnel:

| | evaluations | bound-killed | share |
|---|---:|---:|---|
| device, stride 0 | 3,762,296 | 922 | **0.025%** |
| device, stride 1 | 7,708,066 | 853 | **0.011%** |

So the counters have been comparable all along, and every ratio computed from
them stands.

**The residual itself is gone**, and item A took it. One query per call,
`max_toi = 1`:

| | host | device now | device after the tolerance fix | before it |
|---|---:|---:|---:|---:|
| `w1ee` | 55 | **55** | 274 | 53,015,756 |
| `w2ee` | 543 | **543** | 1,786 | 47,578,924 |
| `w3ee` | 645 | **837** | 1,056 | 2,466,808 |

Exact agreement on two of the three, 1.30× on the third. The 5.0×/3.3×/1.6× that
this item was opened on was the 128-way dice's floor of 128 boxes per query;
running stride 1 one thread per query removed it.

**What remains is not a defect.** At scene level the device still classifies 5.3×
(stride 1) and 5.5× (stride 0) the host's boxes on cloth-funnel. That is the
bound-collapse-rate difference step 0 identified: the host hands out queries in
chunks and a few early chunks collapse the bound for the hundreds of thousands
that follow, while the device starts every query at once. It is a property of
running queries in sequence, not of the search — which is now the same search on
both machines, as the isolated queries show.

## What must not change

- **Conservativeness.** A reported time of impact is at or before the true one,
  and a collision is never missed. Accepting a box is always safe; refining is
  always safe; exhaustion must accept at the box's `t` lower bound. The only
  operation that can lose a root is the origin-containment rejection, and it must
  stay padded by `(vf ? 30 : 28) * eps * min(max_coord, 1)^3`. Any variant that
  reports later than the reference is a defect in that variant, not a trade-off.
- **The root finder computes in double.** Single precision was measured
  non-conservative on armadillo-rollers — nine times of impact later than the true
  root. The interface may stay templated on `T`; the search may not.
- **Quads keep their own path.** They are first class, they have their own
  inclusion function and device kernel, and no variant here may leave them behind.
- **The shipped library keeps zero external dependencies.**

## How to measure

- **One allocation.** Between `srun` allocations this harness varies by ~40%, so a
  difference smaller than that measured across two of them is not a result.
  `benchmark/scripts/assess.sbatch.sh` interleaves variants inside one job.
- **Box counts alongside milliseconds.** `-DSCCD_NP_COUNT_BOXES`. A variant that
  moves time without moving boxes has not done what it claims — that is exactly
  how the batching result fell over.
- **Split by `toi_stride`.** The two device kernels are separate and an aggregate
  over both is an average over populations that differ by 15×.
- **Gate every variant with `ti_oracle`**, which exits non-zero on a missed
  collision or a late time of impact, and with `sccd_narrowphase_cuda_test`.
  Speed that costs the invariant is not a result either.

## Order

Every numbered item below the line has been closed with a measurement. What
remains is listed first.

### Open, largest first

| | | |
|---|---|---|
| ~~**A**~~ | ~~Stride 1 on the thread-per-query kernel~~ | **Done, and it is the default.** 4.47×, 2.76× and 4.43× on the three scenes, identical false positives, `fn=0`, gate green both ways. Detail below. |
| ~~**B**~~ | ~~D2, the capacity-1024 anomaly~~ | **Gone.** It no longer reproduces; the curve is smooth and shallow from 4 to 1024, and no capacity in that range differs from another by more than the run-to-run spread. D1 is unblocked. Detail below. |
| ~~**C**~~ | ~~Refresh the end-to-end table~~ | **Done.** The reduction is `broad_ms + narrow_ms`, no prep, min over modes — recovered the same way as the narrow-phase table. It now has a row per stride, because only one of them moved. |
| ~~**D**~~ | ~~The per-query residual~~ | **Gone.** The device now matches the host box-for-box on two of the three worst queries and is 1.30× on the third. The counting asymmetry it was blamed on is 0.011–0.025%. Detail below. |
| **E** | **A1/A2, per-box cost** | Still last. The stride-1 path is scheduling-bound, so 4× of arithmetic buys nothing there either. |
| — | *Loose end* | The host's mode 0 increments the box counter and nothing prints it — the `fprintf` lives in `narrow_phase_tight_*` and the scalar path has none. Host `Relaxed` box counts are not obtainable from a run. |

### Closed

| | | |
|---|---|---|
| 0 | Establish the cause | The scene-level gap was the bound schedule; the per-query gap was a tolerance on the wrong axis |
| 1 | Check the vertex-face tolerance | Clean. 26,423 vertex-face and 125,525 edge-edge queries against TightInclusion's source; worst relative difference 5.8e-13 and 0 |
| 2 | Re-measure `docs/BENCHMARKS.md` | Done, and it corrected two of my own numbers — the wrong column and an instrumented build |
| 3 | B1, `Relaxed` prepass | The prize is gone. `Tight` costs 1.12–1.53× `Relaxed`, and less on two scenes |
| 4 | The per-query path | Root pre-test landed: 2.9–5.6× less work, no time change. The non-result identified the real bound |
| C1 | Best-first | A no-op on the broad-phase set. Kept behind `SCCD_NP_BEST_FIRST`, off |
| C3 | Per-query bound | Dropped: the device's bound is already tighter than the host's |
| C4 | The seeding dice | Not a work multiplier; it fills the block. Removing it is 1.5–2× slower |

### What the fix put at the top instead

### ~~1. Check the vertex-face tolerance against TightInclusion~~ — clean

The concern was real and the answer is no. The host's *edge-edge* tolerance was
rewritten by hand from TightInclusion; its *vertex-face* tolerance is still the
SymPy generator's output, and so is the device's — so the two machines would agree
just as well if the generated expression were wrong, and every comparison run so
far has been host against device, which cannot see a shared error.

Both were checked against TightInclusion's own source directly
(`_deps/tight_inclusion-src/src/tight_inclusion/ccd.cpp`, `compute_vertex_face_tolerances`
and `compute_edge_edge_tolerances`), evaluated on real dataset queries:

| | queries compared | worst relative difference |
|---|---:|---|
| vertex-face, the generated expression | 26,423 | **5.81e-13** |
| edge-edge, the host's hand-written version | 125,525 | **0.00e+00** |

5.8e-13 is floating-point reassociation — the generated form computes the same
maxima through different intermediates — not a semantic difference. The
vertex-face generator got the bilinear corner right, including
`p011 = v - (f1 + f2 - f0)`, which is the part of the vertex-face formula that is
easy to get wrong and the reason the check was worth running.

So the edge-edge device defect was a transcription failure in one place, not a bad
expression propagated from the generator. Nothing else in the tolerance path is
carrying the same class of error.

### ~~2. Re-measure `docs/BENCHMARKS.md`~~ — done, and it corrects two of my own numbers

The methodology was recovered from `benchmark/assessment/mode-stride-matrix.csv`
by reproducing the published table from it: **sum `narrow_ms_s0` (or
`narrow_ms_s1`) over the scene's 16 cases within one repeat, then take the median
over the three repeats**, vertex-face and edge-edge summed together. That
reproduces all 24 published cells exactly. My first attempt took a median over all
48 rows, which is a different statistic and is where the 2.4× came from.

Two errors of mine fall out of getting this right, and both were in numbers I
reported earlier in this investigation:

- **I was reading the wrong column.** `sccd_bench` emits `narrow_ms` at field 7
  and `query_narrow_ms` at field 8; the stride-1 timing is `narrow_ms_s1` at field
  **13**. `query_narrow_ms` is the accuracy-check pass, not the stride-1 narrow
  phase. Every "stride-1 ms" I quoted was field 8.
- **I was timing an instrumented build.** All of it ran under
  `-DSCCD_NP_COUNT_BOXES`, which this document already warns puts a global atomic
  on the hot path. Measured now: it makes the host 20× slower. So the "172.90 →
  49.69 ms, 3.5×" figure was the wrong column measured on a build whose timings
  are meaningless.

The real table, clean build, correct reduction:

| | was | now | |
|---|---:|---:|---:|
| armadillo-rollers, GPU s1 | 3126.7 | **188.1** | **16.6×** |
| cloth-funnel, GPU s1 | 789.9 | **155.7** | 5.1× |
| armadillo-rollers, GPU s0 | 56.3 | **37.1** | 1.5× |
| cloth-funnel, GPU s0 | 36.5 | **28.3** | 1.3× |
| cloth-ball, GPU s1 | 902.8 | **695.7** | 1.3× |

The CPU columns reproduce the old measurement to within 1–7%, which is the check
that the clean build and the recovered reduction are both right.

Still outstanding: the end-to-end table, whose GPU column sums prep, broad phase
and narrow phase. Refreshing it needs the broad-phase timings from the same run,
which this one did not capture.

### ~~3. Re-price B1 on the post-fix device~~ — closed. The prize is gone

B1 culls queries with a cheap `Relaxed` pass before running `Tight` on the
survivors, so it is worth exactly the ratio between the two modes' cost. Measured
on the device over the real broad-phase candidate sets, post-fix:

| scene | stride | `Relaxed` | `Tight` | `Tight`/`Relaxed` |
|---|---|---:|---:|---:|
| cloth-funnel | 0 | 2,771,806 | 3,758,494 | 1.36× |
| cloth-funnel | 1 | 134,673,458 | 153,543,274 | 1.14× |
| armadillo-rollers | 0 | 11,645,472 | 17,770,366 | 1.53× |
| armadillo-rollers | 1 | 792,280,816 | 532,188,788 | **0.67×** |
| cloth-ball | 0 | 65,927,052 | 73,512,752 | 1.12× |
| cloth-ball | 1 | 6,346,807,782 | 6,042,530,180 | **0.95×** |

`Tight` costs 1.12–1.53× `Relaxed` on the earliest-impact path, and on the
per-query path it is **cheaper than `Relaxed` on two of three scenes**. A prepass
would therefore cost about what it saves and roughly double the work.

**B1's entire economics were the tolerance defect.** Before the fix `Tight` was
doing an order of magnitude more work than it needed on the per-query path, and
that inflated gap is what made a cheap prepass look worth 30×, then 1.6×. With the
defect gone the two modes cost about the same and there is nothing to arbitrage.

### `Relaxed` is no longer the cheaper mode on the device

Worth stating on its own, because the documentation still describes the modes as
a speed-for-tightness trade. On the device that is now false: `Relaxed` costs
**1.5× more** boxes than `Tight` on armadillo-rollers at stride 1, and about the
same on cloth-ball, and the timing table agrees — GPU s1 is 234.5 ms for `Relaxed`
against 188.1 for `Tight`.

The reason is the same asymmetry the whole design rests on. `Relaxed` accepts a
box on a looser test, which ends *that* box sooner but reports a time of impact
further before the true one — and a looser bound prunes less, so the queries that
follow do more work. `Tight` pays more per box and buys a bound that kills more
boxes. On the host, where the bound collapses within a few chunks anyway, the
trade still runs the documented way; on the device it does not.

## A: stride 1 runs one thread per query

The zero-stride body is now templated on `per_query`. With it set, each thread
carries its own query's bound in a register instead of pruning against a
block-shared minimum, adopts the new query's bound when it pops another query's
box off the shared stack, and publishes to `toi[qid]` on the way out and on every
switch. Everything else — the search, the acceptance tests, the stack, the global
queue — is the same code the earliest-impact path runs.

`toi_stride=1`, mode 2, GH200, three repeats, medians:

| scene | block per query | thread per query | |
|---|---:|---:|---:|
| cloth-funnel | 147.2 | **33.0** | **4.47×** |
| armadillo-rollers | 199.8 | **72.4** | **2.76×** |
| cloth-ball | 696.3 | **157.0** | **4.43×** |

False positives 260 / 5,637 / 95,424 with `fn=0` on both variants, identical to
each other and to the host. All 20 configurations conservative either way.
Stride 0 is untouched — 27.3 against 27.6, 37.1 against 37.4, 110.2 against
110.4 — which is the check that the change is confined to the path it targets.

The per-query path now costs about what the earliest-impact path costs: 1.2×,
1.9× and 1.4× of it, against 5.4×, 5.4× and 6.3× before. The host's own ratio is
1.6–4.4×, so the anomaly is gone rather than reduced.

`SCCD_NP_S1_BLOCK_PER_QUERY=1` restores the old kernel. It is kept for two
reasons: the new default has one session's measurement behind it, and a query
heavy enough to genuinely want a whole block is what the global queue exists for
— if such a workload appears, that is the switch to try first.

### 4. The per-query path: its cost is not its work

Post-fix, `toi_stride=1` costs 5.1–6.6× `toi_stride=0` in time on the same
candidate sets. In boxes it was 32–83× — and 128 of the 166–182 boxes per query
were the seeding dice, **70–77% of everything that path did**. (The "3.2%" figure
recorded earlier was measured before the tolerance fix, against a total inflated
by it.)

The block-per-query kernel never tested the root box: it went straight to dicing
the root into 128 subboxes. The multilinear hull over a box is exact, so if the
root does not contain the origin then no subbox of it does either — and most of a
broad-phase candidate set is pairs that never come close. Testing the root first,
**reject only** so no reported value can change, removes those 127 evaluations.

| scene | s1 boxes/query before | after | | GPU s1 ms before | after |
|---|---:|---:|---:|---:|---:|
| cloth-funnel | 182.0 | **63.4** | 2.9× | 155.7 | 151.0 |
| armadillo-rollers | 165.7 | **56.9** | 2.9× | 188.1 | 204.2 |
| cloth-ball | 180.8 | **32.3** | 5.6× | 695.7 | 694.3 |

**The work falls 2.9–5.6× and the time does not move.** False positives are
260 / 5,637 / 95,424 with `fn=0`, identical to the host and to the previous build.

That is the finding, and it is worth more than the optimisation. The dice's 128
evaluations are done by 128 threads in one block-wide step, so removing them saves
one step of latency, not 128 units of work. **The block-per-query kernel is not
work-bound; it is bound by having one block per query** — 843,414 blocks on
cloth-funnel, each doing 63 boxes spread over 128 threads, which is half a box per
thread. The cost is scheduling them.

Which explains the whole stride-0-versus-stride-1 gap without any reference to the
search: the same queries cost 28.3 ms through a kernel that runs one *thread* per
query and 151.0 ms through one that runs one *block* per query.

**So the direction for stride 1 is fewer blocks, not less work.** Run the
per-query path on the thread-per-query kernel — it already exists, it already
carries a per-query `toi` index (`toi_idx = qid * toi_stride`), and its only
missing piece is that it currently writes a single shared minimum — and hand a
query to a whole block only once it has proven itself hard. That is the next
piece of work and it is a larger change than anything else in this document.

The root pre-test is kept. It buys no time today, but it costs one evaluation on
a surviving query, it removes 70% of the path's arithmetic, and it makes the
kernel's cost model legible: with it gone, the box count no longer flatters the
kernel's real bottleneck.
