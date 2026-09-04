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

### The target, finally specific

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

**D2. Understand the capacity-1024 anomaly first.** Under the new queue,
`SCCD_NP_SHARED_STACK_CAP=1024` measures 4,059 ms on cloth-funnel against the old
queue's 634 ms at the same capacity — a 6.4× regression at a setting that is no
longer the default. Every capacity at or below 256 is faster than the old code on
every scene, so the shipped configuration is strictly better, but something about
the drain loop is not understood, and D1 tunes exactly that machinery.

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

0. **Done.** The scene-level gap is the bound schedule; the per-query search is
   2.3×.
1. **B1.** It manufactures the collapsed bound in one cheap pass, with no
   sequence to wait for. Directly on the mechanism.
2. **C1.** The same mechanism from the other side: best-first makes the bound
   arrive early *within* the launch.
3. **The 2.3× residual**, starting with whether the two counters count the same
   thing. Small and contained.
4. D2, because D1 cannot be tuned honestly until it is understood.
5. **A1/A2 are demoted.** 4× of arithmetic behind a 2.3× box factor is not where
   the time is.

Also outstanding, and cheap: **the host box counter only covers the `Tight`
kernel.** `SCCD_NP_HOST_BOX_TICK` lives in `sccd_narrowphase_tight.hpp`, so mode 0
reports zero boxes and `Relaxed`'s per-query cost — the thing B1's economics turn
on — cannot currently be measured on the host at all.
