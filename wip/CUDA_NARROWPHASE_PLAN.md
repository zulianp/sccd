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

**B1. `Relaxed` prepass to cull queries.** *Now the highest-value item, given the
step 0 result: it manufactures the collapsed bound the host gets from sequencing,
in one cheap pass, with no sequence.* `Relaxed` averages 3.6 boxes per query
against `Tight`'s 944, and its answer is *at or before* the true time of impact.
So for the earliest-impact query: if a query's `Relaxed` time of impact is later
than the best `Tight` answer found so far, its true time of impact is later too,
and the query cannot hold the earliest impact — it can be discarded outright.
That is sound, exploits the invariant's asymmetry rather than fighting it, and
costs one cheap pass over the candidate set.

The prize depends on how concentrated the earliest impact is. Worth a cheap
offline estimate from existing benchmark output before writing any kernel: how
many queries have `relaxed_toi` below the scene's global earliest impact?

## Family C — reduce the size of the search

**C1. Best-first instead of depth-first.** The goal on the headline path is the
*earliest* time of impact, and DFS finds it last. A frontier ordered by `tlower`
expands the box that could contain the earliest contact first, so the bound
arrives early and prunes the rest — which is precisely what the host appears to
be getting from its ordering and the device is not. On a GPU this wants a binned
or bucketed priority structure rather than a heap: bucket by `tlower`, drain the
lowest non-empty bucket. The double-buffered queue already gives the launch
structure this needs; the buckets replace the single FIFO.

This is the strongest candidate in the document, and it is also the one step 0
can most cleanly confirm or kill.

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
