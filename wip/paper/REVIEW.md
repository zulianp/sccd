# Referee report

**Manuscript:** *A Conservative Continuous Collision Detection Pipeline for
Triangle and Quad Meshes on CPU and GPU*
**Recommendation:** Major revision.

## Summary

The paper presents a two-phase CCD pipeline with host and device implementations,
argues conservativeness structurally, and evaluates it against TightInclusion on a
dataset with exact roots. The central observation — that a branch-and-bound search
which keeps the minimum over accepted leaves is order-independent, and that this
licenses filling vector lanes and thread blocks with work from unrelated queries —
is a good one, cleanly stated, and the two implementations that follow from it are
genuinely different from prior GPU work. The accuracy result (host output identical
to the reference on all twelve scene-phases) is strong and well supported. The
candour of Section 9 is unusual and welcome.

I nevertheless have four objections I consider blocking, and a number of smaller
ones.

## Major

**M1. A conclusion in Section 8.6 is not supported by its own table.**
The paper states "Neither wins. The sweep takes armadillo-rollers, cloth-ball,
n-body and rod-twist by 1.03x to 1.88x". Table 9 has four numeric columns: *prep*
and *broad* for each strategy. Adding the two, as one must to compare the cost of
running a broad phase, the cell list wins **five of six scenes** and ties the
sixth: armadillo-rollers 11,420 against 29,868 ms, cloth-ball 2,175 against 3,761,
cloth-funnel 7,847 against 21,800, puffer-ball 46,233 against 137,801, rod-twist
95,445 against 190,556, and n-body 9,832 against 9,735. The sweep's sort dominates
its cost on every real scene, and the quoted ratios evidently rank the *traversal*
column alone. The table's own caption compounds this by asserting that *broad* is
"the whole broad phase including it", which the numbers contradict. Either the
ratios or the caption is wrong, and in both readings the section's headline
conclusion does not follow. This must be resolved before the paper can stand.

**M2. The closest prior work is cited but never compared against.**
Reference [2] is a GPU CCD pipeline evaluated on the same six scenes with the same
ground truth, and several of this paper's authors are its authors. A reader's first
question is how the present pipeline compares to it, and the paper never answers.
TightInclusion is a CPU reference; measuring against it establishes accuracy but
says little about whether this parallel design improves on the state of the art in
parallel CCD. If a comparison is impractical, the paper must say so explicitly and
say what would be required.

**M3. Proposition 1 is applied beyond its hypotheses.**
The proposition fixes an acceptance predicate and a splitter and concludes
order-independence. Section 8.4 then reports that the host and device answers
differ, and correctly attributes this to the two evaluating the inclusion function
through different floating-point operation sequences. But that means host and
device do *not* share the accepted set A, so the proposition never applied across
them. This is not an error in the result, but the paper currently invites the
reader to expect exact host-device agreement and then explains its absence
defensively. State the scope once, up front.

The corollary — that TightInclusion's answer is also the minimum over the accepted
set — is asserted rather than argued. TightInclusion returns the first box its
priority queue accepts; one needs the queue to be ordered by the same key that the
minimum is taken over, and one needs its own pruning not to remove the minimiser.
A sentence of justification is needed.

**M4. Quads are claimed as a first-class contribution and are never evaluated.**
Section 5.4 describes a separate inclusion function, host root finder and device
kernel for bilinear quads, and Section 3.3 introduces an uncertified constant for
them. Not one number in Section 8 concerns a quad mesh: all six benchmark scenes
are triangle meshes. Either evaluate the quad path or stop presenting it as a
result; at minimum, Section 9 must record that it is unevaluated, since an
uncertified error constant with no measurements behind it is exactly the kind of
claim this paper is otherwise careful about.

## Minor

**m1.** Section 8.1 states that the harness varies by "about 40%" and that "a
difference smaller than that is not reported as a result". Taken literally this
invalidates almost every ratio in the paper, including the 1.64x-3.00x of
Section 8.5. The 40% is evidently between job allocations, not within one. Separate
the two and state the within-allocation figure, since that is the floor that
actually applies.

**m2.** The claim that the two broad phases return identical candidate sets is
presented as an empirical observation over 25,538 combinations. It looks provable:
both evaluate the same exact floating-point predicate on the same boxes, both
exclude vertex-sharing pairs, and the canonical-cell rule is stated to select
exactly one cell. If it is provable, prove it — a measured coincidence is much
weaker than an identity, and the whole racing argument rests on it.

**m3.** Theorem 1 hypothesis (i) requires the rejection test to be padded *by* the
certified bound, but Section 5.1 says the device pads by max(tau, epsilon), which
is wider. Widening is sound; the hypothesis should say "by at least".

**m4.** Section 4.4 lists "five cheap predictors", but item 5 is a fixed choice,
which is not a predictor. Also, item 5's numbers ("the sweep wins two of three")
contradict Table 9 under either reading. If they come from a different run, say so
and say what differed.

**m5.** Table 1's device column is disclaimed in four separate ways in its own
caption. A number that requires that much hedging is not evidence. Either drop the
column or explain in the text what it is doing there.

**m6.** Nine bibliography entries are never cited.

**m7.** Figure 1 labels both queue halves "buffer A", which is correct but reads as
a mistake on first encounter.

**m8.** The Relaxed mode produces 243,717 false positives on rod-twist against
1,245 for Tight, on 283,655 compared queries. That is a large fraction of all
queries, and it deserves more than the half-sentence it gets, since it bears
directly on whether the mode should ship.

## Things done well

The signed accuracy criterion, and the insistence throughout that a late time of
impact is as bad as a miss, is the right framing and is applied consistently. The
negative results — five failed predictors, the chord splitter that lost, the
hoisting that was slower end to end — are more useful than most positive results in
this literature. The reproducibility apparatus, with every table traceable to a
committed CSV, is exemplary and I would like to see it become normal.

---

# Response to the report

Every point was acted on. What changed, in the referee's numbering.

**M1 — the Section 8.6 conclusion.** The referee is right, and the error was
inherited: the `faster` column ranks the *traversal* column alone, while the
table's upstream caption claims it covers the whole broad phase. Section 8.6 is
rewritten to give both readings and to say that they disagree — on the traversal
the two trade wins, on preparation plus traversal the cell list wins five scenes
of six and ties n-body, with the sums spelled out so a reader can add the columns
themselves. The table's caption and note are corrected for this article (the
generator that produced them belongs to the library's documentation and is not
changed here; the override and its reason are recorded in
`scripts/make_inputs.py`). The abstract, introduction, Section 4.4 and the
conclusion were realigned. The disagreement now carries more weight than the old
claim did: it is a direct argument for measuring rather than choosing.

**M2 — no comparison with [2]. NOT SATISFIED.** The manuscript now states plainly
that the comparison is absent, what TightInclusion does and does not establish,
and what making the comparison would require, and it qualifies the speedups as
margins over a host reference rather than over prior parallel work. The comparison
itself has not been made. See *Outstanding* below.

**M3 — Proposition 1 applied beyond its hypotheses.** A new paragraph after the
proposition (now Proposition 2) states that it fixes the accepted set, and that
two searches agree only if they compute the same values bit for bit — which the
host and the reference do by construction and the device does not. The corollary
about TightInclusion is now Remark 1, with the argument the referee asked for.

**M4 — quads unevaluated. PARTIALLY SATISFIED.** The claim has been withdrawn:
the quad path is now described as a design and marked unevaluated where it is
introduced and again in Section 9, together with the fact that its error constant
is uncertified, and the framing "first-class input" is gone. The referee's first
alternative --- evaluate it --- has not been taken. See *Outstanding* below.

**m1 — the noise floor.** Section 8.1 now separates between-allocation variance
(~40%, the reason comparisons stay inside one allocation) from within-allocation
repeat spread (0.3–3.4%, measured, the floor that actually applies).

**m2 — identical candidate sets.** Now Proposition 1, with a proof. The measured
agreement is kept as a check on the correspondence between algorithm and code,
with a pointer to the sweep bug where that correspondence failed.

**m3 — Theorem 1(i).** Now reads "widened by *at least*", and names the device's
`max(tau, epsilon)`. A fourth hypothesis was also added: the proof uses that the
children of Split tile their parent, which was not previously stated.

**m4 — "five predictors".** Now "four predictors and one constant". Item 5 cites
both runs and says they disagree, rather than quoting one of them as fact.

**m5 — Table 1's device column. PARTIALLY SATISFIED.** A paragraph now says what
the column is for (the size of the tail a work list must survive) and what it is
not evidence of. The referee's other alternative --- dropping it --- was not taken,
on the grounds that the tail size is what motivates Section 6.2 and no other
measurement in the paper gives it. See *Outstanding* below.

**m6 — uncited entries.** All 36 entries are now cited.

**m7 — Figure 1.** Relabelled: buffer A is written by launch k and read by launch
k+1, which also writes buffer B for launch k+2.

**m8 — Relaxed false positives.** Expanded with the rates: 243,717 of 283,655
compared queries on rod-twist, 27,374 against Tight's 558 on puffer-ball, and the
conclusion that on those scenes the mode is not a speed-for-accuracy trade at all.

---

# Outstanding

Three of the referee's points are not, or not fully, met. They are recorded here
because the manuscript states their consequences but cannot close them, and a
later revision should either close them or inherit this list.

**M2 — no numerical comparison against Belgrod et al. (open).**
The manuscript documents the gap rather than filling it, and qualifies every
speedup accordingly. Closing it requires work not done here: building that
pipeline on the same \textsc{GH200} node; matching tolerance, depth cap and
output mode so the two searches do comparable work; and timing both inside one
job allocation, since between-allocation variance on this system is about 40%.
Until that is done the paper establishes accuracy against a certified reference
and speed against a host implementation, and says nothing about its standing
against prior parallel CCD.

**M4 — the quad path is unevaluated (open).**
The referee offered two remedies and the manuscript takes the weaker one: the
claim is withdrawn rather than the evidence supplied. What is missing is a quad
benchmark — a scene set of quad meshes with exact roots, which the dataset used
here does not contain, and which would have to be generated. The uncertified
constant $\kappa = 38$ of Section 3.3 is a second open item of the same kind: it
is chosen to exceed both certified constants and is therefore safe, but a
derivation for the bilinear form does not exist and no measurement here bounds
how loose it is.

**m5 — Table 1's device column (partially addressed).**
The column is kept, with its provenance and its limits stated in the text. It
remains a measurement taken before the redesign it motivated, aggregating two code
paths of which 274 queries on one account for 88% of the boxes counted. A clean
replacement — the same distribution measured on the current kernel, split by path
— would let the table carry the argument without the caveats, and has not been
taken.

**A further gap, now closed: memory and host counters.**
An earlier revision reported what the device kernel is limited by using only the
kernel's own instrumentation, and could say nothing about memory or about Grace.
Both have since been measured on a GH200 node and are reported in
\cref{sec:parallel:limits}, with the raw profiler output published under
`benchmark/results/profile/` and the table generated from it.

The finding changed the section's conclusion. The device kernel is not
bandwidth-bound either: DRAM use runs between 0.08% and 1.2% of sustained peak
across the six scenes, and warps stalled on an execution dependency stay between
1.52 and 1.65 per issue-active cycle on every one of them while the grid size
varies by three orders of magnitude. Grace is the opposite case, retiring 2.8-3.8
instructions per cycle at an L1 miss rate of 0.2-2.1%, with the narrow phase
alone missing L1 on 0.01-0.65% of loads.

Two numbers in earlier drafts of that section were wrong and the measurements
corrected them, which is the argument for running the profiler rather than
reasoning about it. DRAM use is a range and not the flat 0.08% that two scenes
suggested. And a libm cost first inferred from unresolved sample addresses is
14.1% of cycles on one scene and 2.1% on another, not the uniform figure the
inference implied; attributing samples by shared object rather than guessing the
library from the address is what settled it.

What remains open is scope: one narrow-phase mode, the earliest-impact kernel
only, and host counters that are per-process rather than per-phase. That is
recorded in Section 9 of the manuscript.

**A further item the referee did not raise.**
The report of the harness overhead in Section 9 (roughly 0.9 s of measured work
per case inside about 21 s of wall clock on the largest scene) means the timings
are sound as timings but are not a model of the pipeline's cost inside a solver.
The paper is explicit about this, and an in-solver measurement remains the
strongest missing evidence for its practical claims.

---

# Figures from the reference study

Belgrod et al. report fourteen figures. Ours now cover the same ground except
where the data cannot be produced here:

| Theirs | Ours |
|---|---|
| 6.1 row 1, per-frame performance | `per-frame.pdf` |
| 6.1 row 2, timing distributions | `results-grid.pdf` |
| 6.1 row 3, false positives | not drawn: with one method and a distribution that is zero on most steps, the box plot shows less than the conservativeness table |
| 6.1 row 5, memory | not measured; the harness records no peak memory |
| 7.1, runtime against query count | `narrow-per-case.pdf`, against candidate pairs |
| 7.2, strong scaling | `strong-scaling.pdf` |
| 7.3, time-of-impact error | `toi-error.pdf` |
| 7.4, several architectures | one node was available, so no |
| 7.5, runtime breakdown | `runtime-breakdown.pdf` |
| 7.6, batching against memory | no batching sweep; the driver has no memory-limiting knob |
| A.1, voxel-size sensitivity | not applicable: the cell size is derived rather than chosen (Section 4.3), so there is no parameter to sweep |
| 1.1, 4.1, 8.x, simulation renders | no renderer in the repository |

---

# Scientific-writing pass

The manuscript was revised against the criteria of the `scientific-writing`
skill, first as a checklist pass and then as a full rewrite of the prose in
every section. Content was held fixed throughout: no claim, number, citation,
cross-reference, table or figure changed, and `make check` passes against
`docs/BENCHMARKS.md` before and after.

| Criterion | Before | After |
|---|---|---|
| Abstract 100-250 words, one flowing block, no labelled sections | 324 words, 35.9 words/sentence | 258 words, 22.7 words/sentence |
| No bullet lists outside Methods | `itemize` of four contributions in Section 1; `enumerate` of five refuted predictors in Section 4 | both converted to connected prose |
| Every display item integrated into the text | five of nine figures never referenced (`fig:bp`, `fig:grid`, `fig:perframe`, `fig:scaling`, `fig:toierr`) | every figure and every included table referenced and discussed; three captions that had been carrying their own discussion were trimmed once the text took it over |
| Abbreviations expanded at first use | SIMD, ALU, DRAM and UNC used without expansion | each expanded at its first occurrence |
| Average sentence length 15-20 words | 25.0 overall; abstract 35.9, conclusion 31.5, availability 30.6, broad phase 30.3 | 18.5 overall; ten of twelve sections inside the band, the widest being the abstract at 22.7 and the broad phase at 21.3 |
| Vary construction; no single phrasing carrying the argument | "rather than" 55 times in 12,100 words, one every 220 words | 3 times |
| One idea per paragraph, topic sentence first | long multi-clause paragraphs in Sections 5, 6 and 9 | split, each opening on its claim |

**Not satisfied, with the reason.**

*IMRAD structure.* The skill's default is Introduction, Methods, Results and
Discussion. This is a numerical-methods paper and it keeps the structure that
family uses: problem statement, algorithms, parallel implementation, a
correctness argument, then results and limitations. The skill lists theoretical
and methods papers as admitting discipline-specific structures, and this is one.
Forcing IMRAD would also put the manuscript at odds with its target venue.

*Average sentence length.* At 18.5 words the manuscript now sits inside the
skill's 15-20 band overall, but two sections remain above it. Sentences in
Section 4 carry algorithm names, quantities with units and cross-references, and
splitting them further would separate clauses that have to be read together.

*Display-item density.* At 22 tables and figures over 11,500 words the paper
runs at one item per 525 words, against the skill's guideline of one per 1000.
Every item is now referenced and discussed, so the density reflects an
evaluation over six scenes, two processors and two modes, and not items left
undiscussed. Three tables in Section 8.4 do present one comparison three ways,
as absolute time, as a rate and as a per-step cost. Each is kept because each
answers a different question a reader brings, and the text now says which.

*The list in Theorem 7.1.* Hypotheses (i)-(iv) remain an `enumerate`. They are
mathematical antecedents, not prose recast as bullets, and running them together
would make the theorem harder to check.
