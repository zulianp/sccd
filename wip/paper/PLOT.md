# PLOT — what every unit of the paper is for

A map of the article, one line per unit, nested by document level. Tick `[x]` what
should stay as it is; for everything else say in the trailing `()` how it should
change. A later pass revises the paper against this file.

Main article only. `supplemental_material.tex` is not covered here.

**Current shape.** 33 pages, 14 section files, 1,850 lines of prose.
13 sections · 21 subsections · 2 subsubsections · 6 named paragraphs ·
121 body paragraphs · 5 algorithms · 8 figures · 5 tables · 6 theorem-family
environments · 7 numbered equations.

**Markers.** `§ §§ §§§` section levels · `¶n` body paragraph, numbered within its
section · `¶⟨name⟩` a named `\paragraph` · `◧` figure · `▦` table · `⟐` algorithm ·
`⊢` theorem, proposition, definition, remark, proof · `≡` numbered equation.
Trailing tags: `·N` states numeric results · `·HAND` no generator ·
`·TIKZ` hand-drawn · `·NOLABEL` cannot be `\cref`-ed.

---

## 01 · Abstract — abstract.tex (28)

- [x] ¶1  frames CCD; conservativeness is asymmetric; names the assessment scale (8) ·N ()
- [x] ¶2  the three design pieces: cell list, SIMD packing, device queue (9) (we are using a custom variant of cell-list)
- [x] ¶3  headline results, including the configuration where SCCD loses (7) ·N (we emphasize the accuracy gain range with at most 7x price on performance, the loss is not really relevant in the abstract)

## 02 · Introduction — intro.tex + related.tex + contribution.tex (173)

- [x] §  Introduction (2) ()
  - [x] ¶1  defines the CCD problem; the two error modes differ in kind (6) ()
  - [x] ¶2  a late or missed time of impact breaks barrier-potential solvers (9) ()
  - [x] ⊢ def:conservative  formal definition of conservative (6) ()
  - [x] ¶3  the spec is stronger than Boolean correctness; accuracy is signed (7) ()
  - [x] ¶4  positions TightInclusion as the certified reference, used unmodified (7) ()
  - [x] §§ Prior work — related.tex (2) ()
    - [x] ¶1  broad/narrow split; separate literatures and failure modes (5) ()
    - [x] §§§ Broad phase (1) ·NOLABEL ()
      - [x] ¶2  taxonomy: sweep-and-prune, grids and hashes, BVH/KD/interval trees (8) ()
      - [x] ¶3  cell size is the grid's weak point; our two fixes remove it (13) ()
      - [x] ¶4  prior broad phases miss collisions through silent buffer overflow (5) ·N (Belgrod et al. cite it only for the dataset, put the missed collision stuff in a macro (TODO for the original author, also fix it and create a patch for their repo))
    - [x] §§§ Narrow phase (1) ·NOLABEL ()
      - [x] ¶5  surveys cubic solvers, conservative advancement, root parity, symbolic (10) ()
      - [x] ¶6  places the work in the inclusion family; TI is the first sound one (11) ()
      - [x] ¶7  inherits TI's mathematics, differs only in realising the search (8) ()
  - [x] §§ Contributions and structure of the paper — contribution.tex (2) ()
    - [x] ¶1  conservativeness is structural; the dataset carries exact roots (5) ()
    - [x] ¶2  the five claimed contributions, as an enumerate (31) ·N ()
    - [x] ¶3  declares the configurations where the pipeline loses (6) ·N ()
    - [x] ¶4  roadmap of the remaining sections (10) ()

## 03 · Problem statement and the inclusion function — prelim.tex (154)

- [x] §  Problem statement and the inclusion function (2) ()
  - [x] ¶1  multilinearity gives an exact range; float perturbation is the risk (9) ()
  - [x] §§ The queries and their difference function (2) ()
    - [x] ¶2  affine motion reduces meshes to two elementary queries (6) ()
    - [x] ¶3  vertex-face and edge-edge difference functions; defines the true ToI (19) ()
    - [x] ≡ eq:fvf, eq:fee  the two difference functions (—) ()
    - [x] ¶4  adds the bilinear-quad query; node order is lexicographic (12) ()
    - [x] ≡ eq:fvq  the quad difference function (—) ()
    - [x] ¶5  multilinearity makes the corner hull the exact range, not an enclosure (17) ()
    - [x] ≡ eq:hull  the corner hull (—) ()
  - [x] §§ What floating point costs: the error bound and the tolerances (2) (Change it to: Error bound and the tolerances)
    - [x] ¶6  rounding of corner values is the sole exposure to a lost root (5) ()
    - [x] ¶7  states the certified forward error bound and its constants (18) ·N ()
    - [x] ≡ eq:err  the error bound (—) ()
    - [ ] ¶8  delta clamped from above by one; clamping below would explode the pad (4) ·N ()
    - [x] ¶9  the bound holds only for TI's expression order; no fast math (6) (Rewrite "refuses to build under -ffast-math" to something more like "no -ffast-math compilation flags are used here", Who is refusing is unclear)
    - [x] ¶10  the quad constant is uncertified and deliberately over-estimated (4) ·N ()
    - [x] ¶11  per-axis domain tolerances from a codomain distance tolerance (15) ()
    - [x] ≡ eq:tol  the tolerance derivation (—) ()
    - [x] ¶12  tolerance quotients must be capped; the caps match the reference (8) ·N ()
    - [x] ¶13  what the family can promise: conservativeness, not exactness (5) ()
    - [x] ¶14  transition: the cheap filter first, then the search (4) ()

## 04 · Broad phase — broadphase.tex (261)

Restructured: the shared construction now has its own subsection, both queries
read it, `fig:cell2d` is gone and `alg:cell` builds the structure instead of
answering one query.

- [x] §  Broad phase (2) ()
  - [x] ¶1  must be conservative and fast; one grid serves two query shapes (9) ()
  - [x] §§ Swept boxes, and the shape of the output (2) ()
    - [x] ¶2  exact swept boxes, ULP outward rounding, count-then-fill, 32-lane mask (20) ·N ()
  - [x] §§ The cell list ·NEW (2) ()
    - [x] ¶1  two axes, mean-extent spacing, 4N cap, no user voxel size (17) ·N ()
    - [x] ¶2  the swept box sets the spacing; largest-box and two-level both lose (20) ·N ()
    - [x] ¶3  one index per box at its minimum corner; two counted allocations (7) ()
    - [x] ¶4  the two cumulative maxima, and what they license (9) ()
    - [x] ⟐ alg:cell  building the cell list, shared by both queries (29) ()
    - [x] ¶5  a footprint walk is incomplete; the row-major order is total (14) ·N ()
    - [x] ¶6  nothing sorted; cost terms with one binary search per row (7) ()
  - [x] §§ The vertex-face query (2) ()
    - [x] ¶1  faces binned, vertices walk; the range and what bounds each side (12) ()
    - [x] ¶2  why binning the faces is the cheap choice; pair set unchanged (8) ·N ()
  - [x] §§ The edge-edge query (2) ()
    - [x] ¶1  forward from its own cell; equal index decides; no row search (12) ·N ()
    - [x] ⟐ alg:cellmin  the edge-edge self-query pass (41) ()
    - [x] ◧ fig:cellmin  worked trace of the walk and its skipped spans (29) ·TIKZ ()
    - [x] ¶2  intra-cell ordering measured and rejected (6) ()

## 05 · Narrow phase — narrowphase.tex (161)

One acceptance test, not two. Relaxed and the adaptive-splitter subsection are
gone, and the `\textsc{Tight}` name with them, since there is nothing left to
distinguish it from. Carried through sections 07, 08 and 09.

- [x] §  Narrow phase (2) ()
  - [x] ¶1  branch and bound over parameter boxes; names the skeleton (6) ()
  - [x] ⟐ alg:np  branch and bound for one query (32) ()
  - [x] ¶2  acceptance reports a lower bound in t; every non-rejection accepts (3) ()
  - [x] §§ The acceptance test (2) ()
    - [x] ¶3  TI's test transcribed, with its split rule (11) ()
    - [x] ≡ eq:accept-tight  the acceptance test (—) ()
    - [x] ¶4  the test refuses a contact at exactly t=0, at no measured cost (4) ()
  - [x] §§ Order independence (2) ()
    - [x] ¶5  TI needs a priority queue; this search needs none (4) ()
    - [x] ⊢ prop:order  order-independence of the accepted minimum (9) ()
    - [x] ⊢ proof  proof of prop:order (11) ()
    - [x] ⊢ rem:ti  TI's answer is the same minimum, so the two agree exactly (10) ()
    - [x] ¶6  the hypotheses bind: host matches TI bit for bit, the device does not (14) ()
    - [x] ¶7  the freedom granted: any traversal, lanes and blocks mixed across queries (5) ()
  - [x] §§ Quads (2) ()
    - [x] ¶8  a separate quad kernel, unevaluated; its multi-way split; node order (21) ()
    - [x] ¶9  hoisting the split-plane evaluation is faster per query, slower overall (7) ·N ()

## 06 · Parallel implementation — parallel.tex (377)

Restructured and re-measured. The prose now exposes the two kernels and leaves
the discussion to the end; the constants moved into a subsection of their own;
`tab:difficulty` was re-measured on the current kernels, generated rather than
typed, and moved to the end with the discussion. Section 07 was folded into
this one and into 05, as directed.

Both algorithms were checked line by line against the shipped code. Four
corrections: `alg:host` was missing the prune and the domain test applied to each
child before it is written; `alg:device` was missing the clamp of the parent's
$t$ upper bound to the running bound and the acceptance of a child; the drain
loop was shown as bounded at 32 rounds when 32 bounds the grow-and-retry and the
drain loop runs to exhaustion; and the thread-per-query margin was stated as
`2.76x--4.43x` where the three measurements are 4.47, 2.76 and 4.43.

**`tab:difficulty` was stale and the conclusion it supported has changed.** The
old numbers were the pre-redesign diagnostic that motivated the device redesign,
on 881,694 queries. Re-measured through the pipeline on all 25,192,698 candidate
pairs of cloth-funnel, both processors, same candidates, same counters: the
device examines 4.8x as many boxes in total, but the excess is in the middle of
the distribution, and it is now the **host** that owns the tail --- six host
queries above a million boxes and a worst of 4,591,271, against none and 101,246
for the device. The discussion was rewritten accordingly: the old text claimed
the device had the stretched tail.

- [x] Update the algorithms and text to the current version of the code (the optimal only)

- [x] §  Parallel implementation (2) (Improve the prose, make it clearer and more descriptive of the method and algorithms, leave the dicussion for the end of the section)
  - [x] ¶1  the freedom of prop:order, spent in opposite directions (14) ()
  - [x] ¶2  heavy-tailed query cost is what both designs answer (10) ()
  - [x] §§ Host: lanes packed across queries (2) (More descriptive, exposition of the algorithm, less about the mirco-optimizations, which can have their own concise dedicated subsection)
    - [x] ¶3  the box, not the query, is the unit of work across lanes (14) ()
    - [x] ¶4  classify: every fate decided from stored corner values (11) ()
    - [x] ¶5  evaluate: one flat sweep over four mid-face corners (10) ()
    - [x] ¶6  assemble: the lower child stays, the upper is pushed (12) ()
    - [x] ⟐ alg:host  the host SIMD narrow phase, three phases (49) ()
  - [x] §§ Device: a small stack, and a queue that redistributes (2) (Less dictionary style, more descriptive paragraph for a scientific readership audience, clear desriptive and concise technical references)
    - [x] ¶7  one thread per query, one launch per round of work (11) ()
    - [x] ¶8  the small shared stack and why spilling early is the point (11) ·N ()
    - [x] ¶9  the double-buffered queue; a single pool deadlocks (11) ()
    - [x] ¶10  overflow is counted and dropped; why that is sound; grow and retry (14) ()
    - [x] ¶11  sizing the queue from the depth cap (6) ·N ()
    - [x] ¶12  the device broad phase is a conventional port (5) ()
    - [x] ◧ fig:gpu-workflow  redrawn: where a box lives, no prose inside (12) ·TIKZ (Make the figure look more modern, less text overall, clear and concise. )
    - [x] ⟐ alg:device  the device kernel and its relaunching host loop (46) ()
  - [x] §§ Precision (2) ()
    - [x] ¶13  double internally is a requirement; single violates the definition (6) ·N ()
    - [x] ¶14  templated storage; narrowing the output rounds toward minus infinity (6) ()
  - [x] §§ Parameters settled by measurement (2) (NEW: the micro-optimisations, in one concise subsection)
    - [x] ¶15  vector width per query type, and the cross-process timing caveat (16) ·N ()
    - [x] ¶16  shared-stack capacity: 4059.2 to 37.2 ms, identical answers (9) ·N ()
    - [x] ¶17  bound refresh every iteration (7) ·N ()
    - [x] ¶18  thread per query against block per query (6) ·N ()
    - [x] ¶19  work-list capacity is an accuracy limit; the quad example (from 07) (12) ·N ()
  - [x] §§ Why the work divides the way it does (2) (moved to the end, with the table)
    - [x] ¶20  both tails; where the weight sits on each processor (8) ·N ()
    - [x] ¶21  the device's 4.8x total: eight corners against four inherited (10) ·N ()
    - [x] ¶22  the host owns the extreme tail: the bound schedule, not the search (12) ·N ()
    - [x] ¶23  the broad phase inverts both arguments (5) ·N ()
    - [x] ▦ tab:difficulty  boxes per query on cloth-funnel, host against device (generated) (Are these numbers actual? Is the device queries still this amount?, this is not the right place it should be at the end of the section with the discussion)

## 07 · Conservativeness by construction — REMOVED

Removed as directed, and the material that earns its place distributed:

- the safety argument itself is now four sentences in 05, made by inheritance
  from TightInclusion's certified bound rather than as a theorem of our own;
- the two consequences used elsewhere (a looser accept costs accuracy only; a
  cap is an accuracy limit) follow it there, since 08 cites both;
- the two unsound paddings and the depth-cap-that-drops failure join the
  acceptance test in 05, which is where the rejection test is defined;
- the work-list capacity, its derivation from the depth cap, and the vertex-quad
  example are in 06's new parameter subsection and its device subsection.

Dropped in the move: `thm:conservative` and its proof. Nothing else cites them
any more; say the word if the formal statement should come back somewhere.

## 08 · Numerical experiments — results.tex (454)

Six changes, all as marked. The between-allocation variance paragraph is gone,
EToI is expanded where it is defined, the supplement is no longer named in the
setup or in the thread-count subsection, `fig:reference` is split in two, the
TightInclusion subsection is rebuilt around the two questions it answers, and
the strong-scaling plot is back in the article (removed from the supplement, so
the label is not defined twice).

- [x] §  Numerical experiments (2) ()
  - [x] ¶1  three properties evaluated; the losses are included (5) ()
  - [x] §§ Experimental setup (2) (Do not mention the supplemental material)
    - [x] ¶2  one GH200 module per measurement; threads bound with numactl (8) ·N ()
    - [x] ▦ tab:hardware  host and device specifications of one module; 11 rows (27) ·HAND ·N ()
    - [x] ¶3  derives the 8.5x/8.7x ratio speedups should be read against (9) ·N ()
  - [x] §§ Software (2) ()
    - [x] ¶4  build environment and the default search parameters (6) ·N ()
    - [x] ¶5  competitors used as released, same flags, same affinity (7) ()
  - [x] §§ Timing protocol and dataset (2) (EToI expand the acronym explicitly)
    - [x] ¶6  defines the four timing tags (12) ()
    - [x] ¶7  device clocks synchronised; warm-up; repeats in separate processes (7) ·N ()
    - [x] ¶8  the repeat spread, as the floor for claiming a difference (3) ·N (Remove these type of sentences: "Between job allocations this harness varies by
about 40%, so every quantity compared against another is measured inside a single allocation
and no number reported here is placed beside one from a different allocation. ...")
    - [x] ¶9  describes the dataset and defines a case (7) ·N ()
    - [x] ▦ tab:dataset  scene sizes: cases, pairs per step, curated queries; 6 rows (1) ()
  - [x] §§ Conservativeness (2) ()
    - [x] ¶10  the central claim: nothing missed, nothing late (7) ·N ()
    - [x] ▦ tab:conservativeness  signed per-scene conservativeness; 6 rows (1) ()
    - [x] ¶11  why compare against exact roots, not TI; reports false positives (7) ·N ()
  - [x] §§ Against TightInclusion (2) (Improve the explanation, target the reader)
    - [x] ¶12  protocol: narrow phases only, over the same query sets (8) ()
    - [x] ◧ fig:reference  SCCD against TightInclusion per scene (11) (Split into two figures: 1. Against TI and 2. against scalable and ACCD)
    - [x] ◧ fig:competitors  SCCD against Scalable CCD and additive CCD (11) (NEW: the second half of the split)
    - [x] ¶13  speed against TI per scene; rod-twist extreme, cloth-funnel the exception (7) ·N ()
    - [x] ¶14  defines earliness; the host matches TI on all twelve scene-phases (9) ·N ()
    - [x] ¶15  the device agrees closely, not exactly (4) ·N ()
    - [x] ¶16  explains the device gap as a different float operation sequence (11) ()
    - [x] ¶17  worst cases are intrinsic: grazing and already-touching configurations (5) ·N ()
  - [x] §§ CPU against GPU (2) ()
    - [x] ▦ tab:processor  per-scene totals and three phase ratios; 6 rows (1) ()
    - [x] ¶18  per-phase ratios; the margin rests on the narrow phase (8) ·N ()
    - [x] ¶19  the end-to-end winner per scene, including the two the host takes (12) ·N ()
    - [x] ¶20  supplement qualifiers: a standing gap, and drift through rod-twist (6) ·N ()
    - [x] ¶21  phase shares; broad phase tight around its median, narrow phase skewed (10) ·N ()
    - [x] ◧ fig:breakdown  phase shares per processor as stacked bars (9) ()
    - [x] ◧ fig:toierr  earliness distributions, SCCD above, additive CCD below (14) ·N ()
  - [x] §§ Against other implementations (2) ()
    - [x] ¶22  each competitor compared on the question it answers (8) ·N ()
    - [x] ¶23  the device beats Scalable CCD end to end on every scene (8) ·N ()
    - [x] ¶24  the broad-phase margin is standing; separates algorithm from implementation (10) ·N ()
    - [x] ◧ fig:bpvs  device broad phase against Scalable CCD's, step by step (11) ()
    - [x] ¶25  Scalable CCD reports late times of impact; traced to a buffer race (9) ·N ()
    - [x] ¶26  additive CCD is cheaper per candidate; the margin closes with size (7) ·N ()
    - [x] ¶27  that speed costs the answer: earliness and false positives (10) ·N ()
  - [x] §§ Scaling with element count (2) ()
    - [x] ¶28  refinement study design: only the processor varies, no contact (5) ·N ()
    - [x] ◧ fig:scaling  cost against element count, log-log, with fitted exponents (8) ()
    - [x] ¶29  the crossover point and the fitted exponents (9) ·N ()
    - [x] ¶30  structure against queries splits differently per processor (8) ·N ()
    - [x] ¶31  the narrow-phase column prices the rejection test alone (8) ·N ()
  - [x] §§ Scaling with thread count (2) (Make it self contained and more concise, do not cite the suplemental material, reintroduce the plot instead)
    - [x] ¶32  the structure build is the Amdahl bottleneck at 72 cores, with the plot (30) ·N ()
    - [x] ◧ fig:strong  strong scaling of the host pipeline by phase (8) (moved up from the supplement)
  - [x] §§ Trading accuracy for speed --- REMOVED with 05's second acceptance
        test, which it was the evaluation of. Two sentences elsewhere still
        described that mode without introducing it; they are gone too.

## 09 · Limitations and threats to validity — REMOVED

Removed as directed, with `sec:limitations` and its three unlabelled
subsections. Nothing else referenced it apart from the roadmap paragraph in 02,
which no longer names it.

Dropped with it: the single-precision mesh-container bound on accuracy claims,
the disclosure of the device Relaxed defect found during evaluation, the note
that the quad path is evaluated nowhere, and the statement that the speedups are
over our own host reference. The last two are still said elsewhere --- 05 says
the quad path is not evaluated, and 08 names each comparison's reference --- but
the first two are now said nowhere.

## 10 · Software and data availability — availability.tex (35)

- [ ] §  Software and data availability (2) ()
  - [ ] ¶1  BSD 3-Clause; builds with no dependencies or configuration (3) ()
  - [ ] ¶2  the two build and test commands, verbatim (6) ()
  - [ ] ¶3  opt-in CMake options; the oracle is validation, not performance (4) ()
  - [ ] ¶4  every table and figure regenerable from published data (8) ()
  - [ ] ¶5  dataset provenance; TI never modified, the harness adapted instead (7) ()

## 11 · Conclusion — conclusion.tex (33)

- [ ] §  Conclusion (2) ()
  - [ ] ¶1  restates the headline accuracy-and-speed result (7) ·N ()
  - [ ] ¶2  order-independence generalises beyond CCD; both kernels follow from it (8) ()
  - [ ] ¶3  two negative results: statistics fail to predict, principled splitter loses (7) ()
  - [ ] ¶4  three extensions; in-solver evaluation the most important (5) ()

## 12 · Acknowledgements — acknowledgements.tex (27)

- [ ] §  Acknowledgements (2) ()
  - [ ] ¶1  Swiss funding: SNSF, PASC, USI (6) ()
  - [ ] ¶2  co-authors' NSF, NSERC, KAUST and EU grants (4) ()
  - [ ] ¶3  CSCS compute; the dataset authors (5) ()
- [ ] §  Generative AI disclosure (1) ·NOLABEL ()
  - [ ] ¶4  AI used for readability and harness drafting; measurements are real (4) ()

---

## Notes the inventory turned up

- `Prior work` and `Contributions` are subsections **under** Introduction, which is
  what `STRUCTURE_TEMPLATE.md` asks for. Confirm before anything is renumbered.
- Five units carry no label and so cannot be `\cref`-ed: both subsubsections of
  `related.tex`, all three subsections of `limitations.tex`, and the
  Generative AI disclosure heading.
- ~~Three tables are hand-typed with no generator and no `Source:` line.~~ Done:
  `tab:difficulty` and `tab:hardware` are generated now, by
  `scripts/make_difficulty.py` and `scripts/make_hardware.py`, and `tab:stackcap`
  is gone --- the article states the stack-capacity result in prose, so the
  supplement no longer repeats it. No table in either document is typed by hand.
- 18 of 30 generated tables and 2 figures are committed but included nowhere.
  `tab-reference` and `tab-earliness-ref` are no longer among them: they are in
  the supplement, under a section of their own.

  The abstract's `26.6x` was the per scene-phase span while the results section
  quoted the per-scene one (`20.9x`), so the article stated two ranges for one
  comparison. Both are correct aggregations; the article now uses per-scene
  totals throughout, the supplement gives the per scene-phase form, and
  `make check` recomputes the span from `tab-reference` and fails if the
  abstract, the conclusion or the results section stops stating it.

  The `-relaxed` family of tables is orphaned for good now: the subsection that
  would have shown them went with 05's second acceptance test.
- `results.tex` is 437 lines, a quarter of the article's prose, and holds 34 of its
  121 body paragraphs.
