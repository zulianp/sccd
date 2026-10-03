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
- [ ] ¶2  the three design pieces: cell list, SIMD packing, device queue (9) (we are using a custom variant of cell-list)
- [ ] ¶3  headline results, including the configuration where SCCD loses (7) ·N (we emphasize the accuracy gain range with at most 7x price on performance, the loss is not really relevant in the abstract)

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

## 05 · Narrow phase — narrowphase.tex (189)

- [ ] §  Narrow phase (2) ()
  - [ ] ¶1  branch and bound over parameter boxes; names the skeleton (6) ()
  - [ ] ⟐ alg:np  branch and bound for one query (32) ()
  - [ ] ¶2  acceptance reports a lower bound in t; every non-rejection accepts (3) ()
  - [ ] §§ Two acceptance tests (2) ()
    - [ ] ¶3  defines Tight: TI's test transcribed, with its split rule (13) ()
    - [ ] ≡ eq:accept-tight  the Tight acceptance test (—) ()
    - [ ] ¶4  defines Relaxed; both refuse t=0 contacts at no measured cost (9) ()
  - [ ] §§ The adaptive splitter, and what measuring it showed (2) ()
    - [ ] ¶5  the clamp does the work, not Newton; the principled rival loses (20) ·N ()
  - [ ] §§ Order independence (2) ()
    - [ ] ¶6  TI needs a priority queue; this search needs none (4) ()
    - [ ] ⊢ prop:order  order-independence of the accepted minimum (9) ()
    - [ ] ⊢ proof  proof of prop:order (11) ()
    - [ ] ⊢ rem:ti  TI's answer is the same minimum, so the two agree exactly (10) ()
    - [ ] ¶7  the hypotheses bind: host matches TI bit for bit, the device does not (14) ()
    - [ ] ¶8  the freedom granted: any traversal, lanes and blocks mixed across queries (5) ()
  - [ ] §§ Quads (2) ()
    - [ ] ¶9  a separate quad kernel, unevaluated; lexicographic order required (18) ()
    - [ ] ¶10  hoisting the split-plane evaluation is faster per query, slower overall (7) ·N ()

## 06 · Parallel implementation — parallel.tex (294)

- [ ] §  Parallel implementation (2) ()
  - [ ] ¶1  both kernels answer heavy-tailed query cost, in opposite directions (16) ()
  - [ ] ▦ tab:difficulty  boxes per query on cloth-funnel, host against device; 7 rows (25) ·HAND ·N ()
  - [ ] §§ Host: packing lanes across queries (2) ()
    - [ ] ⟐ alg:host  the host SIMD narrow phase, three phases (53) ()
    - [ ] ¶2  the box, not the query, is the unit of work across lanes (7) ·N ()
    - [ ] ¶3  storing eight corner values gives a fourfold reduction by inheritance (5) ·N ()
    - [ ] ¶4  replicating geometry per lane avoids gathers (4) ·N ()
    - [ ] ¶5  classify, evaluate, assemble; the lower child stays for depth-first in t (9) ·N ()
    - [ ] ¶6  query types peak at different lane widths; cross-process timing caveat (13) ·N ()
  - [ ] §§ Device: a small stack, and a queue that redistributes (2) ()
    - [ ] ◧ fig:gpu-workflow  shared stack spilling into a double-buffered queue (11) ·TIKZ ()
    - [ ] ⟐ alg:device  the device kernel and its relaunching host loop (49) ·N ()
    - [ ] ¶7  one thread per query; the obvious alternative fails in three ways (3) ()
    - [ ] ¶⟨The global queue⟩  double buffering avoids a handshake; one buffer hangs (5) ()
    - [ ] ¶⟨The shared stack⟩  a small stack spills heavy queries out for redistribution (9) ·N ()
    - [ ] ¶⟨Overflow handling⟩  dropping overflow is safe; bounded retry after growing (7) ·N ()
    - [ ] ¶8  two measured decisions; the device broad phase is a conventional port (11) ·N ()
  - [ ] §§ Precision (2) ()
    - [ ] ¶9  double internally is a requirement; single violates the definition (6) ·N ()
    - [ ] ¶10  templated storage; narrowing the output rounds toward minus infinity (6) ()
  - [ ] §§ Why the work divides the way it does (2) ()
    - [ ] ¶11  first cause: the device re-evaluates all eight corners per child (7) ·N ()
    - [ ] ¶12  second cause, offered as a candidate: bound sharpness (10) ()
    - [ ] ¶13  the broad phase inverts both arguments; the GPU wins it on all six (4) ·N ()

## 07 · Conservativeness by construction — guarantee.tex (119)

- [ ] §  Conservativeness by construction (2) ()
  - [ ] ¶1  measurement is not proof; promises the structural argument (5) ·N ()
  - [ ] ⊢ thm:conservative  four conditions imply the reported ToI is never late (18) ()
  - [ ] ⊢ proof  proof of thm:conservative (18) ()
  - [ ] ¶2  announces three consequences (1) ()
  - [ ] ¶⟨Accepting is always safe⟩  a looser accept costs accuracy only (9) ()
  - [ ] ¶⟨A cap is an accuracy limit⟩  exhaustion accepts, so a small cap cannot miss (4) ()
  - [ ] ¶3  an undersized work list costs four orders of accuracy; shows three values (13) ·N ()
  - [ ] ¶4  all three values land early; misdiagnosis sends the reader elsewhere (5) ·N ()
  - [ ] ¶5  derive capacity from the depth limit; headroom is nearly free (6) ·N ()
  - [ ] ¶⟨Only an unsound rejection loses a root⟩  condition (i) is the whole exposure (4) ()
  - [ ] ¶6  padding by the distance tolerance or by machine epsilon (8) ·N ()
  - [ ] ¶7  a third failure: a depth cap that drops, visible only in double (6) ()
  - [ ] ¶8  reduces conservativeness to a three-question checklist (8) ()

## 08 · Numerical experiments — results.tex (437)

- [ ] §  Numerical experiments (2) ()
  - [ ] ¶1  three properties evaluated; the losses are included (5) ()
  - [ ] §§ Experimental setup (2) ()
    - [ ] ¶2  one GH200 module per measurement; threads bound with numactl (8) ·N ()
    - [ ] ▦ tab:hardware  host and device specifications of one module; 11 rows (27) ·HAND ·N ()
    - [ ] ¶3  derives the 8.5x/8.7x ratio speedups should be read against (9) ·N ()
  - [ ] §§ Software (2) ()
    - [ ] ¶4  build environment and the default search parameters (6) ·N ()
    - [ ] ¶5  competitors used as released, same flags, same affinity (7) ()
  - [ ] §§ Timing protocol and dataset (2) ()
    - [ ] ¶6  defines the four timing tags (12) ()
    - [ ] ¶7  device clocks synchronised; warm-up; repeats in separate processes (7) ·N ()
    - [ ] ¶8  between- against within-allocation variance; the significance floor (7) ·N ()
    - [ ] ¶9  describes the dataset and defines a case (7) ·N ()
    - [ ] ▦ tab:dataset  scene sizes: cases, pairs per step, curated queries; 6 rows (1) ()
  - [ ] §§ Conservativeness (2) ()
    - [ ] ¶10  the central claim: nothing missed, nothing late (7) ·N ()
    - [ ] ▦ tab:conservativeness  signed per-scene conservativeness; 6 rows (1) ()
    - [ ] ¶11  why compare against exact roots, not TI; reports false positives (7) ·N ()
  - [ ] §§ Against TightInclusion (2) ()
    - [ ] ¶12  protocol: narrow phases only, over the same query sets (8) ()
    - [ ] ◧ fig:reference  SCCD against TI, Scalable CCD and additive CCD per scene (17) ()
    - [ ] ¶13  speed against TI per scene; rod-twist extreme, cloth-funnel the exception (7) ·N ()
    - [ ] ¶14  defines earliness; the host matches TI on all twelve scene-phases (9) ·N ()
    - [ ] ¶15  the device agrees closely, not exactly (4) ·N ()
    - [ ] ¶16  explains the device gap as a different float operation sequence (11) ()
    - [ ] ¶17  worst cases are intrinsic: grazing and already-touching configurations (5) ·N ()
  - [ ] §§ CPU against GPU (2) ()
    - [ ] ▦ tab:processor  per-scene totals and three phase ratios; 6 rows (1) ()
    - [ ] ¶18  per-phase ratios; the margin rests on the narrow phase (8) ·N ()
    - [ ] ¶19  the end-to-end winner per scene, including the two the host takes (12) ·N ()
    - [ ] ¶20  supplement qualifiers: a standing gap, and drift through rod-twist (6) ·N ()
    - [ ] ¶21  phase shares; broad phase tight around its median, narrow phase skewed (10) ·N ()
    - [ ] ◧ fig:breakdown  phase shares per processor as stacked bars (9) ()
    - [ ] ◧ fig:toierr  earliness distributions, SCCD above, additive CCD below (14) ·N ()
  - [ ] §§ Against other implementations (2) ()
    - [ ] ¶22  each competitor compared on the question it answers (8) ·N ()
    - [ ] ¶23  the device beats Scalable CCD end to end on every scene (8) ·N ()
    - [ ] ¶24  the broad-phase margin is standing; separates algorithm from implementation (10) ·N ()
    - [ ] ◧ fig:bpvs  device broad phase against Scalable CCD's, step by step (11) ()
    - [ ] ¶25  Scalable CCD reports late times of impact; traced to a buffer race (9) ·N ()
    - [ ] ¶26  additive CCD is cheaper per candidate; the margin closes with size (7) ·N ()
    - [ ] ¶27  that speed costs the answer: earliness and false positives (10) ·N ()
  - [ ] §§ Scaling with element count (2) ()
    - [ ] ¶28  refinement study design: only the processor varies, no contact (5) ·N ()
    - [ ] ◧ fig:scaling  cost against element count, log-log, with fitted exponents (8) ()
    - [ ] ¶29  the crossover point and the fitted exponents (9) ·N ()
    - [ ] ¶30  structure against queries splits differently per processor (8) ·N ()
    - [ ] ¶31  the narrow-phase column prices the rejection test alone (8) ·N ()
  - [ ] §§ Scaling with thread count (2) ()
    - [ ] ¶32  the structure build is the Amdahl bottleneck at 72 cores (13) ·N ()
  - [ ] §§ Trading accuracy for speed (2) ()
    - [ ] ¶33  Relaxed is safe and looser; the false-positive cost is substantial (11) ·N ()
    - [ ] ¶34  the accuracy cost is unconditional, the speed gain conditional (5) ()

## 09 · Limitations and threats to validity — limitations.tex (102)

- [ ] §  Limitations and threats to validity (2) ()
  - [ ] ¶1  announces three classes of limitation (6) ()
  - [ ] §§ What the ground truth reaches (1) ·NOLABEL ()
    - [ ] ¶2  ground truth covers the query sets; the mesh path is a tripwire (7) ()
    - [ ] ¶3  the mesh container is single precision; bounds any accuracy claim (6) ·N ()
    - [ ] ¶4  inflicted-rounding experiment and control locate the fault in the input (10) ·N ()
  - [ ] §§ What a measurement establishes (1) ·NOLABEL ()
    - [ ] ¶5  shared-minimum pruning makes the last digits non-deterministic (6) ·N ()
    - [ ] ¶6  discloses a device Relaxed defect found and fixed during evaluation (10) ·N ()
    - [ ] ¶7  what a check can establish: signed comparison over every query (4) ·N ()
    - [ ] ¶8  driver overhead excluded; the harness is no model of in-solver cost (5) ·N ()
  - [ ] §§ What was not measured (1) ·NOLABEL ()
    - [ ] ¶9  announces three unmeasured parts (2) ()
    - [ ] ¶10  intra-cell ordering is host-only; open on the device (4) ()
    - [ ] ¶11  the quad path is evaluated nowhere; its constant uncertified (5) ()
    - [ ] ¶12  speedups are over a host reference, not prior parallel work (5) ()
    - [ ] ¶13  counter measurements support ratios only, not absolute time (6) ()
    - [ ] ¶14  no minimum separation, no multi-GPU, no in-solver measurement (4) ()

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
- Three tables are hand-typed with no generator and no `Source:` line:
  `tab:hardware`, `tab:difficulty`, and `tab:stackcap` in the supplement.
  `tab:difficulty`'s own caption flags it as pre-redesign and as aggregating two
  device code paths.
- 18 of 30 generated tables and 2 figures are committed but included nowhere. The
  abstract's `26.6x` is backed only by `tab-reference`, an orphan.
  `Trading accuracy for speed` argues "each mode is therefore reported in its own
  table" and inputs no table at all; the whole `-relaxed` family is orphaned.
- `results.tex` is 437 lines, a quarter of the article's prose, and holds 34 of its
  121 body paragraphs.
