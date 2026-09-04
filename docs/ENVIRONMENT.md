# Environment variables

Every environment variable SCCD reads. Several of these were read but documented
nowhere, which is how `SCCD_ADAPTIVE_SPLIT` came to be the only way to reach a
whole splitting strategy and how an unprefixed lowercase `alpha` survived in a
library.

They are all diagnostics or overrides. **None of them is needed to use SCCD**,
and none of them changes an answer from wrong to right — the defaults are the
supported configuration.

## Choosing an implementation

| Variable | Values | Default | Effect |
|---|---|---|---|
| `SCCD_NARROWPHASE_MODE` | `0`, `1`, `2`, `3` | `0` | Narrow-phase kernel: `0` `Fast`, `1` `FastVectorized`, `2` `Tight`, `3` `TightInclusionCorrected`. Modes `1` and `3` are **validation-only** and fall back to `0` unless the build has TightInclusion — see below. **Ignored for quads**, which have one root-finder variant. |
| `SCCD_BROADPHASE` | `sweep`, `cell2d` | auto | Forces a broad phase instead of letting `choose_broadphase_strategy` decide. Both produce identical pair sets, so this only changes speed. |
| `SCCD_USE_VNARROW_PHASE` | `0`, `1` | build option | Older, narrower switch for the vectorised vertex-face kernel. `SCCD_NARROWPHASE_MODE` supersedes it. |
| `SCCD_VNARROWPHASE_TI_COMPAT` | `0`, `1` | build option | Corrects vectorised vertex-face output against TightInclusion. Requires a TightInclusion build. |
| `SCCD_USE_TI` | `0`, `1` | `0` | Calls TightInclusion directly. Requires a TightInclusion build. Oracle use only. |

### Why modes 1 and 3 need TightInclusion

Both enter the same vectorised kernel. Mode 3 then corrects its answers with the
external library, which makes it an oracle rather than a code path; mode 1 is a
duplicate that lost every scene the assessment measured, and is vertex-face only.
Asking for either without the library gets mode 0 rather than an error.

### The three cases that used to be silent

Falling back is fine; doing it without saying so is not, because a caller who
measures a kernel they did not ask for draws a conclusion about the wrong thing.
Each of these now prints one line to stderr, once per process, and changes
nothing else:

| What you set | What runs | What it says |
|---|---|---|
| `SCCD_NARROWPHASE_MODE=1` or `3`, no TightInclusion | mode 0 | that the mode is validation-only and needs `SCCD_ENABLE_TIGHT_INCLUSION=ON` |
| `SCCD_NARROWPHASE_MODE=20`, or any non-number | the build's default | that the value was ignored, and what the valid ones are |
| `SCCD_NARROWPHASE_MODE=1`, `2` or `3` on a quad mesh | the one quad kernel | that the mode does not reach the quad path |

An unset variable, or `SCCD_NARROWPHASE_MODE=0`, says nothing: there is no
surprise to report. The mode is re-read from the environment on every call, so a
caller that switches it between calls keeps working -- the warnings are one-shot
flags rather than a resolved-once decision.

**Removed variables.** `SCCD_ADAPTIVE_SPLIT` (uniform splitting, demoted to a
spike) and `SCCD_GSTACK_CAP_INIT` (never read). Setting either does nothing.

## Search parameters

| Variable | Type | Default | Effect |
|---|---|---|---|
| `SCCD_MAX_DEPTH` | int | `69` | Maximum subdivision depth. At the cap a box is **accepted** at its `t` lower bound, never dropped — that is what keeps a depth limit from costing a collision. |
| `SCCD_TOL` | float | `3e-8` | Codomain tolerance for the acceptance test. |
| `SCCD_REFINE` | int | `0` | Extra refinement passes in the vertex-quad search. |

Lowering `SCCD_TOL` or raising `SCCD_MAX_DEPTH` tightens the reported time of
impact and costs time. Neither can make the result unsafe: the rejection test
pads by the certified numerical error bound, which these do not touch.

## CUDA tuning

| Variable | Type | Default | Effect |
|---|---|---|---|
| `SCCD_BLOCKS_PER_SM` | int | occupancy API | Blocks resident per SM. |
| `SCCD_BATCH_SIZE` | int | all candidates | Candidates per outer iteration. |
| `SCCD_GSTACK_CAP_MAX` | int | `INT_MAX` | Soft cap on a single growth step of the global stack. The stack starts empty and grows from the deficit the kernel reports. |
| `SCCD_NP_ALPHA` | float | `0.5` | Splitter blending factor. Read as a bare lowercase `alpha` until it was renamed — an unprefixed environment variable in a library's process. |

There is no `SCCD_GSTACK_CAP_INIT`. It was documented in the source for a while
and never read anywhere.

## Diagnostics

`SCCD_NP_COUNT_BOXES` is a **compile-time** define, not an environment variable:
`-DSCCD_NP_COUNT_BOXES` in `CMAKE_CUDA_FLAGS` and `CMAKE_CXX_FLAGS` makes the
device narrow phase and the host TightInclusion kernel each count the boxes they
classify and print one line per call to stderr. Both count the same unit, so the
two numbers are directly comparable — which is what showed that the device
conservative search does 535× the host's work on identical queries. It is a
global atomic on the hot path, so an instrumented build's *timings* mean nothing;
only the counts do.

| Variable | Type | Effect |
|---|---|---|
| `SCCD_BROADPHASE_VERBOSE` | set | Logs the chosen broad phase and the shape statistics behind the choice. |
| `SCCD_MISSING_PAIRS_CSV` | path | `sccd_bench` writes broad-phase pairs the reference found and it did not. |
| `SCCD_BENCH_EXECUTION_SPACE` | `host`, `device` | Where the benchmark drivers run. |
| `SCCD_TOPOLOGY` | `quad` | Makes `refine_scaling` generate a hexahedral cube, so its skin is `QUADSHELL4` rather than `TRISHELL3`. |
| `SCCD_SCALE` | float | Scale factor for `refine_scaling`'s synthesised motion. Negative reflects the mesh through its centre, which is what makes a convex mesh self-intersect. |
| `SCCD_BENCH_MAX_CASES` | int | Evenly spread subsample of a scene's cases, for variant sweeps that must run every variant on the same cases. |
| `SCCD_DB_TO_RAW` | path | smesh's `db_to_raw`, for benchmark mesh conversion. |

## A note on reading these

`SCCD_READ_ENV(name, conversion)` stringifies the variable's own name, so the C++
identifier *is* the environment variable. That is convenient and it is also the
trap that produced `getenv("alpha")`: the identifier has to be named as the
environment variable should be, not as the code would like to read.
