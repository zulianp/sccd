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
| `SCCD_NARROWPHASE_MODE` | `0`, `1`, `2`, `3` | `0` | Narrow-phase kernel: `0` scalar, `1` fast-vector, `2` TightInclusion-exact, `3` TightInclusion-compat. See `src/sccd_narrowphase_mode.hpp`. **Ignored for quads**, which have one root-finder variant. |
| `SCCD_BROADPHASE` | `sweep`, `cell2d` | auto | Forces a broad phase instead of letting `choose_broadphase_strategy` decide. Both produce identical pair sets, so this only changes speed. |
| `SCCD_ADAPTIVE_SPLIT` | `0`, `1` | `1` | `0` splits each interval uniformly; `1` places splitters from a Gauss–Newton step. The assessment measured adaptive ahead on all three scenes, so `0` exists to reproduce that comparison. |
| `SCCD_USE_VNARROW_PHASE` | `0`, `1` | build option | Older, narrower switch for the vectorised vertex-face kernel. `SCCD_NARROWPHASE_MODE` supersedes it. |
| `SCCD_VNARROWPHASE_TI_COMPAT` | `0`, `1` | build option | Corrects vectorised vertex-face output against TightInclusion. Requires a TightInclusion build. |
| `SCCD_USE_TI` | `0`, `1` | `0` | Calls TightInclusion directly. Requires a TightInclusion build. Oracle use only. |

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
