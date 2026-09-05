# Environment variables

Every environment variable SCCD reads. They are all diagnostics or overrides. **None of them is needed to use SCCD**,
and none of them changes an answer from wrong to right — the defaults are the
supported configuration.

## Choosing an implementation

| Variable | Values | Default | Effect |
|---|---|---|---|
| `SCCD_NARROWPHASE_MODE` | `0`, `2` | `0` | Narrow-phase kernel: `0` `Relaxed`, `2` `Tight`. Any other value warns and runs `0`. **Ignored for quads**, which have one root-finder variant. |
| `SCCD_BROADPHASE` | `sweep`, `cell2d` | auto | Forces a broad phase instead of letting `BroadPhaseAutoTuner` race the two and keep the winner. Both produce identical pair sets, so this only changes speed. |
| `SCCD_USE_TI` | `0`, `1` | `0` | Calls TightInclusion directly. Requires a TightInclusion build. Oracle use only. |

### Getting TightInclusion's answer

`SCCD_USE_TI=1` dispatches every query straight to the external library, for both
vertex-face and edge-edge. Needs `SCCD_ENABLE_TIGHT_INCLUSION=ON`. This is the
only way to get the reference's answer; SCCD ships no kernel that is part its own
search and part the reference.

### When SCCD does not run what you asked for

Falling back is fine; doing it without saying so is not, because a caller who
measures a kernel they did not ask for draws a conclusion about the wrong thing.
Each of these prints one line to stderr, once per process, and changes nothing
else:

| What you set | What runs | What it says |
|---|---|---|
| `SCCD_NARROWPHASE_MODE=1` or `3` | `Relaxed` | that the value selects no kernel, and to use `SCCD_USE_TI=1` for the reference |
| `SCCD_NARROWPHASE_MODE=20`, or any non-number | `Relaxed` | that the value was ignored, and what the valid ones are |
| any mode on a quad mesh | the one quad kernel | that the mode does not reach the quad path |

An unset variable, or `SCCD_NARROWPHASE_MODE=0`, says nothing: there is no
surprise to report. The mode is re-read from the environment on every call, so a
caller that switches it between calls keeps working -- the warnings are one-shot
flags rather than a resolved-once decision.

**Not read.** `SCCD_ADAPTIVE_SPLIT` and `SCCD_GSTACK_CAP_INIT` are not read
anywhere; setting them does nothing.

## Search parameters

**The library takes these as parameters, not from the environment.**
`narrow_phase_vf`, `narrow_phase_ee` and `narrow_phase_vq` all require
`max_depth` and `tol` as arguments; nothing inside the library calls `getenv` for
them. `SCCD_MAX_DEPTH` and `SCCD_TOL` are read by the *drivers* -- `sccd_bench`,
`sccd_refine_scaling` and the mesh demos -- which pass the values in. Setting
them changes what those programs ask for; it does not reach into a library call
made by your own code.

| Variable | Type | Default | Effect |
|---|---|---|---|
| `SCCD_MAX_DEPTH` | int | `69` | Maximum subdivision depth. At the cap a box is **accepted** at its `t` lower bound, never dropped — that is what keeps a depth limit from costing a collision. The device vertex-quad kernel holds a stack for depth 128 (`-DSCCD_VQ_MAX_DEPTH`) and reports on stderr if asked for more. |
| `SCCD_TOL` | float | `3e-8` | Codomain tolerance for the acceptance test. |
| `SCCD_REFINE` | int | `0` | Newton polish on an accepted vertex-face box: if it converges to an earlier time inside the box, that time is reported instead. Vertex-face only — the edge-edge and vertex-quad searches have no polish step and ignore it. |

Lowering the tolerance or raising the depth tightens the reported time of impact
and costs time. Neither can make the result unsafe: the rejection test pads by
the certified numerical error bound, which neither touches.

`SCCD_REFINE` is the one of the three the library does read, in
`narrow_phase_vf` and `narrow_phase_ee`.

## CUDA tuning

| Variable | Type | Default | Effect |
|---|---|---|---|
| `SCCD_BLOCKS_PER_SM` | int | occupancy API | Blocks resident per SM. |
| `SCCD_BATCH_SIZE` | int | all candidates | Candidates per outer iteration. |
| `SCCD_NP_S1_BLOCK_PER_QUERY` | `0`, `1` | `0` | Runs the `ToiOutput::PerPair` path as one block per query instead of one thread per query. The default is one thread, which measures 2.8–4.5× faster for identical results; this is here to reproduce the comparison. |
| `SCCD_GSTACK_CAP_MAX` | int | `INT_MAX` | Soft cap on a single growth step of the global stack. The stack starts empty and grows from the deficit the kernel reports. |
| `SCCD_NP_ALPHA` | float | `0.5` | Splitter blending factor. |

## Diagnostics

`SCCD_NP_COUNT_BOXES` is a **compile-time** define, not an environment variable:
`-DSCCD_NP_COUNT_BOXES` in `CMAKE_CUDA_FLAGS` and `CMAKE_CXX_FLAGS` makes the
device narrow phase and the host TightInclusion kernel each count the boxes they
classify and print one line per call to stderr. Both count the same unit, so the
two numbers are directly comparable: on cloth-funnel at `ToiOutput::Earliest` the device
classifies 94× the boxes the host does for the same queries. It is a
global atomic on the hot path, so an instrumented build's *timings* mean nothing;
only the counts do.

| Variable | Type | Effect |
|---|---|---|
| `SCCD_BROADPHASE_VERBOSE` | set | Logs the chosen broad phase and the shape statistics behind the choice. |
| `SCCD_NP_WORST_CSV` | path | The device narrow phase writes its costliest queries here, in the query-CSV format `sccd_np_trace` reads, so the same query can be run on both machines and diffed box by box. Needs a `SCCD_NP_COUNT_BOXES` build. |
| `SCCD_NP_WORST_N` | int (`4`) | How many queries `SCCD_NP_WORST_CSV` writes. |
| `SCCD_NP_NO_GLOBAL_SEED` | `0`, `1` | Seeds every host block at `max_toi` instead of the running global minimum, which removes the host's pruning advantage when comparing its box count against the device's. Needs a `SCCD_NP_COUNT_BOXES` build; it does not exist otherwise. |
| `SCCD_MISSING_PAIRS_CSV` | path | `sccd_bench` writes broad-phase pairs the reference found and it did not. |
| `SCCD_BENCH_EXECUTION_SPACE` | `host`, `device` | Where the benchmark drivers run. |
| `SCCD_TOPOLOGY` | `quad` | Makes `refine_scaling` generate a hexahedral cube, so its skin is `QUADSHELL4` rather than `TRISHELL3`. |
| `SCCD_SCALE` | float | Scale factor for `refine_scaling`'s synthesised motion. Negative reflects the mesh through its centre, which is what makes a convex mesh self-intersect. |
| `SCCD_BENCH_MAX_CASES` | int | Evenly spread subsample of a scene's cases, for variant sweeps that must run every variant on the same cases. |
| `SCCD_DB_TO_RAW` | path | smesh's `db_to_raw`, for benchmark mesh conversion. |
| `SCCD_EXECUTION_SPACE` | `host`, `device` | Fallback `sccd_bench` reads when `SCCD_BENCH_EXECUTION_SPACE` is unset. Prefer the `BENCH` form. |
| `SCCD_USE_FIND_EARLIEST_IMPACT_TIME` | `0`, `1` | Demo only (`mesh_sccd`, `mesh_sccd_cuda`): picks `find_earliest_impact_time` over `find_impact_times`. |
| `SCCD_EXPORT_COLLISIONS` | `0`, `1` | Demo only (`mesh_sccd_cuda`): writes the colliding pairs it found. |

## A note on reading these

`SCCD_READ_ENV(name, conversion)` stringifies the variable's own name, so the C++
identifier *is* the environment variable. Anything read this way must therefore be
named as the environment variable should be, `SCCD_`-prefixed and all, not as the
surrounding code would prefer to read.
