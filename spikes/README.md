# Spikes

Code that does not meet the bar for the shipped library, kept because it may
still be useful to read or revive.

A spike is:

- **not built by default** — everything here is behind `SCCD_ENABLE_SPIKES`,
  which is `OFF`;
- **not installed** — no spike header reaches the install tree;
- **not covered by the correctness gate** — `ctest` does not run it, and nothing
  in `wip/ASSESSMENT.md` depends on it;
- **deletable without notice** — nothing in the shipped library may include a
  spike header or call a spike symbol. That is the property that makes this
  directory worth having, and it is the one to check when adding to it.

The bar it failed is in `wip/ASSESSMENT.md`: a component ships if it
satisfies the conservativeness invariant **and** it is the only implementation of
its job. Where several compete, the best one stays and the rest come here.

## What is here, and why

| Component | Why it was demoted |
|---|---|
| `hacks/` | An out-of-tree adapter for Scalable-CCD that cannot configure: it references a `tests/` directory that does not exist, and needs scalable_ccd, Catch2 and libigl. It is also the sole consumer of the three headers below, which is why they moved with it. |
| `src/cell_broadphase.hpp` | A 3D cell list. Superseded by `src/broadphase/sccd_broadphase_cell2d.hpp`, which is two-dimensional on purpose — a surface does not fill a volume, so the third axis buys cells rather than selectivity. |
| `src/cellmin_left_bound_probe.probe.cpp` | Whether the edge-edge walk has to look left of the querying box's own column on later rows. Restricting it loses 4,488 of 26,950 pairs, which is the footprint rectangle's failure again. |
| `src/cellmin_scene_probe.probe.cpp` | The scene the paper's edge-edge figure draws, run through the shipped bounds, printing each box's cell, walk range and skipped columns so the drawing can be checked against the code. |
| `src/min_corner_self_probe.probe.cpp` | Not a demoted component but a standing measurement, and the reason two entries in `wip/DECISIONS.md` had to be withdrawn. It bins self-query boxes at their minimum corner and reads forward two ways -- over the footprint rectangle, which loses the pairs whose minimum cells are incomparable, and over the row-major linear index with the cell bounds pruning the walk, which loses nothing and tests a half to a fifth of the candidates the shipped scheme does. |
| `src/broadphase_hgrid2.hpp` | A two-level cell list: boxes split by size across a fine grid and a coarse one, each binned at its centroid so a fixed `3x3` stencil finds every partner and no duplicate test is needed. Correct -- it agreed with the sweep pair for pair on every case, with zero missed collisions and zero late times of impact on three scenes -- but 1.3x to 1.5x slower than `src/broadphase/sccd_broadphase_cell2d.hpp` on the host, and two levels is the wrong number for these size ratios. See `wip/DECISIONS.md`. |
| `src/broadphase_lb.hpp` | A lower-bound broad phase with no caller outside `hacks/`, yet installed as public API by the old flat header glob. |
| `src/sccd.hpp` | Aggregate header used only by `hacks/`. |
| `src/cuda/sccd_broadphase_warp.*` | `count_overlaps_warp` is defined and explicitly instantiated, and called by nothing repository-wide. |
| `src/cuda/sccd_lower_bound_all_to_all.*` | Reachable only from its own demo; no path from the public API. Its demo, `demo/cuda/mesh_lower_bound_cuda.exe.cpp`, is run with `scripts/mesh_sccd.sh <scene> <frame0> <frame1> --binary mesh_lower_bound_cuda`; the shell script that used to drive it lived in `scripts/` and carried hard-coded paths into one person's scratch directory. |
| `src/dead.hpp` (narrow phase) | `find_root_newton` and its helpers `norm_diff_vf` and `project_uv_simplex`, `Box::bisect_ee` (its sibling `bisect_vf` is used), `compute_edge_edge_codomain_widths`, `toi_output_name`, and the TightInclusion-guarded `barycentric_triangle_3d` and `isInsideTriangle`. The Newton step was reachable only through an undocumented `SCCD_REFINE` environment variable defaulting to off, and was unsound when set: it reported `min(upper, max(lower, 0.99 * t_approx))` rather than the box's `t` lower bound, which is the only value guaranteed to be at or before a root in the box. It was also half-wired — the edge-edge path took the same flag and ignored it. |
| `src/dead.hpp` (rest of `sccd_objective.hpp`) | `vf_objective` and `vf_objective_dir`, the last two of the ten generated functions. The header kept them while `find_root_newton` used `vf_objective_dir`; with that gone the header held no live code and was deleted. |
| `src/dead.cuh` (`warp_max_32`, `pow2`) | `warp_max_32`'s only caller is the block reduction in `dead.cuh` itself; `pow2` in `sccd_cuda_base.cuh` duplicated `sccd::pow2` and was called by nothing. |
| `src/spikes_only.hpp` (sweep helpers) | `prepare_B_block`, `tail_fill_B`, `scalar_count_range_two_lists` and `scalar_collect_range_two_lists`, whose only callers are the spikes here. They sat in an installed public header, so a consumer compiled them for nothing. |
| ~~`external/json/`~~ | **Returned to `benchmark/json/`.** Demoting it was wrong: it is not in the main `CMakeLists.txt`, which is what I checked, but `benchmark/scripts/bench.sh` configures it as a standalone project, and `sccd_bench` cannot run without what it produces. See `wip/DECISIONS.md`. |
| `research/` | Profiling artefacts and rendered notes, formerly the top-level `research/`. `cuda/` holds Nsight Compute profiles and the tables for `sccd_lower_bound_all_to_all`, itself a spike; `compare/` renders `wip/COMPARE.md`. Nothing in the shipped tree reads any of it. |
| `python/` | The SymPy code generators and the demoted analysis scripts. See [`python/README.md`](python/README.md) — the generators no longer reproduce the headers they seeded, which is why they are here rather than in `python/`. |
| `src/uniform_split.hpp` | Uniform interval splitting, extracted from `sccd_rootfinder.hpp`. A complete second implementation of the splitting job for both vertex-face and edge-edge, ~550 lines, reachable only through `SCCD_ADAPTIVE_SPLIT=0`, and ahead of the adaptive splitter on no real scene. |

## Building them

```sh
cmake -S . -B build -DSCCD_ENABLE_SPIKES=ON
```

Several of these do not build even then — `hacks/` in particular needs
dependencies this repository does not fetch. That is the point: the switch exists
so the code can be found and read, not so it can be relied on.
