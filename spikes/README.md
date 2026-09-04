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
| `src/broadphase_lb.hpp` | A lower-bound broad phase with no caller outside `hacks/`, yet installed as public API by the old flat header glob. |
| `src/sccd.hpp` | Aggregate header used only by `hacks/`. |
| `src/cuda/sccd_broadphase_warp.*` | `count_overlaps_warp` is defined and explicitly instantiated, and called by nothing repository-wide. |
| `src/cuda/sccd_lower_bound_all_to_all.*` | Reachable only from its own demo; no path from the public API. Its demo, `demo/cuda/mesh_lower_bound_cuda.exe.cpp`, is run with `scripts/mesh_sccd.sh <scene> <frame0> <frame1> --binary mesh_lower_bound_cuda`; the shell script that used to drive it lived in `scripts/` and carried hard-coded paths into one person's scratch directory. |
| ~~`external/json/`~~ | **Returned to `benchmark/json/`.** Demoting it was wrong: it is not in the main `CMakeLists.txt`, which is what I checked, but `benchmark/bench.sh` configures it as a standalone project, and `sccd_bench` cannot run without what it produces. See `wip/DECISIONS.md`. |
| `python/` | The SymPy code generators and the demoted analysis scripts. See [`python/README.md`](python/README.md) — the generators no longer reproduce the headers they seeded, which is why they are here rather than in `python/`. |
| `src/uniform_split.hpp` | Uniform interval splitting, extracted from `sccd_rootfinder.hpp`. A complete second implementation of the splitting job for both vertex-face and edge-edge, ~550 lines, reachable only through `SCCD_ADAPTIVE_SPLIT=0`, and ahead of the adaptive splitter on no real scene. |

## Building them

```sh
cmake -S . -B build -DSCCD_ENABLE_SPIKES=ON
```

Several of these do not build even then — `hacks/` in particular needs
dependencies this repository does not fetch. That is the point: the switch exists
so the code can be found and read, not so it can be relied on.
