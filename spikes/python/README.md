# Demoted Python

Same bar as the rest of `spikes/`: nothing here is run by the build, by `ctest`,
or by any script in `benchmark/` or `scripts/`, and nothing in the shipped tree
imports it. Each file keeps its own name rather than being concatenated into one
`dead.py` the way `spikes/src/dead.{hpp,cpp,cu}` collects dead C++, because these
are standalone programs and the concatenation would not run.

Run them against the virtual environment described in
[`python/README.md`](../../python/README.md), with `python/` on `PYTHONPATH` for
the ones that import from it:

```sh
PYTHONPATH=python python3 spikes/python/ccd2d.py
```

Being demoted, they are not maintained and several do not work. `ccd2d.py` and
`ccd3d.py` fail under Numba — one on a type-inference pass, the other on an
`import` inside a jitted function — and did so before they were demoted; without
Numba installed, `myjit` is the identity and both run.

## The code generators

`codegen_*.py` emit C++ from SymPy expressions. They seeded four headers years
ago and were driven by `genmicrokernels.sh`. **They no longer reproduce those
headers**, which have been corrected and extended by hand since, so the headers
are now the source of truth and these are a historical derivation — useful for
seeing where an expression came from, not for regenerating anything.

The divergences, measured by running each generator and diffing against the
header it names:

| Generator | Header it seeded | How they differ now |
|---|---|---|
| `codegen_numerical_error.py` | `src/narrowphase/sccd_numerical_error.hpp` | The generator emits `max(max_abs, 1)^3`; the header computes `min(max_abs, 1)^3`, which is the certified TightInclusion bound the narrow phase's conservativeness rests on. The generator also carries a `use_ms` flag, and a C-style `sccd_get_numerical_error_vf_soa` interface the header replaced with `sccd::numerical_error_bound`. |
| `codegen_tolerance.py` | `src/narrowphase/sccd_tolerance.hpp` | The generator emits the two SoA functions. The header has five: both SoA forms, both array forms, and `clamp_domain_tol`, together with the note on why TightInclusion caps the per-axis tolerances. |
| `codegen_objective.py` | `src/narrowphase/sccd_objective.hpp` | The generator emits all ten objective, gradient and Hessian functions. Eight of them have no caller and now live in `spikes/src/dead.hpp`; the header keeps `vf_objective` and `vf_objective_dir`. The generator also emits `#include "smath.hpp"`, a path that no longer exists. |
| `codegen_interval.py` | `src/stuv.hpp` | The header was deleted outright in `2c7e657` (2025-12-08). The generator has had no target since. |
| `codegen_bisect.py` | none | Emits three midpoint splitters, `mid = (l + u) / 2` per axis. It was never wired to an output and no shipped code resembles it. |

`genmicrokernels.sh` was the driver. In the shipped tree it wrote straight into
`src/`, so once those paths were renamed it would have created three stray
headers that trip the `SCCD_PUBLIC_HEADERS` guard — and had the paths been
repointed instead, it would have silently reverted the `min` clamp above and
resurrected the eight dead objective functions. It now writes to a directory
given on the command line and refuses to write into `src/`.

## The rest

| File | Why it was demoted |
|---|---|
| `exact_roots_ee.py` | A SymPy solver for exact edge-edge roots, the Python half of the pair below. Superseded by the ground truth committed as `roots/<key>/toi.float64`, which the benchmark reads directly. No caller. |
| `exact_roots_vf.wls` | The vertex-face half, and the original: it is Mathematica, not Python. It was tracked as `python/root_vf.py`, where it raised `SyntaxError` on import for as long as it existed. Renamed to the extension that matches its language. |
| `ccd3d.py` | A grid-and-refine vertex-face search with two visualisations. Its only importer bound `find_root_vf` from here and then never called it -- that path goes through the C ABI binding. The implementations came with it out of `python/sccd_reference.py`, where they were reachable only from here. |
| `ccd2d.py` | Two-dimensional CCD, point against segment, with residual plots. SCCD is three-dimensional and nothing imports this. As above, the implementations moved out of `python/sccd_reference.py` with the script that used them -- together they were 958 of that file's 1344 lines. |
| `swept_axis_intersection.py` | A SymPy derivation that aligns a swept triangle to the x-axis and solves for intersection roots. A one-off, referenced by nothing. |
| `plot_ee_splitters.py` | Plots the adaptive bisection splitters for an edge-edge case. It reimplements `sccd::srootfinder::normal_equation_axis_splitters_ee` from `src/srootfinder.hpp` — a namespace and a file that no longer exist — and nothing calls it. Its two rendered figures are alongside it. |
