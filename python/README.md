# Python

Two things live here: the binding that ships with the library, and the tools that
measure it. Nothing here is required to build or use SCCD from C++.

| | |
|---|---|
| `sccd.py` | **Shipped.** `ctypes` binding over the C ABI. Installed to `<prefix>/python`. Standard library only. |
| `sccd_binding_test.py` | The `sccd_binding_test` ctest. Checks the binding against the library it wraps. Standard library only. |
| `sccd_reference.py` | A NumPy reference search, used to cross-check the C ABI's answers. Needs NumPy; Numba is optional and makes it about 100x faster. |
| `sccd_crosscheck.py` | Runs a whole dataset through the C ABI and compares against its exact roots. |
| `sccd_read_queries.py` | Reads a dataset's query CSV. |
| `sccd_read_roots.py` | Reads a dataset's exact-root archive (Wolfram `.wxf` inside a tarball). |
| `sccd_read_mma.py` | Reads a dataset's hit/miss booleans. |
| `sccd_plot_toi_error.py` | Time-of-impact error, from a cross-check table. |
| `sccd_plot_fp_fn.py` | False positives and negatives per scene, from the same tables. |
| `sccd_strip_nonascii.py` | Repairs PLY headers whose comments are not ASCII, which the mesh reader rejects. |

`benchmark/*.py` post-processes the benchmark driver's output and is documented
where it is used. Demoted scripts, including the SymPy code generators, are in
[`spikes/python/`](../spikes/python/README.md).

## Using the binding

```python
import sccd

hit = sccd.find_root_vf(max_iter=69, tol=1e-8,
                        sv=(0.25, 0.25,  1.0),                    # vertex, start of step
                        s1=(0, 0, 0), s2=(1, 0, 0), s3=(0, 1, 0), # triangle, start
                        ev=(0.25, 0.25, -1.0),                    # vertex, end of step
                        e1=(0, 0, 0), e2=(1, 0, 0), e3=(0, 1, 0)) # triangle, end
if hit:
    print(hit.t, hit.u, hit.v)          # 0.499999999 0.25 0.25
```

`find_root_ee`, `find_root_bisection_vf` and the two `find_root_tight_inclusion_*`
oracles have the same shape. The result is a `Root(hit, t, u, v)` named tuple,
falsey on a miss. `dtype="float"` selects the single-precision export where one
exists. `available()` reports which of the C ABI's seven exports the loaded
library has, which is how to tell whether the build has TightInclusion.

The reported time of impact is conservative: at or before the true one, never
after. See [`../docs/API.md`](../docs/API.md) for the C ABI itself.

One query at a time, with no broad phase and no batching — that is what the C ABI
offers. For a pipeline, call the C++ headers.

### Finding the library

`SCCD_LIB_PATH` wins, and may name the library file or a directory holding it.
Otherwise `sccd.py` looks beside itself, then asks the system loader. The library
is loaded on first call, not on import, so importing never fails for want of a
build.

```sh
cmake -S . -B build && cmake --build build -j
SCCD_LIB_PATH=build python3 python/sccd.py     # prints the export list and a smoke test
```

## Running the tools

`sccd.py` and `sccd_binding_test.py` need nothing but the standard library. The
rest need NumPy, SymPy, pandas and Matplotlib:

```sh
python3 -m venv .venv
.venv/bin/pip install numpy sympy pandas matplotlib
```

`scripts/venv.sh`, which `scripts/crosscheck.sh` and `scripts/plot.sh` source,
activates `$SCCD_VENV` if set, then `.venv`, then `data/venv`, and falls through
to the system `python3` if it finds none.

The cross-check needs the benchmark datasets under `data/`, which are gigabytes
and are not in the repository; `benchmark/download_datasets.sh` fetches them.

```sh
scripts/crosscheck.sh          # every query through the C ABI, against exact roots
scripts/plot.sh                # figures from the tables it wrote
```

`sccd_crosscheck.py` and `sccd_plot_toi_error.py` exit non-zero on a missed
collision or a time of impact after the exact root. False positives, and times of
impact earlier than the exact root, are reported and are not failures: they cost
work and never safety.
