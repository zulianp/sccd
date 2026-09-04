# benchmark

The measurement harness. None of it is needed to build, install or test SCCD —
that is `cmake` and `ctest`, with no options and no data. The numbers it produces
are distilled in [`../docs/BENCHMARKS.md`](../docs/BENCHMARKS.md).

## Layout

| | |
|---|---|
| `*.exe.cpp` | The drivers, built as `sccd_bench`, `ti_oracle`, `sccd_refine_scaling`. |
| `scripts/` | Everything that runs them. See below. |
| `*.py` | Post-processing: aggregation, the HTML report, dataset conversion. |
| `json/` | Two converters for the datasets' JSON, as a **separate** CMake project. |
| `assessment/` | Committed measurement tables. Evidence for the keep-or-demote calls in `../wip/ASSESSMENT.md`. |
| `oracle/` | The accuracy gate's query sets, results, and its correctness contract. |
| `out/` | Everything a run produces. Ignored, and safe to delete. |

## The drivers

| Target | What it measures |
|---|---|
| `sccd_bench` | The full pipeline per case: prep, broad phase, narrow phase, false positives and negatives. |
| `ti_oracle` | Accuracy against the datasets' exact roots, and against TightInclusion. The gate; see [`oracle/README.md`](oracle/README.md). |
| `sccd_refine_scaling` | Scaling with element count, on refined triangle and quad meshes. |

`ti_oracle` needs `-DSCCD_ENABLE_TIGHT_INCLUSION=ON`; the other two need
`-DSCCD_ENABLE_SMESH=ON`.

## scripts/

| | |
|---|---|
| `bench.sh` | The sweep. Downloads the data, builds, runs every mode over every scene, post-processes, writes the report. |
| `download_datasets.sh` | Fetches the NYU CCD benchmark sets into `../data`. |
| `crosscheck.sh` | Every query of every dataset through the C ABI, against its exact roots. Writes `out/crosscheck/`. |
| `plot.sh` | Figures from the tables `crosscheck.sh` writes, into `out/crosscheck/figures/`. |
| `assess.sbatch.sh` | The variant sweep behind `../wip/ASSESSMENT.md`, as one Slurm job. |
| `launch.sbatch.sh` | `bench.sh` as a Slurm job. |
| `armadillo_ee_full.sh` | Every armadillo-rollers case, mode 0 against mode 2, per case. |
| `venv.sh` | Sourced by the Python-using scripts to find a virtual environment. |

Each explains itself in the comment block at the top.

## Why `json/` is a separate CMake project

Its two programs turn the datasets' `boxes/*.json` and `mma_bool/*.json` into the
raw arrays `bench.exe.cpp` reads. They need simdjson, fetched at configure time —
and the shipped library must configure and build with nothing fetched from the
network, so this cannot join the main project. `bench.sh` configures it
separately, into `../build_json`.

## Two things that will bite

**Modes are `0` and `2`.** `SCCD_NARROWPHASE_MODE=1` and `=3` no longer exist;
setting either warns and runs `Relaxed`, so a sweep that includes them measures mode
0 twice under two names. Older CSVs under `oracle/` carry the retired names
(`scalar`, `vector`, `ti-vec`, `ti-compat`) and the report renderer maps them.

**Measure inside one allocation.** Between `srun` allocations this harness varies
by about 40%, and several conclusions in this project were retracted after being
read off differences smaller than that. `assess.sbatch.sh` interleaves variants
inside one job for exactly this reason, and says so at length.
