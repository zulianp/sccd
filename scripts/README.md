# scripts

Developer helpers. None is needed to build, install or test SCCD — that is
`cmake` and `ctest`, with no options and no data.

| | |
|---|---|
| `mesh_sccd.sh` | Run the mesh demo on two frames of a scene, on the CPU or the GPU. |
| `crosscheck.sh` | Every query of every dataset through the C ABI, against its exact roots. |
| `plot.sh` | Figures from the tables `crosscheck.sh` writes. |
| `sync_alps.sh` | Push this working tree to the CSCS Alps cluster, for building CUDA. |
| `venv.sh` | Sourced by the above; finds a virtual environment, or falls through to the system `python3`. |

Each takes `--help`, or explains itself in the comment block at the top.

## What they need

`mesh_sccd.sh` needs a build with `-DSCCD_ENABLE_SMESH=ON` and the scene frames
under `data/`, which `benchmark/download_datasets.sh` fetches. Everything else in
`data/` is fetched the same way.

`crosscheck.sh` and `plot.sh` need the Python packages listed in
[`../python/README.md`](../python/README.md); `venv.sh` finds them.

## Two things worth knowing

**These are drivers, not tests.** Nothing here asserts anything about the
library except `crosscheck.sh`, which gates on the conservativeness invariant and
exits non-zero on a missed collision or a time of impact after the exact root.
The test suite is `ctest --test-dir build`.

**A build directory can hold a binary its cache says is disabled.** Turning
`SCCD_ENABLE_SMESH` off does not delete `mesh_sccd`, so a stale one stays
runnable and answers differently from the current source. `mesh_sccd.sh` prints
the binary it runs and takes `--build DIR` for that reason.

## Benchmarks and the cluster

The benchmark harnesses are in [`../benchmark/`](../benchmark/) — `bench.sh` for
the sweep, `download_datasets.sh` for the data, `assess.sbatch.sh` and
`launch.sbatch.sh` for Slurm. Building and running on Alps is in
[`../wip/ALPS.md`](../wip/ALPS.md), including why `ctest` does not work there
under `srun`.
