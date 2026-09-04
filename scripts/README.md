# scripts

Two developer helpers that are not about benchmarking. Neither is needed to
build, install or test SCCD — that is `cmake` and `ctest`, with no options and no
data.

| | |
|---|---|
| `mesh_sccd.sh` | Run the mesh demo on two frames of a scene, on the CPU or the GPU. |
| `sync_alps.sh` | Push this working tree to the CSCS Alps cluster, for building CUDA. |

Both take `--help`, or explain themselves in the comment block at the top.

`mesh_sccd.sh` needs a build with `-DSCCD_ENABLE_SMESH=ON` and the scene frames
under `data/`, which `../benchmark/scripts/download_datasets.sh` fetches.

## Everything else moved

The benchmark harness and its drivers are in
[`../benchmark/scripts/`](../benchmark/README.md) — `bench.sh`, `crosscheck.sh`,
`plot.sh`, the Slurm jobs and the dataset download. They live with the benchmark
because that is what they measure and where their output belongs.

Building and running on Alps is in [`../wip/ALPS.md`](../wip/ALPS.md), including
why `ctest` does not work there under `srun`.

## One thing worth knowing

A build directory keeps a binary after the option that produced it is turned off.
Setting `SCCD_ENABLE_SMESH=OFF` does not delete `mesh_sccd`, so a stale one stays
runnable and answers differently from the current source. `mesh_sccd.sh` prints
the binary it runs and takes `--build DIR` for that reason.
