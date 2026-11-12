# Building and running SCCD on Nibi's AMD GPUs

Nibi (SHARCNET, Digital Research Alliance of Canada) is where SCCD's HIP build is
developed and checked: it has 5 nodes with 4 AMD Instinct MI300A APUs each
(gfx942, 64-lane wavefronts). Fir, the Alliance's other new cluster, has NVIDIA
H100s only.

## Getting on

`ssh nibi` is an alias in `~/.ssh/config` (`nibi.alliancecan.ca`, user
`zulianp`) with a four-hour `ControlMaster`. Multifactor authentication is
mandatory, so the first connection has to be made interactively and approved in
Duo; for the next four hours every `ssh -o BatchMode=yes nibi '<cmd>'` reuses that
connection without asking. From inside the cluster, `ssh l1` (or `l3`, `l4`)
reaches the other login nodes without another approval.

A connection that hangs instead of failing is usually a stale control socket
whose underlying connection has died: `ssh -O exit nibi` and log in again.

## Toolchain

```sh
module load StdEnv/2023 rocm-compilers/19.0.0 hip/6.4.1 rocprim/3.4.0 \
            hipcub/3.3.0 rocthrust/3.3.0 cmake/3.31.0
```

That is HIP 6.4.1 on ROCm's LLVM 19, with g++ 14.3 as the host compiler. `hip`
needs `rocm-compilers` loaded first, and loading `rocm-compilers` replaces the
default `gcc/12.3`. `$SCRATCH/hip-env.sh` on Nibi holds exactly this line.

**A login node can lose its software filesystem.** On 2026-11-12 every `cmake`
under `/cvmfs` returned `Input/output error` on `l2` while `l1`, `l3` and `l4`
read it normally. When a module's binary fails like that, `ssh l1` and work there;
`$SCRATCH` is the same on every node.

## Building

Build in `$SCRATCH` (`/scratch/zulianp`). The login nodes have no GPU but
compile for one, so building there is fine:

```sh
scripts/sync_alps.sh nibi:/scratch/zulianp/sccd      # from the Mac; the script takes any destination
ssh nibi
ssh l1                                               # if l2's /cvmfs is failing
cd $SCRATCH/sccd && source $SCRATCH/hip-env.sh
cmake -S . -B build-hip -DCMAKE_BUILD_TYPE=Release \
      -DSCCD_ENABLE_CUDA=ON -DSCCD_ENABLE_HIP=ON \
      -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_HIP_COMPILER=$(which clang++)
cmake --build build-hip -j 16
```

`SCCD_HIP_ARCHITECTURES` defaults to `gfx942` and `SCCD_HIP_WAVEFRONT_SIZE` to
64, which is MI300A.

### With smesh

SCCD's mesh path (`SCCD_ENABLE_SMESH`) needs smesh built for the same backend:
`-DSMESH_ENABLE_CUDA=ON -DSMESH_ENABLE_HIP=ON`, from smesh's `hip-port` branch
until it is merged. Configuring SCCD against an smesh of the other backend stops
with an error. Benchmarks and accuracy checks need it at `geom_t=double`:

```sh
cmake -S smesh -B smesh/build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=$(which g++) -DCMAKE_C_COMPILER=$(which gcc) \
      -DCMAKE_HIP_COMPILER=$(which clang++) \
      -DSMESH_ENABLE_CUDA=ON -DSMESH_ENABLE_HIP=ON \
      -DSMESH_GEOM_TYPE=float64 -DSMESH_REAL_TYPE=float64 \
      -DSMESH_SCALAR_TYPE=float64 -DSMESH_ACCUMULATOR_TYPE=float64 \
      -DSMESH_ENABLE_OPENMP=ON -DSMESH_ENABLE_MPI=OFF -DSMESH_ENABLE_MPISORT=OFF \
      -DSMESH_ENABLE_DEMO=OFF -DSMESH_ENABLE_RYAML=OFF -DSMESH_ENABLE_DEV_MODE=OFF \
      -DCMAKE_INSTALL_PREFIX=$SCRATCH/smesh-f64-hip
cmake --build smesh/build -j 16 && cmake --install smesh/build
```

then add `-DSCCD_ENABLE_SMESH=ON -Dsmesh_DIR=$SCRATCH/smesh-f64-hip/lib64/cmake/smesh`
to the SCCD configure line above. `$SCRATCH/build-smesh-hip.sh` on Nibi does both.

## Running

The MI300A nodes sit in the general GPU partitions with GRES `gpu:mi300a:4`
(`gpubase_bygpu_b1` at 3 h through `_b5` at 7 days, the `bynode` equivalents,
`gpubase_interac` at 8 h). Name the account and the GPU and leave the partition
out: the scheduler picks one from the time limit, and naming one is rejected.

```sh
srun --account=def-teseo_gpu --gpus-per-node=mi300a:1 --cpus-per-task=8 \
     --mem=32G --time=00:20:00 <binary>
```

`$SCRATCH/run-hip-checks.sh` on Nibi is the batch form of the whole check: every
test binary, `sccd_minimal_cuda`, and `ti_oracle` on cloth-funnel, whose queries
and roots are in `$SCRATCH/sccd-data/` (copied from the Mac's `data/`; nothing
else of the datasets is on Nibi). Submit it with `sbatch $SCRATCH/run-hip-checks.sh`.

The Slurm association with `def-teseo_gpu` did not exist when this account was
first used on 2026-11-12 -- submissions failed with "Invalid account or
account/partition combination" -- and appeared later the same day. If that
message comes back, it is the association, not the job.

## Status (2026-11-12)

Built with `-DSCCD_ENABLE_TIGHT_INCLUSION=ON`, on one MI300A:

- all ten test binaries and `sccd_minimal_cuda` pass;
- `ti_oracle` on cloth-funnel: no missed collisions and no late times of impact
  over 27,092 gated queries, and hits, false positives and the exact-root counts
  identical to the GH200's in every device row. `lateTI` and the error medians
  differ, as they do between two runs on the GH200.

The other five scenes have not been run on AMD.

With smesh (`hip-port` branch, `geom_t=double`), on one MI300A:

- smesh's 22 tests pass. They run its host paths only: `smesh_hierarchies_test`
  reaches the device with `SMESH_EXECUTION_SPACE=device`, but that branch is
  unfinished (a `TODO` where the prolongation belongs, then host functions handed
  device pointers) and fails on the GH200 as well, by segmentation fault. So no
  test covers smesh's own device kernels on either vendor.
- SCCD's tests pass, and `sccd_mesh_cuda_test` on cloth-funnel frames 100-101
  and 220-221 (converted to float64) prints the same as on the GH200: identical
  broad-phase pairs on CPU and GPU, Tight times of impact within 1e-2 of the CPU,
  and in Relaxed two edge-edge pairs where the GPU is earlier (0.64 against
  0.852, 0 against 0.191). `$SCRATCH/run-hip-smesh-checks.sh` is the batch job.
