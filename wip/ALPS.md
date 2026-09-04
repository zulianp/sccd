# Running the tests on Alps

`ctest` does not work under `srun` on Alps right now. Run the binaries directly.

## What happens

```sh
srun --account=c40 --partition=debug --nodes=1 --ntasks=1 --gpus-per-task=1 \
     --time=00:20:00 --uenv=prgenv-gnu/24.11:v2 --view=default \
     ctest --test-dir build-hopper
```

prints `Start 1: <first test>`, that test's own output, and then nothing. It sits
there until the allocation hits its wall clock and Slurm reports `TIMEOUT`.
`stdbuf -oL -eL` does not help, and neither does `--timeout 120`, which should
have killed a stuck *test* and moved on -- so it is `ctest` itself that is stuck,
not any test under it.

## It is not the repository

Measured 2026-09-04. The same command stalls identically on
`$SCRATCH/sccd-npgap`, a tree synced that morning and verified green before the
day's changes. Two different source trees, two different build directories, same
stall after the first test. Do not go looking for the cause in SCCD.

## What to do instead

One `srun`, running each test binary in turn. The whole suite takes about eight
seconds on a GH200 node, so there is nothing to be gained from `ctest` here
anyway:

```sh
srun --account=c40 --partition=debug --nodes=1 --ntasks=1 --gpus-per-task=1 \
     --time=00:25:00 --uenv=prgenv-gnu/24.11:v2 --view=default bash -c '
  export SCCD_LIB_PATH=$PWD/build-hopper
  fail=0
  cd build-hopper
  for t in sccd_minimal sccd_broadphase_cell2d_test sccd_rootfinder_quad_test \
           sccd_c_abi_test sccd_aabb_test sccd_broadphase_cell2d_cuda_test \
           sccd_narrowphase_cuda_test; do
    ./$t || fail=1
  done
  python3 ../python/sccd_binding_test.py || fail=1
  echo "OVERALL=$fail"
'
```

`SCCD_LIB_PATH` is what `ctest` would have set for `sccd_binding_test`; the
binding needs it to find `libsccd.so`.

Redirect to a file and poll it rather than reading the stream over `ssh`: output
through a pipe is block-buffered, and a run that has already finished can look
like a hang for want of a flush. That is a separate trap from the one above, and
it is what made the first `ctest` attempt look like it had got further than it
had.
