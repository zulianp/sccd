# Competitor harnesses

Third-party continuous collision detection run over the same datasets, the same
cases and the same phase split as `benchmark/bench.exe.cpp`, writing the same CSV
schema, so a row from a competitor and a row from SCCD can sit in one file and be
read by `benchmark/report/` without special-casing.

Two are set up here:

| directory | library | what it is |
| --- | --- | --- |
| `scalable_ccd/` | [Scalable CCD](https://github.com/continuous-collision-detection/scalable-ccd) | Sweep and Tiniest Queue broad phase with Tight-Inclusion narrow phase, CPU and CUDA |
| `accd/` | [Additive CCD](https://github.com/ipc-sim/ipc-toolkit) | the conservative advancement method of the IPC line of work, as `ipc::AdditiveCCD` |

Each is compared with SCCD on the question it answers:

* **Scalable CCD: the earliest time of impact.** Both libraries run their
  earliest-time-of-impact path per case -- one query kind of one simulation
  step -- over identical broad-phase candidates. A step's answer is the minimum
  over its cases and is scored, signed, against the step's earliest exact root.
* **Additive CCD: per collision pair.** Both narrow phases are timed on SCCD's
  broad-phase candidates, and scored on every curated query at the coordinates
  its exact root was computed for.

## The rules this directory follows

**Competitor source is never modified.** The same rule the repository already
applies to TightInclusion: a comparison is only worth reporting if the other
implementation is the one its authors published. Every adaptation happens on this
side of the boundary. Where a competitor is fetched, it is pinned to a released
tag, so a rerun a year from now measures the same code.

**Both sides get the same compiler flags.** `cmake/SCCDDependencies.cmake`
already does this for TightInclusion, and the reason is the same here: a
competitor built at the toolchain's defaults while SCCD is built for the host CPU
is not a measurement of either library. The competitor target inherits
`CMAKE_CXX_FLAGS` and the configuration's release flags.

**A competitor is built the way its authors ship it.** Compile-time options are
part of the algorithm being measured. Where a measurement needs a non-default
build -- Scalable CCD's per-query times of impact -- that build is a separate
binary, its rows carry their own mode name, and nothing timed is read from it.

**The phase split is the one the article uses.** The broad phase includes
building whatever acceleration structure it needs — for Scalable CCD that is
constructing the boxes, copying them to the device, and `BroadPhase::build`, not
only the sweep itself. The mesh upload is outside every timed region, as it is in
`sccd_bench`. The narrow phase is the root finding. A library that interleaves
the two, as Scalable CCD does when it batches overlaps under a memory limit, has
each call timed and accumulated into its own total.

**Both sides get the same cores.** A library that exposes only a per-pair kernel
is driven by the parallel loop its own high-level entry point uses, so a timing
comparison is threads against threads rather than one thread against 72.

**The parameters are each library's own.** A knob with a genuine counterpart is
matched; a knob without one takes the competitor's documented default and is
recorded in this file. Names that look alike are not counterparts: Scalable CCD's
`max_iterations` is a budget on box checks whose exhaustion drops the box, not a
subdivision depth.

**Indices are the benchmark's.** Accuracy is scored per query against the same
ground truth, which means a competitor's collision pairs have to be expressed in
the same numbering the curated queries use. Vertices and faces come from the mesh
in its own order; edges come from `benchmark_ordered_edges`, the deterministic
ordering `bench.exe.cpp` derives from the face list. Passing exactly those arrays
to a competitor makes the pairs it returns directly comparable, with no mapping
table to get wrong.

## Building

Off by default, because the shipped library builds with no external dependencies
and these pull several.

```sh
cmake -S . -B build -DSCCD_ENABLE_COMPETITORS=ON -DSCCD_ENABLE_SMESH=ON \
      -Dsmesh_DIR=<prefix>/lib/cmake/smesh
cmake --build build --target scalable_ccd_bench accd_bench
```

Add `-DSCCD_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90` for the device path.
On Alps, pass `-DCMAKE_C_COMPILER` as well as `-DCMAKE_CXX_COMPILER`: left to
itself CMake finds the system GCC 7 for C, OpenMP resolves to its `libgomp`, and
the link fails on `__aarch64_cas4_sync`.

**smesh must be built with double `geom_t`**, and the prepared frames converted
to match (`SCCD_GEOM_PRECISION=float64`). The dataset's exact roots belong to the
rational query coordinates; a float32 mesh answers about coordinates that differ
at about `1e-7`, which near a grazing contact moves the true time of impact by up
to `1e-4` and makes a correct kernel look late. The float64 explicit template
instantiations live on smesh's `sideset` branch -- `master` instantiates `f32`
only, and a `SMESH_GEOM_TYPE=float64` build of it leaves
`mesh_from_folder<int, double>` undefined inside `libsmesh.a`. On GCC that build
also wants `SMESH_ENABLE_DEV_MODE=OFF`, because `-Werror=restrict` rejects
existing semistructured code.

`SCCD_SCALABLE_CCD_TOI_PER_QUERY` (default `OFF`) builds Scalable CCD with
per-query times of impact. It is a diagnostic build and not part of the
comparison.

## Running

The same shape as `sccd_bench`, and the same environment variables select the
scene range, so `benchmark/scripts/sweep.sh` can drive either:

```sh
SCCD_BENCH_EXECUTION_SPACE=device SCCD_BENCH_CASE_BEGIN=0 SCCD_BENCH_CASE_END=200 \
    ./build/benchmark/competitors/scalable_ccd_bench <data-dir> armadillo-rollers
```

`--header` prints the CSV schema, exactly as `sccd_bench --header` does. The two
must agree; `ctest -R competitor_schema` checks that they do.

### The comparison on Alps

`sweep_comparison.sh` runs the comparison over whole scenes as resumable Slurm
chunks, one (scene, repeat, case range) per job on one GH200 module: one Hopper
and one Grace bound to 72 CPUs. Each chunk runs `run_comparison.sh`, which runs
every library over that range in its own process and writes one CSV through a
temporary file. A chunk whose CSV exists is skipped, so an interrupted sweep
resumes.

```sh
# on the Alps login node
bash benchmark/competitors/sweep_comparison.sh --dry-run
nohup bash benchmark/competitors/sweep_comparison.sh > sweep_comparison.log 2>&1 &
bash benchmark/competitors/sweep_comparison.sh --merge
python3 benchmark/competitors/compare_table.py $SCRATCH/sccd-compare/compare.csv
```

The defaults are all six scenes, three repeats and 400 cases per chunk, which is
66 chunks; `--jobs N` keeps N of them in flight at once. The binaries come from
`$SCCD_COMPETITOR_BUILD` (default `$SCRATCH/sccd/build-comp-f64`) and the data
from `$SCCD_DATA_DIR`. The
binding matters for Scalable CCD's host broad phase, which is oneTBB and follows
the CPU affinity mask, not `OMP_NUM_THREADS`.

`compare_table.py` prints two tables per scene, one per comparison:

| table | SCCD narrow phase | competitor | accuracy |
| --- | --- | --- | --- |
| earliest time of impact | `narrow_ms`, `ToiOutput::Earliest` | Scalable CCD device, as shipped | per step: earliest reported minus earliest exact root, every repeat |
| per collision pair | `narrow_ms_s1`, `ToiOutput::PerPair` | ACCD `narrow_ms` over SCCD's host candidates | per curated query: `fp`, `fn`, late, early margins |

## Scalable CCD

**Pinned commit.** `ccf4bc19`, double precision (`SCALABLE_CCD_USE_DOUBLE`, its
default). Parameters: `max_iterations = -1` (unlimited; any finite value makes
exhaustion drop boxes), `tolerance = 1e-6` (a co-domain tolerance, the value its
tests use), `min_distance = 0`, `allow_zero_toi = true`. Thread counts are those
of `cuda::ccd`.

**Its narrow phase is CUDA only.** A host build measures the broad phase alone.

**Earliest time of impact, per case.** Built as shipped, its narrow phase prunes
every box against the running earliest time of impact. `sccd_bench` runs SCCD's
earliest-time-of-impact path once per case from a bound of 1, and this harness
runs Scalable CCD's the same way, so both do the same work: each case starts from
a fresh bound and `s0_toi` is the case's earliest time of impact. Their
`cuda::ccd` carries one bound from the vertex-face pass into the edge-edge pass;
for a deterministic search its answer for a step equals the minimum of the
harness's two case rows, which `--their-ccd` checks. At the pinned commit the two
differ between runs because of the defect below.

| mode | what the row holds |
| --- | --- |
| `scalable-ccd-device` | prep, broad and narrow timings, and the case's earliest time of impact |
| `scalable-ccd-host` | the broad phase on the host |

**Its broad phase agrees with ours to the pair.** On armadillo-rollers `5ee` both
produce `127099` candidates, and `broad_fn` is zero on every case. That is what
two correct sweep-and-prune broad phases give, since the overlap set is a
property of the boxes.

**Its narrow phase loses boxes when its subdivision buffer overflows.**
`CCDBuffer::push` tests `is_full()` and increments the tail as two separate
operations across device threads. When a level of the search needs more room
than the buffer holds, the tail can wrap onto the head with no overflow flagged,
and unprocessed domains are overwritten. The buffer holds `2 * n` domains for `n`
queries, so this happens whenever a list is short relative to the search it
starts. The consequences are missed collisions, times of impact after the true
one, and answers that differ between identical runs, in both builds and through
`cuda::ccd` itself. The comparison reports the library as published; the
diagnosis, and the patched copy used to confirm it, are in
`wip/DECISIONS.md` section 8.

**Checking it.** Two modes separate the library's answer from the harness's
composition:

```sh
# their own entry point on one step
scalable_ccd_bench --their-ccd    <data-dir> armadillo-rollers 151
# the step's curated queries, bare and followed by <pad> pairs rejected on
# their first test; a correct search reports the same collisions both ways
scalable_ccd_bench --buffer-probe <data-dir> armadillo-rollers 107 1000000
```

## Additive CCD

`ipc::AdditiveCCD` from the IPC toolkit `v1.6.0`, with its default
`conservative_rescaling = 0.9`. It is host code only; the toolkit's CUDA files
are all broad phase.

It answers one pair at a time and has no parallelism of its own, so the harness
loops over the candidates with the same `tbb::parallel_for` the toolkit's own
`Candidates::compute_collision_free_stepsize` uses. Both narrow phases then run
on the same cores: a serial loop here would time one thread against SCCD's 72.

It is a narrow phase, so its cost is measured on SCCD's sweep broad-phase
candidates: `prep_ms` and `broad_ms` are SCCD's, and `narrow_ms` is additive CCD
over that list, beside SCCD's per-pair `narrow_ms_s1` over the same list. The
gather from element indices to the eight points it takes sits between the two
timed regions. Its accuracy is scored on the curated queries at their own
coordinates (`queries_raw/<key>.f64`), which is where `sccd_bench` scores SCCD;
`query_narrow_ms` is that pass.

The method advances conservatively and cannot overshoot, so a non-zero `toi_late`
means the harness is wrong and the binary exits non-zero. The columns that
separate it from the others are `toi_med_early` and `toi_max_early`: how far
before contact it stops.
