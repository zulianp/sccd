# CUDA

SCCD runs its whole pipeline on the GPU: swept AABBs, broad phase, narrow phase.
The device API mirrors the host API name for name — every entry point is the
host one under `sccd::device`, with the same argument order and the same
meaning — so porting a working host pipeline is a matter of qualifying the calls
and moving the data.

The guarantee does not change. A time of impact from a device kernel is at or
before the true one, exactly as on the host, and a collision that exists is
reported. That is checked on every scene and every mode in
[`BENCHMARKS.md`](BENCHMARKS.md).

## Build

CUDA is off by default. Turn it on and name the architecture:

```sh
cmake -S . -B build -DSCCD_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90
cmake --build build -j
```

`90` is Hopper; use `80` for Ampere, `70` for Volta. Nothing else is needed —
the device path has no dependency the host path does not have.

A complete, runnable example builds with it:

```sh
cmake --build build -j --target sccd_minimal_cuda
./build/sccd_minimal_cuda
```

`demo/cuda/sccd_minimal_cuda.exe.cpp` is the file to read. It is
`demo/sccd_minimal.exe.cpp` on the device, deliberately: the same two triangles,
the same stages, the same exact answer of `t = 0.5`. Reading them side by side
isolates what is actually different about CUDA. It exits 0 with a message when
no device is present.

## The one rule that is easy to get wrong

**A `T**` argument is a device array of device pointers.**

Geometry is structure-of-arrays on both processors: `T** points` is three row
pointers, one per axis. On the device, both levels must live in device memory —
the rows *and* the array holding them.

This compiles and is wrong:

```cpp
T* rows[3] = {d_x, d_y, d_z};   // rows are device pointers, good
sccd::device::narrow_phase_vf<T, I>(n, ..., rows, ...);   // WRONG: rows is a host array
```

The type is right, so nothing complains; the kernel then dereferences a host
address. Upload the table too:

```cpp
T** d_points = nullptr;
cudaMalloc(&d_points, sizeof(T*) * 3);
cudaMemcpy(d_points, rows, sizeof(T*) * 3, cudaMemcpyHostToDevice);
sccd::device::narrow_phase_vf<T, I>(n, ..., d_points, ...);   // right
```

The same applies to element tables (`I** faces`, `I** edges`) and to the
six-row AABB blocks the broad phase reads and writes.

To read one row back out of such a table — the cumulative-maximum step needs
this — use `sccd::device::soa_device_row(table, row)` rather than indexing the
table from host code.

## Entry points

| Stage | Call | Header |
|---|---|---|
| Swept AABBs | `sccd::device::compute_aabbs` | `sccd_vaabb.cuh` |
| Sort axis | `sccd::device::choose_axis`, `sort_along_axis`, `enumerate` | `sccd_broadphase.cuh` |
| Candidate pairs, two lists | `sccd::device::count_overlaps`, `collect_overlaps` | `sccd_broadphase.cuh` |
| Candidate pairs, one list | `sccd::device::count_self_overlaps`, `collect_self_overlaps` | `sccd_broadphase.cuh` |
| Cell-list broad phase | `sccd::device::cell2d_*` | `sccd_cell2d_broadphase.cuh` |
| Narrow phase, triangles | `sccd::device::narrow_phase_vf`, `narrow_phase_ee` | `sccd_narrowphase.cuh` |
| Narrow phase, quads | `sccd::device::narrow_phase_vq` | `sccd_narrowphase_vq.cuh` |

The broad phase is count-then-fill on both processors: the count pass writes
prefix offsets, so the exact output size is known before a pair is written. Read
the total back from the last offset:

```cpp
std::ptrdiff_t total = 0;
cudaMemcpy(&total, d_offsets + n_elements, sizeof(total), cudaMemcpyDeviceToHost);
```

That is what makes the fill pass a plain parallel write with no atomics, no
compaction and no per-thread growable buffer.

## Output

The time of impact comes back in device memory, so reading it takes a copy.

`sccd::ToiOutput::Earliest` writes a single scalar — allocate one element — and
lets every query prune against the shared running minimum.
`sccd::ToiOutput::PerPair` writes one element per candidate pair and does
correspondingly more work, because there is no shared bound.

The kernels initialise the output themselves, with one exception worth
programming around: **given zero candidate pairs they return immediately and
leave the buffer untouched.** Seed it with `max_toi` so an empty candidate list
reads as "no collision" rather than as whatever was in the allocation.

The return value is a status, not a result: `0` means the call succeeded.

## Choosing a processor

The GPU is not uniformly faster, and which phase wins is not the same on every
scene. Measured over six scenes on a GH200 node ([`BENCHMARKS.md`](BENCHMARKS.md)
has the tables):

- The **broad phase** is the stage whose shape suits a GPU — count, prefix sum,
  scatter, with no sequential window walk. It runs 2.2× to 4.2× faster than
  Grace on five of the six scenes, and 1.2× slower on the sixth.
- The **narrow phase** is usually the slower one on the GPU: 0.27× to 0.76× of
  Grace on five scenes, and 1.26× on the sixth.
- **End to end** the GPU wins five of six by 1.5× to 2.9×, and loses one at
  0.62×.

The two exceptions are different scenes, which is the point: the phase that
wins is a property of the geometry, not a fixed property of the processor. And
the scene the GPU loses outright is the largest in the set by candidate pairs —
30.5 million per step — so "bigger problem, use the GPU" is not a safe rule.
Measure the scene you have.

One trap when comparing: `Tight` is the host's faster mode and the device's
slower one, so pinning the same `SCCD_NARROWPHASE_MODE` on both sides races the
host's best kernel against the device's worst. Compare each processor at its own
best mode unless you specifically mean otherwise.

## Precision

Root finding computes in **double** whatever the storage type. In single
precision the certified numerical error bound and the tolerances that terminate
the search are too close together for the guarantee to survive.

Storing geometry as `float` is fine and is where single precision earns its
throughput; the search widens to double internally, which is exact, and narrows
the result toward negative infinity on the way out. Rounding a time of impact
*down* is always safe; rounding to nearest could round it up, past the root the
search proved.

## Through smesh

If you are using the optional smesh integration, `CCD<T>::create(mesh,
smesh::EXECUTION_SPACE_DEVICE)` selects the device path and the buffers are
managed for you — the pointer-table rule above is handled inside. The staged
interface is the same on both processors:

```cpp
auto ccd = sccd::CCD<double>::create(mesh, smesh::EXECUTION_SPACE_DEVICE);
ccd->broad_phase_prep(points_t0, points_t1);
ccd->broad_phase_fv_step(v_overlap, f_overlap);
ccd->broad_phase_ee_step(e0_overlap, e1_overlap);
ccd->narrow_phase(max_toi, vf_tois, ee_tois, max_depth, tol);
```

`max_toi` is in/out: it bounds the search on the way in and, for the default
`ToiOutput::Earliest`, comes back holding the earliest time of impact. Take it
by reference in any wrapper you write around this.

Quads are handled on the device as well as triangles; `QUADSHELL4` meshes
dispatch to the vertex-quad kernel automatically.
