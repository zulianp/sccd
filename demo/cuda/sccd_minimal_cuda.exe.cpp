// SCCD end to end on the GPU, with nothing but CUDA and the standard library.
//
// This is `demo/sccd_minimal.exe.cpp` run on the device, deliberately: the same
// scene, the same stages, the same exact answer. Read them side by side and the
// only differences are the ones that are actually about CUDA.
//
// Two disconnected triangles, one of which passes through the other during the
// step. Both stages run: swept AABBs and the sweep-and-prune broad phase to get
// candidate pairs, then the narrow phase to turn those into a time of impact.
// The scene is built so the exact answer is known, so the result is verified
// rather than merely printed.
//
// Build (CUDA is off by default):
//
//     cmake -S . -B build -DSCCD_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=<arch>
//     cmake --build build -j --target sccd_minimal_cuda
//     ./build/sccd_minimal_cuda
//
// ## The three things that differ from the host
//
// 1. **Every entry point is the host one, under `sccd::device`.** Same names,
//    same argument order, same semantics. Porting a host pipeline is a matter
//    of qualifying the calls and moving the data.
//
// 2. **A `T**` argument is a *device* array of *device* pointers.** Not a host
//    array. Passing `T* rows[3]` from the stack compiles -- the type is right --
//    and then the kernel dereferences a host address. `upload_table` below is
//    the whole of what is needed: copy each row to the device, then copy the
//    array of those device pointers to the device too.
//
// 3. **Outputs live in device memory.** The time of impact comes back through a
//    device buffer, so it takes a copy to read.
//
// Everything else -- the count-then-fill broad phase, the one-sided
// conservativeness guarantee, the meaning of every parameter -- is unchanged.

#include "sccd_base.hpp"
#include "sccd_broadphase.cuh"
#include "sccd_narrowphase.cuh"
#include "sccd_narrowphase_mode.hpp"
#include "sccd_vaabb.cuh"

#include <cuda_runtime.h>

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <vector>

namespace {

using T = double;
using I = int;

constexpr int kDim = 3;

// Fail loudly and immediately: a CUDA error left unchecked shows up later as a
// wrong answer, which is the one outcome a conservative library must not have.
#define CHECK(call)                                                                       \
    do {                                                                                  \
        const cudaError_t err_ = (call);                                                  \
        if (err_ != cudaSuccess) {                                                        \
            std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(err_), __FILE__,    \
                        __LINE__);                                                        \
            std::exit(1);                                                                 \
        }                                                                                 \
    } while (0)

// The scene, identical to the host example.
//
//   Triangle A, stationary, in the plane z = 0:   nodes 0, 1, 2
//   Triangle B, descending by 2 over the step:    nodes 3, 4, 5
//
// The two triangles share no vertex, so every vertex-face and edge-edge pair
// between them is a genuine candidate. B is tilted, so node 3 is its lowest
// point and leads. The motion is purely vertical, so node 3 stays at (0.5, 0.5)
// in x-y for the whole step, and that point is strictly inside triangle A. Its
// height runs 1.0 -> -1.0, reaching z = 0 at exactly t = 1/2, and every other
// point of B is above it throughout. The first contact is therefore node 3
// against face A at
//
//     t = 0.5, exactly.
//
// All coordinates are exactly representable, so "exactly" is meant literally.
constexpr T kExactToi = 0.5;

struct Scene {
    std::vector<T> x0, y0, z0;  // start of the step, one array per axis
    std::vector<T> x1, y1, z1;  // end of the step
    std::vector<I> f0, f1, f2;  // faces: one row per vertex slot
    std::vector<I> e0, e1;      // edges: one row per endpoint
};

Scene make_scene() {
    Scene s;
    //          A                B
    s.x0 = {0.0, 2.0, 0.0, 0.5, 2.5, 0.5};
    s.y0 = {0.0, 0.0, 2.0, 0.5, 0.5, 2.5};
    s.z0 = {0.0, 0.0, 0.0, 1.0, 2.0, 2.0};

    s.x1 = s.x0;  // B translates straight down by 2
    s.y1 = s.y0;
    s.z1 = {0.0, 0.0, 0.0, -1.0, 0.0, 0.0};

    s.f0 = {0, 3};  // face 0 = (0,1,2), face 1 = (3,4,5)
    s.f1 = {1, 4};
    s.f2 = {2, 5};

    s.e0 = {0, 1, 2, 3, 4, 5};  // three edges per triangle
    s.e1 = {1, 2, 0, 4, 5, 3};
    return s;
}

// --- the only CUDA-specific plumbing in the file ---------------------------

template <typename E>
E* upload(const std::vector<E>& host) {
    E* d = nullptr;
    CHECK(cudaMalloc(&d, sizeof(E) * host.size()));
    CHECK(cudaMemcpy(d, host.data(), sizeof(E) * host.size(), cudaMemcpyHostToDevice));
    return d;
}

template <typename E>
E* device_zeros(const std::size_t n) {
    E* d = nullptr;
    CHECK(cudaMalloc(&d, sizeof(E) * n));
    CHECK(cudaMemset(d, 0, sizeof(E) * n));
    return d;
}

// A `T**` the library can dereference on the device: `rows` are already device
// pointers, and this copies that array of pointers to the device as well. This
// is the step that is easy to skip, because skipping it still compiles.
template <typename E>
E** upload_table(E* const* rows, const int n) {
    E** d = nullptr;
    CHECK(cudaMalloc(&d, sizeof(E*) * n));
    CHECK(cudaMemcpy(d, rows, sizeof(E*) * n, cudaMemcpyHostToDevice));
    return d;
}

// Six device rows -- 0..2 the minimum corner, 3..5 the maximum -- plus the
// device table pointing at them, and the index array the sort permutes.
struct DeviceAabbs {
    T* rows[2 * kDim] = {};
    T** table = nullptr;
    I* idx = nullptr;

    void allocate(const std::ptrdiff_t n) {
        for (int d = 0; d < 2 * kDim; ++d) rows[d] = device_zeros<T>((std::size_t)n);
        table = upload_table(rows, 2 * kDim);
        CHECK(cudaMalloc(&idx, sizeof(I) * (std::size_t)n));
        sccd::device::enumerate<I>(0, n, idx);
    }
};

std::ptrdiff_t read_total(const std::ptrdiff_t* d_offsets, const std::ptrdiff_t n) {
    std::ptrdiff_t total = 0;
    CHECK(cudaMemcpy(&total, d_offsets + n, sizeof(total), cudaMemcpyDeviceToHost));
    return total;
}

}  // namespace

int main() {
    int device_count = 0;
    if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count == 0) {
        std::printf("no CUDA device available; nothing to demonstrate\n");
        return 0;
    }

    const Scene s = make_scene();

    const std::ptrdiff_t n_nodes = (std::ptrdiff_t)s.x0.size();
    const std::ptrdiff_t n_faces = (std::ptrdiff_t)s.f0.size();
    const std::ptrdiff_t n_edges = (std::ptrdiff_t)s.e0.size();

    // ---- upload -------------------------------------------------------------
    // Coordinates are structure-of-arrays: one row per axis, exactly as on the
    // host. The library has no container and no mesh type on either side.
    T* p0_rows[kDim] = {upload(s.x0), upload(s.y0), upload(s.z0)};
    T* p1_rows[kDim] = {upload(s.x1), upload(s.y1), upload(s.z1)};
    I* face_rows[3] = {upload(s.f0), upload(s.f1), upload(s.f2)};
    I* edge_rows[2] = {upload(s.e0), upload(s.e1)};

    T** d_p0 = upload_table(p0_rows, kDim);
    T** d_p1 = upload_table(p1_rows, kDim);
    I** d_faces = upload_table(face_rows, 3);
    I** d_edges = upload_table(edge_rows, 2);

    // ---- broad phase, step 1: swept AABBs -----------------------------------
    // One box per node, face and edge, enclosing the whole trajectory over the
    // step: start position and end position together.
    DeviceAabbs vaabb, faabb, eaabb;
    vaabb.allocate(n_nodes);
    faabb.allocate(n_faces);
    eaabb.allocate(n_edges);

    // BoxRounding::OutwardUlp rounds each bound outward by one ULP. It is a
    // no-op here, where the boxes have the same type as the coordinates, and it
    // is what keeps the broad phase from dropping a candidate when they do not.
    const auto rounding = sccd::BoxRounding::OutwardUlp;
    sccd::device::compute_aabbs(n_nodes, d_p0, d_p1, vaabb.table, rounding);
    sccd::device::compute_aabbs(3, n_faces, d_faces, d_p0, d_p1, faabb.table, rounding);
    sccd::device::compute_aabbs(2, n_edges, d_edges, d_p0, d_p1, eaabb.table, rounding);

    // ---- broad phase, step 2: sort along one axis ---------------------------
    // choose_axis picks the axis with the largest spread of box centres, which
    // is where a sweep prunes best. It returns a host int, having reduced on the
    // device.
    const int sort_axis = sccd::device::choose_axis<T>(n_nodes, vaabb.table);
    T* scratch = device_zeros<T>((std::size_t)std::max(n_nodes, std::max(n_faces, n_edges)));
    sccd::device::sort_along_axis(n_nodes, sort_axis, vaabb.table, vaabb.idx, scratch);
    sccd::device::sort_along_axis(n_faces, sort_axis, faabb.table, faabb.idx, scratch);
    sccd::device::sort_along_axis(n_edges, sort_axis, eaabb.table, eaabb.idx, scratch);

    // ---- broad phase, step 3: candidate pairs -------------------------------
    // Vertex-face and edge-edge are different queries and get separate sweeps.
    // Both are count-then-fill: the count pass writes prefix offsets, so the
    // exact output size is known before a single pair is written and the fill
    // pass needs no reallocation and no per-thread growable buffer. On a GPU
    // that matters more than on a CPU -- it is what lets the fill pass be a
    // plain parallel write with no atomics and no compaction.
    //
    // A pair sharing a vertex is not a candidate, and the sweep masks those out
    // itself given the connectivity, which is why the element arrays go to both
    // passes.

    // Faces against vertices: the second list is points, hence <3, 1> and a null
    // element array with stride 0. This form needs the running maximum of the
    // second list's upper bounds. `soa_device_row` reads one row pointer out of
    // a device table -- the table lives on the device, so its entries cannot
    // simply be indexed from host code.
    T* d_cummax = device_zeros<T>((std::size_t)n_nodes);
    T* vaabb_max_axis = sccd::device::soa_device_row<T>(vaabb.table, kDim + sort_axis);
    sccd::device::cummax(n_nodes, vaabb_max_axis, d_cummax);

    std::ptrdiff_t* d_f_offsets = device_zeros<std::ptrdiff_t>((std::size_t)n_faces + 1);
    sccd::device::count_overlaps<3, 1, T, I>(sort_axis, n_faces, faabb.table, faabb.idx, 1, d_faces,
                                             n_nodes, vaabb.table, vaabb.idx, 0, (I**)nullptr,
                                             d_f_offsets, d_cummax);
    const std::ptrdiff_t n_vf = read_total(d_f_offsets, n_faces);

    I* d_vf_face = nullptr;
    I* d_vf_vertex = nullptr;
    if (n_vf > 0) {
        CHECK(cudaMalloc(&d_vf_face, sizeof(I) * (std::size_t)n_vf));
        CHECK(cudaMalloc(&d_vf_vertex, sizeof(I) * (std::size_t)n_vf));
        sccd::device::collect_overlaps<3, 1, T, I>(sort_axis, n_faces, faabb.table, faabb.idx, 1,
                                                   d_faces, n_nodes, vaabb.table, vaabb.idx, 0,
                                                   (I**)nullptr, d_f_offsets, d_cummax, d_vf_face,
                                                   d_vf_vertex);
    }

    // Edges against edges: one list against itself, so the self variant, which
    // visits each unordered pair once.
    std::ptrdiff_t* d_e_offsets = device_zeros<std::ptrdiff_t>((std::size_t)n_edges + 1);
    sccd::device::count_self_overlaps<2, T, I>(sort_axis, n_edges, eaabb.table, eaabb.idx, 1,
                                               d_edges, d_e_offsets);
    const std::ptrdiff_t n_ee = read_total(d_e_offsets, n_edges);

    I* d_ee_a = nullptr;
    I* d_ee_b = nullptr;
    if (n_ee > 0) {
        CHECK(cudaMalloc(&d_ee_a, sizeof(I) * (std::size_t)n_ee));
        CHECK(cudaMalloc(&d_ee_b, sizeof(I) * (std::size_t)n_ee));
        sccd::device::collect_self_overlaps<2, T, I>(sort_axis, n_edges, eaabb.table, eaabb.idx, 1,
                                                     d_edges, d_e_offsets, d_ee_a, d_ee_b);
    }

    std::printf("broad phase: %td vertex-face and %td edge-edge candidate pairs\n", n_vf, n_ee);

    // The broad phase may over-report -- a candidate is a pair whose swept boxes
    // overlap, not a contact -- but it must not drop the pair that collides.
    // Node 3 against face 0 is the one that does.
    std::vector<I> vf_face((std::size_t)n_vf), vf_vertex((std::size_t)n_vf);
    if (n_vf > 0) {
        CHECK(cudaMemcpy(vf_face.data(), d_vf_face, sizeof(I) * (std::size_t)n_vf,
                         cudaMemcpyDeviceToHost));
        CHECK(cudaMemcpy(vf_vertex.data(), d_vf_vertex, sizeof(I) * (std::size_t)n_vf,
                         cudaMemcpyDeviceToHost));
    }
    bool found = false;
    for (std::size_t i = 0; i < vf_vertex.size(); ++i) {
        found |= (vf_vertex[i] == 3 && vf_face[i] == 0);
    }
    if (!found) {
        std::printf("\nFAIL: broad phase dropped the colliding vertex-face pair\n");
        return 1;
    }

    // ---- narrow phase -------------------------------------------------------
    // ToiOutput::Earliest asks for a single earliest time of impact over all
    // pairs, which lets every query prune against the running minimum -- shared
    // here through one device scalar, updated atomically. Pass ToiOutput::PerPair
    // instead, with a device array of one element per pair, to get them
    // individually.
    //
    // The kernel initialises the output itself, except when there are no
    // candidates at all: it returns early, leaving the buffer untouched. Seed it
    // so an empty candidate list reads as "no collision" rather than as garbage.
    const T max_toi = 1.0;
    const int max_depth = 69;
    const T tol = 3e-8;

    T* d_vf_toi = nullptr;
    T* d_ee_toi = nullptr;
    CHECK(cudaMalloc(&d_vf_toi, sizeof(T)));
    CHECK(cudaMalloc(&d_ee_toi, sizeof(T)));
    CHECK(cudaMemcpy(d_vf_toi, &max_toi, sizeof(T), cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_ee_toi, &max_toi, sizeof(T), cudaMemcpyHostToDevice));

    sccd::device::narrow_phase_vf<T, I>((std::size_t)n_vf, d_vf_vertex, d_vf_face, d_p0, d_p1,
                                        /*element_stride=*/1, d_faces, max_toi, d_vf_toi, max_depth,
                                        tol, sccd::ToiOutput::Earliest);

    sccd::device::narrow_phase_ee<T, I>((std::size_t)n_ee, d_ee_a, d_ee_b, d_p0, d_p1,
                                        /*element_stride=*/1, d_edges, max_toi, d_ee_toi, max_depth,
                                        tol, sccd::ToiOutput::Earliest);

    CHECK(cudaDeviceSynchronize());

    T vf_toi = max_toi;
    T ee_toi = max_toi;
    CHECK(cudaMemcpy(&vf_toi, d_vf_toi, sizeof(T), cudaMemcpyDeviceToHost));
    CHECK(cudaMemcpy(&ee_toi, d_ee_toi, sizeof(T), cudaMemcpyDeviceToHost));

    const T toi = std::min(vf_toi, ee_toi);

    // ---- verify -------------------------------------------------------------
    std::printf("narrow phase: vertex-face %.9f, edge-edge %.9f\n", vf_toi, ee_toi);

    if (toi >= max_toi) {
        std::printf("\nFAIL: no collision reported, but this scene has one\n");
        return 1;
    }

    std::printf("\nearliest time of impact  %.9f\n", toi);
    std::printf("exact answer             %.9f\n", kExactToi);

    // The guarantee is one-sided and is the same guarantee on both processors:
    // the reported time may be at or before the true one, never after. Late
    // would let a solver step through the contact.
    if (toi > kExactToi) {
        std::printf("\nFAIL: later than the true contact by %.3e\n", toi - kExactToi);
        return 1;
    }

    std::printf("early by                 %.3e  (the safe direction, as guaranteed)\n",
                kExactToi - toi);
    return 0;
}
