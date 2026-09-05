// SCCD end to end with nothing but the standard library.
//
// Two disconnected triangles, one of which passes through the other during the
// step. Both stages run: swept AABBs and the sweep-and-prune broad phase to get
// candidate pairs, then the narrow phase to turn those into a time of impact.
// The scene is built so the exact answer is known, so the result is verified
// rather than merely printed.
//
// The library has no container type and no mesh type. Geometry is
// structure-of-arrays -- `T** points` is three row pointers, one per axis -- so
// std::vector is all that is needed.
//
// Build (the default configuration; no options, no dependencies):
//
//     cmake -S . -B build && cmake --build build -j --target sccd_minimal
//     ./build/sccd_minimal

#include "sccd_aabb.hpp"
#include "sccd_broadphase_sweep.hpp"
#include "sccd_narrowphase.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <vector>

namespace {

using T = double;
using I = int;

constexpr int kDim = 3;

// The scene.
//
//   Triangle A, stationary, in the plane z = 0:   nodes 0, 1, 2
//   Triangle B, descending by 2 over the step:    nodes 3, 4, 5
//
// The two triangles share no vertex, so every vertex-face and edge-edge pair
// between them is a genuine candidate.
//
// B is tilted, so node 3 is its lowest point and leads. The motion is purely
// vertical, so node 3 stays at (0.5, 0.5) in x-y for the whole step, and that
// point is strictly inside triangle A. Its height runs 1.0 -> -1.0, reaching
// z = 0 at exactly t = 1/2, and every other point of B is above it throughout.
// The first contact is therefore node 3 against face A at
//
//     t = 0.5, exactly.
//
// All coordinates are exactly representable, so "exactly" is meant literally:
// the true root is the binary floating-point number 0.5, not a rounding of it.
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

// A six-row AABB block: rows 0..2 hold the minimum corner, rows 3..5 the
// maximum. That is the layout every broad-phase entry point expects. `idx` maps
// a sorted position back to the element it came from; the sort permutes both.
struct Aabbs {
    std::vector<T> rows_[6];
    T* rows[6];
    std::vector<I> idx;

    void resize(const std::ptrdiff_t n) {
        for (int d = 0; d < 6; ++d) {
            rows_[d].assign(n, T(0));
            rows[d] = rows_[d].data();
        }
        idx.resize(n);
        for (std::ptrdiff_t i = 0; i < n; ++i) idx[i] = static_cast<I>(i);
    }
};

}  // namespace

int main() {
    Scene s = make_scene();

    T* p0[kDim] = {s.x0.data(), s.y0.data(), s.z0.data()};
    T* p1[kDim] = {s.x1.data(), s.y1.data(), s.z1.data()};
    I* faces[3] = {s.f0.data(), s.f1.data(), s.f2.data()};
    I* edges[2] = {s.e0.data(), s.e1.data()};

    const std::ptrdiff_t n_nodes = static_cast<std::ptrdiff_t>(s.x0.size());
    const std::ptrdiff_t n_faces = static_cast<std::ptrdiff_t>(s.f0.size());
    const std::ptrdiff_t n_edges = static_cast<std::ptrdiff_t>(s.e0.size());

    // ---- broad phase, step 1: swept AABBs -----------------------------------
    // One box per node, face and edge, enclosing the whole trajectory over the
    // step: start position and end position together.
    Aabbs vaabb, faabb, eaabb;
    vaabb.resize(n_nodes);
    faabb.resize(n_faces);
    eaabb.resize(n_edges);

    // `safe_inflate` rounds each bound outward by one ULP. It is a no-op here,
    // where the boxes have the same type as the coordinates and the bounds are
    // therefore exact, and it is what keeps the broad phase from dropping a
    // candidate when they do not -- float boxes over double geometry, say.
    const bool safe_inflate = true;
    sccd::compute_aabbs(kDim, n_nodes, p0, p1, vaabb.rows, vaabb.rows + kDim, safe_inflate);
    sccd::compute_aabbs(3, n_faces, faces, kDim, p0, p1, faabb.rows, faabb.rows + kDim, safe_inflate);
    sccd::compute_aabbs(2, n_edges, edges, kDim, p0, p1, eaabb.rows, eaabb.rows + kDim, safe_inflate);

    // ---- broad phase, step 2: sort along one axis ---------------------------
    // The sweep needs the lists ordered on a common axis. choose_axis picks the
    // one with the largest spread of box centres, which is where a sweep prunes
    // best. `scratch` is sized for the longest list.
    const int sort_axis = sccd::choose_axis(n_nodes, vaabb.rows);
    std::vector<T> scratch(std::max(n_nodes, std::max(n_faces, n_edges)));
    sccd::sort_along_axis(n_nodes, sort_axis, vaabb.rows, vaabb.idx.data(), scratch.data());
    sccd::sort_along_axis(n_faces, sort_axis, faabb.rows, faabb.idx.data(), scratch.data());
    sccd::sort_along_axis(n_edges, sort_axis, eaabb.rows, eaabb.idx.data(), scratch.data());

    // ---- broad phase, step 3: candidate pairs -------------------------------
    // Vertex-face and edge-edge are different queries and get separate sweeps.
    // Both are count-then-fill: the count pass writes prefix offsets, so the
    // exact output size is known before a single pair is written and the fill
    // pass needs no reallocation and no per-thread growable buffer.
    //
    // A pair sharing a vertex is not a candidate -- a vertex cannot collide
    // with a face it belongs to -- and the sweep masks those out itself, given
    // the connectivity. That is why the element arrays are passed to both
    // passes.

    // Faces against vertices. The second list is points, hence <3, 1> and a
    // null element array with stride 0. This form needs the running maximum of
    // the second list's upper bounds to find each face's candidate window.
    std::vector<T> cmax(n_nodes);
    sccd::cummax(n_nodes, vaabb.rows[kDim + sort_axis], cmax.data());

    std::vector<std::ptrdiff_t> f_offsets(n_faces + 1, 0);
    sccd::count_overlaps<3, 1, T, I>(sort_axis, n_faces, faabb.rows, faabb.idx.data(), 1, faces,
                                     n_nodes, vaabb.rows, vaabb.idx.data(), 0, nullptr,
                                     f_offsets.data(), cmax.data());

    std::vector<I> vf_face(f_offsets[n_faces]), vf_vertex(f_offsets[n_faces]);
    if (!vf_face.empty()) {
        sccd::collect_overlaps<3, 1, T, I>(sort_axis, n_faces, faabb.rows, faabb.idx.data(), 1, faces,
                                           n_nodes, vaabb.rows, vaabb.idx.data(), 0, nullptr,
                                           f_offsets.data(), cmax.data(), vf_face.data(),
                                           vf_vertex.data());
    }

    // Edges against edges. One list against itself, so the self variant, which
    // visits each unordered pair once.
    std::vector<std::ptrdiff_t> e_offsets(n_edges + 1, 0);
    sccd::count_self_overlaps<2>(sort_axis, n_edges, eaabb.rows, eaabb.idx.data(), 1, edges,
                                 e_offsets.data());

    std::vector<I> ee_a(e_offsets[n_edges]), ee_b(e_offsets[n_edges]);
    if (!ee_a.empty()) {
        sccd::collect_self_overlaps<2>(sort_axis, n_edges, eaabb.rows, eaabb.idx.data(), 1, edges,
                                       e_offsets.data(), ee_a.data(), ee_b.data());
    }

    std::printf("broad phase: %zu vertex-face and %zu edge-edge candidate pairs\n",
                vf_vertex.size(), ee_a.size());

    // The broad phase may over-report -- a candidate is a pair whose swept boxes
    // overlap, not a contact -- but it must not drop the pair that collides. Node
    // 3 against face 0 is the one that does.
    bool found = false;
    for (std::size_t i = 0; i < vf_vertex.size(); ++i) {
        found |= (vf_vertex[i] == 3 && vf_face[i] == 0);
    }
    if (!found) {
        std::printf("\nFAIL: broad phase dropped the colliding vertex-face pair\n");
        return 1;
    }

    // ---- narrow phase -------------------------------------------------------
    // toi_stride = 0 asks for a single earliest time of impact over all pairs,
    // which lets every query prune against the running minimum. Pass 1 instead,
    // with an output array of one element per pair, to get them individually.
    const T max_toi = 1.0;
    const int max_depth = 69;
    const T tol = 3e-8;

    T vf_toi = max_toi;
    sccd::narrow_phase_vf<T, I>(vf_vertex.size(), vf_vertex.data(), vf_face.data(), p0, p1,
                                   /*face_stride=*/1, faces, max_toi, &vf_toi, max_depth, tol,
                                   /*toi_stride=*/0);

    T ee_toi = max_toi;
    sccd::narrow_phase_ee<T, I>(ee_a.size(), ee_a.data(), ee_b.data(), p0, p1,
                                /*edge_stride=*/1, edges, max_toi, &ee_toi, max_depth, tol,
                                /*toi_stride=*/0);

    const T toi = std::min(vf_toi, ee_toi);

    // ---- verify -------------------------------------------------------------
    std::printf("narrow phase: vertex-face %.9f, edge-edge %.9f\n", vf_toi, ee_toi);

    if (toi >= max_toi) {
        std::printf("\nFAIL: no collision reported, but this scene has one\n");
        return 1;
    }

    std::printf("\nearliest time of impact  %.9f\n", toi);
    std::printf("exact answer             %.9f\n", kExactToi);

    // The guarantee is one-sided: the reported time may be at or before the true
    // one, never after. Late would let a solver step through the contact.
    if (toi > kExactToi) {
        std::printf("\nFAIL: later than the true contact by %.3e\n", toi - kExactToi);
        return 1;
    }

    std::printf("early by                 %.3e  (the safe direction, as guaranteed)\n",
                kExactToi - toi);
    return 0;
}
