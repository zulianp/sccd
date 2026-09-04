// SCCD with nothing but the standard library.
//
// The narrow phase takes structure-of-arrays geometry: one contiguous array per
// axis, and one per face vertex slot. That is all `T** v0` means -- three row
// pointers, not a matrix type -- so std::vector is a perfectly good container
// for it and no mesh library is involved.
//
// Build (the default configuration; no options, no dependencies):
//
//     cmake -S . -B build && cmake --build build -j --target sccd_minimal
//     ./build/sccd_minimal
//
// The scene: one triangle sitting in the z = 0 plane, and two vertices moving
// straight down through the step. The first passes through the middle of the
// triangle at t = 0.25; the second travels parallel to the plane and never
// touches it. The exact time of impact of the first is 0.25 by construction, so
// the printed result can be checked against it.

#include "sccd_narrowphase.hpp"

#include <cstdio>
#include <vector>

int main() {
    using T = double;
    using I = int;

    // ---- geometry, start of step -------------------------------------------
    // Nodes 0..2 are the triangle; node 3 and node 4 are the moving vertices.
    // One vector per axis: this is the structure-of-arrays layout the kernels
    // want, and the reason they take T** rather than a point type.
    std::vector<T> x0 = {0.0, 1.0, 0.0, 0.25, 5.0};
    std::vector<T> y0 = {0.0, 0.0, 1.0, 0.25, 5.0};
    std::vector<T> z0 = {0.0, 0.0, 0.0, 1.00, 1.0};

    // ---- geometry, end of step ---------------------------------------------
    // The triangle does not move. Vertex 3 descends from z=1 to z=-3, crossing
    // z=0 at t = 1/4. Vertex 4 slides sideways at constant height and misses.
    std::vector<T> x1 = {0.0, 1.0, 0.0, 0.25, 6.0};
    std::vector<T> y1 = {0.0, 0.0, 1.0, 0.25, 5.0};
    std::vector<T> z1 = {0.0, 0.0, 0.0, -3.0, 1.0};

    T* v0[3] = {x0.data(), y0.data(), z0.data()};
    T* v1[3] = {x1.data(), y1.data(), z1.data()};

    // ---- connectivity ------------------------------------------------------
    // One row per vertex slot, indexed by face. Face 0 is nodes 0, 1, 2.
    std::vector<I> f0 = {0};
    std::vector<I> f1 = {1};
    std::vector<I> f2 = {2};
    I* faces[3] = {f0.data(), f1.data(), f2.data()};

    // ---- candidate pairs ---------------------------------------------------
    // Normally the broad phase produces these. Two queries: vertex 3 against
    // face 0, and vertex 4 against face 0.
    std::vector<I> vertex_of_pair = {3, 4};
    std::vector<I> face_of_pair = {0, 0};

    // ---- run ---------------------------------------------------------------
    // toi_stride = 1 asks for one time of impact per pair. Pass 0 instead and
    // toi[0] comes back holding the single earliest time over all of them.
    const std::size_t n_pairs = vertex_of_pair.size();
    std::vector<T> toi(n_pairs, T(1));

    const T max_toi = 1.0;   // search the whole step
    const int max_depth = 69;
    const T tol = 3e-8;

    sccd::narrow_phase_vf<3, T, I>(n_pairs,
                                   vertex_of_pair.data(),
                                   face_of_pair.data(),
                                   v0,
                                   v1,
                                   /*face_stride=*/1,
                                   faces,
                                   max_toi,
                                   toi.data(),
                                   max_depth,
                                   tol,
                                   /*toi_stride=*/1);

    // ---- results -----------------------------------------------------------
    // A time of impact of 1 means "no collision in this step". Anything below is
    // a contact, reported at or before the true time -- never after.
    const T expected = 0.25;
    std::printf("query  vertex  face  time of impact\n");
    for (std::size_t i = 0; i < n_pairs; ++i) {
        std::printf("  %zu       %d      %d    ", i, vertex_of_pair[i], face_of_pair[i]);
        if (toi[i] < max_toi) {
            std::printf("%.9f\n", toi[i]);
        } else {
            std::printf("no collision\n");
        }
    }

    const bool hit_found = toi[0] < max_toi;
    const bool conservative = hit_found && toi[0] <= expected;
    std::printf("\nexact time of impact for query 0 is %.9f\n", expected);
    std::printf("reported %.9f, which is %s (early by %.3e)\n",
                toi[0],
                conservative ? "at or before it, as guaranteed" : "LATER -- that is a bug",
                expected - toi[0]);

    if (!hit_found || !conservative || toi[1] < max_toi) {
        std::printf("\nunexpected result\n");
        return 1;
    }
    return 0;
}
