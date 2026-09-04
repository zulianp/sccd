// The device narrow phase shipped without a test. src/tests/cuda/ held a broad
// phase test and a mesh test, so the three device root finders -- vertex-face,
// edge-edge and vertex-quad -- were covered by nothing that ctest runs.
//
// The first run of this file found that `narrow_phase_vq` dereferenced its
// `v0`, `v1` and `quads` arguments on the host, though they are device arrays of
// device pointers -- the convention `narrow_phase_vf` and `narrow_phase_ee` take,
// and what the CCD integration hands all three. It segfaulted on the first call.
// Nothing had reached that entry point before: the quad device kernel is new, and
// the assessment's quad device rows predate it.
//
// ## Why the queries are built rather than sampled
//
// Comparing the device against the host would only be a parity check, and parity
// is the wrong property: two kernels that both drop the same root agree. The
// invariant is signed and absolute -- a reported time of impact must be at or
// before the true one, never after -- so the test needs the true one.
//
// Every query here therefore has a closed-form time of impact. The primitive
// stays in the z = 0 plane for the whole step and the vertex (or the second
// edge) descends through that plane, monotonically, arriving exactly on the
// primitive at a chosen t*. Since z is linear in t and strictly decreasing,
// z(t) = 0 has one root, so t* is not merely *a* contact time, it is *the* first
// one, and "earlier than t*" and "later than t*" both mean something exact.
//
// ## Why the coordinates are not round numbers
//
// The first version of this test built everything from dyadic values so that the
// reference time of impact was exact. It passed against a kernel with a known
// unsound rejection, and the reason is worth recording: with dyadic inputs of
// small mantissa the corner evaluations are themselves exact, so a box holding a
// root has fmin <= 0 exactly and the width of the rejection pad never decides
// anything. An exact test case cannot exercise a numerical guard.
//
// Coordinates therefore use full mantissas, and the reference comes from solving
// the crossing in `long double` and rounding the result *up* to double. The
// primitive lies in z = 0 and the vertex is strictly interior laterally, so
// contact happens exactly when z reaches zero and the crossing time is
// `pz / (pz - qz)`, one division in a wider format. Rounding up keeps the
// reference an upper bound on the true root, which is the direction that makes
// `toi <= reference` a sound statement of the invariant.
//
// A quarter of the cases are grazing: the vertex descends with a very small
// vertical velocity, so the root is ill-conditioned and the function values near
// it are the size of the rounding error rather than far above it. That is the
// regime where the pad is the deciding quantity.
//
// ## What is asserted
//
//   * a collision must be reported -- a miss is a false negative, illegal;
//   * the reported time of impact must satisfy toi <= t*.
//
// Nothing is asserted about how much earlier than t* a kernel reports. That is
// accuracy, not safety, and the two are printed separately so a regression in
// one is not mistaken for the other. False positives are not failures either,
// for the same reason.
//
// ## The tight-tolerance pass
//
// The second pass runs every query again at tol = 1e-16, below the certified
// numerical error bound, which in double is at most
// (vf ? 30 : 28) * eps * min(max_coord, 1)^3 ~ 6.7e-15. That is the regime the
// device's mode-0 rejection used to get wrong: it padded the origin-containment
// test with the caller's tolerance instead of the bound, so a pad of 1e-16 was
// narrower than the error in the corner values it was testing and a box holding
// a root could be discarded. At the usual 3e-8 the pad was four hundred thousand
// times *wider* than it needed to be, which is why no scene ever showed it.
//
// Be clear about what this pass does and does not establish: the pre-fix kernel
// **passes** it. Reverting the pad to `tol` alone and re-running changes nothing,
// here or with the grazing cases, so this is coverage of the regime rather than a
// regression test for that fix. Losing a root needs the computed fmin to land
// above the pad while the true one is at or below zero, which needs a box small
// enough that F is within rounding error of zero across it -- and the acceptance
// conditions fire well before that. Recorded so the next person does not read a
// green run as proof the narrow pad was harmless, or spend another afternoon
// trying to construct the case.

#include "sccd_narrowphase.cuh"
#include "sccd_narrowphase_vq.cuh"

#include "sccd_narrowphase.hpp"
#include "sccd_narrowphase_quad.hpp"
#include "sccd_narrowphase_mode.hpp"

#include <cuda_runtime.h>

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

using scalar_t = double;
using idx_t = std::int32_t;

namespace {

#define SCCD_CUDA_CHECK(call)                                                                                  \
    do {                                                                                             \
        const cudaError_t err_ = (call);                                                             \
        if (err_ != cudaSuccess) {                                                                   \
            std::fprintf(                                                                            \
                stderr, "CUDA %s at %s:%d\n", cudaGetErrorString(err_), __FILE__, __LINE__);         \
            std::exit(2);                                                                            \
        }                                                                                            \
    } while (0)

    constexpr scalar_t kMaxToi = scalar_t(1);
    constexpr int kMaxDepth = 69;

    /**
     * \brief Deterministic full-mantissa source.
     *
     * A fixed sequence rather than a seeded RNG so a failure is reproducible from
     * the case index alone. `real` fills the mantissa deliberately: round numbers
     * make the kernels' corner evaluations exact, and an exact case cannot
     * exercise a numerical guard.
     */
    struct Rng {
        std::uint64_t s = 0x9e3779b97f4a7c15ull;
        std::uint64_t next() {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            return s * 0x2545f4914f6cdd1dull;
        }
        int next_int(const int lo, const int hi) {
            return lo + (int)((next() >> 24) % (std::uint64_t)(hi - lo + 1));
        }
        /** Uniform in [lo, hi) with all 53 mantissa bits populated. */
        scalar_t real(const scalar_t lo, const scalar_t hi) {
            const scalar_t u = (scalar_t)(next() >> 11) * (scalar_t(1) / scalar_t(9007199254740992.0));
            return lo + u * (hi - lo);
        }
    };

    /**
     * \brief The crossing time of a linear z, solved wider and rounded up.
     *
     * The primitive is in z = 0 for the whole step, so contact is exactly when the
     * moving point's z reaches zero: t = pz / (pz - qz). Solving in `long double`
     * and rounding the result up to double keeps it an upper bound on the true
     * root, which is what makes `toi <= reference` a sound reading of the
     * invariant rather than a coin flip on the last bit.
     */
    scalar_t crossing_upper_bound(const scalar_t pz, const scalar_t qz) {
        const long double t = (long double)pz / ((long double)pz - (long double)qz);
        scalar_t r = (scalar_t)t;
        if ((long double)r < t) r = std::nextafter(r, std::numeric_limits<scalar_t>::infinity());
        return r;
    }

    /**
     * \brief A set of queries in SoA layout, with the exact time of impact of each.
     *
     * `nxe` is 3 for a triangle, 4 for a quad, and 2 for an edge (edge-edge uses
     * two rows and indexes both operands into the same element array).
     */
    struct Scene {
        std::vector<scalar_t> p0[3];
        std::vector<scalar_t> p1[3];
        std::vector<idx_t> elem[4];
        int nxe = 3;
        std::vector<idx_t> q0, q1;
        /** Upper bound on each query's true time of impact, from the wider solve. */
        std::vector<scalar_t> t_star;

        idx_t add_vertex(const scalar_t a[3], const scalar_t b[3]) {
            const idx_t i = (idx_t)p0[0].size();
            for (int d = 0; d < 3; ++d) {
                p0[d].push_back(a[d]);
                p1[d].push_back(b[d]);
            }
            return i;
        }
        idx_t add_element(const idx_t* nodes) {
            const idx_t e = (idx_t)elem[0].size();
            for (int v = 0; v < nxe; ++v) elem[v].push_back(nodes[v]);
            return e;
        }
        size_t nq() const { return q0.size(); }
    };

    /**
     * \brief Aim a descending point at \p x, arriving near \p t_star.
     *
     * The primitive lies in z = 0, so the point starts above it and ends below:
     * `P1 = P0 + (X - P0) / t_star`. z is then linear and strictly decreasing
     * with a single zero, which is the contact. Where that zero falls is read
     * back from the constructed endpoints by `crossing_upper_bound` rather than
     * assumed to be t_star -- the division above rounds, and a reference the
     * construction merely intended is not a reference.
     */
    void descend_onto(const scalar_t x[3],
                      const scalar_t p0_in[3],
                      const scalar_t t_star,
                      scalar_t p0_out[3],
                      scalar_t p1_out[3]) {
        const scalar_t inv = scalar_t(1) / t_star;
        for (int d = 0; d < 3; ++d) {
            p0_out[d] = p0_in[d];
            p1_out[d] = p0_in[d] + (x[d] - p0_in[d]) * inv;
        }
    }

    const scalar_t kTStars[] = {
        scalar_t(0.5), scalar_t(0.25), scalar_t(0.125), scalar_t(0.0625), scalar_t(0.03125)};
    constexpr int kNTStars = (int)(sizeof(kTStars) / sizeof(kTStars[0]));

    /** \brief Vertex-face: a triangle translating inside z = 0, a vertex dropping
     *         through a strictly interior point of it. */
    Scene make_vf(const int ncases) {
        Scene s;
        s.nxe = 3;
        Rng r;
        for (int c = 0; c < ncases; ++c) {
            const scalar_t t_star = kTStars[c % kNTStars];
            const bool grazing = (c % 4) == 0;

            // A non-degenerate triangle in the plane, plus a rigid in-plane
            // translation over the step.
            const scalar_t dx = r.real(-1, 1);
            const scalar_t dyy = r.real(-1, 1);
            scalar_t tri0[3][3], tri1[3][3];
            tri0[0][0] = r.real(-1, 1);
            tri0[0][1] = r.real(-1, 1);
            tri0[1][0] = tri0[0][0] + r.real(0.3, 1.0);
            tri0[1][1] = tri0[0][1] - r.real(0.3, 1.0);
            tri0[2][0] = tri0[0][0] + r.real(0.3, 1.0);
            tri0[2][1] = tri0[0][1] + r.real(0.3, 1.0);
            for (int k = 0; k < 3; ++k) {
                tri0[k][2] = 0;
                tri1[k][0] = tri0[k][0] + dx;
                tri1[k][1] = tri0[k][1] + dyy;
                tri1[k][2] = 0;
            }

            // A strictly interior barycentric point, well away from the edges so
            // a last-bit wobble cannot move the contact off the triangle.
            const scalar_t w0 = r.real(0.2, 0.4);
            const scalar_t w1 = r.real(0.2, 0.4);
            const scalar_t w2 = scalar_t(1) - w0 - w1;

            scalar_t x[3];
            for (int d = 0; d < 3; ++d) {
                const scalar_t a = tri0[0][d] + t_star * (tri1[0][d] - tri0[0][d]);
                const scalar_t b = tri0[1][d] + t_star * (tri1[1][d] - tri0[1][d]);
                const scalar_t cc = tri0[2][d] + t_star * (tri1[2][d] - tri0[2][d]);
                x[d] = w0 * a + w1 * b + w2 * cc;
            }

            // Grazing cases descend with a small vertical velocity, so the root
            // is ill-conditioned and the corner values near it are far closer to
            // the rounding error than in the well-conditioned ones. The starting
            // height stays above the largest distance tolerance used below: a
            // point that begins within tol of the plane is already touching, and
            // a kernel reporting t = 0 for it is right rather than inaccurate.
            const scalar_t h0 = grazing ? r.real(1e-6, 1e-4) : r.real(0.05, 1.0);
            scalar_t seed[3] = {x[0] + r.real(-0.3, 0.3), x[1] + r.real(-0.3, 0.3), h0};
            scalar_t v0[3], v1[3];
            descend_onto(x, seed, t_star, v0, v1);

            idx_t nodes[3];
            for (int k = 0; k < 3; ++k) nodes[k] = s.add_vertex(tri0[k], tri1[k]);
            const idx_t f = s.add_element(nodes);
            const idx_t v = s.add_vertex(v0, v1);

            s.q0.push_back(v);
            s.q1.push_back(f);
            s.t_star.push_back(crossing_upper_bound(v0[2], v1[2]));
        }
        return s;
    }

    /** \brief Edge-edge: one edge along x in z = 0, the other along y descending
     *         through it. They are coplanar only at the instant z reaches 0. */
    Scene make_ee(const int ncases) {
        Scene s;
        s.nxe = 2;
        Rng r;
        for (int c = 0; c < ncases; ++c) {
            const scalar_t t_star = kTStars[c % kNTStars];
            const bool grazing = (c % 4) == 0;

            const scalar_t cx = r.real(-1, 1);
            const scalar_t cy = r.real(-1, 1);
            const scalar_t half_a = r.real(0.5, 1.5);
            const scalar_t half_b = r.real(0.5, 1.5);

            // Edge A: along x, stationary in the plane.
            const scalar_t a0[3] = {cx - half_a, cy, 0};
            const scalar_t a1[3] = {cx + half_a, cy, 0};

            // Edge B: along y, sliding a little in x and dropping from h0 > 0
            // through the plane. The two are coplanar only when z reaches zero.
            const scalar_t slide = r.real(-0.2, 0.2);
            const scalar_t h0 = grazing ? r.real(1e-6, 1e-4) : r.real(0.05, 1.0);
            const scalar_t h1 = h0 + (scalar_t(0) - h0) / t_star;

            const scalar_t b0_start[3] = {cx, cy - half_b, h0};
            const scalar_t b1_start[3] = {cx, cy + half_b, h0};
            const scalar_t b0_end[3] = {cx + slide, cy - half_b, h1};
            const scalar_t b1_end[3] = {cx + slide, cy + half_b, h1};

            idx_t na[2] = {s.add_vertex(a0, a0), s.add_vertex(a1, a1)};
            const idx_t ea = s.add_element(na);
            idx_t nb[2] = {s.add_vertex(b0_start, b0_end), s.add_vertex(b1_start, b1_end)};
            const idx_t eb = s.add_element(nb);

            s.q0.push_back(ea);
            s.q1.push_back(eb);
            s.t_star.push_back(crossing_upper_bound(h0, h1));
        }
        return s;
    }

    /** \brief Vertex-quad: a sheared but planar quad in z = 0, a vertex dropping
     *         through a strictly interior bilinear point of it. */
    Scene make_vq(const int ncases) {
        Scene s;
        s.nxe = 4;
        Rng r;
        for (int c = 0; c < ncases; ++c) {
            const scalar_t t_star = kTStars[c % kNTStars];
            const bool grazing = (c % 4) == 0;

            const scalar_t ox = r.real(-1, 1);
            const scalar_t oy = r.real(-1, 1);
            const scalar_t side = r.real(0.5, 1.5);
            const scalar_t shear = r.real(-0.25, 0.25);
            const scalar_t dx = r.real(-1, 1);
            const scalar_t dyy = r.real(-1, 1);

            // Lexicographic order, not cyclic: the vertex-quad form weights node
            // k by w1=(1-u)(1-v), w2=u(1-v), w3=(1-u)v, w4=uv, so node 3 is the
            // (0,1) corner and node 4 is (1,1). Winding the nodes round the quad
            // instead puts the surface somewhere the kernel is not looking, and
            // the queries then legitimately have no contact. All in z = 0.
            scalar_t quad0[4][3] = {{ox, oy, 0},
                                    {ox + side, oy, 0},
                                    {ox + shear, oy + side, 0},
                                    {ox + side + shear, oy + side, 0}};
            scalar_t quad1[4][3];
            for (int k = 0; k < 4; ++k) {
                quad1[k][0] = quad0[k][0] + dx;
                quad1[k][1] = quad0[k][1] + dyy;
                quad1[k][2] = 0;
            }

            // Strictly interior bilinear parameters, away from the boundary.
            const scalar_t u = r.real(0.25, 0.75);
            const scalar_t v = r.real(0.25, 0.75);

            scalar_t x[3];
            for (int d = 0; d < 3; ++d) {
                scalar_t q[4];
                for (int k = 0; k < 4; ++k) q[k] = quad0[k][d] + t_star * (quad1[k][d] - quad0[k][d]);
                const scalar_t w1 = (scalar_t(1) - u) * (scalar_t(1) - v);
                const scalar_t w2 = u * (scalar_t(1) - v);
                const scalar_t w3 = (scalar_t(1) - u) * v;
                const scalar_t w4 = u * v;
                x[d] = w1 * q[0] + w2 * q[1] + w3 * q[2] + w4 * q[3];
            }

            const scalar_t h0 = grazing ? r.real(1e-6, 1e-4) : r.real(0.05, 1.0);
            scalar_t seed[3] = {x[0] + r.real(-0.3, 0.3), x[1] + r.real(-0.3, 0.3), h0};
            scalar_t p0[3], p1[3];
            descend_onto(x, seed, t_star, p0, p1);

            idx_t nodes[4];
            for (int k = 0; k < 4; ++k) nodes[k] = s.add_vertex(quad0[k], quad1[k]);
            const idx_t q = s.add_element(nodes);
            const idx_t vtx = s.add_vertex(p0, p1);

            s.q0.push_back(vtx);
            s.q1.push_back(q);
            s.t_star.push_back(crossing_upper_bound(p0[2], p1[2]));
        }
        return s;
    }

    /**
     * \brief The first \p n queries of \p s, sharing its geometry.
     *
     * Dropping queries cannot invalidate anything: the vertex and element arrays
     * are indexed by the queries, not the other way round. Used to keep the
     * tight-tolerance pass to a test's runtime -- the search depth needed to
     * satisfy a 1e-16 tolerance is what makes it expensive, and that cost is per
     * query, so fewer queries buys the whole saving without weakening what each
     * one checks.
     */
    Scene head(const Scene& s, const size_t n) {
        Scene o = s;
        const size_t k = n < s.nq() ? n : s.nq();
        o.q0.resize(k);
        o.q1.resize(k);
        o.t_star.resize(k);
        return o;
    }

    // -----------------------------------------------------------------------
    // Device buffers
    // -----------------------------------------------------------------------

    template <typename T>
    T* dup(const std::vector<T>& h) {
        T* d = nullptr;
        SCCD_CUDA_CHECK(cudaMalloc(&d, sizeof(T) * (h.empty() ? 1 : h.size())));
        if (!h.empty()) SCCD_CUDA_CHECK(cudaMemcpy(d, h.data(), sizeof(T) * h.size(), cudaMemcpyHostToDevice));
        return d;
    }

    struct DeviceScene {
        scalar_t* rows0[3] = {};
        scalar_t* rows1[3] = {};
        idx_t* elem_rows[4] = {};
        scalar_t** d_p0 = nullptr;
        scalar_t** d_p1 = nullptr;
        idx_t** d_elem = nullptr;
        idx_t* d_q0 = nullptr;
        idx_t* d_q1 = nullptr;
        scalar_t* d_toi = nullptr;
        int nxe = 0;

        void upload(const Scene& s) {
            nxe = s.nxe;
            for (int d = 0; d < 3; ++d) {
                rows0[d] = dup(s.p0[d]);
                rows1[d] = dup(s.p1[d]);
            }
            SCCD_CUDA_CHECK(cudaMalloc(&d_p0, sizeof(scalar_t*) * 3));
            SCCD_CUDA_CHECK(cudaMalloc(&d_p1, sizeof(scalar_t*) * 3));
            SCCD_CUDA_CHECK(cudaMemcpy(d_p0, rows0, sizeof(scalar_t*) * 3, cudaMemcpyHostToDevice));
            SCCD_CUDA_CHECK(cudaMemcpy(d_p1, rows1, sizeof(scalar_t*) * 3, cudaMemcpyHostToDevice));

            for (int v = 0; v < nxe; ++v) elem_rows[v] = dup(s.elem[v]);
            SCCD_CUDA_CHECK(cudaMalloc(&d_elem, sizeof(idx_t*) * 4));
            SCCD_CUDA_CHECK(cudaMemcpy(d_elem, elem_rows, sizeof(idx_t*) * (size_t)nxe, cudaMemcpyHostToDevice));

            d_q0 = dup(s.q0);
            d_q1 = dup(s.q1);
            SCCD_CUDA_CHECK(cudaMalloc(&d_toi, sizeof(scalar_t) * s.nq()));
        }

        void free_all() {
            for (int d = 0; d < 3; ++d) {
                cudaFree(rows0[d]);
                cudaFree(rows1[d]);
            }
            for (int v = 0; v < nxe; ++v) cudaFree(elem_rows[v]);
            cudaFree(d_p0);
            cudaFree(d_p1);
            cudaFree(d_elem);
            cudaFree(d_q0);
            cudaFree(d_q1);
            cudaFree(d_toi);
        }
    };

    // -----------------------------------------------------------------------
    // Checking
    // -----------------------------------------------------------------------

    struct Verdict {
        int missed = 0;
        int late = 0;
        int hit = 0;
        double worst_late = 0.0;   // largest (toi - t*) over the late queries
        double worst_early = 0.0;  // largest (t* - toi) over the on-time queries
    };

    Verdict check(const Scene& s, const std::vector<scalar_t>& toi) {
        Verdict v;
        for (size_t i = 0; i < s.nq(); ++i) {
            const double t = (double)toi[i];
            const double ref = (double)s.t_star[i];  // upper bound on the true root
            if (!(t < (double)kMaxToi)) {
                ++v.missed;
                continue;
            }
            ++v.hit;
            if (t > ref) {
                ++v.late;
                if (t - ref > v.worst_late) v.worst_late = t - ref;
            } else if (ref - t > v.worst_early) {
                v.worst_early = ref - t;
            }
        }
        return v;
    }

    int report(const char* label, const Scene& s, const Verdict& v) {
        const bool bad = (v.missed != 0) || (v.late != 0);
        std::printf("  %-34s queries=%-5zu hit=%-5d missed=%-4d late=%-4d worst_late=%-10.3e max_early=%.3e  %s\n",
                    label,
                    s.nq(),
                    v.hit,
                    v.missed,
                    v.late,
                    v.worst_late,
                    v.worst_early,
                    bad ? "FAIL" : "ok");
        return bad ? 1 : 0;
    }

    void set_mode(const int mode) {
        const std::string m = std::to_string(mode);
        setenv("SCCD_NARROWPHASE_MODE", m.c_str(), 1);
    }

    // -----------------------------------------------------------------------
    // Runners
    // -----------------------------------------------------------------------

    std::vector<scalar_t> host_vf(const Scene& s, const scalar_t tol) {
        std::vector<scalar_t> toi(s.nq(), kMaxToi);
        scalar_t* p0[3] = {(scalar_t*)s.p0[0].data(), (scalar_t*)s.p0[1].data(), (scalar_t*)s.p0[2].data()};
        scalar_t* p1[3] = {(scalar_t*)s.p1[0].data(), (scalar_t*)s.p1[1].data(), (scalar_t*)s.p1[2].data()};
        idx_t* el[3] = {(idx_t*)s.elem[0].data(), (idx_t*)s.elem[1].data(), (idx_t*)s.elem[2].data()};
        sccd::narrow_phase_vf<3, scalar_t, idx_t>(
            s.nq(), s.q0.data(), s.q1.data(), p0, p1, 1, el, kMaxToi, toi.data(), kMaxDepth, tol, 1);
        return toi;
    }

    std::vector<scalar_t> host_ee(const Scene& s, const scalar_t tol) {
        std::vector<scalar_t> toi(s.nq(), kMaxToi);
        scalar_t* p0[3] = {(scalar_t*)s.p0[0].data(), (scalar_t*)s.p0[1].data(), (scalar_t*)s.p0[2].data()};
        scalar_t* p1[3] = {(scalar_t*)s.p1[0].data(), (scalar_t*)s.p1[1].data(), (scalar_t*)s.p1[2].data()};
        idx_t* el[2] = {(idx_t*)s.elem[0].data(), (idx_t*)s.elem[1].data()};
        sccd::narrow_phase_ee<scalar_t, idx_t>(
            s.nq(), s.q0.data(), s.q1.data(), p0, p1, 1, el, kMaxToi, toi.data(), kMaxDepth, tol, 1);
        return toi;
    }

    std::vector<scalar_t> host_vq(const Scene& s, const scalar_t tol) {
        std::vector<scalar_t> toi(s.nq(), kMaxToi);
        scalar_t* p0[3] = {(scalar_t*)s.p0[0].data(), (scalar_t*)s.p0[1].data(), (scalar_t*)s.p0[2].data()};
        scalar_t* p1[3] = {(scalar_t*)s.p1[0].data(), (scalar_t*)s.p1[1].data(), (scalar_t*)s.p1[2].data()};
        idx_t* el[4] = {(idx_t*)s.elem[0].data(),
                        (idx_t*)s.elem[1].data(),
                        (idx_t*)s.elem[2].data(),
                        (idx_t*)s.elem[3].data()};
        sccd::narrow_phase_vq<4, scalar_t, idx_t>(
            s.nq(), s.q0.data(), s.q1.data(), p0, p1, 1, el, kMaxToi, toi.data(), kMaxDepth, tol, 1);
        return toi;
    }

    enum class Kind { VF, EE, VQ };

    std::vector<scalar_t> device_run(const Scene& s, const Kind kind, const scalar_t tol) {
        DeviceScene d;
        d.upload(s);
        switch (kind) {
            case Kind::VF:
                sccd::device::narrow_phase_vf<3, scalar_t, idx_t>(s.nq(),
                                                                  d.d_q0,
                                                                  d.d_q1,
                                                                  d.d_p0,
                                                                  d.d_p1,
                                                                  1,
                                                                  d.d_elem,
                                                                  kMaxToi,
                                                                  d.d_toi,
                                                                  kMaxDepth,
                                                                  tol,
                                                                  1);
                break;
            case Kind::EE:
                sccd::device::narrow_phase_ee<scalar_t, idx_t>(s.nq(),
                                                               d.d_q0,
                                                               d.d_q1,
                                                               d.d_p0,
                                                               d.d_p1,
                                                               1,
                                                               d.d_elem,
                                                               kMaxToi,
                                                               d.d_toi,
                                                               kMaxDepth,
                                                               tol,
                                                               1);
                break;
            case Kind::VQ:
                sccd::device::narrow_phase_vq<scalar_t, idx_t>(s.nq(),
                                                               d.d_q0,
                                                               d.d_q1,
                                                               d.d_p0,
                                                               d.d_p1,
                                                               1,
                                                               d.d_elem,
                                                               kMaxToi,
                                                               d.d_toi,
                                                               kMaxDepth,
                                                               tol,
                                                               1);
                break;
        }
        SCCD_CUDA_CHECK(cudaDeviceSynchronize());
        std::vector<scalar_t> toi(s.nq());
        SCCD_CUDA_CHECK(cudaMemcpy(toi.data(), d.d_toi, sizeof(scalar_t) * s.nq(), cudaMemcpyDeviceToHost));
        d.free_all();
        return toi;
    }

}  // namespace

int main() {
    // Unbuffered: this runs under ctest with stdout redirected, and a crash with
    // a full buffer loses exactly the line that says where it happened.
    setvbuf(stdout, nullptr, _IONBF, 0);

    int devices = 0;
    if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) {
        std::printf("no CUDA device available\n");
        return 0;
    }

    const Scene vf_all = make_vf(600);
    const Scene ee_all = make_ee(600);
    const Scene vq_all = make_vq(400);

    std::printf("vertex-face %zu queries, edge-edge %zu, vertex-quad %zu\n\n",
                vf_all.nq(),
                ee_all.nq(),
                vq_all.nq());

    int bad = 0;

    // 3e-8 is what the benchmarks and the CCD interface use. 1e-16 is below the
    // certified error bound in double (<= 6.7e-15), which is the regime an
    // origin-containment test padded with the caller's tolerance gets wrong.
    for (const scalar_t tol : {scalar_t(3e-8), scalar_t(1e-16)}) {
        const bool tight = tol < scalar_t(1e-14);
        std::printf("tol = %.1e%s\n", (double)tol, tight ? "   (below the certified error bound)" : "");

        const Scene vf = tight ? head(vf_all, 100) : vf_all;
        const Scene ee = tight ? head(ee_all, 100) : ee_all;
        const Scene vq = tight ? head(vq_all, 100) : vq_all;

        for (const int mode : {0, 2}) {
            set_mode(mode);
            char label[128];

            std::snprintf(label, sizeof(label), "host   vertex-face  mode %d", mode);
            std::fflush(stdout);
            bad += report(label, vf, check(vf, host_vf(vf, tol)));
            std::snprintf(label, sizeof(label), "device vertex-face  mode %d", mode);
            bad += report(label, vf, check(vf, device_run(vf, Kind::VF, tol)));

            std::snprintf(label, sizeof(label), "host   edge-edge    mode %d", mode);
            bad += report(label, ee, check(ee, host_ee(ee, tol)));
            std::snprintf(label, sizeof(label), "device edge-edge    mode %d", mode);
            bad += report(label, ee, check(ee, device_run(ee, Kind::EE, tol)));
        }

        // Quads have one root-finder variant on each side and never consult the
        // mode enum, so running them per mode would measure the same code twice.
        unsetenv("SCCD_NARROWPHASE_MODE");
        bad += report("host   vertex-quad", vq, check(vq, host_vq(vq, tol)));
        bad += report("device vertex-quad", vq, check(vq, device_run(vq, Kind::VQ, tol)));

        std::printf("\n");
    }

    if (bad) {
        std::printf("FAILED: %d configuration(s) missed a collision or reported a late time of impact\n",
                    bad);
        return 1;
    }
    std::printf("all configurations conservative\n");
    return 0;
}
