#ifndef SCCD_NARROWPHASE_TIGHT_HPP
#define SCCD_NARROWPHASE_TIGHT_HPP

// A vectorized narrow phase that reproduces TightInclusion's semantics exactly,
// for both vertex-face and edge-edge queries.
//
// Why a second vectorized kernel rather than a fix to the first: the fast kernel
// in sccd_vnarrowphase.hpp trades TightInclusion's predicates for cheaper ones
// (it compares codomain widths against domain tolerances, and has two extra
// acceptance paths), which makes it report times of impact slightly *after* the
// true one. That breaks conservativeness -- see benchmark/oracle/README.md. This
// kernel keeps what actually made the fast kernel fast, the lane-packed work
// list, and replaces the mathematics with TightInclusion's:
//
//   * corner arithmetic in TightInclusion's exact expression form, so the
//     30*eps*delta^3 / 28*eps*delta^3 error certificate stays valid;
//   * acceptance on *domain* widths against the domain tolerances, plus the
//     "function box is inside the error box" case, and nothing else;
//   * bisection on the axis with the largest width/tolerance ratio, with
//     TightInclusion's u+v<=1 and t<toi pruning;
//   * exact dyadic box endpoints (see tight_bisect below).
//
// Ordering: TightInclusion pops boxes from a priority queue ordered by ascending
// t lower bound and returns the first box it accepts, which is therefore the
// accepted box with the globally smallest t. This kernel instead explores the
// same tree in any order, keeps the minimum t over all accepted boxes, and
// prunes every box whose t lower bound already exceeds that minimum. The two
// agree: any box TightInclusion would have examined before its answer is also
// examined here and judged by the same predicate, and every box reached only
// afterwards has a t lower bound at or above the answer and cannot lower the
// minimum. Decoupling the result from the traversal order is what allows lanes
// to be packed with boxes from different queries.

#include "sccd_base.hpp"
#include "sccd_narrowphase_mode.hpp"
#include "sccd_math.hpp"
#include "sccd_numerical_error.hpp"
#include "sccd_tolerance.hpp"
#include "sccd_parallel.hpp"

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

// How many boxes are carried in lanes at once.
//
// Not simply "as wide as possible": half of this kernel's time goes on the
// per-lane scalar bookkeeping between vector sweeps (classify, choose a split
// axis, assemble two children), and that cost is linear in the lane count while
// the vector work is not. Measured on NEON, where a vector holds 2 doubles, 4
// lanes beat 8 by 19% on vertex-face and 27% on edge-edge; 16 and 32 were worse
// still. The rule below is "two vector operations per corner evaluation".
//
// The x86 values are the prior, not a measurement -- re-run
// benchmark/ti_oracle and the cloth-funnel probe there before trusting them.
#ifndef SCCD_NARROWPHASE_TIGHT_VSIZE
#if defined(__AVX512F__) || defined(__AVX2__)
#define SCCD_NARROWPHASE_TIGHT_VSIZE 8
#else
#define SCCD_NARROWPHASE_TIGHT_VSIZE 4
#endif
#endif

// Vertex-face and edge-edge do not want the same width. A vertex-face corner is
// the more expensive of the two (an extra multiply-add and one more subtraction
// per corner), so its vector sweeps take longer to amortize the per-lane scalar
// bookkeeping, and it peaks narrower.
//
// Measured on NEON over armadillo-rollers (20k queries), each phase timed in
// its own process: vertex-face runs 297.7 -> 246.0 ms going from 4 lanes to 2
// (17% faster), while edge-edge runs 599.9 -> 494.6 ms going from 2 to 4 (18%
// faster). 8 and 16 were far worse for both.
//
// Time each phase separately when redoing this. Running both in one process
// makes the numbers useless: the edge-edge time tracks the *vertex-face*
// instantiation's lane width by around 20% even under --phase ee, which never
// executes vertex-face code at all. That is code layout, not work, and it is
// large enough to invent or hide a result of this size.
//
// x86 is left on the single value: it has never been measured, and guessing a
// split there would repeat the mistake this comment block already warns about.
#ifndef SCCD_NARROWPHASE_TIGHT_VSIZE_VF
#if defined(__ARM_NEON) && SCCD_NARROWPHASE_TIGHT_VSIZE > 2
#define SCCD_NARROWPHASE_TIGHT_VSIZE_VF (SCCD_NARROWPHASE_TIGHT_VSIZE / 2)
#else
#define SCCD_NARROWPHASE_TIGHT_VSIZE_VF SCCD_NARROWPHASE_TIGHT_VSIZE
#endif
#endif

#ifndef SCCD_NARROWPHASE_TIGHT_VSIZE_EE
#define SCCD_NARROWPHASE_TIGHT_VSIZE_EE SCCD_NARROWPHASE_TIGHT_VSIZE
#endif

#include <cstdio>
#include <cstdlib>

namespace sccd {

#ifdef SCCD_NP_COUNT_BOXES
    // Boxes classified by the host TightInclusion search, for comparison against
    // the device counter of the same name in src/cuda/sccd_narrowphase.cu. Both
    // count the same unit -- one box whose inclusion function was evaluated and
    // whose fate was decided -- so the two numbers are directly comparable, which
    // is the whole point of counting them. Not thread-safe by design: a relaxed
    // atomic on the hot path of a parallel loop would change what it measures, and
    // an approximate count is enough to tell 10 million from 800 million.
    extern unsigned long long g_np_host_boxes;
#ifndef SCCD_NP_HOST_BOX_TICK
#define SCCD_NP_HOST_BOX_TICK() (++::sccd::g_np_host_boxes)
#endif

    // Per-query counts, bucketed by log2 on the way out. 843k queries is too many
    // to print and the mean is the statistic that hides the answer: a mean of 944
    // is a different defect depending on whether every query costs 944 boxes or
    // three queries cost a quarter of a billion each.
    extern unsigned long long g_np_host_hist[24];
    inline void np_hist_add(unsigned long long* hist, unsigned long long n) {
        int b = 0;
        while (n > 1 && b < 23) {
            n >>= 1;
            ++b;
        }
        ++hist[b];
    }
#define SCCD_NP_HOST_PERQ_TICK(q) (++perq_[(q)])
    // Matching the device's g_np_level / g_np_depth_accept, same units, so the two
    // level distributions can be laid side by side.
    extern unsigned long long g_np_host_level[80];
    extern unsigned long long g_np_host_depth_accept;
#else
#define SCCD_NP_HOST_BOX_TICK() ((void)0)
#endif

    SCCD_FP_STRICT_BEGIN

    /**
     * \brief One corner evaluation, in TightInclusion's expression form.
     *
     * For vertex-face, (a0..a3) are the vertex and the three face corners at t0
     * and (b0..b3) the same at t1; for edge-edge they are (ea0, ea1, eb0, eb1).
     * The two share a signature so the traversal can be written once.
     *
     * The expression order matters and is not free to rearrange: TightInclusion's
     * numerical error bound is derived for exactly this sequence of operations,
     * and the acceptance predicate pads by that bound. An algebraically identical
     * but differently rounded form would invalidate the certificate and can lose
     * a root, which is a missed collision.
     */
    /**
     * \brief The (u, v) half of a corner evaluation, given the four points at t.
     *
     * Split out so a caller evaluating several corners that share a t value can
     * compute the points once. The expressions are the ones tight_corner always
     * used, in the same order, so the values are bit-identical either way and
     * the error certificate is unaffected.
     */
    template <bool IsVertexFace, typename T>
    static SCCD_ALWAYS_INLINE T tight_corner_at(const T p0, const T p1, const T p2, const T p3, const T u, const T v) {
        if constexpr (IsVertexFace) {
            // vertex minus the point (u, v) of the triangle
            const T face = (p2 - p1) * u + (p3 - p1) * v + p1;
            return p0 - face;
        } else {
            // point u of edge a minus point v of edge b
            const T va = (p1 - p0) * u + p0;
            const T vb = (p3 - p2) * v + p2;
            return va - vb;
        }
    }

    template <bool IsVertexFace, typename T>
    static SCCD_ALWAYS_INLINE T tight_corner(const T a0,
                                          const T a1,
                                          const T a2,
                                          const T a3,
                                          const T d0,
                                          const T d1,
                                          const T d2,
                                          const T d3,
                                          const T t,
                                          const T u,
                                          const T v) {
        // d_k is (end - start), exactly the subexpression TightInclusion forms
        // as (x_t1 - x_t0); hoisting it out of the corner loop does not change
        // any rounding, since it is the same single subtraction either way.
        const T p0 = d0 * t + a0;
        const T p1 = d1 * t + a1;
        const T p2 = d2 * t + a2;
        const T p3 = d3 * t + a3;

        return tight_corner_at<IsVertexFace, T>(p0, p1, p2, p3, u, v);
    }

    /**
     * \brief Bisect [lo, hi] at its midpoint, reporting whether it was exact.
     *
     * TightInclusion carries endpoints as dyadic rationals (numerator over a
     * power of two) so a midpoint is always exact. Starting from [0, 1], every
     * bisection midpoint is m/2^k, which a double represents exactly while
     * k <= 52 -- deeper than any box these queries reach in practice. Rather
     * than assume that, the split is verified: when the midpoint fails to
     * separate the endpoints the interval has hit the resolution of the format,
     * and the caller must accept the box rather than split it. Accepting reports
     * the box's t lower bound, which is the conservative direction; dropping it
     * would risk losing a root.
     */
    template <typename T>
    static SCCD_ALWAYS_INLINE bool tight_bisect(const T lo, const T hi, T& mid) {
        mid = (lo + hi) * T(0.5);
        return lo < mid && mid < hi;
    }

    /** \brief Geometry for a block of queries, as structure-of-arrays. */
    template <typename T, int VSIZE>
    struct TightQueryBlock {
        // [point][dim][lane]: points are (vertex, f0, f1, f2) or (ea0, ea1, eb0, eb1)
        T start[4][3][VSIZE];
        T delta[4][3][VSIZE];  // end - start, hoisted out of the corner loop
        T tol[3][VSIZE];       // per-axis domain tolerance
        T err[3][VSIZE];       // per-axis numerical error bound
    };

    /**
     * \brief A box on the work list.
     *
     * Carries the function value at each of its 8 corners, per dimension.
     * Bisection makes the two children share the face at the midpoint, and each
     * child inherits the parent's face on its outer side, so a split only has to
     * evaluate the 4 corners on the mid face rather than 8 per child. That is a
     * fourfold cut in the evaluation that dominates this kernel. The inherited
     * values come from the same expression at the same arguments, so nothing
     * about the result changes.
     *
     * Corner c holds the upper bound of axis d when bit (1 << (2 - d)) is set.
     */
    template <typename I>
    struct TightBox {
        I query;
        int depth;
        double lo[3];
        double hi[3];
        double corner[3][8];
    };

    /** \brief Bit of a corner index selecting the upper bound of axis \p d. */
    static SCCD_ALWAYS_INLINE int tight_corner_bit(const int d) { return 1 << (2 - d); }

    /**
     * \brief Corner indices of the mid face, per split axis.
     *
     * tight_mid_corner[axis][k] is the k-th of the 4 corners that have `axis` at its
     * upper bound; k enumerates the other two axes' lo/hi combinations. Tabulated
     * because deriving it inline costs a nested loop with a branch, per lane, per
     * corner, on the hot path.
     */
    static constexpr int tight_mid_corner[3][4] = {{4, 6, 5, 7}, {2, 6, 3, 7}, {1, 5, 3, 7}};

    namespace detail {

        /**
         * \brief Componentwise bounds of the function over a box, per lane.
         *
         * The function is multilinear in (t, u, v), so its exact range over the
         * box is the hull of the eight corner values; the only slack is the
         * floating-point rounding that the error bound covers.
         */
        /**
         * \brief Geometry replicated per lane rather than per query.
         *
         * Lanes hold boxes from different queries, so indexing the block's
         * geometry by the lane's query is a gather on every corner of every
         * dimension. Copying the 30 values a lane needs when it is refilled
         * turns the inner loop into contiguous loads, which is what lets it
         * vectorize; a refill costs far less than the 24 corner evaluations it
         * feeds.
         */
        template <typename T, int VSIZE>
        struct TightLaneGeometry {
            T start[4][3][VSIZE];
            T delta[4][3][VSIZE];
            T tol[3][VSIZE];
            T err[3][VSIZE];
        };

        template <typename T, int VSIZE>
        static SCCD_ALWAYS_INLINE void tight_load_lane(const TightQueryBlock<T, VSIZE>& block,
                                                    const int q,
                                                    const int l,
                                                    TightLaneGeometry<T, VSIZE>& lane) {
            for (int k = 0; k < 4; ++k) {
                for (int d = 0; d < 3; ++d) {
                    lane.start[k][d][l] = block.start[k][d][q];
                    lane.delta[k][d][l] = block.delta[k][d][q];
                }
            }
            for (int d = 0; d < 3; ++d) {
                lane.tol[d][l] = block.tol[d][q];
                lane.err[d][l] = block.err[d][q];
            }
        }

        /** \brief Per-dimension min/max over a lane's 8 stored corner values. */
        template <typename T, int VSIZE>
        static void tight_bounds_from_corners(const T (&corner)[3][8][VSIZE],
                                           T (&fmin)[3][VSIZE],
                                           T (&fmax)[3][VSIZE]) {
            for (int d = 0; d < 3; ++d) {
                T* const SCCD_RESTRICT lo_out = fmin[d];
                T* const SCCD_RESTRICT hi_out = fmax[d];
                SCCD_VECTORIZE_LOOP
                for (int l = 0; l < VSIZE; ++l) {
                    T lo = corner[d][0][l];
                    T hi = corner[d][0][l];
                    for (int c = 1; c < 8; ++c) {
                        const T v = corner[d][c][l];
                        lo = sccd::min<T>(lo, v);
                        hi = sccd::max<T>(hi, v);
                    }
                    lo_out[l] = lo;
                    hi_out[l] = hi;
                }
            }
        }

        /**
         * \brief Evaluate the 4 corners of each lane's mid face, across lanes.
         *
         * The parameter values are precomputed by the caller so that this stays
         * one flat SIMD loop even though each lane may be splitting a different
         * axis.
         */
        template <bool IsVertexFace, typename T, int VSIZE>
        static void tight_eval_mid_face(const TightLaneGeometry<T, VSIZE>& lane,
                                     const T (&mt)[4][VSIZE],
                                     const T (&mu)[4][VSIZE],
                                     const T (&mv)[4][VSIZE],
                                     const uint8_t* const SCCD_RESTRICT wanted,
                                     T (&out)[3][4][VSIZE]) {
            for (int d = 0; d < 3; ++d) {
                // Hoisted into restrict-qualified locals. Without this the
                // compiler cannot prove `out` does not alias the geometry it
                // reads, and it emits the whole loop scalar -- which it did:
                // this kernel's hottest loop had no vector instructions in it at
                // all. #pragma omp simd does not help, since it is inert unless
                // the consumer happens to compile with OpenMP enabled, and the
                // kernel must not depend on that.
                const T* const SCCD_RESTRICT a0 = lane.start[0][d];
                const T* const SCCD_RESTRICT a1 = lane.start[1][d];
                const T* const SCCD_RESTRICT a2 = lane.start[2][d];
                const T* const SCCD_RESTRICT a3 = lane.start[3][d];
                const T* const SCCD_RESTRICT d0 = lane.delta[0][d];
                const T* const SCCD_RESTRICT d1 = lane.delta[1][d];
                const T* const SCCD_RESTRICT d2 = lane.delta[2][d];
                const T* const SCCD_RESTRICT d3 = lane.delta[3][d];

                // Corners k and k+2 of a mid face always share a t value, for
                // every split axis: tight_mid_corner's four entries pair up as
                // (lo, hi, lo, hi) in the t bit, and a t-split pins all four to
                // the midpoint. So the four points at t are computed twice per
                // lane rather than four times, halving the FMA count of this
                // loop. The corner expressions themselves are untouched, so the
                // results stay bit-identical -- see tight_corner_at.
                for (int g = 0; g < 2; ++g) {
                    const int k0 = g;
                    const int k1 = g + 2;
                    const T* const SCCD_RESTRICT tk = mt[k0];
                    const T* const SCCD_RESTRICT u0 = mu[k0];
                    const T* const SCCD_RESTRICT v0 = mv[k0];
                    const T* const SCCD_RESTRICT u1 = mu[k1];
                    const T* const SCCD_RESTRICT v1 = mv[k1];
                    T* const SCCD_RESTRICT o0 = out[d][k0];
                    T* const SCCD_RESTRICT o1 = out[d][k1];
                    // Only about half the lanes split on any given pass -- the rest
                    // were accepted or rejected and their mid face is never read.
                    // Measured at 51%, so skipping them removes roughly a quarter
                    // of this kernel's total work.
                    //
                    // The test is per lane rather than a compaction because the
                    // geometry is indexed by lane: compacting would trade these
                    // branches for a gather of the same width. On a target where
                    // the loop below does vectorize, prefer predication here.
                    for (int l = 0; l < VSIZE; ++l) {
                        if (!wanted[l]) {
                            continue;
                        }
                        const T t = tk[l];
                        assert(t == mt[k1][l] && "mid-face corners k and k+2 must share a t value");
                        const T p0 = d0[l] * t + a0[l];
                        const T p1 = d1[l] * t + a1[l];
                        const T p2 = d2[l] * t + a2[l];
                        const T p3 = d3[l] * t + a3[l];
                        o0[l] = tight_corner_at<IsVertexFace, T>(p0, p1, p2, p3, u0[l], v0[l]);
                        o1[l] = tight_corner_at<IsVertexFace, T>(p0, p1, p2, p3, u1[l], v1[l]);
                    }
                }
            }
        }

        /** \brief All 8 corners of one box, scalar; used only to seed a root. */
        template <bool IsVertexFace, typename T>
        static void tight_eval_all_corners(const T start[4][3],
                                        const T delta[4][3],
                                        const T lo[3],
                                        const T hi[3],
                                        T corner[3][8]) {
            for (int d = 0; d < 3; ++d) {
                for (int c = 0; c < 8; ++c) {
                    const T t = (c & 4) ? hi[0] : lo[0];
                    const T u = (c & 2) ? hi[1] : lo[1];
                    const T v = (c & 1) ? hi[2] : lo[2];
                    corner[d][c] = tight_corner<IsVertexFace, T>(
                        start[0][d], start[1][d], start[2][d], start[3][d],
                        delta[0][d], delta[1][d], delta[2][d], delta[3][d], t, u, v);
                }
            }
        }

        /**
         * \brief TightInclusion's acceptance test for one lane.
         * \return False when the box provably contains no root and can be dropped.
         */
        template <typename T, int VSIZE>
        static SCCD_ALWAYS_INLINE bool tight_classify(const T (&fmin)[3][VSIZE],
                                                   const T (&fmax)[3][VSIZE],
                                                   const T (&box_lo)[3][VSIZE],
                                                   const T (&box_hi)[3][VSIZE],
                                                   const TightLaneGeometry<T, VSIZE>& lane,
                                                   const int l,
                                                   bool& accept) {
            bool box_in = true;
            bool within_tol = true;

            for (int d = 0; d < 3; ++d) {
                const T eps = lane.err[d][l];
                if (fmin[d][l] > eps || fmax[d][l] < -eps) {
                    accept = false;
                    return false;  // origin is outside the function's bounding box
                }
                box_in = box_in && (fmin[d][l] >= -eps) && (fmax[d][l] <= eps);
                // Note: the width tested here is the *domain* width, against the
                // domain tolerance. Comparing codomain widths against these is
                // what made the fast kernel report late times of impact.
                within_tol = within_tol && ((box_hi[d][l] - box_lo[d][l]) <= lane.tol[d][l]);
            }

            accept = box_in || within_tol;
            return true;
        }

        /** \brief Axis with the largest width/tolerance ratio (TightInclusion's rule). */
        template <typename T, int VSIZE>
        static SCCD_ALWAYS_INLINE int tight_split_axis(const T (&box_lo)[3][VSIZE],
                                                    const T (&box_hi)[3][VSIZE],
                                                    const TightLaneGeometry<T, VSIZE>& lane,
                                                    const int l) {
            int best = 0;
            T best_ratio = -std::numeric_limits<T>::max();
            for (int d = 0; d < 3; ++d) {
                const T width = box_hi[d][l] - box_lo[d][l];
                const T tol = lane.tol[d][l];
                if (width > tol) {
                    const T ratio = width / tol;
                    if (ratio > best_ratio) {
                        best_ratio = ratio;
                        best = d;
                    }
                }
            }
            return best;
        }

    }  // namespace detail

    /**
     * \brief TightInclusion-equivalent vectorized narrow phase.
     *
     * \tparam IsVertexFace true for vertex-face queries, false for edge-edge.
     * \param overlap0 vertex ids (VF) or first-edge ids (EE).
     * \param overlap1 face ids (VF) or second-edge ids (EE).
     * \param prims    faces[3] (VF) or edges[2] (EE).
     * \param toi_output 0 writes one global minimum, 1 writes one value per query.
     */
    template <bool IsVertexFace, typename T, typename I>
    static int narrow_phase_tight_impl(const size_t noverlaps,
                                      const I* const SCCD_RESTRICT overlap0,
                                      const I* const SCCD_RESTRICT overlap1,
                                      T** const SCCD_RESTRICT v0,
                                      T** const SCCD_RESTRICT v1,
                                      const size_t element_stride,
                                      I** const SCCD_RESTRICT prims,
                                      const T max_toi,
                                      T* const SCCD_RESTRICT toi,
                                      const int max_depth,
                                      const T tol,
                                      const ToiOutput toi_output) {
        using T_HP = double;
        constexpr int VSIZE =
            IsVertexFace ? SCCD_NARROWPHASE_TIGHT_VSIZE_VF : SCCD_NARROWPHASE_TIGHT_VSIZE_EE;
        if (noverlaps == 0) {
            if (toi != nullptr && toi_output == ToiOutput::Earliest) {
                toi[0] = sccd::min<T>(T(1), max_toi);
            }
            return 0;
        }
        assert(toi != nullptr);

        const T_HP domain_toi = sccd::min<T_HP>(T_HP(1), T_HP(max_toi));
        const ptrdiff_t nblocks = static_cast<ptrdiff_t>((noverlaps + VSIZE - 1) / VSIZE);
        std::atomic<T_HP> global_min{domain_toi};

#ifdef SCCD_NP_COUNT_BOXES
        // Each block owns VSIZE consecutive queries and no two blocks share one,
        // so the counters are disjoint per thread and need no synchronisation.
        std::vector<unsigned long long> perq_(noverlaps, 0ull);
#endif
        sccd::parallel_for_br_dynamic(0, nblocks, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            // Deliberately a plain local, not thread_local. The per-chunk
            // reserve looks like allocation churn, but measured on NEON over
            // armadillo-rollers it is worth nothing at all (-1.5% VF, -0.2% EE,
            // i.e. inside the noise): the scheduler's chunks are large enough to
            // amortize it and the allocator reuses the block. thread_local would
            // trade that for a permanent ~250 KB per thread and an indirection
            // on every access.
            std::vector<TightBox<int>> stack;
            stack.reserve(1024);

            for (ptrdiff_t ib = rbegin; ib < rend; ++ib) {
                const ptrdiff_t block_begin = ib * VSIZE;
                (void)block_begin;
                const int block_size =
                    static_cast<int>(sccd::min<ptrdiff_t>(VSIZE, static_cast<ptrdiff_t>(noverlaps) - block_begin));

                TightQueryBlock<T_HP, VSIZE> block;
                T_HP toi_q[VSIZE];

#ifdef SCCD_NP_COUNT_BOXES
                // Measurement knob, present only in an instrumented build.
                //
                // For toi_output == ToiOutput::Earliest the host seeds each block of queries with
                // the global minimum as it stands when the block starts, so a
                // block scheduled late searches a t-window that earlier blocks
                // have already collapsed. The device has no such sequence: one
                // launch starts every query at once against max_toi and the bound
                // only tightens as the launch runs.
                //
                // SCCD_NP_NO_GLOBAL_SEED=1 removes the host's advantage, seeding
                // every block at max_toi. If the host's box count then approaches
                // the device's, the 94x gap is the bound schedule and not the
                // kernel -- which is the question this knob exists to answer.
                static const bool no_global_seed_ = [] {
                    const char* e = getenv("SCCD_NP_NO_GLOBAL_SEED");
                    return e != nullptr && atoi(e) != 0;
                }();
                const T_HP seed = (toi_output == ToiOutput::Earliest && !no_global_seed_)
                                      ? global_min.load(std::memory_order_relaxed)
                                      : domain_toi;
#else
                const T_HP seed = toi_output == ToiOutput::Earliest ? global_min.load(std::memory_order_relaxed) : domain_toi;
#endif

                for (int q = 0; q < VSIZE; ++q) {
                    toi_q[q] = seed;
                    for (int k = 0; k < 4; ++k) {
                        for (int d = 0; d < 3; ++d) {
                            block.start[k][d][q] = T_HP(0);
                            block.delta[k][d][q] = T_HP(0);
                        }
                    }
                    for (int d = 0; d < 3; ++d) {
                        block.tol[d][q] = T_HP(1);
                        block.err[d][q] = T_HP(0);
                    }
                }

                for (int q = 0; q < block_size; ++q) {
                    const ptrdiff_t i = block_begin + q;
                    I node[4];
                    if constexpr (IsVertexFace) {
                        const size_t f = static_cast<size_t>(overlap1[i]) * element_stride;
                        node[0] = overlap0[i];
                        node[1] = prims[0][f];
                        node[2] = prims[1][f];
                        node[3] = prims[2][f];
                    } else {
                        const size_t ea = static_cast<size_t>(overlap0[i]) * element_stride;
                        const size_t eb = static_cast<size_t>(overlap1[i]) * element_stride;
                        node[0] = prims[0][ea];
                        node[1] = prims[1][ea];
                        node[2] = prims[0][eb];
                        node[3] = prims[1][eb];
                    }

                    T_HP ps[4][3];
                    T_HP pe[4][3];
                    for (int k = 0; k < 4; ++k) {
                        for (int d = 0; d < 3; ++d) {
                            ps[k][d] = T_HP(v0[d][node[k]]);
                            pe[k][d] = T_HP(v1[d][node[k]]);
                            block.start[k][d][q] = ps[k][d];
                            block.delta[k][d][q] = pe[k][d] - ps[k][d];
                        }
                    }

                    T_HP qt[3];
                    if constexpr (IsVertexFace) {
                        compute_face_vertex_tolerance<T_HP>(
                            T_HP(tol), ps[0], ps[1], ps[2], ps[3], pe[0], pe[1], pe[2], pe[3], qt);
                    } else {
                        compute_edge_edge_tolerance<T_HP>(
                            T_HP(tol), ps[0], ps[1], ps[2], ps[3], pe[0], pe[1], pe[2], pe[3], qt);
                    }
                    T_HP qe[3];
                    numerical_error_bound<IsVertexFace, T_HP>(
                        ps[0], ps[1], ps[2], ps[3], pe[0], pe[1], pe[2], pe[3], qe);

                    for (int d = 0; d < 3; ++d) {
                        block.tol[d][q] = qt[d];
                        block.err[d][q] = qe[d];
                    }
                }

                stack.clear();
                for (int q = block_size - 1; q >= 0; --q) {
                    TightBox<int> root;
                    root.query = q;
                    root.depth = 0;
                    root.lo[0] = T_HP(0);
                    root.hi[0] = toi_q[q];
                    root.lo[1] = T_HP(0);
                    root.hi[1] = T_HP(1);
                    root.lo[2] = T_HP(0);
                    root.hi[2] = T_HP(1);
                    T_HP rs[4][3];
                    T_HP rd[4][3];
                    for (int k = 0; k < 4; ++k) {
                        for (int d = 0; d < 3; ++d) {
                            rs[k][d] = block.start[k][d][q];
                            rd[k][d] = block.delta[k][d][q];
                        }
                    }
                    detail::tight_eval_all_corners<IsVertexFace, T_HP>(rs, rd, root.lo, root.hi, root.corner);

                    // Pushed even when the time budget is already zero: a box of
                    // zero extent in t can still hold a contact at t == 0, and
                    // dropping it would be a missed collision.
                    stack.push_back(root);
                }

                alignas(64) T_HP box_lo[3][VSIZE];
                alignas(64) T_HP box_hi[3][VSIZE];
                alignas(64) T_HP corner[3][8][VSIZE];
                alignas(64) T_HP fmin[3][VSIZE];
                alignas(64) T_HP fmax[3][VSIZE];
                alignas(64) T_HP mt[4][VSIZE];
                alignas(64) T_HP mu[4][VSIZE];
                alignas(64) T_HP mv[4][VSIZE];
                alignas(64) T_HP mid_face[3][4][VSIZE];
                int split_axis[VSIZE];
                T_HP split_mid[VSIZE];
                uint8_t will_split[VSIZE];
                detail::TightLaneGeometry<T_HP, VSIZE> lane;
                int lane_query[VSIZE];
                int lane_depth[VSIZE];
                uint8_t active[VSIZE];

                for (int l = 0; l < VSIZE; ++l) {
                    active[l] = 0;
                    lane_query[l] = 0;
                    lane_depth[l] = 0;
                    for (int d = 0; d < 3; ++d) {
                        box_lo[d][l] = T_HP(0);
                        box_hi[d][l] = T_HP(0);
                        for (int c = 0; c < 8; ++c) {
                            corner[d][c][l] = T_HP(0);
                        }
                    }
                    split_axis[l] = 0;
                    split_mid[l] = T_HP(0);
                    will_split[l] = 0;
                    // Idle lanes evaluate a zero box over query 0's geometry;
                    // the result is finite and discarded.
                    detail::tight_load_lane<T_HP, VSIZE>(block, 0, l, lane);
                }

                while (true) {
                    // Refill idle lanes from the shared work list. Boxes that can
                    // no longer improve their query's minimum are dropped here.
                    int filled = 0;
                    for (int l = 0; l < VSIZE; ++l) {
                        if (active[l]) {
                            ++filled;
                            continue;
                        }
                        while (!stack.empty()) {
                            const TightBox<int> b = stack.back();
                            stack.pop_back();
                            if (b.lo[0] >= toi_q[b.query]) {
                                continue;
                            }
                            if (lane_query[l] != b.query) {
                                detail::tight_load_lane<T_HP, VSIZE>(block, b.query, l, lane);
                            }
                            lane_query[l] = b.query;
                            lane_depth[l] = b.depth;
                            for (int d = 0; d < 3; ++d) {
                                box_lo[d][l] = b.lo[d];
                                box_hi[d][l] = b.hi[d];
                                for (int c = 0; c < 8; ++c) {
                                    corner[d][c][l] = b.corner[d][c];
                                }
                            }
                            active[l] = 1;
                            ++filled;
                            break;
                        }
                        // Deliberately no early exit when the work list is
                        // empty: lanes now persist across iterations (a split
                        // keeps its lower half here), so a later lane can still
                        // be live even when this one cannot be refilled.
                    }

                    if (filled == 0) {
                        break;
                    }

                    detail::tight_bounds_from_corners<T_HP, VSIZE>(corner, fmin, fmax);

                    // Phase 1: decide each lane's fate from the bounds it already
                    // carries. No function evaluation happens here.
                    for (int l = 0; l < VSIZE; ++l) {
                        will_split[l] = 0;
                        // Idle lanes take part in phase 2 with a degenerate face.
                        split_axis[l] = 0;
                        split_mid[l] = box_lo[0][l];

                        if (!active[l]) {
                            continue;
                        }
                        const int q = lane_query[l];
                        active[l] = 0;  // this box is consumed either way

                        bool accept = false;
                        SCCD_NP_HOST_BOX_TICK();
#ifdef SCCD_NP_COUNT_BOXES
                        SCCD_NP_HOST_PERQ_TICK(block_begin + lane_query[l]);
                        {
                            const int lv = lane_depth[l];
                            ++g_np_host_level[lv < 80 ? lv : 79];
                        }
#endif
                        if (!detail::tight_classify<T_HP, VSIZE>(
                                fmin, fmax, box_lo, box_hi, lane, l, accept)) {
                            continue;  // no root in this box
                        }

#ifdef SCCD_NP_COUNT_BOXES
                        if (!accept && lane_depth[l] >= max_depth) ++g_np_host_depth_accept;
#endif
                        if (accept || lane_depth[l] >= max_depth) {
                            if (box_lo[0][l] < toi_q[q]) {
                                toi_q[q] = box_lo[0][l];
                            }
                            continue;
                        }

                        const int axis = detail::tight_split_axis<T_HP, VSIZE>(box_lo, box_hi, lane, l);
                        T_HP mid;
                        if (!tight_bisect<T_HP>(box_lo[axis][l], box_hi[axis][l], mid)) {
                            // Cannot subdivide further; accepting is conservative.
                            if (box_lo[0][l] < toi_q[q]) {
                                toi_q[q] = box_lo[0][l];
                            }
                            continue;
                        }

                        split_axis[l] = axis;
                        split_mid[l] = mid;
                        will_split[l] = 1;
                    }

                    // Phase 2: lay out the 4 mid-face corners of every lane and
                    // evaluate them in one vectorized sweep. Lanes that are not
                    // splitting ride along; discarding their result is cheaper
                    // than branching inside the SIMD loop.
                    for (int l = 0; l < VSIZE; ++l) {
                        if (!will_split[l]) {
                            continue;
                        }
                        const int axis = split_axis[l];
                        const T_HP mid = split_mid[l];
                        for (int k = 0; k < 4; ++k) {
                            const int c = tight_mid_corner[axis][k];
                            mt[k][l] = (axis == 0) ? mid : ((c & 4) ? box_hi[0][l] : box_lo[0][l]);
                            mu[k][l] = (axis == 1) ? mid : ((c & 2) ? box_hi[1][l] : box_lo[1][l]);
                            mv[k][l] = (axis == 2) ? mid : ((c & 1) ? box_hi[2][l] : box_lo[2][l]);
                        }
                    }

                    detail::tight_eval_mid_face<IsVertexFace, T_HP, VSIZE>(
                        lane, mt, mu, mv, will_split, mid_face);

                    // Phase 3: assemble both children from the parent's corners
                    // plus the shared mid face, and push them.
                    for (int l = 0; l < VSIZE; ++l) {
                        if (!will_split[l]) {
                            continue;
                        }
                        const int q = lane_query[l];
                        const int axis = split_axis[l];
                        const int bit = tight_corner_bit(axis);
                        const T_HP mid = split_mid[l];
                        const T_HP other_lo = (axis == 1) ? box_lo[2][l] : box_lo[1][l];

                        // Upper half [mid, hi]: outer face inherited from the
                        // parent, inner face is the newly evaluated mid face.
                        // The prune is tested before the slot is filled so a
                        // rejected child never costs a 200-byte write.
                        const bool keep_upper =
                            axis == 0 ? (mid < toi_q[q])
                                      : (!IsVertexFace || (mid + other_lo <= T_HP(1)));
                        if (keep_upper) {
                            stack.emplace_back();
                            TightBox<int>& child = stack.back();
                            child.query = q;
                            child.depth = lane_depth[l] + 1;
                            for (int d = 0; d < 3; ++d) {
                                child.lo[d] = box_lo[d][l];
                                child.hi[d] = box_hi[d][l];
                            }
                            child.lo[axis] = mid;
                            child.hi[axis] = box_hi[axis][l];
                            for (int d = 0; d < 3; ++d) {
                                for (int k = 0; k < 4; ++k) {
                                    const int c = tight_mid_corner[axis][k];
                                    child.corner[d][c] = corner[d][c][l];
                                    child.corner[d][c & ~bit] = mid_face[d][k][l];
                                }
                            }
                        }

                        // Lower half [lo, mid] stays in this lane rather than
                        // going through the work list. It is the half we want to
                        // descend into first anyway, so keeping it here makes the
                        // traversal depth-first per lane -- the minimum converges
                        // sooner and prunes harder -- while halving the stack
                        // traffic these fat boxes would otherwise cost and
                        // skipping a refill (the lane's query is unchanged, so
                        // its geometry does not need reloading).
                        //
                        // Safe to overwrite the lane in place only because the
                        // upper child above has already read the parent's corners.
                        const bool keep_lower =
                            axis == 0 ? (box_lo[0][l] < toi_q[q])
                                      : (!IsVertexFace || (box_lo[axis][l] + other_lo <= T_HP(1)));
                        if (keep_lower) {
                            for (int d = 0; d < 3; ++d) {
                                for (int k = 0; k < 4; ++k) {
                                    corner[d][tight_mid_corner[axis][k]][l] = mid_face[d][k][l];
                                }
                            }
                            box_hi[axis][l] = mid;
                            lane_depth[l] = lane_depth[l] + 1;
                            active[l] = 1;
                        }
                    }
                }

                for (int q = 0; q < block_size; ++q) {
                    if (toi_output == ToiOutput::PerPair) {
                        toi[block_begin + q] = static_cast<T>(toi_q[q]);
                    } else {
                        T_HP current = global_min.load(std::memory_order_relaxed);
                        while (toi_q[q] < current &&
                               !global_min.compare_exchange_weak(
                                   current, toi_q[q], std::memory_order_relaxed, std::memory_order_relaxed)) {
                        }
                    }
                }
            }
        });

        if (toi_output == ToiOutput::Earliest) {
            toi[0] = static_cast<T>(global_min.load(std::memory_order_relaxed));
        }
#ifdef SCCD_NP_COUNT_BOXES
        {
            unsigned long long hist[24] = {0};
            unsigned long long worst = 0;
            size_t worst_q = 0;
            for (size_t q = 0; q < noverlaps; ++q) {
                np_hist_add(hist, perq_[q]);
                if (perq_[q] > worst) {
                    worst = perq_[q];
                    worst_q = q;
                }
            }
            fprintf(stderr, "sccd-np-hist host stride=%d queries=%zu worst=%llu at=%zu hist=",
                    toi_output, noverlaps, worst, worst_q);
            for (int b = 0; b < 24; ++b) fprintf(stderr, "%llu%s", hist[b], b == 23 ? "\n" : ",");
            fprintf(stderr, "sccd-np-level host depth_accept=%llu levels=", g_np_host_depth_accept);
            for (int b = 0; b < 80; ++b) fprintf(stderr, "%llu%s", g_np_host_level[b], b == 79 ? "\n" : ",");
            for (int b = 0; b < 80; ++b) g_np_host_level[b] = 0;
            g_np_host_depth_accept = 0;
        }
#endif
        return 0;
    }

    template <typename T, typename I>
    int narrow_phase_tight_vf(const size_t noverlaps,
                             const I* const SCCD_RESTRICT voverlap,
                             const I* const SCCD_RESTRICT first_out,
                             T** const SCCD_RESTRICT v0,
                             T** const SCCD_RESTRICT v1,
                             const size_t element_stride,
                             I** const SCCD_RESTRICT faces,
                             const T max_toi,
                             T* const SCCD_RESTRICT toi,
                             const int max_depth,
                             const T tol,
                             const ToiOutput toi_output) {
#ifdef SCCD_NP_COUNT_BOXES
        const unsigned long long before_ = g_np_host_boxes;
        const int rc_ = narrow_phase_tight_impl<true, T, I>(
            noverlaps, voverlap, first_out, v0, v1, element_stride, faces, max_toi, toi, max_depth, tol, toi_output);
        fprintf(stderr, "sccd-np-count vf host-conservative stride=%d queries=%zu boxes=%llu per_query=%.1f\n",
                toi_output, noverlaps, g_np_host_boxes - before_,
                noverlaps ? (double)(g_np_host_boxes - before_) / (double)noverlaps : 0.0);
        return rc_;
#else
        return narrow_phase_tight_impl<true, T, I>(
            noverlaps, voverlap, first_out, v0, v1, element_stride, faces, max_toi, toi, max_depth, tol, toi_output);
#endif
    }

    template <typename T, typename I>
    int narrow_phase_tight_ee(const size_t noverlaps,
                             const I* const SCCD_RESTRICT e0overlap,
                             const I* const SCCD_RESTRICT e1overlap,
                             T** const SCCD_RESTRICT v0,
                             T** const SCCD_RESTRICT v1,
                             const size_t element_stride,
                             I** const SCCD_RESTRICT edges,
                             const T max_toi,
                             T* const SCCD_RESTRICT toi,
                             const int max_depth,
                             const T tol,
                             const ToiOutput toi_output) {
#ifdef SCCD_NP_COUNT_BOXES
        const unsigned long long before_ = g_np_host_boxes;
        const int rc_ = narrow_phase_tight_impl<false, T, I>(
            noverlaps, e0overlap, e1overlap, v0, v1, element_stride, edges, max_toi, toi, max_depth, tol, toi_output);
        fprintf(stderr, "sccd-np-count ee host-conservative stride=%d queries=%zu boxes=%llu per_query=%.1f\n",
                toi_output, noverlaps, g_np_host_boxes - before_,
                noverlaps ? (double)(g_np_host_boxes - before_) / (double)noverlaps : 0.0);
        return rc_;
#else
        return narrow_phase_tight_impl<false, T, I>(
            noverlaps, e0overlap, e1overlap, v0, v1, element_stride, edges, max_toi, toi, max_depth, tol, toi_output);
#endif
    }

    SCCD_FP_STRICT_END

}  // namespace sccd

#endif  // SCCD_NARROWPHASE_TIGHT_HPP
