#ifndef SCCD_ROOTFINDER_QUAD_HPP
#define SCCD_ROOTFINDER_QUAD_HPP

#include "sccd_rootfinder.hpp"

namespace sccd {

    template <typename T>
    inline void diff_vq(const T sv[3],
                        const T s1[3],
                        const T s2[3],
                        const T s3[3],
                        const T s4[3],
                        const T ev[3],
                        const T e1[3],
                        const T e2[3],
                        const T e3[3],
                        const T e4[3],
                        const T t,
                        const T u,
                        const T v,
                        T *const SCCD_RESTRICT diff) {
        const T omt = T(1) - t;
        const T omu = T(1) - u;
        const T omv = T(1) - v;
        const T w1 = omu * omv;
        const T w2 = u * omv;
        const T w3 = omu * v;
        const T w4 = u * v;

        for (int d = 0; d < 3; ++d) {
            const T vertex = omt * sv[d] + t * ev[d];
            const T quad_s = w1 * s1[d] + w2 * s2[d] + w3 * s3[d] + w4 * s4[d];
            const T quad_e = w1 * e1[d] + w2 * e2[d] + w3 * e3[d] + w4 * e4[d];
            diff[d] = vertex - (omt * quad_s + t * quad_e);
        }
    }

    /**
     * \brief The vertex-quad inclusion function, factored over the quad's own
     *        parameter domain.
     *
     * `diff_vq` above is the reference form and stays that way: it is what the
     * equivalence test checks against. This is the same function regrouped so
     * that the work which does not depend on (u, v) can be hoisted out of the
     * corner loop.
     *
     * Write the quad's bilinear blend with weights w1..w4 summing to one. The
     * reference evaluates
     *
     *     F = V(t) - [ (1-t) * SUM_k w_k s_k  +  t * SUM_k w_k e_k ]
     *
     * blending in (u, v) first and in t second. Since the weights do not depend
     * on t, that equals
     *
     *     F = V(t) - SUM_k w_k * [ (1-t) s_k + t e_k ]  =  V(t) - SUM_k w_k P_k(t)
     *
     * which blends in t first. The point of the second form is that P_k(t) and
     * V(t) depend on nothing but t: a box has two t bounds, and a split along u
     * or v holds both of them fixed across every corner it evaluates. So the
     * frame is computed twice for the whole split instead of once per corner.
     *
     * The two forms are equal in exact arithmetic but not bit-for-bit in
     * floating point, since the multiplications associate differently. Both are
     * trilinear evaluations of the same depth -- each term is a product of three
     * rounded factors summed with the others -- so the certified bound that makes
     * the rejection sound covers this form as it covers the reference. That is an
     * argument, not a proof, which is why `vertex_quad_root_test` checks the two
     * against each other and against a densely sampled reference root rather than
     * taking it on trust.
     */
    template <typename T>
    struct VQFrame {
        T V[3];
        T P[4][3];
    };

    /** \brief Evaluate the t-dependent part of the inclusion function. */
    template <typename T>
    inline void vq_frame_at(const T sv[3],
                            const T s1[3],
                            const T s2[3],
                            const T s3[3],
                            const T s4[3],
                            const T ev[3],
                            const T e1[3],
                            const T e2[3],
                            const T e3[3],
                            const T e4[3],
                            const T t,
                            VQFrame<T> &frame) {
        const T omt = T(1) - t;
        const T *const sk[4] = {s1, s2, s3, s4};
        const T *const ek[4] = {e1, e2, e3, e4};
        for (int d = 0; d < 3; ++d) {
            frame.V[d] = omt * sv[d] + t * ev[d];
        }
        for (int k = 0; k < 4; ++k) {
            for (int d = 0; d < 3; ++d) {
                frame.P[k][d] = omt * sk[k][d] + t * ek[k][d];
            }
        }
    }

    /** \brief Blend a frame at one (u, v) corner. */
    template <typename T>
    inline void vq_eval_frame(const VQFrame<T> &frame, const T u, const T v, T *const SCCD_RESTRICT F) {
        const T omu = T(1) - u;
        const T omv = T(1) - v;
        const T w1 = omu * omv;
        const T w2 = u * omv;
        const T w3 = omu * v;
        const T w4 = u * v;
        for (int d = 0; d < 3; ++d) {
            F[d] = frame.V[d] -
                   (w1 * frame.P[0][d] + w2 * frame.P[1][d] + w3 * frame.P[2][d] + w4 * frame.P[3][d]);
        }
    }

    /**
     * \brief The per-query reductions the search needs, in one pass.
     *
     * Three things are derived from a vertex-quad query before the search
     * starts, and they were computed by three separate walks over the same
     * thirty coordinates:
     *
     *  - the per-axis Lipschitz constants, which are the codomain widths;
     *  - the per-axis tolerances, which are a constant divided by those same
     *    Lipschitz constants;
     *  - the certified numerical error bound, from the largest coordinate.
     *
     * The first two were computing the identical reduction twice --
     * `compute_vertex_quad_tolerance` and `compute_vertex_quad_codomain_widths`
     * have the same loop body, one dividing at the end and the other not.
     *
     * That matters more than it sounds. In the regime the mesh pipeline actually
     * runs -- a time of impact shared across queries, so the bound tightens and
     * prunes -- the search explores about 1.2 boxes per query. The preamble is
     * not amortised over a deep search; it is most of the query.
     */
    template <typename T>
    struct VQBounds {
        T lipschitz[3];
        T max_coord[3];
    };

    template <typename T>
    inline void vq_bounds(const T sv[3],
                          const T s1[3],
                          const T s2[3],
                          const T s3[3],
                          const T s4[3],
                          const T ev[3],
                          const T e1[3],
                          const T e2[3],
                          const T e3[3],
                          const T e4[3],
                          VQBounds<T> &out) {
        for (int i = 0; i < 3; ++i) out.lipschitz[i] = T(0);

        for (int d = 0; d < 3; ++d) {
            const T vt = ev[d] - sv[d];
            out.lipschitz[0] = sccd::max<T>(
                out.lipschitz[0],
                sccd::max<T>(sccd::max<T>(sccd::abs<T>(vt - (e1[d] - s1[d])), sccd::abs<T>(vt - (e2[d] - s2[d]))),
                             sccd::max<T>(sccd::abs<T>(vt - (e3[d] - s3[d])), sccd::abs<T>(vt - (e4[d] - s4[d])))));
            out.lipschitz[1] = sccd::max<T>(
                out.lipschitz[1],
                sccd::max<T>(sccd::max<T>(sccd::abs<T>(s2[d] - s1[d]), sccd::abs<T>(s4[d] - s3[d])),
                             sccd::max<T>(sccd::abs<T>(e2[d] - e1[d]), sccd::abs<T>(e4[d] - e3[d]))));
            out.lipschitz[2] = sccd::max<T>(
                out.lipschitz[2],
                sccd::max<T>(sccd::max<T>(sccd::abs<T>(s3[d] - s1[d]), sccd::abs<T>(s4[d] - s2[d])),
                             sccd::max<T>(sccd::abs<T>(e3[d] - e1[d]), sccd::abs<T>(e4[d] - e2[d]))));

            // Same walk, so the error bound's reduction rides along with it.
            const T m = sccd::max<T>(
                sccd::max<T>(sccd::max<T>(sccd::abs<T>(sv[d]), sccd::abs<T>(s1[d])),
                             sccd::max<T>(sccd::abs<T>(s2[d]), sccd::abs<T>(s3[d]))),
                sccd::max<T>(sccd::max<T>(sccd::abs<T>(s4[d]), sccd::abs<T>(ev[d])),
                             sccd::max<T>(sccd::max<T>(sccd::abs<T>(e1[d]), sccd::abs<T>(e2[d])),
                                          sccd::max<T>(sccd::abs<T>(e3[d]), sccd::abs<T>(e4[d])))));
            out.max_coord[d] = sccd::max<T>(T(1), m);
        }
    }

    template <typename T>
    inline void compute_vertex_quad_tolerance(const T codomain_tol,
                                              const T sv[3],
                                              const T s1[3],
                                              const T s2[3],
                                              const T s3[3],
                                              const T s4[3],
                                              const T ev[3],
                                              const T e1[3],
                                              const T e2[3],
                                              const T e3[3],
                                              const T e4[3],
                                              T *const SCCD_RESTRICT tol) {
        VQBounds<T> b;
        vq_bounds<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, b);
        const T *const lipschitz = b.lipschitz;

        const T axis_tol = codomain_tol / T(3);
        const T eps = std::numeric_limits<T>::epsilon();
        for (int d = 0; d < 3; ++d) {
            tol[d] = lipschitz[d] > eps ? axis_tol / lipschitz[d] : T(1);
        }
    }

    template <typename T>
    inline void compute_vertex_quad_codomain_widths(const T sv[3],
                                                    const T s1[3],
                                                    const T s2[3],
                                                    const T s3[3],
                                                    const T s4[3],
                                                    const T ev[3],
                                                    const T e1[3],
                                                    const T e2[3],
                                                    const T e3[3],
                                                    const T e4[3],
                                                    T widths[3]) {
        VQBounds<T> b;
        vq_bounds<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, b);
        widths[0] = b.lipschitz[0];
        widths[1] = b.lipschitz[1];
        widths[2] = b.lipschitz[2];
    }

    template <typename T>
    inline void normalize_vertex_quad_codomain_widths(T widths[3]) {
        const T total = widths[0] + widths[1] + widths[2] + std::numeric_limits<T>::epsilon();
        for (int d = 0; d < 3; ++d) {
            widths[d] = sccd::max<T>(T(1e-4), widths[d] / total);
        }
    }

    template <typename T>
    static inline void sccd_get_numerical_error_vq_soa(const int use_ms,
                                                       const T v0sx,
                                                       const T v0sy,
                                                       const T v0sz,
                                                       const T v1sx,
                                                       const T v1sy,
                                                       const T v1sz,
                                                       const T v2sx,
                                                       const T v2sy,
                                                       const T v2sz,
                                                       const T v3sx,
                                                       const T v3sy,
                                                       const T v3sz,
                                                       const T v4sx,
                                                       const T v4sy,
                                                       const T v4sz,
                                                       const T v0ex,
                                                       const T v0ey,
                                                       const T v0ez,
                                                       const T v1ex,
                                                       const T v1ey,
                                                       const T v1ez,
                                                       const T v2ex,
                                                       const T v2ey,
                                                       const T v2ez,
                                                       const T v3ex,
                                                       const T v3ey,
                                                       const T v3ez,
                                                       const T v4ex,
                                                       const T v4ey,
                                                       const T v4ez,
                                                       T *const SCCD_RESTRICT errx,
                                                       T *const SCCD_RESTRICT erry,
                                                       T *const SCCD_RESTRICT errz) {
        const T kFilter = (use_ms ? T(42) : T(38)) * std::numeric_limits<T>::epsilon();
        const T maxx =
            sccd::max<T>(T(1),
                         sccd::max<T>(
                             sccd::max<T>(sccd::max<T>(sccd::abs<T>(v0sx), sccd::abs<T>(v1sx)),
                                          sccd::max<T>(sccd::abs<T>(v2sx), sccd::abs<T>(v3sx))),
                             sccd::max<T>(sccd::max<T>(sccd::abs<T>(v4sx), sccd::abs<T>(v0ex)),
                                          sccd::max<T>(
                                              sccd::max<T>(sccd::abs<T>(v1ex), sccd::abs<T>(v2ex)),
                                              sccd::max<T>(sccd::abs<T>(v3ex), sccd::abs<T>(v4ex))))));
        const T maxy =
            sccd::max<T>(T(1),
                         sccd::max<T>(
                             sccd::max<T>(sccd::max<T>(sccd::abs<T>(v0sy), sccd::abs<T>(v1sy)),
                                          sccd::max<T>(sccd::abs<T>(v2sy), sccd::abs<T>(v3sy))),
                             sccd::max<T>(sccd::max<T>(sccd::abs<T>(v4sy), sccd::abs<T>(v0ey)),
                                          sccd::max<T>(
                                              sccd::max<T>(sccd::abs<T>(v1ey), sccd::abs<T>(v2ey)),
                                              sccd::max<T>(sccd::abs<T>(v3ey), sccd::abs<T>(v4ey))))));
        const T maxz =
            sccd::max<T>(T(1),
                         sccd::max<T>(
                             sccd::max<T>(sccd::max<T>(sccd::abs<T>(v0sz), sccd::abs<T>(v1sz)),
                                          sccd::max<T>(sccd::abs<T>(v2sz), sccd::abs<T>(v3sz))),
                             sccd::max<T>(sccd::max<T>(sccd::abs<T>(v4sz), sccd::abs<T>(v0ez)),
                                          sccd::max<T>(
                                              sccd::max<T>(sccd::abs<T>(v1ez), sccd::abs<T>(v2ez)),
                                              sccd::max<T>(sccd::abs<T>(v3ez), sccd::abs<T>(v4ez))))));
        *errx = sccd::pow3<T>(maxx) * kFilter;
        *erry = sccd::pow3<T>(maxy) * kFilter;
        *errz = sccd::pow3<T>(maxz) * kFilter;
    }

    /**
     * \brief The certified numerical error bound for one vertex-quad query.
     *
     * A thin array-taking wrapper over the 30-scalar form above, so that callers
     * can compute it once per query and hand the same three numbers to every
     * box. It depends only on the query's coordinates, so recomputing it per box
     * -- which is what the search used to do -- is pure waste.
     */
    template <typename T>
    inline void vq_numerical_error(const T sv[3],
                                   const T s1[3],
                                   const T s2[3],
                                   const T s3[3],
                                   const T s4[3],
                                   const T ev[3],
                                   const T e1[3],
                                   const T e2[3],
                                   const T e3[3],
                                   const T e4[3],
                                   T out[3]) {
        sccd_get_numerical_error_vq_soa<T>(/*use_ms=*/0,
                                           sv[0], sv[1], sv[2],
                                           s1[0], s1[1], s1[2],
                                           s2[0], s2[1], s2[2],
                                           s3[0], s3[1], s3[2],
                                           s4[0], s4[1], s4[2],
                                           ev[0], ev[1], ev[2],
                                           e1[0], e1[1], e1[2],
                                           e2[0], e2[1], e2[2],
                                           e3[0], e3[1], e3[2],
                                           e4[0], e4[1], e4[2],
                                           &out[0], &out[1], &out[2]);
    }

    /**
     * \brief Everything the search needs from a query, from one pass.
     *
     * The hot path calls this instead of the three functions above, which now
     * share its reduction and remain for callers that want one of the three on
     * its own.
     */
    template <typename T>
    inline void vq_prepare(const T codomain_tol,
                           const T sv[3],
                           const T s1[3],
                           const T s2[3],
                           const T s3[3],
                           const T s4[3],
                           const T ev[3],
                           const T e1[3],
                           const T e2[3],
                           const T e3[3],
                           const T e4[3],
                           T tols[3],
                           T widths[3],
                           T numerical_error[3]) {
        VQBounds<T> b;
        vq_bounds<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, b);

        const T axis_tol = codomain_tol / T(3);
        const T eps = std::numeric_limits<T>::epsilon();
        const T kFilter = T(38) * eps;
        for (int d = 0; d < 3; ++d) {
            tols[d] = b.lipschitz[d] > eps ? axis_tol / b.lipschitz[d] : T(1);
            widths[d] = b.lipschitz[d];
            numerical_error[d] = sccd::pow3<T>(b.max_coord[d]) * kFilter;
        }
        normalize_vertex_quad_codomain_widths<T>(widths);
    }

    template <typename T>
    inline bool accept_grid_root_vq(const Box<T> &box, T &toi, T &u, T &v) {
        const T t_approx = box.tuv[0].lower;
        if (t_approx < toi) {
            toi = t_approx;
            u = box.tuv[1].lower;
            v = box.tuv[2].lower;
            return true;
        }
        return false;
    }

    template <int SplitDim, int N, typename T>
    inline static void normal_equation_axis_splitters_vq(const Box<T> &domain,
                                                         const T sv[3],
                                                         const T s1[3],
                                                         const T s2[3],
                                                         const T s3[3],
                                                         const T s4[3],
                                                         const T ev[3],
                                                         const T e1[3],
                                                         const T e2[3],
                                                         const T e3[3],
                                                         const T e4[3],
                                                         T *const SCCD_RESTRICT splitters) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        static_assert(N > 0);

        const T lo = domain.tuv[SplitDim].lower;
        const T hi = domain.tuv[SplitDim].upper;
        const T h = (hi - lo) / T(N + 1);

        T damping = T(0.45);
        if constexpr (N == 1) {
            damping = T(0.6);
        }

        const T radius = h * damping;
        const T mid_t = (domain.tuv[0].lower + domain.tuv[0].upper) * T(0.5);
        const T mid_u = (domain.tuv[1].lower + domain.tuv[1].upper) * T(0.5);
        const T mid_v = (domain.tuv[2].lower + domain.tuv[2].upper) * T(0.5);
        const T eps = std::numeric_limits<T>::epsilon();

        T F_base[3];
        T J_axis[3];
        T H_axis = T(0);

        const T base_t = SplitDim == 0 ? T(0) : mid_t;
        const T base_u = SplitDim == 1 ? T(0) : mid_u;
        const T base_v = SplitDim == 2 ? T(0) : mid_v;
        diff_vq<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, base_t, base_u, base_v, F_base);

        const T omt = T(1) - mid_t;
        const T omu = T(1) - mid_u;
        const T omv = T(1) - mid_v;
        for (int d = 0; d < 3; ++d) {
            if constexpr (SplitDim == 0) {
                const T sq = omu * omv * s1[d] + mid_u * omv * s2[d] + omu * mid_v * s3[d] + mid_u * mid_v * s4[d];
                const T eq = omu * omv * e1[d] + mid_u * omv * e2[d] + omu * mid_v * e3[d] + mid_u * mid_v * e4[d];
                J_axis[d] = (ev[d] - sv[d]) - (eq - sq);
            } else if constexpr (SplitDim == 1) {
                const T du_s = omv * (s2[d] - s1[d]) + mid_v * (s4[d] - s3[d]);
                const T du_e = omv * (e2[d] - e1[d]) + mid_v * (e4[d] - e3[d]);
                J_axis[d] = -(omt * du_s + mid_t * du_e);
            } else {
                const T dv_s = omu * (s3[d] - s1[d]) + mid_u * (s4[d] - s2[d]);
                const T dv_e = omu * (e3[d] - e1[d]) + mid_u * (e4[d] - e2[d]);
                J_axis[d] = -(omt * dv_s + mid_t * dv_e);
            }
            H_axis += J_axis[d] * J_axis[d];
        }

        const T step_scale = H_axis > eps ? T(1) / H_axis : T(0.00001);

#pragma omp simd aligned(splitters)
        for (int i = 0; i < N; ++i) {
            const T x0 = lo + h * T(i + 1);
            const T xmin = sccd::max<T>(lo, x0 - radius);
            const T xmax = sccd::min<T>(hi, x0 + radius);

            T g = T(0);
            for (int d = 0; d < 3; ++d) {
                const T J = J_axis[d];
                g += (F_base[d] + x0 * J) * J;
            }

            const T step = g * step_scale;
            splitters[i] = sccd::min<T>(xmax, sccd::max<T>(xmin, x0 - step));
        }
    }

    template <int SplitDim, int N, typename T>
    inline bool grid_search_adaptive_split_vq_axis(const sccd::Box<T> &domain,
                                                   const int max_iter,
                                                   const T tol,
                                                   const T tols[3],
                                                   const T numerical_error[3],
                                                   const T sv[3],
                                                   const T s1[3],
                                                   const T s2[3],
                                                   const T s3[3],
                                                   const T s4[3],
                                                   const T ev[3],
                                                   const T e1[3],
                                                   const T e2[3],
                                                   const T e3[3],
                                                   const T e4[3],
                                                   T &toi,
                                                   T &u,
                                                   T &v,
                                                   std::vector<sccd::Box<T>> &stack,
                                                   const bool refine) {
        // The certified numerical error bound is what makes the rejection below
        // sound, and it depends only on the query's coordinates -- not on the
        // box. It used to be recomputed here, on every box popped from the
        // stack: 30 absolute values, 30 maxima and three cubes per split, to
        // arrive at the same three numbers every time. It is now computed once
        // per query and passed in.
        //
        // (It also used to not be called at all. The acceptance test padded with
        // machine epsilon instead, roughly 30x too small for unit-scale
        // geometry, which let it discard boxes that contained a root.)

        (void)refine;

        const T lo = domain.tuv[SplitDim].lower;
        const T hi = domain.tuv[SplitDim].upper;

        alignas(64) T splitters[N];
        normal_equation_axis_splitters_vq<SplitDim, N, T>(
            domain, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, splitters);

        alignas(64) T samples[N + 2];
        samples[0] = lo;
        samples[N + 1] = hi;
#pragma omp simd aligned(splitters, samples)
        for (int i = 0; i < N; ++i) {
            samples[i + 1] = splitters[i];
        }

        // Evaluate each sample plane at most once, and only when a sub-box that
        // reads it survives the cutoff.
        //
        // Two separate savings, and they pull against each other if done naively.
        //
        // The first is sharing. Splitting into N+1 sub-boxes along one axis, the
        // naive loop evaluates eight corners each, 8(N+1) in total; but
        // consecutive sub-boxes share a face, so only 4(N+2) distinct corner
        // values exist. Carrying the previous plane forward gets that for free.
        //
        // The second is the frame. A corner costs a blend rather than a full
        // interpolation once the t-dependent part is precomputed, and a split
        // along u or v holds t at the domain's two bounds for every corner in the
        // split, so two frames serve all of them.
        //
        // The trap is that both invite hoisting the evaluation out of the
        // sub-box loop, and that loses more than it gains. A sub-box whose t
        // lower bound has reached the running time of impact is discarded
        // unlooked at, and the running value gets sharper *during* this loop
        // whenever a sub-box accepts. Evaluating all the planes up front measured
        // 1.3x faster per isolated query and 1.35x SLOWER end to end, where a
        // shared global minimum prunes hard; computing the cutoff once on entry
        // instead of re-reading it each iteration still left it slower than the
        // naive loop. So the planes are produced lazily, one per iteration, and
        // the cutoff is re-read exactly where the naive loop read it.
        T plane_a[4][3];
        T plane_b[4][3];
        bool have_a = false;

        const T dom_lo[3] = {domain.tuv[0].lower, domain.tuv[1].lower, domain.tuv[2].lower};
        const T dom_hi[3] = {domain.tuv[0].upper, domain.tuv[1].upper, domain.tuv[2].upper};

        // The two axes held fixed while SplitDim sweeps.
        constexpr int A = (SplitDim == 0) ? 1 : 0;
        constexpr int B = (SplitDim == 2) ? 1 : 2;

        VQFrame<T> frame_lo;
        VQFrame<T> frame_hi;
        if constexpr (SplitDim != 0) {
            vq_frame_at<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, dom_lo[0], frame_lo);
            vq_frame_at<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, dom_hi[0], frame_hi);
        }

        const auto eval_plane = [&](const int j, T dst[4][3]) {
            if constexpr (SplitDim == 0) {
                VQFrame<T> frame;
                vq_frame_at<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, samples[j], frame);
                for (int c = 0; c < 4; ++c) {
                    vq_eval_frame<T>(frame,
                                     (c & 1) ? dom_hi[A] : dom_lo[A],
                                     (c & 2) ? dom_hi[B] : dom_lo[B],
                                     dst[c]);
                }
            } else {
                for (int c = 0; c < 4; ++c) {
                    // Bit 0 selects the t bound, bit 1 the remaining fixed axis.
                    const VQFrame<T> &frame = (c & 1) ? frame_hi : frame_lo;
                    const T other = (c & 2) ? dom_hi[B] : dom_lo[B];
                    const T uu = (SplitDim == 1) ? samples[j] : other;
                    const T vv = (SplitDim == 1) ? other : samples[j];
                    vq_eval_frame<T>(frame, uu, vv, dst[c]);
                }
            }
        };

        const auto stack_size = stack.size();
        bool found = false;
        for (int i = 0; i < N + 1; ++i) {
            const T sample_min = samples[i];
            const T sample_max = samples[i + 1];
            const T tt_min = SplitDim == 0 ? sample_min : domain.tuv[0].lower;

            if (tt_min >= toi) {
                have_a = false;
                continue;
            }

            if (!have_a) {
                eval_plane(i, plane_a);
                have_a = true;
            }
            eval_plane(i + 1, plane_b);

            T fmin[3];
            T fmax[3];
            init_codomain_bounds<T>(fmin, fmax);
            for (int c = 0; c < 4; ++c) {
                update_codomain_bounds<T>(plane_a[c], fmin, fmax);
                update_codomain_bounds<T>(plane_b[c], fmin, fmax);
            }

            // The upper face of this sub-box is the lower face of the next.
            for (int c = 0; c < 4; ++c) {
                for (int d = 0; d < 3; ++d) plane_a[c][d] = plane_b[c][d];
            }

            bool accepted = false;
            if (!codomain_acceptance_vq<T>(fmin, fmax, tol, tols, numerical_error, accepted)) {
                continue;
            }

            // Reject a box whose t lower bound is zero. This is TightInclusion's
            // no_zero_toi policy: a contact reported at exactly t == 0 stalls a
            // solver, which cannot advance the step at all.
            //
            // It does not cost a collision. Measured against the datasets' exact
            // roots, every mode reports gt_missed == 0 and gt_late == 0 with this
            // in place -- see benchmark/oracle/README.md. TightInclusion does
            // report a hit at t == 0 on some of these queries, but its answer is
            // a conservative lower bound and over-reports; the exact roots agree
            // with rejecting them.
            accepted = accepted && (tt_min > 0);

            Box<T> box = split_axis_box<SplitDim, T>(domain, sample_min, sample_max);
            if (accepted || box.depth > max_iter) {
                found |= accept_grid_root_vq<T>(box, toi, u, v);
                continue;
            }

            stack.push_back(box);
        }

        if constexpr (SplitDim == 0) {
            std::reverse(stack.begin() + stack_size, stack.end());
        }

        return found;
    }

    template <int N, typename T>
    inline bool grid_search_adaptive_split_vq(const sccd::Box<T> &domain,
                                              const int max_iter,
                                              const T tol,
                                              const T tols[3],
                                              const T numerical_error[3],
                                              const T codomain_widths[3],
                                              const T sv[3],
                                              const T s1[3],
                                              const T s2[3],
                                              const T s3[3],
                                              const T s4[3],
                                              const T ev[3],
                                              const T e1[3],
                                              const T e2[3],
                                              const T e3[3],
                                              const T e4[3],
                                              T &toi,
                                              T &u,
                                              T &v,
                                              std::vector<sccd::Box<T>> &stack,
                                              const bool refine) {
        const int split_dim = domain.widest_dimension(codomain_widths);
        if (split_dim == 0) {
            return grid_search_adaptive_split_vq_axis<0, N, T>(
                domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_adaptive_split_vq_axis<1, N, T>(
                domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        return grid_search_adaptive_split_vq_axis<2, N, T>(
            domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, toi, u, v, stack, refine);
    }

    template <int N, typename T>
    inline bool grid_search_adaptive_split_vq(const sccd::Box<T> &domain,
                                              const int max_iter,
                                              const T tol,
                                              const T tols[3],
                                              const T sv[3],
                                              const T s1[3],
                                              const T s2[3],
                                              const T s3[3],
                                              const T s4[3],
                                              const T ev[3],
                                              const T e1[3],
                                              const T e2[3],
                                              const T e3[3],
                                              const T e4[3],
                                              T &toi,
                                              T &u,
                                              T &v,
                                              std::vector<sccd::Box<T>> &stack,
                                              const bool refine) {
        const T codomain_widths[3] = {T(1), T(1), T(1)};
        T numerical_error[3];
        vq_numerical_error<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, numerical_error);
        return grid_search_adaptive_split_vq<N, T>(
            domain, max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vq(const int max_iter,
                                          const T tol,
                                          const T tols[3],
                                          const T numerical_error[3],
                                          const T codomain_widths[3],
                                          const T sv[3],
                                          const T s1[3],
                                          const T s2[3],
                                          const T s3[3],
                                          const T s4[3],
                                          const T ev[3],
                                          const T e1[3],
                                          const T e2[3],
                                          const T e3[3],
                                          const T e4[3],
                                          const Box<T> &initial_domain,
                                          T &t,
                                          T &u,
                                          T &v,
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        if (initial_domain.tuv[0].lower >= t) {
            return false;
        }

        return grid_search_adaptive_split_vq<SCCD_ADAPTIVE_NUM_SPLITS, T>(
            initial_domain, max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vq(const int max_iter,
                                          const T tol,
                                          const T tols[3],
                                          const T sv[3],
                                          const T s1[3],
                                          const T s2[3],
                                          const T s3[3],
                                          const T s4[3],
                                          const T ev[3],
                                          const T e1[3],
                                          const T e2[3],
                                          const T e3[3],
                                          const T e4[3],
                                          const Box<T> &initial_domain,
                                          T &t,
                                          T &u,
                                          T &v,
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        const T codomain_widths[3] = {T(1), T(1), T(1)};
        T numerical_error[3];
        vq_numerical_error<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, numerical_error);
        return find_root_grid_adaptive_split_vq<T>(
            max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vq(const int max_iter,
                                          const T tol,
                                          const T sv[3],
                                          const T s1[3],
                                          const T s2[3],
                                          const T s3[3],
                                          const T s4[3],
                                          const T ev[3],
                                          const T e1[3],
                                          const T e2[3],
                                          const T e3[3],
                                          const T e4[3],
                                          const Box<T> &initial_domain,
                                          T &t,
                                          T &u,
                                          T &v,
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        T tols[3];
        T codomain_widths[3];
        compute_vertex_quad_tolerance<T>(tol, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, tols);
        compute_vertex_quad_codomain_widths<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, codomain_widths);
        normalize_vertex_quad_codomain_widths<T>(codomain_widths);
        T numerical_error[3];
        vq_numerical_error<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, numerical_error);
        return find_root_grid_adaptive_split_vq<T>(
            max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vq(const int max_iter,
                                          const T tol,
                                          const T sv[3],
                                          const T s1[3],
                                          const T s2[3],
                                          const T s3[3],
                                          const T s4[3],
                                          const T ev[3],
                                          const T e1[3],
                                          const T e2[3],
                                          const T e3[3],
                                          const T e4[3],
                                          T &t,
                                          T &u,
                                          T &v,
                                          std::vector<Box<T>> &stack,
                                          const bool refine = false) {
        return find_root_grid_adaptive_split_vq<T>(
            max_iter, tol, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, unit_domain_box<T>(), t, u, v, stack, refine);
    }

} // namespace sccd

#endif
