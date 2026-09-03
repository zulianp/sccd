#ifndef S_ROOT_FINDER_VERTEX_QUAD_HPP
#define S_ROOT_FINDER_VERTEX_QUAD_HPP

#include "srootfinder.hpp"

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
        T lipschitz[3] = {T(0), T(0), T(0)};
        for (int d = 0; d < 3; ++d) {
            const T vt = ev[d] - sv[d];
            lipschitz[0] =
                sccd::max<T>(lipschitz[0],
                              sccd::max<T>(
                                  sccd::max<T>(sccd::abs<T>(vt - (e1[d] - s1[d])),
                                               sccd::abs<T>(vt - (e2[d] - s2[d]))),
                                  sccd::max<T>(sccd::abs<T>(vt - (e3[d] - s3[d])),
                                               sccd::abs<T>(vt - (e4[d] - s4[d])))));

            lipschitz[1] =
                sccd::max<T>(lipschitz[1],
                              sccd::max<T>(
                                  sccd::max<T>(sccd::abs<T>(s2[d] - s1[d]), sccd::abs<T>(s4[d] - s3[d])),
                                  sccd::max<T>(sccd::abs<T>(e2[d] - e1[d]), sccd::abs<T>(e4[d] - e3[d]))));

            lipschitz[2] =
                sccd::max<T>(lipschitz[2],
                              sccd::max<T>(
                                  sccd::max<T>(sccd::abs<T>(s3[d] - s1[d]), sccd::abs<T>(s4[d] - s2[d])),
                                  sccd::max<T>(sccd::abs<T>(e3[d] - e1[d]), sccd::abs<T>(e4[d] - e2[d]))));
        }

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
        T wt = T(0);
        T wu = T(0);
        T wv = T(0);
        for (int d = 0; d < 3; ++d) {
            const T vt = ev[d] - sv[d];
            wt = sccd::max<T>(
                wt,
                sccd::max<T>(sccd::max<T>(sccd::abs<T>(vt - (e1[d] - s1[d])), sccd::abs<T>(vt - (e2[d] - s2[d]))),
                              sccd::max<T>(sccd::abs<T>(vt - (e3[d] - s3[d])), sccd::abs<T>(vt - (e4[d] - s4[d])))));
            wu = sccd::max<T>(wu,
                               sccd::max<T>(
                                   sccd::max<T>(sccd::abs<T>(s2[d] - s1[d]), sccd::abs<T>(s4[d] - s3[d])),
                                   sccd::max<T>(sccd::abs<T>(e2[d] - e1[d]), sccd::abs<T>(e4[d] - e3[d]))));
            wv = sccd::max<T>(wv,
                               sccd::max<T>(
                                   sccd::max<T>(sccd::abs<T>(s3[d] - s1[d]), sccd::abs<T>(s4[d] - s2[d])),
                                   sccd::max<T>(sccd::abs<T>(e3[d] - e1[d]), sccd::abs<T>(e4[d] - e2[d]))));
        }
        widths[0] = wt;
        widths[1] = wu;
        widths[2] = wv;
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
        // The certified numerical error bound, which is what makes the rejection
        // below sound. sccd_get_numerical_error_vq_soa already existed in this
        // file but was never called: the acceptance test padded with machine
        // epsilon instead, roughly 30x too small for unit-scale geometry, which
        // let it discard boxes that contained a root.
        T numerical_error[3];
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
                                           &numerical_error[0],
                                           &numerical_error[1],
                                           &numerical_error[2]);

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

        const auto stack_size = stack.size();
        bool found = false;
        for (int i = 0; i < N + 1; ++i) {
            const T sample_min = samples[i];
            const T sample_max = samples[i + 1];
            const T tt_min = SplitDim == 0 ? sample_min : domain.tuv[0].lower;
            const T tt_max = SplitDim == 0 ? sample_max : domain.tuv[0].upper;
            const T uu_min = SplitDim == 1 ? sample_min : domain.tuv[1].lower;
            const T uu_max = SplitDim == 1 ? sample_max : domain.tuv[1].upper;
            const T vv_min = SplitDim == 2 ? sample_min : domain.tuv[2].lower;
            const T vv_max = SplitDim == 2 ? sample_max : domain.tuv[2].upper;

            if (tt_min >= toi) {
                continue;
            }

            T fmin[3];
            T fmax[3];
            init_codomain_bounds<T>(fmin, fmax);

            for (int mask = 0; mask < 8; ++mask) {
                const T ct = (mask & 1) ? tt_max : tt_min;
                const T cu = (mask & 2) ? uu_max : uu_min;
                const T cv = (mask & 4) ? vv_max : vv_min;
                T F[3];
                diff_vq<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, ct, cu, cv, F);
                update_codomain_bounds<T>(F, fmin, fmax);
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
                domain, max_iter, tol, tols, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_adaptive_split_vq_axis<1, N, T>(
                domain, max_iter, tol, tols, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        return grid_search_adaptive_split_vq_axis<2, N, T>(
            domain, max_iter, tol, tols, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, toi, u, v, stack, refine);
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
        return grid_search_adaptive_split_vq<N, T>(
            domain, max_iter, tol, tols, codomain_widths, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_adaptive_split_vq(const int max_iter,
                                          const T tol,
                                          const T tols[3],
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

        return grid_search_adaptive_split_vq<ADAPTIVE_NUM_SPLITS, T>(
            initial_domain, max_iter, tol, tols, codomain_widths, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, stack, refine);
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
        return find_root_grid_adaptive_split_vq<T>(
            max_iter, tol, tols, codomain_widths, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
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
        return find_root_grid_adaptive_split_vq<T>(
            max_iter, tol, tols, codomain_widths, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
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
