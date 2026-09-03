#ifndef SCCD_SPIKE_UNIFORM_SPLIT_HPP
#define SCCD_SPIKE_UNIFORM_SPLIT_HPP

// Uniform interval splitting: a spike. See spikes/README.md for the rules.
//
// This is the branch-and-bound search splitting each interval at its midpoints
// rather than placing splitters from a Gauss-Newton step. It was a complete
// second implementation of the splitting job, for both vertex-face and
// edge-edge, reachable only through SCCD_ADAPTIVE_SPLIT=0.
//
// It was demoted because it never won. Measured end to end on three real scenes
// (benchmark/ASSESSMENT.md), adaptive is ahead on all three, and on cloth-funnel
// by 19% against a 6-10% run-to-run spread. The other two margins are inside the
// noise and prove nothing on their own -- but a duplicate that never wins does
// not meet the keep bar, and this one carried about 550 lines.
//
// It is kept because it is the alternative the adaptive splitter is measured
// against, and reading it is the quickest way to see what "adaptive" is adaptive
// relative to. It still compiles against the shipped headers -- spikes/ builds a
// translation unit that includes it, so it fails loudly rather than rotting --
// but nothing calls it: the dispatch that did is gone, and reviving it means
// putting that branch back.

#include "srootfinder.hpp"

namespace sccd {

    template <int SplitDim, int N, typename T>
    inline bool grid_search_uniform_split_vf_axis(const sccd::Box<T> &domain,
                                                  const int max_iter,
                                                  const T tol,
                                                  const T tols[3],
                                                  const T numerical_error[3],
                                                  const T sv[3],
                                                  const T s1[3],
                                                  const T s2[3],
                                                  const T s3[3],
                                                  const T ev[3],
                                                  const T e1[3],
                                                  const T e2[3],
                                                  const T e3[3],
                                                  T &toi,
                                                  T &u,
                                                  T &v,
                                                  std::vector<sccd::Box<T>> &stack,
                                                  const bool refine) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        static_assert(N > 0);

        const T lo = domain.tuv[SplitDim].lower;
        const T h = (domain.tuv[SplitDim].upper - lo) / T(N);

        alignas(64) T sample_min[N];
        alignas(64) T sample_max[N];
        alignas(64) T fmin[3][N];
        alignas(64) T fmax[3][N];
        alignas(64) bool contains_zero[N];
        alignas(64) bool accept[N];

#pragma omp simd aligned(sample_min, sample_max)
        for (int i = 0; i < N; ++i) {
            sample_min[i] = lo + h * T(i);
            sample_max[i] = lo + h * T(i + 1);
        }

        for (int d = 0; d < 3; ++d) {
#pragma omp simd aligned(fmin, fmax)
            for (int i = 0; i < N; ++i) {
                fmin[d][i] = std::numeric_limits<T>::max();
                fmax[d][i] = std::numeric_limits<T>::lowest();
            }
        }

        for (int mask = 0; mask < 8; ++mask) {
            const bool mt = mask & 1;
            const bool mu = mask & 2;
            const bool mv = mask & 4;
#pragma omp simd aligned(sample_min, sample_max, fmin, fmax)
            for (int i = 0; i < N; ++i) {
                const T tt = SplitDim == 0 ? (mt ? sample_max[i] : sample_min[i])
                                           : (mt ? domain.tuv[0].upper : domain.tuv[0].lower);
                const T uu = SplitDim == 1 ? (mu ? sample_max[i] : sample_min[i])
                                           : (mu ? domain.tuv[1].upper : domain.tuv[1].lower);
                const T vv = SplitDim == 2 ? (mv ? sample_max[i] : sample_min[i])
                                           : (mv ? domain.tuv[2].upper : domain.tuv[2].lower);
                T F[3];
                diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, tt, uu, vv, F);
                for (int d = 0; d < 3; ++d) {
                    fmin[d][i] = sccd::min<T>(fmin[d][i], F[d]);
                    fmax[d][i] = sccd::max<T>(fmax[d][i], F[d]);
                }
            }
        }

#pragma omp simd aligned(fmin, fmax, contains_zero, accept)
        for (int i = 0; i < N; ++i) {
            const T fmin_i[3] = {fmin[0][i], fmin[1][i], fmin[2][i]};
            const T fmax_i[3] = {fmax[0][i], fmax[1][i], fmax[2][i]};
            contains_zero[i] = codomain_acceptance<T>(fmin_i, fmax_i, tol, tols, numerical_error, accept[i]);
        }

        auto stack_size = stack.size();

        bool found = false;
        for (int i = 0; i < N; ++i) {
            if (!contains_zero[i]) {
                continue;
            }

            const T sample_lo = sample_min[i];
            const T sample_hi = sample_max[i];
            const T tt_min = SplitDim == 0 ? sample_lo : domain.tuv[0].lower;
            const T tt_max = SplitDim == 0 ? sample_hi : domain.tuv[0].upper;
            const T uu_min = SplitDim == 1 ? sample_lo : domain.tuv[1].lower;
            const T uu_max = SplitDim == 1 ? sample_hi : domain.tuv[1].upper;
            const T vv_min = SplitDim == 2 ? sample_lo : domain.tuv[2].lower;
            const T vv_max = SplitDim == 2 ? sample_hi : domain.tuv[2].upper;

            if (tt_min >= toi || uu_min + vv_min >= T(1) + tol) {
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
            accept[i] = accept[i] && (tt_min > 0);

            Box<T> box = split_axis_box<SplitDim, T>(domain, sample_lo, sample_hi);
            if (accept[i] || box.depth > max_iter) {
                found |= accept_grid_root_vf<T>(box, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, refine);
                continue;
            }

            // box.tuv[0].lower = std::min(box.tuv[0].lower, toi);
            stack.push_back(box);
        }

        // Make sure the tmin is on top of the stack
        if constexpr (SplitDim == 0) {
            std::reverse(stack.begin() + stack_size, stack.end());
        }

        return found;
    }

    template <int N, typename T>
    inline bool grid_search_uniform_split_vf(const sccd::Box<T> &domain,
                                             const int max_iter,
                                             const T tol,
                                             const T tols[3],
                                             const T numerical_error[3],
                                             const T codomain_widths[3],
                                             const T sv[3],
                                             const T s1[3],
                                             const T s2[3],
                                             const T s3[3],
                                             const T ev[3],
                                             const T e1[3],
                                             const T e2[3],
                                             const T e3[3],
                                             T &toi,
                                             T &u,
                                             T &v,
                                             std::vector<sccd::Box<T>> &stack,
                                             const bool refine) {
        const int split_dim = domain.widest_dimension(codomain_widths);
        if (split_dim == 0) {
            return grid_search_uniform_split_vf_axis<0, N, T>(
                domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_uniform_split_vf_axis<1, N, T>(
                domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
        }
        return grid_search_uniform_split_vf_axis<2, N, T>(
            domain, max_iter, tol, tols, numerical_error, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
    }

    template <int N, typename T>
    inline bool grid_search_uniform_split_vf(const sccd::Box<T> &domain,
                                             const int max_iter,
                                             const T tol,
                                             const T tols[3],
                                             const T sv[3],
                                             const T s1[3],
                                             const T s2[3],
                                             const T s3[3],
                                             const T ev[3],
                                             const T e1[3],
                                             const T e2[3],
                                             const T e3[3],
                                             T &toi,
                                             T &u,
                                             T &v,
                                             std::vector<sccd::Box<T>> &stack,
                                             const bool refine) {
        const T codomain_widths[3] = {T(1), T(1), T(1)};
        T numerical_error[3];
        numerical_error_bound<true, T>(sv, s1, s2, s3, ev, e1, e2, e3, numerical_error);
        return grid_search_uniform_split_vf<N, T>(
            domain, max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, ev, e1, e2, e3, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_vf(const int max_iter,
                                         const T tol,
                                         const T tols[3],
                                         const T numerical_error[3],
                                         const T codomain_widths[3],
                                         const T sv[3],
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T ev[3],
                                         const T e1[3],
                                         const T e2[3],
                                         const T e3[3],
                                         const Box<T> &initial_domain,
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        if (initial_domain.tuv[0].lower >= t) {
            return false;
        }

        return grid_search_uniform_split_vf<UNIFORM_NUM_SPLITS, T>(
            initial_domain, max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_vf(const int max_iter,
                                         const T tol,
                                         const T tols[3],
                                         const T sv[3],
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T ev[3],
                                         const T e1[3],
                                         const T e2[3],
                                         const T e3[3],
                                         const Box<T> &initial_domain,
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        const T codomain_widths[3] = {T(1), T(1), T(1)};
        T numerical_error[3];
        numerical_error_bound<true, T>(sv, s1, s2, s3, ev, e1, e2, e3, numerical_error);
        return find_root_grid_uniform_split_vf<T>(
            max_iter, tol, tols, numerical_error, codomain_widths, sv, s1, s2, s3, ev, e1, e2, e3, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_vf(const int max_iter,
                                         const T tol,
                                         const T sv[3],
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T ev[3],
                                         const T e1[3],
                                         const T e2[3],
                                         const T e3[3],
                                         const Box<T> &initial_domain,
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        T tols[3];
        compute_face_vertex_tolerance<T>(tol, sv, s1, s2, s3, ev, e1, e2, e3, tols);
        return find_root_grid_uniform_split_vf<T>(
            max_iter, tol, tols, sv, s1, s2, s3, ev, e1, e2, e3, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_vf(const int max_iter,
                                         const T tol,
                                         const T sv[3],
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T ev[3],
                                         const T e1[3],
                                         const T e2[3],
                                         const T e3[3],
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        return find_root_grid_uniform_split_vf<T>(
            max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, unit_domain_box<T>(), t, u, v, stack, refine);
    }

    template <int SplitDim, int N, typename T>
    inline bool grid_search_uniform_split_ee_axis(const sccd::Box<T> &domain,
                                                  const int max_iter,
                                                  const T tol,
                                                  const T tols[3],
                                                  const T numerical_error[3],
                                                  const T s1[3],
                                                  const T s2[3],
                                                  const T s3[3],
                                                  const T s4[3],
                                                  const T e1[3],
                                                  const T e2[3],
                                                  const T e3[3],
                                                  const T e4[3],
                                                  T &toi,
                                                  T &u,
                                                  T &v,
                                                  std::vector<sccd::Box<T>> &stack,
                                                  const bool refine) {
        static_assert(SplitDim >= 0 && SplitDim < 3);
        static_assert(N > 0);
        (void)refine;

        const T lo = domain.tuv[SplitDim].lower;
        const T h = (domain.tuv[SplitDim].upper - lo) / T(N);

        alignas(64) T sample_min[N];
        alignas(64) T sample_max[N];
        alignas(64) T fmin[3][N];
        alignas(64) T fmax[3][N];
        alignas(64) bool contains_zero[N];
        alignas(64) bool accept[N];

#pragma omp simd
        for (int i = 0; i < N; ++i) {
            sample_min[i] = lo + h * T(i);
            sample_max[i] = lo + h * T(i + 1);
        }

        for (int d = 0; d < 3; ++d) {
#pragma omp simd
            for (int i = 0; i < N; ++i) {
                fmin[d][i] = std::numeric_limits<T>::max();
                fmax[d][i] = std::numeric_limits<T>::lowest();
            }
        }

        for (int mask = 0; mask < 8; ++mask) {
            const bool mt = mask & 1;
            const bool mu = mask & 2;
            const bool mv = mask & 4;
#pragma omp simd
            for (int i = 0; i < N; ++i) {
                const T tt = SplitDim == 0 ? (mt ? sample_max[i] : sample_min[i])
                                           : (mt ? domain.tuv[0].upper : domain.tuv[0].lower);
                const T uu = SplitDim == 1 ? (mu ? sample_max[i] : sample_min[i])
                                           : (mu ? domain.tuv[1].upper : domain.tuv[1].lower);
                const T vv = SplitDim == 2 ? (mv ? sample_max[i] : sample_min[i])
                                           : (mv ? domain.tuv[2].upper : domain.tuv[2].lower);
                T F[3];
                diff_ee<T>(s1, s2, s3, s4, e1, e2, e3, e4, tt, uu, vv, F);
                for (int d = 0; d < 3; ++d) {
                    fmin[d][i] = sccd::min<T>(fmin[d][i], F[d]);
                    fmax[d][i] = sccd::max<T>(fmax[d][i], F[d]);
                }
            }
        }

#pragma omp simd
        for (int i = 0; i < N; ++i) {
            const T fmin_i[3] = {fmin[0][i], fmin[1][i], fmin[2][i]};
            const T fmax_i[3] = {fmax[0][i], fmax[1][i], fmax[2][i]};
            contains_zero[i] = codomain_acceptance<T>(fmin_i, fmax_i, tol, tols, numerical_error, accept[i]);
        }

        bool found = false;
        for (int i = 0; i < N; ++i) {
            if (!contains_zero[i]) {
                continue;
            }

            const T sample_lo = sample_min[i];
            const T sample_hi = sample_max[i];
            const T tt_min = SplitDim == 0 ? sample_lo : domain.tuv[0].lower;
            const T tt_max = SplitDim == 0 ? sample_hi : domain.tuv[0].upper;
            const T uu_min = SplitDim == 1 ? sample_lo : domain.tuv[1].lower;
            const T uu_max = SplitDim == 1 ? sample_hi : domain.tuv[1].upper;
            const T vv_min = SplitDim == 2 ? sample_lo : domain.tuv[2].lower;
            const T vv_max = SplitDim == 2 ? sample_hi : domain.tuv[2].upper;

            if (tt_min >= toi) {
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
            accept[i] = accept[i] && (tt_min > 0);

            Box<T> box = split_axis_box<SplitDim, T>(domain, sample_lo, sample_hi);
            if (accept[i] || box.depth > max_iter) {
                found |= accept_grid_root_ee<T>(box, toi, u, v);
                continue;
            }

            // box.tuv[0].lower = std::min(box.tuv[0].lower, toi);
            stack.push_back(box);
        }

        return found;
    }

    template <int N, typename T>
    inline bool grid_search_uniform_split_ee(const sccd::Box<T> &domain,
                                             const int max_iter,
                                             const T tol,
                                             const T tols[3],
                                             const T numerical_error[3],
                                             const T codomain_widths[3],
                                             const T s1[3],
                                             const T s2[3],
                                             const T s3[3],
                                             const T s4[3],
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
            return grid_search_uniform_split_ee_axis<0, N, T>(
                domain, max_iter, tol, tols, numerical_error, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        if (split_dim == 1) {
            return grid_search_uniform_split_ee_axis<1, N, T>(
                domain, max_iter, tol, tols, numerical_error, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
        }
        return grid_search_uniform_split_ee_axis<2, N, T>(
            domain, max_iter, tol, tols, numerical_error, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
    }

    template <int N, typename T>
    inline bool grid_search_uniform_split_ee(const sccd::Box<T> &domain,
                                             const int max_iter,
                                             const T tol,
                                             const T tols[3],
                                             const T s1[3],
                                             const T s2[3],
                                             const T s3[3],
                                             const T s4[3],
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
        numerical_error_bound<false, T>(s1, s2, s3, s4, e1, e2, e3, e4, numerical_error);
        return grid_search_uniform_split_ee<N, T>(
            domain, max_iter, tol, tols, numerical_error, codomain_widths, s1, s2, s3, s4, e1, e2, e3, e4, toi, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_ee(const int max_iter,
                                         const T tol,
                                         const T tols[3],
                                         const T numerical_error[3],
                                         const T codomain_widths[3],
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T s4[3],
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

        return grid_search_uniform_split_ee<UNIFORM_NUM_SPLITS, T>(
            initial_domain, max_iter, tol, tols, numerical_error, codomain_widths, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_ee(const int max_iter,
                                         const T tol,
                                         const T tols[3],
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T s4[3],
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
        numerical_error_bound<false, T>(s1, s2, s3, s4, e1, e2, e3, e4, numerical_error);
        return find_root_grid_uniform_split_ee<T>(
            max_iter, tol, tols, numerical_error, codomain_widths, s1, s2, s3, s4, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_ee(const int max_iter,
                                         const T tol,
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T s4[3],
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
        compute_edge_edge_tolerance<T>(tol, s1, s2, s3, s4, e1, e2, e3, e4, tols);
        return find_root_grid_uniform_split_ee<T>(
            max_iter, tol, tols, s1, s2, s3, s4, e1, e2, e3, e4, initial_domain, t, u, v, stack, refine);
    }

    template <typename T>
    bool find_root_grid_uniform_split_ee(const int max_iter,
                                         const T tol,
                                         const T s1[3],
                                         const T s2[3],
                                         const T s3[3],
                                         const T s4[3],
                                         const T e1[3],
                                         const T e2[3],
                                         const T e3[3],
                                         const T e4[3],
                                         T &t,
                                         T &u,
                                         T &v,
                                         std::vector<Box<T>> &stack,
                                         const bool refine = false) {
        return find_root_grid_uniform_split_ee<T>(
            max_iter, tol, s1, s2, s3, s4, e1, e2, e3, e4, unit_domain_box<T>(), t, u, v, stack, refine);
    }

}  // namespace sccd

#endif  // SCCD_SPIKE_UNIFORM_SPLIT_HPP
