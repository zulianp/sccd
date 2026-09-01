#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP
// #include <tbb/blocked_range.h>
// #include <tbb/parallel_for.h>

#include "assert.h"

#include "roots.hpp"
#include "sccd_base.hpp"
#include "sccd_vnarrowphase.hpp"
#include "srootfinder.hpp"
#include "vaabb.hpp"

#include <atomic>

#define SCCD_ENABLE_CODOMAIN_SCALING 1

namespace sccd {

    template <typename T>
    T atomic_min(std::atomic<T>& atm, T val) {
        T expected = atm.load();
        while (expected > val) {
            // If 'atm' is still equal to 'expected', set it to 'val'.
            // If not, 'expected' is updated to the actual current value of 'atm', and the loop retries.
            if (atm.compare_exchange_strong(expected, val)) {
                break;
            }
            // No need for a separate load() here; compare_exchange_strong updates 'expected' on failure.
        }
        return expected;  // Returns the value 'expected' held when the operation succeeded (the prior minimum).
    }

    template <int nxe, typename T, typename I>
    int narrow_phase_newton_pass_vf(const size_t noverlaps,
                                    const I* const SCCD_RESTRICT voveralp,
                                    const I* const SCCD_RESTRICT foveralp,
                                    T** const SCCD_RESTRICT v0,
                                    T** const SCCD_RESTRICT v1,
                                    const size_t face_stride,
                                    I** const SCCD_RESTRICT faces,
                                    const T max_toi,
                                    T* const SCCD_RESTRICT toi,
                                    const T tol,
                                    const int toi_stride = 0);

    template <int nxe, typename T, typename I>
    int narrow_phase_vf(const size_t noverlaps,
                        const I* const SCCD_RESTRICT voveralp,
                        const I* const SCCD_RESTRICT foveralp,
                        // Geometric data
                        T** const SCCD_RESTRICT v0,
                        T** const SCCD_RESTRICT v1,
                        const size_t face_stride,
                        I** const SCCD_RESTRICT faces,
                        const T max_toi,
                        T* const SCCD_RESTRICT toi,
                        const int max_depth,
                        const T tol,
                        const int toi_stride = 0) {
        using T_HP = double;
        const T infty = 100000;

        assert(toi_stride == 0 || toi_stride == 1);
        if (noverlaps == 0) {
            if (toi != nullptr && toi_stride == 0) toi[0] = max_toi;
            return 0;
        }
        assert(toi != nullptr);

        int SCCD_USE_VNARROW_PHASE = 0;
        SCCD_READ_ENV(SCCD_USE_VNARROW_PHASE, atoi);

        if (SCCD_USE_VNARROW_PHASE) {
            return v_narrow_phase_vf<nxe, T, I>(
                noverlaps, voveralp, foveralp, v0, v1, face_stride, faces, max_toi, toi, max_depth, tol, toi_stride);
        }

        int SCCD_USE_TI = 0;
        SCCD_READ_ENV(SCCD_USE_TI, atoi);

        int SCCD_REFINE = 0;
        SCCD_READ_ENV(SCCD_REFINE, atoi);

        int SCCD_ADAPTIVE_SPLIT = 1;
        SCCD_READ_ENV(SCCD_ADAPTIVE_SPLIT, atoi);

        std::atomic<T> min_t = max_toi;

        if (toi_stride == 1) {
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < noverlaps; i++) {
                toi[i] = max_toi;
            }
        }

        if (toi_stride == 0 && min_t == 0) {
            toi[0] = 0;
            // printf("[after newton pass vf] min_t is 0 returning\n");
            return 0;
        }

        if (toi_stride == 0) toi[0] = max_toi;
        // sccd::parallel_for_br_dynamic(0, noverlaps, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
        // std::vector<Box<T_HP>> stack;

#pragma omp parallel
        {
            std::vector<Box<T_HP>> stack;
            stack.reserve(64);

#pragma omp for schedule(dynamic, 64) nowait
            for (ptrdiff_t i = 0; i < noverlaps; i++) {
                if (toi_stride == 1) toi[i] = max_toi;

                const I vi = voveralp[i];
                const I fi = foveralp[i];

                I nodes[3] = {faces[0][fi * face_stride], faces[1][fi * face_stride], faces[2][fi * face_stride]};

                const T_HP sv[3] = {v0[0][vi], v0[1][vi], v0[2][vi]};
                const T_HP ev[3] = {v1[0][vi], v1[1][vi], v1[2][vi]};

                const T_HP s1[3] = {v0[0][nodes[0]], v0[1][nodes[0]], v0[2][nodes[0]]};
                const T_HP s2[3] = {v0[0][nodes[1]], v0[1][nodes[1]], v0[2][nodes[1]]};
                const T_HP s3[3] = {v0[0][nodes[2]], v0[1][nodes[2]], v0[2][nodes[2]]};

                const T_HP e1[3] = {v1[0][nodes[0]], v1[1][nodes[0]], v1[2][nodes[0]]};
                const T_HP e2[3] = {v1[0][nodes[1]], v1[1][nodes[1]], v1[2][nodes[1]]};
                const T_HP e3[3] = {v1[0][nodes[2]], v1[1][nodes[2]], v1[2][nodes[2]]};

                // Iteration variables
                T_HP t = toi_stride == 0 ? T_HP(min_t.load()) : T_HP(toi[i]);
                T_HP u = 0;
                T_HP v = 0;

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#warning "SCCD_ENABLE_TIGHT_INCLUSION"
                if (SCCD_USE_TI) {
                    if (find_root_tight_inclusion_vf<T_HP>(max_depth, tol, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
                        if (toi_stride == 0) {
                            atomic_min<T>(min_t, t);
                        } else {
                            toi[i] = t;
                        }
                    }
                    continue;
                }
#endif

                const Box<T_HP> initial_domain(Interval<T_HP>{T_HP(0), sccd::min<T_HP>(t, T_HP(1))},
                                               Interval<T_HP>{T_HP(0), T_HP(1)},
                                               Interval<T_HP>{T_HP(0), T_HP(1)},
                                               0);

#ifndef SCCD_ENABLE_CODOMAIN_SCALING
                T_HP codomain_widths[3] = {1, 1, 1};
#else
                T_HP codomain_widths[3];
                codomain_widths[0] = T_HP(0);
                codomain_widths[1] = T_HP(0);
                codomain_widths[2] = T_HP(0);
                for (int d = 0; d < 3; ++d) {
                    const T_HP vv = ev[d] - sv[d];
                    const T_HP p1 = e1[d] - s1[d];
                    const T_HP p2 = e2[d] - s2[d];
                    const T_HP p3 = e3[d] - s3[d];

                    codomain_widths[0] = sccd::max<T_HP>(
                        codomain_widths[0],
                        sccd::max<T_HP>(sccd::max<T_HP>(sccd::abs<T_HP>(vv - p1), sccd::abs<T_HP>(vv - p2)),
                                        sccd::max<T_HP>(sccd::abs<T_HP>(vv - p3), sccd::abs<T_HP>(vv + p1 - p2 - p3))));

                    const T_HP su = s2[d] - s1[d];
                    const T_HP eu = e2[d] - e1[d];
                    const T_HP u_upper = su + (eu - su) * t;
                    codomain_widths[1] = sccd::max<T_HP>(
                        codomain_widths[1], sccd::max<T_HP>(sccd::abs<T_HP>(su), sccd::abs<T_HP>(u_upper)));

                    const T_HP sv0 = s3[d] - s1[d];
                    const T_HP ev0 = e3[d] - e1[d];
                    const T_HP v_upper = sv0 + (ev0 - sv0) * t;
                    codomain_widths[2] = sccd::max<T_HP>(
                        codomain_widths[2], sccd::max<T_HP>(sccd::abs<T_HP>(sv0), sccd::abs<T_HP>(v_upper)));
                }

                T_HP tot_widths = codomain_widths[0] + codomain_widths[1] + codomain_widths[2] + 1e-16;

                for (int d = 0; d < 3; d++) {
                    codomain_widths[d] /= tot_widths;
                    codomain_widths[d] = sccd::max<T_HP>(1e-4, codomain_widths[d]);
                }
#endif

                T_HP tols[3];
                compute_face_vertex_tolerance<T_HP>(tol, sv, s1, s2, s3, ev, e1, e2, e3, tols);
                T_HP numerical_error[3];
                numerical_error_bound<true, T_HP>(sv, s1, s2, s3, ev, e1, e2, e3, numerical_error);

                stack.clear();
                stack.push_back(initial_domain);

                while (!stack.empty()) {
                    Box<T_HP> box = stack.back();
                    stack.pop_back();

                    if (box.tuv[0].lower >= t) {
                        continue;
                    }

                    box.tuv[0].upper = sccd::min<T_HP>(box.tuv[0].upper, t);

                    bool found = false;
                    if (SCCD_ADAPTIVE_SPLIT) {
                        found = find_root_grid_adaptive_split_vf<T_HP>(max_depth,
                                                                       tol,
                                                                       tols,
                                                                       numerical_error,
                                                                       codomain_widths,
                                                                       sv,
                                                                       s1,
                                                                       s2,
                                                                       s3,
                                                                       ev,
                                                                       e1,
                                                                       e2,
                                                                       e3,
                                                                       box,
                                                                       t,
                                                                       u,
                                                                       v,
                                                                       stack,
                                                                       SCCD_REFINE);
                    } else {
                        found = find_root_grid_uniform_split_vf<T_HP>(max_depth,
                                                                      tol,
                                                                      tols,
                                                                      numerical_error,
                                                                      codomain_widths,
                                                                      sv,
                                                                      s1,
                                                                      s2,
                                                                      s3,
                                                                      ev,
                                                                      e1,
                                                                      e2,
                                                                      e3,
                                                                      box,
                                                                      t,
                                                                      u,
                                                                      v,
                                                                      stack,
                                                                      SCCD_REFINE);
                    }

                    if (found) {
                        if (toi_stride == 0) {
                            auto ret = atomic_min<T>(min_t, t);
                            if (ret < t) {
                                t = ret;
                            }

                        } else {
                            toi[i] = t;
                        }
                    } else if (!stack.empty()) {
                        if (toi_stride == 0) {
                            t = sccd::min<T_HP>(t, min_t.load());
                        }
                    }
                }
            }

            // printf("VF max capacity: %zu\n", stack.capacity());
        }
        // );

        if (toi_stride == 0) toi[0] = min_t;
        return 0;
    }

    // Legacy direct-call overload: the old API only exposed toi_stride here.
    template <int nxe, typename T, typename I>
    int narrow_phase_vf(const size_t noverlaps,
                        const I* const SCCD_RESTRICT voveralp,
                        const I* const SCCD_RESTRICT foveralp,
                        T** const SCCD_RESTRICT v0,
                        T** const SCCD_RESTRICT v1,
                        const size_t face_stride,
                        I** const SCCD_RESTRICT faces,
                        const T max_toi,
                        T* const SCCD_RESTRICT toi,
                        const int toi_stride = 0) {
        int SCCD_MAX_DEPTH = 69;
        SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);

        T SCCD_TOL = T(3e-8);
        SCCD_READ_ENV(SCCD_TOL, atof);

        return narrow_phase_vf<nxe, T, I>(noverlaps,
                                          voveralp,
                                          foveralp,
                                          v0,
                                          v1,
                                          face_stride,
                                          faces,
                                          max_toi,
                                          toi,
                                          SCCD_MAX_DEPTH,
                                          SCCD_TOL,
                                          toi_stride);
    }

    template <typename T, typename I>
    int narrow_phase_ee(const size_t noverlaps,
                        const I* const SCCD_RESTRICT e0overalp,
                        const I* const SCCD_RESTRICT e1overalp,
                        // Geometric data
                        T** const SCCD_RESTRICT v0,
                        T** const SCCD_RESTRICT v1,
                        const size_t edge_stride,
                        I** const SCCD_RESTRICT edges,
                        // Output
                        const T max_toi,
                        T* const SCCD_RESTRICT toi,
                        const int max_depth,
                        const T tol,
                        const int toi_stride = 0) {
        using T_HP = double;
        const T infty = 100000;

        if (toi_stride == 0 && max_toi == 0) {
            toi[0] = 0;
            // printf("max_toi is 0 returning\n");
            return 0;
        }

        if (toi_stride == 1) {
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < noverlaps; i++) {
                toi[i] = max_toi;
            }
        }

        assert(toi_stride == 0 || toi_stride == 1);
        if (noverlaps == 0) {
            if (toi != nullptr && toi_stride == 0) toi[0] = max_toi;
            return 0;
        }
        assert(toi != nullptr);

        int SCCD_USE_TI = 0;
        SCCD_READ_ENV(SCCD_USE_TI, atoi);

        int SCCD_REFINE = 0;
        SCCD_READ_ENV(SCCD_REFINE, atoi);

        int SCCD_ADAPTIVE_SPLIT = 1;
        SCCD_READ_ENV(SCCD_ADAPTIVE_SPLIT, atoi);

        std::atomic<T> min_t = max_toi;
        if (toi_stride == 0) toi[0] = max_toi;

        // sccd::parallel_for_br_dynamic(0, noverlaps, [&](const ptrdiff_t rbegin, const ptrdiff_t rend)

#pragma omp parallel
        {
            std::vector<Box<T_HP>> stack;
            stack.reserve(64);

#pragma omp for schedule(dynamic, 64) nowait
            for (ptrdiff_t i = 0; i < noverlaps; i++) {
                if (toi_stride == 1) toi[i] = max_toi;

                const I i0 = e0overalp[i];
                const I i1 = e1overalp[i];

                I nodes0[2] = {edges[0][i0 * edge_stride], edges[1][i0 * edge_stride]};
                I nodes1[2] = {edges[0][i1 * edge_stride], edges[1][i1 * edge_stride]};

                const T_HP s1[3] = {v0[0][nodes0[0]], v0[1][nodes0[0]], v0[2][nodes0[0]]};
                const T_HP s2[3] = {v0[0][nodes0[1]], v0[1][nodes0[1]], v0[2][nodes0[1]]};

                const T_HP s3[3] = {v0[0][nodes1[0]], v0[1][nodes1[0]], v0[2][nodes1[0]]};
                const T_HP s4[3] = {v0[0][nodes1[1]], v0[1][nodes1[1]], v0[2][nodes1[1]]};

                const T_HP e1[3] = {v1[0][nodes0[0]], v1[1][nodes0[0]], v1[2][nodes0[0]]};
                const T_HP e2[3] = {v1[0][nodes0[1]], v1[1][nodes0[1]], v1[2][nodes0[1]]};

                const T_HP e3[3] = {v1[0][nodes1[0]], v1[1][nodes1[0]], v1[2][nodes1[0]]};
                const T_HP e4[3] = {v1[0][nodes1[1]], v1[1][nodes1[1]], v1[2][nodes1[1]]};

                // Iteration variables
                T_HP t = toi_stride == 0 ? T_HP(min_t.load()) : T_HP(toi[i]);
                T_HP u = 0;
                T_HP v = 0;
                const T_HP t_upper = sccd::min<T_HP>(t, T_HP(1));
                const Box<T_HP> initial_domain(Interval<T_HP>{T_HP(0), t_upper},
                                               Interval<T_HP>{T_HP(0), T_HP(1)},
                                               Interval<T_HP>{T_HP(0), T_HP(1)},
                                               0);

#ifndef SCCD_ENABLE_CODOMAIN_SCALING
                T_HP codomain_widths[3] = {1, 1, 1};
#else
                T_HP codomain_widths[3];
                codomain_widths[0] = T_HP(0);
                codomain_widths[1] = T_HP(0);
                codomain_widths[2] = T_HP(0);
                for (int d = 0; d < 3; ++d) {
                    const T_HP a0 = e1[d] - s1[d];
                    const T_HP a1 = e2[d] - s2[d];
                    const T_HP b0 = e3[d] - s3[d];
                    const T_HP b1 = e4[d] - s4[d];

                    codomain_widths[0] = sccd::max<T_HP>(
                        codomain_widths[0],
                        sccd::max<T_HP>(sccd::max<T_HP>(sccd::abs<T_HP>(a0 - b0), sccd::abs<T_HP>(a0 - b1)),
                                        sccd::max<T_HP>(sccd::abs<T_HP>(a1 - b0), sccd::abs<T_HP>(a1 - b1))));

                    const T_HP su = s2[d] - s1[d];
                    const T_HP eu = e2[d] - e1[d];
                    const T_HP u_upper = su + (eu - su) * t_upper;
                    codomain_widths[1] = sccd::max<T_HP>(
                        codomain_widths[1], sccd::max<T_HP>(sccd::abs<T_HP>(su), sccd::abs<T_HP>(u_upper)));

                    const T_HP sv = s4[d] - s3[d];
                    const T_HP ev = e4[d] - e3[d];
                    const T_HP v_upper = sv + (ev - sv) * t_upper;
                    codomain_widths[2] = sccd::max<T_HP>(
                        codomain_widths[2], sccd::max<T_HP>(sccd::abs<T_HP>(sv), sccd::abs<T_HP>(v_upper)));
                }

                T_HP tot_widths = codomain_widths[0] + codomain_widths[1] + codomain_widths[2] + 1e-16;

                for (int d = 0; d < 3; d++) {
                    codomain_widths[d] /= tot_widths;
                    codomain_widths[d] = sccd::max<T_HP>(1e-4, codomain_widths[d]);
                }

                // printf("%g %g %g\n", codomain_widths[0], codomain_widths[1], codomain_widths[2]);
#endif

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#warning "SCCD_ENABLE_TIGHT_INCLUSION"
                if (SCCD_USE_TI) {
                    if (find_root_tight_inclusion_ee<T_HP>(max_depth, tol, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v)) {
                        if (toi_stride == 0) {
                            atomic_min<T>(min_t, t);
                        } else {
                            toi[i] = t;
                        }
                    }
                    continue;
                }
#endif
                T_HP tols[3];
                compute_edge_edge_tolerance<T_HP>(tol, s1, s2, s3, s4, e1, e2, e3, e4, tols);
                T_HP numerical_error[3];
                numerical_error_bound<false, T_HP>(s1, s2, s3, s4, e1, e2, e3, e4, numerical_error);

                stack.clear();
                stack.push_back(initial_domain);
                while (!stack.empty()) {
                    Box<T_HP> box = stack.back();
                    stack.pop_back();

                    if (box.tuv[0].lower >= t) {
                        continue;
                    }

                    box.tuv[0].upper = sccd::min<T_HP>(box.tuv[0].upper, t);

                    bool found = false;
                    if (SCCD_ADAPTIVE_SPLIT) {
                        found = find_root_grid_adaptive_split_ee<T_HP>(max_depth,
                                                                       tol,
                                                                       tols,
                                                                       numerical_error,
                                                                       codomain_widths,
                                                                       s1,
                                                                       s2,
                                                                       s3,
                                                                       s4,
                                                                       e1,
                                                                       e2,
                                                                       e3,
                                                                       e4,
                                                                       box,
                                                                       t,
                                                                       u,
                                                                       v,
                                                                       stack,
                                                                       SCCD_REFINE);
                    } else {
                        found = find_root_grid_uniform_split_ee<T_HP>(max_depth,
                                                                      tol,
                                                                      tols,
                                                                      numerical_error,
                                                                      codomain_widths,
                                                                      s1,
                                                                      s2,
                                                                      s3,
                                                                      s4,
                                                                      e1,
                                                                      e2,
                                                                      e3,
                                                                      e4,
                                                                      box,
                                                                      t,
                                                                      u,
                                                                      v,
                                                                      stack,
                                                                      SCCD_REFINE);
                    }

                    if (found) {
                        if (toi_stride == 0) {
                            auto ret = atomic_min<T>(min_t, t);
                            if (ret < t) {
                                t = ret;
                            }

                        } else {
                            toi[i] = t;
                        }
                    } else if (!stack.empty()) {
                        if (toi_stride == 0) {
                            t = sccd::min<T_HP>(t, min_t.load());
                        }
                    }
                }
            }

            // printf("EE max capacity: %zu\n", stack.capacity());
        }
        // );

        if (toi_stride == 0) toi[0] = min_t;
        return 0;
    }

    // Legacy direct-call overload: the old API only exposed toi_stride here.
    template <typename T, typename I>
    int narrow_phase_ee(const size_t noverlaps,
                        const I* const SCCD_RESTRICT e0overalp,
                        const I* const SCCD_RESTRICT e1overalp,
                        T** const SCCD_RESTRICT v0,
                        T** const SCCD_RESTRICT v1,
                        const size_t edge_stride,
                        I** const SCCD_RESTRICT edges,
                        const T max_toi,
                        T* const SCCD_RESTRICT toi,
                        const int toi_stride = 0) {
        int SCCD_MAX_DEPTH = 69;
        SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);

        T SCCD_TOL = T(3e-8);
        SCCD_READ_ENV(SCCD_TOL, atof);

        return narrow_phase_ee<T, I>(noverlaps,
                                     e0overalp,
                                     e1overalp,
                                     v0,
                                     v1,
                                     edge_stride,
                                     edges,
                                     max_toi,
                                     toi,
                                     SCCD_MAX_DEPTH,
                                     SCCD_TOL,
                                     toi_stride);
    }

}  // namespace sccd

#endif
