#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP
// #include <tbb/blocked_range.h>
// #include <tbb/parallel_for.h>

#include "assert.h"

#include "roots.hpp"
#include "sccd_base.hpp"
#include "srootfinder.hpp"
#include "vaabb.hpp"

#include <atomic>

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
                        const int toi_stride = 0) {
        using T_HP = double;
        const T infty = 100000;

        assert(toi_stride == 0 || toi_stride == 1);
        if (noverlaps == 0) {
            if (toi != nullptr && toi_stride == 0) toi[0] = max_toi;
            return 0;
        }
        assert(toi != nullptr);

        int SCCD_USE_TI = 0;
        SCCD_READ_ENV(SCCD_USE_TI, atoi);

        int SCCD_MAX_DEPTH = 32;
        SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);

        T_HP SCCD_TOL = std::is_same_v<float, T_HP> ? T(1e-8) : T(1e-14);
        SCCD_READ_ENV(SCCD_TOL, atof);

        int SCCD_REFINE = 0;
        SCCD_READ_ENV(SCCD_REFINE, atoi);

        int SCCD_ADAPTIVE_SPLIT = 0;
        SCCD_READ_ENV(SCCD_ADAPTIVE_SPLIT, atoi);

        std::atomic<T> min_t = max_toi;

        int SCCD_USE_NEWTON_PASS = 0;
        SCCD_READ_ENV(SCCD_USE_NEWTON_PASS, atoi);
        if (SCCD_USE_NEWTON_PASS) {
            narrow_phase_newton_pass_vf<nxe, T, I>(
                noverlaps, voveralp, foveralp, v0, v1, face_stride, faces, max_toi, toi, toi_stride);

            if (toi_stride == 0) {
                min_t = sccd::min<T>(min_t, T(toi[0]));
            }
        }

        if (toi_stride == 0 && min_t == 0) {
            toi[0] = 0;
            printf("[after newton pass vf] min_t is 0 returning\n");
            return 0;
        }

        if (toi_stride == 0) toi[0] = max_toi;
        sccd::parallel_for_br(0, noverlaps, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            std::vector<Box<T_HP>> stack;

            for (ptrdiff_t i = rbegin; i < rend; i++) {
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
                T_HP t = min_t;
                T_HP u = 0;
                T_HP v = 0;

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#warning "SCCD_ENABLE_TIGHT_INCLUSION"
                if (SCCD_USE_TI) {
                    if (find_root_tight_inclusion_vf<T_HP>(
                            SCCD_MAX_DEPTH, SCCD_TOL, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
                        if (toi_stride == 0) {
                            atomic_min<T>(min_t, t);
                        } else {
                            toi[i] = t;
                        }
                    }
                    continue;
                }
#endif
                bool found = false;
                if (SCCD_ADAPTIVE_SPLIT) {
                    found = find_root_grid_adaptive_split_vf<T_HP>(
                        SCCD_MAX_DEPTH, SCCD_TOL, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack, SCCD_REFINE);
                } else {
                    found = find_root_grid_uniform_split_vf<T_HP>(
                        SCCD_MAX_DEPTH, SCCD_TOL, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v, stack, SCCD_REFINE);
                }

                if (found) {
                    if (toi_stride == 0) {
                        atomic_min<T>(min_t, t);
                    } else {
                        toi[i] = t;
                    }
                }
            }
        });

        if (toi_stride == 0) toi[0] = min_t;
        return 0;
    }

    template <typename T, typename I>
    int narrow_phase_newton_pass_ee(const size_t noverlaps,
                                    const I* const SCCD_RESTRICT e0overalp,
                                    const I* const SCCD_RESTRICT e1overalp,
                                    T** const SCCD_RESTRICT v0,
                                    T** const SCCD_RESTRICT v1,
                                    const size_t edge_stride,
                                    I** const SCCD_RESTRICT edges,
                                    const T max_toi,
                                    T* const SCCD_RESTRICT toi,
                                    const int toi_stride = 0);

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
                        const int toi_stride = 0) {
        using T_HP = double;
        const T infty = 100000;

        if (toi_stride == 0 && max_toi == 0) {
            toi[0] = 0;
            printf("max_toi is 0 returning\n");
            return 0;
        }

        assert(toi_stride == 0 || toi_stride == 1);
        if (noverlaps == 0) {
            if (toi != nullptr && toi_stride == 0) toi[0] = max_toi;
            return 0;
        }
        assert(toi != nullptr);

        int SCCD_USE_TI = 0;
        SCCD_READ_ENV(SCCD_USE_TI, atoi);

        int SCCD_MAX_DEPTH = 32;
        SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);

        T_HP SCCD_TOL = std::is_same_v<float, T_HP> ? T(1e-8) : T(1e-14);
        SCCD_READ_ENV(SCCD_TOL, atof);

        int SCCD_REFINE = 0;
        SCCD_READ_ENV(SCCD_REFINE, atoi);

        int SCCD_ADAPTIVE_SPLIT = 0;
        SCCD_READ_ENV(SCCD_ADAPTIVE_SPLIT, atoi);

        std::atomic<T> min_t = max_toi;

        // int SCCD_USE_NEWTON_PASS = 0;
        // SCCD_READ_ENV(SCCD_USE_NEWTON_PASS, atoi);
        // if (SCCD_USE_NEWTON_PASS) {
        //     narrow_phase_newton_pass_ee<T, I>(
        //         noverlaps, e0overalp, e1overalp, v0, v1, edge_stride, edges, max_toi, toi, toi_stride);

        //     if (toi_stride == 0) {
        //         min_t = sccd::min<T>(min_t, T(toi[0]));
        //     }
        // }

        // if (toi_stride == 0 && min_t == 0) {
        //     toi[0] = 0;
        //     printf("[after newton pass] min_t is 0 returning\n");
        //     return 0;
        // }

        if (toi_stride == 0) toi[0] = max_toi;
        sccd::parallel_for_br(0, noverlaps, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            std::vector<Box<T_HP>> stack;

            for (ptrdiff_t i = rbegin; i < rend; i++) {
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
                T_HP t = min_t;
                T_HP u = 0;
                T_HP v = 0;

#ifdef SCCD_ENABLE_TIGHT_INCLUSION
#warning "SCCD_ENABLE_TIGHT_INCLUSION"
                if (SCCD_USE_TI) {
                    if (find_root_tight_inclusion_ee<T_HP>(
                            SCCD_MAX_DEPTH, SCCD_TOL, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v)) {
                        if (toi_stride == 0) {
                            atomic_min<T>(min_t, t);
                        } else {
                            toi[i] = t;
                        }
                    }
                    continue;
                }
#endif
                bool found = false;
                if (SCCD_ADAPTIVE_SPLIT) {
                    found = find_root_grid_adaptive_split_ee<T_HP>(
                        SCCD_MAX_DEPTH, SCCD_TOL, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v, stack, SCCD_REFINE);
                } else {
                    found = find_root_grid_uniform_split_ee<T_HP>(
                        SCCD_MAX_DEPTH, SCCD_TOL, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v, stack, SCCD_REFINE);
                }

                if (found) {
                    if (toi_stride == 0) {
                        atomic_min<T>(min_t, t);
                    } else {
                        toi[i] = t;
                    }
                }
            }
        });

        if (toi_stride == 0) toi[0] = min_t;
        return 0;
    }

    template <typename T>
    inline T clamp_unit(const T x) {
        return sccd::min<T>(T(1), sccd::max<T>(T(0), x));
    }

    template <typename T>
    inline bool find_root_newton_ee(const int max_iter,
                                    const T atol,
                                    const T s1[3],
                                    const T s2[3],
                                    const T s3[3],
                                    const T s4[3],
                                    const T e1[3],
                                    const T e2[3],
                                    const T e3[3],
                                    const T e4[3],
                                    T& t,
                                    T& u,
                                    T& v) {
        t = clamp_unit<T>(t);
        u = clamp_unit<T>(u);
        v = clamp_unit<T>(v);

        T f = 0;
        ee_objective<T>(s1, s1, s2, s3, s4, e1, e1, e2, e3, e4, t, u, v, &f);
        const T ftol = T(0.5) * atol * atol;
        if (f <= ftol) {
            return true;
        }

        for (int k = 0; k < max_iter; ++k) {
            T p[3] = {0, 0, 0};
            ee_objective_dir<T>(s1, s1, s2, s3, s4, e1, e1, e2, e3, e4, t, u, v, &f, p);
            if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
                break;
            }

            T best_t = t;
            T best_u = u;
            T best_v = v;
            T best_f = f;
            T alpha = 1;
            bool improved = false;

            for (int j = 0; j < 12; ++j) {
                const T cand_t = clamp_unit<T>(t - alpha * p[0]);
                const T cand_u = clamp_unit<T>(u - alpha * p[1]);
                const T cand_v = clamp_unit<T>(v - alpha * p[2]);

                T fnext = 0;
                ee_objective<T>(s1, s1, s2, s3, s4, e1, e1, e2, e3, e4, cand_t, cand_u, cand_v, &fnext);
                if (fnext < best_f) {
                    best_t = cand_t;
                    best_u = cand_u;
                    best_v = cand_v;
                    best_f = fnext;
                    improved = true;
                    break;
                }

                alpha *= T(0.5);
            }

            t = best_t;
            u = best_u;
            v = best_v;
            f = best_f;

            if (f <= ftol) {
                return true;
            }

            if (!improved) {
                break;
            }
        }

        return false;
    }

    template <typename T>
    inline bool find_root_newton_vf(const int max_iter,
                                    const T atol,
                                    const T sv[3],
                                    const T s1[3],
                                    const T s2[3],
                                    const T s3[3],
                                    const T ev[3],
                                    const T e1[3],
                                    const T e2[3],
                                    const T e3[3],
                                    T& t,
                                    T& u,
                                    T& v) {
        // t = clamp_unit<T>(t);
        // project_uv_simplex<T>(u, v);

        T s4[3] = {0, 0, 0};
        T e4[3] = {0, 0, 0};

        T f = 0;
        // vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f);
        const T ftol = T(0.5) * atol * atol;
        // if (f <= ftol) {
        //     return true;
        // }

        for (int k = 0; k < max_iter; ++k) {
            T p[3] = {0, 0, 0};
            vf_objective_dir<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f, p);
            // if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2])) {
            //     break;
            // }

            t -= p[0];
            u -= p[1];
            v -= p[2];

            // T best_t = t;
            // T best_u = u;
            // T best_v = v;
            // T best_f = f;
            // T alpha = 1;
            // bool improved = false;

            // for (int j = 0; j < 2; ++j) {
            //     const T cand_t = clamp_unit<T>(t - alpha * p[0]);
            //     T cand_u = u - alpha * p[1];
            //     T cand_v = v - alpha * p[2];
            //     project_uv_simplex<T>(cand_u, cand_v);

            //     T fnext = 0;
            //     vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, cand_t, cand_u, cand_v, &fnext);
            //     if (fnext < best_f) {
            //         best_t = cand_t;
            //         best_u = cand_u;
            //         best_v = cand_v;
            //         best_f = fnext;
            //         improved = true;
            //         break;
            //     }

            //     alpha *= T(0.5);
            // }

            // t = best_t;
            // u = best_u;
            // v = best_v;
            // f = best_f;

            // if (!improved) {
            //     break;
            // }
        }

        if (f <= ftol) {
            return true;
        }

        return false;
    }

    template <typename T>
    inline bool newton_epsilon_box_accept_ee(const T tol,
                                             const T s1[3],
                                             const T s2[3],
                                             const T s3[3],
                                             const T s4[3],
                                             const T e1[3],
                                             const T e2[3],
                                             const T e3[3],
                                             const T e4[3],
                                             const T t,
                                             const T u,
                                             const T v) {
        T tols[3];
        compute_edge_edge_tolerance<T>(tol, s1, s2, s3, s4, e1, e2, e3, e4, tols);

        const T radius = sccd::max<T>(tol, std::numeric_limits<T>::epsilon() * T(64));
        const T t0 = sccd::max<T>(T(0), t - radius);
        const T t1 = sccd::min<T>(T(1), t + radius);
        const T u0 = sccd::max<T>(T(0), u - radius);
        const T u1 = sccd::min<T>(T(1), u + radius);
        const T v0 = sccd::max<T>(T(0), v - radius);
        const T v1 = sccd::min<T>(T(1), v + radius);

        T fmin[3];
        T fmax[3];
        init_codomain_bounds<T>(fmin, fmax);

        for (int mask = 0; mask < 8; ++mask) {
            const T ct = (mask & 1) ? t1 : t0;
            const T cu = (mask & 2) ? u1 : u0;
            const T cv = (mask & 4) ? v1 : v0;
            T F[3];
            diff_ee<T>(s1, s2, s3, s4, e1, e2, e3, e4, ct, cu, cv, F);
            update_codomain_bounds<T>(F, fmin, fmax);
        }

        bool accepted = false;
        return codomain_acceptance<T>(fmin, fmax, tol, tols, accepted) && accepted;
    }

    template <typename T>
    inline bool newton_epsilon_box_accept_vf(const T tol,
                                             const T sv[3],
                                             const T s1[3],
                                             const T s2[3],
                                             const T s3[3],
                                             const T ev[3],
                                             const T e1[3],
                                             const T e2[3],
                                             const T e3[3],
                                             const T t,
                                             const T u,
                                             const T v) {
        T tols[3];
        compute_face_vertex_tolerance<T>(tol, sv, s1, s2, s3, ev, e1, e2, e3, tols);

        const T radius = sccd::max<T>(tol, std::numeric_limits<T>::epsilon() * T(64));
        const T t0 = sccd::max<T>(T(0), t - radius);
        const T t1 = sccd::min<T>(T(1), t + radius);
        const T u0 = sccd::max<T>(T(0), u - radius);
        const T u1 = sccd::min<T>(T(1), u + radius);
        const T v0 = sccd::max<T>(T(0), v - radius);
        const T v1 = sccd::min<T>(T(1), v + radius);

        T fmin[3];
        T fmax[3];
        init_codomain_bounds<T>(fmin, fmax);

        for (int mask = 0; mask < 8; ++mask) {
            const T ct = (mask & 1) ? t1 : t0;
            const T cu = (mask & 2) ? u1 : u0;
            const T cv = (mask & 4) ? v1 : v0;
            T F[3];
            diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, ct, cu, cv, F);
            update_codomain_bounds<T>(F, fmin, fmax);
        }

        bool accepted = false;
        return codomain_acceptance<T>(fmin, fmax, tol, tols, accepted) && accepted;
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
                                    const int toi_stride) {
        using T_HP = double;

        assert(toi_stride == 0 || toi_stride == 1);
        if (noverlaps == 0) {
            if (toi != nullptr && toi_stride == 0) toi[0] = max_toi;
            return 0;
        }
        assert(toi != nullptr);

        int SCCD_MAX_DEPTH = 32;
        SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);

        T_HP SCCD_TOL = std::is_same_v<float, T_HP> ? T(1e-8) : T(1e-14);
        SCCD_READ_ENV(SCCD_TOL, atof);

        std::atomic<T> min_t = max_toi;
        if (toi_stride == 0) toi[0] = max_toi;

        ptrdiff_t accepted = 0;

        sccd::parallel_for_br(0, noverlaps, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t i = rbegin; i < rend; ++i) {
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

                T_HP t = T_HP(0.5) * sccd::min<T_HP>(T_HP(max_toi), T_HP(1));
                // T_HP t = 0;
                T_HP u = T_HP(1.0 / 3.0);
                T_HP v = T_HP(1.0 / 3.0);

                if (!find_root_newton_vf<T_HP>(SCCD_MAX_DEPTH, SCCD_TOL, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
                    continue;
                }

                if (t < T_HP(0) || t > T_HP(max_toi) || u < T_HP(0) || v < T_HP(0) || u + v > T_HP(1)) {
                    continue;
                }

                if (!newton_epsilon_box_accept_vf<T_HP>(SCCD_TOL, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
                    continue;
                }

#pragma omp atomic update
                accepted++;

                if (toi_stride == 0) {
                    atomic_min<T>(min_t, T(t));
                } else {
                    toi[i] = T(t);
                }
            }
        });

        if (toi_stride == 0) toi[0] = min_t;

        if (toi_stride == 0) printf("[vf] min_t: %g (max_toi=%g, accepted=%ld)\n", T(min_t), T(max_toi), accepted);
        return 0;
    }

    template <typename T, typename I>
    int narrow_phase_newton_pass_ee(const size_t noverlaps,
                                    const I* const SCCD_RESTRICT e0overalp,
                                    const I* const SCCD_RESTRICT e1overalp,
                                    T** const SCCD_RESTRICT v0,
                                    T** const SCCD_RESTRICT v1,
                                    const size_t edge_stride,
                                    I** const SCCD_RESTRICT edges,
                                    const T max_toi,
                                    T* const SCCD_RESTRICT toi,
                                    const int toi_stride) {
        using T_HP = double;

        assert(toi_stride == 0 || toi_stride == 1);
        if (noverlaps == 0) {
            if (toi != nullptr && toi_stride == 0) toi[0] = max_toi;
            return 0;
        }
        assert(toi != nullptr);

        int SCCD_MAX_DEPTH = 32;
        SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);

        T_HP SCCD_TOL = std::is_same_v<float, T_HP> ? T(1e-8) : T(1e-14);
        SCCD_READ_ENV(SCCD_TOL, atof);

        std::atomic<T> min_t = max_toi;
        if (toi_stride == 0) toi[0] = max_toi;

        ptrdiff_t accepted = 0;

        sccd::parallel_for_br(0, noverlaps, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t i = rbegin; i < rend; ++i) {
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

                T_HP t = T_HP(0.5) * sccd::min<T_HP>(T_HP(max_toi), T_HP(1));
                T_HP u = T_HP(0.5);
                T_HP v = T_HP(0.5);

                if (!find_root_newton_ee<T_HP>(SCCD_MAX_DEPTH, SCCD_TOL, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v)) {
                    continue;
                }

                if (t < T_HP(0) || t > T_HP(max_toi) || u < T_HP(0) || u > T_HP(1) || v < T_HP(0) || v > T_HP(1)) {
                    continue;
                }

                if (!newton_epsilon_box_accept_ee<T_HP>(SCCD_TOL, s1, s2, s3, s4, e1, e2, e3, e4, t, u, v)) {
                    continue;
                }

                // #ifdef SCCD_ENABLE_TIGHT_INCLUSION
                //                 T_HP ti_t = t;
                //                 T_HP ti_u = u;
                //                 T_HP ti_v = v;
                //                 if (!find_root_tight_inclusion_ee<T_HP>(
                //                         SCCD_MAX_DEPTH, SCCD_TOL, s1, s2, s3, s4, e1, e2, e3, e4, ti_t, ti_u, ti_v))
                //                         {
                //                     continue;
                //                 }
                //                 if (ti_t < T_HP(0) || ti_t > T_HP(max_toi)) {
                //                     continue;
                //                 }
                //                 t = ti_t;
                // #endif

#pragma omp atomic update
                accepted++;

                if (toi_stride == 0) {
                    atomic_min<T>(min_t, T(t));
                } else {
                    toi[i] = T(t);
                }
            }
        });

        if (toi_stride == 0) toi[0] = min_t;
        if (toi_stride == 0) printf("[ee] min_t: %g (max_toi=%g, accepted=%ld)\n", T(min_t), T(max_toi), accepted);
        return 0;
    }

}  // namespace sccd

#endif
