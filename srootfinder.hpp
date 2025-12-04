#ifndef S_ROOT_FINDER_HPP
#define S_ROOT_FINDER_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "vaabb.h"

#include "roots.hpp"
#include "snumtol.hpp"

namespace sccd {

    template <typename T>
    inline void project_uv_simplex(T &u, T &v) {
        u = sccd::max<T>(u, 0);
        v = sccd::max<T>(v, 0);
        const T s = u + v;
        if (s <= static_cast<T>(1)) {
            return;
        }

        T u_proj = static_cast<T>(0.5) * (u - v + 1);
        u_proj = sccd::min<T>(static_cast<T>(1), sccd::max<T>(0, u_proj));
        v = static_cast<T>(1) - u_proj;
        u = u_proj;
    }

    template <typename T>
    inline void diff_vf(const T sv[3],
                        const T s1[3],
                        const T s2[3],
                        const T s3[3],
                        const T ev[3],
                        const T e1[3],
                        const T e2[3],
                        const T e3[3],
                        const T &t,
                        const T &u,
                        const T &v,
                        T *const SFEM_RESTRICT diff) {
        T t0 = (1 - t);
        T t1 = t;
        T o = (1 - u - v);
        for (int d = 0; d < 3; d++) {
            diff[d] = t0 * sv[d] + t1 * ev[d] - (o * s1[d] + u * s2[d] + v * s3[d]);
        }
    }

    template <typename T>
    inline T norm_diff_vf(const T sv[3],
                          const T s1[3],
                          const T s2[3],
                          const T s3[3],
                          const T ev[3],
                          const T e1[3],
                          const T e2[3],
                          const T e3[3],
                          T &t,
                          T &u,
                          T &v) {
        T t0 = (1 - t);
        T t1 = t;
        T o = (1 - u - v);
        T diff[3];
        for (int d = 0; d < 3; d++) {
            diff[d] = t0 * sv[d] + t1 * ev[d] - (o * s1[d] + u * s2[d] + v * s3[d]);
        }
        return sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
    }

    template <typename T>
    bool find_root_newton(const int max_iter,
                          const T atol,
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
                          T &v) {
        project_uv_simplex<T>(u, v);
        t = sccd::min<T>(static_cast<T>(1), sccd::max<T>(0, t));

        T s4[3] = {0, 0, 0};
        T e4[3] = {0, 0, 0};

        T f = 0;
        vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f);

        for (int k = 0; k < max_iter; k++) {
            T p[3] = {0, 0, 0};
            vf_objective_dir<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f, p);

            T best_t = t;
            T best_u = u;
            T best_v = v;
            T best_f = f;

            T alpha = 1;
            bool improved = false;
            for (int j = 0; j < 12; j++) {
                T cand_t = t - alpha * p[0];
                T cand_u = u - alpha * p[1];
                T cand_v = v - alpha * p[2];

                cand_t = sccd::min<T>(static_cast<T>(1), sccd::max<T>(0, cand_t));
                project_uv_simplex<T>(cand_u, cand_v);

                T fnext = 0;
                vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, cand_t, cand_u, cand_v, &fnext);

                if (fnext < best_f) {
                    best_t = cand_t;
                    best_u = cand_u;
                    best_v = cand_v;
                    best_f = fnext;
                    improved = true;
                    break;
                }

                alpha *= static_cast<T>(0.5);
            }

            t = best_t;
            u = best_u;
            v = best_v;
            f = best_f;

            const T norm_diff = norm_diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, t, u, v);
            if (norm_diff < atol) {
                return (u >= -atol && v >= -atol && u + v <= 1 + atol && t >= 0 && t <= 1);
            }

            if (!improved) {
                break;
            }
        }

        return false;
    }

    template <typename T>
    inline static void detect_zero_soa(const int n,
                                       const T *const SFEM_RESTRICT F000,
                                       const T *const SFEM_RESTRICT F001,
                                       const T *const SFEM_RESTRICT F010,
                                       const T *const SFEM_RESTRICT F011,
                                       const T *const SFEM_RESTRICT F100,
                                       const T *const SFEM_RESTRICT F101,
                                       const T *const SFEM_RESTRICT F110,
                                       const T *const SFEM_RESTRICT F111,
                                       const T tol,
                                       int *const SFEM_RESTRICT contains_zero) {
        for (int k = 0; k < n; k++) {
            T fmin = sccd::min(sccd::min(sccd::min(F000[k], F001[k]), sccd::min(F010[k], F011[k])),
                               sccd::min(sccd::min(F100[k], F101[k]), sccd::min(F110[k], F111[k])));

            T fmax = sccd::max(sccd::max(sccd::max(F000[k], F001[k]), sccd::max(F010[k], F011[k])),
                               sccd::max(sccd::max(F100[k], F101[k]), sccd::max(F110[k], F111[k])));
            // 0 in [fmin, fmax] => fmin <= 0 <= fmax
            contains_zero[k] &= fmin <= tol & fmax >= -tol;
        }
    }

    template <typename T>
    inline static void detect_zero(const int n_a,
                                   const int n_b,
                                   const int n_c,
                                   const int stride_a,
                                   const int stride_b, /*stride_c = 1*/
                                   T *const SFEM_RESTRICT F,
                                   const T tol,
                                   int *const SFEM_RESTRICT contains_zero) {
        for (int a = 0; a < n_a; a++) {
            for (int b = 0; b < n_b; b++) {
                const int i0 = a * n_b * n_c + b * n_c;
                const int i1 = a * n_b * n_c + b * n_c + 1;
                const int i2 = a * n_b * n_c + (b + 1) * n_c;
                const int i3 = a * n_b * n_c + (b + 1) * n_c + 1;
                const int i4 = (a + 1) * n_b * n_c + b * n_c;
                const int i5 = (a + 1) * n_b * n_c + b * n_c + 1;
                const int i6 = (a + 1) * n_b * n_c + (b + 1) * n_c;
                const int i7 = (a + 1) * n_b * n_c + (b + 1) * n_c + 1;

                detect_zero_soa(
                    n_c, &F[i0], &F[i1], &F[i2], &F[i3], &F[i4], &F[i5], &F[i6], &F[i7], tol, &contains_zero[i0]);
            }
        }
    }

    template <typename T>
    inline static void sample_Fvf(const int n_a,
                                  const int n_b,
                                  const int n_c,
                                  const int stride_a,
                                  const int stride_b,
                                  const T start_a,
                                  const T start_b,
                                  const T start_c,
                                  const T ha,
                                  const T hb,
                                  const T hc,
                                  const T sv,
                                  const T ev,
                                  const T s1,
                                  const T s2,
                                  const T s3,
                                  const T e1,
                                  const T e2,
                                  const T e3,
                                  T *const SFEM_RESTRICT F) {
        for (int a = 0; a < n_a; a++) {
            for (int b = 0; b < n_b; b++) {
                for (int c = 0; c < n_c; c++) {
                    const int idx = a * stride_a + b * stride_b + c;
                    const T t = start_a + a * ha;
                    const T u = start_b + b * hb;
                    const T v = start_c + c * hc;
                    const T t0 = (1 - t);
                    const T t1 = t;
                    const T o = (1 - u - v);
                    F[idx] = t0 * sv + t1 * ev - (o * s1 + u * s2 + v * s3);
                }
            }
        }
    }

#define ROOT_FINDING_CHUNK_SIZE 729

    template <typename T>
    bool find_root(const int max_iter,
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
                   T &v) {
        T Fx[ROOT_FINDING_CHUNK_SIZE];
        T Fy[ROOT_FINDING_CHUNK_SIZE];
        T Fz[ROOT_FINDING_CHUNK_SIZE];
        int contains_zero[ROOT_FINDING_CHUNK_SIZE];

        sample_Fvf<T>(2, 2, 2, 4, 2, 0, 0, 0, 1, 1, 1, sv[0], s1[0], s2[0], s3[0], ev[0], e1[0], e2[0], e3[0], Fx);
        sample_Fvf<T>(2, 2, 2, 4, 2, 0, 0, 0, 1, 1, 1, sv[1], s1[1], s2[1], s3[1], ev[1], e1[1], e2[1], e3[1], Fy);
        sample_Fvf<T>(2, 2, 2, 4, 2, 0, 0, 0, 1, 1, 1, sv[2], s1[2], s2[2], s3[2], ev[2], e1[2], e2[2], e3[2], Fz);

        const T Fx_aabb[2] = {sccd::array_min<T>(8, Fx), sccd::array_max<T>(8, Fx)};
        const T Fy_aabb[2] = {sccd::array_min<T>(8, Fy), sccd::array_max<T>(8, Fy)};
        const T Fz_aabb[2] = {sccd::array_min<T>(8, Fz), sccd::array_max<T>(8, Fz)};

        bool intersects = Fx_aabb[0] <= tol & Fx_aabb[1] >= -tol && Fy_aabb[0] <= tol & Fy_aabb[1] >= -tol &&
                          Fz_aabb[0] <= tol & Fz_aabb[1] >= -tol;

        if (!intersects) {
            return false;
        }

        T x_inv = sccd::max<int>(2, 1. / sccd::max<T>(1e-5, (Fx_aabb[1] - Fx_aabb[0])));
        T y_inv = sccd::max<int>(2, 1. / sccd::max<T>(1e-5, (Fy_aabb[1] - Fy_aabb[0])));
        T z_inv = sccd::max<int>(2, 1. / sccd::max<T>(1e-5, (Fz_aabb[1] - Fz_aabb[0])));
        T aabb_inv[3] = {x_inv, y_inv, z_inv};
        T tot_inv = nextafter_up(x_inv * y_inv * z_inv + 1e-5);

        int args[3] = {0, 1, 2};
        std::sort(args, args + 3, [&](int a, int b) { return aabb_inv[a] < aabb_inv[b]; });

        int t_n = 8;
        int u_n = 8;
        int v_n = 8;
        int total = (t_n + 1) * (u_n + 1) * (v_n + 1);

        const T t_min = 0;
        const T u_min = 0;
        const T v_min = 0;
        const T t_max = 1;
        const T u_max = 1;
        const T v_max = 1;
        const T t_h = (t_max - t_min) / t_n;
        const T u_h = (u_max - u_min) / u_n;
        const T v_h = (v_max - v_min) / v_n;

        int t_stride = (u_n + 1) * (v_n + 1);
        int u_stride = v_n + 1;

        sample_Fvf<T>(t_n + 1,
                      u_n + 1,
                      v_n + 1,
                      t_stride,
                      u_stride,
                      0,
                      0,
                      0,
                      t_h,
                      u_h,
                      v_h,
                      sv[0],
                      s1[0],
                      s2[0],
                      s3[0],
                      ev[0],
                      e1[0],
                      e2[0],
                      e3[0],
                      Fx);

        sample_Fvf<T>(t_n + 1,
                      u_n + 1,
                      v_n + 1,
                      t_stride,
                      u_stride,
                      0,
                      0,
                      0,
                      t_h,
                      u_h,
                      v_h,
                      sv[1],
                      s1[1],
                      s2[1],
                      s3[1],
                      ev[1],
                      e1[1],
                      e2[1],
                      e3[1],
                      Fy);

        sample_Fvf<T>(t_n + 1,
                      u_n + 1,
                      v_n + 1,
                      t_stride,
                      u_stride,
                      0,
                      0,
                      0,
                      t_h,
                      u_h,
                      v_h,
                      sv[2],
                      s1[2],
                      s2[2],
                      s3[2],
                      ev[2],
                      e1[2],
                      e2[2],
                      e3[2],
                      Fz);

        for (int i = 0; i < t_n * u_n * v_n; i++) {
            contains_zero[i] = 1;
        }

        detect_zero(t_n, u_n, v_n, t_stride, u_stride, Fx, tol, contains_zero);
        detect_zero(t_n, u_n, v_n, t_stride, u_stride, Fy, tol, contains_zero);
        detect_zero(t_n, u_n, v_n, t_stride, u_stride, Fz, tol, contains_zero);

        t = t_max + 1;
        u = u_max + 1;
        v = v_max + 1;
        bool found = false;
        int count = 0;
        for (int t_i = 0; t_i < t_n; t_i++) {
            if (found) break;
            for (int u_i = 0; u_i < u_n; u_i++) {
                for (int v_i = 0; v_i < v_n; v_i++) {
                    const int idx = t_i * t_stride + u_i * u_stride + v_i;
                    if (contains_zero[idx]) {
                        T t_current = t_h * t_i + t_min;
                        T u_current = u_h * u_i + u_min;
                        T v_current = v_h * v_i + v_min;

                        if (u_h + v_h > 1 + tol) continue;
                        bool ok = find_root_newton<T>(
                            max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, t_current, u_current, v_current);

                        count += ok;
                        if (ok && t_current < t) {
                            t = t_current;
                            u = u_current;
                            v = v_current;
                            found = true;
                        }
                    }
                }
            }
        }

        return found;
    }

    template <typename T>
    bool find_root_bisection(const int max_iter,
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
                             T &v) {
        struct Box {
            T tl, tu;
            T ul, uu;
            T vl, vu;
        };

        // Compute per-axis tolerances (matching snumtol.hpp signature)
        T tol0 = tol, tol1 = tol, tol2 = tol;
        compute_face_vertex_tolerance_soa<T>(tol,
                                             sv[0],
                                             sv[1],
                                             sv[2],
                                             s1[0],
                                             s1[1],
                                             s1[2],
                                             s2[0],
                                             s2[1],
                                             s2[2],
                                             s3[0],
                                             s3[1],
                                             s3[2],
                                             ev[0],
                                             ev[1],
                                             ev[2],
                                             e1[0],
                                             e1[1],
                                             e1[2],
                                             e2[0],
                                             e2[1],
                                             e2[2],
                                             e3[0],
                                             e3[1],
                                             e3[2],
                                             &tol0,
                                             &tol1,
                                             &tol2);

        printf("tol %f -> tol0: %f, tol1: %f, tol2: %f\n", tol, tol0, tol1, tol2);

        // Depth-first stack
        std::vector<Box> stack;
        stack.reserve(max_iter * 2 + 8);
        stack.push_back({T(0), T(1), T(0), T(1), T(0), T(1)});

        int iter = 0;
        while (!stack.empty() && iter < max_iter) {
            Box box = stack.back();
            stack.pop_back();
            ++iter;

            // Quick rectangle-simplex disjoint test for the (u,v) slice
            const bool disjoint_simplex = (box.ul >= T(1) && box.uu >= T(1)) || (box.vl >= T(1) && box.vu >= T(1)) ||
                                          ((box.ul + box.vl > T(1)) && (box.ul + box.vu > T(1)) &&
                                           (box.uu + box.vl > T(1)) && (box.uu + box.vu > T(1)));
            if (disjoint_simplex) {
                continue;
            }

            T fmin[3] = {std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
            T fmax[3] = {
                std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};

            auto inside_simplex = [](const T u, const T v) -> bool {
                return (u >= T(0) && v >= T(0) && u + v <= T(1));
            };
            auto inside_rect = [&](T u, T v) { return (u >= box.ul && u <= box.uu && v >= box.vl && v <= box.vu); };

            std::vector<std::pair<T, T>> uv_samples;
            uv_samples.reserve(12);
            auto add_uv = [&](T u, T v) {
                const T eps = T(1e-12);
                for (const auto &p : uv_samples) {
                    if (sccd::abs<T>(p.first - u) <= eps && sccd::abs<T>(p.second - v) <= eps) return;
                }
                uv_samples.emplace_back(u, v);
            };

            const T uvals[2] = {box.ul, box.uu};
            const T vvals[2] = {box.vl, box.vu};
            // Rectangle corners inside simplex
            for (int iu = 0; iu < 2; ++iu) {
                for (int iv = 0; iv < 2; ++iv) {
                    const T u = uvals[iu];
                    const T v = vvals[iv];
                    if (inside_simplex(u, v)) add_uv(u, v);
                }
            }
            // Intersections with u+v=1 on edges
            for (int iu = 0; iu < 2; ++iu) {
                const T u = uvals[iu];
                const T v_int = T(1) - u;
                if (v_int >= box.vl && v_int <= box.vu) add_uv(u, v_int);
            }
            for (int iv = 0; iv < 2; ++iv) {
                const T v = vvals[iv];
                const T u_int = T(1) - v;
                if (u_int >= box.ul && u_int <= box.uu) add_uv(u_int, v);
            }
            // Triangle vertices inside rectangle
            if (inside_rect(T(0), T(0))) add_uv(T(0), T(0));
            if (inside_rect(T(1), T(0))) add_uv(T(1), T(0));
            if (inside_rect(T(0), T(1))) add_uv(T(0), T(1));
            // Fallback: project rectangle corners
            if (uv_samples.empty()) {
                for (int iu = 0; iu < 2; ++iu) {
                    for (int iv = 0; iv < 2; ++iv) {
                        T u = uvals[iu];
                        T v = vvals[iv];
                        project_uv_simplex<T>(u, v);
                        add_uv(u, v);
                    }
                }
            }

            auto eval_at_t = [&](const T tt) {
                for (const auto &uv : uv_samples) {
                    const T u = uv.first;
                    const T v = uv.second;
                    T Fv[3];
                    diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, tt, u, v, Fv);
                    fmin[0] = sccd::min(fmin[0], Fv[0]);
                    fmin[1] = sccd::min(fmin[1], Fv[1]);
                    fmin[2] = sccd::min(fmin[2], Fv[2]);
                    fmax[0] = sccd::max(fmax[0], Fv[0]);
                    fmax[1] = sccd::max(fmax[1], Fv[1]);
                    fmax[2] = sccd::max(fmax[2], Fv[2]);
                }
            };
            eval_at_t(box.tl);
            eval_at_t(box.tu);

            const bool intersects = (fmin[0] <= tol0 && fmax[0] >= -tol0) && (fmin[1] <= tol1 && fmax[1] >= -tol1) &&
                                    (fmin[2] <= tol2 && fmax[2] >= -tol2);

            if (!intersects) {
                // printf("NO INTERSECTION: Iteration %d\n", iter);
                // printf("BOX: (%f, %f, %f, %f, %f, %f)\n", box.tl, box.tu, box.ul, box.uu, box.vl, box.vu);
                // printf("fmin: (%f, %f, %f), fmax: (%f, %f, %f)\n", fmin[0], fmin[1], fmin[2], fmax[0], fmax[1],
                // fmax[2]); printf("--------------------------------\n");
                continue;
            } else {
                printf("INTERSECTION: Iteration %d/%d\n", iter, max_iter);
                printf("BOX: (%f, %f, %f, %f, %f, %f)\n", box.tl, box.tu, box.ul, box.uu, box.vl, box.vu);
                printf(
                    "fmin: (%f, %f, %f), fmax: (%f, %f, %f)\n", fmin[0], fmin[1], fmin[2], fmax[0], fmax[1], fmax[2]);
                printf("--------------------------------\n");
            }

            const T dt = box.tu - box.tl;
            const T du = box.uu - box.ul;
            const T dv = box.vu - box.vl;

            const bool small_enough = dt <= 10 * tol0 && du <= 10 * tol1 && dv <= 10 * tol2;

            if (small_enough) {
                T tt = (box.tl + box.tu) * T(0.5);
                T uu = (box.ul + box.uu) * T(0.5);
                T vv = (box.vl + box.vu) * T(0.5);

                T t_candidate = tt;
                T u_candidate = uu;
                T v_candidate = vv;

                if (find_root_newton<T>(
                        max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3, t_candidate, u_candidate, v_candidate)) {
                    t = t_candidate;
                    u = u_candidate;
                    v = v_candidate;
                    return true;
                }

                // Only accept midpoint if residual is small
                {
                    T Fmid[3];
                    diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, t_candidate, u_candidate, v_candidate, Fmid);
                    const T fn = sccd::abs<T>(Fmid[0]) + sccd::abs<T>(Fmid[1]) + sccd::abs<T>(Fmid[2]);
                    if (fn <= (tol0 + tol1 + tol2)) {
                        t = t_candidate;
                        u = u_candidate;
                        v = v_candidate;
                        return true;
                    }
                }

                // continue to split
            }

            // Split along the widest dimension to reduce volume fast
            int split_dim = 0;
            if (du > dt && du >= dv) {
                split_dim = 1;
            } else if (dv > dt && dv > du) {
                split_dim = 2;
            }

            if (split_dim == 0) {
                const T mid = (box.tl + box.tu) * T(0.5);
                Box upper = box;
                upper.tl = mid;
                Box lower = box;
                lower.tu = mid;

                stack.push_back(upper);
                stack.push_back(lower);
            } else if (split_dim == 1) {
                const T mid = (box.ul + box.uu) * T(0.5);
                Box upper = box;
                upper.ul = mid;
                // if (upper.ul + upper.vl <= tol1 + tol2) {
                stack.push_back(upper);
                // }

                Box lower = box;
                lower.uu = mid;
                stack.push_back(lower);
            } else {
                const T mid = (box.vl + box.vu) * T(0.5);
                Box upper = box;
                upper.vl = mid;
                // if (upper.ul + upper.vl <= tol1 + tol2) {
                stack.push_back(upper);
                // }

                Box lower = box;
                lower.vu = mid;
                stack.push_back(lower);
            }
        }

        return false;
    }

    template <typename T>
    bool find_root_check(const int max_iter,
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
                         T &v) {
        T Fx_coarse[8];
        T Fy_coarse[8];
        T Fz_coarse[8];

        sample_Fvf<T>(
            2, 2, 2, 4, 2, 0, 0, 0, 1, 1, 1, sv[0], s1[0], s2[0], s3[0], ev[0], e1[0], e2[0], e3[0], Fx_coarse);
        sample_Fvf<T>(
            2, 2, 2, 4, 2, 0, 0, 0, 1, 1, 1, sv[1], s1[1], s2[1], s3[1], ev[1], e1[1], e2[1], e3[1], Fy_coarse);
        sample_Fvf<T>(
            2, 2, 2, 4, 2, 0, 0, 0, 1, 1, 1, sv[2], s1[2], s2[2], s3[2], ev[2], e1[2], e2[2], e3[2], Fz_coarse);

        const T Fx_min = sccd::array_min<T>(8, Fx_coarse);
        const T Fx_max = sccd::array_max<T>(8, Fx_coarse);
        const T Fy_min = sccd::array_min<T>(8, Fy_coarse);
        const T Fy_max = sccd::array_max<T>(8, Fy_coarse);
        const T Fz_min = sccd::array_min<T>(8, Fz_coarse);
        const T Fz_max = sccd::array_max<T>(8, Fz_coarse);

        if (!((Fx_min <= tol && Fx_max >= -tol) && (Fy_min <= tol && Fy_max >= -tol) &&
              (Fz_min <= tol && Fz_max >= -tol))) {
            return false;
        }

        auto safe_inv_range = [](const T vmin, const T vmax) -> T {
            const T rng = sccd::max<T>(static_cast<T>(1e-5), vmax - vmin);
            return sccd::max<T>(static_cast<T>(2), static_cast<T>(1) / rng);
        };

        const T invs[3] = {
            safe_inv_range(Fx_min, Fx_max), safe_inv_range(Fy_min, Fy_max), safe_inv_range(Fz_min, Fz_max)};
        int order[3] = {0, 1, 2};
        std::sort(order, order + 3, [&](int a, int b) { return invs[a] < invs[b]; });

        int N[3] = {0, 0, 0};
        int total = 1;
        T tot_inv = invs[0] * invs[1] * invs[2];
        const int base = ROOT_FINDING_CHUNK_SIZE;
        for (int oi = 0; oi < 3; ++oi) {
            const int idx = order[oi];
            const T ratio = invs[idx] / sccd::max<T>(static_cast<T>(1e-12), tot_inv);
            const T exponent = static_cast<T>(1.0) / static_cast<T>(3 - oi);
            const T target = (static_cast<T>(base) / static_cast<T>(total)) * std::pow(ratio, exponent);
            int n_val = static_cast<int>(std::floor(target));
            const int min_val = (idx == 0) ? 8 : 2;
            if (n_val < min_val) n_val = min_val;
            N[idx] = n_val;
            total *= N[idx];
            tot_inv = sccd::max<T>(static_cast<T>(1e-12), tot_inv / invs[idx]);
        }

        while (total > base) {
            int i = (N[0] >= N[1] && N[0] >= N[2]) ? 0 : ((N[1] >= N[2]) ? 1 : 2);
            if (N[i] > 2) {
                total /= N[i];
                --N[i];
                total *= N[i];
            } else {
                break;
            }
        }

        const int Nt = N[0];
        const int Nu = N[1];
        const int Nv = N[2];
        const int t_n = Nt - 1;
        const int u_n = Nu - 1;
        const int v_n = Nv - 1;
        if (t_n <= 0 || u_n <= 0 || v_n <= 0) {
            return false;
        }

        const T t_min = 0;
        const T u_min = 0;
        const T v_min = 0;
        const T t_max = 1;
        const T u_max = 1;
        const T v_max = 1;
        const T t_h = (t_max - t_min) / t_n;
        const T u_h = (u_max - u_min) / u_n;
        const T v_h = (v_max - v_min) / v_n;

        const int t_stride = (u_n + 1) * (v_n + 1);
        const int u_stride = (v_n + 1);

        std::vector<T> Fx((t_n + 1) * (u_n + 1) * (v_n + 1));
        std::vector<T> Fy((t_n + 1) * (u_n + 1) * (v_n + 1));
        std::vector<T> Fz((t_n + 1) * (u_n + 1) * (v_n + 1));

        sample_Fvf<T>(t_n + 1,
                      u_n + 1,
                      v_n + 1,
                      t_stride,
                      u_stride,
                      t_min,
                      u_min,
                      v_min,
                      t_h,
                      u_h,
                      v_h,
                      sv[0],
                      s1[0],
                      s2[0],
                      s3[0],
                      ev[0],
                      e1[0],
                      e2[0],
                      e3[0],
                      Fx.data());
        sample_Fvf<T>(t_n + 1,
                      u_n + 1,
                      v_n + 1,
                      t_stride,
                      u_stride,
                      t_min,
                      u_min,
                      v_min,
                      t_h,
                      u_h,
                      v_h,
                      sv[1],
                      s1[1],
                      s2[1],
                      s3[1],
                      ev[1],
                      e1[1],
                      e2[1],
                      e3[1],
                      Fy.data());
        sample_Fvf<T>(t_n + 1,
                      u_n + 1,
                      v_n + 1,
                      t_stride,
                      u_stride,
                      t_min,
                      u_min,
                      v_min,
                      t_h,
                      u_h,
                      v_h,
                      sv[2],
                      s1[2],
                      s2[2],
                      s3[2],
                      ev[2],
                      e1[2],
                      e2[2],
                      e3[2],
                      Fz.data());

        std::vector<int> contains_zero(t_n * u_n * v_n, 1);
        detect_zero(t_n, u_n, v_n, t_stride, u_stride, Fx.data(), tol, contains_zero.data());
        detect_zero(t_n, u_n, v_n, t_stride, u_stride, Fy.data(), tol, contains_zero.data());
        detect_zero(t_n, u_n, v_n, t_stride, u_stride, Fz.data(), tol, contains_zero.data());

        struct Seed {
            T prio;
            T fn;
            int ti, ui, vi;
            T t, u, v;
        };
        std::vector<Seed> seeds;
        seeds.reserve(t_n * u_n * v_n);

        for (int ti = 0; ti < t_n; ++ti) {
            for (int ui = 0; ui < u_n; ++ui) {
                for (int vi = 0; vi < v_n; ++vi) {
                    const int cell_idx = ti * u_n * v_n + ui * v_n + vi;
                    if (!contains_zero[cell_idx]) continue;
                    const T t_c = t_min + (static_cast<T>(ti) + static_cast<T>(0.5)) * t_h;
                    T u_c = u_min + (static_cast<T>(ui) + static_cast<T>(0.5)) * u_h;
                    T v_c = v_min + (static_cast<T>(vi) + static_cast<T>(0.5)) * v_h;

                    if ((u_c + v_c) > (static_cast<T>(1) + static_cast<T>(1e-8))) {
                        project_uv_simplex<T>(u_c, v_c);
                    }

                    T Fv_eval[3];
                    diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, t_c, u_c, v_c, Fv_eval);
                    const T fn_c =
                        std::sqrt(Fv_eval[0] * Fv_eval[0] + Fv_eval[1] * Fv_eval[1] + Fv_eval[2] * Fv_eval[2]);
                    seeds.push_back({t_c, fn_c, ti, ui, vi, t_c, u_c, v_c});
                }
            }
        }

        std::sort(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b) {
            if (a.prio == b.prio) return a.fn < b.fn;
            return a.prio < b.prio;
        });

        const size_t max_seeds = std::min<size_t>(1024, seeds.size());
        const T refine_tol = sccd::max<T>(static_cast<T>(1e-10), tol * static_cast<T>(1e-2));
        bool found_any = false;
        T best_t = std::numeric_limits<T>::max();
        T best_u = 0;
        T best_v = 0;

        for (size_t k = 0; k < max_seeds; ++k) {
            T t_current = seeds[k].t;
            T u_current = seeds[k].u;
            T v_current = seeds[k].v;
            const bool ok_ref = find_root_newton<T>(
                max_iter, refine_tol, sv, s1, s2, s3, ev, e1, e2, e3, t_current, u_current, v_current);
            if (ok_ref && (t_current >= 0) && (t_current <= 1) && (u_current >= static_cast<T>(-1e-8)) &&
                (v_current >= static_cast<T>(-1e-8)) &&
                (u_current + v_current <= static_cast<T>(1) + static_cast<T>(1e-8))) {
                if (!found_any || t_current < best_t) {
                    found_any = true;
                    best_t = t_current;
                    best_u = u_current;
                    best_v = v_current;
                }
            }
        }

        if (found_any) {
            t = best_t;
            u = best_u;
            v = best_v;
            return true;
        }

        auto dot3 = [](const T *a, const T *b) -> T { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
        auto closest_point_barycentric =
            [&](const T p[3], const T a[3], const T b[3], const T c[3], T &w0, T &w1, T &w2, T closest[3]) {
                T ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
                T ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
                T ap[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};

                const T d1 = dot3(ab, ap);
                const T d2 = dot3(ac, ap);
                if (d1 <= 0 && d2 <= 0) {
                    w0 = 1;
                    w1 = 0;
                    w2 = 0;
                    closest[0] = a[0];
                    closest[1] = a[1];
                    closest[2] = a[2];
                    return;
                }

                T bp[3] = {p[0] - b[0], p[1] - b[1], p[2] - b[2]};
                const T d3 = dot3(ab, bp);
                const T d4 = dot3(ac, bp);
                if (d3 >= 0 && d4 <= d3) {
                    w0 = 0;
                    w1 = 1;
                    w2 = 0;
                    closest[0] = b[0];
                    closest[1] = b[1];
                    closest[2] = b[2];
                    return;
                }

                const T vc = d1 * d4 - d3 * d2;
                if (vc <= 0 && d1 >= 0 && d3 <= 0) {
                    const T v_edge = d1 / (d1 - d3);
                    w0 = 1 - v_edge;
                    w1 = v_edge;
                    w2 = 0;
                    closest[0] = a[0] + v_edge * ab[0];
                    closest[1] = a[1] + v_edge * ab[1];
                    closest[2] = a[2] + v_edge * ab[2];
                    return;
                }

                T cp_vec[3] = {p[0] - c[0], p[1] - c[1], p[2] - c[2]};
                const T d5 = dot3(ab, cp_vec);
                const T d6 = dot3(ac, cp_vec);
                if (d6 >= 0 && d5 <= d6) {
                    w0 = 0;
                    w1 = 0;
                    w2 = 1;
                    closest[0] = c[0];
                    closest[1] = c[1];
                    closest[2] = c[2];
                    return;
                }

                const T vb = d5 * d2 - d1 * d6;
                if (vb <= 0 && d2 >= 0 && d6 <= 0) {
                    const T w_edge = d2 / (d2 - d6);
                    w0 = 1 - w_edge;
                    w1 = 0;
                    w2 = w_edge;
                    closest[0] = a[0] + w_edge * ac[0];
                    closest[1] = a[1] + w_edge * ac[1];
                    closest[2] = a[2] + w_edge * ac[2];
                    return;
                }

                const T va = d3 * d6 - d5 * d4;
                if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
                    const T w_edge = (d4 - d3) / ((d4 - d3) + (d5 - d6));
                    w0 = 0;
                    w1 = 1 - w_edge;
                    w2 = w_edge;
                    closest[0] = b[0] + w_edge * (c[0] - b[0]);
                    closest[1] = b[1] + w_edge * (c[1] - b[1]);
                    closest[2] = b[2] + w_edge * (c[2] - b[2]);
                    return;
                }

                const T denom = static_cast<T>(1) / (va + vb + vc);
                w1 = vb * denom;
                w2 = vc * denom;
                w0 = 1 - w1 - w2;
                closest[0] = w0 * a[0] + w1 * b[0] + w2 * c[0];
                closest[1] = w0 * a[1] + w1 * b[1] + w2 * c[1];
                closest[2] = w0 * a[2] + w1 * b[2] + w2 * c[2];
            };

        const int samples = 201;
        const T tol_g = sccd::max<T>(static_cast<T>(1e-6), tol * static_cast<T>(1e4));
        T best_fn = std::numeric_limits<T>::max();
        T best_t_sw = 0;
        T best_u_sw = 0;
        T best_v_sw = 0;

        for (int i = 0; i < samples; ++i) {
            const T t_s = static_cast<T>(i) / static_cast<T>(samples - 1);
            const T t0 = static_cast<T>(1) - t_s;
            const T t1 = t_s;
            T vp[3] = {t0 * sv[0] + t1 * ev[0], t0 * sv[1] + t1 * ev[1], t0 * sv[2] + t1 * ev[2]};
            T f1[3] = {t0 * s1[0] + t1 * e1[0], t0 * s1[1] + t1 * e1[1], t0 * s1[2] + t1 * e1[2]};
            T f2[3] = {t0 * s2[0] + t1 * e2[0], t0 * s2[1] + t1 * e2[1], t0 * s2[2] + t1 * e2[2]};
            T f3[3] = {t0 * s3[0] + t1 * e3[0], t0 * s3[1] + t1 * e3[1], t0 * s3[2] + t1 * e3[2]};

            T w0, w1, w2b;
            T closest[3];
            closest_point_barycentric(vp, f1, f2, f3, w0, w1, w2b, closest);
            const T rx = vp[0] - closest[0];
            const T ry = vp[1] - closest[1];
            const T rz = vp[2] - closest[2];
            const T fn = std::sqrt(rx * rx + ry * ry + rz * rz);
            if (fn < best_fn) {
                best_fn = fn;
                best_t_sw = t_s;
                best_u_sw = w1;
                best_v_sw = w2b;
            }
        }

        if (best_fn > tol_g) {
            return false;
        }

        {
            T t_candidate = best_t_sw;
            T u_candidate = best_u_sw;
            T v_candidate = best_v_sw;
            const bool ok_ref = find_root_newton<T>(
                max_iter, refine_tol, sv, s1, s2, s3, ev, e1, e2, e3, t_candidate, u_candidate, v_candidate);
            if (ok_ref && (t_candidate >= 0) && (t_candidate <= 1) && (u_candidate >= static_cast<T>(-1e-8)) &&
                (v_candidate >= static_cast<T>(-1e-8)) &&
                (u_candidate + v_candidate <= static_cast<T>(1) + static_cast<T>(1e-8))) {
                t = t_candidate;
                u = u_candidate;
                v = v_candidate;
                return true;
            }
        }

        return false;
    }


    template <typename T>
    bool find_root_dfs(const int max_iter,
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
                         T &v) {
        // Depth-first root search mirroring find_root_dfs_3D in numeric_roots.py
        auto eval_F = [&](const T tt, const T uu, const T vv, T &fx, T &fy, T &fz) {
            const T t0 = static_cast<T>(1) - tt;
            const T t1 = tt;
            const T o = static_cast<T>(1) - uu - vv;

            const T vx = t0 * sv[0] + t1 * ev[0];
            const T vy = t0 * sv[1] + t1 * ev[1];
            const T vz = t0 * sv[2] + t1 * ev[2];

            const T f1x = t0 * s1[0] + t1 * e1[0];
            const T f1y = t0 * s1[1] + t1 * e1[1];
            const T f1z = t0 * s1[2] + t1 * e1[2];
            const T f2x = t0 * s2[0] + t1 * e2[0];
            const T f2y = t0 * s2[1] + t1 * e2[1];
            const T f2z = t0 * s2[2] + t1 * e2[2];
            const T f3x = t0 * s3[0] + t1 * e3[0];
            const T f3y = t0 * s3[1] + t1 * e3[1];
            const T f3z = t0 * s3[2] + t1 * e3[2];

            const T fx_face = o * f1x + uu * f2x + vv * f3x;
            const T fy_face = o * f1y + uu * f2y + vv * f3y;
            const T fz_face = o * f1z + uu * f2z + vv * f3z;

            fx = vx - fx_face;
            fy = vy - fy_face;
            fz = vz - fz_face;
        };

        int side = std::max<int>(4, static_cast<int>(std::round(
                                       std::pow(static_cast<T>(ROOT_FINDING_CHUNK_SIZE), static_cast<T>(1.0 / 3.0)))));
        while ((side * side * side) > ROOT_FINDING_CHUNK_SIZE && side > 2) {
            --side;
        }

        const int t_n = side;
        const int u_n = side;
        const int v_n = side;

        const T t_min = static_cast<T>(0);
        const T u_min = static_cast<T>(0);
        const T v_min = static_cast<T>(0);
        const T t_max = static_cast<T>(1);
        const T u_max = static_cast<T>(1);
        const T v_max = static_cast<T>(1);

        const T t_h = (t_max - t_min) / static_cast<T>(t_n);
        const T u_h = (u_max - u_min) / static_cast<T>(u_n);
        const T v_h = (v_max - v_min) / static_cast<T>(v_n);

        const int t_stride = (u_n + 1) * (v_n + 1);
        const int u_stride = (v_n + 1);
        const int total_samples = (t_n + 1) * (u_n + 1) * (v_n + 1);

        std::vector<T> Fx(total_samples);
        std::vector<T> Fy(total_samples);
        std::vector<T> Fz(total_samples);

        auto refine_root = [&](T &tt, T &uu, T &vv) -> bool {
            auto eval = [&](const T t_eval, const T u_eval, const T v_eval, T out[3]) {
                eval_F(t_eval, u_eval, v_eval, out[0], out[1], out[2]);
            };

            tt = sccd::min<T>(static_cast<T>(1), sccd::max<T>(static_cast<T>(0), tt));
            project_uv_simplex<T>(uu, vv);

            const int refine_max_iter = 100;
            const T tol_f = sccd::max<T>(static_cast<T>(1e-10), tol * static_cast<T>(1e-2));

            T x[3] = {tt, uu, vv};
            for (int iter = 0; iter < refine_max_iter; ++iter) {
                T F0[3];
                eval(x[0], x[1], x[2], F0);
                const T fn = std::sqrt(F0[0] * F0[0] + F0[1] * F0[1] + F0[2] * F0[2]);
                if (fn <= tol_f) {
                    tt = x[0];
                    uu = x[1];
                    vv = x[2];
                    return true;
                }

                const T eps_t = static_cast<T>(1e-7);
                const T eps_u = static_cast<T>(1e-7);
                const T eps_v = static_cast<T>(1e-7);

                T Ft[3], Fu[3], Fv[3];
                eval(sccd::min<T>(static_cast<T>(1), x[0] + eps_t), x[1], x[2], Ft);
                eval(x[0], x[1] + eps_u, x[2], Fu);
                eval(x[0], x[1], x[2] + eps_v, Fv);

                T J[3][3];
                for (int r = 0; r < 3; ++r) {
                    J[r][0] = (Ft[r] - F0[r]) / eps_t;
                    J[r][1] = (Fu[r] - F0[r]) / eps_u;
                    J[r][2] = (Fv[r] - F0[r]) / eps_v;
                }

                // Normal equations: (J^T J) p = -J^T F0
                T JTJ[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                T JTF[3] = {0, 0, 0};
                for (int c = 0; c < 3; ++c) {
                    for (int r = 0; r < 3; ++r) {
                        JTJ[c][r] = J[0][c] * J[0][r] + J[1][c] * J[1][r] + J[2][c] * J[2][r];
                    }
                    JTF[c] = J[0][c] * F0[0] + J[1][c] * F0[1] + J[2][c] * F0[2];
                }

                const T lambda = static_cast<T>(1e-12);
                JTJ[0][0] += lambda;
                JTJ[1][1] += lambda;
                JTJ[2][2] += lambda;

                // Solve JTJ p = -JTF (Gaussian elimination with partial pivoting)
                T A[3][4] = {{JTJ[0][0], JTJ[0][1], JTJ[0][2], -JTF[0]},
                             {JTJ[1][0], JTJ[1][1], JTJ[1][2], -JTF[1]},
                             {JTJ[2][0], JTJ[2][1], JTJ[2][2], -JTF[2]}};

                for (int col = 0; col < 3; ++col) {
                    int pivot = col;
                    T max_abs = sccd::abs(A[col][col]);
                    for (int r = col + 1; r < 3; ++r) {
                        const T a = sccd::abs(A[r][col]);
                        if (a > max_abs) {
                            max_abs = a;
                            pivot = r;
                        }
                    }
                    if (max_abs < static_cast<T>(1e-20)) continue;
                    if (pivot != col) {
                        for (int c = col; c < 4; ++c) std::swap(A[col][c], A[pivot][c]);
                    }
                    const T inv_p = static_cast<T>(1) / A[col][col];
                    for (int c = col; c < 4; ++c) A[col][c] *= inv_p;
                    for (int r = 0; r < 3; ++r) {
                        if (r == col) continue;
                        const T factor = A[r][col];
                        if (factor == static_cast<T>(0)) continue;
                        for (int c = col; c < 4; ++c) {
                            A[r][c] -= factor * A[col][c];
                        }
                    }
                }

                T p[3] = {A[0][3], A[1][3], A[2][3]};

                T alpha = static_cast<T>(1);
                bool improved = false;
                for (int ls = 0; ls < 12; ++ls) {
                    T xn[3] = {x[0] + alpha * p[0], x[1] + alpha * p[1], x[2] + alpha * p[2]};
                    xn[0] = sccd::min<T>(static_cast<T>(1), sccd::max<T>(static_cast<T>(0), xn[0]));
                    project_uv_simplex<T>(xn[1], xn[2]);

                    T F1[3];
                    eval(xn[0], xn[1], xn[2], F1);
                    const T fn1 = std::sqrt(F1[0] * F1[0] + F1[1] * F1[1] + F1[2] * F1[2]);
                    if (fn1 < fn) {
                        x[0] = xn[0];
                        x[1] = xn[1];
                        x[2] = xn[2];
                        improved = true;
                        break;
                    }
                    alpha *= static_cast<T>(0.5);
                }

                if (!improved) {
                    break;
                }
            }

            T F_final[3];
            eval(x[0], x[1], x[2], F_final);
            const T fn_final = std::sqrt(F_final[0] * F_final[0] + F_final[1] * F_final[1] + F_final[2] * F_final[2]);
            tt = x[0];
            uu = x[1];
            vv = x[2];
            return fn_final <= sccd::max<T>(static_cast<T>(1e-10), tol * static_cast<T>(1e-2));
        };

        for (int ti = 0; ti <= t_n; ++ti) {
            const T tt = t_min + static_cast<T>(ti) * t_h;
            const int base_t = ti * t_stride;
            for (int ui = 0; ui <= u_n; ++ui) {
                const T uu = u_min + static_cast<T>(ui) * u_h;
                const int base_u = base_t + ui * u_stride;
                for (int vi = 0; vi <= v_n; ++vi) {
                    const T vv = v_min + static_cast<T>(vi) * v_h;
                    const int idx = base_u + vi;
                    T fx, fy, fz;
                    eval_F(tt, uu, vv, fx, fy, fz);
                    Fx[idx] = fx;
                    Fy[idx] = fy;
                    Fz[idx] = fz;
                }
            }
        }

        std::vector<int> contains_zero(t_n * u_n * v_n, 1);
        for (int ti = 0; ti < t_n; ++ti) {
            for (int ui = 0; ui < u_n; ++ui) {
                for (int vi = 0; vi < v_n; ++vi) {
                    const int cell_idx = ti * u_n * v_n + ui * v_n + vi;
                    const int i0 = ti * t_stride + ui * u_stride + vi;
                    const int i1 = i0 + 1;
                    const int i2 = i0 + u_stride;
                    const int i3 = i2 + 1;
                    const int i4 = i0 + t_stride;
                    const int i5 = i1 + t_stride;
                    const int i6 = i2 + t_stride;
                    const int i7 = i3 + t_stride;

                    const T xs[8] = {Fx[i0], Fx[i1], Fx[i2], Fx[i3], Fx[i4], Fx[i5], Fx[i6], Fx[i7]};
                    const T ys[8] = {Fy[i0], Fy[i1], Fy[i2], Fy[i3], Fy[i4], Fy[i5], Fy[i6], Fy[i7]};
                    const T zs[8] = {Fz[i0], Fz[i1], Fz[i2], Fz[i3], Fz[i4], Fz[i5], Fz[i6], Fz[i7]};

                    const T fx_min = sccd::array_min<T>(8, xs);
                    const T fx_max = sccd::array_max<T>(8, xs);
                    const T fy_min = sccd::array_min<T>(8, ys);
                    const T fy_max = sccd::array_max<T>(8, ys);
                    const T fz_min = sccd::array_min<T>(8, zs);
                    const T fz_max = sccd::array_max<T>(8, zs);

                    const bool has_zero_x = (fx_min <= tol) && (fx_max >= -tol);
                    const bool has_zero_y = (fy_min <= tol) && (fy_max >= -tol);
                    const bool has_zero_z = (fz_min <= tol) && (fz_max >= -tol);
                    contains_zero[cell_idx] = has_zero_x && has_zero_y && has_zero_z;
                }
            }
        }

        auto cell_has_zero = [&](const T t0, const T t1, const T u0, const T u1, const T v0, const T v1) -> bool {
            if ((sccd::max<T>(static_cast<T>(0), u0) + sccd::max<T>(static_cast<T>(0), v0)) >
                static_cast<T>(1) + static_cast<T>(1e-8)) {
                return false;
            }

            T fx_min = std::numeric_limits<T>::max();
            T fy_min = std::numeric_limits<T>::max();
            T fz_min = std::numeric_limits<T>::max();
            T fx_max = std::numeric_limits<T>::lowest();
            T fy_max = std::numeric_limits<T>::lowest();
            T fz_max = std::numeric_limits<T>::lowest();

            for (const T tt : {t0, t1}) {
                for (const T uu : {u0, u1}) {
                    for (const T vv : {v0, v1}) {
                        T fx, fy, fz;
                        eval_F(tt, uu, vv, fx, fy, fz);
                        fx_min = sccd::min<T>(fx_min, fx);
                        fy_min = sccd::min<T>(fy_min, fy);
                        fz_min = sccd::min<T>(fz_min, fz);
                        fx_max = sccd::max<T>(fx_max, fx);
                        fy_max = sccd::max<T>(fy_max, fy);
                        fz_max = sccd::max<T>(fz_max, fz);
                    }
                }
            }

            const bool has_x = (fx_min <= tol) && (fx_max >= -tol);
            const bool has_y = (fy_min <= tol) && (fy_max >= -tol);
            const bool has_z = (fz_min <= tol) && (fz_max >= -tol);
            return has_x && has_y && has_z;
        };

        struct Cell {
            T t0, t1, u0, u1, v0, v1;
            int depth;
        };

        std::vector<Cell> seeds;
        seeds.reserve(static_cast<size_t>(t_n * u_n * v_n));
        for (int ti = 0; ti < t_n; ++ti) {
            for (int ui = 0; ui < u_n; ++ui) {
                const T u0 = u_min + static_cast<T>(ui) * u_h;
                for (int vi = 0; vi < v_n; ++vi) {
                    const int cell_idx = ti * u_n * v_n + ui * v_n + vi;
                    if (!contains_zero[cell_idx]) continue;
                    const T v0 = v_min + static_cast<T>(vi) * v_h;
                    if ((u0 + v0) > static_cast<T>(1) + static_cast<T>(1e-8)) continue;

                    const T t0 = t_min + static_cast<T>(ti) * t_h;
                    seeds.push_back({t0, t0 + t_h, u0, u0 + u_h, v0, v0 + v_h, 0});
                }
            }
        }

        if (seeds.empty()) {
            return false;
        }

        std::sort(seeds.begin(), seeds.end(), [](const Cell &a, const Cell &b) { return a.t0 < b.t0; });
        std::vector<Cell> stack;
        stack.reserve(seeds.size());
        for (auto it = seeds.rbegin(); it != seeds.rend(); ++it) {
            stack.push_back(*it);
        }

        const T min_dim = sccd::max<T>(static_cast<T>(1e-6), tol * static_cast<T>(10));
        const int max_depth = 100;
        const T refine_tol = sccd::max<T>(static_cast<T>(1e-10), tol * static_cast<T>(1e-2));

        while (!stack.empty()) {
            Cell cell = stack.back();
            stack.pop_back();

            if ((sccd::max<T>(static_cast<T>(0), cell.u0) + sccd::max<T>(static_cast<T>(0), cell.v0)) >
                static_cast<T>(1) + static_cast<T>(1e-8)) {
                continue;
            }

            const T tc = static_cast<T>(0.5) * (cell.t0 + cell.t1);
            const T uc = static_cast<T>(0.5) * (cell.u0 + cell.u1);
            const T vc = static_cast<T>(0.5) * (cell.v0 + cell.v1);

                if (uc >= static_cast<T>(-1e-8) && vc >= static_cast<T>(-1e-8) &&
                    (uc + vc) <= static_cast<T>(1) + static_cast<T>(1e-8)) {
                    T fx_c, fy_c, fz_c;
                    eval_F(tc, uc, vc, fx_c, fy_c, fz_c);
                    const T fn_center = std::sqrt(fx_c * fx_c + fy_c * fy_c + fz_c * fz_c);
                    if (fn_center <= tol) {
                    T t_candidate = tc;
                    T u_candidate = uc;
                    T v_candidate = vc;
                    if (refine_root(t_candidate, u_candidate, v_candidate) ||
                        find_root_newton<T>(100, refine_tol, sv, s1, s2, s3, ev, e1, e2, e3, t_candidate, u_candidate, v_candidate)) {
                        t = t_candidate;
                        u = u_candidate;
                        v = v_candidate;
                        return true;
                    } else {
                        continue;
                    }
                }
            }

            const T size = sccd::max<T>(cell.t1 - cell.t0, sccd::max<T>(cell.u1 - cell.u0, cell.v1 - cell.v0));
            if (size <= min_dim || cell.depth >= max_iter) {
                continue;
            }

            const T tm = static_cast<T>(0.5) * (cell.t0 + cell.t1);
            const T um = static_cast<T>(0.5) * (cell.u0 + cell.u1);
            const T vm = static_cast<T>(0.5) * (cell.v0 + cell.v1);

            const Cell subcells[8] = {{cell.t0, tm, cell.u0, um, cell.v0, vm, cell.depth + 1},
                                      {cell.t0, tm, cell.u0, um, vm, cell.v1, cell.depth + 1},
                                      {cell.t0, tm, um, cell.u1, cell.v0, vm, cell.depth + 1},
                                      {cell.t0, tm, um, cell.u1, vm, cell.v1, cell.depth + 1},
                                      {tm, cell.t1, cell.u0, um, cell.v0, vm, cell.depth + 1},
                                      {tm, cell.t1, cell.u0, um, vm, cell.v1, cell.depth + 1},
                                      {tm, cell.t1, um, cell.u1, cell.v0, vm, cell.depth + 1},
                                      {tm, cell.t1, um, cell.u1, vm, cell.v1, cell.depth + 1}};

            if (cell.depth + 1 >= max_depth) {
                break;
            }

            for (const Cell &sub : subcells) {
                if (cell_has_zero(sub.t0, sub.t1, sub.u0, sub.u1, sub.v0, sub.v1)) {
                    stack.push_back(sub);
                }
            }
        }

        return false;
                         }

}  // namespace sccd

#endif  // S_ROOT_FINDER_HPP
