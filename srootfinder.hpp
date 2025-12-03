#ifndef S_ROOT_FINDER_HPP
#define S_ROOT_FINDER_HPP

#include <algorithm>
#include <limits>
#include <type_traits>
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

            T fmin[3] = {std::numeric_limits<T>::max(), std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
            T fmax[3] = {
                std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};

            for (int corner = 0; corner < 8; corner++) {
                const T tt = (corner & 1) ? box.tu : box.tl;
                const T uu_ = (corner & 2) ? box.uu : box.ul;
                const T vv_ = (corner & 4) ? box.vu : box.vl;

                T F[3];
                diff_vf<T>(sv, s1, s2, s3, ev, e1, e2, e3, tt, uu_, vv_, F);

                fmin[0] = sccd::min(fmin[0], F[0]);
                fmin[1] = sccd::min(fmin[1], F[1]);
                fmin[2] = sccd::min(fmin[2], F[2]);

                fmax[0] = sccd::max(fmax[0], F[0]);
                fmax[1] = sccd::max(fmax[1], F[1]);
                fmax[2] = sccd::max(fmax[2], F[2]);
            }

            const bool intersects = (fmin[0] <= tol0 && fmax[0] >= -tol0) && (fmin[1] <= tol1 && fmax[1] >= -tol1) &&
                                    (fmin[2] <= tol2 && fmax[2] >= -tol2);

            if (!intersects) {
                // printf("NO INTERSECTION: Iteration %d\n", iter);
                // printf("BOX: (%f, %f, %f, %f, %f, %f)\n", box.tl, box.tu, box.ul, box.uu, box.vl, box.vu);
                // printf("fmin: (%f, %f, %f), fmax: (%f, %f, %f)\n", fmin[0], fmin[1], fmin[2], fmax[0], fmax[1], fmax[2]);
                // printf("--------------------------------\n");
                continue;
            } else {
              printf("INTERSECTION: Iteration %d/%d\n", iter, max_iter);
              printf("BOX: (%f, %f, %f, %f, %f, %f)\n", box.tl, box.tu, box.ul, box.uu, box.vl, box.vu);
              printf("fmin: (%f, %f, %f), fmax: (%f, %f, %f)\n", fmin[0], fmin[1], fmin[2], fmax[0], fmax[1], fmax[2]);
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

                t = t_candidate;
                u = u_candidate;
                v = v_candidate;

                if (t >= -tol0 && t <= 1 + tol0 && u >= -tol1 && v >= -tol2 && u + v <= 1 + tol1 + tol2) {
                    // t = t_candidate;
                    // u = u_candidate;
                    // v = v_candidate;
                    return true;
                }
                // continue;
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

}  // namespace sccd

#endif  // S_ROOT_FINDER_HPP
