#ifndef S_ROOT_FINDER_HPP
#define S_ROOT_FINDER_HPP


#include <type_traits>
#include <algorithm>

#include "vaabb.h"
#include "roots.hpp"

namespace sccd {

template <typename T>
bool find_root_newton(const int max_iter, const T atol, const T sv[3],
                      const T s1[3], const T s2[3], const T s3[3],
                      const T ev[3], const T e1[3], const T e2[3],
                      const T e3[3], T &t, T &u, T &v) {
  T t0_ = t;
  T u0_ = u;
  T v0_ = v;

  T f0 = 0;
  T p[3] = {0, 0, 0};

  T s4[3] = {0, 0, 0}; // FIXME
  T e4[3] = {0, 0, 0}; // FIXME

  vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f0);

  T f = f0;
  T fnext = 0;
  bool outside = false;
  bool converged = false;
  int k = 0;
  T alpha = 1;
  T norm_diff = 0;
  for (; k < max_iter; k++) {
    vf_objective_dir<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &fnext,
                        p);
    // vf_gradient<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, p);

    //
    alpha = 2;
    for (int j = 0; j < 10 && fnext > f; j++) {
      alpha /= 2;
      p[0] *= alpha;
      p[1] *= alpha;
      p[2] *= alpha;

      T temp_t = t - p[0];
      T temp_u = u - p[1];
      T temp_v = v - p[2];

      vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, temp_t, temp_u,
                      temp_v, &fnext);
    }

    t -= p[0];
    u -= p[1];
    v -= p[2];

    const T t0 = (1 - t);
    const T t1 = t;
    const T o = (1 - u - v);
    T diff[3];
    norm_diff = 0;
    for (int d = 0; d < 3; d++) {
      // Vertex
      T vx = t0 * sv[d] + t1 * ev[d];

      // Face
      T fx = t0 * (o * s1[d] + u * s2[d] + v * s3[d]) +
             t1 * (o * e1[d] + u * e2[d] + v * e3[d]);

      diff[d] = vx - fx;
      norm_diff += diff[d] * diff[d];
    }

    norm_diff = sqrt(norm_diff);
    if (norm_diff < atol) {
      converged = true;
      break;
    }

    f = fnext;

    if ((k + 1) % 10 == 0) {
      outside = !(u >= -1e-4 && v >= -1e-4 && u < 1 + 1e-4 && v < 1 + 1e-4);
      if (outside) {
        break;
      }
    }
  }

  // Check if inside face
  return converged &&
         (u >= atol && v >= atol && u + v <= 1 + atol && t >= 0 && t <= 1);
}

template <typename T>
inline static void detect_zero_soa(
    const int n, const T *const SFEM_RESTRICT F000,
    const T *const SFEM_RESTRICT F001, const T *const SFEM_RESTRICT F010,
    const T *const SFEM_RESTRICT F011, const T *const SFEM_RESTRICT F100,
    const T *const SFEM_RESTRICT F101, const T *const SFEM_RESTRICT F110,
    const T *const SFEM_RESTRICT F111, const T tol,
    int *const SFEM_RESTRICT contains_zero) {
  for (int k = 0; k < n; k++) {
    T fmin = sccd::min(
        sccd::min(sccd::min(F000[k], F001[k]), sccd::min(F010[k], F011[k])),
        sccd::min(sccd::min(F100[k], F101[k]), sccd::min(F110[k], F111[k])));

    T fmax = sccd::max(
        sccd::max(sccd::max(F000[k], F001[k]), sccd::max(F010[k], F011[k])),
        sccd::max(sccd::max(F100[k], F101[k]), sccd::max(F110[k], F111[k])));
    // 0 in [fmin, fmax] => fmin <= 0 <= fmax
    contains_zero[k] &= fmin <= tol & fmax >= -tol;
  }
}

template <typename T>
inline static void detect_zero(const int n_a, const int n_b, const int n_c,
                               const int stride_a,
                               const int stride_b, /*stride_c = 1*/
                               T *const SFEM_RESTRICT F, const T tol,
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

      detect_zero_soa(n_c, &F[i0], &F[i1], &F[i2], &F[i3], &F[i4], &F[i5],
                      &F[i6], &F[i7], tol, &contains_zero[i0]);
    }
  }
}

template <typename T>
inline static void
sample_Fvf(const int n_a, const int n_b,
           const int n_c,                          //
           const int stride_a, const int stride_b, /*stride_c = 1*/
           const T start_a, const T start_b, const T start_c, const T ha,
           const T hb, const T hc, const T sv, const T ev, const T s1,
           const T s2, const T s3, const T e1, const T e2, const T e3,
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

template <typename T>
inline T norm_diff_vf(const T sv[3], const T s1[3], const T s2[3],
                      const T s3[3], const T ev[3], const T e1[3],
                      const T e2[3], const T e3[3], T &t, T &u, T &v) {

  T t0 = (1 - t);
  T t1 = t;
  T o = (1 - u - v);
  T diff[3];
  for (int d = 0; d < 3; d++) {
    diff[d] = t0 * sv[d] + t1 * ev[d] - (o * s1[d] + u * s2[d] + v * s3[d]);
  }
  return sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
}

// #define ROOT_FINDING_TN 32
// #define ROOT_FINDING_UN 32
// #define ROOT_FINDING_VN 32
// #define ROOT_FINDING_CHUNK_SIZE \
//     ((ROOT_FINDING_TN + 1) * (ROOT_FINDING_UN + 1) * (ROOT_FINDING_VN + 1))

#define ROOT_FINDING_CHUNK_SIZE 4096

template <typename T>
bool find_root(const int max_iter, const T tol, const T sv[3], const T s1[3],
               const T s2[3], const T s3[3], const T ev[3], const T e1[3],
               const T e2[3], const T e3[3], T &t, T &u, T &v) {
  T Fx[ROOT_FINDING_CHUNK_SIZE];
  T Fy[ROOT_FINDING_CHUNK_SIZE];
  T Fz[ROOT_FINDING_CHUNK_SIZE];
  int contains_zero[ROOT_FINDING_CHUNK_SIZE];

  sample_Fvf<T>(2, 2, 2, // sizes
                4, 2,    // strides
                0, 0, 0, // start
                1, 1, 1, // step
                sv[0], s1[0], s2[0], s3[0], ev[0], e1[0], e2[0], e3[0], Fx);

  sample_Fvf<T>(2, 2, 2, // sizes
                4, 2,    // strides
                0, 0, 0, // start
                1, 1, 1, // step
                sv[1], s1[1], s2[1], s3[1], ev[1], e1[1], e2[1], e3[1], Fy);

  sample_Fvf<T>(2, 2, 2, // sizes
                4, 2,    // strides
                0, 0, 0, // start
                1, 1, 1, // step
                sv[2], s1[2], s2[2], s3[2], ev[2], e1[2], e2[2], e3[2], Fz);

  const T Fx_aabb[2] = {sccd::array_min<T>(8, Fx), sccd::array_max<T>(8, Fx)};
  const T Fy_aabb[2] = {sccd::array_min<T>(8, Fy), sccd::array_max<T>(8, Fy)};
  const T Fz_aabb[2] = {sccd::array_min<T>(8, Fz), sccd::array_max<T>(8, Fz)};

  bool intersects = Fx_aabb[0] <= tol & Fx_aabb[1] >= -tol &&
                    Fy_aabb[0] <= tol & Fy_aabb[1] >= -tol &&
                    Fz_aabb[0] <= tol & Fz_aabb[1] >= -tol;

  if (!intersects) {
    return false;
  }

  T x_inv =
      sccd::max<int>(2, 1. / sccd::max<T>(1e-5, (Fx_aabb[1] - Fx_aabb[0])));
  T y_inv =
      sccd::max<int>(2, 1. / sccd::max<T>(1e-5, (Fy_aabb[1] - Fy_aabb[0])));
  T z_inv =
      sccd::max<int>(2, 1. / sccd::max<T>(1e-5, (Fz_aabb[1] - Fz_aabb[0])));
  T aabb_inv[3] = {x_inv, y_inv, z_inv};
  T tot_inv = nextafter_up(x_inv * y_inv * z_inv + 1e-5);

  int args[3] = {0, 1, 2};
  std::sort(args, args + 3,
            [&](int a, int b) { return aabb_inv[a] < aabb_inv[b]; });

  int N[3];
  T chunk_size = ROOT_FINDING_CHUNK_SIZE;
  int total = 1;
  for (int i = 0; i < 3; i++) {
    int ii = args[i];
    N[ii] = sccd::max<int>(2, floor((ROOT_FINDING_CHUNK_SIZE / T(total)) *
                                    (aabb_inv[ii] / tot_inv)));
    total *= N[ii];
    tot_inv /= aabb_inv[ii];
    tot_inv = nextafter_up(tot_inv);
  }

  int t_n = N[0] - 1;
  int u_n = N[1] - 1;
  int v_n = N[2] - 1;

  // total = (N[0] + 1) * (N[1] + 1) * (N[2] + 1);

  // printf("t_n: %d, u_n: %d, v_n: %d (Total: %d <= %d)\n", t_n, u_n, v_n,
  // total, ROOT_FINDING_CHUNK_SIZE);
  assert(total <= ROOT_FINDING_CHUNK_SIZE);

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

  sample_Fvf<T>(t_n + 1, u_n + 1, v_n + 1, // sizes
                t_stride, u_stride,        // strides
                0, 0, 0,                   // start
                t_h, u_h, v_h,             // step
                sv[0], s1[0], s2[0], s3[0], ev[0], e1[0], e2[0], e3[0], Fx);

  sample_Fvf<T>(t_n + 1, u_n + 1, v_n + 1, // sizes
                t_stride, u_stride,        // strides
                0, 0, 0,                   // start
                t_h, u_h, v_h,             // step
                sv[1], s1[1], s2[1], s3[1], ev[1], e1[1], e2[1], e3[1], Fy);

  sample_Fvf<T>(t_n + 1, u_n + 1, v_n + 1, // sizes
                t_stride, u_stride,        // strides
                0, 0, 0,                   // start
                t_h, u_h, v_h,             // step
                sv[2], s1[2], s2[2], s3[2], ev[2], e1[2], e2[2], e3[2], Fz);

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
    if (found)
      break;
    for (int u_i = 0; u_i < u_n; u_i++) {
      for (int v_i = 0; v_i < v_n; v_i++) {
        const int idx = t_i * t_stride + u_i * u_stride + v_i;
        if (contains_zero[idx]) {
          T t_current = t_h * t_i + t_min;
          T u_current = u_h * u_i + u_min;
          T v_current = v_h * v_i + v_min;

          if (u_h + v_h > 1 + tol)
            continue;
          bool ok = true;

          // bool ok = find_root_newton<T>(
          //   max_iter, tol, sv, s1, s2, s3, ev, e1, e2, e3,
          //         t_current, u_current, v_current);

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

  // if (count)
  //     printf("Count: %d\n", count);
  return found;
}

} // namespace sccd

#endif // S_ROOT_FINDER_HPP
