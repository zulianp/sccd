#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "vaabb.h"

#include "roots.hpp"

#include "tight_inclusion/ccd.hpp"
#include "tight_inclusion/interval_root_finder.hpp"
#include <Eigen/Dense>

namespace sccd {

bool barycentricTriangle3D(
    const ticcd::Vector3& A,
    const ticcd::Vector3& B,
    const ticcd::Vector3& C,
    const ticcd::Vector3& P,
    double& u,
    double& v)
{
    using std::abs;

    ticcd::Vector3 e1 = B - A;
    ticcd::Vector3 e2 = C - A;
    ticcd::Vector3 n = e1.cross(e2).eval();

    ticcd::Vector3 dir = P - A;
    double dist = n.dot(dir);
    if (dist * dist > 1e-5) {
        return false;
    }

    // Compute local coordinates u and v: (P - A) = u * (B - A) + v * (C - A)
    // Solve: dir = u * e1 + v * e2
    // Using dot product method (more numerically stable):
    // dir · e1 = u * (e1 · e1) + v * (e2 · e1)
    // dir · e2 = u * (e1 · e2) + v * (e2 · e2)
    double d00 = e1.dot(e1);
    double d01 = e1.dot(e2);
    double d11 = e2.dot(e2);
    double d20 = dir.dot(e1);
    double d21 = dir.dot(e2);

    double denom = d00 * d11 - d01 * d01;
    if (abs(denom) < 1e-10) {
        // Degenerate triangle
        return false;
    }

    u = (d11 * d20 - d01 * d21) / denom;
    v = (d00 * d21 - d01 * d20) / denom;

    return true;
}

bool isInsideTriangle(
    const ticcd::Vector3& lambda, ticcd::Scalar tol = ticcd::Scalar(1e-6))
{
    return (lambda.array() >= -tol).all()
        && (lambda.array() <= ticcd::Scalar(1) + tol).all()
        && std::abs(lambda.sum() - ticcd::Scalar(1)) <= ticcd::Scalar(1e-6);
}

template <typename T>
bool find_root_newton(
    const int max_iter,
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
    T& v)
{
    T t0_ = t;
    T u0_ = u;
    T v0_ = v;

    T f0 = 0;
    T p[3] = { 0, 0, 0 };

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
        vf_objective_dir<T>(
            sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &fnext, p);
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

            vf_objective<T>(
                sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, temp_t, temp_u, temp_v,
                &fnext);
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
            T fx = t0 * (o * s1[d] + u * s2[d] + v * s3[d])
                + t1 * (o * e1[d] + u * e2[d] + v * e3[d]);

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
            outside =
                !(u >= -1e-4 && v >= -1e-4 && u < 1 + 1e-4 && v < 1 + 1e-4);
            if (outside) {
                break;
            }
        }
    }

    // Check if inside face
    return converged
        && (u >= atol && v >= atol && u + v <= 1 + atol && t >= 0 && t <= 1);
}

template <typename T>
inline static void detect_zero_soa(
    const int n,
    const T* const SFEM_RESTRICT F000,
    const T* const SFEM_RESTRICT F001,
    const T* const SFEM_RESTRICT F010,
    const T* const SFEM_RESTRICT F011,
    const T* const SFEM_RESTRICT F100,
    const T* const SFEM_RESTRICT F101,
    const T* const SFEM_RESTRICT F110,
    const T* const SFEM_RESTRICT F111,
    const T tol,
    int* const SFEM_RESTRICT contains_zero)
{
    for (int k = 0; k < n; k++) {
        T fmin = sccd::min(
            sccd::min(sccd::min(F000[k], F001[k]), sccd::min(F010[k], F011[k])),
            sccd::min(
                sccd::min(F100[k], F101[k]), sccd::min(F110[k], F111[k])));

        T fmax = sccd::max(
            sccd::max(sccd::max(F000[k], F001[k]), sccd::max(F010[k], F011[k])),
            sccd::max(
                sccd::max(F100[k], F101[k]), sccd::max(F110[k], F111[k])));
        // 0 in [fmin, fmax] => fmin <= 0 <= fmax
        contains_zero[k] &= fmin <= tol & fmax >= -tol;
    }
}

template <typename T>
inline static void detect_zero(
    const int n_a,
    const int n_b,
    const int n_c,
    const int stride_a,
    const int stride_b, /*stride_c = 1*/
    T* const SFEM_RESTRICT F,
    const T tol,
    int* const SFEM_RESTRICT contains_zero)
{

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
                n_c, &F[i0], &F[i1], &F[i2], &F[i3], &F[i4], &F[i5], &F[i6],
                &F[i7], tol, &contains_zero[i0]);
        }
    }
}

template <typename T>
inline static void sample_Fvf(
    const int n_a,
    const int n_b,
    const int n_c, //
    const int stride_a,
    const int stride_b, /*stride_c = 1*/
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
    T* const SFEM_RESTRICT F)
{
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
inline T norm_diff_vf(
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
    T& v)
{

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
// #define ROOT_FINDING_CHUNK_SIZE                                                \
//     ((ROOT_FINDING_TN + 1) * (ROOT_FINDING_UN + 1) * (ROOT_FINDING_VN + 1))


#define ROOT_FINDING_CHUNK_SIZE 4096

template <typename T>
bool find_root(
    const int max_iter,
    const T tol,
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
    T& v)
{
    T Fx[ROOT_FINDING_CHUNK_SIZE];
    T Fy[ROOT_FINDING_CHUNK_SIZE];
    T Fz[ROOT_FINDING_CHUNK_SIZE];
    int contains_zero[ROOT_FINDING_CHUNK_SIZE];

    sample_Fvf<T>(
        2, 2, 2, // sizes
        4, 2,    // strides
        0, 0, 0, // start
        1, 1, 1, // step
        sv[0], s1[0], s2[0], s3[0], ev[0], e1[0], e2[0], e3[0], Fx);

    sample_Fvf<T>(
        2, 2, 2, // sizes
        4, 2,    // strides
        0, 0, 0, // start
        1, 1, 1, // step
        sv[1], s1[1], s2[1], s3[1], ev[1], e1[1], e2[1], e3[1], Fy);

    sample_Fvf<T>(
        2, 2, 2, // sizes
        4, 2,    // strides
        0, 0, 0, // start
        1, 1, 1, // step
        sv[2], s1[2], s2[2], s3[2], ev[2], e1[2], e2[2], e3[2], Fz);

    const T Fx_aabb[2] = { sccd::array_min<T>(8, Fx),
                           sccd::array_max<T>(8, Fx) };
    const T Fy_aabb[2] = { sccd::array_min<T>(8, Fy),
                           sccd::array_max<T>(8, Fy) };
    const T Fz_aabb[2] = { sccd::array_min<T>(8, Fz),
                           sccd::array_max<T>(8, Fz) };

    bool intersects = Fx_aabb[0] <= tol & Fx_aabb[1] >= -tol
        && Fy_aabb[0] <= tol & Fy_aabb[1] >= -tol
        && Fz_aabb[0] <= tol & Fz_aabb[1] >= -tol;

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

    int N[3];
    T chunk_size = ROOT_FINDING_CHUNK_SIZE;
    int total = 1;
    for(int i = 0; i < 3; i++) {
      int ii = args[i];
        N[ii] = sccd::max<int>(2, floor((ROOT_FINDING_CHUNK_SIZE / T(total)) * (aabb_inv[ii]/tot_inv)));
        total *= N[ii];
        tot_inv /= aabb_inv[ii]; 
        tot_inv = nextafter_up(tot_inv);
    }
 
    int t_n = N[0] - 1;
    int u_n = N[1] - 1;
    int v_n = N[2] - 1;

    // total = (N[0] + 1) * (N[1] + 1) * (N[2] + 1);
    
    // printf("t_n: %d, u_n: %d, v_n: %d (Total: %d <= %d)\n", t_n, u_n, v_n, total, ROOT_FINDING_CHUNK_SIZE);
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

    sample_Fvf<T>(
        t_n + 1, u_n + 1, v_n + 1, // sizes
        t_stride, u_stride,        // strides
        0, 0, 0,                   // start
        t_h, u_h, v_h,             // step
        sv[0], s1[0], s2[0], s3[0], ev[0], e1[0], e2[0], e3[0], Fx);

    sample_Fvf<T>(
        t_n + 1, u_n + 1, v_n + 1, // sizes
        t_stride, u_stride,        // strides
        0, 0, 0,                   // start
        t_h, u_h, v_h,             // step
        sv[1], s1[1], s2[1], s3[1], ev[1], e1[1], e2[1], e3[1], Fy);

    sample_Fvf<T>(
        t_n + 1, u_n + 1, v_n + 1, // sizes
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
        if(found) break;
        for (int u_i = 0; u_i < u_n; u_i++) {
            for (int v_i = 0; v_i < v_n; v_i++) {
                const int idx = t_i * t_stride + u_i * u_stride + v_i;
                if (contains_zero[idx]) {
                    T t_current = t_h * t_i + t_min;
                    T u_current = u_h * u_i + u_min;
                    T v_current = v_h * v_i + v_min;

                    if(u_h + v_h > 1 + tol) continue;
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

template <typename T>
bool find_root_oracle(
    const int max_iter,
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
    T& v)
{
    ticcd::Vector3 v_t0(sv[0], sv[1], sv[2]);
    ticcd::Vector3 f0_t0(s1[0], s1[1], s1[2]);
    ticcd::Vector3 f1_t0(s2[0], s2[1], s2[2]);
    ticcd::Vector3 f2_t0(s3[0], s3[1], s3[2]);

    ticcd::Vector3 v_t1(ev[0], ev[1], ev[2]);
    ticcd::Vector3 f0_t1(e1[0], e1[1], e1[2]);
    ticcd::Vector3 f1_t1(e2[0], e2[1], e2[2]);
    ticcd::Vector3 f2_t1(e3[0], e3[1], e3[2]);

    ticcd::Array3 tol(atol, atol, atol);
    ticcd::Array3 err(1e-10, 1e-10, 1e-10);

    ticcd::Scalar ms = 0;
    ticcd::Scalar max_time = 1;
    ticcd::Scalar max_itr = 10000;
    ticcd::Scalar toi = 0;
    ticcd::Scalar output_tolerance = 1e-6;

    bool test_ok = ticcd::vertexFaceCCD(
        v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, err, ms, toi,
        1e-6, 1, 1000, output_tolerance, false,
        ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    double u0 = -1, v0 = -1;
    double discrepancy = -1;
    if (test_ok) {
        auto f0 = f0_t0 * (1 - toi) + toi * (f0_t1);
        auto f1 = f1_t0 * (1 - toi) + toi * (f1_t1);
        auto f2 = f2_t0 * (1 - toi) + toi * (f2_t1);
        auto pt = (1 - toi) * v_t0 + toi * v_t1;

        const bool inplane = barycentricTriangle3D(
            f0.eval(), f1.eval(), f2.eval(), pt.eval(), u0, v0);
        assert(inplane);

        test_ok =
            (u0 >= -1e-8 && v0 >= -1e-8 && u0 + v0 <= 1 + 1e-8 && toi >= -1e-8
             && toi <= 1 + 1e-8);

        auto pt_rec = (1 - u0 - v0) * f0 + u0 * f1 + v0 * f2;
        auto diff = pt_rec - pt;

        discrepancy = diff.dot(diff);
        t = toi;
        u = u0;
        v = v0;
    }
    return test_ok;
}

template <int nxe, typename T>
T narrow_phase_vf(
    const size_t noverlaps,
    const idx_t* const SFEM_RESTRICT voveralp,
    const idx_t* const SFEM_RESTRICT foveralp,
    // Geometric data
    T** const SFEM_RESTRICT v0,
    T** const SFEM_RESTRICT v1,
    const size_t face_stride,
    idx_t** const SFEM_RESTRICT faces,
    // Output
    T* const SFEM_RESTRICT toi)
{
    const T infty = 100000;
    T min_t = infty;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, noverlaps),
                      [&](const tbb::blocked_range<size_t> &r) {
                        for (size_t i = r.begin(); i < r.end(); i++) {
    // for (size_t i = 0; i < noverlaps; i++) {
        const idx_t vi = voveralp[i];
        const idx_t fi = foveralp[i];

        idx_t nodes[3] = { faces[0][fi * face_stride],
                           faces[1][fi * face_stride],
                           faces[2][fi * face_stride] };

        const double sv[3] = { v0[0][vi], v0[1][vi], v0[2][vi] };
        const double ev[3] = { v1[0][vi], v1[1][vi], v1[2][vi] };

        const double s1[3] = { v0[0][nodes[0]], v0[1][nodes[0]],
                               v0[2][nodes[0]] };
        const double s2[3] = { v0[0][nodes[1]], v0[1][nodes[1]],
                               v0[2][nodes[1]] };
        const double s3[3] = { v0[0][nodes[2]], v0[1][nodes[2]],
                               v0[2][nodes[2]] };

        const double e1[3] = { v1[0][nodes[0]], v1[1][nodes[0]],
                               v1[2][nodes[0]] };
        const double e2[3] = { v1[0][nodes[1]], v1[1][nodes[1]],
                               v1[2][nodes[1]] };
        const double e3[3] = { v1[0][nodes[2]], v1[1][nodes[2]],
                               v1[2][nodes[2]] };

        // Iteration variables
        double t = 0;
        double u = 0;
        double v = 0;

#define TEST_ORACLE 0
#if TEST_ORACLE
        double t_oracle = 0;
        double u_oracle = 0;
        double v_oracle = 0;
        bool test_ok = find_root_oracle<double>(
            1000, 1e-10, sv, s1, s2, s3, ev, e1, e2, e3, t_oracle, u_oracle,
            v_oracle);
#endif
        if (find_root<double>(
                100, 1e-10, sv, s1, s2, s3, ev, e1, e2, e3, t, u, v)) {
            toi[i] = t;
            min_t = sccd::min<T>(t, min_t);
#if TEST_ORACLE
            if (!test_ok) {
                printf(
                    "RF(Ours: (%g, %g, %g))\n", double(t), double(u),
                    double(v));
                // assert(false);
            }
#endif

        } else {
            toi[i] = infty;
#if TEST_ORACLE
            if (test_ok) {
                printf(
                    "Missed root: (%g, %g, %g) == (%g, %g, %g)\n", t_oracle,
                    u_oracle, v_oracle, t, u, v);
                assert(false);
            }
#endif
        }
    }

    //
    });

    return min_t;
}

} // namespace sccd

#endif