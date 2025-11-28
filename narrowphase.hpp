#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include "vaabb.h"

#include "roots.hpp"

#include "tight_inclusion/interval_root_finder.hpp"
#include "tight_inclusion/ccd.hpp"
#include <Eigen/Dense>

// // #include "snumerr.hpp"
// // #include "snumtol.hpp"
// // #include "stuv.hpp"
// // #include <limits>

namespace sccd {

// // enum class CCDStepResult {
// //   Continue = 0,
// //   NoInclusion,
// //   DomainTolSatisfied,
// //   BoxInside,
// //   CodomainTolSatisfied,
// //   DegenerateSplit
// // };

// // template <typename T>
// // int ccd_ee(const int max_iter, const T codomain_tol, const T ms,
// //            const int use_ms, const int allow_zero_toi,
// //            // Boxes
// //            // Start
// //            const T v0_sx, const T v0_sy, const T v0_sz,
// //            // V1
// //            const T v1_sx, const T v1_sy, const T v1_sz,
// //            // V2
// //            const T v2_sx, const T v2_sy, const T v2_sz,
// //            // V3
// //            const T v3_sx, const T v3_sy, const T v3_sz,
// //            // End
// //            // V0
// //            const T v0_ex, const T v0_ey, const T v0_ez,
// //            // V1
// //            const T v1_ex, const T v1_ey, const T v1_ez,
// //            // V2
// //            const T v2_ex, const T v2_ey, const T v2_ez,
// //            // V3
// //            const T v3_ex, const T v3_ey, const T v3_ez,
// //            // Output
// //            T *const toi) {

// // 	// This can be vectorized?
// //   int box_in = 0;
// //   const int ok = sccd_origin_in_inclusion_ee_soa<T>(
// //       t_l, t_u, u_l, u_u, v_l, v_u, v0sx, v0sy, v0sz, v1sx, v1sy, v1sz,
// v2sx,
// //       v2sy, v2sz, v3sx, v3sy, v3sz, v0ex, v0ey, v0ez, v1ex, v1ey, v1ez,
// v2ex,
// //       v2ey, v2ez, v3ex, v3ey, v3ez, ms, ms, ms, errx, erry, errz,
// &true_tol,
// //       &box_in);

// //   if (!ok) {
// //     return CCDStepResult::NoInclusion;
// //   }
// // }

//         const DeviceMatrix<Scalar>& d_vertices_t0,
//         const DeviceMatrix<Scalar>& d_vertices_t1,
//         const DeviceMatrix<int>& d_edges,
//         const DeviceMatrix<int>& d_faces,
//         const std::shared_ptr<DeviceAABBs> d_vertex_boxes,
//         const std::shared_ptr<DeviceAABBs> d_edge_boxes,
//         const std::shared_ptr<DeviceAABBs> d_face_boxes,
//         const Scalar min_distance,
//         const int max_iterations,
//         const Scalar tolerance,
//         const bool allow_zero_toi,
// #ifdef SCALABLE_CCD_TOI_PER_QUERY
//         std::vector<std::tuple<int, int, Scalar>>& collisions,

bool barycentricTriangle3D(
    const ticcd::Vector3& A,
    const ticcd::Vector3& B,
    const ticcd::Vector3& C,
    const ticcd::Vector3& P,
    double &u,
    double &v)
{
    using std::abs;

    ticcd::Vector3 e1 = B - A;
    ticcd::Vector3 e2 = C - A;
    ticcd::Vector3 n = e1.cross(e2).eval();

    ticcd::Vector3 dir = P - A;
    double dist = n.dot(dir);
    if(dist*dist > 1e-5) {
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
bool find_root(
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
    T& t,
    T& u,
    T& v)
{

#if 1
    ticcd::Vector3 v_t0(sv[0], sv[1], sv[2]);
    ticcd::Vector3 f0_t0(s1[0], s1[1], s1[2]);
    ticcd::Vector3 f1_t0(s2[0], s2[1], s2[2]);
    ticcd::Vector3 f2_t0(s3[0], s3[1], s3[2]);

    ticcd::Vector3 v_t1(ev[0], ev[1], ev[2]);
    ticcd::Vector3 f0_t1(e1[0], e1[1], e1[2]);
    ticcd::Vector3 f1_t1(e2[0], e2[1], e2[2]);
    ticcd::Vector3 f2_t1(e3[0], e3[1], e3[2]);

    ticcd::Array3 tol(1e-10, 1e-10, 1e-10);
    ticcd::Array3 err(1e-10, 1e-10, 1e-10);

    ticcd::Scalar ms = 0;
    ticcd::Scalar max_time = 1;
    ticcd::Scalar max_itr = 10000;
    ticcd::Scalar toi = 0;
    ticcd::Scalar output_tolerance = 1e-10;
    
    bool test_ok = ticcd::vertexFaceCCD(
        v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1,
        err, ms, toi, 1e-8, 1, 1000, output_tolerance, false,
        ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    double u0 = -1, v0 = -1;
    double discrepancy = -1;
    if (test_ok) {
        auto f0 = f0_t0 * (1 - toi) + toi * (f0_t1);
        auto f1 = f1_t0 * (1 - toi) + toi * (f1_t1);
        auto f2 = f2_t0 * (1 - toi) + toi * (f2_t1);
        auto pt = (1 - toi) * v_t0 + toi * v_t1;

        
        const bool inplane = barycentricTriangle3D(f0.eval(), f1.eval(), f2.eval(), pt.eval(), u0, v0);
        assert(inplane);

   

        test_ok =
        (u0 >= -1e-8 && v0 >= -1e-8 && u0 + v0 <= 1 + 1e-8 && toi >= -1e-8 && toi <= 1 + 1e-8);

        auto pt_rec = (1 - u0 - v0) * f0  + u0 * f1 + v0 * f2;
        auto diff = pt_rec - pt;

        discrepancy = diff.dot(diff);
        // printf("(%g, %g, %g) %g\n", t, u, v, diff.dot(diff));

        // t = toi;
        // u = u0;
        // v = v0;
    }

#endif

    T f0 = 0;
    T p[3] = { 0, 0, 0 };

    vf_objective<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f0);

    T f = f0;
    T fnext = 0;
    bool outside = false;
    bool converged = false;
    int k = 0;
    T alpha = 1;
    T norm_diff = 0;
    for (; k < 10000; k++) {
        // vf_objective_dir<T>(
        //     sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &fnext, p);
        vf_gradient<T>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, p);
    
    //
        alpha = 2;
        for (int j = 0; j < 10 && fnext > f; j++) {
            alpha /= 2;
            p[0] *= alpha;
            p[1] *= alpha;
            p[2] *= alpha;

            vf_objective<T>(
                sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &fnext);
        }

        // if(fnext > f) {
          

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
        if (norm_diff < 1e-14) {
            converged = true;
            break;
        }

        f = fnext;

        // if ((k + 1) % 40 == 0) {
        //     outside =
        //         !(u >= -1e-4 && v >= -1e-4 && u < 1 + 1e-4 && v < 1 + 1e-4);
        //     if (outside) {
        //         break;
        //     }
        // }
    }

#if 1
    bool ok = converged
        && (u >= -1e-8 && v >= -1e-8 && u + v <= 1 + 1e-8 && t >= -1e-8 && t <= 1 + 1e-8);

    if (test_ok != ok) {
        printf(
            "RF(%d) vs Ours(%d) (%g, %g, %g)[%g] == (%g, %g, %g) (k=%d, diff=%g, converged=%d)\n", int(test_ok), int(ok),
            double(toi), double(u0), double(v0),  double(discrepancy), double(t), double(u), double(v), k, double(norm_diff), int(converged));
    }
#endif

    // Check if inside face
    return converged
        && (u >= 1e-14 && v >= 1e-14 && u + v <= 1 + 1e-14 && t >= 0 && t <= 1);
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

    // tbb::parallel_for(tbb::blocked_range<size_t>(0, noverlaps),
    //                   [&](const tbb::blocked_range<size_t> &r) {
    //                     for (size_t i = r.begin(); i < r.end(); i++)
    for (size_t i = 0; i < noverlaps; i++) {
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
        const double s4[3] = { 0, 0, 0 }; // FIXME: remove me

        const double e1[3] = { v1[0][nodes[0]], v1[1][nodes[0]],
                               v1[2][nodes[0]] };
        const double e2[3] = { v1[0][nodes[1]], v1[1][nodes[1]],
                               v1[2][nodes[1]] };
        const double e3[3] = { v1[0][nodes[2]], v1[1][nodes[2]],
                               v1[2][nodes[2]] };
        const double e4[3] = { 0, 0, 0 }; // FIXME: remove me

        // Iteration variables
        double t = 0;
        double u = 0;
        double v = 0;

        if (find_root<double>(
                sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v)) {
            toi[i] = t;
            min_t = sccd::min<T>(t, min_t);
        } else {
            toi[i] = infty;
        }
    }
    //
    // });

    return min_t;
}
} // namespace sccd

#endif