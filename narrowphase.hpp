#ifndef NARROWPHASE_HPP
#define NARROWPHASE_HPP

#include "vaabb.h"

#include "roots.hpp"

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
    T infty = 100000;
    T min_t = infty;

    for (size_t i = 0; i < noverlaps; i++) {
        const idx_t vi = voveralp[i];
        const idx_t fi = foveralp[i];

        idx_t nodes[3] = { faces[0][fi * face_stride], faces[1][fi * face_stride], faces[2][fi * face_stride] };

        const double sv[3] = { v0[0][vi], v0[1][vi], v0[2][vi] };
        const double ev[3] = { v1[0][vi], v1[1][vi], v1[2][vi] };

        const double s1[3] = { v0[0][nodes[0]], v0[1][nodes[0]], v0[2][nodes[0]] };
        const double s2[3] = { v0[0][nodes[1]], v0[1][nodes[1]], v0[2][nodes[1]] };
        const double s3[3] = { v0[0][nodes[2]], v0[1][nodes[2]], v0[2][nodes[2]] };
        const double s4[3] = { 0, 0, 0 }; //FIXME: remove me

        const double e1[3] = { v1[0][nodes[0]], v1[1][nodes[0]], v1[2][nodes[0]] };
        const double e2[3] = { v1[0][nodes[1]], v1[1][nodes[1]], v1[2][nodes[1]] };
        const double e3[3] = { v1[0][nodes[2]], v1[1][nodes[2]], v1[2][nodes[2]] };
        const double e4[3] = { 0, 0, 0 };  //FIXME: remove me

        // Iteration variables
        double t = 0;
        double u = 1./3;
        double v = 1./3;
        
        double f0 = 0;
        double p[3] = { 0, 0, 0 };

        vf_objective_dir<double>(
            sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &f0, p);

        double f = f0;
        t -= p[0];
        u -= p[1];
        v -= p[2];
        double fnext = 0;
		bool outside = false;
		int k = 0; 
        for (;k < 600; k++) {
            vf_objective_dir<double>(
                sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, t, u, v, &fnext, p);

            double alpha = 1;
            for(int j = 0; j < 10 && fnext > f; j++) {
            	alpha /= 2;
            	p[0] *= alpha;
            	p[1] *= alpha;
            	p[2] *= alpha;

            	vf_objective<double>(
            	    sv,
            	    s1,
            	    s2,
            	    s3,
            	    s4,
            	    ev,
            	    e1,
            	    e2,
            	    e3,
            	    e4,
            	    t,
            	    u,
            	    v,
            	    &fnext);
            }

            t -= p[0];
            u -= p[1];
            v -= p[2];

            if (sccd::abs(fnext - f) < 1e-10) {
                break;
            }

            f = fnext;

            if((k + 1) % 10 == 0) {
            	outside = !(u >= -1e-4 && v >= -1e-4 && u < 1 + 1e-4 && v < 1 + 1e-4);
            	if(outside) break;
        	}
        }

        if(outside) {
            // printf("Not converged(%d)\n", k);
            // printf("t: %f, u: %f, v: %f\n", t, u, v);
            // printf("f: %f, fnext: %f\n", f0, f);
            // printf("p: %f, %f, %f\n", p[0], p[1], p[2]);

            // printf("sv: %f, %f, %f\n", sv[0], sv[1], sv[2]);
            // printf("ev: %f, %f, %f\n", ev[0], ev[1], ev[2]);

            // printf("s1: %f, %f, %f\n", s1[0], s1[1], s1[2]);
            // printf("s2: %f, %f, %f\n", s2[0], s2[1], s2[2]);
            // printf("s3: %f, %f, %f\n", s3[0], s3[1], s3[2]);
            // printf("e1: %f, %f, %f\n", e1[0], e1[1], e1[2]);
            // printf("e2: %f, %f, %f\n", e2[0], e2[1], e2[2]);
            // printf("e3: %f, %f, %f\n", e3[0], e3[1], e3[2]);
            // toi[i] = infty;
            // continue;
        }

        // Check if inside face
        if (u >= 1e-10 && v >= 1e-10 && sccd::abs<T>(u + v - 1) < 1e-10 && t > -1) {
            toi[i] = t;
            min_t = sccd::min<T>(t, min_t);
            printf("%f\n", t);
        } else {
            toi[i] = infty;
        }

        
    }

    return min_t;
}
} // namespace sccd

#endif