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
// //       t_l, t_u, u_l, u_u, v_l, v_u, v0sx, v0sy, v0sz, v1sx, v1sy, v1sz, v2sx,
// //       v2sy, v2sz, v3sx, v3sy, v3sz, v0ex, v0ey, v0ez, v1ex, v1ey, v1ez, v2ex,
// //       v2ey, v2ez, v3ex, v3ey, v3ez, ms, ms, ms, errx, erry, errz, &true_tol,
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

	template<int nxe>
	geom_t narrow_phase_vf(
		const size_t noverlaps,
		const idx_t * const SFEM_RESTRICT voveralp,
		const idx_t * const SFEM_RESTRICT foveralp,
		// Geometric data
		geom_t **const SFEM_RESTRICT v0, 
		geom_t **const SFEM_RESTRICT v1,
		idx_t **const SFEM_RESTRICT faces,
		// Output
		geom_t * const SFEM_RESTRICT toi)
	{
		for(size_t i = 0; i < noverlaps; i++) {
			const idx_t vi = voveralp[i];
			const idx_t fi = foveralp[i];


		}
	}

} // namespace sccd

#endif 