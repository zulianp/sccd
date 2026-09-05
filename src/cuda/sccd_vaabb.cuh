#ifndef SCCD_VAABB_CUH
#define SCCD_VAABB_CUH

#include "sccd_base.hpp"

#include <cstddef>

namespace sccd {
    namespace device {

        template <typename idx_t, typename geom_t, typename aabb_t>
        void compute_aabbs(const int nxe,
                           const ptrdiff_t n_elements,
                           const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);

        template <typename idx_t, typename geom_t, typename aabb_t>
        void compute_aabbs(const int nxe,
                           const ptrdiff_t n_elements,
                           const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                           const BoxRounding rounding);

        template <typename geom_t, typename aabb_t>
        void compute_aabbs(                           const ptrdiff_t n_nodes,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);

        template <typename geom_t, typename aabb_t>
        void compute_aabbs(                           const ptrdiff_t n_nodes,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                           const BoxRounding rounding);

    }  // namespace device

}  // namespace sccd

#endif  // SCCD_VAABB_CUH
