#ifndef SCCD_BROADPHASE_WARP_CUH
#define SCCD_BROADPHASE_WARP_CUH

#include "sccd_base.hpp"

#include <cstddef>

namespace sccd {
    namespace device {

        template <int first_nxe, int second_nxe, typename T, typename I>
        void count_overlaps_warp(const int sort_axis,
                                 const ptrdiff_t first_count,
                                 T** const SCCD_RESTRICT first_aabbs,
                                 I* const SCCD_RESTRICT first_idx,
                                 const ptrdiff_t first_element_stride,
                                 I** const SCCD_RESTRICT first_elements,
                                 const ptrdiff_t second_count,
                                 T** const SCCD_RESTRICT second_aabbs,
                                 I* const SCCD_RESTRICT second_idx,
                                 const ptrdiff_t second_element_stride,
                                 I** const SCCD_RESTRICT second_elements,
                                 ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                 const T* const SCCD_RESTRICT cummax);

    }  // namespace device
}  // namespace sccd

#endif  // SCCD_BROADPHASE_WARP_CUH
