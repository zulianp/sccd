#ifndef SCCD_BROADPHASE_CUH
#define SCCD_BROADPHASE_CUH

#include "sccd_base.hpp"

#include <cstddef>

namespace sccd {
    namespace device {
        template <typename T>
        int choose_axis(const int dim, const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);
    }  // namespace device
}  // namespace sccd

#endif  // SCCD_BROADPHASE_CUH