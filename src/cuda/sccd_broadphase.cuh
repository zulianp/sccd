#ifndef SCCD_BROADPHASE_CUH
#define SCCD_BROADPHASE_CUH

#include "sccd_base.hpp"

#include <cstddef>

namespace sccd {
    namespace device {
        template <typename T>
        int choose_axis(const int dim, const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);

        template <typename T>
        void enumerate(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx);

        template <typename T, typename I>
        void sort_along_axis(const int dim,
                             const ptrdiff_t n,
                             const int sort_axis,
                             T** const SCCD_RESTRICT arrays,
                             I* const SCCD_RESTRICT idx,
                             T* const SCCD_RESTRICT scratch);
    }  // namespace device
}  // namespace sccd

#endif  // SCCD_BROADPHASE_CUH