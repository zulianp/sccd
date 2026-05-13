#ifndef SCCD_LOWER_BOUND_ALL_TO_ALL_CUH
#define SCCD_LOWER_BOUND_ALL_TO_ALL_CUH

#include "sccd_base.hpp"

namespace sccd {

    namespace host {
        template <typename T, typename I>
        void lower_bound_all_to_all(const ptrdiff_t count_search_keys,
                                    const T* const SCCD_RESTRICT sorted_search_keys,
                                    const ptrdiff_t count_sorted_keys,
                                    const T* const SCCD_RESTRICT sorted_keys,
                                    I* const SCCD_RESTRICT indices);
    }  // namespace host

    namespace device {

        template <typename T, typename I>
        void lower_bound_all_to_all(const ptrdiff_t count_search_keys,
                                    const T* const SCCD_RESTRICT sorted_search_keys,
                                    const ptrdiff_t count_sorted_keys,
                                    const T* const SCCD_RESTRICT sorted_keys,
                                    I* const SCCD_RESTRICT indices);

    }  // namespace device
}  // namespace sccd

#endif  // SCCD_LOWER_BOUND_ALL_TO_ALL_CUH
