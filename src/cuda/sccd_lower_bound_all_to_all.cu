#include "sccd_lower_bound_all_to_all.cuh"

#include "sccd_cuda_base.cuh"

#include <algorithm>
#include <cassert>
#include <cstddef>

namespace sccd {
    namespace host {
        template <typename T, typename I>
        void lower_bound_all_to_all(const ptrdiff_t count_search_keys,
                                    const T* const SCCD_RESTRICT sorted_search_keys,
                                    const ptrdiff_t count_sorted_keys,
                                    const T* const SCCD_RESTRICT sorted_keys,
                                    I* const SCCD_RESTRICT indices) {
#pragma omp parallel for
            for (ptrdiff_t i = 0; i < count_search_keys; i++) {
                ptrdiff_t start = (i == 0) ? 0 : indices[i - 1];
                indices[i] =
                    std::lower_bound(sorted_keys, sorted_keys + count_sorted_keys, sorted_search_keys[i]) - sorted_keys;
            }
        }
    }  // namespace host

    namespace device {

        template <typename T>
        __device__ __forceinline__ const T* lower_bound(const T* const SCCD_RESTRICT begin,
                                                        const T* const SCCD_RESTRICT end,
                                                        const T val) {
            ptrdiff_t count = end - begin;
            const T* it = begin;

            while (count > 0) {
                const ptrdiff_t step = count >> 1;
                const T* const mid = it + step;

                if (*mid < val) {
                    it = mid + 1;
                    count -= step + 1;
                } else {
                    count = step;
                }
            }

            return it;
        }

        template <typename T, typename I>
        __global__ void lower_bound_all_to_all_baseline_kernel(const ptrdiff_t count_search_keys,
                                                               const T* const SCCD_RESTRICT sorted_search_keys,
                                                               const ptrdiff_t count_sorted_keys,
                                                               const T* const SCCD_RESTRICT sorted_keys,
                                                               I* const SCCD_RESTRICT indices) {
            const ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= count_search_keys) return;
            indices[i] = lower_bound(sorted_keys, sorted_keys + count_sorted_keys, sorted_search_keys[i]) - sorted_keys;
        }

        template <typename T, typename I>
        void lower_bound_all_to_all(const ptrdiff_t count_search_keys,
                                    const T* const SCCD_RESTRICT sorted_search_keys,
                                    const ptrdiff_t count_sorted_keys,
                                    const T* const SCCD_RESTRICT sorted_keys,
                                    I* const SCCD_RESTRICT indices) {
            int SCCD_LB_BASELINE = 1;
            SCCD_READ_ENV(SCCD_LB_BASELINE, atoi);

            if (SCCD_LB_BASELINE) {
                dim3 block(128, 1, 1);
                dim3 grid((count_search_keys + block.x - 1) / block.x, 1, 1);

                lower_bound_all_to_all_baseline_kernel<T, I>
                    <<<grid, block>>>(count_search_keys, sorted_search_keys, count_sorted_keys, sorted_keys, indices);

                SCCD_CUDA_LAST_ERROR();
            } else {
                assert(false);
            }
        }
    }  // namespace device
}  // namespace sccd

#define SCCD_LOWER_BOUND_ALL_TO_ALL_INSTANTIATE(T, I)                                                                 \
    template void sccd::host::lower_bound_all_to_all<T, I>(const ptrdiff_t, const T*, const ptrdiff_t, const T*, I*); \
    template void sccd::device::lower_bound_all_to_all<T, I>(const ptrdiff_t, const T*, const ptrdiff_t, const T*, I*);

SCCD_LOWER_BOUND_ALL_TO_ALL_INSTANTIATE(float, int)
SCCD_LOWER_BOUND_ALL_TO_ALL_INSTANTIATE(double, int)
SCCD_LOWER_BOUND_ALL_TO_ALL_INSTANTIATE(float, ptrdiff_t)
SCCD_LOWER_BOUND_ALL_TO_ALL_INSTANTIATE(double, ptrdiff_t)
