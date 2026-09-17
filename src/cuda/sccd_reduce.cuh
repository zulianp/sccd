#ifndef SCCD_REDUCE_CUH
#define SCCD_REDUCE_CUH

#include "sccd_base.hpp"
#include "sccd_cuda_base.cuh"

#include <cassert>
#include <cstddef>

namespace sccd {
    namespace device {

        inline __device__ unsigned int lane_id(const int thIdx) { return thIdx % SCCD_WARP_SIZE; }
        inline __device__ unsigned int lane_id() { return threadIdx.x % SCCD_WARP_SIZE; }

        template <typename T>
        struct device_max_op {
            __host__ __device__ __forceinline__ T operator()(const T& a, const T& b) const { return a > b ? a : b; }
        };

        template <typename T>
        __device__ T warp_reduce_32(const T in) {
            static_assert(SCCD_WARP_SIZE == 32, "Only implemented for CUDA with warp size of 32!");
            T out = in;
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 16, SCCD_WARP_SIZE);  // 0-16, 1-17, ..., 15-31
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 8, SCCD_WARP_SIZE);   // 0-8, ..., 1-7, ..., 23-31
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 4, SCCD_WARP_SIZE);
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 2, SCCD_WARP_SIZE);
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 1, SCCD_WARP_SIZE);
            return out;
        }

        /**
         * \brief Sum `val` across the block and add the total to `result`.
         *
         * **Every thread of the block must reach this call.** Lane 0 of each warp
         * writes that warp's partial sum to `block_accumulator[warp_id]`, and warp
         * 0 then reads all `n_warps` slots; a warp that returned earlier leaves
         * its slot holding whatever that multiprocessor's shared memory last held,
         * and the sum silently absorbs it. A caller guarding on an element index
         * therefore passes an identity for threads past the end rather than
         * returning, which leaves the sum unchanged and every slot written. The
         * __syncthreads and the full-mask shuffles need the same uniformity.
         *
         * `block_accumulator` is shared memory with at least one entry per warp.
         */
        template <typename T>
        __device__ void block_reduce_to_gmem(const T val,
                                             T* const SCCD_RESTRICT block_accumulator,
                                             T* const SCCD_RESTRICT result) {
            T acc = warp_reduce_32(val);
            const unsigned int warp_id = threadIdx.x / SCCD_WARP_SIZE;
            const unsigned int lid = lane_id();
            const unsigned int n_warps = (blockDim.x + SCCD_WARP_SIZE - 1) / SCCD_WARP_SIZE;

            if (!lid) {
                block_accumulator[warp_id] = acc;
            }

            __syncthreads();

            if (!warp_id) {
                assert(warp_id < n_warps);
                acc = lid < n_warps ? block_accumulator[lid] : 0;
                acc = warp_reduce_32(acc);

                if (!threadIdx.x) {
                    assert(acc == acc);
                    atomicAdd(result, acc);
                }
            }

            __syncthreads();
        }

    }  // namespace device
}  // namespace sccd

#endif  // SCCD_REDUCE_CUH
