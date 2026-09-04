#ifndef SCCD_SPIKES_DEAD_CUH
#define SCCD_SPIKES_DEAD_CUH

// Device code that shipped in the library and was called by nothing.
//
// The *_overlaps_with_starts family is the largest piece: a complete second
// broad-phase collect path, roughly 160 lines of kernel, that was compiled into
// every CUDA build through 8 explicit template instantiations and reached from
// nowhere. The only spikes/ hit for the name is a different function that
// happens to share it.
//
// spikes/README.md applies: not built by default, not installed, not covered by
// the correctness gate, deletable without notice.
//
// This is a header rather than a .cu so it needs no separate compilation unit;
// nothing includes it, which is the point.

#include "sccd_base.hpp"
#include "sccd_cuda_base.cuh"

#include <cstddef>
#include <cstdint>

namespace sccd {
    namespace device {
        namespace dead {

        // ---- count_overlaps_with_starts decl  (from sccd_broadphase.cuh)

        template <int first_nxe, int second_nxe, typename T, typename I>
        void count_overlaps_with_starts(const int sort_axis,
                                        const ptrdiff_t first_count,
                                        T** const SCCD_RESTRICT first_aabbs,
                                        I* const SCCD_RESTRICT first_idx,
                                        const ptrdiff_t first_stride,
                                        I** const SCCD_RESTRICT first_elements,
                                        const ptrdiff_t second_count,
                                        T** const SCCD_RESTRICT second_aabbs,
                                        I* const SCCD_RESTRICT second_idx,
                                        const ptrdiff_t second_stride,
                                        I** const SCCD_RESTRICT second_elements,
                                        ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                        const I* const SCCD_RESTRICT starts);

        // ---- collect_overlaps_with_starts decl  (from sccd_broadphase.cuh)

        template <int first_nxe, int second_nxe, typename T, typename I>
        void collect_overlaps_with_starts(const int sort_axis,
                                          const ptrdiff_t first_count,
                                          T** const SCCD_RESTRICT first_aabbs,
                                          I* const SCCD_RESTRICT first_idx,
                                          const ptrdiff_t first_stride,
                                          I** SCCD_RESTRICT const first_elements,
                                          const ptrdiff_t second_count,
                                          T** const SCCD_RESTRICT second_aabbs,
                                          I* const SCCD_RESTRICT second_idx,
                                          const ptrdiff_t second_stride,
                                          I** SCCD_RESTRICT const second_elements,
                                          const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                          const I* const SCCD_RESTRICT starts,
                                          I* SCCD_RESTRICT foverlap,
                                          I* SCCD_RESTRICT noverlap);

        // ---- count_overlaps_with_starts_kernel  (from sccd_broadphase.cu)

        template <int first_nxe, int second_nxe, typename T, typename I>
        __global__ void count_overlaps_with_starts_kernel(const int sort_axis,
                                                          const ptrdiff_t first_count,
                                                          T** const SCCD_RESTRICT first_aabbs,
                                                          I* const SCCD_RESTRICT first_idx,
                                                          const ptrdiff_t first_stride,
                                                          I** const SCCD_RESTRICT first_elements,
                                                          const ptrdiff_t second_count,
                                                          T** const SCCD_RESTRICT second_aabbs,
                                                          I* const SCCD_RESTRICT second_idx,
                                                          const ptrdiff_t second_stride,
                                                          I** const SCCD_RESTRICT second_elements,
                                                          ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                          const I* const SCCD_RESTRICT starts) {
            ptrdiff_t fi = blockIdx.x * blockDim.x + threadIdx.x;
            if (fi >= first_count) return;
            if (fi == 0) {
                ccdptr[0] = 0;
            }

            ptrdiff_t begin = starts[fi];
            ccdptr[fi + 1] = count_overlaps_kernel_aux<first_nxe, second_nxe, T, I>(fi,
                                                                                    begin,
                                                                                    sort_axis,
                                                                                    first_count,
                                                                                    first_aabbs,
                                                                                    first_idx,
                                                                                    first_stride,
                                                                                    first_elements,
                                                                                    second_count,
                                                                                    second_aabbs,
                                                                                    second_idx,
                                                                                    second_stride,
                                                                                    second_elements);
        }

        // ---- count_overlaps_with_starts  (from sccd_broadphase.cu)

        template <int first_nxe, int second_nxe, typename T, typename I>
        void count_overlaps_with_starts(const int sort_axis,
                                        const ptrdiff_t first_count,
                                        T** const SCCD_RESTRICT first_aabbs,
                                        I* const SCCD_RESTRICT first_idx,
                                        const ptrdiff_t first_stride,
                                        I** const SCCD_RESTRICT first_elements,
                                        const ptrdiff_t second_count,
                                        T** const SCCD_RESTRICT second_aabbs,
                                        I* const SCCD_RESTRICT second_idx,
                                        const ptrdiff_t second_stride,
                                        I** const SCCD_RESTRICT second_elements,
                                        ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                        const I* const SCCD_RESTRICT starts) {
            SCCD_CUDA_LAST_ERROR();

            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((first_count + block.x - 1) / block.x);
            count_overlaps_with_starts_kernel<first_nxe, second_nxe, T, I><<<grid, block>>>(sort_axis,
                                                                                            first_count,
                                                                                            first_aabbs,
                                                                                            first_idx,
                                                                                            first_stride,
                                                                                            first_elements,
                                                                                            second_count,
                                                                                            second_aabbs,
                                                                                            second_idx,
                                                                                            second_stride,
                                                                                            second_elements,
                                                                                            ccdptr,
                                                                                            starts);

            SCCD_CUDA_LAST_ERROR();

            size_t tmp_storage_bytes = 0;
            SCCD_CHECK_CUDA(cub::DeviceScan::InclusiveSum(nullptr, tmp_storage_bytes, ccdptr, ccdptr, first_count + 1));
            void* const tmp_storage = workspace(WorkspaceSlot::TempStorage).get(tmp_storage_bytes);
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveSum(tmp_storage, tmp_storage_bytes, ccdptr, ccdptr, first_count + 1));

            SCCD_CUDA_LAST_ERROR();
        }

        // ---- collect_overlaps_with_starts_kernel  (from sccd_broadphase.cu)

        template <int first_nxe, int second_nxe, typename T, typename I>
        __global__ void collect_overlaps_with_starts_kernel(const int sort_axis,
                                                            const ptrdiff_t first_count,
                                                            T** const SCCD_RESTRICT first_aabbs,
                                                            I* const SCCD_RESTRICT first_idx,
                                                            const ptrdiff_t first_stride,
                                                            I** SCCD_RESTRICT const first_elements,
                                                            const ptrdiff_t second_count,
                                                            T** const SCCD_RESTRICT second_aabbs,
                                                            I* const SCCD_RESTRICT second_idx,
                                                            const ptrdiff_t second_stride,
                                                            I** SCCD_RESTRICT const second_elements,
                                                            const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                            const I* const SCCD_RESTRICT starts,
                                                            I* SCCD_RESTRICT foverlap,
                                                            I* SCCD_RESTRICT noverlap) {
            ptrdiff_t fi = blockIdx.x * blockDim.x + threadIdx.x;
            if (fi >= first_count) return;
            const ptrdiff_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
            if (!expected_count) return;

            ptrdiff_t begin = starts[fi];

            collect_overlaps_kernel_aux<first_nxe, second_nxe, T, I>(fi,
                                                                     begin,
                                                                     sort_axis,
                                                                     first_count,
                                                                     first_aabbs,
                                                                     first_idx,
                                                                     first_stride,
                                                                     first_elements,
                                                                     second_count,
                                                                     second_aabbs,
                                                                     second_idx,
                                                                     second_stride,
                                                                     second_elements,
                                                                     ccdptr,
                                                                     foverlap,
                                                                     noverlap);
        }

        // ---- collect_overlaps_with_starts  (from sccd_broadphase.cu)

        template <int first_nxe, int second_nxe, typename T, typename I>
        void collect_overlaps_with_starts(const int sort_axis,
                                          const ptrdiff_t first_count,
                                          T** const SCCD_RESTRICT first_aabbs,
                                          I* const SCCD_RESTRICT first_idx,
                                          const ptrdiff_t first_stride,
                                          I** SCCD_RESTRICT const first_elements,
                                          const ptrdiff_t second_count,
                                          T** const SCCD_RESTRICT second_aabbs,
                                          I* const SCCD_RESTRICT second_idx,
                                          const ptrdiff_t second_stride,
                                          I** SCCD_RESTRICT const second_elements,
                                          const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                          const I* const SCCD_RESTRICT starts,
                                          I* SCCD_RESTRICT foverlap,
                                          I* SCCD_RESTRICT noverlap) {
            SCCD_CUDA_LAST_ERROR();

            dim3 block(SCCD_BP_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((first_count + block.x - 1) / block.x);
            collect_overlaps_with_starts_kernel<first_nxe, second_nxe, T, I><<<grid, block>>>(sort_axis,
                                                                                              first_count,
                                                                                              first_aabbs,
                                                                                              first_idx,
                                                                                              first_stride,
                                                                                              first_elements,
                                                                                              second_count,
                                                                                              second_aabbs,
                                                                                              second_idx,
                                                                                              second_stride,
                                                                                              second_elements,
                                                                                              ccdptr,
                                                                                              starts,
                                                                                              foverlap,
                                                                                              noverlap);

            SCCD_CUDA_LAST_ERROR();
        }

        // ---- reset_batch_narrow_phase_kernel  (from sccd_narrowphase.cu)
        // a __global__ that is never launched. The driver zeroes the
        // two counters with one cudaMemsetAsync instead.

        __global__ void reset_batch_narrow_phase_kernel(int* SCCD_RESTRICT g_top, int* SCCD_RESTRICT g_request) {
            *g_top = 0;
            *g_request = 0;
        }

        // ---- release_slot  (from sccd_narrowphase.cu)
        // orphaned by the double-buffered work queue in 09757ac3, which
        // replaced the stack discipline with a bump allocator and a claim cursor.

        static inline __device__ int release_slot(int* SCCD_RESTRICT counter) {
            int old = atomicAdd(counter, 0);
            while (true) {
                if (old <= 0) return -1;
                const int prev = atomicCAS(counter, old, old - 1);
                if (prev == old) return old - 1;
                old = prev;
            }
        }

        // ---- block_max_to_root  (from sccd_reduce.cuh)
        // block-wide max reduction into thread 0; no caller.

        template <typename T>
        __device__ void block_max_to_root(const int thIdx,
                                          const int blDim,
                                          const T val,
                                          T* const SCCD_RESTRICT block_accumulator,
                                          T* const SCCD_RESTRICT result) {
            T acc = warp_max_32(val);

            const unsigned int warp_id = thIdx / SCCD_WARP_SIZE;
            const unsigned int lid = lane_id(thIdx);
            const unsigned int n_warps = (blDim + SCCD_WARP_SIZE - 1) / SCCD_WARP_SIZE;

            if (!lid) {
                block_accumulator[warp_id] = acc;
            }

            __syncthreads();

            if (!warp_id) {
                assert(warp_id < n_warps);
                acc = lid < n_warps ? block_accumulator[lid] : block_accumulator[0];
                acc = warp_max_32(acc);

                if (!thIdx) {
                    assert(acc == acc);
                    *result = acc;
                }
            }
        }

        // ---- broadcast_to_block  (from sccd_reduce.cuh)
        // no caller, and internally dead: it computes warp_id, lid,
        // root_warp and root_lid and then uses none of them.

        template <typename T>
        __device__ __forceinline__ T broadcast_to_block(const int root, const T val) {
            const int thIdx = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * threadIdx.z);
            const int warp_id = thIdx / SCCD_WARP_SIZE;
            const int lid = thIdx % SCCD_WARP_SIZE;
            const int root_warp = root / SCCD_WARP_SIZE;
            const int root_lid = root % SCCD_WARP_SIZE;

            __shared__ T warp_value;
            if (root == 0) {
                warp_value = val;
            }
            __syncthreads();

            return warp_value;
        }

        }  // namespace dead
    }  // namespace device
}  // namespace sccd

#endif  // SCCD_SPIKES_DEAD_CUH
