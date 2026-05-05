#include "sccd_broadphase_warp.cuh"

#include <cassert>
#include <cstdint>

#include <cub/device/device_scan.cuh>

#include "sccd_cuda_base.cuh"

#define SCCD_BP_WARP_WARPS_PER_BLOCK 8
#define SCCD_BP_WARP_THREADS (SCCD_BP_WARP_WARPS_PER_BLOCK * SCCD_WARP_SIZE)
#define SCCD_BP_WARP_TILE_SIZE SCCD_BP_WARP_THREADS

namespace sccd {
    namespace device {
        namespace {

            template <typename T>
            __device__ __forceinline__ const T* lower_bound_device(const T* const SCCD_RESTRICT begin,
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

            template <typename T>
            __device__ __forceinline__ const T* upper_bound_device(const T* const SCCD_RESTRICT begin,
                                                                   const T* const SCCD_RESTRICT end,
                                                                   const T val) {
                ptrdiff_t count = end - begin;
                const T* it = begin;

                while (count > 0) {
                    const ptrdiff_t step = count >> 1;
                    const T* const mid = it + step;

                    if (!(val < *mid)) {
                        it = mid + 1;
                        count -= step + 1;
                    } else {
                        count = step;
                    }
                }

                return it;
            }

            template <typename T>
            __device__ __forceinline__ uint32_t disjoint_device(const T aminx,
                                                                const T aminy,
                                                                const T aminz,
                                                                const T amaxx,
                                                                const T amaxy,
                                                                const T amaxz,
                                                                const T bminx,
                                                                const T bminy,
                                                                const T bminz,
                                                                const T bmaxx,
                                                                const T bmaxy,
                                                                const T bmaxz) {
                return aminx > bmaxx | aminy > bmaxy | aminz > bmaxz | bminx > amaxx | bminy > amaxy | bminz > amaxz;
            }

            __device__ __forceinline__ uint64_t warp_broadcast_u64(const uint64_t value) {
                const uint2 parts = make_uint2((uint32_t)value, (uint32_t)(value >> 32));
                const uint32_t lo = __shfl_sync(SCCD_WARP_FULL_MASK, parts.x, 0, SCCD_WARP_SIZE);
                const uint32_t hi = __shfl_sync(SCCD_WARP_FULL_MASK, parts.y, 0, SCCD_WARP_SIZE);
                return ((uint64_t)hi << 32) | lo;
            }

            template <typename I>
            __device__ __forceinline__ I warp_broadcast_index(const I value) {
                if constexpr (sizeof(I) <= sizeof(uint32_t)) {
                    return (I)__shfl_sync(SCCD_WARP_FULL_MASK, (uint32_t)value, 0, SCCD_WARP_SIZE);
                } else {
                    return (I)warp_broadcast_u64((uint64_t)value);
                }
            }

            __device__ __forceinline__ ptrdiff_t warp_reduce_sum_ptrdiff(ptrdiff_t value) {
                uint64_t raw = (uint64_t)value;
                for (int offset = SCCD_WARP_SIZE >> 1; offset > 0; offset >>= 1) {
                    const uint2 parts = make_uint2((uint32_t)raw, (uint32_t)(raw >> 32));
                    const uint32_t lo = __shfl_down_sync(SCCD_WARP_FULL_MASK, parts.x, offset, SCCD_WARP_SIZE);
                    const uint32_t hi = __shfl_down_sync(SCCD_WARP_FULL_MASK, parts.y, offset, SCCD_WARP_SIZE);
                    raw += (((uint64_t)hi << 32) | lo);
                }
                return (ptrdiff_t)raw;
            }

            template <typename T>
            __device__ __forceinline__ T warp_broadcast_value(const T value) {
                return __shfl_sync(SCCD_WARP_FULL_MASK, value, 0, SCCD_WARP_SIZE);
            }

            template <>
            __device__ __forceinline__ double warp_broadcast_value<double>(const double value) {
                return __shfl_sync(SCCD_WARP_FULL_MASK, value, 0, SCCD_WARP_SIZE);
            }

            template <int first_nxe, int second_nxe, typename T, typename I>
            __global__ void count_overlaps_warp_kernel(const int sort_axis,
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
                                                       const T* const SCCD_RESTRICT cummax,
                                                       ptrdiff_t* const SCCD_RESTRICT ccdptr) {
                static_assert(SCCD_BP_WARP_THREADS % SCCD_WARP_SIZE == 0, "Block must contain full warps");

                __shared__ ptrdiff_t warp_begin[SCCD_BP_WARP_WARPS_PER_BLOCK];
                __shared__ ptrdiff_t warp_end[SCCD_BP_WARP_WARPS_PER_BLOCK];
                __shared__ ptrdiff_t block_window[2];
                __shared__ T tile_minx[SCCD_BP_WARP_TILE_SIZE];
                __shared__ T tile_miny[SCCD_BP_WARP_TILE_SIZE];
                __shared__ T tile_minz[SCCD_BP_WARP_TILE_SIZE];
                __shared__ T tile_maxx[SCCD_BP_WARP_TILE_SIZE];
                __shared__ T tile_maxy[SCCD_BP_WARP_TILE_SIZE];
                __shared__ T tile_maxz[SCCD_BP_WARP_TILE_SIZE];

                const int lane = threadIdx.x & (SCCD_WARP_SIZE - 1);
                const int warp_in_block = threadIdx.x / SCCD_WARP_SIZE;
                const ptrdiff_t fi = (ptrdiff_t)blockIdx.x * SCCD_BP_WARP_WARPS_PER_BLOCK + warp_in_block;
                const bool active = fi < first_count;

                if (blockIdx.x == 0 && threadIdx.x == 0) {
                    ccdptr[0] = 0;
                }

                const T* const SCCD_RESTRICT first_minx = first_aabbs[0];
                const T* const SCCD_RESTRICT first_miny = first_aabbs[1];
                const T* const SCCD_RESTRICT first_minz = first_aabbs[2];
                const T* const SCCD_RESTRICT first_maxx = first_aabbs[3];
                const T* const SCCD_RESTRICT first_maxy = first_aabbs[4];
                const T* const SCCD_RESTRICT first_maxz = first_aabbs[5];
                const T* const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
                const T* const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];

                const T* const SCCD_RESTRICT second_minx = second_aabbs[0];
                const T* const SCCD_RESTRICT second_miny = second_aabbs[1];
                const T* const SCCD_RESTRICT second_minz = second_aabbs[2];
                const T* const SCCD_RESTRICT second_maxx = second_aabbs[3];
                const T* const SCCD_RESTRICT second_maxy = second_aabbs[4];
                const T* const SCCD_RESTRICT second_maxz = second_aabbs[5];
                const T* const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
                const T* const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

                T fimin = 0;
                T fimax = 0;
                T aminx = 0;
                T aminy = 0;
                T aminz = 0;
                T amaxx = 0;
                T amaxy = 0;
                T amaxz = 0;
                I first_idxi = 0;
                I ev[first_nxe];
                for (int v = 0; v < first_nxe; ++v) {
                    ev[v] = 0;
                }

                if (active && lane == 0) {
                    fimin = first_xmin[fi];
                    fimax = first_xmax[fi];
                    aminx = first_minx[fi];
                    aminy = first_miny[fi];
                    aminz = first_minz[fi];
                    amaxx = first_maxx[fi];
                    amaxy = first_maxy[fi];
                    amaxz = first_maxz[fi];
                    first_idxi = first_idx[fi];
                    assert(first_idxi >= 0);
                    assert(first_idxi < first_count);

                    for (int v = 0; v < first_nxe; ++v) {
                        ev[v] = first_elements[v][first_idxi * first_stride];
                    }

                    ptrdiff_t begin = lower_bound_device(cummax, cummax + second_count, fimin) - cummax;
                    while (begin < second_count && !(fimin < second_xmax[begin])) {
                        ++begin;
                    }
                    const ptrdiff_t end = upper_bound_device(second_xmin, second_xmin + second_count, fimax) - second_xmin;
                    warp_begin[warp_in_block] = begin;
                    warp_end[warp_in_block] = end;
                } else if (!active && lane == 0) {
                    warp_begin[warp_in_block] = second_count;
                    warp_end[warp_in_block] = 0;
                }

                fimin = warp_broadcast_value(fimin);
                fimax = warp_broadcast_value(fimax);
                aminx = warp_broadcast_value(aminx);
                aminy = warp_broadcast_value(aminy);
                aminz = warp_broadcast_value(aminz);
                amaxx = warp_broadcast_value(amaxx);
                amaxy = warp_broadcast_value(amaxy);
                amaxz = warp_broadcast_value(amaxz);
                first_idxi = warp_broadcast_index(first_idxi);
                for (int v = 0; v < first_nxe; ++v) {
                    ev[v] = warp_broadcast_index(ev[v]);
                }

                __syncthreads();

                if (threadIdx.x == 0) {
                    ptrdiff_t begin = second_count;
                    ptrdiff_t end = 0;
                    for (int w = 0; w < SCCD_BP_WARP_WARPS_PER_BLOCK; ++w) {
                        begin = min(begin, warp_begin[w]);
                        end = max(end, warp_end[w]);
                    }
                    block_window[0] = begin;
                    block_window[1] = end;
                }

                __syncthreads();

                const ptrdiff_t begin = warp_begin[warp_in_block];
                const ptrdiff_t end = warp_end[warp_in_block];

                ptrdiff_t count = 0;
                const ptrdiff_t block_begin = block_window[0];
                const ptrdiff_t block_end = block_window[1];

                for (ptrdiff_t tile = block_begin; tile < block_end; tile += SCCD_BP_WARP_TILE_SIZE) {
                    const ptrdiff_t tile_len = min((ptrdiff_t)SCCD_BP_WARP_TILE_SIZE, block_end - tile);

                    if (threadIdx.x < tile_len) {
                        const ptrdiff_t j = tile + threadIdx.x;
                        tile_minx[threadIdx.x] = second_minx[j];
                        tile_miny[threadIdx.x] = second_miny[j];
                        tile_minz[threadIdx.x] = second_minz[j];
                        tile_maxx[threadIdx.x] = second_maxx[j];
                        tile_maxy[threadIdx.x] = second_maxy[j];
                        tile_maxz[threadIdx.x] = second_maxz[j];
                    }

                    __syncthreads();

                    if (active && begin < end) {
                        const ptrdiff_t first_j = max(begin, tile);
                        const ptrdiff_t last_j = min(end, tile + tile_len);
                        for (ptrdiff_t j = first_j + lane; j < last_j; j += SCCD_WARP_SIZE) {
                            const int tile_lane = (int)(j - tile);
                            if (disjoint_device(aminx,
                                                aminy,
                                                aminz,
                                                amaxx,
                                                amaxy,
                                                amaxz,
                                                tile_minx[tile_lane],
                                                tile_miny[tile_lane],
                                                tile_minz[tile_lane],
                                                tile_maxx[tile_lane],
                                                tile_maxy[tile_lane],
                                                tile_maxz[tile_lane])) {
                                continue;
                            }

                            bool share = false;
                            const I jidx = second_idx[j];
                            assert(jidx >= 0);
                            assert(jidx < second_count);
                            if constexpr (second_nxe > 1) {
                                I sev[second_nxe];
                                for (int v = 0; v < second_nxe; ++v) {
                                    sev[v] = second_elements[v][jidx * second_stride];
                                }
                                for (int a = 0; a < first_nxe; ++a) {
                                    for (int b = 0; b < second_nxe; ++b) {
                                        share |= ev[a] == sev[b];
                                    }
                                }
                            } else {
                                for (int a = 0; a < first_nxe; ++a) {
                                    share |= ev[a] == jidx;
                                }
                            }

                            count += share ? 0 : 1;
                        }
                    }

                    __syncthreads();
                }

                count = warp_reduce_sum_ptrdiff(count);

                if (active && lane == 0) {
                    ccdptr[fi + 1] = count;
                }
            }

        }  // namespace

        template <int first_nxe, int second_nxe, typename T, typename I>
        void count_overlaps_warp(const int sort_axis,
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
                                 const T* const SCCD_RESTRICT cummax) {
            SCCD_CUDA_LAST_ERROR();

            if (first_count <= 0) {
                SCCD_CHECK_CUDA(cudaMemset(ccdptr, 0, sizeof(*ccdptr)));
                return;
            }

            dim3 block(SCCD_BP_WARP_THREADS);
            dim3 grid((first_count + SCCD_BP_WARP_WARPS_PER_BLOCK - 1) / SCCD_BP_WARP_WARPS_PER_BLOCK);
            count_overlaps_warp_kernel<first_nxe, second_nxe, T, I><<<grid, block>>>(sort_axis,
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
                                                                                     cummax,
                                                                                     ccdptr);
            SCCD_CUDA_LAST_ERROR();

            void* tmp_storage = nullptr;
            size_t tmp_storage_bytes = 0;
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveSum(nullptr, tmp_storage_bytes, ccdptr, ccdptr, first_count + 1));
            SCCD_CHECK_CUDA(cudaMalloc(&tmp_storage, tmp_storage_bytes));
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveSum(tmp_storage, tmp_storage_bytes, ccdptr, ccdptr, first_count + 1));
            SCCD_CHECK_CUDA(cudaFree(tmp_storage));

            SCCD_CUDA_LAST_ERROR();
        }

#define INSTANTIATE_COUNT_OVERLAPS_WARP(FIRST_NXE, SECOND_NXE, T, I)                                      \
    template void sccd::device::count_overlaps_warp<FIRST_NXE, SECOND_NXE, T, I>(const int sort_axis,     \
                                                                                  const ptrdiff_t first_count, \
                                                                                  T** const SCCD_RESTRICT first_aabbs, \
                                                                                  I* const SCCD_RESTRICT first_idx, \
                                                                                  const ptrdiff_t first_stride, \
                                                                                  I** const SCCD_RESTRICT first_elements, \
                                                                                  const ptrdiff_t second_count, \
                                                                                  T** const SCCD_RESTRICT second_aabbs, \
                                                                                  I* const SCCD_RESTRICT second_idx, \
                                                                                  const ptrdiff_t second_stride, \
                                                                                  I** const SCCD_RESTRICT second_elements, \
                                                                                  ptrdiff_t* const SCCD_RESTRICT ccdptr, \
                                                                                  const T* const SCCD_RESTRICT cummax)

        INSTANTIATE_COUNT_OVERLAPS_WARP(3, 1, float, int32_t);
        INSTANTIATE_COUNT_OVERLAPS_WARP(3, 1, float, int64_t);
        INSTANTIATE_COUNT_OVERLAPS_WARP(3, 1, double, int32_t);
        INSTANTIATE_COUNT_OVERLAPS_WARP(3, 1, double, int64_t);

#undef INSTANTIATE_COUNT_OVERLAPS_WARP

    }  // namespace device
}  // namespace sccd

#undef SCCD_BP_WARP_WARPS_PER_BLOCK
#undef SCCD_BP_WARP_THREADS
#undef SCCD_BP_WARP_TILE_SIZE
