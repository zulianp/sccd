#include "sccd_broadphase.cuh"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_scan.cuh>

#include "sccd_base.hpp"

#define SCCD_N_WARPS_PER_BLOCK 8
#define SCCD_WARP_SIZE 32
#define SCCD_WARP_FULL_MASK 0xffffffff

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

namespace sccd {
    namespace device {

        void cuda_check(const cudaError_t error) {
            if (error != cudaSuccess) {
                fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
                exit(1);
            }
        }

#define SCCD_CHECK_CUDA(error) \
    cudaDeviceSynchronize();   \
    cuda_check(error)

#define SCCD_CUDA_LAST_ERROR() SCCD_CHECK_CUDA(cudaGetLastError())

        inline __device__ unsigned int lane_id() { return threadIdx.x % SCCD_WARP_SIZE; }

        template <typename T>
        __device__ T warp_reduce_32(const T in) {
            static_assert(SCCD_WARP_SIZE == 32, "Only implemented for CUDA!");
            T out = in;
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 16, SCCD_WARP_SIZE);  // 0-16, 1-17, ..., 15-31
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 8, SCCD_WARP_SIZE);   // 0-8, ..., 1-7, ..., 23-31
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 4, SCCD_WARP_SIZE);
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 2, SCCD_WARP_SIZE);
            out += __shfl_xor_sync(SCCD_WARP_FULL_MASK, out, 1, SCCD_WARP_SIZE);
            return out;
        }

        template <typename T>
        __device__ void t_warp_reduce(const T val,
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
                assert(warp_id < SCCD_N_WARPS_PER_BLOCK);
                acc = lid < n_warps ? block_accumulator[lid] : 0;
                acc = warp_reduce_32(acc);

                if (!threadIdx.x) {
                    assert(acc == acc);
                    atomicAdd(result, acc);
                }
            }
        }

        template <typename T>
        __global__ void choose_axis_mean_kernel(const int dim,
                                                const ptrdiff_t n,
                                                const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                                T* const SCCD_RESTRICT mean) {
            ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= n) return;

            T local_mean[3] = {0};

            for (int d = 0; d < dim; d++) {
                const T p0 = aabbs[d][i];
                const T p1 = aabbs[dim + d][i];
                const T p = (p0 + p1) / 2;
                local_mean[d] += p;
            }

            __shared__ T block_accumulator[SCCD_N_WARPS_PER_BLOCK];
            for (int d = 0; d < dim; d++) {
                t_warp_reduce<T>(local_mean[d], block_accumulator, &mean[d]);
            }
        }

        template <typename T>
        __global__ void choose_axis_var_kernel(const int dim,
                                               const ptrdiff_t n,
                                               const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                               T* const SCCD_RESTRICT mean,
                                               T* const SCCD_RESTRICT var) {
            ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= n) return;

            T local_var[3] = {0};
            for (int d = 0; d < dim; d++) {
                const T m = mean[d] / n;
                const T p = (aabbs[d][i] + aabbs[dim + d][i]) / 2;
                local_var[d] += (p - m) * (p - m);
            }

            __shared__ T block_accumulator[SCCD_N_WARPS_PER_BLOCK];
            for (int d = 0; d < dim; d++) {
                t_warp_reduce<T>(local_var[d], block_accumulator, &var[d]);
            }
        }

        template <typename T>
        int choose_axis(const int dim, const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs) {
            dim3 block(SCCD_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((n + block.x - 1) / block.x);

            T* mean = nullptr;
            cudaMalloc(&mean, dim * sizeof(T));
            T* var = nullptr;
            cudaMalloc(&var, dim * sizeof(T));

            cudaMemset(mean, 0, dim * sizeof(T));
            cudaMemset(var, 0, dim * sizeof(T));

            choose_axis_mean_kernel<T><<<grid, block>>>(dim, n, aabbs, mean);
            choose_axis_var_kernel<T><<<grid, block>>>(dim, n, aabbs, mean, var);

            cudaError_t error = cudaGetLastError();

            T* hvar = (T*)malloc(dim * sizeof(T));
            cudaMemcpy(hvar, var, dim * sizeof(T), cudaMemcpyDeviceToHost);

            int fargmax = 0;
            T fmax = hvar[0];

            for (int d = 1; d < 3; d++) {
                if (fmax < hvar[d]) {
                    fmax = hvar[d];
                    fargmax = d;
                }
            }

            cudaFree(mean);
            cudaFree(var);
            free(hvar);

            if (error != cudaSuccess) {
                fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
                exit(1);
            }

            return fargmax;
        }

        template <typename T>
        __global__ void enumerate_kernel(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx) {
            ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= end - begin) return;
            idx[i] = begin + i;
        }

        template <typename T>
        void enumerate(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx) {
            dim3 block(SCCD_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((end - begin + block.x - 1) / block.x);
            enumerate_kernel<T><<<grid, block>>>(begin, end, idx);
        }

        template <typename T, typename I>
        __global__ void permute_kernel(const ptrdiff_t n,
                                       const I* const SCCD_RESTRICT idx,
                                       const T* const SCCD_RESTRICT src,
                                       T* const SCCD_RESTRICT dst) {
            const ptrdiff_t i = blockIdx.x * blockDim.x + threadIdx.x;
            if (i >= n) return;
            dst[i] = src[idx[i]];
        }

        bool is_ptr_device(const void* ptr) {
            cudaPointerAttributes attributes;
            cudaError_t err = cudaPointerGetAttributes(&attributes, ptr);

            if (err != cudaSuccess) {
                fprintf(stderr, "cudaPointerGetAttributes failed: %s\n", cudaGetErrorString(err));
            }

#if CUDART_VERSION >= 10000
            // CUDA 10.0 and newer
            return attributes.type == cudaMemoryTypeDevice;
#else
            return attributes.memoryType == cudaMemoryTypeDevice;
#endif
        }

        template <typename T, typename I>
        void sort_along_axis(const int dim,
                             const ptrdiff_t n,
                             const int sort_axis,
                             T** const SCCD_RESTRICT arrays,
                             I* const SCCD_RESTRICT idx,
                             T* const SCCD_RESTRICT scratch) {
            if (n <= 0) return;

            SCCD_CUDA_LAST_ERROR();

            if (n == 1) {
                enumerate<I>(0, n, idx);
                SCCD_CUDA_LAST_ERROR();
                return;
            }

            T* host_arrays[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

            if (is_ptr_device(arrays)) {
                SCCD_CHECK_CUDA(cudaMemcpy(host_arrays, arrays, 2 * dim * sizeof(T*), cudaMemcpyDeviceToHost));
            } else {
                for (int d = 0; d < 2 * dim; d++) {
                    host_arrays[d] = arrays[d];
                }
            }

            const T* const SCCD_RESTRICT x = host_arrays[sort_axis];

            I* tmp_idx = nullptr;
            void* tmp_storage = nullptr;
            size_t tmp_storage_bytes = 0;

            SCCD_CHECK_CUDA(cudaMalloc(&tmp_idx, n * sizeof(I)));
            enumerate<I>(0, n, tmp_idx);
            SCCD_CUDA_LAST_ERROR();

            SCCD_CHECK_CUDA(cub::DeviceRadixSort::SortPairs(nullptr, tmp_storage_bytes, x, scratch, tmp_idx, idx, n));
            SCCD_CHECK_CUDA(cudaMalloc(&tmp_storage, tmp_storage_bytes));
            SCCD_CHECK_CUDA(
                cub::DeviceRadixSort::SortPairs(tmp_storage, tmp_storage_bytes, x, scratch, tmp_idx, idx, n));

            SCCD_CHECK_CUDA(cudaMemcpy(host_arrays[sort_axis], scratch, n * sizeof(T), cudaMemcpyDeviceToDevice));

            dim3 block(SCCD_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((n + block.x - 1) / block.x);
            for (int d = 0; d < 2 * dim; d++) {
                if (d == sort_axis) continue;
                permute_kernel<T, I><<<grid, block>>>(n, idx, host_arrays[d], scratch);
                SCCD_CUDA_LAST_ERROR();
                SCCD_CHECK_CUDA(cudaMemcpy(host_arrays[d], scratch, n * sizeof(T), cudaMemcpyDeviceToDevice));
            }

            SCCD_CHECK_CUDA(cudaFree(tmp_storage));
            SCCD_CHECK_CUDA(cudaFree(tmp_idx));

            SCCD_CUDA_LAST_ERROR();
        }

        template <typename T>
        static inline __device__ void compute_candidate_window_progressive(const T fimin,
                                                                           const T fimax,
                                                                           const T* const SCCD_RESTRICT second_xmax,
                                                                           const T* const SCCD_RESTRICT second_xmin,
                                                                           const ptrdiff_t second_count,
                                                                           ptrdiff_t& begin,
                                                                           ptrdiff_t& end) {
            for (; begin < second_count; ++begin) {
                if (fimin < second_xmax[begin]) {
                    break;
                }
            }
            end = begin;
            for (; end < second_count; ++end) {
                if (fimax < second_xmin[end]) {
                    break;
                }
            }
        }

        template <int nxe, typename I>
        static inline __device__ void load_ev(I** const SCCD_RESTRICT elements,
                                              const I elem_idx,
                                              const ptrdiff_t stride,
                                              I (&out)[nxe]) {
            for (int v = 0; v < nxe; ++v) {
                out[v] = elements[v][elem_idx * stride];
            }
        }

        template <int n1, int n2, typename I>
        static inline __device__ bool shares_vertex(const I (&a)[n1], const I (&b)[n2]) {
            for (int i = 0; i < n1; ++i) {
                for (int j = 0; j < n2; ++j) {
                    if (a[i] == b[j]) {
                        return true;
                    }
                }
            }
            return false;
        }

        template <typename T>
        static inline __device__ uint32_t disjoint(const T aminx,
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

        template <int nxe, typename T, typename I>
        __global__ void count_self_overlaps_kernel(const int sort_axis,
                                                   const ptrdiff_t element_count,
                                                   T** const SCCD_RESTRICT aabbs,
                                                   I* const SCCD_RESTRICT idx,
                                                   const ptrdiff_t stride,
                                                   I** const SCCD_RESTRICT elements,
                                                   ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            ptrdiff_t fi = blockIdx.x * blockDim.x + threadIdx.x;
            if (fi >= element_count) return;
            if (fi == 0) {
                ccdptr[0] = 0;
            }

            const T* const SCCD_RESTRICT xmin = aabbs[sort_axis];
            const T* const SCCD_RESTRICT xmax = aabbs[3 + sort_axis];

            const T fimin = xmin[fi];
            const T fimax = xmax[fi];
            const I idxi = idx[fi];

            assert(idxi >= 0);
            assert(idxi < element_count);

            I ev[nxe];
            for (int v = 0; v < nxe; v++) {
                ev[v] = elements[v][idxi * stride];
            }

            ptrdiff_t begin = fi + 1;
            ptrdiff_t end = begin;
            compute_candidate_window_progressive(fimin, fimax, xmax, xmin, element_count, begin, end);

            if (begin >= end) {
                ccdptr[fi + 1] = 0;
                return;
            }

            ptrdiff_t count = 0;
            const T aminx = aabbs[0][fi];
            const T aminy = aabbs[1][fi];
            const T aminz = aabbs[2][fi];
            const T amaxx = aabbs[3][fi];
            const T amaxy = aabbs[4][fi];
            const T amaxz = aabbs[5][fi];

            for (ptrdiff_t j = begin; j < end; ++j) {
                if (disjoint(aminx,
                             aminy,
                             aminz,
                             amaxx,
                             amaxy,
                             amaxz,
                             aabbs[0][j],
                             aabbs[1][j],
                             aabbs[2][j],
                             aabbs[3][j],
                             aabbs[4][j],
                             aabbs[5][j])) {
                    continue;
                }
                const I jidx = idx[j];
                I sev[nxe];
                load_ev<nxe>(elements, jidx, stride, sev);
                const bool share = shares_vertex<nxe, nxe>(ev, sev);
                count += share ? 0 : 1;
            }

            ccdptr[fi + 1] = count;
        }

        template <int nxe, typename T, typename I>
        void count_self_overlaps(const int sort_axis,
                                 const ptrdiff_t element_count,
                                 T** const SCCD_RESTRICT aabbs,
                                 I* const SCCD_RESTRICT idx,
                                 const ptrdiff_t stride,
                                 I** const SCCD_RESTRICT elements,
                                 ptrdiff_t* const SCCD_RESTRICT ccdptr) {
            SCCD_CUDA_LAST_ERROR();

            if (element_count <= 0) {
                SCCD_CHECK_CUDA(cudaMemset(ccdptr, 0, sizeof(*ccdptr)));
                return;
            }

            dim3 block(SCCD_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((element_count + block.x - 1) / block.x);
            count_self_overlaps_kernel<nxe, T, I>
                <<<grid, block>>>(sort_axis, element_count, aabbs, idx, stride, elements, ccdptr);
            SCCD_CUDA_LAST_ERROR();

            void* tmp_storage = nullptr;
            size_t tmp_storage_bytes = 0;
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveSum(nullptr, tmp_storage_bytes, ccdptr, ccdptr, element_count + 1));
            SCCD_CHECK_CUDA(cudaMalloc(&tmp_storage, tmp_storage_bytes));
            SCCD_CHECK_CUDA(
                cub::DeviceScan::InclusiveSum(tmp_storage, tmp_storage_bytes, ccdptr, ccdptr, element_count + 1));
            SCCD_CHECK_CUDA(cudaFree(tmp_storage));

            SCCD_CUDA_LAST_ERROR();
        }

        template <int nxe, typename T, typename I>
        __global__ void collect_self_overlaps_kernel(const int sort_axis,
                                                     const ptrdiff_t element_count,
                                                     T** const SCCD_RESTRICT aabbs,
                                                     I* const SCCD_RESTRICT idx,
                                                     const ptrdiff_t stride,
                                                     I** const elements,
                                                     const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                                     I* SCCD_RESTRICT foverlap,
                                                     I* SCCD_RESTRICT noverlap) {
            ptrdiff_t fi = blockIdx.x * blockDim.x + threadIdx.x;
            if (fi >= element_count) return;

            const ptrdiff_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
            if (!expected_count) return;

            const T fimin = aabbs[sort_axis][fi];
            const T fimax = aabbs[3 + sort_axis][fi];
            const I idxi = idx[fi];

            I ev[nxe];
            for (int v = 0; v < nxe; v++) {
                ev[v] = elements[v][idxi * stride];
            }

            ptrdiff_t begin = fi + 1;
            ptrdiff_t end = begin;
            compute_candidate_window_progressive(
                fimin, fimax, aabbs[3 + sort_axis], aabbs[sort_axis], element_count, begin, end);

            if (begin >= end) {
                return;
            }

            I* SCCD_RESTRICT const first_local_elements = &foverlap[ccdptr[fi]];
            I* SCCD_RESTRICT const second_local_elements = &noverlap[ccdptr[fi]];

            ptrdiff_t count = 0;
            const T aminx = aabbs[0][fi];
            const T aminy = aabbs[1][fi];
            const T aminz = aabbs[2][fi];
            const T amaxx = aabbs[3][fi];
            const T amaxy = aabbs[4][fi];
            const T amaxz = aabbs[5][fi];

            for (ptrdiff_t j = begin; j < end; ++j) {
                if (disjoint(aminx,
                             aminy,
                             aminz,
                             amaxx,
                             amaxy,
                             amaxz,
                             aabbs[0][j],
                             aabbs[1][j],
                             aabbs[2][j],
                             aabbs[3][j],
                             aabbs[4][j],
                             aabbs[5][j])) {
                    continue;
                }
                const I jidx = idx[j];
                I sev[nxe];
                load_ev<nxe>(elements, jidx, stride, sev);
                if (shares_vertex<nxe, nxe>(ev, sev)) {
                    continue;
                }

                first_local_elements[count] = MIN(idxi, jidx);
                second_local_elements[count] = MAX(idxi, jidx);

                count++;
            }

            assert(expected_count == count);
        }

        template <int nxe, typename T, typename I>
        void collect_self_overlaps(const int sort_axis,
                                   const ptrdiff_t element_count,
                                   T** const SCCD_RESTRICT aabbs,
                                   I* const SCCD_RESTRICT idx,
                                   const ptrdiff_t stride,
                                   I** const elements,
                                   const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                   I* SCCD_RESTRICT foverlap,
                                   I* SCCD_RESTRICT noverlap) {
            SCCD_CUDA_LAST_ERROR();

            dim3 block(SCCD_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((element_count + block.x - 1) / block.x);
            collect_self_overlaps_kernel<nxe, T, I>
                <<<grid, block>>>(sort_axis, element_count, aabbs, idx, stride, elements, ccdptr, foverlap, noverlap);

            SCCD_CUDA_LAST_ERROR();
        }

        template <int first_nxe, int second_nxe, typename T, typename I>
        __global__ void count_overlaps_kernel(const int sort_axis,
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
                                              ptrdiff_t* const SCCD_RESTRICT ccdptr) {}

        template <int first_nxe, int second_nxe, typename T, typename I>
        void count_overlaps(const int sort_axis,
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
                            ptrdiff_t* const SCCD_RESTRICT ccdptr) {}

    }  // namespace device
}  // namespace sccd

#define INSTANTIATE_CHOOSE_AXIS(T)             \
    template int sccd::device::choose_axis<T>( \
        const int dim, const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs)

#define INSTANTIATE_ENUMERATE(T) \
    template void sccd::device::enumerate<T>(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx)

#define INSTANTIATE_SORT_ALONG_AXIS(T, I)                                             \
    template void sccd::device::sort_along_axis<T, I>(const int dim,                  \
                                                      const ptrdiff_t n,              \
                                                      const int sort_axis,            \
                                                      T** const SCCD_RESTRICT arrays, \
                                                      I* const SCCD_RESTRICT idx,     \
                                                      T* const SCCD_RESTRICT scratch)

#define INSTANTIATE_COUNT_SELF_OVERLAPS(NXE, T, I)                                               \
    template void sccd::device::count_self_overlaps<NXE, T, I>(const int sort_axis,              \
                                                               const ptrdiff_t element_count,    \
                                                               T** const SCCD_RESTRICT aabbs,    \
                                                               I* const SCCD_RESTRICT idx,       \
                                                               const ptrdiff_t stride,           \
                                                               I** const SCCD_RESTRICT elements, \
                                                               ptrdiff_t* const SCCD_RESTRICT ccdptr)

#define INSTANTIATE_COLLECT_SELF_OVERLAPS(NXE, T, I)                                                          \
    template void sccd::device::collect_self_overlaps<NXE, T, I>(const int sort_axis,                         \
                                                                 const ptrdiff_t element_count,               \
                                                                 T** const SCCD_RESTRICT aabbs,               \
                                                                 I* const SCCD_RESTRICT idx,                  \
                                                                 const ptrdiff_t stride,                      \
                                                                 I** const SCCD_RESTRICT elements,            \
                                                                 const ptrdiff_t* const SCCD_RESTRICT ccdptr, \
                                                                 I* SCCD_RESTRICT foverlap,                   \
                                                                 I* SCCD_RESTRICT noverlap)

#define INSTANTIATE_COUNT_OVERLAPS(FIRST_NXE, SECOND_NXE, T, I)                                                      \
    template void sccd::device::count_overlaps<FIRST_NXE, SECOND_NXE, T, I>(const int sort_axis,                     \
                                                                            const ptrdiff_t first_count,             \
                                                                            T** const SCCD_RESTRICT first_aabbs,     \
                                                                            I* const SCCD_RESTRICT first_idx,        \
                                                                            const ptrdiff_t first_stride,            \
                                                                            I** const SCCD_RESTRICT first_elements,  \
                                                                            const ptrdiff_t second_count,            \
                                                                            T** const SCCD_RESTRICT second_aabbs,    \
                                                                            I* const SCCD_RESTRICT second_idx,       \
                                                                            const ptrdiff_t second_stride,           \
                                                                            I** const SCCD_RESTRICT second_elements, \
                                                                            const T* const SCCD_RESTRICT cummax,     \
                                                                            ptrdiff_t* const SCCD_RESTRICT ccdptr)

INSTANTIATE_CHOOSE_AXIS(float);
INSTANTIATE_CHOOSE_AXIS(double);
INSTANTIATE_ENUMERATE(int32_t);
INSTANTIATE_ENUMERATE(int64_t);
INSTANTIATE_SORT_ALONG_AXIS(float, int32_t);
INSTANTIATE_SORT_ALONG_AXIS(float, int64_t);
INSTANTIATE_SORT_ALONG_AXIS(double, int32_t);
INSTANTIATE_SORT_ALONG_AXIS(double, int64_t);

INSTANTIATE_COUNT_SELF_OVERLAPS(2, float, int32_t);
INSTANTIATE_COUNT_SELF_OVERLAPS(2, float, int64_t);
INSTANTIATE_COUNT_SELF_OVERLAPS(2, double, int32_t);
INSTANTIATE_COUNT_SELF_OVERLAPS(2, double, int64_t);

INSTANTIATE_COLLECT_SELF_OVERLAPS(2, float, int32_t);
INSTANTIATE_COLLECT_SELF_OVERLAPS(2, float, int64_t);
INSTANTIATE_COLLECT_SELF_OVERLAPS(2, double, int32_t);
INSTANTIATE_COLLECT_SELF_OVERLAPS(2, double, int64_t);

INSTANTIATE_COUNT_OVERLAPS(3, 1, float, int32_t);
INSTANTIATE_COUNT_OVERLAPS(3, 1, float, int64_t);
INSTANTIATE_COUNT_OVERLAPS(3, 1, double, int32_t);
INSTANTIATE_COUNT_OVERLAPS(3, 1, double, int64_t);

#undef INSTANTIATE_CHOOSE_AXIS
#undef INSTANTIATE_ENUMERATE
#undef INSTANTIATE_SORT_ALONG_AXIS
#undef INSTANTIATE_COUNT_SELF_OVERLAPS
#undef INSTANTIATE_COLLECT_SELF_OVERLAPS
// Clean-up kernel macros
#undef SCCD_N_WARPS_PER_BLOCK
#undef SCCD_WARP_SIZE
#undef SCCD_WARP_FULL_MASK
