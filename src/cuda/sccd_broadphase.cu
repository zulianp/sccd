#include "sccd_broadphase.cuh"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <cub/device/device_radix_sort.cuh>

#include "sccd_base.hpp"

#define SCCD_N_WARPS_PER_BLOCK 8
#define SCCD_WARP_SIZE 32
#define SCCD_WARP_FULL_MASK 0xffffffff

namespace sccd {
    namespace device {

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

            auto check_cuda = [](const cudaError_t error) {
                if (error != cudaSuccess) {
                    fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
                    exit(1);
                }
            };

            if (n == 1) {
                enumerate<I>(0, n, idx);
                check_cuda(cudaGetLastError());
                return;
            }

            T* host_arrays[6] = {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};

            if (is_ptr_device(arrays)) {
                check_cuda(cudaMemcpy(host_arrays, arrays, 2 * dim * sizeof(T*), cudaMemcpyDeviceToHost));
            } else {
                for (int d = 0; d < 2 * dim; d++) {
                    host_arrays[d] = arrays[d];
                }
            }

            const T* const SCCD_RESTRICT x = host_arrays[sort_axis];

            I* tmp_idx = nullptr;
            void* tmp_storage = nullptr;
            size_t tmp_storage_bytes = 0;

            check_cuda(cudaMalloc(&tmp_idx, n * sizeof(I)));
            enumerate<I>(0, n, tmp_idx);
            check_cuda(cudaGetLastError());

            check_cuda(cub::DeviceRadixSort::SortPairs(nullptr, tmp_storage_bytes, x, scratch, tmp_idx, idx, n));
            check_cuda(cudaMalloc(&tmp_storage, tmp_storage_bytes));
            check_cuda(cub::DeviceRadixSort::SortPairs(tmp_storage, tmp_storage_bytes, x, scratch, tmp_idx, idx, n));

            check_cuda(cudaMemcpy(host_arrays[sort_axis], scratch, n * sizeof(T), cudaMemcpyDeviceToDevice));

            dim3 block(SCCD_N_WARPS_PER_BLOCK * SCCD_WARP_SIZE);
            dim3 grid((n + block.x - 1) / block.x);
            for (int d = 0; d < 2 * dim; d++) {
                if (d == sort_axis) continue;
                permute_kernel<T, I><<<grid, block>>>(n, idx, host_arrays[d], scratch);
                check_cuda(cudaGetLastError());
                check_cuda(cudaMemcpy(host_arrays[d], scratch, n * sizeof(T), cudaMemcpyDeviceToDevice));
            }

            check_cuda(cudaFree(tmp_storage));
            check_cuda(cudaFree(tmp_idx));
        }

    }  // namespace device
}  // namespace sccd

#define INSTANTIATE_CHOOSE_AXIS(T)             \
    template int sccd::device::choose_axis<T>( \
        const int dim, const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);

#define INSTANTIATE_ENUMERATE(T) \
    template void sccd::device::enumerate<T>(const ptrdiff_t begin, const ptrdiff_t end, T* const SCCD_RESTRICT idx);

#define INSTANTIATE_SORT_ALONG_AXIS(T, I)                                             \
    template void sccd::device::sort_along_axis<T, I>(const int dim,                  \
                                                      const ptrdiff_t n,              \
                                                      const int sort_axis,            \
                                                      T** const SCCD_RESTRICT arrays, \
                                                      I* const SCCD_RESTRICT idx,     \
                                                      T* const SCCD_RESTRICT scratch);

INSTANTIATE_CHOOSE_AXIS(float);
INSTANTIATE_CHOOSE_AXIS(double);
INSTANTIATE_ENUMERATE(int32_t);
INSTANTIATE_ENUMERATE(int64_t);
INSTANTIATE_SORT_ALONG_AXIS(float, int32_t);
INSTANTIATE_SORT_ALONG_AXIS(float, int64_t);
INSTANTIATE_SORT_ALONG_AXIS(double, int32_t);
INSTANTIATE_SORT_ALONG_AXIS(double, int64_t);

#undef INSTANTIATE_CHOOSE_AXIS
#undef INSTANTIATE_ENUMERATE
#undef INSTANTIATE_SORT_ALONG_AXIS

// Clean-up kernel macros
#undef SCCD_N_WARPS_PER_BLOCK
#undef SCCD_WARP_SIZE
#undef SCCD_WARP_FULL_MASK
