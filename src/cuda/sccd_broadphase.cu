#include "sccd_broadphase.cuh"

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include "sccd_base.hpp"

//  template <typename T>
//     static int choose_axis(const ptrdiff_t n, T** const SCCD_RESTRICT aabb) {
//     T mean[3] = {0};
//     T var[3] = {0};
//     for (int d = 0; d < 3; d++) {
//         for (ptrdiff_t i = 0; i < n; i++) {
//             const T c = (aabb[d + 3][i] + aabb[d][i]) / 2;
//             mean[d] += c;
//         }

//         mean[d] /= n;
//         for (ptrdiff_t i = 0; i < n; i++) {
//             const T c = (aabb[d + 3][i] + aabb[d][i]) / 2;
//             var[d] += (c - mean[d]) * (c - mean[d]);
//         }
//     }

//     int fargmax = 0;
//     T fmax = var[0];

//     for (int d = 1; d < 3; d++) {
//         if (fmax < var[d]) {
//             fmax = var[d];
//             fargmax = d;
//         }
//     }

//     return fargmax;
// }

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
        __device__ void t_warp_reduce(const T val, T* block_accumulator, T* result) {
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
                                                T* const SCCD_RESTRICT* const SCCD_RESTRICT mean) {
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
                t_warp_reduce(local_mean[d], block_accumulator, mean[d]);
            }
        }

        template <typename T>
        __global__ void choose_axis_var_kernel(const int dim,
                                               const ptrdiff_t n,
                                               const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                               T* const SCCD_RESTRICT mean,
                                               T* const SCCD_RESTRICT* const SCCD_RESTRICT var) {
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
                t_warp_reduce(local_var[d], block_accumulator, var[d]);
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
    }  // namespace device
}  // namespace sccd

#define INSTANTIATE_CHOOSE_AXIS(T)             \
    template int sccd::device::choose_axis<T>( \
        const int dim, const ptrdiff_t n, const T* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);

INSTANTIATE_CHOOSE_AXIS(float);
INSTANTIATE_CHOOSE_AXIS(double);

#undef INSTANTIATE_CHOOSE_AXIS