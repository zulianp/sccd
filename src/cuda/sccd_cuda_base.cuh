#ifndef SCCD_CUDA_BASE_CUH
#define SCCD_CUDA_BASE_CUH

#include "sccd_base.hpp"

#include <cstddef>
#include <cstdio>

#include <cuda_runtime.h>

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#define SCCD_WARP_SIZE 32
#define SCCD_WARP_FULL_MASK 0xffffffff

namespace sccd {
    namespace device {

        static void cuda_check(const cudaError_t error) {
            if (error != cudaSuccess) {
                fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
                exit(1);
            }
        }

        template <typename T>
        static inline __device__ T max(const T a, const T b) {
            return (a > b) ? a : b;
        }

        template <typename T>
        static inline __device__ T min(const T a, const T b) {
            return (a < b) ? a : b;
        }

        template <typename T>
        static inline __device__ T abs(const T x) {
            return x < 0 ? -x : x;
        }

        template <typename T>
        static inline __device__ T pow2(const T x) {
            return x * x;
        }
        template <typename T>
        static inline __device__ T pow3(const T x) {
            return x * x * x;
        }

        template <typename T>
        struct Vec4Type {};

        template <>
        struct Vec4Type<float> {
            using type = float4;
        };

        template <>
        struct Vec4Type<double> {
            using type = double4;
        };

        inline __device__ float4 operator-(const float4& a, const float4& b) {
            return make_float4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
        }

        inline __device__ double4 operator-(const double4& a, const double4& b) {
            return make_double4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
        }

        inline __device__ float4 operator*(const float4& a, const float4& b) {
            return make_float4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
        }

        inline __device__ double4 operator*(const double4& a, const double4& b) {
            return make_double4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
        }

        inline __device__ float4 operator+(const float4& a, const float4& b) {
            return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
        }

        inline __device__ double4 operator+(const double4& a, const double4& b) {
            return make_double4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
        }

        inline __device__ double4 operator/(const double4& a, const double4& b) {
            return make_double4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
        }

        inline __device__ float4 operator/(const float4& a, const float4& b) {
            return make_float4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
        }

        inline __device__ float4 operator*(const float alpha, const float4& b) {
            return make_float4(alpha * b.x, alpha * b.y, alpha * b.z, alpha * b.w);
        }

        inline __device__ double4 operator*(const double alpha, const double4& b) {
            return make_double4(alpha * b.x, alpha * b.y, alpha * b.z, alpha * b.w);
        }

    }  // namespace device
}  // namespace sccd

#define SCCD_CHECK_CUDA(error) \
    cudaDeviceSynchronize();   \
    cuda_check(error)

#define SCCD_CUDA_LAST_ERROR() SCCD_CHECK_CUDA(cudaGetLastError())

#endif  // SCCD_CUDA_BASE_CUH
