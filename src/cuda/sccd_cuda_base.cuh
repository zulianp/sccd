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

        static void cuda_check(const cudaError_t error, const char* file, const int line) {
            if (error != cudaSuccess) {
                fprintf(stderr, "CUDA error: %s in %s:%d\n", cudaGetErrorString(error), file, line);
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

        // CAS-based atomic-min for float/double, works on both global and shared addresses.
        static inline __device__ float atomic_min(float* address, float val) {
            int* a = reinterpret_cast<int*>(address);
            int old = __float_as_int(*address), assumed;
            do {
                if (__int_as_float(old) <= val) return __int_as_float(old);
                assumed = old;
                old = atomicCAS(a, assumed, __float_as_int(val));
            } while (assumed != old);
            return __int_as_float(old);
        }

        static inline __device__ double atomic_min(double* address, double val) {
            unsigned long long* a = reinterpret_cast<unsigned long long*>(address);
            unsigned long long old = __double_as_longlong(*address), assumed;
            do {
                if (__longlong_as_double(old) <= val) return __longlong_as_double(old);
                assumed = old;
                old = atomicCAS(a, assumed, __double_as_longlong(val));
            } while (assumed != old);
            return __longlong_as_double(old);
        }

    }  // namespace device
}  // namespace sccd

#define SCCD_CHECK_CUDA(error) \
    cudaDeviceSynchronize();   \
    cuda_check(error, __FILE__, __LINE__)

#define SCCD_CUDA_LAST_ERROR() SCCD_CHECK_CUDA(cudaGetLastError())

#endif  // SCCD_CUDA_BASE_CUH
