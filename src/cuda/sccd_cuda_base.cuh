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
    }  // namespace device
}  // namespace sccd

#define SCCD_CHECK_CUDA(error) \
    cudaDeviceSynchronize();   \
    cuda_check(error)

#define SCCD_CUDA_LAST_ERROR() SCCD_CHECK_CUDA(cudaGetLastError())

#endif  // SCCD_CUDA_BASE_CUH
