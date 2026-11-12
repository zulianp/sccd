#ifndef SCCD_GPU_CUB_CUH
#define SCCD_GPU_CUB_CUH

// The device-wide primitives the broad phase uses: CUB on CUDA, hipCUB on HIP.
// hipCUB mirrors CUB's interface in its own namespace, so under HIP `cub`
// names it and the call sites are the same on both.

#if defined(SCCD_ENABLE_HIP)

#include <hipcub/hipcub.hpp>

namespace cub = hipcub;

#else

#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_scan.cuh>

#endif

#endif  // SCCD_GPU_CUB_CUH
