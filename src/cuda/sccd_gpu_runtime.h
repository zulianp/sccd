#ifndef SCCD_GPU_RUNTIME_H
#define SCCD_GPU_RUNTIME_H

// The GPU runtime the device path is built against: CUDA, or HIP for AMD GPUs.
//
// SCCD's device code is written against the CUDA runtime API and compiles
// unchanged under HIP, whose runtime mirrors it name for name. Under
// SCCD_ENABLE_HIP this header includes HIP and maps the cuda* names SCCD uses
// onto their hip* counterparts; otherwise it includes the CUDA runtime. It is
// the one place that knows which runtime is in use, so every translation unit
// that calls the runtime -- kernels, host code, tests, benchmarks -- includes it
// instead of <cuda_runtime.h>. Only the names SCCD calls are mapped.

#if defined(SCCD_ENABLE_HIP)

#include <hip/hip_runtime.h>

#define cudaError_t hipError_t
#define cudaSuccess hipSuccess
#define cudaGetErrorString hipGetErrorString
#define cudaGetLastError hipGetLastError

#define cudaMalloc hipMalloc
#define cudaFree hipFree
#define cudaHostAlloc hipHostMalloc
#define cudaHostAllocDefault hipHostMallocDefault
#define cudaFreeHost hipHostFree
#define cudaMemGetInfo hipMemGetInfo

#define cudaMemcpy hipMemcpy
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemcpyToSymbol hipMemcpyToSymbol
#define cudaMemcpyFromSymbol hipMemcpyFromSymbol
#define cudaMemset hipMemset
#define cudaMemsetAsync hipMemsetAsync

#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaGetDevice hipGetDevice
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaDeviceProp hipDeviceProp_t
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaOccupancyMaxActiveBlocksPerMultiprocessor hipOccupancyMaxActiveBlocksPerMultiprocessor

#define cudaPointerAttributes hipPointerAttribute_t
#define cudaPointerGetAttributes hipPointerGetAttributes
#define cudaMemoryTypeDevice hipMemoryTypeDevice

#define cudaEvent_t hipEvent_t
#define cudaEventCreate hipEventCreate
#define cudaEventDestroy hipEventDestroy
#define cudaEventRecord hipEventRecord
#define cudaEventSynchronize hipEventSynchronize
#define cudaEventElapsedTime hipEventElapsedTime

#else

#include <cuda_runtime.h>

#endif

#endif  // SCCD_GPU_RUNTIME_H
