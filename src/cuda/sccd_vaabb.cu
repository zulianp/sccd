#include "sccd_vaabb.cuh"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace sccd {
    namespace device {
        template <typename T>
        static __device__ __forceinline__ T dmin(const T a, const T b) {
            return b < a ? b : a;
        }

        template <typename T>
        static __device__ __forceinline__ T dmax(const T a, const T b) {
            return a < b ? b : a;
        }

        static __device__ __forceinline__ float dnextafter_down(const float x) {
            return nextafterf(x, -FLT_MAX);
        }

        static __device__ __forceinline__ float dnextafter_up(const float x) {
            return nextafterf(x, FLT_MAX);
        }

        static __device__ __forceinline__ double dnextafter_down(const double x) {
            return nextafter(x, -DBL_MAX);
        }

        static __device__ __forceinline__ double dnextafter_up(const double x) {
            return nextafter(x, DBL_MAX);
        }

        template <typename idx_t, typename geom_t, typename aabb_t>
        __global__ void compute_aabbs_face_kernel(const int nxe,
                                                  const ptrdiff_t n_elements,
                                                  const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                                                  const int dim,
                                                  const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                                                  const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                                                  aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                                  const bool safe_inflate) {
            ptrdiff_t e = blockIdx.x * blockDim.x + threadIdx.x;
            if (e >= n_elements) return;

            for (int d = 0; d < dim; d++) {
                const idx_t ii = elements[0][e];
                const geom_t p0 = points0[d][ii];
                const geom_t p1 = points1[d][ii];

                if (safe_inflate) {
                    aabbs[d][e] = dnextafter_down(dmin<geom_t>(p0, p1));
                    aabbs[dim + d][e] = dnextafter_up(dmax<geom_t>(p0, p1));
                } else {
                    aabbs[d][e] = dmin<geom_t>(p0, p1);
                    aabbs[dim + d][e] = dmax<geom_t>(p0, p1);
                }

                for (int v = 1; v < nxe; v++) {
                    const idx_t ii = elements[v][e];
                    const geom_t p0 = points0[d][ii];
                    const geom_t p1 = points1[d][ii];
                    const aabb_t p_min =
                        safe_inflate ? dnextafter_down(dmin<geom_t>(p0, p1)) : dmin<geom_t>(p0, p1);
                    const aabb_t p_max =
                        safe_inflate ? dnextafter_up(dmax<geom_t>(p0, p1)) : dmax<geom_t>(p0, p1);
                    aabbs[d][e] = dmin<aabb_t>(aabbs[d][e], p_min);
                    aabbs[dim + d][e] = dmax<aabb_t>(aabbs[dim + d][e], p_max);
                }
            }
        }

        template <typename idx_t, typename geom_t, typename aabb_t>
        void compute_aabbs(const int nxe,
                           const ptrdiff_t n_elements,
                           const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                           const int dim,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs) {
            compute_aabbs(nxe, n_elements, elements, dim, points0, points1, aabbs, false);
        }

        template <typename idx_t, typename geom_t, typename aabb_t>
        void compute_aabbs(const int nxe,
                           const ptrdiff_t n_elements,
                           const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                           const int dim,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                           const bool safe_inflate) {
            dim3 block(128);
            dim3 grid((n_elements + block.x - 1) / block.x);

            compute_aabbs_face_kernel<idx_t, geom_t, aabb_t>
                <<<grid, block>>>(nxe, n_elements, elements, dim, points0, points1, aabbs, safe_inflate);
            cudaError_t error = cudaGetLastError();

            if (error != cudaSuccess) {
                fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
                exit(1);
            }
        }

        //

        template <typename geom_t, typename aabb_t>
        __global__ void compute_aabbs_node_kernel(const int dim,
                                                  const ptrdiff_t n_nodes,
                                                  const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                                                  const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                                                  aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                                  const bool safe_inflate) {
            ptrdiff_t n = blockIdx.x * blockDim.x + threadIdx.x;
            if (n >= n_nodes) return;

            for (int d = 0; d < dim; d++) {
                const geom_t p0 = points0[d][n];
                const geom_t p1 = points1[d][n];
                if (safe_inflate) {
                    aabbs[d][n] = dnextafter_down(dmin<geom_t>(p0, p1));
                    aabbs[dim + d][n] = dnextafter_up(dmax<geom_t>(p0, p1));
                } else {
                    aabbs[d][n] = dmin<geom_t>(p0, p1);
                    aabbs[dim + d][n] = dmax<geom_t>(p0, p1);
                }
            }
        }

        template <typename geom_t, typename aabb_t>
        void compute_aabbs(const int dim,
                           const ptrdiff_t n_nodes,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs) {
            compute_aabbs(dim, n_nodes, points0, points1, aabbs, false);
        }

        template <typename geom_t, typename aabb_t>
        void compute_aabbs(const int dim,
                           const ptrdiff_t n_nodes,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                           const bool safe_inflate) {
            dim3 block(128);
            dim3 grid((n_nodes + block.x - 1) / block.x);

            compute_aabbs_node_kernel<geom_t, aabb_t>
                <<<grid, block>>>(dim, n_nodes, points0, points1, aabbs, safe_inflate);
            cudaError_t error = cudaGetLastError();

            if (error != cudaSuccess) {
                fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
                exit(1);
            }
        }
    }  // namespace device
}  // namespace sccd

#define SCCD_VAABB_INSTANTIATE(idx_t, geom_t, aabb_t)                    \
    template void sccd::device::compute_aabbs<idx_t, geom_t, aabb_t>(   \
        const int nxe,                                                  \
        const ptrdiff_t n_elements,                                     \
        const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements, \
        const int dim,                                                  \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1, \
        aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);        \
    template void sccd::device::compute_aabbs<idx_t, geom_t, aabb_t>(   \
        const int nxe,                                                  \
        const ptrdiff_t n_elements,                                     \
        const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements, \
        const int dim,                                                  \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1, \
        aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,         \
        const bool safe_inflate);                                       \
    template void sccd::device::compute_aabbs<geom_t, aabb_t>(          \
        const int dim,                                                  \
        const ptrdiff_t n_nodes,                                        \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1, \
        aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);        \
    template void sccd::device::compute_aabbs<geom_t, aabb_t>(          \
        const int dim,                                                  \
        const ptrdiff_t n_nodes,                                        \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1, \
        aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,         \
        const bool safe_inflate)

SCCD_VAABB_INSTANTIATE(int, float, float);
SCCD_VAABB_INSTANTIATE(int, double, double);

#undef SCCD_VAABB_INSTANTIATE
