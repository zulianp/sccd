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
                                                                         const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                                                  const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                                                  aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                                  const BoxRounding rounding) {
            ptrdiff_t e = blockIdx.x * blockDim.x + threadIdx.x;
            if (e >= n_elements) return;

            // Reduce over the element's vertices in registers and write each
            // bound once. Accumulating through aabbs[d][e] instead costs two
            // global reads and two global writes per vertex per dimension --
            // about thirty memory operations for a triangle where six writes
            // do -- and the row pointers themselves live in device memory, so
            // aabbs[d] is another load every time it is named.
            //
            // Rounding outward once at the end gives the same box as rounding
            // every vertex: dnextafter_down and dnextafter_up are monotone, so
            // the minimum of the rounded values is the rounded minimum. The box
            // is therefore identical, not merely close, which is what the
            // conservativeness argument needs.
            for (int d = 0; d < SCCD_DIM; d++) {
                const geom_t* const SCCD_RESTRICT p0d = points0[d];
                const geom_t* const SCCD_RESTRICT p1d = points1[d];

                const idx_t i0 = elements[0][e];
                aabb_t e_min = aabb_t(dmin<geom_t>(p0d[i0], p1d[i0]));
                aabb_t e_max = aabb_t(dmax<geom_t>(p0d[i0], p1d[i0]));

                for (int v = 1; v < nxe; v++) {
                    const idx_t ii = elements[v][e];
                    const geom_t q0 = p0d[ii];
                    const geom_t q1 = p1d[ii];
                    e_min = dmin<aabb_t>(e_min, aabb_t(dmin<geom_t>(q0, q1)));
                    e_max = dmax<aabb_t>(e_max, aabb_t(dmax<geom_t>(q0, q1)));
                }

                if (rounding == BoxRounding::OutwardUlp) {
                    e_min = dnextafter_down(e_min);
                    e_max = dnextafter_up(e_max);
                }

                aabbs[d][e] = e_min;
                aabbs[SCCD_DIM + d][e] = e_max;
            }
        }

        template <typename idx_t, typename geom_t, typename aabb_t>
        void compute_aabbs(const int nxe,
                           const ptrdiff_t n_elements,
                           const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs) {
            compute_aabbs(nxe, n_elements, elements, points0, points1, aabbs, sccd::BoxRounding::Exact);
        }

        template <typename idx_t, typename geom_t, typename aabb_t>
        void compute_aabbs(const int nxe,
                           const ptrdiff_t n_elements,
                           const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                           const BoxRounding rounding) {
            dim3 block(128);
            dim3 grid((n_elements + block.x - 1) / block.x);

            compute_aabbs_face_kernel<idx_t, geom_t, aabb_t>
                <<<grid, block>>>(nxe, n_elements, elements, points0, points1, aabbs, rounding);
            cudaError_t error = cudaGetLastError();

            if (error != cudaSuccess) {
                fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(error));
                exit(1);
            }
        }

        //

        template <typename geom_t, typename aabb_t>
        __global__ void compute_aabbs_node_kernel(                                                  const ptrdiff_t n_nodes,
                                                  const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                                                  const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                                                  aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                                                  const BoxRounding rounding) {
            ptrdiff_t n = blockIdx.x * blockDim.x + threadIdx.x;
            if (n >= n_nodes) return;

            for (int d = 0; d < SCCD_DIM; d++) {
                const geom_t p0 = points0[d][n];
                const geom_t p1 = points1[d][n];
                aabb_t p_min = aabb_t(dmin<geom_t>(p0, p1));
                aabb_t p_max = aabb_t(dmax<geom_t>(p0, p1));
                if (rounding == BoxRounding::OutwardUlp) {
                    p_min = dnextafter_down(p_min);
                    p_max = dnextafter_up(p_max);
                }
                aabbs[d][n] = p_min;
                aabbs[SCCD_DIM + d][n] = p_max;
            }
        }

        template <typename geom_t, typename aabb_t>
        void compute_aabbs(                           const ptrdiff_t n_nodes,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs) {
            compute_aabbs(n_nodes, points0, points1, aabbs, sccd::BoxRounding::Exact);
        }

        template <typename geom_t, typename aabb_t>
        void compute_aabbs(                           const ptrdiff_t n_nodes,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0,
                           const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1,
                           aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,
                           const BoxRounding rounding) {
            dim3 block(128);
            dim3 grid((n_nodes + block.x - 1) / block.x);

            compute_aabbs_node_kernel<geom_t, aabb_t>
                <<<grid, block>>>(n_nodes, points0, points1, aabbs, rounding);
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
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1, \
        aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);        \
    template void sccd::device::compute_aabbs<idx_t, geom_t, aabb_t>(   \
        const int nxe,                                                  \
        const ptrdiff_t n_elements,                                     \
        const idx_t* const SCCD_RESTRICT* const SCCD_RESTRICT elements, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1, \
        aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,         \
        const BoxRounding rounding);                                       \
    template void sccd::device::compute_aabbs<geom_t, aabb_t>(          \
        const ptrdiff_t n_nodes,                                        \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1, \
        aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs);        \
    template void sccd::device::compute_aabbs<geom_t, aabb_t>(          \
        const ptrdiff_t n_nodes,                                        \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points0, \
        const geom_t* const SCCD_RESTRICT* const SCCD_RESTRICT points1, \
        aabb_t* const SCCD_RESTRICT* const SCCD_RESTRICT aabbs,         \
        const BoxRounding rounding)

// Only the diagonals, geom_t == aabb_t.
//
// The two are deliberately separate parameters: boxes narrower than the
// geometry are a supported configuration, and `rounding` exists for exactly
// that case -- it rounds each bound outward by one ULP so a narrower box cannot
// exclude a candidate the wider geometry would have found. Nothing in the
// repository instantiates a mixed pair today, so none is listed; a caller who
// wants float boxes over double geometry adds
// SCCD_VAABB_INSTANTIATE(int, double, float) here and gets a link error until
// they do, rather than a silently wrong answer.
SCCD_VAABB_INSTANTIATE(int, float, float);
SCCD_VAABB_INSTANTIATE(int, double, double);

#undef SCCD_VAABB_INSTANTIATE
