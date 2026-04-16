#include "sccd_narrowphase.cuh"

#include <stdint.h>

#include "sccd_cuda_base.cuh"
#include "sccd_reduce.cuh"

namespace sccd {
    namespace device {

        template <typename T, typename Vec4>
        static inline __device__ void compute_edge_edge_tolerance_soa(const T codomain_tol,
                                                                      const Vec4 sx,
                                                                      const Vec4 sy,
                                                                      const Vec4 sz,
                                                                      const Vec4 ex,
                                                                      const Vec4 ey,
                                                                      const Vec4 ez,
                                                                      T* const SCCD_RESTRICT tol0,
                                                                      T* const SCCD_RESTRICT tol1,
                                                                      T* const SCCD_RESTRICT tol2) {
            const T v0sx = sx.x;
            const T v0sy = sy.x;
            const T v0sz = sz.x;
            const T v1sx = sx.y;
            const T v1sy = sy.y;
            const T v1sz = sz.y;
            const T v2sx = sx.z;
            const T v2sy = sy.z;
            const T v2sz = sz.z;
            const T v3sx = sx.w;
            const T v3sy = sy.w;
            const T v3sz = sz.w;
            const T v0ex = ex.x;
            const T v0ey = ey.x;
            const T v0ez = ez.x;
            const T v1ex = ex.y;
            const T v1ey = ey.y;
            const T v1ez = ez.y;
            const T v2ex = ex.z;
            const T v2ey = ey.z;
            const T v2ez = ez.z;
            const T v3ex = ex.w;
            const T v3ey = ey.w;
            const T v3ez = ez.w;

            const T ssa0 = v0ex - v0sx;
            const T ssa1 = -v2ex + v2sx;
            const T ssa2 = -v3ex + v3sx;
            const T ssa3 = v0ey - v0sy;
            const T ssa4 = -v2ey + v2sy;
            const T ssa5 = -v3ey + v3sy;
            const T ssa6 = v0ez - v0sz;
            const T ssa7 = -v2ez + v2sz;
            const T ssa8 = -v3ez + v3sz;
            const T ssa9 = -v1sx;
            const T ssa10 = ssa9 + v1ex;
            const T ssa11 = -v1sy;
            const T ssa12 = ssa11 + v1ey;
            const T ssa13 = -v1sz;
            const T ssa14 = ssa13 + v1ez;
            const T ssa15 = (1.0 / 3.0) * codomain_tol;
            const T ssa16 =
                ssa15 / device::max<T>(
                            device::abs<T>(ssa0 + ssa1),
                            device::max<T>(
                                device::abs<T>(ssa0 + ssa2),
                                device::max<T>(
                                    device::abs<T>(ssa1 + ssa10),
                                    device::max<T>(
                                        device::abs<T>(ssa10 + ssa2),
                                        device::max<T>(
                                            device::abs<T>(ssa12 + ssa4),
                                            device::max<T>(
                                                device::abs<T>(ssa12 + ssa5),
                                                device::max<T>(
                                                    device::abs<T>(ssa14 + ssa7),
                                                    device::max<T>(
                                                        device::abs<T>(ssa14 + ssa8),
                                                        device::max<T>(
                                                            device::abs<T>(ssa3 + ssa4),
                                                            device::max<T>(
                                                                device::abs<T>(ssa3 + ssa5),
                                                                device::max<T>(device::abs<T>(ssa6 + ssa7),
                                                                               device::abs<T>(ssa6 + ssa8))))))))))));
            *tol0 = ssa16;
            *tol1 = ssa16;
            *tol2 = ssa15 /
                    device::max<T>(
                        device::abs<T>(ssa11 + v0sy),
                        device::max<T>(device::abs<T>(ssa13 + v0sz),
                                       device::max<T>(device::abs<T>(ssa9 + v0sx),
                                                      device::max<T>(device::abs<T>(v0ex - v1ex),
                                                                     device::max<T>(device::abs<T>(v0ey - v1ey),
                                                                                    device::abs<T>(v0ez - v1ez))))));
        }

        template <typename T, typename Vec4>
        __device__ void sample_f_ee(const T tl,
                                    const T tu,
                                    const T ul,
                                    const T uu,
                                    const T vl,
                                    const T vu,
                                    const Vec4 sx,
                                    const Vec4 ex,
                                    T* const f) {
            // Compute spatial displacements range
            Vec4 dx = ex - sx;

            // Compute temporal displacements for lower bound
            {
                Vec4 xt = tl * dx + sx;
                f[0] = ((xt.y - xt.x) * ul + xt.x - (xt.w - xt.y) * vl + xt.y);
                f[1] = ((xt.y - xt.x) * ul + xt.x - (xt.w - xt.y) * vu + xt.y);
                f[2] = ((xt.y - xt.x) * uu + xt.x - (xt.w - xt.y) * vl + xt.y);
                f[3] = ((xt.y - xt.x) * uu + xt.x - (xt.w - xt.y) * vu + xt.y);
            }

            // Compute temporal displacements for upper bound
            {
                Vec4 xt = tu * dx + sx;
                f[4] = ((xt.y - xt.x) * ul + xt.x - (xt.w - xt.y) * vl + xt.y);
                f[5] = ((xt.y - xt.x) * ul + xt.x - (xt.w - xt.y) * vu + xt.y);
                f[6] = ((xt.y - xt.x) * uu + xt.x - (xt.w - xt.y) * vl + xt.y);
                f[7] = ((xt.y - xt.x) * uu + xt.x - (xt.w - xt.y) * vu + xt.y);
            }
        }

        template <typename T>
        __device__ void fminmax(const T* const SCCD_RESTRICT f, T& fmin, T& fmax) {
            fmin = f[0];
            fmax = f[0];
            for (int i = 1; i < 8; i++) {
                fmin = MIN(fmin, f[i]);
                fmax = MAX(fmax, f[i]);
            }
        }

#define MASK_FULL 0xf
#define MASK_DOMAIN_SMALLER_THAN_TOL 1
#define MASK_BOX_INSIDE_EPSILON_BOX 2
#define MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE 4
#define MASK_INTERVAL_TERMINAL 8

        template <typename T>
        __device__ uint8_t cond_mask(const T fmin, const T fmax, const T tol, const T adaptive_tol) {
            bool cond1 = (fmax - fmin <= adaptive_tol);
            bool cond2 = !((fmin < tol) | (fmax > -tol));
            bool cond3 = (fmax - fmin < tol);
            bool cond4 = (fmin >= fmax);

            uint8_t cond_mask = (cond1 ? MASK_DOMAIN_SMALLER_THAN_TOL : 0);
            cond_mask |= (cond2 ? MASK_BOX_INSIDE_EPSILON_BOX : 0);
            cond_mask |= (cond3 ? MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE : 0);
            cond_mask |= (cond4 ? MASK_INTERVAL_TERMINAL : 0);
            cond_mask &= ((fmin <= tol) & ((fmax >= -tol) ? MASK_FULL : 0));
            return cond_mask;
        }

        template <typename T, typename Vec4>
        __device__ void contains_origin_ee(const int nx,
                                           const int ny,
                                           const int x,
                                           const int y,
                                           const T tlower,
                                           const T tupper,
                                           const T ulower,
                                           const T uupper,
                                           const T vlower,
                                           const T vupper,
                                           const Vec4 sx,
                                           const Vec4 sy,
                                           const Vec4 sz,
                                           const Vec4 ex,
                                           const Vec4 ey,
                                           const Vec4 ez,
                                           const T tol,
                                           const T* const SCCD_RESTRICT adaptive_tol,
                                           int* const SCCD_RESTRICT contains_origin,
                                           int* const SCCD_RESTRICT accept) {
            const T u_h = (uupper - ulower) / nx;
            const T v_h = (vupper - vlower) / ny;

            const T ul = ulower + x * u_h;
            const T vl = vlower + y * v_h;
            const T uu = ulower + (x + 1) * u_h;
            const T vu = vlower + (y + 1) * v_h;

            T fmin, fmax;
            T f[8];

            sample_f_ee(tlower, tupper, ul, uu, vl, vu, sx, ex, f);
            fminmax(f, fmin, fmax);
            const uint8_t x_mask = cond_mask(fmin, fmax, tol, adaptive_tol[0]);
            *contains_origin = (fmin <= tol) & (fmax >= -tol);

            sample_f_ee(tlower, tupper, ul, uu, vl, vu, sy, ey, f);
            fminmax(f, fmin, fmax);
            const uint8_t y_mask = cond_mask(fmin, fmax, tol, adaptive_tol[1]);
            *contains_origin &= (fmin <= tol) & (fmax >= -tol);

            sample_f_ee(tlower, tupper, ul, uu, vl, vu, sz, ez, f);
            fminmax(f, fmin, fmax);
            const uint8_t z_mask = cond_mask(fmin, fmax, tol, adaptive_tol[2]);
            *contains_origin &= (fmin <= tol) & (fmax >= -tol);

            const uint8_t and_mask = (x_mask & y_mask & z_mask);
            const uint8_t or_mask = (x_mask | y_mask | z_mask);

            const bool cond1 = and_mask & MASK_DOMAIN_SMALLER_THAN_TOL;
            const bool cond2 = and_mask & MASK_BOX_INSIDE_EPSILON_BOX;
            const bool cond3 = or_mask & MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE;
            const bool cond4 = and_mask & MASK_INTERVAL_TERMINAL;

            *accept = *contains_origin && (cond1 || cond2 || cond3 || cond4);
            // If contains origin and does not accept, it is a nutcase (should we handle it separately?)
        }

#define SCCD_MAX_STACK_SIZE 32

        template <typename T, typename I>
        __global__ void narrow_phase_ee_kernel(const size_t noverlaps,
                                               const I* const SCCD_RESTRICT e0overalp,
                                               const I* const SCCD_RESTRICT e1overalp,
                                               // Geometric data
                                               T** const SCCD_RESTRICT sp,
                                               T** const SCCD_RESTRICT ep,
                                               const size_t edge_stride,
                                               I** const SCCD_RESTRICT edges,
                                               // Output
                                               const T max_toi,
                                               const T tol,
                                               T* SCCD_RESTRICT toi) {
            __shared__ T stack_tlower[SCCD_MAX_STACK_SIZE];
            __shared__ T stack_tupper[SCCD_MAX_STACK_SIZE];
            __shared__ int stack_level[SCCD_MAX_STACK_SIZE];
            __shared__ T block_accumulator[1024];

            using Vec4 = typename device::Vec4Type<T>::type;

            const int thIdx = threadIdx.x + threadIdx.y * blockDim.x;
            const int blDim = blockDim.x * blockDim.y;
            const bool root = (thIdx == 0);

            int stack_size = 0;
            if (root) {
                stack_tlower[stack_size] = 0;
                stack_tupper[stack_size] = 1;
                stack_level[stack_size] = 0;
                stack_size++;
            }

            const I ea = e0overalp[blockIdx.x];
            const I eb = e1overalp[blockIdx.x];

            Vec4 sx, ex;
            Vec4 sy, ey;
            Vec4 sz, ez;

            const ptrdiff_t idxa0 = edges[0][ea * edge_stride];
            const ptrdiff_t idxa1 = edges[1][ea * edge_stride];

            const ptrdiff_t idxb0 = edges[0][eb * edge_stride];
            const ptrdiff_t idxb1 = edges[1][eb * edge_stride];

            // All x-coordinates
            sx.x = sp[0][idxa0];
            sx.y = sp[0][idxa1];
            sx.z = sp[0][idxb0];
            sx.w = sp[0][idxb1];

            ex.x = ep[0][idxa0];
            ex.y = ep[0][idxa1];
            ex.z = ep[0][idxb0];
            ex.w = ep[0][idxb1];

            // All y-coordinates
            sy.x = sp[1][idxa0];
            sy.y = sp[1][idxa1];
            sy.z = sp[1][idxb0];
            sy.w = sp[1][idxb1];

            ey.x = ep[1][idxa0];
            ey.y = ep[1][idxa1];
            ey.z = ep[1][idxb0];
            ey.w = ep[1][idxb1];

            // All z-coordinates
            sz.x = sp[2][idxa0];
            sz.y = sp[2][idxa1];
            sz.z = sp[2][idxb0];
            sz.w = sp[2][idxb1];

            ez.x = ep[2][idxa0];
            ez.y = ep[2][idxa1];
            ez.z = ep[2][idxb0];
            ez.w = ep[2][idxb1];

            T adaptive_tol[3];
            compute_edge_edge_tolerance_soa<T, Vec4>(
                tol, sx, sy, sz, ex, ey, ez, &adaptive_tol[0], &adaptive_tol[1], &adaptive_tol[2]);

            T current_toi = max_toi;
            while (stack_size < SCCD_MAX_STACK_SIZE && stack_size > 0) {
                __syncthreads();

                --stack_size;  // pop the top of the stack
                T tlower = stack_tlower[stack_size];
                T tupper = stack_tupper[stack_size];
                int level = stack_level[stack_size];

                T tmid = (tlower + tupper) * 0.5;

                int left_contains_origin = 0;
                int left_accept = 0;

                contains_origin_ee<T, Vec4>(blockDim.x,
                                            blockDim.y,
                                            threadIdx.x,
                                            threadIdx.y,
                                            tlower,
                                            tmid,
                                            0,
                                            1,
                                            0,
                                            1,
                                            sx,
                                            sy,
                                            sz,
                                            ex,
                                            ey,
                                            ez,
                                            tol,
                                            adaptive_tol,
                                            &left_contains_origin,
                                            &left_accept);

                if (left_accept) {
                    current_toi = MIN(current_toi, tlower);
                }

                // Reduce left_accept and left_contains_origin
                int or_left_accept = 0;
                device::block_max_to_root(thIdx, blDim, left_accept, (int*)block_accumulator, &or_left_accept);
                or_left_accept = device::broadcast_to_block(0, or_left_accept);

                if (or_left_accept) {
                    // We have found the toi, so we are done
                    break;
                }

                int or_left_contains_origin = 0;
                device::block_max_to_root(
                    thIdx, blDim, left_contains_origin, (int*)block_accumulator, &or_left_contains_origin);

                if (root) {
                    if (or_left_contains_origin) {
                        if (level < SCCD_MAX_STACK_SIZE && tlower < current_toi) {
                            stack_tlower[stack_size] = tlower;
                            stack_tupper[stack_size] = tmid;
                            stack_level[stack_size] = level + 1;
                            stack_size++;
                        } else {
                            // TODO
                        }
                    }
                }

                int right_contains_origin = 0;
                int right_accept = 0;

                contains_origin_ee<T, Vec4>(blockDim.x,
                                            blockDim.y,
                                            threadIdx.x,
                                            threadIdx.y,
                                            tmid,
                                            tupper,
                                            0,
                                            1,
                                            0,
                                            1,
                                            sx,
                                            sy,
                                            sz,
                                            ex,
                                            ey,
                                            ez,
                                            tol,
                                            adaptive_tol,
                                            &right_contains_origin,
                                            &right_accept);

                if (right_accept) {
                    current_toi = MIN(current_toi, tmid);
                }

                // Reduce left_accept and left_contains_origin
                int or_right_accept = 0;
                device::block_max_to_root(thIdx, blDim, right_accept, (int*)block_accumulator, &or_right_accept);
                or_right_accept = device::broadcast_to_block(0, or_right_accept);

                if (or_right_accept) {
                    // We have found the toi, so we are done
                    break;
                }

                int or_right_contains_origin = 0;
                device::block_max_to_root(
                    thIdx, blDim, right_contains_origin, (int*)block_accumulator, &or_right_contains_origin);

                if (root) {
                    if (or_right_contains_origin) {
                        if (level < SCCD_MAX_STACK_SIZE && tmid < current_toi) {
                            stack_tlower[stack_size] = tmid;
                            stack_tupper[stack_size] = tupper;
                            stack_level[stack_size] = level + 1;
                            stack_size++;
                        } else {
                            // TODO
                        }
                    }
                }
            }
        }

        template <typename T, typename I>
        T narrow_phase_ee(const size_t noverlaps,
                          const I* const SCCD_RESTRICT e0overalp,
                          const I* const SCCD_RESTRICT e1overalp,
                          // Geometric data
                          T** const SCCD_RESTRICT v0,
                          T** const SCCD_RESTRICT v1,
                          const size_t edge_stride,
                          I** const SCCD_RESTRICT edges,
                          // Output
                          const T max_toi) {
            SCCD_CUDA_LAST_ERROR();

            dim3 block(32, 32, 1);
            dim3 grid(noverlaps, 1, 1);

            T* d_toi = nullptr;
            cudaMalloc(&d_toi, noverlaps * sizeof(T));

            narrow_phase_ee_kernel<T, I>
                <<<grid, block>>>(noverlaps, e0overalp, e1overalp, v0, v1, edge_stride, edges, max_toi, 1e-8, d_toi);

            SCCD_CUDA_LAST_ERROR();

            T h_toi;
            cudaMemcpy(&h_toi, d_toi, sizeof(T), cudaMemcpyDeviceToHost);
            cudaFree(d_toi);

            SCCD_CUDA_LAST_ERROR();
            return h_toi;
        }

        template <int nxe, typename T, typename I>
        T narrow_phase_vf(const size_t noverlaps,
                          const I* const SCCD_RESTRICT voveralp,
                          const I* const SCCD_RESTRICT foveralp,
                          // Geometric data
                          T** const SCCD_RESTRICT v0,
                          T** const SCCD_RESTRICT v1,
                          const size_t face_stride,
                          I** const SCCD_RESTRICT faces,
                          const T max_toi) {
            return max_toi;
        }

    }  // namespace device
}  // namespace sccd

#define INSTANTIATE_NARROW_PHASE_EE(T, I)                                                  \
    template T sccd::device::narrow_phase_ee<T, I>(const size_t noverlaps,                 \
                                                   const I* const SCCD_RESTRICT e0overalp, \
                                                   const I* const SCCD_RESTRICT e1overalp, \
                                                   T** const SCCD_RESTRICT v0,             \
                                                   T** const SCCD_RESTRICT v1,             \
                                                   const size_t edge_stride,               \
                                                   I** const SCCD_RESTRICT edges,          \
                                                   const T max_toi);

#define INSTANTIATE_NARROW_PHASE_VF(NXE, T, I)                                                 \
    template T sccd::device::narrow_phase_vf<NXE, T, I>(const size_t noverlaps,                \
                                                        const I* const SCCD_RESTRICT voveralp, \
                                                        const I* const SCCD_RESTRICT foveralp, \
                                                        T** const SCCD_RESTRICT v0,            \
                                                        T** const SCCD_RESTRICT v1,            \
                                                        const size_t face_stride,              \
                                                        I** const SCCD_RESTRICT faces,         \
                                                        const T max_toi);

INSTANTIATE_NARROW_PHASE_EE(float, int32_t);
INSTANTIATE_NARROW_PHASE_EE(float, int64_t);
INSTANTIATE_NARROW_PHASE_EE(double, int32_t);
INSTANTIATE_NARROW_PHASE_EE(double, int64_t);

INSTANTIATE_NARROW_PHASE_VF(3, float, int32_t);
INSTANTIATE_NARROW_PHASE_VF(3, float, int64_t);
INSTANTIATE_NARROW_PHASE_VF(3, double, int32_t);
INSTANTIATE_NARROW_PHASE_VF(3, double, int64_t);

#undef INSTANTIATE_NARROW_PHASE_EE
#undef INSTANTIATE_NARROW_PHASE_VF