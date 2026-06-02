#include "sccd_narrowphase.cuh"

#include <stdint.h>
#include <stdlib.h>

#include <cassert>
#include <climits>
#include <limits>

#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <thrust/functional.h>
#include <thrust/memory.h>
#include <thrust/reduce.h>

#include "sccd_cuda_base.cuh"
#include "sccd_reduce.cuh"

#ifndef SCCD_NP_SHARED_STACK_CAP
#define SCCD_NP_SHARED_STACK_CAP 1024
#endif

#ifndef SCCD_NP_THREADS_PER_BLOCK
#define SCCD_NP_THREADS_PER_BLOCK 128
#endif

#ifndef SCCD_CUDA_ADAPTIVE_SPLIT
#define SCCD_CUDA_ADAPTIVE_SPLIT 1
#endif

namespace sccd {
    namespace device {

        // Axis-aligned subdomain in (t, u, v) parameter space.
        template <typename T>
        struct Domain {
            T tlower;
            T tupper;
            T ulower;
            T uupper;
            T vlower;
            T vupper;
        };

        template <typename T>
        struct Stack {
            T* tlower;
            T* tupper;
            T* ulower;
            T* uupper;
            T* vlower;
            T* vupper;
            int* level;
            int* qid;
            int* top;
            int* request;
            int capacity;
        };

        template <bool is_vf, typename T>
        static inline __device__ bool is_domain_valid(const Domain<T>& domain,
                                                      const T toi,
                                                      const T* const SCCD_RESTRICT tols) {
            if constexpr (is_vf) {
                return domain.tlower < toi && (domain.ulower + domain.vlower < 1 + tols[1] + tols[2]);
            } else {
                return domain.tlower < toi;
            }
        }

        template <typename T, typename Vec4>
        static inline __device__ void compute_edge_edge_tolerance(const T codomain_tol,
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
            const T ssa15 = T(1.0 / 3.0) * codomain_tol;
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
        static inline __device__ void compute_face_vertex_tolerance(const T codomain_tol,
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
            const T ssa1 = -v2ex;
            const T ssa2 = ssa0 + v2sx;
            const T ssa3 = -v3ex;
            const T ssa4 = ssa3 + v3sx;
            const T ssa5 = v0ey - v0sy;
            const T ssa6 = -v2ey;
            const T ssa7 = ssa5 + v2sy;
            const T ssa8 = -v3ey;
            const T ssa9 = ssa8 + v3sy;
            const T ssa10 = v0ez - v0sz;
            const T ssa11 = -v2ez;
            const T ssa12 = ssa10 + v2sz;
            const T ssa13 = -v3ez;
            const T ssa14 = ssa13 + v3sz;
            const T ssa15 = ssa1 + v1ex;
            const T ssa16 = ssa6 + v1ey;
            const T ssa17 = ssa11 + v1ez;
            const T ssa18 = T(1.0 / 3.0) * codomain_tol;
            *tol0 = ssa18 /
                    device::max<T>(
                        device::abs<T>(ssa0 + ssa4),
                        device::max<T>(
                            device::abs<T>(ssa1 + ssa2),
                            device::max<T>(
                                device::abs<T>(ssa10 + ssa14),
                                device::max<T>(
                                    device::abs<T>(ssa11 + ssa12),
                                    device::max<T>(
                                        device::abs<T>(ssa5 + ssa9),
                                        device::max<T>(
                                            device::abs<T>(ssa6 + ssa7),
                                            device::max<T>(
                                                device::abs<T>(ssa0 - v1ex + v1sx),
                                                device::max<T>(
                                                    device::abs<T>(ssa10 - v1ez + v1sz),
                                                    device::max<T>(
                                                        device::abs<T>(ssa5 - v1ey + v1sy),
                                                        device::max<T>(
                                                            device::abs<T>(ssa12 + ssa14 + ssa17 - v1sz),
                                                            device::max<T>(
                                                                device::abs<T>(ssa15 + ssa2 + ssa4 - v1sx),
                                                                device::abs<T>(ssa16 + ssa7 + ssa9 - v1sy))))))))))));
            *tol1 = ssa18 /
                    device::max<T>(
                        device::abs<T>(ssa15),
                        device::max<T>(device::abs<T>(ssa16),
                                       device::max<T>(device::abs<T>(ssa17),
                                                      device::max<T>(device::abs<T>(v1sx - v2sx),
                                                                     device::max<T>(device::abs<T>(v1sy - v2sy),
                                                                                    device::abs<T>(v1sz - v2sz))))));
            *tol2 = ssa18 /
                    device::max<T>(
                        device::abs<T>(ssa13 + v1ez),
                        device::max<T>(device::abs<T>(ssa3 + v1ex),
                                       device::max<T>(device::abs<T>(ssa8 + v1ey),
                                                      device::max<T>(device::abs<T>(v1sx - v3sx),
                                                                     device::max<T>(device::abs<T>(v1sy - v3sy),
                                                                                    device::abs<T>(v1sz - v3sz))))));
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

            {
                Vec4 xt = tl * dx + sx;
                const T pa_l = (xt.y - xt.x) * ul + xt.x;
                const T pa_u = (xt.y - xt.x) * uu + xt.x;
                const T pb_l = (xt.w - xt.z) * vl + xt.z;
                const T pb_u = (xt.w - xt.z) * vu + xt.z;
                f[0] = pa_l - pb_l;
                f[1] = pa_l - pb_u;
                f[2] = pa_u - pb_l;
                f[3] = pa_u - pb_u;
            }

            {
                Vec4 xt = tu * dx + sx;
                const T pa_l = (xt.y - xt.x) * ul + xt.x;
                const T pa_u = (xt.y - xt.x) * uu + xt.x;
                const T pb_l = (xt.w - xt.z) * vl + xt.z;
                const T pb_u = (xt.w - xt.z) * vu + xt.z;
                f[4] = pa_l - pb_l;
                f[5] = pa_l - pb_u;
                f[6] = pa_u - pb_l;
                f[7] = pa_u - pb_u;
            }
        }

        template <typename T, typename Vec4>
        __device__ void sample_f_vf(const T tl,
                                    const T tu,
                                    const T ul,
                                    const T uu,
                                    const T vl,
                                    const T vu,
                                    const Vec4 sx,
                                    const Vec4 ex,
                                    T* const SCCD_RESTRICT F) {
            Vec4 dx = ex - sx;

            {
                Vec4 xt = tl * dx + sx;
                const T vertex = xt.x;

                const T face0 = (xt.z - xt.y) * ul + (xt.w - xt.y) * vl + xt.y;
                const T face1 = (xt.z - xt.y) * ul + (xt.w - xt.y) * vu + xt.y;
                const T face2 = (xt.z - xt.y) * uu + (xt.w - xt.y) * vl + xt.y;
                const T face3 = (xt.z - xt.y) * uu + (xt.w - xt.y) * vu + xt.y;

                F[0] = vertex - face0;
                F[1] = vertex - face1;
                F[2] = vertex - face2;
                F[3] = vertex - face3;
            }

            {
                Vec4 xt = tu * dx + sx;
                const T vertex = xt.x;

                const T face0 = (xt.z - xt.y) * ul + (xt.w - xt.y) * vl + xt.y;
                const T face1 = (xt.z - xt.y) * ul + (xt.w - xt.y) * vu + xt.y;
                const T face2 = (xt.z - xt.y) * uu + (xt.w - xt.y) * vl + xt.y;
                const T face3 = (xt.z - xt.y) * uu + (xt.w - xt.y) * vu + xt.y;

                F[4] = vertex - face0;
                F[5] = vertex - face1;
                F[6] = vertex - face2;
                F[7] = vertex - face3;
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

        template <bool is_vf, typename T, typename Vec4>
        static inline __device__ void sample_f(const T tl,
                                               const T tu,
                                               const T ul,
                                               const T uu,
                                               const T vl,
                                               const T vu,
                                               const Vec4 s,
                                               const Vec4 e,
                                               T* const SCCD_RESTRICT f) {
            if constexpr (is_vf) {
                sample_f_vf<T, Vec4>(tl, tu, ul, uu, vl, vu, s, e, f);
            } else {
                sample_f_ee<T, Vec4>(tl, tu, ul, uu, vl, vu, s, e, f);
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
            cond_mask &= ((fmin <= tol) && (fmax >= -tol)) ? MASK_FULL : 0;
            return cond_mask;
        }

        template <bool is_vf, typename T, typename Vec4>
        static inline __device__ void evaluate_cell_3d(const Domain<T>& cell,
                                                       const Vec4 sx,
                                                       const Vec4 sy,
                                                       const Vec4 sz,
                                                       const Vec4 ex,
                                                       const Vec4 ey,
                                                       const Vec4 ez,
                                                       const T tol,
                                                       const T* const SCCD_RESTRICT adaptive_tol,
                                                       int& contains_origin,
                                                       int& accept) {
            const T tl = cell.tlower, tu = cell.tupper;
            const T ul = cell.ulower, uu = cell.uupper;
            const T vl = cell.vlower, vu = cell.vupper;

            T fmin, fmax;
            T f[8];

            sample_f<is_vf, T, Vec4>(tl, tu, ul, uu, vl, vu, sx, ex, f);

            fminmax<T>(f, fmin, fmax);
            const T x_width = fmax - fmin;
            const uint8_t x_mask = cond_mask<T>(fmin, fmax, tol, adaptive_tol[0]);
            int co = (fmin <= tol) & (fmax >= -tol);

            sample_f<is_vf, T, Vec4>(tl, tu, ul, uu, vl, vu, sy, ey, f);

            fminmax<T>(f, fmin, fmax);
            const T y_width = fmax - fmin;
            const uint8_t y_mask = cond_mask<T>(fmin, fmax, tol, adaptive_tol[1]);
            co &= (fmin <= tol) & (fmax >= -tol);

            sample_f<is_vf, T, Vec4>(tl, tu, ul, uu, vl, vu, sz, ez, f);

            fminmax<T>(f, fmin, fmax);
            const T z_width = fmax - fmin;
            const uint8_t z_mask = cond_mask<T>(fmin, fmax, tol, adaptive_tol[2]);
            co &= (fmin <= tol) & (fmax >= -tol);

            const uint8_t and_mask = (x_mask & y_mask & z_mask);
            const T true_tol = device::max<T>(x_width, device::max<T>(y_width, z_width));

            const bool cond1 = and_mask & MASK_DOMAIN_SMALLER_THAN_TOL;
            const bool cond2 = and_mask & MASK_BOX_INSIDE_EPSILON_BOX;
            const bool cond3 = (tl > T(0)) && (true_tol < tol);
            const bool cond4 = and_mask & MASK_INTERVAL_TERMINAL;

            contains_origin = co;
            accept = co && (cond1 || cond2 || cond3 || cond4);
        }

        template <bool is_vf, typename T, typename Vec4>
        static inline __device__ void sample_cell_3d(const int ti,
                                                     const int ui,
                                                     const int vi,
                                                     const int nt,
                                                     const int nu,
                                                     const int nv,
                                                     const Domain<T>& parent,
                                                     const Vec4 sx,
                                                     const Vec4 sy,
                                                     const Vec4 sz,
                                                     const Vec4 ex,
                                                     const Vec4 ey,
                                                     const Vec4 ez,
                                                     const T tol,
                                                     const T* const SCCD_RESTRICT adaptive_tol,
                                                     int& contains_origin,
                                                     int& accept,
                                                     Domain<T>& cell) {
            const T t_h = (parent.tupper - parent.tlower) / nt;
            const T u_h = (parent.uupper - parent.ulower) / nu;
            const T v_h = (parent.vupper - parent.vlower) / nv;

            cell.tlower = parent.tlower + ti * t_h;
            cell.tupper = cell.tlower + t_h;
            cell.ulower = parent.ulower + ui * u_h;
            cell.uupper = cell.ulower + u_h;
            cell.vlower = parent.vlower + vi * v_h;
            cell.vupper = cell.vlower + v_h;

            evaluate_cell_3d<is_vf, T, Vec4>(cell, sx, sy, sz, ex, ey, ez, tol, adaptive_tol, contains_origin, accept);
        }

        template <typename T, typename Vec4, typename I>
        static inline __device__ void load_query_ee(const int qid,
                                                    const I* const SCCD_RESTRICT overlap0,
                                                    const I* const SCCD_RESTRICT overlap1,
                                                    T** const SCCD_RESTRICT sp,
                                                    T** const SCCD_RESTRICT ep,
                                                    const size_t element_stride,
                                                    I** const SCCD_RESTRICT elements,
                                                    Vec4& sx,
                                                    Vec4& sy,
                                                    Vec4& sz,
                                                    Vec4& ex,
                                                    Vec4& ey,
                                                    Vec4& ez) {
            const I ea = overlap0[qid];
            const I eb = overlap1[qid];

            const ptrdiff_t idxa0 = elements[0][ea * element_stride];
            const ptrdiff_t idxa1 = elements[1][ea * element_stride];
            const ptrdiff_t idxb0 = elements[0][eb * element_stride];
            const ptrdiff_t idxb1 = elements[1][eb * element_stride];

            sx.x = sp[0][idxa0];
            sx.y = sp[0][idxa1];
            sx.z = sp[0][idxb0];
            sx.w = sp[0][idxb1];
            sy.x = sp[1][idxa0];
            sy.y = sp[1][idxa1];
            sy.z = sp[1][idxb0];
            sy.w = sp[1][idxb1];
            sz.x = sp[2][idxa0];
            sz.y = sp[2][idxa1];
            sz.z = sp[2][idxb0];
            sz.w = sp[2][idxb1];

            ex.x = ep[0][idxa0];
            ex.y = ep[0][idxa1];
            ex.z = ep[0][idxb0];
            ex.w = ep[0][idxb1];
            ey.x = ep[1][idxa0];
            ey.y = ep[1][idxa1];
            ey.z = ep[1][idxb0];
            ey.w = ep[1][idxb1];
            ez.x = ep[2][idxa0];
            ez.y = ep[2][idxa1];
            ez.z = ep[2][idxb0];
            ez.w = ep[2][idxb1];
        }

        template <typename T, typename Vec4, typename I>
        static inline __device__ void load_query_vf(const int qid,
                                                    const I* const SCCD_RESTRICT voverlap,
                                                    const I* const SCCD_RESTRICT foverlap,
                                                    T** const SCCD_RESTRICT sp,
                                                    T** const SCCD_RESTRICT ep,
                                                    const size_t element_stride,
                                                    I** const SCCD_RESTRICT elements,
                                                    Vec4& sx,
                                                    Vec4& sy,
                                                    Vec4& sz,
                                                    Vec4& ex,
                                                    Vec4& ey,
                                                    Vec4& ez) {
            const I va = voverlap[qid];
            const I vb = foverlap[qid];

            const auto i0 = elements[0][vb * element_stride];
            const auto i1 = elements[1][vb * element_stride];
            const auto i2 = elements[2][vb * element_stride];

            // ---------------
            // Start
            // ---------------
            sx.x = sp[0][va];
            sx.y = sp[0][i0];
            sx.z = sp[0][i1];
            sx.w = sp[0][i2];

            sy.x = sp[1][va];
            sy.y = sp[1][i0];
            sy.z = sp[1][i1];
            sy.w = sp[1][i2];

            sz.x = sp[2][va];
            sz.y = sp[2][i0];
            sz.z = sp[2][i1];
            sz.w = sp[2][i2];

            // ---------------
            // End
            // ---------------
            ex.x = ep[0][va];
            ex.y = ep[0][i0];
            ex.z = ep[0][i1];
            ex.w = ep[0][i2];

            ey.x = ep[1][va];
            ey.y = ep[1][i0];
            ey.z = ep[1][i1];
            ey.w = ep[1][i2];

            ez.x = ep[2][va];
            ez.y = ep[2][i0];
            ez.z = ep[2][i1];
            ez.w = ep[2][i2];
        }

        template <bool is_vf, typename T, typename Vec4, typename I>
        static inline __device__ void load_query_and_tol(const int qid,
                                                         const I* const SCCD_RESTRICT overlap0,
                                                         const I* const SCCD_RESTRICT overlap1,
                                                         T** const SCCD_RESTRICT sp,
                                                         T** const SCCD_RESTRICT ep,
                                                         const size_t element_stride,
                                                         I** const SCCD_RESTRICT elements,
                                                         const T tol,
                                                         Vec4& sx,
                                                         Vec4& sy,
                                                         Vec4& sz,
                                                         Vec4& ex,
                                                         Vec4& ey,
                                                         Vec4& ez,
                                                         T* const SCCD_RESTRICT atol) {
            if constexpr (is_vf) {
                load_query_vf<T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, sx, sy, sz, ex, ey, ez);
                compute_face_vertex_tolerance<T, Vec4>(tol, sx, sy, sz, ex, ey, ez, &atol[0], &atol[1], &atol[2]);
            } else {
                load_query_ee<T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, sx, sy, sz, ex, ey, ez);
                compute_edge_edge_tolerance<T, Vec4>(tol, sx, sy, sz, ex, ey, ez, &atol[0], &atol[1], &atol[2]);
            }
        }

        // Per-slot validity tag stored in the qid array.  Three states:
        //   SCCD_QID_EMPTY   (-1): free, available to be claimed by a writer
        //   SCCD_QID_WRITING (-2): claimed by a writer, fields being filled
        //   qid >= 0             : committed, fields are safe to read
        static constexpr int SCCD_QID_EMPTY = -1;
        static constexpr int SCCD_QID_WRITING = -2;

        template <typename T>
        __global__ void init_narrow_phase_kernel(const size_t toi_n, const T max_toi, T* SCCD_RESTRICT toi) {
            const size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
            if (i < toi_n) toi[i] = max_toi;
        }

        __global__ void reset_batch_narrow_phase_kernel(int* SCCD_RESTRICT g_top, int* SCCD_RESTRICT g_request) {
            *g_top = 0;
            *g_request = 0;
        }

        static inline __device__ int reserve_slots(int* SCCD_RESTRICT counter, int k, int capacity) {
            int old = atomicAdd(counter, 0);
            while (true) {
                if (old + k > capacity) return -1;
                const int prev = atomicCAS(counter, old, old + k);
                if (prev == old) return old;
                old = prev;
            }
        }

        static inline __device__ int release_slot(int* SCCD_RESTRICT counter) {
            int old = atomicAdd(counter, 0);
            while (true) {
                if (old <= 0) return -1;
                const int prev = atomicCAS(counter, old, old - 1);
                if (prev == old) return old - 1;
                old = prev;
            }
        }

        template <bool Shared, typename T>
        static inline __device__ int try_pop(const Stack<T>& stk, Domain<T>& d, int& level, int& qid) {
            const int slot = release_slot(stk.top);
            if (slot < 0) return 0;

            int q;
            do {
                q = atomicAdd(&stk.qid[slot], 0);
            } while (q < 0);
            if (Shared) {
                __threadfence_block();
            } else {
                __threadfence();
            }

            d.tlower = stk.tlower[slot];
            d.tupper = stk.tupper[slot];
            d.ulower = stk.ulower[slot];
            d.uupper = stk.uupper[slot];
            d.vlower = stk.vlower[slot];
            d.vupper = stk.vupper[slot];
            level = stk.level[slot];
            qid = q;

            atomicExch(&stk.qid[slot], SCCD_QID_EMPTY);
            return 1;
        }

        template <int N>
        struct DfsSplit;

        template <>
        struct DfsSplit<64> {
            static constexpr int NT = 4;
            static constexpr int NU = 4;
            static constexpr int NV = 4;
        };

        template <>
        struct DfsSplit<128> {
            static constexpr int NT = 4;
            static constexpr int NU = 4;
            static constexpr int NV = 8;
        };

        template <>
        struct DfsSplit<256> {
            static constexpr int NT = 4;
            static constexpr int NU = 8;
            static constexpr int NV = 8;
        };

        template <int N>
        static inline __device__ int block_popc(const int pred, int* SCCD_RESTRICT warp_sums) {
            const int lane = threadIdx.x & 31;
            const int warp = threadIdx.x >> 5;
            int v = pred;
            for (int o = 16; o > 0; o >>= 1) {
                v += __shfl_xor_sync(0xffffffffu, v, o);
            }
            if (lane == 0) warp_sums[warp] = v;
            __syncthreads();
            if (warp == 0) {
                v = (lane < (N >> 5)) ? warp_sums[lane] : 0;
                for (int o = 16; o > 0; o >>= 1) {
                    v += __shfl_xor_sync(0xffffffffu, v, o);
                }
                if (lane == 0) warp_sums[0] = v;
            }
            __syncthreads();
            return warp_sums[0];
        }

        template <typename T>
        static inline __device__ void bisect_longest_axis(const Domain<T>& in,
                                                          const T* const SCCD_RESTRICT atol,
                                                          Domain<T>& left,
                                                          Domain<T>& right) {
            left = in;
            right = in;

            const T dt = (in.tupper - in.tlower);  // / atol[0];
            const T du = (in.uupper - in.ulower);  // / atol[1];
            const T dv = (in.vupper - in.vlower);  // / atol[2];

            if (dt >= du && dt >= dv) {
                const T m = (in.tlower + in.tupper) * T(0.5);
                left.tupper = m;
                right.tlower = m;
            } else if (du >= dv) {
                const T m = (in.ulower + in.uupper) * T(0.5);
                left.uupper = m;
                right.ulower = m;
            } else {
                const T m = (in.vlower + in.vupper) * T(0.5);
                left.vupper = m;
                right.vlower = m;
            }
        }

        template <typename T>
        static inline __device__ int widest_dimension(const Domain<T>& in) {
            const T dt = in.tupper - in.tlower;
            const T du = in.uupper - in.ulower;
            const T dv = in.vupper - in.vlower;
            if (du > dt && du >= dv) return 1;
            if (dv > dt && dv > du) return 2;
            return 0;
        }

        template <bool is_vf, typename T, typename Vec4>
        static inline __device__ void adaptive_axis_terms_component(const int split_dim,
                                                                    const Vec4 s,
                                                                    const Vec4 e,
                                                                    const T mid_t,
                                                                    const T mid_u,
                                                                    const T mid_v,
                                                                    T& F_base,
                                                                    T& J_axis) {
            const T base_t = split_dim == 0 ? T(0) : mid_t;
            const T base_u = split_dim == 1 ? T(0) : mid_u;
            const T base_v = split_dim == 2 ? T(0) : mid_v;

            if constexpr (is_vf) {
                const T base_omt = T(1) - base_t;
                const T base_o = T(1) - base_u - base_v;
                const T vertex = base_omt * s.x + base_t * e.x;
                const T face = base_omt * (base_o * s.y + base_u * s.z + base_v * s.w) +
                               base_t * (base_o * e.y + base_u * e.z + base_v * e.w);
                F_base = vertex - face;

                if (split_dim == 0) {
                    const T o = T(1) - mid_u - mid_v;
                    J_axis = (e.x - s.x) - (o * (e.y - s.y) + mid_u * (e.z - s.z) + mid_v * (e.w - s.w));
                } else if (split_dim == 1) {
                    J_axis = -((T(1) - mid_t) * (s.z - s.y) + mid_t * (e.z - e.y));
                } else {
                    J_axis = -((T(1) - mid_t) * (s.w - s.y) + mid_t * (e.w - e.y));
                }
            } else {
                const T ea0 = (e.x - s.x) * base_t + s.x;
                const T ea1 = (e.y - s.y) * base_t + s.y;
                const T eb0 = (e.z - s.z) * base_t + s.z;
                const T eb1 = (e.w - s.w) * base_t + s.w;
                F_base = ((ea1 - ea0) * base_u + ea0) - ((eb1 - eb0) * base_v + eb0);

                if (split_dim == 0) {
                    J_axis = (T(1) - mid_u) * (e.x - s.x) + mid_u * (e.y - s.y) - (T(1) - mid_v) * (e.z - s.z) -
                             mid_v * (e.w - s.w);
                } else if (split_dim == 1) {
                    J_axis = (T(1) - mid_t) * (s.y - s.x) + mid_t * (e.y - e.x);
                } else {
                    J_axis = -((T(1) - mid_t) * (s.w - s.z) + mid_t * (e.w - e.z));
                }
            }
        }

        template <bool is_vf, typename T, typename Vec4>
        static inline __device__ void adaptive_split_longest_axis(const Domain<T>& in,
                                                                  const Vec4 sx,
                                                                  const Vec4 sy,
                                                                  const Vec4 sz,
                                                                  const Vec4 ex,
                                                                  const Vec4 ey,
                                                                  const Vec4 ez,
                                                                  Domain<T>& left,
                                                                  Domain<T>& right) {
            left = in;
            right = in;

            const int split_dim = widest_dimension<T>(in);
            const T lo = split_dim == 0 ? in.tlower : (split_dim == 1 ? in.ulower : in.vlower);
            const T hi = split_dim == 0 ? in.tupper : (split_dim == 1 ? in.uupper : in.vupper);
            const T h = (hi - lo) * T(0.5);
            const T radius = h * T(0.6);
            const T x0 = lo + h;
            const T mid_t = (in.tlower + in.tupper) * T(0.5);
            const T mid_u = (in.ulower + in.uupper) * T(0.5);
            const T mid_v = (in.vlower + in.vupper) * T(0.5);

            T Fx;
            T Fy;
            T Fz;
            T Jx;
            T Jy;
            T Jz;
            adaptive_axis_terms_component<is_vf, T, Vec4>(split_dim, sx, ex, mid_t, mid_u, mid_v, Fx, Jx);
            adaptive_axis_terms_component<is_vf, T, Vec4>(split_dim, sy, ey, mid_t, mid_u, mid_v, Fy, Jy);
            adaptive_axis_terms_component<is_vf, T, Vec4>(split_dim, sz, ez, mid_t, mid_u, mid_v, Fz, Jz);

            const T H_axis = Jx * Jx + Jy * Jy + Jz * Jz;
            const T eps = sizeof(T) == sizeof(float) ? T(1.1920928955078125e-7) : T(2.2204460492503131e-16);
            const T step_scale = H_axis > eps ? T(1) / H_axis : T(0.00001);
            const T g = (Fx + x0 * Jx) * Jx + (Fy + x0 * Jy) * Jy + (Fz + x0 * Jz) * Jz;
            const T xmin = device::max<T>(lo, x0 - radius);
            const T xmax = device::min<T>(hi, x0 + radius);
            const T m = device::min<T>(xmax, device::max<T>(xmin, x0 - g * step_scale));

            if (split_dim == 0) {
                left.tupper = m;
                right.tlower = m;
            } else if (split_dim == 1) {
                left.uupper = m;
                right.ulower = m;
            } else {
                left.vupper = m;
                right.vlower = m;
            }
        }

        template <bool is_vf, int N, typename T, typename I>
        static __device__ void narrow_phase_dfs_zero_stride_body(const I* const SCCD_RESTRICT overlap0,
                                                                 const I* const SCCD_RESTRICT overlap1,
                                                                 T** const SCCD_RESTRICT sp,
                                                                 T** const SCCD_RESTRICT ep,
                                                                 const size_t element_stride,
                                                                 I** const SCCD_RESTRICT elements,
                                                                 const T tol,
                                                                 const int max_depth,
                                                                 T* SCCD_RESTRICT toi,
                                                                 Stack<T> g_stack,
                                                                 int qid_in,
                                                                 Domain<T> cur_in,
                                                                 int level_in,
                                                                 int active_in) {
            static_assert(N == 64 || N == 128 || N == 256, "SCCD_NP_THREADS_PER_BLOCK must be one of 64/128/256");
            using Vec4 = typename device::Vec4Type<T>::type;
            constexpr int S_CAP = SCCD_NP_SHARED_STACK_CAP;
            __shared__ T s_tlower[S_CAP];
            __shared__ T s_tupper[S_CAP];
            __shared__ T s_ulower[S_CAP];
            __shared__ T s_uupper[S_CAP];
            __shared__ T s_vlower[S_CAP];
            __shared__ T s_vupper[S_CAP];
            __shared__ int s_level[S_CAP];
            __shared__ int s_qid[S_CAP];
            __shared__ int s_top;
            __shared__ T s_toi;
            __shared__ int warp_sums[N >> 5];

            const int tid = threadIdx.x;

            if (tid == 0) {
                s_top = 0;
                s_toi = toi[0];
            }
            __syncthreads();

            Stack<T> s_stack = {s_tlower,
                                s_tupper,
                                s_ulower,
                                s_uupper,
                                s_vlower,
                                s_vupper,
                                s_level,
                                s_qid,
                                &s_top,
                                /*request=*/(int*)nullptr,
                                S_CAP};

            int qid = qid_in;
            Domain<T> cur = cur_in;
            int level = level_in;
            int active = active_in;
            Vec4 sx, sy, sz, ex, ey, ez;
            T atol[3] = {T(0), T(0), T(0)};

            if (active) {
                load_query_and_tol<is_vf, T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol);
            }

            while (true) {
                if (tid == 0) {
                    const T g = device::atomic_min(&toi[0], s_toi);
                    if (g < s_toi) s_toi = g;
                }
                __syncthreads();

                if (active && cur.tlower >= s_toi) active = 0;

                if (active && level >= max_depth) {
                    device::atomic_min(&s_toi, cur.tlower);
                    active = 0;
                }

                Domain<T> push_box;
                int push_level = 0;
                int will_push = 0;

                if (active) {
                    Domain<T> left, right;
                    cur.tupper = device::min<T>(cur.tupper, s_toi);

                    if (SCCD_CUDA_ADAPTIVE_SPLIT) {
                        adaptive_split_longest_axis<is_vf, T, Vec4>(cur, sx, sy, sz, ex, ey, ez, left, right);
                    } else {
                        bisect_longest_axis<T>(cur, atol, left, right);
                    }
                    int cl = 0, cr = 0, al = 0, ar = 0;
                    evaluate_cell_3d<is_vf, T, Vec4>(left, sx, sy, sz, ex, ey, ez, tol, atol, cl, al);
                    evaluate_cell_3d<is_vf, T, Vec4>(right, sx, sy, sz, ex, ey, ez, tol, atol, cr, ar);

                    if (!is_domain_valid<is_vf>(left, s_toi, atol)) {
                        cl = 0;
                        al = 0;
                    }
                    if (!is_domain_valid<is_vf>(right, s_toi, atol)) {
                        cr = 0;
                        ar = 0;
                    }

                    if (al) {
                        device::atomic_min(&s_toi, left.tlower);
                        cl = 0;
                    }
                    if (ar) {
                        device::atomic_min(&s_toi, right.tlower);
                        cr = 0;
                    }

                    if (cl && cr) {
                        const bool left_first = left.tlower <= right.tlower;
                        cur = left_first ? left : right;
                        push_box = left_first ? right : left;
                        push_level = level + 1;
                        will_push = 1;
                        level += 1;
                    } else if (cl) {
                        cur = left;
                        level += 1;
                    } else if (cr) {
                        cur = right;
                        level += 1;
                    } else {
                        active = 0;
                    }
                }

                if (will_push) {
                    const int slot = reserve_slots(&s_top, 1, S_CAP);
                    if (slot >= 0 && slot < S_CAP) {
                        s_stack.tlower[slot] = push_box.tlower;
                        s_stack.tupper[slot] = push_box.tupper;
                        s_stack.ulower[slot] = push_box.ulower;
                        s_stack.uupper[slot] = push_box.uupper;
                        s_stack.vlower[slot] = push_box.vlower;
                        s_stack.vupper[slot] = push_box.vupper;
                        s_stack.level[slot] = push_level;
                        s_stack.qid[slot] = qid;
                    } else {
                        const int g_slot = reserve_slots(g_stack.top, 1, g_stack.capacity);
                        if (g_slot < 0) {
                            atomicAdd(g_stack.request, 1);
                        } else {
                            while (atomicCAS(&g_stack.qid[g_slot], SCCD_QID_EMPTY, SCCD_QID_WRITING) !=
                                   SCCD_QID_EMPTY) {
                                // busy-wait
                            }
                            g_stack.tlower[g_slot] = push_box.tlower;
                            g_stack.tupper[g_slot] = push_box.tupper;
                            g_stack.ulower[g_slot] = push_box.ulower;
                            g_stack.uupper[g_slot] = push_box.uupper;
                            g_stack.vlower[g_slot] = push_box.vlower;
                            g_stack.vupper[g_slot] = push_box.vupper;
                            g_stack.level[g_slot] = push_level;
                            __threadfence();
                            atomicExch(&g_stack.qid[g_slot], qid);
                        }
                    }
                }

                __syncthreads();

                if (!active) {
                    int t = atomicAdd(&s_top, 0);
                    while (t > 0) {
                        const int prev = atomicCAS(&s_top, t, t - 1);
                        if (prev == t) {
                            const int slot = t - 1;
                            if (slot >= 0 && slot < S_CAP) {
                                const int new_qid = s_stack.qid[slot];
                                cur.tlower = s_stack.tlower[slot];
                                cur.tupper = s_stack.tupper[slot];
                                cur.ulower = s_stack.ulower[slot];
                                cur.uupper = s_stack.uupper[slot];
                                cur.vlower = s_stack.vlower[slot];
                                cur.vupper = s_stack.vupper[slot];
                                level = s_stack.level[slot];
                                if (new_qid != qid) {
                                    qid = new_qid;
                                    load_query_and_tol<is_vf, T, Vec4, I>(qid,
                                                                          overlap0,
                                                                          overlap1,
                                                                          sp,
                                                                          ep,
                                                                          element_stride,
                                                                          elements,
                                                                          tol,
                                                                          sx,
                                                                          sy,
                                                                          sz,
                                                                          ex,
                                                                          ey,
                                                                          ez,
                                                                          atol);
                                }
                                active = 1;
                            }
                            break;
                        }
                        t = prev;
                    }
                }

                __syncthreads();

                const int n_active = block_popc<N>(active, warp_sums);
                if (n_active == 0) {
                    const int s_now = atomicAdd(&s_top, 0);
                    if (s_now <= 0) break;
                }
            }

            if (tid == 0) device::atomic_min(&toi[0], s_toi);
        }

        template <bool is_vf, int N, typename T, typename I>
        __global__ void narrow_phase_dfs_zero_stride_kernel(const I* const SCCD_RESTRICT overlap0,
                                                            const I* const SCCD_RESTRICT overlap1,
                                                            T** const SCCD_RESTRICT sp,
                                                            T** const SCCD_RESTRICT ep,
                                                            const size_t element_stride,
                                                            I** const SCCD_RESTRICT elements,
                                                            const T tol,
                                                            const int max_depth,
                                                            T* SCCD_RESTRICT toi,
                                                            Stack<T> g_stack,
                                                            const int seed_begin,
                                                            const int seed_end) {
            using Vec4 = typename device::Vec4Type<T>::type;

            const int tid = threadIdx.x;
            const int my_seed = seed_begin + (int)blockIdx.x * N + tid;
            const bool has_seed = my_seed < seed_end;

            int qid = -1;
            Domain<T> cur = {T(0), T(0), T(0), T(0), T(0), T(0)};
            int level = 0;
            int active = 0;

            if (has_seed) {
                qid = my_seed;
                Vec4 sx, sy, sz, ex, ey, ez;
                T atol[3];
                load_query_and_tol<is_vf, T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol);

                Domain<T> root = {T(0), T(1), T(0), T(1), T(0), T(1)};
                int contains = 0;
                int accept = 0;
                evaluate_cell_3d<is_vf, T, Vec4>(root, sx, sy, sz, ex, ey, ez, tol, atol, contains, accept);

                if (contains && is_domain_valid<is_vf>(root, toi[0], atol)) {
                    cur = root;
                    level = 0;
                    active = 1;
                }
            }

            narrow_phase_dfs_zero_stride_body<is_vf, N, T, I>(overlap0,
                                                              overlap1,
                                                              sp,
                                                              ep,
                                                              element_stride,
                                                              elements,
                                                              tol,
                                                              max_depth,
                                                              toi,
                                                              g_stack,
                                                              qid,
                                                              cur,
                                                              level,
                                                              active);
        }

        template <bool is_vf, int N, typename T, typename I>
        __global__ void narrow_phase_dfs_zero_stride_from_stack_kernel(const I* const SCCD_RESTRICT overlap0,
                                                                       const I* const SCCD_RESTRICT overlap1,
                                                                       T** const SCCD_RESTRICT sp,
                                                                       T** const SCCD_RESTRICT ep,
                                                                       const size_t element_stride,
                                                                       I** const SCCD_RESTRICT elements,
                                                                       const T tol,
                                                                       const int max_depth,
                                                                       T* SCCD_RESTRICT toi,
                                                                       Stack<T> g_stack) {
            int qid = -1;
            Domain<T> cur = {T(0), T(0), T(0), T(0), T(0), T(0)};
            int level = 0;
            int active = 0;

            if (try_pop<false, T>(g_stack, cur, level, qid)) {
                active = 1;
            }

            narrow_phase_dfs_zero_stride_body<is_vf, N, T, I>(overlap0,
                                                              overlap1,
                                                              sp,
                                                              ep,
                                                              element_stride,
                                                              elements,
                                                              tol,
                                                              max_depth,
                                                              toi,
                                                              g_stack,
                                                              qid,
                                                              cur,
                                                              level,
                                                              active);
        }

        template <bool is_vf, int N, typename T, typename I>
        static __device__ void narrow_phase_dfs_body(const I* const SCCD_RESTRICT overlap0,
                                                     const I* const SCCD_RESTRICT overlap1,
                                                     T** const SCCD_RESTRICT sp,
                                                     T** const SCCD_RESTRICT ep,
                                                     const size_t element_stride,
                                                     I** const SCCD_RESTRICT elements,
                                                     const T tol,
                                                     const int max_depth,
                                                     T* SCCD_RESTRICT toi,
                                                     const int toi_stride,
                                                     Stack<T> g_stack,
                                                     const T alpha,
                                                     const int qid,
                                                     Domain<T> sampling_root,
                                                     const int initial_level,
                                                     const bool do_hard_defer) {
            static_assert(N == 64 || N == 128 || N == 256, "SCCD_NP_THREADS_PER_BLOCK must be one of 64/128/256");
            using Vec4 = typename device::Vec4Type<T>::type;
            constexpr int NT = DfsSplit<N>::NT;
            constexpr int NU = DfsSplit<N>::NU;
            constexpr int NV = DfsSplit<N>::NV;
            constexpr int S_CAP = SCCD_NP_SHARED_STACK_CAP;
            const int tid = threadIdx.x;

            __shared__ T s_tlower[S_CAP];
            __shared__ T s_tupper[S_CAP];
            __shared__ T s_ulower[S_CAP];
            __shared__ T s_uupper[S_CAP];
            __shared__ T s_vlower[S_CAP];
            __shared__ T s_vupper[S_CAP];
            __shared__ int s_level[S_CAP];
            __shared__ int s_qid[S_CAP];
            __shared__ int s_top;
            __shared__ T s_toi;
            __shared__ int s_hard;
            __shared__ int s_defer_base;
            __shared__ int s_defer_cursor;
            __shared__ int warp_sums[N >> 5];

            const int toi_idx = qid * toi_stride;

            for (int i = tid; i < S_CAP; i += N) s_qid[i] = SCCD_QID_EMPTY;
            if (tid == 0) {
                s_top = 0;
                s_toi = toi[toi_idx];
                s_hard = 0;
                s_defer_base = -1;
                s_defer_cursor = 0;
            }
            __syncthreads();

            Stack<T> s_stack = {s_tlower,
                                s_tupper,
                                s_ulower,
                                s_uupper,
                                s_vlower,
                                s_vupper,
                                s_level,
                                s_qid,
                                &s_top,
                                /*request=*/(int*)nullptr,
                                S_CAP};

            T atol[3];
            Vec4 sx, sy, sz, ex, ey, ez;

            load_query_and_tol<is_vf, T, Vec4, I>(
                qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol);

            const int ti = tid / (NU * NV);
            const int rem = tid % (NU * NV);
            const int ui = rem / NV;
            const int vi = rem % NV;

            Domain<T> cur;
            int contains = 0;
            int accept = 0;
            sample_cell_3d<is_vf, T, Vec4>(
                ti, ui, vi, NT, NU, NV, sampling_root, sx, sy, sz, ex, ey, ez, tol, atol, contains, accept, cur);

            if (accept && is_domain_valid<is_vf>(cur, s_toi, atol)) {
                device::atomic_min(&s_toi, cur.tlower);
            }
            __syncthreads();

            const int active_seed = contains && !accept && (cur.tlower < s_toi);

            const int co_count = block_popc<N>(active_seed, warp_sums);
            if (!co_count) {
                if (tid == 0) device::atomic_min(&toi[toi_idx], s_toi);
                return;
            }

            if (do_hard_defer && tid == 0) {
                if ((T)co_count > alpha * (T)N) {
                    s_hard = 1;
                    const int base = reserve_slots(g_stack.top, co_count, g_stack.capacity);
                    if (base < 0) {
                        // Deficit will be allocated on the host's retry pass.
                        atomicAdd(g_stack.request, co_count);
                    }
                    s_defer_base = base;
                    s_defer_cursor = 0;
                }
            }
            __syncthreads();

            if (s_hard) {
                if (active_seed && s_defer_base >= 0) {
                    const int rank = atomicAdd(&s_defer_cursor, 1);
                    const int slot = s_defer_base + rank;
                    if (slot >= 0 && slot < g_stack.capacity) {
                        g_stack.tlower[slot] = cur.tlower;
                        g_stack.tupper[slot] = sccd::min<T>(cur.tupper, s_toi);
                        g_stack.ulower[slot] = cur.ulower;
                        g_stack.uupper[slot] = cur.uupper;
                        g_stack.vlower[slot] = cur.vlower;
                        g_stack.vupper[slot] = cur.vupper;
                        g_stack.level[slot] = initial_level;
                        __threadfence();
                        atomicExch(&g_stack.qid[slot], qid);
                    }
                }
                if (tid == 0) device::atomic_min(&toi[toi_idx], s_toi);
                return;
            }

            int active = active_seed;
            int level = initial_level + 1;
            while (true) {
                if (toi_stride == 0) {
                    if (tid == 0) {
                        const T g = device::atomic_min(&toi[toi_idx], s_toi);
                        if (g < s_toi) s_toi = g;
                    }
                    __syncthreads();
                }

                if (active && cur.tlower >= s_toi) active = 0;

                if (active && level >= max_depth) {
                    if (is_domain_valid<is_vf>(cur, s_toi, atol)) {
                        device::atomic_min(&s_toi, cur.tlower);
                    }
                    active = 0;
                }

                Domain<T> push_box;
                int push_level = 0;
                int will_push = 0;

                if (active) {
                    Domain<T> left, right;

                    if (SCCD_CUDA_ADAPTIVE_SPLIT) {
                        adaptive_split_longest_axis<is_vf, T, Vec4>(cur, sx, sy, sz, ex, ey, ez, left, right);
                    } else {
                        bisect_longest_axis<T>(cur, atol, left, right);
                    }

                    int cl = 0, cr = 0, al = 0, ar = 0;
                    evaluate_cell_3d<is_vf, T, Vec4>(left, sx, sy, sz, ex, ey, ez, tol, atol, cl, al);
                    evaluate_cell_3d<is_vf, T, Vec4>(right, sx, sy, sz, ex, ey, ez, tol, atol, cr, ar);

                    if (!is_domain_valid<is_vf>(left, s_toi, atol)) {
                        cl = 0;
                        al = 0;
                    }
                    if (!is_domain_valid<is_vf>(right, s_toi, atol)) {
                        cr = 0;
                        ar = 0;
                    }

                    if (al) {
                        device::atomic_min(&s_toi, left.tlower);
                        cl = 0;
                    }
                    if (ar) {
                        device::atomic_min(&s_toi, right.tlower);
                        cr = 0;
                    }

                    if (cl && cr) {
                        const bool left_first = left.tlower <= right.tlower;
                        cur = left_first ? left : right;
                        push_box = left_first ? right : left;
                        push_level = level + 1;
                        will_push = 1;
                        level += 1;
                    } else if (cl) {
                        cur = left;
                        level += 1;
                    } else if (cr) {
                        cur = right;
                        level += 1;
                    } else {
                        active = 0;
                    }
                }

                if (will_push) {
                    const int slot = reserve_slots(&s_top, 1, S_CAP);
                    if (slot >= 0 && slot < S_CAP) {
                        s_stack.tlower[slot] = push_box.tlower;
                        s_stack.tupper[slot] = push_box.tupper;
                        s_stack.ulower[slot] = push_box.ulower;
                        s_stack.uupper[slot] = push_box.uupper;
                        s_stack.vlower[slot] = push_box.vlower;
                        s_stack.vupper[slot] = push_box.vupper;
                        s_stack.level[slot] = push_level;
                        s_stack.qid[slot] = qid;
                    } else {
                        const int g_slot = reserve_slots(g_stack.top, 1, g_stack.capacity);
                        if (g_slot < 0) {
                            atomicAdd(g_stack.request, 1);
                        } else {
                            while (atomicCAS(&g_stack.qid[g_slot], SCCD_QID_EMPTY, SCCD_QID_WRITING) !=
                                   SCCD_QID_EMPTY) {
                                // busy-wait
                            }
                            g_stack.tlower[g_slot] = push_box.tlower;
                            g_stack.tupper[g_slot] = push_box.tupper;
                            g_stack.ulower[g_slot] = push_box.ulower;
                            g_stack.uupper[g_slot] = push_box.uupper;
                            g_stack.vlower[g_slot] = push_box.vlower;
                            g_stack.vupper[g_slot] = push_box.vupper;
                            g_stack.level[g_slot] = push_level;
                            __threadfence();
                            atomicExch(&g_stack.qid[g_slot], qid);
                        }
                    }
                }

                __syncthreads();

                if (!active) {
                    int t = atomicAdd(&s_top, 0);
                    while (t > 0) {
                        const int prev = atomicCAS(&s_top, t, t - 1);
                        if (prev == t) {
                            const int slot = t - 1;
                            if (slot >= 0 && slot < S_CAP) {
                                cur.tlower = s_stack.tlower[slot];
                                cur.tupper = s_stack.tupper[slot];
                                cur.ulower = s_stack.ulower[slot];
                                cur.uupper = s_stack.uupper[slot];
                                cur.vlower = s_stack.vlower[slot];
                                cur.vupper = s_stack.vupper[slot];
                                level = s_stack.level[slot];
                                active = 1;
                            }
                            break;
                        }
                        t = prev;
                    }
                }

                __syncthreads();

                const int n_active = block_popc<N>(active, warp_sums);
                if (n_active == 0 && s_top == 0) break;
            }

            if (tid == 0) device::atomic_min(&toi[toi_idx], s_toi);
        }

        template <bool is_vf, int N, typename T, typename I>
        __global__ void narrow_phase_dfs_kernel(const I* const SCCD_RESTRICT overlap0,
                                                const I* const SCCD_RESTRICT overlap1,
                                                T** const SCCD_RESTRICT sp,
                                                T** const SCCD_RESTRICT ep,
                                                const size_t element_stride,
                                                I** const SCCD_RESTRICT elements,
                                                const T tol,
                                                const int max_depth,
                                                T* SCCD_RESTRICT toi,
                                                const int toi_stride,
                                                Stack<T> g_stack,
                                                const T alpha,
                                                const int seed_begin,
                                                const int seed_end) {
            const int qid = seed_begin + (int)blockIdx.x;
            if (qid >= seed_end) return;

            Domain<T> root = {T(0), T(1), T(0), T(1), T(0), T(1)};
            narrow_phase_dfs_body<is_vf, N, T, I>(overlap0,
                                                  overlap1,
                                                  sp,
                                                  ep,
                                                  element_stride,
                                                  elements,
                                                  tol,
                                                  max_depth,
                                                  toi,
                                                  toi_stride,
                                                  g_stack,
                                                  alpha,
                                                  qid,
                                                  root,
                                                  /*initial_level=*/0,
                                                  /*do_hard_defer=*/true);
        }

        template <bool is_vf, int N, typename T, typename I>
        __global__ void narrow_phase_dfs_from_stack_kernel(const I* const SCCD_RESTRICT overlap0,
                                                           const I* const SCCD_RESTRICT overlap1,
                                                           T** const SCCD_RESTRICT sp,
                                                           T** const SCCD_RESTRICT ep,
                                                           const size_t element_stride,
                                                           I** const SCCD_RESTRICT elements,
                                                           const T tol,
                                                           const int max_depth,
                                                           T* SCCD_RESTRICT toi,
                                                           const int toi_stride,
                                                           Stack<T> g_stack) {
            __shared__ int b_qid;
            __shared__ int b_level;
            __shared__ int b_have_work;
            __shared__ Domain<T> b_cur;

            if (threadIdx.x == 0) {
                Domain<T> popped_cur;
                int popped_level = 0;
                int popped_qid = -1;
                if (try_pop<false, T>(g_stack, popped_cur, popped_level, popped_qid)) {
                    b_qid = popped_qid;
                    b_level = popped_level;
                    b_cur = popped_cur;
                    b_have_work = 1;
                } else {
                    b_have_work = 0;
                }
            }
            __syncthreads();

            if (!b_have_work) return;

            narrow_phase_dfs_body<is_vf, N, T, I>(overlap0,
                                                  overlap1,
                                                  sp,
                                                  ep,
                                                  element_stride,
                                                  elements,
                                                  tol,
                                                  max_depth,
                                                  toi,
                                                  toi_stride,
                                                  g_stack,
                                                  T(0),
                                                  b_qid,
                                                  b_cur,
                                                  b_level,
                                                  /*do_hard_defer=*/false);
        }

        template <bool is_vf, typename T, typename I>
        int narrow_phase_generic(const size_t noverlaps,
                                 const I* const SCCD_RESTRICT overlap0,
                                 const I* const SCCD_RESTRICT overlap1,
                                 // Geometric data
                                 T** const SCCD_RESTRICT v0,
                                 T** const SCCD_RESTRICT v1,
                                 const size_t element_stride,
                                 I** const SCCD_RESTRICT elements,
                                 // Output
                                 const T max_toi,
                                 T* const SCCD_RESTRICT d_toi,
                                 const int max_depth,
                                 const T tol,
                                 const int toi_stride) {
            SCCD_CUDA_LAST_ERROR();

            if (noverlaps == 0) return max_toi;
            assert(d_toi != nullptr);
            assert(toi_stride == 0 || toi_stride == 1);

            // toi length: 1 when stride==0 (all candidates share toi[0]),
            //             noverlaps when stride==1 (one toi per candidate).
            const size_t toi_n = (toi_stride == 0) ? 1 : noverlaps;

            T SCCD_NP_ALPHA = T(0.5);
            {
                double alpha = (double)SCCD_NP_ALPHA;
                SCCD_READ_ENV(alpha, atof);
                SCCD_NP_ALPHA = (T)alpha;
            }

            // ----------------------------------------------------------------
            // Auto-sized hyperparameters.
            //
            //   SCCD_BLOCKS_PER_SM    -> from CUDA occupancy API
            //   SCCD_GSTACK_CAP_INIT  -> initial global-stack capacity
            //                            (default 0; grown on demand from
            //                            kernel-reported deficit)
            //   SCCD_GSTACK_CAP_MAX   -> soft cap on a single grow step
            //   SCCD_BATCH_SIZE       -> candidates per outer iteration
            //                            (default: noverlaps)
            //
            // Any of these can be overridden via the matching env var.
            // ----------------------------------------------------------------
            int dev = 0;
            cudaGetDevice(&dev);
            cudaDeviceProp prop{};
            cudaGetDeviceProperties(&prop, dev);

            constexpr int N = SCCD_NP_THREADS_PER_BLOCK;

            int SCCD_BLOCKS_PER_SM = 4;
            {
                int occ = 0;
                if (cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                        &occ, (const void*)narrow_phase_dfs_zero_stride_from_stack_kernel<is_vf, N, T, I>, N, 0) ==
                        cudaSuccess &&
                    occ > 0) {
                    SCCD_BLOCKS_PER_SM = occ;
                }
            }

            SCCD_READ_ENV(SCCD_BLOCKS_PER_SM, atoi);
            // printf("SCCD_BLOCKS_PER_SM: %d, SCCD_NP_THREADS_PER_BLOCK: %d\n",
            //        SCCD_BLOCKS_PER_SM,
            //        SCCD_NP_THREADS_PER_BLOCK);

            int base_grid_blocks = prop.multiProcessorCount * SCCD_BLOCKS_PER_SM;
            if (base_grid_blocks <= 0) base_grid_blocks = 1;

            // Batch size: candidates processed per outer iteration.
            int SCCD_BATCH_SIZE = 0;
            SCCD_READ_ENV(SCCD_BATCH_SIZE, atoi);
            const size_t batch_size = (SCCD_BATCH_SIZE > 0) ? (size_t)SCCD_BATCH_SIZE : noverlaps;

            int gstack_cap = 0;

            int SCCD_GSTACK_CAP_MAX = INT_MAX;
            SCCD_READ_ENV(SCCD_GSTACK_CAP_MAX, atoi);

            T* g_tlower = nullptr;
            T* g_tupper = nullptr;
            T* g_ulower = nullptr;
            T* g_uupper = nullptr;
            T* g_vlower = nullptr;
            T* g_vupper = nullptr;
            int* g_level = nullptr;
            int* g_qid = nullptr;
            int* g_top = nullptr;
            int* g_request = nullptr;
            cudaMalloc(&g_top, sizeof(int));
            cudaMalloc(&g_request, sizeof(int));

            auto grow_stack = [&](int new_cap) {
                if (new_cap <= gstack_cap) return;
                cudaFree(g_tlower);
                cudaFree(g_tupper);
                cudaFree(g_ulower);
                cudaFree(g_uupper);
                cudaFree(g_vlower);
                cudaFree(g_vupper);
                cudaFree(g_level);
                cudaFree(g_qid);
                g_tlower = nullptr;
                g_tupper = nullptr;
                g_ulower = nullptr;
                g_uupper = nullptr;
                g_vlower = nullptr;
                g_vupper = nullptr;
                g_level = nullptr;
                g_qid = nullptr;
                cudaMalloc(&g_tlower, new_cap * sizeof(T));
                cudaMalloc(&g_tupper, new_cap * sizeof(T));
                cudaMalloc(&g_ulower, new_cap * sizeof(T));
                cudaMalloc(&g_uupper, new_cap * sizeof(T));
                cudaMalloc(&g_vlower, new_cap * sizeof(T));
                cudaMalloc(&g_vupper, new_cap * sizeof(T));
                cudaMalloc(&g_level, new_cap * sizeof(int));
                cudaMalloc(&g_qid, new_cap * sizeof(int));
                // SCCD_QID_EMPTY == -1, so 0xFF byte pattern initialises
                // every int to -1 without launching the init kernel.
                cudaMemsetAsync(g_qid, 0xFF, new_cap * sizeof(int));
                gstack_cap = new_cap;
            };

            if (gstack_cap > 0) grow_stack(gstack_cap);

            // One-time toi init.  Stack arrays are sized on demand below.
            {
                const int block = 256;
                const size_t grid_sz = (toi_n + block - 1) / block;
                const int grid = (grid_sz > (size_t)INT_MAX) ? INT_MAX : (int)grid_sz;
                init_narrow_phase_kernel<T><<<grid, block>>>(toi_n, max_toi, d_toi);
                SCCD_CUDA_LAST_ERROR();
            }

            dim3 block_pass1(SCCD_NP_THREADS_PER_BLOCK, 1, 1);

            // Batch loop.  Each batch runs Pass 1 (seed-driven) then a
            // Pass 2 drain loop on whatever spilled to g_stack.  If any
            // push overflowed (g_request > 0) the entire batch is
            // retried after growing the stack -- TOIs are preserved
            // across retries, so subsequent attempts prune more
            // aggressively.
            for (size_t begin = 0; begin < noverlaps; begin += batch_size) {
                const size_t end = (begin + batch_size < noverlaps) ? (begin + batch_size) : noverlaps;
                const int this_batch = (int)(end - begin);

                while (true) {
                    reset_batch_narrow_phase_kernel<<<1, 1>>>(g_top, g_request);
                    SCCD_CUDA_LAST_ERROR();

                    Stack<T> g_stack = {g_tlower,
                                        g_tupper,
                                        g_ulower,
                                        g_uupper,
                                        g_vlower,
                                        g_vupper,
                                        g_level,
                                        g_qid,
                                        g_top,
                                        g_request,
                                        gstack_cap};

                    // Pass 1: seed-driven.
                    if (toi_stride == 0) {
                        const int grid_blocks_zs = (this_batch + N - 1) / N;
                        dim3 grid_pass1_zs(grid_blocks_zs, 1, 1);
                        narrow_phase_dfs_zero_stride_kernel<is_vf, N, T, I>
                            <<<grid_pass1_zs, block_pass1>>>(overlap0,
                                                             overlap1,
                                                             v0,
                                                             v1,
                                                             element_stride,
                                                             elements,
                                                             tol,
                                                             max_depth,
                                                             d_toi,
                                                             g_stack,
                                                             (int)begin,
                                                             (int)end);
                    } else {
                        dim3 grid_pass1(this_batch, 1, 1);
                        narrow_phase_dfs_kernel<is_vf, SCCD_NP_THREADS_PER_BLOCK, T, I>
                            <<<grid_pass1, block_pass1>>>(overlap0,
                                                          overlap1,
                                                          v0,
                                                          v1,
                                                          element_stride,
                                                          elements,
                                                          tol,
                                                          max_depth,
                                                          d_toi,
                                                          toi_stride,
                                                          g_stack,
                                                          SCCD_NP_ALPHA,
                                                          (int)begin,
                                                          (int)end);
                    }
                    SCCD_CUDA_LAST_ERROR();

                    int h_g_top = 0;
                    SCCD_CHECK_CUDA(cudaMemcpy(&h_g_top, g_top, sizeof(int), cudaMemcpyDeviceToHost));

                    // Pass 2: drain whatever made it onto g_stack.  Each
                    // launch consumes up to base_grid_blocks worth of
                    // entries; relaunch until empty.  Spillover during
                    // drain is recorded in g_request and handled by the
                    // outer retry below.
                    while (h_g_top > 0) {
                        // printf("Draining g_stack with from-stack kernel (%d)\n", h_g_top);
                        int grid_blocks = (toi_stride == 0) ? (h_g_top + N - 1) / N : h_g_top;
                        if (grid_blocks > base_grid_blocks) grid_blocks = base_grid_blocks;
                        if (grid_blocks < 1) grid_blocks = 1;
                        dim3 grid_pass2(grid_blocks, 1, 1);

                        if (toi_stride == 0) {
                            narrow_phase_dfs_zero_stride_from_stack_kernel<is_vf, N, T, I><<<grid_pass2, block_pass1>>>(
                                overlap0, overlap1, v0, v1, element_stride, elements, tol, max_depth, d_toi, g_stack);
                        } else {
                            narrow_phase_dfs_from_stack_kernel<is_vf, N, T, I>
                                <<<grid_pass2, block_pass1>>>(overlap0,
                                                              overlap1,
                                                              v0,
                                                              v1,
                                                              element_stride,
                                                              elements,
                                                              tol,
                                                              max_depth,
                                                              d_toi,
                                                              toi_stride,
                                                              g_stack);
                        }
                        SCCD_CUDA_LAST_ERROR();

                        cudaMemcpy(&h_g_top, g_top, sizeof(int), cudaMemcpyDeviceToHost);
                    }

                    // Did anything overflow during Pass 1 or Pass 2?
                    // Grow and retry the entire batch from seeds; the
                    // tighter TOIs from this attempt cut the work tree
                    // on the next pass.
                    int h_g_request = 0;
                    cudaMemcpy(&h_g_request, g_request, sizeof(int), cudaMemcpyDeviceToHost);
                    if (h_g_request <= 0) break;

                    // printf("Overflowed: %d\n", h_g_request);

                    int grow_by = h_g_request;
                    if (grow_by > SCCD_GSTACK_CAP_MAX) grow_by = SCCD_GSTACK_CAP_MAX;
                    const long long target_ll = (long long)gstack_cap + (long long)grow_by;
                    const int new_cap = (target_ll > (long long)INT_MAX) ? INT_MAX : (int)target_ll;
                    grow_stack(new_cap);
                }
            }

            cudaFree(g_tlower);
            cudaFree(g_tupper);
            cudaFree(g_ulower);
            cudaFree(g_uupper);
            cudaFree(g_vlower);
            cudaFree(g_vupper);
            cudaFree(g_level);
            cudaFree(g_qid);
            cudaFree(g_top);
            cudaFree(g_request);

            SCCD_CUDA_LAST_ERROR();
            return 0;
        }

        template <typename T, typename I>
        int narrow_phase_ee(const size_t noverlaps,
                            const I* const SCCD_RESTRICT overlap0,
                            const I* const SCCD_RESTRICT overlap1,
                            // Geometric data
                            T** const SCCD_RESTRICT v0,
                            T** const SCCD_RESTRICT v1,
                            const size_t edge_stride,
                            I** const SCCD_RESTRICT edges,
                            // Output
                            const T max_toi,
                            T* const SCCD_RESTRICT toi,
                            const int max_depth,
                            const T tol,
                            const int toi_stride) {
            return narrow_phase_generic<false, T, I>(
                noverlaps, overlap0, overlap1, v0, v1, edge_stride, edges, max_toi, toi, max_depth, tol, toi_stride);
        }

        template <int nxe, typename T, typename I>
        int narrow_phase_vf(const size_t noverlaps,
                            const I* const SCCD_RESTRICT voveralp,
                            const I* const SCCD_RESTRICT foveralp,
                            // Geometric data
                            T** const SCCD_RESTRICT v0,
                            T** const SCCD_RESTRICT v1,
                            const size_t face_stride,
                            I** const SCCD_RESTRICT faces,
                            const T max_toi,
                            T* const SCCD_RESTRICT toi,
                            const int max_depth,
                            const T tol,
                            const int toi_stride) {
            return narrow_phase_generic<true, T, I>(
                noverlaps, voveralp, foveralp, v0, v1, face_stride, faces, max_toi, toi, max_depth, tol, toi_stride);
        }

        template <typename T>
        int minmax(const T* const SCCD_RESTRICT data, const size_t n, T* const h_min, T* const h_max) {
            if (n == 0) {
                *h_min = std::numeric_limits<T>::max();
                *h_max = std::numeric_limits<T>::lowest();
                return 0;
            }

            auto begin = thrust::device_pointer_cast(data);
            auto minmax = thrust::minmax_element(begin, begin + n);
            SCCD_CHECK_CUDA(
                cudaMemcpy(h_min, thrust::raw_pointer_cast(minmax.first), sizeof(*h_min), cudaMemcpyDeviceToHost));
            SCCD_CHECK_CUDA(
                cudaMemcpy(h_max, thrust::raw_pointer_cast(minmax.second), sizeof(*h_max), cudaMemcpyDeviceToHost));
            return 0;
        }

    }  // namespace device
}  // namespace sccd

#define INSTANTIATE_NARROW_PHASE_EE(T, I)                                                   \
    template int sccd::device::narrow_phase_ee<T, I>(const size_t noverlaps,                \
                                                     const I* const SCCD_RESTRICT overlap0, \
                                                     const I* const SCCD_RESTRICT overlap1, \
                                                     T** const SCCD_RESTRICT v0,            \
                                                     T** const SCCD_RESTRICT v1,            \
                                                     const size_t element_stride,           \
                                                     I** const SCCD_RESTRICT elements,      \
                                                     const T max_toi,                       \
                                                     T* const SCCD_RESTRICT toi,            \
                                                     const int max_depth,                   \
                                                     const T tol,                           \
                                                     const int toi_stride);

#define INSTANTIATE_NARROW_PHASE_VF(NXE, T, I)                                                   \
    template int sccd::device::narrow_phase_vf<NXE, T, I>(const size_t noverlaps,                \
                                                          const I* const SCCD_RESTRICT voveralp, \
                                                          const I* const SCCD_RESTRICT foveralp, \
                                                          T** const SCCD_RESTRICT v0,            \
                                                          T** const SCCD_RESTRICT v1,            \
                                                          const size_t element_stride,           \
                                                          I** const SCCD_RESTRICT elements,      \
                                                          const T max_toi,                       \
                                                          T* const SCCD_RESTRICT toi,            \
                                                          const int max_depth,                   \
                                                          const T tol,                           \
                                                          const int toi_stride);

INSTANTIATE_NARROW_PHASE_EE(float, int32_t);
INSTANTIATE_NARROW_PHASE_EE(float, int64_t);
INSTANTIATE_NARROW_PHASE_EE(double, int32_t);
INSTANTIATE_NARROW_PHASE_EE(double, int64_t);

INSTANTIATE_NARROW_PHASE_VF(3, float, int32_t);
INSTANTIATE_NARROW_PHASE_VF(3, float, int64_t);
INSTANTIATE_NARROW_PHASE_VF(3, double, int32_t);
INSTANTIATE_NARROW_PHASE_VF(3, double, int64_t);

template int sccd::device::minmax<float>(const float* const SCCD_RESTRICT data,
                                         const size_t n,
                                         float* const h_min,
                                         float* const h_max);
template int sccd::device::minmax<double>(const double* const SCCD_RESTRICT data,
                                          const size_t n,
                                          double* const h_min,
                                          double* const h_max);

#undef INSTANTIATE_NARROW_PHASE_EE
#undef INSTANTIATE_NARROW_PHASE_VF
