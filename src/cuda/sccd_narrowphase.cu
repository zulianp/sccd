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

// Compile-time configuration of the warp-per-query narrowphase EE kernel.
#ifndef SCCD_NP_WARPS_PER_BLOCK
#define SCCD_NP_WARPS_PER_BLOCK 4
#endif

#ifndef SCCD_NP_SHARED_STACK_CAP
#define SCCD_NP_SHARED_STACK_CAP 2048
#endif

#ifndef SCCD_NP_THREADS_PER_BLOCK
#define SCCD_NP_THREADS_PER_BLOCK 128
#endif

// Maximum number of bisections (per-thread DFS depth) before a candidate
// subdomain is "accepted" with cur.tlower as a conservative TOI for that
// cell.  This bounds the depth of the depth-first traversal in
// narrow_phase_dfs_kernel and prevents pathologically deep recursion.
#ifndef SCCD_NP_MAX_BISECTIONS
#define SCCD_NP_MAX_BISECTIONS 69
// #define SCCD_NP_MAX_BISECTIONS 128
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

        // A LIFO of subdomains.  `top` is a POINTER to the counter so the
        // same struct can represent either a shared-memory stack (top
        // lives in shared) or a device-global stack (top lives in
        // global), and atomics on *top hit the appropriate memory
        // scope.  All array fields use an SoA layout for coalesced
        // access during push/pop by peer warps.
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
        static inline __device__ void compute_face_vertex_tolerance_soa(const T codomain_tol,
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

            // f(t,u,v) = position_on_edge_a(t,u) - position_on_edge_b(t,v)
            // Vec4 layout: (.x,.y) = edge a vertices (P0_a, P1_a)
            //              (.z,.w) = edge b vertices (P0_b, P1_b).
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
            // Gate the flags on the origin-containment test in this axis.
            // Without MASK_FULL the bitwise-AND would clip bits 1..3 (see
            // operator precedence: the ternary binds tighter than &).
            cond_mask &= ((fmin <= tol) && (fmax >= -tol)) ? MASK_FULL : 0;
            return cond_mask;
        }

        // ----------------------------------------------------------------
        // Pure cell evaluation: tests an already-built (t, u, v) domain
        // for origin containment and the four accept conditions.  Use
        // this when the caller already has the cell box in hand (e.g.
        // after a bisection) and no subdivision is needed.
        // ----------------------------------------------------------------
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

        // ----------------------------------------------------------------
        // Per-cell evaluation with subdivision: each lane processes one
        // (ti, ui, vi) cell of a 2 x 4 x 4 decomposition (or any other
        // nt x nu x nv) of the parent (t, u, v) subdomain.  The built
        // sub-cell is returned via `cell`; the rest delegates to
        // evaluate_cell_3d.
        // ----------------------------------------------------------------
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

        // Load the 4 endpoint coordinates (sx,sy,sz,ex,ey,ez) for a
        // given query (overlap candidate) into registers.
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
            // TODO: Implement this
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
                compute_face_vertex_tolerance_soa<T, Vec4>(tol, sx, sy, sz, ex, ey, ez, &atol[0], &atol[1], &atol[2]);
            } else {
                load_query_ee<T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, sx, sy, sz, ex, ey, ez);
                compute_edge_edge_tolerance_soa<T, Vec4>(tol, sx, sy, sz, ex, ey, ez, &atol[0], &atol[1], &atol[2]);
            }
        }

        // Backoff used when both stacks are empty but other warps are
        // still in-flight.  Falls back to a busy spin on pre-Volta GPUs.
        static inline __device__ void backoff_ns(const unsigned int ns) {
#if (__CUDA_ARCH__ >= 700)
            __nanosleep(ns);
#else
            (void)ns;
#endif
        }

        // Per-slot validity tag stored in the qid array.  Three states:
        //   SCCD_QID_EMPTY   (-1): free, available to be claimed by a writer
        //   SCCD_QID_WRITING (-2): claimed by a writer, fields being filled
        //   qid >= 0              : committed, fields are safe to read
        //
        // Writer:
        //   1. reserve_slots bumps `top` (exclusive ownership of the slot index)
        //   2. CAS qid[slot] EMPTY -> WRITING (spin-wait until the previous
        //      popper, if any, clears the slot)
        //   3. write fields, fence, atomicExch qid[slot] -> real qid
        //
        // Popper:
        //   1. release_slot decrements `top`
        //   2. spin until qid[slot] >= 0 (committed)
        //   3. fence, read fields, atomicExch qid[slot] -> EMPTY
        //
        // The WRITING state is essential: without it, a pusher that reserves
        // a slot just after a popper decremented `top` would clobber the
        // still-in-use fields the popper is about to read.
        static constexpr int SCCD_QID_EMPTY = -1;
        static constexpr int SCCD_QID_WRITING = -2;

        // ----------------------------------------------------------------
        // One-time init
        //  - toi[i] = max_toi for i in [0, noverlaps)
        //  - g_qid[i] = SCCD_QID_EMPTY for i in [0, g_capacity)
        //
        // The old seed kernel used to push a full-domain entry per query
        // onto the global LIFO, which persistent workers would then pop
        // back out.  That is pure bookkeeping: we know a priori that
        // every seed is the full unit cube at level 0.  Workers now pull
        // seeds lazily via an atomic counter (see seed_cursor), so the
        // global stack starts empty and we avoid O(noverlaps) wasted
        // writes-then-reads.
        // ----------------------------------------------------------------
        template <typename T>
        __global__ void init_narrow_phase_kernel(const size_t toi_n,
                                                 const int g_capacity,
                                                 const T max_toi,
                                                 T* SCCD_RESTRICT toi,
                                                 int* SCCD_RESTRICT g_qid) {
            const size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
            if (i < toi_n) toi[i] = max_toi;
            if (i < (size_t)g_capacity) g_qid[i] = SCCD_QID_EMPTY;
        }

        // Per-batch reset.  Launched as <<<1,1>>> between batches; much
        // cheaper than a memset+memcpy dance.  Stream-ordered against
        // the preceding worker launch, so no explicit sync is needed.
        __global__ void reset_batch_narrow_phase_kernel(int* SCCD_RESTRICT g_top,
                                                        int* SCCD_RESTRICT inflight,
                                                        int* SCCD_RESTRICT seed_cursor,
                                                        int* SCCD_RESTRICT halt,
                                                        const int seed_begin) {
            *g_top = 0;
            *inflight = 0;
            *seed_cursor = seed_begin;
            *halt = 0;
        }

        // Reserve k contiguous slots in a stack via bounded CAS.
        // Returns the base index, or -1 when there is not enough room.
        // Bounded so the counter never exceeds 'capacity'.
        static inline __device__ int reserve_slots(int* SCCD_RESTRICT counter, int k, int capacity) {
            int old = atomicAdd(counter, 0);
            while (true) {
                if (old + k > capacity) return -1;
                const int prev = atomicCAS(counter, old, old + k);
                if (prev == old) return old;
                old = prev;
            }
        }

        // Pop one slot from a stack via bounded CAS.
        // Returns the slot index to read (= old_top - 1), or -1 if empty.
        // Bounded so the counter never goes below 0.
        static inline __device__ int release_slot(int* SCCD_RESTRICT counter) {
            int old = atomicAdd(counter, 0);
            while (true) {
                if (old <= 0) return -1;
                const int prev = atomicCAS(counter, old, old - 1);
                if (prev == old) return old - 1;
                old = prev;
            }
        }

        // Pop one entry from a stack into Domain + level + qid.
        // `Shared` selects the memory scope of the fence that orders
        // the validity-tag read against the field reads.  Both variants
        // share identical logic otherwise; template dispatch lets the
        // compiler pick the cheaper `__threadfence_block` for the
        // intra-block case.
        template <bool Shared, typename T>
        static inline __device__ int try_pop(const Stack<T>& stk, Domain<T>& d, int& level, int& qid) {
            const int slot = release_slot(stk.top);
            if (slot < 0) return 0;
            // Spin until the writer publishes a committed (non-negative)
            // qid.  EMPTY (-1) and WRITING (-2) both mean "not yet
            // published"; any q >= 0 means the writer's field writes
            // retired before the qid publish (they fenced first).
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
            // Reset valid tag so the slot is reusable.
            atomicExch(&stk.qid[slot], SCCD_QID_EMPTY);
            return 1;
        }

        // Warp-cooperative push.  Each lane holds one candidate cell
        // (its own slice of the parent's 2 x 4 x 4 subdivision); the
        // `keep` bit indicates whether this lane's cell should be
        // enqueued for further subdivision.  All keeping lanes push
        // in a single contiguous reservation, written in parallel.
        //
        // Memory discipline (same as before): the global stack is split
        // into [0, g_normal_cap) for normal pushes and
        // [g_normal_cap, g_capacity) for halt-flush.  If both the
        // shared stack and the normal global zone can't fit the group
        // of `n = __popc(keep_mask)` children, halt is raised and the
        // push goes into the reserve.  Host sizing (R) guarantees the
        // reserve can always absorb at most 32 children per warp plus
        // the block-flush.
        template <typename T>
        static inline __device__ void push_warp_cells(const Stack<T>& s_stack,
                                                      const Stack<T>& g_stack,
                                                      const int g_normal_cap,
                                                      int* SCCD_RESTRICT halt,
                                                      const int lane,
                                                      const bool keep,
                                                      const Domain<T>& cell,
                                                      const int level,
                                                      const int qid) {
            const unsigned keep_mask = __ballot_sync(0xffffffffu, keep);
            const int n = __popc(keep_mask);
            if (n == 0) return;

            // Lane 0 reserves contiguous slots: shared first, then
            // global normal, then global reserve (under halt).
            int slot_base = -1;
            int is_shared = 0;
            if (lane == 0) {
                slot_base = reserve_slots(s_stack.top, n, s_stack.capacity);
                if (slot_base >= 0) {
                    is_shared = 1;
                } else {
                    slot_base = reserve_slots(g_stack.top, n, g_normal_cap);
                    if (slot_base < 0) {
                        atomicExch(halt, 1);
                        slot_base = reserve_slots(g_stack.top, n, g_stack.capacity);
                        assert(slot_base >= 0);
                    }
                }
            }
            slot_base = __shfl_sync(0xffffffffu, slot_base, 0);
            is_shared = __shfl_sync(0xffffffffu, is_shared, 0);
            // If even the reserve zone overflowed (host under-sized R, or
            // NDEBUG disabled the assert above), bail out before issuing
            // out-of-bounds writes.  halt is already raised so the host
            // will observe the failure and can relaunch.
            if (slot_base < 0) return;

            // Each keeping lane writes into slot_base + (its rank among
            // keepers).  `keep_mask & ((1u << lane) - 1u)` is the
            // prefix, popc gives the rank.  `is_shared` is uniform
            // across the warp (it came from a single-source shfl), so
            // the ternary below never diverges.
            const int rank = __popc(keep_mask & ((1u << lane) - 1u));
            const int slot = slot_base + rank;

            if (keep) {
                const Stack<T>& dst = is_shared ? s_stack : g_stack;
                // Claim the slot by transitioning qid EMPTY -> WRITING.
                // A concurrent popper that decremented `top` past this
                // index still holds an outstanding read on the old
                // fields; spin until it finishes (atomicExch to EMPTY)
                // before we overwrite them.
                while (atomicCAS(&dst.qid[slot], SCCD_QID_EMPTY, SCCD_QID_WRITING) != SCCD_QID_EMPTY) {
                    // busy-wait
                }
                dst.tlower[slot] = cell.tlower;
                dst.tupper[slot] = cell.tupper;
                dst.ulower[slot] = cell.ulower;
                dst.uupper[slot] = cell.uupper;
                dst.vlower[slot] = cell.vlower;
                dst.vupper[slot] = cell.vupper;
                dst.level[slot] = level;
            }
            // Make sure all field writes retire before any qid publish.
            __syncwarp();
            if (is_shared) {
                __threadfence_block();
                if (keep) atomicExch(&s_stack.qid[slot], qid);
                __threadfence_block();
            } else {
                __threadfence();
                if (keep) atomicExch(&g_stack.qid[slot], qid);
                __threadfence();
            }
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

        // Per-thread DFS variant for the toi_stride==0 case (one global toi
        // shared across all candidates).  Each thread owns one initial seed
        // and its own DFS register state; the block shares a stack used as
        // a load-balancing reservoir, with overflow spilling to the global
        // stack so other blocks can pick it up.
        //
        // Layout: blockDim.x == N, qid = seed_begin + blockIdx.x * N + tid
        // (static assignment).  Termination is best-effort: a block exits
        // when it observes its own threads idle, the shared stack empty,
        // and the global stack empty for several consecutive iterations.
        // Any work still in g_stack at exit is drained by the existing
        // warp-cooperative kernel from the host driver.
        template <bool is_vf, int N, typename T, typename I>
        __global__ void narrow_phase_dfs_zero_stride_kernel(const I* const SCCD_RESTRICT overlap0,
                                                            const I* const SCCD_RESTRICT overlap1,
                                                            T** const SCCD_RESTRICT sp,
                                                            T** const SCCD_RESTRICT ep,
                                                            const size_t element_stride,
                                                            I** const SCCD_RESTRICT elements,
                                                            const T tol,
                                                            //    One global toi for all candidates
                                                            T* SCCD_RESTRICT toi,
                                                            Stack<T> g_stack,
                                                            const int g_normal_cap,
                                                            int* SCCD_RESTRICT halt,
                                                            const T alpha,
                                                            const int seed_begin,
                                                            const int seed_end) {
            static_assert(N == 64 || N == 128 || N == 256, "SCCD_NP_THREADS_PER_BLOCK must be one of 64/128/256");
            (void)alpha;  // unused: this kernel does not hard-defer.
            using Vec4 = typename device::Vec4Type<T>::type;
            constexpr int S_CAP = SCCD_NP_SHARED_STACK_CAP;
            constexpr int max_bisections = SCCD_NP_MAX_BISECTIONS;
            // Number of consecutive iterations the block must observe an
            // empty workload before declaring termination.  Higher values
            // reduce the chance of exiting just before another block
            // pushes new work into g_stack.
            constexpr int MAX_IDLE_ITERS = 8;

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
            const int my_seed = seed_begin + (int)blockIdx.x * N + tid;
            const bool has_seed = my_seed < seed_end;

            // Init shared validity tags (only needed for global stack -- the
            // shared stack uses Phase A/B sync separation), s_top, s_toi.
            for (int i = tid; i < S_CAP; i += N) s_qid[i] = SCCD_QID_EMPTY;
            if (tid == 0) {
                s_top = 0;
                s_toi = toi[0];
            }
            __syncthreads();

            Stack<T> s_stack = {
                s_tlower, s_tupper, s_ulower, s_uupper, s_vlower, s_vupper, s_level, s_qid, &s_top, S_CAP};

            // Per-thread DFS state.  qid identifies the candidate currently
            // owning `cur`; it changes when the thread pops a slot tagged
            // with a different candidate, in which case geometry/atol must
            // be reloaded.
            int qid = -1;
            Domain<T> cur = {T(0), T(0), T(0), T(0), T(0), T(0)};
            int level = 0;
            int active = 0;
            Vec4 sx, sy, sz, ex, ey, ez;
            T atol[3] = {T(0), T(0), T(0)};

            // Initial seed: load geometry and use the root only as a
            // refinement seed.  Do not accept the full [0,1]^3 root domain:
            // its tlower is exactly 0 and creates a false zero TOI.
            if (has_seed) {
                qid = my_seed;
                load_query_and_tol<is_vf, T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol);

                Domain<T> root = {T(0), T(1), T(0), T(1), T(0), T(1)};
                int contains = 0;
                int accept = 0;
                evaluate_cell_3d<is_vf, T, Vec4>(root, sx, sy, sz, ex, ey, ez, tol, atol, contains, accept);
                if (contains && is_domain_valid<is_vf>(root, s_toi, atol)) {
                    cur = root;
                    level = 0;
                    active = 1;
                }
            }

            int idle_iters = 0;

            while (true) {
                // Cross-block prune: fold s_toi into the global toi and
                // refresh s_toi with whichever side is tighter.  The block
                // can then prune dominated cells before any further work.
                if (tid == 0) {
                    const T g = device::atomic_min(&toi[0], s_toi);
                    if (g < s_toi) s_toi = g;
                }
                __syncthreads();

                if (active && cur.tlower >= s_toi) active = 0;

                // Depth cap: cur is origin-containing by construction; at
                // max depth take its tlower as a conservative TOI for this
                // cell instead of bisecting further.
                if (active && level >= max_bisections) {
                    device::atomic_min(&s_toi, cur.tlower);
                    active = 0;
                }

                // ----------- Phase A: bisect & decide what to push -----------
                Domain<T> push_box;
                int push_level = 0;
                int will_push = 0;

                if (active) {
                    Domain<T> left, right;
                    bisect_longest_axis<T>(cur, atol, left, right);
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

                // Push: try shared first (single-writer slot, no tag needed
                // because the __syncthreads below separates writers from
                // readers).  Overflow into global uses the EMPTY/WRITING
                // tag protocol because cross-block pop/push interleave.
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
                        int g_slot = reserve_slots(g_stack.top, 1, g_normal_cap);
                        if (g_slot < 0) {
                            atomicExch(halt, 1);
                            g_slot = reserve_slots(g_stack.top, 1, g_stack.capacity);
                        }
                        if (g_slot >= 0 && g_slot < g_stack.capacity) {
                            while (atomicCAS(&g_stack.qid[g_slot], SCCD_QID_EMPTY, SCCD_QID_WRITING) !=
                                   SCCD_QID_EMPTY) {
                                // Slot still being read by a popper that
                                // released the index but has not yet
                                // exchanged the qid back to EMPTY.
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

                // ----------- Phase B: idle threads pull more work -----------
                if (!active) {
                    // Shared pop: bounded CAS decrement, then plain reads.
                    // Pushes from Phase A are already visible across the
                    // barrier above.
                    int t = atomicAdd(&s_top, 0);
                    int got = 0;
                    while (t > 0 && !got) {
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
                                got = 1;
                            }
                            break;
                        }
                        t = prev;
                    }

                    // Global pop: tag protocol.  Spin until the writer
                    // commits a non-negative qid, then drain fields and
                    // recycle the slot back to EMPTY.
                    if (!got) {
                        const int g_slot = release_slot(g_stack.top);
                        if (g_slot >= 0 && g_slot < g_stack.capacity) {
                            int q;
                            do {
                                q = atomicAdd(&g_stack.qid[g_slot], 0);
                            } while (q < 0);
                            __threadfence();
                            cur.tlower = g_stack.tlower[g_slot];
                            cur.tupper = g_stack.tupper[g_slot];
                            cur.ulower = g_stack.ulower[g_slot];
                            cur.uupper = g_stack.uupper[g_slot];
                            cur.vlower = g_stack.vlower[g_slot];
                            cur.vupper = g_stack.vupper[g_slot];
                            level = g_stack.level[g_slot];
                            atomicExch(&g_stack.qid[g_slot], SCCD_QID_EMPTY);
                            if (q != qid) {
                                qid = q;
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
                    }
                }

                __syncthreads();

                // Termination (best-effort): block exits once it has seen
                // no active thread and both stacks empty for MAX_IDLE_ITERS
                // consecutive iterations.  A halt flag also forces exit so
                // the host can react to a full g_stack.
                const int n_active = block_popc<N>(active, warp_sums);
                if (n_active == 0) {
                    const int s_now = atomicAdd(&s_top, 0);
                    const int g_now = atomicAdd(g_stack.top, 0);
                    const int h_now = atomicAdd(halt, 0);
                    if ((s_now <= 0 && g_now <= 0) || h_now != 0) {
                        idle_iters += 1;
                        if (idle_iters >= MAX_IDLE_ITERS) break;
                    } else {
                        idle_iters = 0;
                    }
                } else {
                    idle_iters = 0;
                }
            }

            // Final flush: ensure this block's tightest s_toi is folded
            // into the global toi before exit.
            if (tid == 0) device::atomic_min(&toi[0], s_toi);
        }

        template <bool is_vf, int N, typename T, typename I>
        __global__ void narrow_phase_dfs_kernel(const I* const SCCD_RESTRICT overlap0,
                                                const I* const SCCD_RESTRICT overlap1,
                                                T** const SCCD_RESTRICT sp,
                                                T** const SCCD_RESTRICT ep,
                                                const size_t element_stride,
                                                I** const SCCD_RESTRICT elements,
                                                const T tol,
                                                T* SCCD_RESTRICT toi,
                                                const int toi_stride,
                                                Stack<T> g_stack,
                                                const int g_normal_cap,
                                                int* SCCD_RESTRICT halt,
                                                const T alpha,
                                                const int seed_begin,
                                                const int seed_end) {
            static_assert(N == 64 || N == 128 || N == 256, "SCCD_NP_THREADS_PER_BLOCK must be one of 64/128/256");
            using Vec4 = typename device::Vec4Type<T>::type;
            constexpr int NT = DfsSplit<N>::NT;
            constexpr int NU = DfsSplit<N>::NU;
            constexpr int NV = DfsSplit<N>::NV;
            constexpr int S_CAP = SCCD_NP_SHARED_STACK_CAP;

            const int tid = threadIdx.x;
            const int qid = seed_begin + (int)blockIdx.x;
            if (qid >= seed_end) return;

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

            Stack<T> s_stack = {
                s_tlower, s_tupper, s_ulower, s_uupper, s_vlower, s_vupper, s_level, s_qid, &s_top, S_CAP};

            T atol[3];
            Vec4 sx, sy, sz, ex, ey, ez;

            load_query_and_tol<is_vf, T, Vec4, I>(
                qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol);

            const int ti = tid / (NU * NV);
            const int rem = tid % (NU * NV);
            const int ui = rem / NV;
            const int vi = rem % NV;

            Domain<T> root = {T(0), T(1), T(0), T(1), T(0), T(1)};
            Domain<T> cur;
            int contains = 0;
            int accept = 0;
            sample_cell_3d<is_vf, T, Vec4>(
                ti, ui, vi, NT, NU, NV, root, sx, sy, sz, ex, ey, ez, tol, atol, contains, accept, cur);

            // Fold accepted cells into s_toi first, so the prune below can
            // use the tightest bound available within this block.
            if (accept && is_domain_valid<is_vf>(cur, s_toi, atol)) {
                device::atomic_min(&s_toi, cur.tlower);
            }
            __syncthreads();

            // Drop any seed cell whose time interval is already dominated
            // by s_toi: the smallest TOI it could ever produce is tlower,
            // which is >= the current best, so further bisection is wasted.
            const int active_seed = contains && !accept && (cur.tlower < s_toi);

            const int co_count = block_popc<N>(active_seed, warp_sums);
            if (!co_count) {
                // Block has no work left, but s_toi may have been improved
                // by accepts in the initial seed -- flush it to the global
                // toi before exiting so the resu   lt is not lost.
                if (tid == 0) device::atomic_min(&toi[toi_idx], s_toi);
                return;
            }

            if (tid == 0) {
                if ((T)co_count > alpha * (T)N) {
                    s_hard = 1;
                    int base = reserve_slots(g_stack.top, co_count, g_normal_cap);
                    if (base < 0) {
                        atomicExch(halt, 1);
                        base = reserve_slots(g_stack.top, co_count, g_stack.capacity);
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
                        g_stack.tupper[slot] = cur.tupper;
                        g_stack.ulower[slot] = cur.ulower;
                        g_stack.uupper[slot] = cur.uupper;
                        g_stack.vlower[slot] = cur.vlower;
                        g_stack.vupper[slot] = cur.vupper;
                        g_stack.level[slot] = 0;
                        g_stack.qid[slot] = qid;
                    }
                }
                if (tid == 0) device::atomic_min(&toi[toi_idx], s_toi);
                return;
            }

            constexpr int max_bisections = SCCD_NP_MAX_BISECTIONS;

            int active = active_seed;
            int level = 1;
            while (true) {
                // Cross-block pruning: when the toi is shared across all
                // candidates (toi_stride == 0), fold s_toi into the
                // global toi and refresh s_toi with whichever side was
                // tighter.  device::atomic_min returns the value of
                // *address before the (potential) update, so the value
                // returned is min(prev_global, s_toi) modulo concurrent
                // updates, and writing it back to s_toi keeps the block
                // in sync with the global minimum.
                if (toi_stride == 0) {
                    if (tid == 0) {
                        const T g = device::atomic_min(&toi[toi_idx], s_toi);
                        if (g < s_toi) s_toi = g;
                    }
                    __syncthreads();
                }

                if (active && cur.tlower >= s_toi) active = 0;

                // Depth cap: cur is origin-containing by construction (we
                // only keep origin-containing halves during DFS).  At max
                // depth, accept its tlower as a conservative TOI for this
                // cell instead of bisecting further.
                if (active && level >= max_bisections) {
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
                    bisect_longest_axis<T>(cur, atol, left, right);
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

                // Phase A: pushes. Bounded-CAS reserve, no spin against
                // any other thread's progress. Overflow spills to global
                // stack (Pass 2 will pick those up).  No qid-tag protocol
                // is needed: each reserved slot is touched by a single
                // thread, and the producer/consumer ordering is provided
                // by the __syncthreads at the end of this iteration and
                // the kernel-launch boundary between Pass 1 and Pass 2.
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
                        int g_slot = reserve_slots(g_stack.top, 1, g_normal_cap);
                        if (g_slot < 0) {
                            atomicExch(halt, 1);
                            g_slot = reserve_slots(g_stack.top, 1, g_stack.capacity);
                        }
                        if (g_slot >= 0 && g_slot < g_stack.capacity) {
                            g_stack.tlower[g_slot] = push_box.tlower;
                            g_stack.tupper[g_slot] = push_box.tupper;
                            g_stack.ulower[g_slot] = push_box.ulower;
                            g_stack.uupper[g_slot] = push_box.uupper;
                            g_stack.vlower[g_slot] = push_box.vlower;
                            g_stack.vupper[g_slot] = push_box.vupper;
                            g_stack.level[g_slot] = push_level;
                            g_stack.qid[g_slot] = qid;
                        }
                    }
                }

                __syncthreads();

                // Phase B: pops. Bounded-CAS decrement, no spin against
                // any other thread's progress. Each popper gets a unique
                // slot in [0, s_top); pushes from Phase A are visible.
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
                // if (!tid) printf("n_active: %d, s_top: %d\n", n_active, s_top);
                if (n_active == 0 && s_top == 0) break;
            }

            if (tid == 0) device::atomic_min(&toi[toi_idx], s_toi);
        }

        // ----------------------------------------------------------------
        // narrow_phase_kernel (warp-per-query, persistent worker)
        //
        // Block layout: (32, W, 1).  Each warp consumes one subdomain at
        // a time; lane 0 manages stack ops on behalf of the warp.  All
        // 32 lanes evaluate the 2 x 4 x 4 decomposition (one cell per
        // lane) in parallel via sample_cell_3d.
        //
        // Work sources, tried in order when lane 0 is hunting for work:
        //   1. block-local shared stack (hottest, cheapest)
        //   2. device-global LIFO stack (overflow from (1), cross-block)
        //   3. the seed cursor: an atomic counter in [seed_begin,
        //      seed_end).  Claiming seed i synthesizes the root entry
        //      (tl=0, tu=1, ul=0, uu=1, vl=0, vu=1, level=0, qid=i) in
        //      registers -- no global-memory dance.
        //
        // Termination: warp exits once (1) shared stack is empty, (2)
        // global stack is empty, (3) no warp is mid-processing
        // (inflight == 0), AND (4) seed_cursor >= seed_end.  Inflight
        // is incremented before any work-claim attempt to prevent
        // termination from racing ahead.
        // ----------------------------------------------------------------
        template <bool is_vf, typename T, typename I>
        __global__ void narrow_phase_kernel(const I* const SCCD_RESTRICT overlap0,
                                            const I* const SCCD_RESTRICT overlap1,
                                            T** const SCCD_RESTRICT sp,
                                            T** const SCCD_RESTRICT ep,
                                            const size_t element_stride,
                                            I** const SCCD_RESTRICT elements,
                                            const T tol,
                                            T* SCCD_RESTRICT toi,
                                            const int toi_stride,
                                            Stack<T> g_stack,
                                            const int g_normal_cap,
                                            int* SCCD_RESTRICT inflight,
                                            int* SCCD_RESTRICT seed_cursor,
                                            const int seed_end,
                                            int* SCCD_RESTRICT halt) {
            using Vec4 = typename device::Vec4Type<T>::type;

            constexpr int S_CAP = SCCD_NP_SHARED_STACK_CAP;
            constexpr int W = SCCD_NP_WARPS_PER_BLOCK;

            // Per-block shared stack storage.
            __shared__ T s_tlower[S_CAP];
            __shared__ T s_tupper[S_CAP];
            __shared__ T s_ulower[S_CAP];
            __shared__ T s_uupper[S_CAP];
            __shared__ T s_vlower[S_CAP];
            __shared__ T s_vupper[S_CAP];
            __shared__ int s_level[S_CAP];
            __shared__ int s_qid[S_CAP];
            __shared__ int s_top;

            // Per-warp slot used by lane 0 to publish the popped entry to
            // the rest of the warp.
            __shared__ int w_status[W];  // 0=work, 1=backoff, 2=terminate
            __shared__ Domain<T> w_domain[W];
            __shared__ int w_level[W];
            __shared__ int w_qid[W];

            const int lane = threadIdx.x;
            const int warp = threadIdx.y;
            const int tid_in_block = lane + warp * 32;
            const int threads_in_block = 32 * W;

            // Initialize shared-stack validity tags and counter, then
            // assemble a Stack<T> view over the shared arrays.
            for (int i = tid_in_block; i < S_CAP; i += threads_in_block) {
                s_qid[i] = SCCD_QID_EMPTY;
            }
            if (tid_in_block == 0) s_top = 0;
            __syncthreads();

            Stack<T> s_stack = {
                s_tlower, s_tupper, s_ulower, s_uupper, s_vlower, s_vupper, s_level, s_qid, &s_top, S_CAP};

            unsigned int backoff = 64;

            while (true) {
                if (lane == 0) {
                    // Fast-path abort if another warp raised the halt
                    // flag (global normal zone was full on a push).
                    // In-flight warps still finish their current
                    // iteration; they'll exit here next time around.
                    // The block-flush below drains the shared stack
                    // into global so the host can relaunch and pick
                    // up any residual work.
                    w_status[warp] = (atomicAdd(halt, 0) != 0) ? 2 : 0;
                }
                __syncwarp();
                if (w_status[warp] == 2) break;

                if (lane == 0) {
                    // Optimistically claim a worker slot before touching
                    // any stack counter so that termination cannot race
                    // ahead of an in-progress pop.
                    atomicAdd(inflight, 1);

                    Domain<T> d;
                    int level, qid;
                    int got = try_pop<true, T>(s_stack, d, level, qid);
                    if (!got) got = try_pop<false, T>(g_stack, d, level, qid);

                    // Source 3: lazy seed claim.  Synthesizes the root
                    // entry for query `old` in registers -- no global
                    // push/pop roundtrip.  A seed_cursor that overshoots
                    // seed_end is harmless: peers observing the larger
                    // counter simply read "no seeds left".
                    if (!got) {
                        const int old = atomicAdd(seed_cursor, 1);
                        if (old < seed_end) {
                            d = {T(0), T(1), T(0), T(1), T(0), T(1)};
                            level = 0;
                            qid = old;
                            got = 1;
                        }
                    }

                    if (got) {
                        w_domain[warp] = d;
                        w_level[warp] = level;
                        w_qid[warp] = qid;
                        w_status[warp] = 0;
                    } else {
                        // Release the optimistic claim and assess
                        // termination.  Ordering matters: only after
                        // releasing can other warps observe a true
                        // 'no inflight' state.
                        atomicSub(inflight, 1);
                        const int gtop_now = atomicAdd(g_stack.top, 0);
                        const int stop_now = atomicAdd(s_stack.top, 0);
                        const int infl_now = atomicAdd(inflight, 0);
                        const int seed_now = atomicAdd(seed_cursor, 0);
                        if (gtop_now <= 0 && stop_now <= 0 && infl_now <= 0 && seed_now >= seed_end) {
                            w_status[warp] = 2;
                        } else {
                            w_status[warp] = 1;
                        }
                    }
                }
                __syncwarp();

                const int status = w_status[warp];
                if (status == 2) break;
                if (status == 1) {
                    backoff_ns(backoff);
                    if (backoff < 4096) backoff <<= 1;
                    continue;
                }
                backoff = 64;

                const Domain<T> parent = w_domain[warp];
                const int level = w_level[warp];
                const int qid = w_qid[warp];

                // Early cull if a tighter toi has already been recorded
                // for this query.  Read once on lane 0 and broadcast so
                // all 32 lanes agree; otherwise a concurrent atomic_min
                // by another warp could make the plain load diverge
                // across lanes and break the __syncwarp() below.
                const int toi_idx = qid * toi_stride;

                T cur_toi = T(0);
                if (lane == 0) cur_toi = toi[toi_idx];
                cur_toi = __shfl_sync(0xffffffffu, cur_toi, 0);
                if (parent.tlower >= cur_toi) {
                    if (lane == 0) atomicSub(inflight, 1);
                    __syncwarp();
                    continue;
                }

                Vec4 sx, sy, sz, ex, ey, ez;
                T atol[3];

                load_query_and_tol<is_vf, T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol);

                // 2 x 4 x 4 lane decomposition: each lane evaluates one
                // cell of the parent subdomain.
                const int ti = lane >> 4;
                const int ui = (lane >> 2) & 0x3;
                const int vi = lane & 0x3;

                int contains = 0;
                int accept = 0;
                Domain<T> cell;
                sample_cell_3d<is_vf, T, Vec4>(
                    ti, ui, vi, 2, 4, 4, parent, sx, sy, sz, ex, ey, ez, tol, atol, contains, accept, cell);

                const bool accept_valid = accept && is_domain_valid<is_vf>(cell, cur_toi, atol);
                const unsigned acc_mask = __ballot_sync(0xffffffffu, accept_valid);

                if (acc_mask) {
                    // Reduce min(cell.tlower) over accepting lanes.  Use
                    // the parent's tupper as a sentinel so non-accepting
                    // lanes never win the reduction.
                    T best = accept_valid ? cell.tlower : parent.tupper;
                    for (int o = 16; o > 0; o >>= 1) {
                        T other = __shfl_xor_sync(0xffffffffu, best, o);
                        best = (best < other) ? best : other;
                    }
                    if (lane == 0) device::atomic_min(&toi[toi_idx], best);
                }

                // Refine every lane whose cell still contains the origin
                // AND has not yet been accepted AND is still competitive
                // with the current best toi[qid].  One 2x4x4 sampling
                // produces up to 32 children (one per lane), fully
                // using the work done by all 32 lanes in sample_cell_3d.
                T cur_toi_now = T(0);
                if (lane == 0) cur_toi_now = toi[toi_idx];
                cur_toi_now = __shfl_sync(0xffffffffu, cur_toi_now, 0);
                const bool keep_cell = contains && !accept && (cell.tlower < cur_toi_now);

                push_warp_cells<T>(s_stack, g_stack, g_normal_cap, halt, lane, keep_cell, cell, level + 1, qid);

                if (lane == 0) atomicSub(inflight, 1);
                __syncwarp();
            }

            // ------------------------------------------------------------
            // Block-flush: drain any residual shared-stack entries into
            // the global stack so the host can relaunch and finish them.
            //
            // Under normal termination s_top is 0 and this is a no-op;
            // under halt termination s_top may still hold valid entries
            // that must not be lost.  A single team-wide reserve on the
            // global stack amortizes the CAS cost; individual threads
            // then copy one entry each.
            // ------------------------------------------------------------
            __syncthreads();

            __shared__ int flush_n;
            __shared__ int flush_base;
            if (tid_in_block == 0) {
                flush_n = s_top;
                flush_base = (flush_n > 0) ? reserve_slots(g_stack.top, flush_n, g_stack.capacity) : 0;
                // By construction (host sizes g_capacity with reserve
                // R >= grid_blocks * (S_CAP + 32*W)) this must succeed.
                if (flush_n > 0) assert(flush_base >= 0);
            }
            __syncthreads();

            for (int i = tid_in_block; i < flush_n; i += threads_in_block) {
                const int src = i;
                const int dst = flush_base + i;
                // Same claim protocol as push_warp_cells: other blocks
                // may still be popping/pushing on g_stack while this
                // block flushes.
                while (atomicCAS(&g_stack.qid[dst], SCCD_QID_EMPTY, SCCD_QID_WRITING) != SCCD_QID_EMPTY) {
                    // busy-wait
                }
                g_stack.tlower[dst] = s_stack.tlower[src];
                g_stack.tupper[dst] = s_stack.tupper[src];
                g_stack.ulower[dst] = s_stack.ulower[src];
                g_stack.uupper[dst] = s_stack.uupper[src];
                g_stack.vlower[dst] = s_stack.vlower[src];
                g_stack.vupper[dst] = s_stack.vupper[src];
                g_stack.level[dst] = s_stack.level[src];
                // s_qid[src] must be a published (non-negative) tag
                // because the invariant is that slots in [0, s_top) are
                // fully committed.  Republish on the global side last,
                // after a fence, so a future popper sees fields first.
                const int q = s_stack.qid[src];
                __threadfence();
                atomicExch(&g_stack.qid[dst], q);
            }
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
                                 const int toi_stride) {
            SCCD_CUDA_LAST_ERROR();

            if (noverlaps == 0) return max_toi;
            assert(d_toi != nullptr);
            assert(toi_stride == 0 || toi_stride == 1);

            // toi length: 1 when stride==0 (all candidates share toi[0]),
            //             noverlaps when stride==1 (one toi per candidate).
            const size_t toi_n = (toi_stride == 0) ? 1 : noverlaps;

            T tol = std::is_same_v<T, float> ? T(1e-8) : T(1e-12);
            {
                double SCCD_TOL = (double)tol;
                SCCD_READ_ENV(SCCD_TOL, atof);
                tol = (T)SCCD_TOL;
            }
            T SCCD_NP_ALPHA = T(0.5);
            {
                double alpha = (double)SCCD_NP_ALPHA;
                SCCD_READ_ENV(alpha, atof);
                SCCD_NP_ALPHA = (T)alpha;
            }

            // ----------------------------------------------------------------
            // Auto-sized hyperparameters.
            //
            //   SCCD_BLOCKS_PER_SM  -> from CUDA occupancy API
            //   SCCD_GSTACK_CAP     -> from cudaMemGetInfo, capped by a
            //                          heuristic to avoid over-allocation
            //                          when noverlaps is small
            //   SCCD_BATCH_SIZE     -> from gstack headroom, so a single
            //                          launch is likely to finish a batch
            //                          without halt-flush
            //
            // Any of these can be overridden via the matching env var.
            // ----------------------------------------------------------------
            int dev = 0;
            cudaGetDevice(&dev);
            cudaDeviceProp prop{};
            cudaGetDeviceProperties(&prop, dev);

            const int W = SCCD_NP_WARPS_PER_BLOCK;
            const int block_threads = 32 * W;

            // Occupancy-derived blocks-per-SM (0 dynamic shared mem).
            int SCCD_BLOCKS_PER_SM = 4;
            {
                int occ = 0;
                if (cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                        &occ, (const void*)narrow_phase_kernel<is_vf, T, I>, block_threads, 0) == cudaSuccess &&
                    occ > 0) {
                    SCCD_BLOCKS_PER_SM = occ;
                }
            }
            SCCD_READ_ENV(SCCD_BLOCKS_PER_SM, atoi);

            int base_grid_blocks = prop.multiProcessorCount * SCCD_BLOCKS_PER_SM;
            if (base_grid_blocks <= 0) base_grid_blocks = 1;

            // Reserve count R sized so that at halt time every block can
            // flush its entire shared stack (S_CAP) and every warp can
            // complete an in-flight push of up to 32 children (one per
            // lane, from the 2 x 4 x 4 cell subdivision) without
            // running out of global-stack slots.
            const int reserve_R = base_grid_blocks * (SCCD_NP_SHARED_STACK_CAP + 32 * W);

            // Memory-derived default global-stack cap.  We budget a
            // fraction (default 0.25) of currently-free device memory
            // for the narrowphase scratch buffers, subtract the toi
            // array, and convert the remainder into stack slots.  The
            // effective cap is further clamped below by 2*R (needed
            // for the halt reserve) and above by a heuristic tied to
            // noverlaps so we don't reserve gigabytes for tiny inputs.
            size_t free_bytes = 0, total_bytes = 0;
            cudaMemGetInfo(&free_bytes, &total_bytes);

            double SCCD_MEM_FRACTION = 0.25;
            SCCD_READ_ENV(SCCD_MEM_FRACTION, atof);
            if (SCCD_MEM_FRACTION <= 0.0 || SCCD_MEM_FRACTION > 1.0) SCCD_MEM_FRACTION = 0.25;

            const size_t per_slot_bytes = 6 * sizeof(T) + 2 * sizeof(int);
            const size_t scratch_overhead = 4 * sizeof(int);  // counters
            const size_t budget = (size_t)((double)free_bytes * SCCD_MEM_FRACTION);

            size_t slots_from_mem = 0;
            if (budget > scratch_overhead) {
                slots_from_mem = (budget - scratch_overhead) / per_slot_bytes;
            }
            const size_t slots_heuristic = noverlaps * 32 + 4096;
            size_t auto_slots = (slots_from_mem < slots_heuristic) ? slots_from_mem : slots_heuristic;
            if (auto_slots < (size_t)(2 * reserve_R)) auto_slots = (size_t)(2 * reserve_R);
            if (auto_slots > (size_t)INT_MAX) auto_slots = (size_t)INT_MAX;

            int SCCD_GSTACK_CAP = (int)auto_slots;
            SCCD_READ_ENV(SCCD_GSTACK_CAP, atoi);
            int gstack_cap = SCCD_GSTACK_CAP;
            if (gstack_cap < 2 * reserve_R) gstack_cap = 2 * reserve_R;
            const int g_normal_cap = gstack_cap - reserve_R;

            // Candidates-per-kernel cap.  Default: sized so the normal
            // zone can (on average) hold a modest expansion per seed
            // without triggering a halt-flush.  We budget ~8 slots of
            // normal headroom per seed, so batch_size = normal_cap/8
            // -- clamped to [1, noverlaps].  Override via SCCD_BATCH_SIZE.
            int SCCD_BATCH_SIZE = 0;
            {
                size_t auto_batch = (size_t)g_normal_cap / 8;
                if (auto_batch < 1) auto_batch = 1;
                if (auto_batch > noverlaps) auto_batch = noverlaps;
                if (auto_batch > (size_t)INT_MAX) auto_batch = (size_t)INT_MAX;
                SCCD_BATCH_SIZE = (int)auto_batch;
            }
            SCCD_READ_ENV(SCCD_BATCH_SIZE, atoi);
            const size_t batch_size = (SCCD_BATCH_SIZE > 0) ? (size_t)SCCD_BATCH_SIZE : noverlaps;

            // Global LIFO stack arrays (shared across batches).
            T* g_tlower = nullptr;
            T* g_tupper = nullptr;
            T* g_ulower = nullptr;
            T* g_uupper = nullptr;
            T* g_vlower = nullptr;
            T* g_vupper = nullptr;
            int* g_level = nullptr;
            int* g_qid = nullptr;
            int* g_top = nullptr;
            int* inflight = nullptr;
            int* seed_cursor = nullptr;
            int* halt = nullptr;
            cudaMalloc(&g_tlower, gstack_cap * sizeof(T));
            cudaMalloc(&g_tupper, gstack_cap * sizeof(T));
            cudaMalloc(&g_ulower, gstack_cap * sizeof(T));
            cudaMalloc(&g_uupper, gstack_cap * sizeof(T));
            cudaMalloc(&g_vlower, gstack_cap * sizeof(T));
            cudaMalloc(&g_vupper, gstack_cap * sizeof(T));
            cudaMalloc(&g_level, gstack_cap * sizeof(int));
            cudaMalloc(&g_qid, gstack_cap * sizeof(int));
            cudaMalloc(&g_top, sizeof(int));
            cudaMalloc(&inflight, sizeof(int));
            cudaMalloc(&seed_cursor, sizeof(int));
            cudaMalloc(&halt, sizeof(int));

            // One-time init: toi[0..toi_n) = max_toi, g_qid[] = SCCD_QID_EMPTY.
            // Workers now pull seeds lazily, so the global stack is left
            // empty here (no per-query push-then-pop bookkeeping).
            {
                const int block = 256;
                // Compare in size_t; (int)noverlaps would be UB when
                // noverlaps > INT_MAX.
                const size_t work = ((size_t)gstack_cap > toi_n) ? (size_t)gstack_cap : toi_n;
                const size_t grid_sz = (work + block - 1) / block;
                const int grid = (grid_sz > (size_t)INT_MAX) ? INT_MAX : (int)grid_sz;
                init_narrow_phase_kernel<T><<<grid, block>>>(toi_n, gstack_cap, max_toi, d_toi, g_qid);
                SCCD_CUDA_LAST_ERROR();
            }

            dim3 block_pass2(32, SCCD_NP_WARPS_PER_BLOCK, 1);
            dim3 block_pass1(SCCD_NP_THREADS_PER_BLOCK, 1, 1);

            // Batch loop.  Each iteration claims [begin, end) of the
            // input.  Inside, we may relaunch the worker multiple times
            // until that range's work fully drains; this is how we cap
            // peak device memory -- when the global stack's normal zone
            // fills, workers halt + flush + exit and the host relaunches
            // with `g_top` and `seed_cursor` preserved.
            for (size_t begin = 0; begin < noverlaps; begin += batch_size) {
                const size_t end = (begin + batch_size < noverlaps) ? (begin + batch_size) : noverlaps;
                const int this_batch = (int)(end - begin);

                reset_batch_narrow_phase_kernel<<<1, 1>>>(g_top, inflight, seed_cursor, halt, (int)begin);
                SCCD_CUDA_LAST_ERROR();

                Stack<T> g_stack = {
                    g_tlower, g_tupper, g_ulower, g_uupper, g_vlower, g_vupper, g_level, g_qid, g_top, gstack_cap};

                // Pass 1 has two flavors selected by toi_stride:
                //   stride == 0 : per-thread DFS, one global toi shared
                //                 across candidates.  Grid is sized to
                //                 cover this_batch threads (N per block).
                //   stride == 1 : per-block DFS, one toi per candidate.
                //                 Grid is one block per candidate.
                if (toi_stride == 0) {
                    printf("Using narrow_phase_dfs_zero_stride_kernel (%d)\n", this_batch);
                    constexpr int N = SCCD_NP_THREADS_PER_BLOCK;
                    const int grid_blocks_zs = (this_batch + N - 1) / N;
                    dim3 grid_pass1_zs(grid_blocks_zs, 1, 1);
                    narrow_phase_dfs_zero_stride_kernel<is_vf, N, T, I><<<grid_pass1_zs, block_pass1>>>(overlap0,
                                                                                                        overlap1,
                                                                                                        v0,
                                                                                                        v1,
                                                                                                        element_stride,
                                                                                                        elements,
                                                                                                        tol,
                                                                                                        d_toi,
                                                                                                        g_stack,
                                                                                                        g_normal_cap,
                                                                                                        halt,
                                                                                                        SCCD_NP_ALPHA,
                                                                                                        (int)begin,
                                                                                                        (int)end);
                } else {
                    printf("Using narrow_phase_dfs_kernel (%d)\n", this_batch);
                    dim3 grid_pass1(this_batch, 1, 1);
                    narrow_phase_dfs_kernel<is_vf, SCCD_NP_THREADS_PER_BLOCK, T, I>
                        <<<grid_pass1, block_pass1>>>(overlap0,
                                                      overlap1,
                                                      v0,
                                                      v1,
                                                      element_stride,
                                                      elements,
                                                      tol,
                                                      d_toi,
                                                      toi_stride,
                                                      g_stack,
                                                      g_normal_cap,
                                                      halt,
                                                      SCCD_NP_ALPHA,
                                                      (int)begin,
                                                      (int)end);
                }
                SCCD_CUDA_LAST_ERROR();

                int h_g_top = 0;
                SCCD_CHECK_CUDA(cudaMemcpy(&h_g_top, g_top, sizeof(int), cudaMemcpyDeviceToHost));
                if (h_g_top <= 0) continue;

                printf("gtop = %d\n", h_g_top);

                const int pass2_seed = (int)end;
                cudaMemcpy(seed_cursor, &pass2_seed, sizeof(int), cudaMemcpyHostToDevice);
                cudaMemsetAsync(inflight, 0, sizeof(int));
                cudaMemsetAsync(halt, 0, sizeof(int));

                while (true) {
                    int grid_blocks = base_grid_blocks;
                    if (grid_blocks > this_batch) grid_blocks = this_batch;
                    dim3 grid_pass2(grid_blocks, 1, 1);

                    narrow_phase_kernel<is_vf, T, I><<<grid_pass2, block_pass2>>>(overlap0,
                                                                                  overlap1,
                                                                                  v0,
                                                                                  v1,
                                                                                  element_stride,
                                                                                  elements,
                                                                                  tol,
                                                                                  d_toi,
                                                                                  toi_stride,
                                                                                  g_stack,
                                                                                  g_normal_cap,
                                                                                  inflight,
                                                                                  seed_cursor,
                                                                                  (int)end,
                                                                                  halt);
                    SCCD_CUDA_LAST_ERROR();

                    // Check whether the batch is fully drained.  A D2H
                    // copy per relaunch is cheap compared to the kernel
                    // runtime and is only paid once per overflow event.
                    int h_seed_cursor = 0;
                    cudaMemcpy(&h_g_top, g_top, sizeof(int), cudaMemcpyDeviceToHost);
                    cudaMemcpy(&h_seed_cursor, seed_cursor, sizeof(int), cudaMemcpyDeviceToHost);
                    if (h_g_top <= 0 && h_seed_cursor >= (int)end) break;

                    // Work remains; reset transient counters and relaunch.
                    // g_top, g_qid[], seed_cursor, d_toi are preserved.
                    cudaMemsetAsync(halt, 0, sizeof(int));
                    cudaMemsetAsync(inflight, 0, sizeof(int));
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
            cudaFree(inflight);
            cudaFree(seed_cursor);
            cudaFree(halt);

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
                            const int toi_stride) {
            return narrow_phase_generic<false, T, I>(
                noverlaps, overlap0, overlap1, v0, v1, edge_stride, edges, max_toi, toi, toi_stride);
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
                            const int toi_stride) {
            return narrow_phase_generic<true, T, I>(
                noverlaps, voveralp, foveralp, v0, v1, face_stride, faces, max_toi, toi, toi_stride);
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
