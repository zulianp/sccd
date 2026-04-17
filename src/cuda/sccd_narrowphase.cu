#include "sccd_narrowphase.cuh"

#include <stdint.h>
#include <stdlib.h>

#include <cassert>

#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>

#include "sccd_cuda_base.cuh"
#include "sccd_reduce.cuh"

// Compile-time configuration of the warp-per-query narrowphase EE kernel.
#ifndef SCCD_NP_WARPS_PER_BLOCK
#define SCCD_NP_WARPS_PER_BLOCK 4
#endif

#ifndef SCCD_NP_SHARED_STACK_CAP
#define SCCD_NP_SHARED_STACK_CAP 128
#endif

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

        // ----------------------------------------------------------------
        // Per-cell evaluation: each lane processes one (ti, ui, vi) cell
        // of a 2 x 4 x 4 decomposition (or any other nt x nu x nv) of the
        // current (t, u, v) subdomain and returns whether the cell
        // contains the origin and whether it satisfies the accept
        // termination criteria.
        // ----------------------------------------------------------------
        template <typename T, typename Vec4>
        static inline __device__ void sample_cell_3d(const int ti,
                                                     const int ui,
                                                     const int vi,
                                                     const int nt,
                                                     const int nu,
                                                     const int nv,
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
                                                     int& contains_origin,
                                                     int& accept,
                                                     T& cell_tlower) {
            const T t_h = (tupper - tlower) / nt;
            const T u_h = (uupper - ulower) / nu;
            const T v_h = (vupper - vlower) / nv;

            const T tl = tlower + ti * t_h;
            const T tu = tl + t_h;
            const T ul = ulower + ui * u_h;
            const T uu = ul + u_h;
            const T vl = vlower + vi * v_h;
            const T vu = vl + v_h;

            cell_tlower = tl;

            T fmin, fmax;
            T f[8];

            sample_f_ee<T, Vec4>(tl, tu, ul, uu, vl, vu, sx, ex, f);
            fminmax<T>(f, fmin, fmax);
            const uint8_t x_mask = cond_mask<T>(fmin, fmax, tol, adaptive_tol[0]);
            int co = (fmin <= tol) & (fmax >= -tol);

            sample_f_ee<T, Vec4>(tl, tu, ul, uu, vl, vu, sy, ey, f);
            fminmax<T>(f, fmin, fmax);
            const uint8_t y_mask = cond_mask<T>(fmin, fmax, tol, adaptive_tol[1]);
            co &= (fmin <= tol) & (fmax >= -tol);

            sample_f_ee<T, Vec4>(tl, tu, ul, uu, vl, vu, sz, ez, f);
            fminmax<T>(f, fmin, fmax);
            const uint8_t z_mask = cond_mask<T>(fmin, fmax, tol, adaptive_tol[2]);
            co &= (fmin <= tol) & (fmax >= -tol);

            const uint8_t and_mask = (x_mask & y_mask & z_mask);
            const uint8_t or_mask = (x_mask | y_mask | z_mask);

            const bool cond1 = and_mask & MASK_DOMAIN_SMALLER_THAN_TOL;
            const bool cond2 = and_mask & MASK_BOX_INSIDE_EPSILON_BOX;
            const bool cond3 = or_mask & MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE;
            const bool cond4 = and_mask & MASK_INTERVAL_TERMINAL;

            contains_origin = co;
            accept = co && (cond1 || cond2 || cond3 || cond4);
        }

        // Load the 4 endpoint coordinates (sx,sy,sz,ex,ey,ez) for a
        // given query (overlap candidate) into registers.
        template <typename T, typename Vec4, typename I>
        static inline __device__ void load_query_ee(const int qid,
                                                    const I* const SCCD_RESTRICT e0overlap,
                                                    const I* const SCCD_RESTRICT e1overlap,
                                                    T** const SCCD_RESTRICT sp,
                                                    T** const SCCD_RESTRICT ep,
                                                    const size_t edge_stride,
                                                    I** const SCCD_RESTRICT edges,
                                                    Vec4& sx,
                                                    Vec4& sy,
                                                    Vec4& sz,
                                                    Vec4& ex,
                                                    Vec4& ey,
                                                    Vec4& ez) {
            const I ea = e0overlap[qid];
            const I eb = e1overlap[qid];

            const ptrdiff_t idxa0 = edges[0][ea * edge_stride];
            const ptrdiff_t idxa1 = edges[1][ea * edge_stride];
            const ptrdiff_t idxb0 = edges[0][eb * edge_stride];
            const ptrdiff_t idxb1 = edges[1][eb * edge_stride];

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

        // Backoff used when both stacks are empty but other warps are
        // still in-flight.  Falls back to a busy spin on pre-Volta GPUs.
        static inline __device__ void backoff_ns(const unsigned int ns) {
#if (__CUDA_ARCH__ >= 700)
            __nanosleep(ns);
#else
            (void)ns;
#endif
        }

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
        __global__ void init_narrow_phase_ee_kernel(const size_t noverlaps,
                                                    const int g_capacity,
                                                    const T max_toi,
                                                    T* SCCD_RESTRICT toi,
                                                    int* SCCD_RESTRICT g_qid) {
            const size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
            if (i < noverlaps) toi[i] = max_toi;
            if (i < (size_t)g_capacity) g_qid[i] = SCCD_QID_EMPTY;
        }

        // Per-batch reset.  Launched as <<<1,1>>> between batches; much
        // cheaper than a memset+memcpy dance.  Stream-ordered against
        // the preceding worker launch, so no explicit sync is needed.
        __global__ void reset_batch_narrow_phase_ee_kernel(int* SCCD_RESTRICT g_top,
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

        // We use the qid array of each stack as a per-slot "valid" tag:
        //   qid == SCCD_QID_EMPTY  --> slot empty / being written
        //   qid != SCCD_QID_EMPTY  --> fields fully published
        // Writers store qid LAST (after a fence); readers spin on qid
        // becoming non-empty after their CAS pop, then clear it.
        // This avoids the race where a popper observes a counter
        // increment ahead of the corresponding field writes.
        static constexpr int SCCD_QID_EMPTY = -1;

        // Pop one entry from the shared stack into out params (lane 0 only).
        template <typename T>
        static inline __device__ int try_pop_shared(int* SCCD_RESTRICT s_top,
                                                    T* SCCD_RESTRICT s_tlower,
                                                    T* SCCD_RESTRICT s_tupper,
                                                    T* SCCD_RESTRICT s_ulower,
                                                    T* SCCD_RESTRICT s_uupper,
                                                    T* SCCD_RESTRICT s_vlower,
                                                    T* SCCD_RESTRICT s_vupper,
                                                    int* SCCD_RESTRICT s_level,
                                                    int* SCCD_RESTRICT s_qid,
                                                    T& tl,
                                                    T& tu,
                                                    T& ul,
                                                    T& uu,
                                                    T& vl,
                                                    T& vu,
                                                    int& level,
                                                    int& qid) {
            const int slot = release_slot(s_top);
            if (slot < 0) return 0;
            // Spin until the writer publishes qid (they fenced first).
            int q;
            do {
                q = atomicAdd(&s_qid[slot], 0);
            } while (q == SCCD_QID_EMPTY);
            __threadfence_block();
            tl = s_tlower[slot];
            tu = s_tupper[slot];
            ul = s_ulower[slot];
            uu = s_uupper[slot];
            vl = s_vlower[slot];
            vu = s_vupper[slot];
            level = s_level[slot];
            qid = q;
            // Reset valid tag so the slot is reusable.
            atomicExch(&s_qid[slot], SCCD_QID_EMPTY);
            return 1;
        }

        // Pop one entry from the global stack into out params (lane 0 only).
        template <typename T>
        static inline __device__ int try_pop_global(int* SCCD_RESTRICT g_top,
                                                    T* SCCD_RESTRICT g_tlower,
                                                    T* SCCD_RESTRICT g_tupper,
                                                    T* SCCD_RESTRICT g_ulower,
                                                    T* SCCD_RESTRICT g_uupper,
                                                    T* SCCD_RESTRICT g_vlower,
                                                    T* SCCD_RESTRICT g_vupper,
                                                    int* SCCD_RESTRICT g_level,
                                                    int* SCCD_RESTRICT g_qid,
                                                    T& tl,
                                                    T& tu,
                                                    T& ul,
                                                    T& uu,
                                                    T& vl,
                                                    T& vu,
                                                    int& level,
                                                    int& qid) {
            const int slot = release_slot(g_top);
            if (slot < 0) return 0;
            int q;
            do {
                q = atomicAdd(&g_qid[slot], 0);
            } while (q == SCCD_QID_EMPTY);
            __threadfence();
            tl = g_tlower[slot];
            tu = g_tupper[slot];
            ul = g_ulower[slot];
            uu = g_uupper[slot];
            vl = g_vlower[slot];
            vu = g_vupper[slot];
            level = g_level[slot];
            qid = q;
            atomicExch(&g_qid[slot], SCCD_QID_EMPTY);
            return 1;
        }

        // Push up to two child entries (lane 0 only).  Any child whose
        // tlower is already >= the current toi[qid] is culled; this
        // rereads toi[qid] just before reserving slots so pushes see the
        // freshest tight bound established by peer warps.
        //
        // Memory discipline.  The global stack is split into two zones:
        //   [0, g_normal_cap)              : normal pushes
        //   [g_normal_cap, g_capacity)     : reserve for halt-flush
        //
        // If the shared stack is full AND the normal zone is full, the
        // push raises the halt flag and uses the reserve zone.  The
        // reserve is sized by the host (R = grid_blocks * (S_CAP + 2*W))
        // to guarantee that every block can later flush its shared
        // stack and every in-flight warp can finalize its push, so this
        // fallback path always succeeds.
        template <typename T>
        static inline __device__ void push_children(int* SCCD_RESTRICT s_top,
                                                    T* SCCD_RESTRICT s_tlower,
                                                    T* SCCD_RESTRICT s_tupper,
                                                    T* SCCD_RESTRICT s_ulower,
                                                    T* SCCD_RESTRICT s_uupper,
                                                    T* SCCD_RESTRICT s_vlower,
                                                    T* SCCD_RESTRICT s_vupper,
                                                    int* SCCD_RESTRICT s_level,
                                                    int* SCCD_RESTRICT s_qid,
                                                    int s_capacity,
                                                    int* SCCD_RESTRICT g_top,
                                                    T* SCCD_RESTRICT g_tlower,
                                                    T* SCCD_RESTRICT g_tupper,
                                                    T* SCCD_RESTRICT g_ulower,
                                                    T* SCCD_RESTRICT g_uupper,
                                                    T* SCCD_RESTRICT g_vlower,
                                                    T* SCCD_RESTRICT g_vupper,
                                                    int* SCCD_RESTRICT g_level,
                                                    int* SCCD_RESTRICT g_qid,
                                                    int g_capacity,
                                                    int g_normal_cap,
                                                    int* SCCD_RESTRICT halt,
                                                    const T* SCCD_RESTRICT toi_qid,
                                                    const T tl0,
                                                    const T tu0,
                                                    const T ul0,
                                                    const T uu0,
                                                    const T vl0,
                                                    const T vu0,
                                                    const T tl1,
                                                    const T tu1,
                                                    const T ul1,
                                                    const T uu1,
                                                    const T vl1,
                                                    const T vu1,
                                                    const int level,
                                                    const int qid) {
            // Relaxed load is fine: missing a freshly-tightened toi only
            // results in an extra push that a later pop will cull.
            const T cur = *toi_qid;
            const bool keep0 = tl0 < cur;
            const bool keep1 = tl1 < cur;
            const int n = (keep0 ? 1 : 0) + (keep1 ? 1 : 0);
            if (n == 0) return;

            int slot = reserve_slots(s_top, n, s_capacity);
            if (slot >= 0) {
                int k = 0;
                if (keep0) {
                    s_tlower[slot + k] = tl0;
                    s_tupper[slot + k] = tu0;
                    s_ulower[slot + k] = ul0;
                    s_uupper[slot + k] = uu0;
                    s_vlower[slot + k] = vl0;
                    s_vupper[slot + k] = vu0;
                    s_level[slot + k] = level;
                    ++k;
                }
                if (keep1) {
                    s_tlower[slot + k] = tl1;
                    s_tupper[slot + k] = tu1;
                    s_ulower[slot + k] = ul1;
                    s_uupper[slot + k] = uu1;
                    s_vlower[slot + k] = vl1;
                    s_vupper[slot + k] = vu1;
                    s_level[slot + k] = level;
                }

                __threadfence_block();
                // Publish qid LAST so a popper that observes a non-empty
                // tag is guaranteed to see committed field writes.
                atomicExch(&s_qid[slot], qid);
                if (n == 2) atomicExch(&s_qid[slot + 1], qid);
                // Fence again so the subsequent atomicSub(inflight, 1)
                // in the caller is ordered AFTER this push is visible.
                // Without this, a peer warp observing inflight==0 could
                // still see a stale s_top and terminate prematurely.
                __threadfence_block();
                return;
            }

            // Normal zone first.  If it's full, raise halt and spill
            // into the reserve zone -- which, by construction, always
            // has room for in-flight pushes plus the eventual shared-
            // stack flushes.
            slot = reserve_slots(g_top, n, g_normal_cap);
            if (slot < 0) {
                atomicExch(halt, 1);
                slot = reserve_slots(g_top, n, g_capacity);
                assert(slot >= 0);
            }
            int k = 0;
            if (keep0) {
                g_tlower[slot + k] = tl0;
                g_tupper[slot + k] = tu0;
                g_ulower[slot + k] = ul0;
                g_uupper[slot + k] = uu0;
                g_vlower[slot + k] = vl0;
                g_vupper[slot + k] = vu0;
                g_level[slot + k] = level;
                ++k;
            }
            if (keep1) {
                g_tlower[slot + k] = tl1;
                g_tupper[slot + k] = tu1;
                g_ulower[slot + k] = ul1;
                g_uupper[slot + k] = uu1;
                g_vlower[slot + k] = vl1;
                g_vupper[slot + k] = vu1;
                g_level[slot + k] = level;
            }

            __threadfence();
            atomicExch(&g_qid[slot], qid);
            if (n == 2) atomicExch(&g_qid[slot + 1], qid);
            // Device-scope fence so any peer block observing the
            // caller's upcoming atomicSub(inflight, 1) also sees g_top
            // and the published qids.
            __threadfence();
        }

        // ----------------------------------------------------------------
        // narrow_phase_ee_kernel (warp-per-query, persistent worker)
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
        template <typename T, typename I>
        __global__ void narrow_phase_ee_kernel(const I* const SCCD_RESTRICT e0overlap,
                                               const I* const SCCD_RESTRICT e1overlap,
                                               T** const SCCD_RESTRICT sp,
                                               T** const SCCD_RESTRICT ep,
                                               const size_t edge_stride,
                                               I** const SCCD_RESTRICT edges,
                                               const T tol,
                                               T* SCCD_RESTRICT toi,
                                               T* SCCD_RESTRICT g_tlower,
                                               T* SCCD_RESTRICT g_tupper,
                                               T* SCCD_RESTRICT g_ulower,
                                               T* SCCD_RESTRICT g_uupper,
                                               T* SCCD_RESTRICT g_vlower,
                                               T* SCCD_RESTRICT g_vupper,
                                               int* SCCD_RESTRICT g_level,
                                               int* SCCD_RESTRICT g_qid,
                                               int* SCCD_RESTRICT g_top,
                                               int g_capacity,
                                               int g_normal_cap,
                                               int* SCCD_RESTRICT inflight,
                                               int* SCCD_RESTRICT seed_cursor,
                                               const int seed_end,
                                               int* SCCD_RESTRICT halt) {
            using Vec4 = typename device::Vec4Type<T>::type;

            constexpr int S_CAP = SCCD_NP_SHARED_STACK_CAP;
            constexpr int W = SCCD_NP_WARPS_PER_BLOCK;

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
            // the rest of the warp (avoids needing a stream of shfls).
            __shared__ int w_status[W];  // 0=work, 1=backoff, 2=terminate
            __shared__ T w_tl[W], w_tu[W], w_ul[W], w_uu[W], w_vl[W], w_vu[W];
            __shared__ int w_level[W];
            __shared__ int w_qid[W];

            const int lane = threadIdx.x;
            const int warp = threadIdx.y;

            // Initialize shared stack header and tag all slots as empty
            // before any warp can touch them.
            const int tid_in_block = lane + warp * 32;
            const int threads_in_block = 32 * SCCD_NP_WARPS_PER_BLOCK;
            for (int i = tid_in_block; i < S_CAP; i += threads_in_block) {
                s_qid[i] = SCCD_QID_EMPTY;
            }
            if (tid_in_block == 0) {
                s_top = 0;
            }
            __syncthreads();

            unsigned int backoff = 64;

            while (true) {
                if (lane == 0) {
                    // Fast-path abort if another warp raised the halt
                    // flag (global stack hit its normal-zone cap).  We
                    // leave in-flight warps to finish their current
                    // iteration; at next top-of-loop they'll also exit
                    // here.  The block-flush below moves remaining
                    // shared-stack work to global so the host can
                    // relaunch and pick it up.
                    w_status[warp] = (atomicAdd(halt, 0) != 0) ? 2 : 0;
                }
                __syncwarp();
                if (w_status[warp] == 2) break;

                if (lane == 0) {
                    // Optimistically claim a worker slot before touching
                    // any stack counter so that termination cannot race
                    // ahead of an in-progress pop.
                    atomicAdd(inflight, 1);

                    T tl, tu, ul, uu, vl, vu;
                    int level, qid;
                    int got = try_pop_shared<T>(&s_top,
                                                s_tlower,
                                                s_tupper,
                                                s_ulower,
                                                s_uupper,
                                                s_vlower,
                                                s_vupper,
                                                s_level,
                                                s_qid,
                                                tl,
                                                tu,
                                                ul,
                                                uu,
                                                vl,
                                                vu,
                                                level,
                                                qid);
                    if (!got) {
                        got = try_pop_global<T>(g_top,
                                                g_tlower,
                                                g_tupper,
                                                g_ulower,
                                                g_uupper,
                                                g_vlower,
                                                g_vupper,
                                                g_level,
                                                g_qid,
                                                tl,
                                                tu,
                                                ul,
                                                uu,
                                                vl,
                                                vu,
                                                level,
                                                qid);
                    }

                    // Source 3: lazy seed claim.  Synthesizes the root
                    // entry for query `old` in registers -- no global
                    // push/pop roundtrip.  A seed_cursor that overshoots
                    // seed_end is harmless: peers observing the larger
                    // counter simply read "no seeds left".
                    if (!got) {
                        const int old = atomicAdd(seed_cursor, 1);
                        if (old < seed_end) {
                            tl = T(0);
                            tu = T(1);
                            ul = T(0);
                            uu = T(1);
                            vl = T(0);
                            vu = T(1);
                            level = 0;
                            qid = old;
                            got = 1;
                        }
                    }

                    if (got) {
                        w_tl[warp] = tl;
                        w_tu[warp] = tu;
                        w_ul[warp] = ul;
                        w_uu[warp] = uu;
                        w_vl[warp] = vl;
                        w_vu[warp] = vu;
                        w_level[warp] = level;
                        w_qid[warp] = qid;
                        w_status[warp] = 0;
                    } else {
                        // Release the optimistic claim and assess
                        // termination.  Ordering matters: only after
                        // releasing can other warps observe a true
                        // 'no inflight' state.
                        atomicSub(inflight, 1);
                        const int gtop_now = atomicAdd(g_top, 0);
                        const int stop_now = atomicAdd(&s_top, 0);
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
                if (status == 2) {
                    break;
                }
                if (status == 1) {
                    backoff_ns(backoff);
                    if (backoff < 4096) backoff <<= 1;
                    continue;
                }
                backoff = 64;

                const T tl = w_tl[warp];
                const T tu = w_tu[warp];
                const T ul = w_ul[warp];
                const T uu = w_uu[warp];
                const T vl = w_vl[warp];
                const T vu = w_vu[warp];
                const int level = w_level[warp];
                const int qid = w_qid[warp];

                // Early cull if a tighter toi has already been recorded
                // for this query.  Read once on lane 0 and broadcast so
                // all 32 lanes agree on the value; otherwise a concurrent
                // atomic_min by another warp could make the plain load
                // diverge across lanes and break the __syncwarp() below.
                T cur_toi = T(0);
                if (lane == 0) cur_toi = toi[qid];
                cur_toi = __shfl_sync(0xffffffffu, cur_toi, 0);
                if (tl >= cur_toi) {
                    if (lane == 0) atomicSub(inflight, 1);
                    __syncwarp();
                    continue;
                }

                Vec4 sx, sy, sz, ex, ey, ez;
                load_query_ee<T, Vec4, I>(
                    qid, e0overlap, e1overlap, sp, ep, edge_stride, edges, sx, sy, sz, ex, ey, ez);

                T atol[3];
                compute_edge_edge_tolerance_soa<T, Vec4>(tol, sx, sy, sz, ex, ey, ez, &atol[0], &atol[1], &atol[2]);

                // 2 x 4 x 4 lane decomposition.
                const int ti = lane >> 4;
                const int ui = (lane >> 2) & 0x3;
                const int vi = lane & 0x3;

                int contains = 0;
                int accept = 0;
                T cell_tl = tl;
                sample_cell_3d<T, Vec4>(ti,
                                        ui,
                                        vi,
                                        2,
                                        4,
                                        4,
                                        tl,
                                        tu,
                                        ul,
                                        uu,
                                        vl,
                                        vu,
                                        sx,
                                        sy,
                                        sz,
                                        ex,
                                        ey,
                                        ez,
                                        tol,
                                        atol,
                                        contains,
                                        accept,
                                        cell_tl);

                const unsigned acc_mask = __ballot_sync(0xffffffffu, accept);
                const unsigned cont_mask = __ballot_sync(0xffffffffu, contains);

                if (acc_mask) {
                    // Reduce min(cell_tl) over accepting lanes.  Use the
                    // parent subdomain upper-bound as a sentinel so that
                    // non-accepting lanes never win the reduction.
                    T best = accept ? cell_tl : tu;
                    for (int o = 16; o > 0; o >>= 1) {
                        T other = __shfl_xor_sync(0xffffffffu, best, o);
                        best = (best < other) ? best : other;
                    }
                    if (lane == 0) {
                        device::atomic_min(&toi[qid], best);
                        atomicSub(inflight, 1);
                    }
                    __syncwarp();
                    continue;
                }

                if (cont_mask) {
                    // Bisect the widest of (t, u, v) at the parent level.
                    const T tw = tu - tl;
                    const T uw = uu - ul;
                    const T vw = vu - vl;
                    int axis = 0;
                    T width = tw;
                    if (uw > width) {
                        axis = 1;
                        width = uw;
                    }
                    if (vw > width) {
                        axis = 2;
                    }

                    T tl0 = tl, tu0 = tu, ul0 = ul, uu0 = uu, vl0 = vl, vu0 = vu;
                    T tl1 = tl, tu1 = tu, ul1 = ul, uu1 = uu, vl1 = vl, vu1 = vu;
                    if (axis == 0) {
                        const T m = (tl + tu) * T(0.5);
                        tu0 = m;
                        tl1 = m;
                    } else if (axis == 1) {
                        const T m = (ul + uu) * T(0.5);
                        uu0 = m;
                        ul1 = m;
                    } else {
                        const T m = (vl + vu) * T(0.5);
                        vu0 = m;
                        vl1 = m;
                    }

                    if (lane == 0) {
                        push_children<T>(&s_top,
                                         s_tlower,
                                         s_tupper,
                                         s_ulower,
                                         s_uupper,
                                         s_vlower,
                                         s_vupper,
                                         s_level,
                                         s_qid,
                                         S_CAP,
                                         g_top,
                                         g_tlower,
                                         g_tupper,
                                         g_ulower,
                                         g_uupper,
                                         g_vlower,
                                         g_vupper,
                                         g_level,
                                         g_qid,
                                         g_capacity,
                                         g_normal_cap,
                                         halt,
                                         &toi[qid],
                                         tl0,
                                         tu0,
                                         ul0,
                                         uu0,
                                         vl0,
                                         vu0,
                                         tl1,
                                         tu1,
                                         ul1,
                                         uu1,
                                         vl1,
                                         vu1,
                                         level + 1,
                                         qid);
                        atomicSub(inflight, 1);
                    }
                    __syncwarp();
                    continue;
                }

                // Neither accept nor contains-origin: discard.
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
                flush_base = (flush_n > 0) ? reserve_slots(g_top, flush_n, g_capacity) : 0;
                // By construction (host sizes g_capacity with reserve
                // >= grid_blocks * S_CAP) this reserve must succeed.
                if (flush_n > 0) assert(flush_base >= 0);
            }
            __syncthreads();

            for (int i = tid_in_block; i < flush_n; i += threads_in_block) {
                const int src = i;
                const int dst = flush_base + i;
                g_tlower[dst] = s_tlower[src];
                g_tupper[dst] = s_tupper[src];
                g_ulower[dst] = s_ulower[src];
                g_uupper[dst] = s_uupper[src];
                g_vlower[dst] = s_vlower[src];
                g_vupper[dst] = s_vupper[src];
                g_level[dst] = s_level[src];
                // s_qid[src] must be a published (non-EMPTY) tag because
                // the invariant is that slots in [0, s_top) are fully
                // committed.  Republish on the global side last, after a
                // fence, so future popper sees fields first.
                const int q = s_qid[src];
                __threadfence();
                atomicExch(&g_qid[dst], q);
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

            if (noverlaps == 0) return max_toi;

            T tol = T(1e-8);
            {
                double SCCD_TOL = (double)tol;
                SCCD_READ_ENV(SCCD_TOL, atof);
                tol = (T)SCCD_TOL;
            }

            // Candidates-per-kernel cap.  Default 0 == unbounded (all
            // overlaps in one worker launch).  When positive, the host
            // splits the input into contiguous [begin, end) batches and
            // relaunches the worker for each one, reusing all device
            // buffers across batches.
            int SCCD_BATCH_SIZE = 0;
            SCCD_READ_ENV(SCCD_BATCH_SIZE, atoi);
            const size_t batch_size = (SCCD_BATCH_SIZE > 0) ? (size_t)SCCD_BATCH_SIZE : noverlaps;

            // Grid sizing up front: we need it to bound the reserve R.
            int dev = 0;
            cudaGetDevice(&dev);
            cudaDeviceProp prop{};
            cudaGetDeviceProperties(&prop, dev);

            int SCCD_BLOCKS_PER_SM = 4;
            SCCD_READ_ENV(SCCD_BLOCKS_PER_SM, atoi);

            int base_grid_blocks = prop.multiProcessorCount * SCCD_BLOCKS_PER_SM;
            if (base_grid_blocks <= 0) base_grid_blocks = 1;

            // Reserve count R sized so that at halt time every block can
            // flush its entire shared stack (S_CAP) and every warp can
            // complete an in-flight push (2 children) without running
            // out of global-stack slots.
            const int W = SCCD_NP_WARPS_PER_BLOCK;
            const int reserve_R = base_grid_blocks * (SCCD_NP_SHARED_STACK_CAP + 2 * W);

            // Global stack cap.  Caller's env value is a starting point;
            // we enforce gstack_cap >= 2*R so there's a non-trivial
            // normal zone on top of the reserve.  With tiny caps the
            // host simply relaunches more times.
            int SCCD_GSTACK_CAP = (int)(noverlaps * 16) + 1024;
            SCCD_READ_ENV(SCCD_GSTACK_CAP, atoi);
            int gstack_cap = SCCD_GSTACK_CAP;
            if (gstack_cap < 2 * reserve_R) gstack_cap = 2 * reserve_R;
            const int g_normal_cap = gstack_cap - reserve_R;

            // Per-query toi output (allocated once across all batches).
            T* d_toi = nullptr;
            cudaMalloc(&d_toi, noverlaps * sizeof(T));

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

            // One-time init: toi[] = max_toi, g_qid[] = SCCD_QID_EMPTY.
            // Workers now pull seeds lazily, so the global stack is left
            // empty here (no per-query push-then-pop bookkeeping).
            {
                const int block = 256;
                const size_t work = (gstack_cap > (int)noverlaps) ? (size_t)gstack_cap : noverlaps;
                const int grid = (int)((work + block - 1) / block);
                init_narrow_phase_ee_kernel<T><<<grid, block>>>(noverlaps,
                                                                gstack_cap,
                                                                max_toi,
                                                                d_toi,
                                                                g_qid);
                SCCD_CUDA_LAST_ERROR();
            }

            dim3 block(32, SCCD_NP_WARPS_PER_BLOCK, 1);

            // Batch loop.  Each iteration claims [begin, end) of the
            // input.  Inside, we may relaunch the worker multiple times
            // until that range's work fully drains; this is how we cap
            // peak device memory -- when the global stack's normal zone
            // fills, workers halt + flush + exit and the host relaunches
            // with `g_top` and `seed_cursor` preserved.
            for (size_t begin = 0; begin < noverlaps; begin += batch_size) {
                const size_t end = (begin + batch_size < noverlaps) ? (begin + batch_size) : noverlaps;
                const int this_batch = (int)(end - begin);

                reset_batch_narrow_phase_ee_kernel<<<1, 1>>>(g_top, inflight, seed_cursor, halt, (int)begin);
                SCCD_CUDA_LAST_ERROR();

                while (true) {
                    int grid_blocks = base_grid_blocks;
                    if (grid_blocks > this_batch) grid_blocks = this_batch;
                    dim3 grid(grid_blocks, 1, 1);

                    narrow_phase_ee_kernel<T, I><<<grid, block>>>(e0overalp,
                                                                  e1overalp,
                                                                  v0,
                                                                  v1,
                                                                  edge_stride,
                                                                  edges,
                                                                  tol,
                                                                  d_toi,
                                                                  g_tlower,
                                                                  g_tupper,
                                                                  g_ulower,
                                                                  g_uupper,
                                                                  g_vlower,
                                                                  g_vupper,
                                                                  g_level,
                                                                  g_qid,
                                                                  g_top,
                                                                  gstack_cap,
                                                                  g_normal_cap,
                                                                  inflight,
                                                                  seed_cursor,
                                                                  (int)end,
                                                                  halt);
                    SCCD_CUDA_LAST_ERROR();

                    // Check whether the batch is fully drained.  A D2H
                    // copy per relaunch is cheap compared to the kernel
                    // runtime and is only paid once per overflow event.
                    int h_g_top = 0;
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

            // Reduce per-query toi to a single value on device, then copy.
            T h_toi = thrust::reduce(thrust::device_pointer_cast(d_toi),
                                     thrust::device_pointer_cast(d_toi) + noverlaps,
                                     max_toi,
                                     thrust::minimum<T>());

            cudaFree(d_toi);
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