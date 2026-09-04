#include "sccd_narrowphase.cuh"

#include <stdint.h>
#include <stdlib.h>

#include <cassert>
#include <cstdio>
#include <climits>
#include <limits>
#include <vector>

#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <thrust/functional.h>
#include <thrust/memory.h>
#include <thrust/reduce.h>

#include "sccd_cuda_base.cuh"
#include "sccd_narrowphase_mode.hpp"
#include "sccd_reduce.cuh"

// Entries in the per-block shared stack. An entry is six scalars plus two ints,
// so 1024 entries is 32,800 bytes of static __shared__ for float and 57,376 for
// double. See SharedStackCap below for why that is legal and what it costs.
#ifndef SCCD_NP_SHARED_STACK_CAP
#define SCCD_NP_SHARED_STACK_CAP 64
#endif

#ifndef SCCD_NP_THREADS_PER_BLOCK
#define SCCD_NP_THREADS_PER_BLOCK 128
#endif

#ifndef SCCD_CUDA_ADAPTIVE_SPLIT
#define SCCD_CUDA_ADAPTIVE_SPLIT 1
#endif

// How many boxes one block takes off the global stack per launch. Bounded only
// so a pathological run cannot spin a block forever; the host loop relaunches
// while anything remains, so lowering this costs rounds but never correctness.
#ifndef SCCD_NP_DRAIN_PER_BLOCK
#define SCCD_NP_DRAIN_PER_BLOCK 4096
#endif

// Levels of breadth-first bisection used to seed a block on the conservative
// path. Enough to fill 128 threads when the tree branches (2^7), with headroom;
// a tree too narrow to fill the block in this many levels is handed to the DFS
// rather than pursued single-file here.
// Levels of breadth-first bisection used to seed a block on the conservative
// path. 7 levels is enough to fill 128 threads when the tree branches; the extra
// is headroom for levels that prune. A tree too narrow to fill the block within
// this many levels is handed to the DFS rather than pursued single-file here.
#ifndef SCCD_NP_RAMP_LEVELS
#define SCCD_NP_RAMP_LEVELS 12
#endif

// How many times a batch may be regrown and rerun before the driver gives up.
// See the comment at the retry site: giving up keeps the result conservative,
// where spinning does not keep it finite.
#ifndef SCCD_NP_MAX_RETRY_ROUNDS
#define SCCD_NP_MAX_RETRY_ROUNDS 32
#endif

// Without a bound, nvcc budgets registers for the largest block it must support
// rather than the 128 threads these kernels are always launched with, and it is
// free to let the register count drift between toolkit versions. Stating the
// block size costs nothing and pins that down.
//
// The second argument -- a minimum number of resident blocks per SM -- is the
// occupancy/spill tuning knob and is deliberately NOT set by default: these
// kernels inline adaptive_split_longest_axis, which is large, and forcing the
// budget down would trade a visible occupancy number for invisible local-memory
// spills. Set SCCD_NP_MIN_BLOCKS_PER_SM only against a measured -Xptxas -v run
// on the target architecture (see SCCD_CUDA_VERBOSE_PTXAS).
#ifdef SCCD_NP_MIN_BLOCKS_PER_SM
#define SCCD_NP_LAUNCH_BOUNDS(N) __launch_bounds__(N, SCCD_NP_MIN_BLOCKS_PER_SM)
#else
#define SCCD_NP_LAUNCH_BOUNDS(N) __launch_bounds__(N)
#endif

namespace sccd {
    namespace device {

        /**
         * \brief atomicMin on a T-typed time of impact, from a double the search
         *        computed, rounding the narrowing DOWN.
         *
         * Round-to-nearest is wrong here. Narrowing a double time of impact to
         * float can round *up*, which would publish a value later than the one
         * the search actually proved -- and a time of impact later than the truth
         * is the failure this whole kernel exists to prevent. Rounding toward
         * negative infinity can only report earlier, which is always safe.
         *
         * Returns the running minimum, widened back to double.
         */
        template <typename T>
        static inline __device__ double atomic_min_toi(T* const address, const double value);

        template <>
        inline __device__ double atomic_min_toi<double>(double* const address, const double value) {
            return device::atomic_min(address, value);
        }

        template <>
        inline __device__ double atomic_min_toi<float>(float* const address, const double value) {
            return (double)device::atomic_min(address, __double2float_rd(value));
        }

        /**
         * \brief The scalar type the narrow-phase search computes in: always double.
         *
         * The conservativeness argument rests on TightInclusion's numerical error
         * bound and on tolerances derived from it, and single precision is not
         * accurate enough to carry them: measured on GH200, the float search
         * reported nine times of impact later than the true root on
         * armadillo-rollers, which is an invariant violation.
         *
         * The template parameter T stays the *storage* type. Geometry is read as
         * T and widened on load, which is exact, so a float caller keeps float
         * buffers and gets a float time of impact -- only the search itself is
         * double. The host kernel does the same thing:
         * sccd_vnarrowphase_ti.hpp opens with `using T_HP = double`.
         */
        using TC = double;

        // Corner-evaluation counter, for answering "is the gap work or cost?".
        //
        // Compiled out unless SCCD_NP_COUNT_BOXES is defined, so the shipped
        // kernel is untouched -- it is a diagnostic, and one that changes what it
        // measures: a global atomicAdd on every evaluation serialises the kernel,
        // so an instrumented build's *timings* mean nothing. Only the counts do.
#ifdef SCCD_NP_COUNT_BOXES
        __device__ unsigned long long g_np_evals = 0;
        // Per-query counts, so the mean can be taken apart. Null unless the
        // driver allocated it for this call.
        __device__ unsigned long long* g_np_perq = nullptr;
        // Boxes by DFS level, and how many were accepted only because the depth
        // cap forced it. Together these say whether the search converges: a
        // converging search empties out well before the cap and almost never
        // accepts at it.
        __device__ unsigned long long g_np_level[80] = {0};
        __device__ unsigned long long g_np_depth_accept = 0;
        // Where pushed boxes land: the block's shared stack, the shared global
        // queue when that is full, or nowhere when the queue is full too. Says
        // whether the flush-and-drain machinery is being exercised at all.
        __device__ unsigned long long g_np_push_shared = 0;
        __device__ unsigned long long g_np_push_global = 0;
        __device__ unsigned long long g_np_push_lost = 0;
#define SCCD_NP_EVAL_TICK() atomicAdd(&g_np_evals, 1ull)
#define SCCD_NP_PERQ_TICK(qid, n)                                             \
    do {                                                                      \
        if (g_np_perq != nullptr && (qid) >= 0) {                             \
            atomicAdd(&g_np_perq[(qid)], (unsigned long long)(n));            \
        }                                                                     \
    } while (0)
#else
#define SCCD_NP_EVAL_TICK() ((void)0)
#define SCCD_NP_PERQ_TICK(qid, n) ((void)0)
#endif

        // Shared-stack capacity for scalar type T.
        //
        // Small on purpose. The block-local stack is where a query's subtree
        // lives, so a large one keeps a heavy query inside the block that seeded
        // it: measured per query on cloth-funnel, the conservative search needs
        // one box for most queries and over a million for 208 of them, and at
        // 1024 entries per block 99% of pushes never reach the global queue at
        // all. One block then grinds through a 19.7-million-box query while its
        // neighbours finish and idle.
        //
        // Spilling early hands those boxes to the global queue, which the next
        // drain round redistributes across every block. Measured on GH200, mode 2,
        // 16 cases (narrow-phase milliseconds):
        //
        //                        cap 1024   cap 256   cap 64
        //   cloth-funnel           4059.2     198.2     37.2
        //   armadillo-rollers      3133.7     277.2     75.8
        //   cloth-ball              129.8     110.2    105.9
        //
        // False positives and negatives are identical at every capacity, and mode
        // 0 -- which averages under four boxes per query and so never fills even a
        // small stack -- is unaffected.
        //
        // This only became available once the global queue stopped deadlocking
        // under sustained use; see the Stack comment above. Before that, shrinking
        // this number hung the kernel, which is why it was 1024.
        //
        // 64 entries is 3,584 bytes for double against 57,376 at 1024. That does
        // not buy occupancy -- at 238 registers per thread the kernel is capped at
        // two blocks per SM either way -- it buys balance.
        template <typename T>
        struct SharedStackCap {
            static constexpr int value = SCCD_NP_SHARED_STACK_CAP;
        };

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

        /**
         * \brief One half of the double-buffered global work queue.
         *
         * A launch reads from one buffer and writes to the other, never both, so
         * there is no producer/consumer handshake here and no spinning. Reading
         * uses `top` as a claim cursor over the `count` entries a previous launch
         * left; writing uses `top` as a bump allocator and `count` is unused.
         *
         * This replaces a single shared buffer in which producers and consumers
         * ran concurrently and coordinated through the `qid` array: a writer
         * spun on `atomicCAS(&qid[slot], EMPTY, WRITING)` until a reader released
         * the slot, and a reader spun until the writer committed. With the queue
         * lightly used that handshake was invisible; driven hard -- which is
         * exactly what using it to balance the load requires -- it hung the
         * kernel. Separating the buffers removes the interaction rather than
         * tuning around it.
         */
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
            int* top;      // write cursor when writing, claim cursor when reading
            int* request;  // boxes that did not fit, for the host to size the retry
            int capacity;
            int count;     // entries available to readers; 0 for a write buffer
        };

        /**
         * \brief Claim one slot and write a box to it. No waiting, ever.
         *
         * The index comes from a bump allocator, so it belongs to this thread
         * alone and nobody reads the buffer until the launch ends. A box that does
         * not fit is counted in `request` and dropped, which is safe -- dropping
         * can only leave the time of impact earlier than the truth, never later --
         * and the host grows the queue and retries.
         */
        template <typename T>
        static inline __device__ void push_global(const Stack<T>& out,
                                                  const Domain<T>& box,
                                                  const int level,
                                                  const int qid) {
            const int slot = atomicAdd(out.top, 1);
            if (slot >= out.capacity) {
                atomicAdd(out.request, 1);
                return;
            }
            out.tlower[slot] = box.tlower;
            out.tupper[slot] = box.tupper;
            out.ulower[slot] = box.ulower;
            out.uupper[slot] = box.uupper;
            out.vlower[slot] = box.vlower;
            out.vupper[slot] = box.vupper;
            out.level[slot] = level;
            out.qid[slot] = qid;
        }

        /**
         * \brief Take the next box from the read buffer, or return 0 when empty.
         *
         * A plain claim on a cursor over entries a previous launch committed, so
         * the loads need no ordering against a concurrent writer -- there is none.
         */
        template <typename T>
        static inline __device__ int pop_global(const Stack<T>& in,
                                                Domain<T>& d,
                                                int& level,
                                                int& qid) {
            if (in.count <= 0) return 0;
            const int slot = atomicAdd(in.top, 1);
            if (slot >= in.count) return 0;

            d.tlower = in.tlower[slot];
            d.tupper = in.tupper[slot];
            d.ulower = in.ulower[slot];
            d.uupper = in.uupper[slot];
            d.vlower = in.vlower[slot];
            d.vupper = in.vupper[slot];
            level = in.level[slot];
            qid = in.qid[slot];
            return 1;
        }

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
                fmin = SCCD_MIN(fmin, f[i]);
                fmax = SCCD_MAX(fmax, f[i]);
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

#define SCCD_MASK_FULL 0xf
#define SCCD_MASK_DOMAIN_SMALLER_THAN_TOL 1
#define SCCD_MASK_BOX_INSIDE_EPSILON_BOX 2
#define SCCD_MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE 4
#define SCCD_MASK_INTERVAL_TERMINAL 8

        template <typename T>
        __device__ uint8_t cond_mask(const T fmin, const T fmax, const T tol, const T adaptive_tol) {
            bool cond1 = (fmax - fmin <= adaptive_tol);
            bool cond2 = (fmin > -tol) & (fmax < tol);
            bool cond3 = (fmax - fmin < tol);
            bool cond4 = (fmin >= fmax);

            uint8_t cond_mask = (cond1 ? SCCD_MASK_DOMAIN_SMALLER_THAN_TOL : 0);
            cond_mask |= (cond2 ? SCCD_MASK_BOX_INSIDE_EPSILON_BOX : 0);
            cond_mask |= (cond3 ? SCCD_MASK_REAL_TOLERANCE_SMALLER_THAN_INT_TOLERANCE : 0);
            cond_mask |= (cond4 ? SCCD_MASK_INTERVAL_TERMINAL : 0);
            cond_mask &= ((fmin <= tol) && (fmax >= -tol)) ? SCCD_MASK_FULL : 0;
            return cond_mask;
        }

        /**
         * \brief Mode-0 acceptance test.
         *
         * The origin-containment test -- the only test here whose failure
         * *rejects* a box -- is padded by `max(tol, aerr[d])`. It used to be
         * padded by `tol` alone, and that is unsound whenever the caller asks for
         * a distance tolerance below the certified numerical error bound: the pad
         * is then narrower than the error in the corner values, and a box holding
         * a root can be discarded. In double, with the bound at
         * `(vf ? 30 : 28) * eps * min(max_coord, 1)^3 <= 6.7e-15`, a typical
         * tolerance of 3e-8 is far wider and nothing changes -- which is why the
         * scenes measured never saw it. Taking the max makes that structural
         * instead of a property of the tolerance the caller happened to pass.
         *
         * The acceptance conditions keep using `tol`. Accepting a box is always
         * safe however loose the test, so those are an accuracy question, not a
         * soundness one.
         */
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
                                                       const T* const SCCD_RESTRICT aerr,
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
            const T x_pad = device::max<T>(tol, aerr[0]);
            int co = (fmin <= x_pad) & (fmax >= -x_pad);

            sample_f<is_vf, T, Vec4>(tl, tu, ul, uu, vl, vu, sy, ey, f);

            fminmax<T>(f, fmin, fmax);
            const T y_width = fmax - fmin;
            const uint8_t y_mask = cond_mask<T>(fmin, fmax, tol, adaptive_tol[1]);
            const T y_pad = device::max<T>(tol, aerr[1]);
            co &= (fmin <= y_pad) & (fmax >= -y_pad);

            sample_f<is_vf, T, Vec4>(tl, tu, ul, uu, vl, vu, sz, ez, f);

            fminmax<T>(f, fmin, fmax);
            const T z_width = fmax - fmin;
            const uint8_t z_mask = cond_mask<T>(fmin, fmax, tol, adaptive_tol[2]);
            const T z_pad = device::max<T>(tol, aerr[2]);
            co &= (fmin <= z_pad) & (fmax >= -z_pad);

            const uint8_t and_mask = (x_mask & y_mask & z_mask);
            const T true_tol = device::max<T>(x_width, device::max<T>(y_width, z_width));

            const bool cond1 = and_mask & SCCD_MASK_DOMAIN_SMALLER_THAN_TOL;
            const bool cond2 = and_mask & SCCD_MASK_BOX_INSIDE_EPSILON_BOX;
            const bool cond3 = (tl > T(0)) && (true_tol < tol);
            const bool cond4 = and_mask & SCCD_MASK_INTERVAL_TERMINAL;

            contains_origin = co;
            accept = co && (cond1 || cond2 || cond3 || cond4);
        }

        // ---------------------------------------------------------------------
        // TightInclusion-equivalent predicate and split.
        //
        // The kernels above implement the host's mode-0 acceptance test. It is
        // looser than TightInclusion's -- it compares a *codomain* width against a
        // *domain* tolerance -- but looseness is on the accepting side and costs
        // accuracy, not safety. Its rejection is now padded by the same certified
        // error bound as the kernels below (see evaluate_cell_3d); it used to be
        // padded by the caller's distance tolerance alone, which is the unsound
        // rejection that cost late times of impact in single precision on GH200
        // (benchmark/oracle/README.md).
        //
        // What follows is the device twin of src/sccd_vnarrowphase_ti.hpp, which
        // reproduces TightInclusion exactly and is the host's mode 2.
        // ---------------------------------------------------------------------

        template <typename T>
        struct FpTraits;

        template <>
        struct FpTraits<float> {
            static inline __device__ float epsilon() { return 1.1920928955078125e-7f; }
        };

        template <>
        struct FpTraits<double> {
            static inline __device__ double epsilon() { return 2.220446049250313e-16; }
        };

        /**
         * \brief Per-axis numerical error bound, the device twin of
         *        snumerical_error.hpp's numerical_error_bound_component.
         *
         * The host spells the cube as std::pow(delta, 3) to stay bit-identical
         * with TightInclusion. There is no reason to pay a libm call here: a
         * plain cube is within two roundings of it, and scaling by (1 + 4 eps)
         * puts the result at or above the correctly rounded value. Erring large
         * is the conservative direction -- a wider pad keeps boxes that would
         * otherwise be rejected, and a rejected box that held a root is a missed
         * collision.
         */
        template <bool is_vf, typename T, typename Vec4>
        static inline __device__ T numerical_error_bound_component(const Vec4 s, const Vec4 e) {
            const T max_abs =
                device::max<T>(device::max<T>(device::max<T>(device::abs<T>(s.x), device::abs<T>(s.y)),
                                              device::max<T>(device::abs<T>(s.z), device::abs<T>(s.w))),
                               device::max<T>(device::max<T>(device::abs<T>(e.x), device::abs<T>(e.y)),
                                              device::max<T>(device::abs<T>(e.z), device::abs<T>(e.w))));
            const T delta = device::min<T>(max_abs, T(1));
            const T eps = FpTraits<T>::epsilon();
            const T filter = T(is_vf ? 30 : 28) * eps;
            const T cube = delta * delta * delta;
            return filter * cube * (T(1) + T(4) * eps);
        }

        /**
         * \brief TightInclusion's acceptance test on one box.
         *
         * Three conditions and nothing else, transcribed from tight_classify in
         * src/sccd_vnarrowphase_ti.hpp:
         *   reject  if the origin is outside the padded range on any axis;
         *   accept  if the whole range is inside the error box on every axis;
         *   accept  if the DOMAIN width is within the domain tolerance on every
         *           axis -- a domain width against a domain tolerance, which is
         *           the comparison the mode-0 test gets wrong.
         */
        template <bool is_vf, typename T, typename Vec4>
        static inline __device__ void evaluate_cell_3d_tight(const Domain<T>& cell,
                                                          const Vec4 sx,
                                                          const Vec4 sy,
                                                          const Vec4 sz,
                                                          const Vec4 ex,
                                                          const Vec4 ey,
                                                          const Vec4 ez,
                                                          const T* const SCCD_RESTRICT atol,
                                                          const T* const SCCD_RESTRICT aerr,
                                                          int& contains_origin,
                                                          int& accept) {
            const T tl = cell.tlower, tu = cell.tupper;
            const T ul = cell.ulower, uu = cell.uupper;
            const T vl = cell.vlower, vu = cell.vupper;

            T f[8];
            T fmin, fmax;

            int co = 1;
            int box_in = 1;

            sample_f<is_vf, T, Vec4>(tl, tu, ul, uu, vl, vu, sx, ex, f);
            fminmax<T>(f, fmin, fmax);
            co &= (fmin <= aerr[0]) & (fmax >= -aerr[0]);
            box_in &= (fmin >= -aerr[0]) & (fmax <= aerr[0]);

            sample_f<is_vf, T, Vec4>(tl, tu, ul, uu, vl, vu, sy, ey, f);
            fminmax<T>(f, fmin, fmax);
            co &= (fmin <= aerr[1]) & (fmax >= -aerr[1]);
            box_in &= (fmin >= -aerr[1]) & (fmax <= aerr[1]);

            sample_f<is_vf, T, Vec4>(tl, tu, ul, uu, vl, vu, sz, ez, f);
            fminmax<T>(f, fmin, fmax);
            co &= (fmin <= aerr[2]) & (fmax >= -aerr[2]);
            box_in &= (fmin >= -aerr[2]) & (fmax <= aerr[2]);

            const int within_tol = ((tu - tl) <= atol[0]) & ((uu - ul) <= atol[1]) & ((vu - vl) <= atol[2]);

            contains_origin = co;
            accept = co & (box_in | within_tol);
        }

        /**
         * \brief Whether TightInclusion's split can still make progress on this box.
         *
         * False once the chosen axis' midpoint no longer separates its endpoints,
         * i.e. the interval has reached the resolution of the format. Checked
         * before splitting rather than reported from inside it, so the caller
         * handles it next to the depth cutoff, which is the same situation: the
         * box must be accepted at its t lower bound, never dropped.
         */
        template <typename T>
        static inline __device__ bool tight_can_split(const Domain<T>& in, const T* const SCCD_RESTRICT atol) {
            const T width[3] = {in.tupper - in.tlower, in.uupper - in.ulower, in.vupper - in.vlower};
            int axis = 0;
            T best = -T(1);
            for (int d = 0; d < 3; ++d) {
                if (width[d] > atol[d]) {
                    const T ratio = width[d] / atol[d];
                    if (ratio > best) {
                        best = ratio;
                        axis = d;
                    }
                }
            }
            const T lo = (axis == 0) ? in.tlower : ((axis == 1) ? in.ulower : in.vlower);
            const T hi = (axis == 0) ? in.tupper : ((axis == 1) ? in.uupper : in.vupper);
            const T mid = (lo + hi) * T(0.5);
            return lo < mid && mid < hi;
        }

        /**
         * \brief TightInclusion's split: the axis furthest past its tolerance,
         *        halved at the midpoint.
         *
         * Returns false when the midpoint fails to separate the endpoints, i.e.
         * the interval has reached the resolution of the format. The caller must
         * then accept the box rather than split it: accepting reports the box's t
         * lower bound, which is early and therefore safe, while dropping it could
         * lose a root.
         */
        template <typename T>
        static inline __device__ bool bisect_ti_axis(const Domain<T>& in,
                                                     const T* const SCCD_RESTRICT atol,
                                                     Domain<T>& left,
                                                     Domain<T>& right) {
            const T width[3] = {in.tupper - in.tlower, in.uupper - in.ulower, in.vupper - in.vlower};

            int axis = 0;
            T best = -T(1);
            for (int d = 0; d < 3; ++d) {
                if (width[d] > atol[d]) {
                    const T ratio = width[d] / atol[d];
                    if (ratio > best) {
                        best = ratio;
                        axis = d;
                    }
                }
            }

            const T lo = (axis == 0) ? in.tlower : ((axis == 1) ? in.ulower : in.vlower);
            const T hi = (axis == 0) ? in.tupper : ((axis == 1) ? in.uupper : in.vupper);
            const T mid = (lo + hi) * T(0.5);
            if (!(lo < mid && mid < hi)) {
                return false;
            }

            left = in;
            right = in;
            if (axis == 0) {
                left.tupper = mid;
                right.tlower = mid;
            } else if (axis == 1) {
                left.uupper = mid;
                right.ulower = mid;
            } else {
                left.vupper = mid;
                right.vlower = mid;
            }
            return true;
        }

        template <bool is_vf, bool conservative, typename T, typename Vec4>
        static inline __device__ void evaluate_cell_3d_policy(const Domain<T>& cell,
                                                              const Vec4 sx,
                                                              const Vec4 sy,
                                                              const Vec4 sz,
                                                              const Vec4 ex,
                                                              const Vec4 ey,
                                                              const Vec4 ez,
                                                              const T tol,
                                                              const T* const SCCD_RESTRICT atol,
                                                              const T* const SCCD_RESTRICT aerr,
                                                              int& contains_origin,
                                                              int& accept);

        template <bool is_vf, bool conservative, typename T, typename Vec4>
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
                                                     const T* const SCCD_RESTRICT aerr,
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

            evaluate_cell_3d_policy<is_vf, conservative, TC, Vec4>(
                cell, sx, sy, sz, ex, ey, ez, tol, adaptive_tol, aerr, contains_origin, accept);
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

        /**
         * \brief Load a query's geometry and its per-axis domain tolerances.
         *
         * With \p conservative the per-axis numerical error bound is filled in
         * too. It is templated rather than always computed because these kernels
         * are register-bound -- the double instantiations sit at 164-220
         * registers with no spills -- so three extra live values are not free for
         * a path that never reads them.
         */
        template <bool is_vf, bool conservative, typename TS, typename Vec4, typename I>
        static inline __device__ void load_query_and_tol(const int qid,
                                                         const I* const SCCD_RESTRICT overlap0,
                                                         const I* const SCCD_RESTRICT overlap1,
                                                         TS** const SCCD_RESTRICT sp,
                                                         TS** const SCCD_RESTRICT ep,
                                                         const size_t element_stride,
                                                         I** const SCCD_RESTRICT elements,
                                                         const TC tol,
                                                         Vec4& sx,
                                                         Vec4& sy,
                                                         Vec4& sz,
                                                         Vec4& ex,
                                                         Vec4& ey,
                                                         Vec4& ez,
                                                         TC* const SCCD_RESTRICT atol,
                                                         TC* const SCCD_RESTRICT aerr) {
            // TS is the storage type and Vec4 the double vector: the loads below
            // widen float geometry on the way into registers, which is exact.
            if constexpr (is_vf) {
                load_query_vf<TS, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, sx, sy, sz, ex, ey, ez);
                compute_face_vertex_tolerance<TC, Vec4>(tol, sx, sy, sz, ex, ey, ez, &atol[0], &atol[1], &atol[2]);
            } else {
                load_query_ee<TS, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, sx, sy, sz, ex, ey, ez);
                compute_edge_edge_tolerance<TC, Vec4>(tol, sx, sy, sz, ex, ey, ez, &atol[0], &atol[1], &atol[2]);
            }

            // Both paths need this, not just the conservative one. Mode 0's
            // origin-containment test pads with the caller's distance tolerance,
            // which is only wide enough when that tolerance happens to exceed the
            // certified bound; see evaluate_cell_3d.
            aerr[0] = numerical_error_bound_component<is_vf, TC, Vec4>(sx, ex);
            aerr[1] = numerical_error_bound_component<is_vf, TC, Vec4>(sy, ey);
            aerr[2] = numerical_error_bound_component<is_vf, TC, Vec4>(sz, ez);
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

        /**
         * \brief DFS stack buffers that persist between narrow-phase calls.
         *
         * cudaMalloc and cudaFree synchronize the whole device. Allocating and
         * releasing the eight stack arrays plus the two counters on every
         * invocation put twenty such stalls into each CCD step, and because the
         * capacity was a local it always restarted from zero and re-grew.
         * Keeping them alive across calls also preserves the capacity the
         * previous step converged on, so steady state does no allocation at all.
         *
         * One instance per scalar type; like the rest of the device pipeline it
         * assumes a single host thread drives it.
         */
        template <typename T>
        struct PersistentDfsStack {
            T* tlower{nullptr};
            T* tupper{nullptr};
            T* ulower{nullptr};
            T* uupper{nullptr};
            T* vlower{nullptr};
            T* vupper{nullptr};
            int* level{nullptr};
            int* qid{nullptr};
            int* counters{nullptr};       // device, [0]/[1] = write cursors for the
                                          // two buffers, [2] = g_request
            int* host_counters{nullptr};  // pinned staging for reading them back
            int cap{0};

            ~PersistentDfsStack() { release(); }

            void release() {
                cudaFree(tlower);
                cudaFree(tupper);
                cudaFree(ulower);
                cudaFree(uupper);
                cudaFree(vlower);
                cudaFree(vupper);
                cudaFree(level);
                cudaFree(qid);
                cudaFree(counters);
                cudaFreeHost(host_counters);
                tlower = tupper = ulower = uupper = vlower = vupper = nullptr;
                level = qid = counters = nullptr;
                host_counters = nullptr;
                cap = 0;
            }
        };

        template <typename T>
        static PersistentDfsStack<T>& persistent_dfs_stack() {
            static PersistentDfsStack<T> stack;
            return stack;
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



        template <int N>
        struct DfsSplit;

        // A single warp per block. The DFS loop carries a __syncthreads() per
        // iteration, so a block runs at the pace of the deepest subtree any of its
        // threads holds; at N == 32 that barrier degenerates to a warp and the
        // imbalance is bounded by 32 rather than by the block. It also quadruples
        // the number of resident blocks, which matters when 220-240 registers
        // already cap a 128-thread block at 2 per SM.
        template <>
        struct DfsSplit<32> {
            static constexpr int NT = 2;
            static constexpr int NU = 4;
            static constexpr int NV = 4;
        };

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

        /**
         * \brief Move the block-shared stack's smallest-t box to the top.
         *
         * The search wants the *earliest* time of impact, and the box that could
         * contain it is the one with the smallest `tlower`. A LIFO refill hands an
         * idle thread the most recently pushed box instead, which is whatever the
         * last split produced -- often a deep box of some other query whose t
         * interval starts far above the answer. Expanding that box cannot tighten
         * the bound, and until the bound tightens nothing else prunes.
         *
         * One argmin over at most SCCD_NP_SHARED_STACK_CAP entries, then one swap,
         * makes the next claim take the minimum. It does not make the block fully
         * best-first -- the second and later claimers in the same round still take
         * LIFO order -- but it puts the most promising box in front on every
         * refill, which is the part that moves the bound.
         *
         * Costs one block reduction per refill round against eight corner
         * evaluations per box, so it is cheap in the units that matter here.
         */
        template <int N, typename TC>
        static inline __device__ void promote_min_tlower(TC* SCCD_RESTRICT s_tlower,
                                                         TC* SCCD_RESTRICT s_tupper,
                                                         TC* SCCD_RESTRICT s_ulower,
                                                         TC* SCCD_RESTRICT s_uupper,
                                                         TC* SCCD_RESTRICT s_vlower,
                                                         TC* SCCD_RESTRICT s_vupper,
                                                         int* SCCD_RESTRICT s_level,
                                                         int* SCCD_RESTRICT s_qid,
                                                         const int top,
                                                         int* SCCD_RESTRICT warp_sums) {
            if (top <= 1) return;

            // Argmin over the live entries, one per thread, reduced across the
            // block. Ties break on the lower index so the choice is deterministic.
            int best = -1;
            TC best_t = TC(0);
            for (int i = (int)threadIdx.x; i < top; i += N) {
                const TC t = s_tlower[i];
                if (best < 0 || t < best_t) {
                    best = i;
                    best_t = t;
                }
            }

            const int lane = threadIdx.x & 31;
            for (int o = 16; o > 0; o >>= 1) {
                const TC other_t = __shfl_xor_sync(0xffffffffu, best_t, o);
                const int other = __shfl_xor_sync(0xffffffffu, best, o);
                if (other >= 0 && (best < 0 || other_t < best_t || (other_t == best_t && other < best))) {
                    best = other;
                    best_t = other_t;
                }
            }
            // warp_sums doubles as the cross-warp scratch; it is sized N/32 and is
            // free at this point in the loop.
            const int warp = threadIdx.x >> 5;
            if (lane == 0) warp_sums[warp] = best;
            __syncthreads();
            if (threadIdx.x == 0) {
                int win = -1;
                TC win_t = TC(0);
                for (int w = 0; w < (N >> 5); ++w) {
                    const int cand = warp_sums[w];
                    if (cand < 0) continue;
                    const TC t = s_tlower[cand];
                    if (win < 0 || t < win_t || (t == win_t && cand < win)) {
                        win = cand;
                        win_t = t;
                    }
                }
                const int last = top - 1;
                if (win >= 0 && win != last) {
                    TC tmp;
                    int itmp;
                    tmp = s_tlower[win]; s_tlower[win] = s_tlower[last]; s_tlower[last] = tmp;
                    tmp = s_tupper[win]; s_tupper[win] = s_tupper[last]; s_tupper[last] = tmp;
                    tmp = s_ulower[win]; s_ulower[win] = s_ulower[last]; s_ulower[last] = tmp;
                    tmp = s_uupper[win]; s_uupper[win] = s_uupper[last]; s_uupper[last] = tmp;
                    tmp = s_vlower[win]; s_vlower[win] = s_vlower[last]; s_vlower[last] = tmp;
                    tmp = s_vupper[win]; s_vupper[win] = s_vupper[last]; s_vupper[last] = tmp;
                    itmp = s_level[win]; s_level[win] = s_level[last]; s_level[last] = itmp;
                    itmp = s_qid[win];   s_qid[win]   = s_qid[last];   s_qid[last]   = itmp;
                }
            }
            __syncthreads();
        }

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
#define SCCD_ENABLE_CODOMAIN_SCALING
#ifdef SCCD_ENABLE_CODOMAIN_SCALING
            const T t0 = in.tlower;
            const T t1 = in.tupper;
            const T u0 = in.ulower;
            const T u1 = in.uupper;
            const T v0 = in.vlower;
            const T v1 = in.vupper;
            T wt = T(0);
            T wu = T(0);
            T wv = T(0);

            if constexpr (is_vf) {
                const T tvx = ex.x - sx.x;
                const T tp1x = ex.y - sx.y;
                const T tp2x = ex.z - sx.z;
                const T tp3x = ex.w - sx.w;
                const T jtx0 = tvx - tp1x;
                const T jtxu = tp1x - tp2x;
                const T jtxv = tp1x - tp3x;
                wt = device::max<T>(
                    wt,
                    device::max<T>(
                        device::max<T>(device::abs<T>(jtx0 + jtxu * u0 + jtxv * v0),
                                       device::abs<T>(jtx0 + jtxu * u1 + jtxv * v0)),
                        device::max<T>(device::abs<T>(jtx0 + jtxu * u0 + jtxv * v1),
                                       device::abs<T>(jtx0 + jtxu * u1 + jtxv * v1))));
                const T sux = sx.z - sx.y;
                const T eux = ex.z - ex.y;
                const T svx = sx.w - sx.y;
                const T evx = ex.w - ex.y;
                wu = device::max<T>(wu, device::max<T>(device::abs<T>(sux + (eux - sux) * t0),
                                                       device::abs<T>(sux + (eux - sux) * t1)));
                wv = device::max<T>(wv, device::max<T>(device::abs<T>(svx + (evx - svx) * t0),
                                                       device::abs<T>(svx + (evx - svx) * t1)));

                const T tvy = ey.x - sy.x;
                const T tp1y = ey.y - sy.y;
                const T tp2y = ey.z - sy.z;
                const T tp3y = ey.w - sy.w;
                const T jty0 = tvy - tp1y;
                const T jtyu = tp1y - tp2y;
                const T jtyv = tp1y - tp3y;
                wt = device::max<T>(
                    wt,
                    device::max<T>(
                        device::max<T>(device::abs<T>(jty0 + jtyu * u0 + jtyv * v0),
                                       device::abs<T>(jty0 + jtyu * u1 + jtyv * v0)),
                        device::max<T>(device::abs<T>(jty0 + jtyu * u0 + jtyv * v1),
                                       device::abs<T>(jty0 + jtyu * u1 + jtyv * v1))));
                const T suy = sy.z - sy.y;
                const T euy = ey.z - ey.y;
                const T svy = sy.w - sy.y;
                const T evy = ey.w - ey.y;
                wu = device::max<T>(wu, device::max<T>(device::abs<T>(suy + (euy - suy) * t0),
                                                       device::abs<T>(suy + (euy - suy) * t1)));
                wv = device::max<T>(wv, device::max<T>(device::abs<T>(svy + (evy - svy) * t0),
                                                       device::abs<T>(svy + (evy - svy) * t1)));

                const T tvz = ez.x - sz.x;
                const T tp1z = ez.y - sz.y;
                const T tp2z = ez.z - sz.z;
                const T tp3z = ez.w - sz.w;
                const T jtz0 = tvz - tp1z;
                const T jtzu = tp1z - tp2z;
                const T jtzv = tp1z - tp3z;
                wt = device::max<T>(
                    wt,
                    device::max<T>(
                        device::max<T>(device::abs<T>(jtz0 + jtzu * u0 + jtzv * v0),
                                       device::abs<T>(jtz0 + jtzu * u1 + jtzv * v0)),
                        device::max<T>(device::abs<T>(jtz0 + jtzu * u0 + jtzv * v1),
                                       device::abs<T>(jtz0 + jtzu * u1 + jtzv * v1))));
                const T suz = sz.z - sz.y;
                const T euz = ez.z - ez.y;
                const T svz = sz.w - sz.y;
                const T evz = ez.w - ez.y;
                wu = device::max<T>(wu, device::max<T>(device::abs<T>(suz + (euz - suz) * t0),
                                                       device::abs<T>(suz + (euz - suz) * t1)));
                wv = device::max<T>(wv, device::max<T>(device::abs<T>(svz + (evz - svz) * t0),
                                                       device::abs<T>(svz + (evz - svz) * t1)));
            } else {
                const T a0x = ex.x - sx.x;
                const T a1x = ex.y - sx.y;
                const T b0x = ex.z - sx.z;
                const T b1x = ex.w - sx.w;
                const T jtx0 = a0x - b0x;
                const T jtxu = a1x - a0x;
                const T jtxv = b0x - b1x;
                wt = device::max<T>(
                    wt,
                    device::max<T>(
                        device::max<T>(device::abs<T>(jtx0 + jtxu * u0 + jtxv * v0),
                                       device::abs<T>(jtx0 + jtxu * u1 + jtxv * v0)),
                        device::max<T>(device::abs<T>(jtx0 + jtxu * u0 + jtxv * v1),
                                       device::abs<T>(jtx0 + jtxu * u1 + jtxv * v1))));
                const T sux = sx.y - sx.x;
                const T eux = ex.y - ex.x;
                const T svx = sx.w - sx.z;
                const T evx = ex.w - ex.z;
                wu = device::max<T>(wu, device::max<T>(device::abs<T>(sux + (eux - sux) * t0),
                                                       device::abs<T>(sux + (eux - sux) * t1)));
                wv = device::max<T>(wv, device::max<T>(device::abs<T>(svx + (evx - svx) * t0),
                                                       device::abs<T>(svx + (evx - svx) * t1)));

                const T a0y = ey.x - sy.x;
                const T a1y = ey.y - sy.y;
                const T b0y = ey.z - sy.z;
                const T b1y = ey.w - sy.w;
                const T jty0 = a0y - b0y;
                const T jtyu = a1y - a0y;
                const T jtyv = b0y - b1y;
                wt = device::max<T>(
                    wt,
                    device::max<T>(
                        device::max<T>(device::abs<T>(jty0 + jtyu * u0 + jtyv * v0),
                                       device::abs<T>(jty0 + jtyu * u1 + jtyv * v0)),
                        device::max<T>(device::abs<T>(jty0 + jtyu * u0 + jtyv * v1),
                                       device::abs<T>(jty0 + jtyu * u1 + jtyv * v1))));
                const T suy = sy.y - sy.x;
                const T euy = ey.y - ey.x;
                const T svy = sy.w - sy.z;
                const T evy = ey.w - ey.z;
                wu = device::max<T>(wu, device::max<T>(device::abs<T>(suy + (euy - suy) * t0),
                                                       device::abs<T>(suy + (euy - suy) * t1)));
                wv = device::max<T>(wv, device::max<T>(device::abs<T>(svy + (evy - svy) * t0),
                                                       device::abs<T>(svy + (evy - svy) * t1)));

                const T a0z = ez.x - sz.x;
                const T a1z = ez.y - sz.y;
                const T b0z = ez.z - sz.z;
                const T b1z = ez.w - sz.w;
                const T jtz0 = a0z - b0z;
                const T jtzu = a1z - a0z;
                const T jtzv = b0z - b1z;
                wt = device::max<T>(
                    wt,
                    device::max<T>(
                        device::max<T>(device::abs<T>(jtz0 + jtzu * u0 + jtzv * v0),
                                       device::abs<T>(jtz0 + jtzu * u1 + jtzv * v0)),
                        device::max<T>(device::abs<T>(jtz0 + jtzu * u0 + jtzv * v1),
                                       device::abs<T>(jtz0 + jtzu * u1 + jtzv * v1))));
                const T suz = sz.y - sz.x;
                const T euz = ez.y - ez.x;
                const T svz = sz.w - sz.z;
                const T evz = ez.w - ez.z;
                wu = device::max<T>(wu, device::max<T>(device::abs<T>(suz + (euz - suz) * t0),
                                                       device::abs<T>(suz + (euz - suz) * t1)));
                wv = device::max<T>(wv, device::max<T>(device::abs<T>(svz + (evz - svz) * t0),
                                                       device::abs<T>(svz + (evz - svz) * t1)));
            }

            const T inv_total = T(1) / (wt + wu + wv + T(1e-16));
            wt = device::max<T>(T(1e-4), wt * inv_total);
            wu = device::max<T>(T(1e-4), wu * inv_total);
            wv = device::max<T>(T(1e-4), wv * inv_total);

            const T dt = (in.tupper - in.tlower) * wt;
            const T du = (in.uupper - in.ulower) * wu;
            const T dv = (in.vupper - in.vlower) * wv;
            const int split_dim = (du > dt && du >= dv) ? 1 : ((dv > dt && dv > du) ? 2 : 0);
#else
            const int split_dim = widest_dimension<T>(in);
#endif
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

        /**
         * \brief Evaluate one box under whichever acceptance test the kernel was
         *        instantiated for.
         *
         * One seam rather than an `if constexpr` at each of the five call sites,
         * so the two predicates cannot drift apart in how they are invoked.
         */
        template <bool is_vf, bool conservative, typename T, typename Vec4>
        static inline __device__ void evaluate_cell_3d_policy(const Domain<T>& cell,
                                                              const Vec4 sx,
                                                              const Vec4 sy,
                                                              const Vec4 sz,
                                                              const Vec4 ex,
                                                              const Vec4 ey,
                                                              const Vec4 ez,
                                                              const T tol,
                                                              const T* const SCCD_RESTRICT atol,
                                                              const T* const SCCD_RESTRICT aerr,
                                                              int& contains_origin,
                                                              int& accept) {
            SCCD_NP_EVAL_TICK();
            if constexpr (conservative) {
                evaluate_cell_3d_tight<is_vf, T, Vec4>(
                    cell, sx, sy, sz, ex, ey, ez, atol, aerr, contains_origin, accept);
            } else {
                evaluate_cell_3d<is_vf, T, Vec4>(
                    cell, sx, sy, sz, ex, ey, ez, tol, atol, aerr, contains_origin, accept);
            }
        }

        /**
         * \brief Split one box under whichever rule the kernel was instantiated for.
         *
         * Returns false only in the conservative case, when the interval has hit
         * the resolution of the format and the caller must accept rather than
         * split. The mode-0 splitters always produce two children.
         */
        template <bool is_vf, bool conservative, typename T, typename Vec4>
        static inline __device__ bool split_cell_policy(const Domain<T>& cur,
                                                        const Vec4 sx,
                                                        const Vec4 sy,
                                                        const Vec4 sz,
                                                        const Vec4 ex,
                                                        const Vec4 ey,
                                                        const Vec4 ez,
                                                        const T* const SCCD_RESTRICT atol,
                                                        Domain<T>& left,
                                                        Domain<T>& right) {
            if constexpr (conservative) {
                return bisect_ti_axis<T>(cur, atol, left, right);
            } else {
                if (SCCD_CUDA_ADAPTIVE_SPLIT) {
                    adaptive_split_longest_axis<is_vf, T, Vec4>(cur, sx, sy, sz, ex, ey, ez, left, right);
                } else {
                    bisect_longest_axis<T>(cur, atol, left, right);
                }
                return true;
            }
        }

        template <bool is_vf, bool conservative, int N, typename T, typename I>
        static __device__ void narrow_phase_dfs_zero_stride_body(const I* const SCCD_RESTRICT overlap0,
                                                                 const I* const SCCD_RESTRICT overlap1,
                                                                 T** const SCCD_RESTRICT sp,
                                                                 T** const SCCD_RESTRICT ep,
                                                                 const size_t element_stride,
                                                                 I** const SCCD_RESTRICT elements,
                                                                 const T tol,
                                                                 const int max_depth,
                                                                 T* SCCD_RESTRICT toi,
                                                                 Stack<TC> g_out,
                                                                 int qid_in,
                                                                 Domain<TC> cur_in,
                                                                 int level_in,
                                                                 int active_in) {
            static_assert(N == 32 || N == 64 || N == 128 || N == 256,
                          "SCCD_NP_THREADS_PER_BLOCK must be one of 32/64/128/256");
            using Vec4 = typename device::Vec4Type<TC>::type;
            constexpr int S_CAP = SharedStackCap<T>::value;
            __shared__ TC s_tlower[S_CAP];
            __shared__ TC s_tupper[S_CAP];
            __shared__ TC s_ulower[S_CAP];
            __shared__ TC s_uupper[S_CAP];
            __shared__ TC s_vlower[S_CAP];
            __shared__ TC s_vupper[S_CAP];
            __shared__ int s_level[S_CAP];
            __shared__ int s_qid[S_CAP];
            __shared__ int s_top;
            __shared__ TC s_toi;
            __shared__ int warp_sums[N >> 5];

            const int tid = threadIdx.x;

            if (tid == 0) {
                s_top = 0;
                s_toi = (TC)toi[0];
            }
            __syncthreads();

            Stack<TC> s_stack = {s_tlower,
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
            Domain<TC> cur = cur_in;
            int level = level_in;
            int active = active_in;
            Vec4 sx, sy, sz, ex, ey, ez;
            TC atol[3] = {TC(0), TC(0), TC(0)};
            // The certified error bound, per axis. Both acceptance policies read
            // it now: the conservative one for all three of its conditions, mode 0
            // as the floor under its rejection pad.
            TC aerr[3] = {TC(0), TC(0), TC(0)};

            if (active) {
                load_query_and_tol<is_vf, conservative, T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol, aerr);
            }

            while (true) {
                // Refreshed every iteration, deliberately. Every resident block
                // hits this one word atomically, which looks like the kind of
                // contention worth batching away, and is not: refreshing every
                // 4th, 16th, 64th and 256th iteration instead was measured
                // monotonically worse on all three scenes -- 638 -> 3377 ms on
                // cloth-funnel at 256 -- because the pruning a fresh bound buys is
                // worth more than the atomic costs. See wip/ASSESSMENT.md.
                if (tid == 0) {
                    const TC g = atomic_min_toi<T>(&toi[0], s_toi);
                    if (g < s_toi) s_toi = g;
                }
                __syncthreads();

                if (active && cur.tlower >= s_toi) active = 0;

#ifdef SCCD_NP_COUNT_BOXES
                if (active) atomicAdd(&g_np_level[level < 80 ? level : 79], 1ull);
#endif

                if (active && level >= max_depth) {
#ifdef SCCD_NP_COUNT_BOXES
                    atomicAdd(&g_np_depth_accept, 1ull);
#endif
                    // Same guard as the block-per-query body: a box outside the
                    // barycentric simplex holds no contact, so accepting it would
                    // report a time of impact for a point that is not on the
                    // triangle. The test is already padded by the u and v
                    // tolerances, so it cannot reject a box that straddles the
                    // edge.
                    if (is_domain_valid<is_vf>(cur, s_toi, atol)) {
                        device::atomic_min(&s_toi, cur.tlower);
                    }
                    active = 0;
                }

                // The conservative split halves an interval, so it eventually runs
                // out of representable midpoints. That is the same situation as
                // the depth cutoff and gets the same answer: accept the box at its
                // t lower bound, which is early and therefore safe.
                if (conservative && active && !tight_can_split<TC>(cur, atol)) {
                    if (is_domain_valid<is_vf>(cur, s_toi, atol)) {
                        device::atomic_min(&s_toi, cur.tlower);
                    }
                    active = 0;
                }

                Domain<TC> push_box;
                int push_level = 0;
                int will_push = 0;

                if (active) {
                    Domain<TC> left, right;
                    SCCD_NP_PERQ_TICK(qid, 2);
                    cur.tupper = device::min<TC>(cur.tupper, s_toi);

                    // Guarded by the tight_can_split check above, so this cannot fail.
                    split_cell_policy<is_vf, conservative, TC, Vec4>(
                        cur, sx, sy, sz, ex, ey, ez, atol, left, right);
                    int cl = 0, cr = 0, al = 0, ar = 0;
                    evaluate_cell_3d_policy<is_vf, conservative, TC, Vec4>(
                        left, sx, sy, sz, ex, ey, ez, tol, atol, aerr, cl, al);
                    evaluate_cell_3d_policy<is_vf, conservative, TC, Vec4>(
                        right, sx, sy, sz, ex, ey, ez, tol, atol, aerr, cr, ar);

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
#ifdef SCCD_NP_COUNT_BOXES
                    atomicAdd((slot >= 0 && slot < S_CAP) ? &g_np_push_shared : &g_np_push_global, 1ull);
#endif
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
                        push_global<TC>(g_out, push_box, push_level, qid);
                    }
                }

                __syncthreads();

#ifdef SCCD_NP_BEST_FIRST
                // Best-first refill: put the smallest-t box in front before the
                // claims run. See promote_min_tlower.
                {
                    const int n_idle = block_popc<N>(!active, warp_sums);
                    if (n_idle > 0) {
                        promote_min_tlower<N, TC>(s_tlower, s_tupper, s_ulower, s_uupper,
                                                  s_vlower, s_vupper, s_level, s_qid,
                                                  atomicAdd(&s_top, 0), warp_sums);
                    }
                }
#endif

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
                                    load_query_and_tol<is_vf, conservative, T, Vec4, I>(qid,
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
                                                                          atol,
                                                                          aerr);
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

            if (tid == 0) atomic_min_toi<T>(&toi[0], s_toi);
        }

        template <bool is_vf, bool conservative, int N, typename T, typename I>
        __global__ SCCD_NP_LAUNCH_BOUNDS(N)
        void narrow_phase_dfs_zero_stride_kernel(const I* const SCCD_RESTRICT overlap0,
                                                 const I* const SCCD_RESTRICT overlap1,
                                                 T** const SCCD_RESTRICT sp,
                                                 T** const SCCD_RESTRICT ep,
                                                 const size_t element_stride,
                                                 I** const SCCD_RESTRICT elements,
                                                 const T tol,
                                                 const int max_depth,
                                                 T* SCCD_RESTRICT toi,
                                                 Stack<TC> g_out,
                                                 const int seed_begin,
                                                 const int seed_end) {
            using Vec4 = typename device::Vec4Type<TC>::type;

            const int tid = threadIdx.x;
            const int my_seed = seed_begin + (int)blockIdx.x * N + tid;
            const bool has_seed = my_seed < seed_end;

            int qid = -1;
            Domain<TC> cur = {TC(0), TC(0), TC(0), TC(0), TC(0), TC(0)};
            int level = 0;
            int active = 0;

            if (has_seed) {
                qid = my_seed;
                Vec4 sx, sy, sz, ex, ey, ez;
                TC atol[3];
                TC aerr[3] = {TC(0), TC(0), TC(0)};
                load_query_and_tol<is_vf, conservative, T, Vec4, I>(
                    qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol, aerr);

                Domain<TC> root = {TC(0), TC(1), TC(0), TC(1), TC(0), TC(1)};
                int contains = 0;
                int accept = 0;
                SCCD_NP_PERQ_TICK(qid, 1);
                evaluate_cell_3d_policy<is_vf, conservative, TC, Vec4>(
                    root, sx, sy, sz, ex, ey, ez, tol, atol, aerr, contains, accept);

                if (contains && is_domain_valid<is_vf>(root, (TC)toi[0], atol)) {
                    cur = root;
                    level = 0;
                    active = 1;
                }
            }

            narrow_phase_dfs_zero_stride_body<is_vf, conservative, N, T, I>(overlap0,
                                                              overlap1,
                                                              sp,
                                                              ep,
                                                              element_stride,
                                                              elements,
                                                              tol,
                                                              max_depth,
                                                              toi,
                                                              g_out,
                                                              qid,
                                                              cur,
                                                              level,
                                                              active);
        }

        template <bool is_vf, bool conservative, int N, typename T, typename I>
        __global__ SCCD_NP_LAUNCH_BOUNDS(N)
        void narrow_phase_dfs_zero_stride_from_stack_kernel(const I* const SCCD_RESTRICT overlap0,
                                                            const I* const SCCD_RESTRICT overlap1,
                                                            T** const SCCD_RESTRICT sp,
                                                            T** const SCCD_RESTRICT ep,
                                                            const size_t element_stride,
                                                            I** const SCCD_RESTRICT elements,
                                                            const T tol,
                                                            const int max_depth,
                                                            T* SCCD_RESTRICT toi,
                                                            Stack<TC> g_in,
                                                            Stack<TC> g_out) {
            int qid = -1;
            Domain<TC> cur = {TC(0), TC(0), TC(0), TC(0), TC(0), TC(0)};
            int level = 0;
            int active = 0;

            if (pop_global<TC>(g_in, cur, level, qid)) {
                active = 1;
            }

            narrow_phase_dfs_zero_stride_body<is_vf, conservative, N, T, I>(overlap0,
                                                              overlap1,
                                                              sp,
                                                              ep,
                                                              element_stride,
                                                              elements,
                                                              tol,
                                                              max_depth,
                                                              toi,
                                                              g_out,
                                                              qid,
                                                              cur,
                                                              level,
                                                              active);
        }

        template <bool is_vf, bool conservative, int N, typename T, typename I>
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
                                                     Stack<TC> g_out,
                                                     const T alpha,
                                                     const int qid,
                                                     Domain<TC> sampling_root,
                                                     const int initial_level,
                                                     const bool do_hard_defer) {
            static_assert(N == 32 || N == 64 || N == 128 || N == 256,
                          "SCCD_NP_THREADS_PER_BLOCK must be one of 32/64/128/256");
            using Vec4 = typename device::Vec4Type<TC>::type;
            constexpr int NT = DfsSplit<N>::NT;
            constexpr int NU = DfsSplit<N>::NU;
            constexpr int NV = DfsSplit<N>::NV;
            constexpr int S_CAP = SharedStackCap<T>::value;
            const int tid = threadIdx.x;

            __shared__ TC s_tlower[S_CAP];
            __shared__ TC s_tupper[S_CAP];
            __shared__ TC s_ulower[S_CAP];
            __shared__ TC s_uupper[S_CAP];
            __shared__ TC s_vlower[S_CAP];
            __shared__ TC s_vupper[S_CAP];
            __shared__ int s_level[S_CAP];
            __shared__ int s_qid[S_CAP];
            __shared__ int s_top;
            __shared__ TC s_toi;
            __shared__ int s_hard;
            __shared__ int s_defer_base;
            __shared__ int s_defer_cursor;
            __shared__ int warp_sums[N >> 5];

            const int toi_idx = qid * toi_stride;

            for (int i = tid; i < S_CAP; i += N) s_qid[i] = SCCD_QID_EMPTY;
            if (tid == 0) {
                s_top = 0;
                s_toi = (TC)toi[toi_idx];
                s_hard = 0;
                s_defer_base = -1;
                s_defer_cursor = 0;
            }
            __syncthreads();

            Stack<TC> s_stack = {s_tlower,
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

            TC atol[3];
            TC aerr[3] = {TC(0), TC(0), TC(0)};
            Vec4 sx, sy, sz, ex, ey, ez;

            load_query_and_tol<is_vf, conservative, T, Vec4, I>(
                qid, overlap0, overlap1, sp, ep, element_stride, elements, tol, sx, sy, sz, ex, ey, ez, atol, aerr);

            const int ti = tid / (NU * NV);
            const int rem = tid % (NU * NV);
            const int ui = rem / NV;
            const int vi = rem % NV;

            // Seeding is the uniform NT*NU*NV dice on both paths.
            //
            // The dice is the wrong shape for the conservative search in
            // principle -- it commits to splitting u and v whether those axes
            // need it, where TightInclusion splits the one axis furthest past its
            // tolerance -- but measurement says that is not what costs. Replacing
            // it with a single root seed, which removes the shape mismatch
            // entirely, moved armadillo-rollers edge-edge from 803 ms to 871 ms:
            // no better, and the opposite sign from a work multiplier. A
            // breadth-first bisection ramp, which fixes the shape *and* fills the
            // block, was worse still and needed a larger shared stack merely to
            // finish.
            //
            // So the cost is in the search loop below, not in how it is seeded.
            // See benchmark/oracle/README.md for where the evidence points next.
            Domain<TC> cur;
            int contains = 0;
            int accept = 0;
            SCCD_NP_PERQ_TICK(qid, 1);
            sample_cell_3d<is_vf, conservative, TC, Vec4>(ti,
                                                          ui,
                                                          vi,
                                                          NT,
                                                          NU,
                                                          NV,
                                                          sampling_root,
                                                          sx,
                                                          sy,
                                                          sz,
                                                          ex,
                                                          ey,
                                                          ez,
                                                          tol,
                                                          atol,
                                                          aerr,
                                                          contains,
                                                          accept,
                                                          cur);

            if (accept && is_domain_valid<is_vf>(cur, s_toi, atol)) {
                device::atomic_min(&s_toi, cur.tlower);
            }
            __syncthreads();

            const int active_seed = contains && !accept && (cur.tlower < s_toi);

            const int co_count = block_popc<N>(active_seed, warp_sums);
            if (!co_count) {
                if (tid == 0) atomic_min_toi<T>(&toi[toi_idx], s_toi);
                return;
            }

            if (do_hard_defer && tid == 0) {
                if ((T)co_count > alpha * (T)N) {
                    s_hard = 1;
                    const int base = atomicAdd(g_out.top, co_count);
                    if (base + co_count > g_out.capacity) {
                        // Deficit will be allocated on the host's retry pass.
                        atomicAdd(g_out.request, co_count);
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
                    if (slot >= 0 && slot < g_out.capacity) {
                        g_out.tlower[slot] = cur.tlower;
                        g_out.tupper[slot] = device::min<TC>(cur.tupper, s_toi);
                        g_out.ulower[slot] = cur.ulower;
                        g_out.uupper[slot] = cur.uupper;
                        g_out.vlower[slot] = cur.vlower;
                        g_out.vupper[slot] = cur.vupper;
                        g_out.level[slot] = initial_level;
                        g_out.qid[slot] = qid;
                    }
                }
                if (tid == 0) atomic_min_toi<T>(&toi[toi_idx], s_toi);
                return;
            }

            int active = active_seed;
            int level = initial_level + 1;
            while (true) {
                if (toi_stride == 0) {
                    if (tid == 0) {
                        const TC g = atomic_min_toi<T>(&toi[toi_idx], s_toi);
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

                // The conservative split halves an interval, so it eventually runs
                // out of representable midpoints. That is the same situation as
                // the depth cutoff and gets the same answer: accept the box at its
                // t lower bound, which is early and therefore safe.
                if (conservative && active && !tight_can_split<TC>(cur, atol)) {
                    if (is_domain_valid<is_vf>(cur, s_toi, atol)) {
                        device::atomic_min(&s_toi, cur.tlower);
                    }
                    active = 0;
                }

                Domain<TC> push_box;
                int push_level = 0;
                int will_push = 0;

                if (active) {
                    Domain<TC> left, right;

                    cur.tupper = device::min<TC>(cur.tupper, s_toi);

                    // Guarded by the tight_can_split check above, so this cannot fail.
                    split_cell_policy<is_vf, conservative, TC, Vec4>(
                        cur, sx, sy, sz, ex, ey, ez, atol, left, right);

                    SCCD_NP_PERQ_TICK(qid, 2);
                    int cl = 0, cr = 0, al = 0, ar = 0;
                    evaluate_cell_3d_policy<is_vf, conservative, TC, Vec4>(
                        left, sx, sy, sz, ex, ey, ez, tol, atol, aerr, cl, al);
                    evaluate_cell_3d_policy<is_vf, conservative, TC, Vec4>(
                        right, sx, sy, sz, ex, ey, ez, tol, atol, aerr, cr, ar);

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
                        push_global<TC>(g_out, push_box, push_level, qid);
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

            if (tid == 0) atomic_min_toi<T>(&toi[toi_idx], s_toi);
        }

        template <bool is_vf, bool conservative, int N, typename T, typename I>
        __global__ SCCD_NP_LAUNCH_BOUNDS(N)
        void narrow_phase_dfs_kernel(const I* const SCCD_RESTRICT overlap0,
                                     const I* const SCCD_RESTRICT overlap1,
                                     T** const SCCD_RESTRICT sp,
                                     T** const SCCD_RESTRICT ep,
                                     const size_t element_stride,
                                     I** const SCCD_RESTRICT elements,
                                     const T tol,
                                     const int max_depth,
                                     T* SCCD_RESTRICT toi,
                                     const int toi_stride,
                                     Stack<TC> g_out,
                                     const T alpha,
                                     const int seed_begin,
                                     const int seed_end) {
            const int qid = seed_begin + (int)blockIdx.x;
            if (qid >= seed_end) return;

            Domain<TC> root = {TC(0), TC(1), TC(0), TC(1), TC(0), TC(1)};
            narrow_phase_dfs_body<is_vf, conservative, N, T, I>(overlap0,
                                                  overlap1,
                                                  sp,
                                                  ep,
                                                  element_stride,
                                                  elements,
                                                  tol,
                                                  max_depth,
                                                  toi,
                                                  toi_stride,
                                                  g_out,
                                                  alpha,
                                                  qid,
                                                  root,
                                                  /*initial_level=*/0,
                                                  /*do_hard_defer=*/true);
        }

        template <bool is_vf, bool conservative, int N, typename T, typename I>
        __global__ SCCD_NP_LAUNCH_BOUNDS(N)
        void narrow_phase_dfs_from_stack_kernel(const I* const SCCD_RESTRICT overlap0,
                                                const I* const SCCD_RESTRICT overlap1,
                                                T** const SCCD_RESTRICT sp,
                                                T** const SCCD_RESTRICT ep,
                                                const size_t element_stride,
                                                I** const SCCD_RESTRICT elements,
                                                const T tol,
                                                const int max_depth,
                                                T* SCCD_RESTRICT toi,
                                                const int toi_stride,
                                                Stack<TC> g_in,
                                                Stack<TC> g_out) {
            __shared__ int b_qid;
            __shared__ int b_level;
            __shared__ int b_have_work;
            __shared__ Domain<TC> b_cur;

            // Drain boxes until the global stack is empty, rather than one per
            // launch.
            //
            // Every round of the host's drain loop costs a kernel launch, a
            // device sync and a blocking readback, and consumes at most
            // base_grid_blocks boxes -- 264 on this part, because occupancy is
            // register-bound at 2 blocks/SM. A conservative search produces far
            // more boxes than that on hard geometry, so the run degenerated into
            // tens of thousands of near-empty rounds: measured at 11.5 s for 39k
            // edge-edge queries on armadillo-rollers, against 26 ms for the
            // mode-0 kernel.
            //
            // This is purely an optimization. A block stops as soon as the stack
            // looks empty, and the host loop still relaunches while g_top > 0, so
            // work pushed by another block after this one gave up is not lost.
            for (int drained = 0; drained < SCCD_NP_DRAIN_PER_BLOCK; ++drained) {
                if (threadIdx.x == 0) {
                    Domain<TC> popped_cur;
                    int popped_level = 0;
                    int popped_qid = -1;
                    if (pop_global<TC>(g_in, popped_cur, popped_level, popped_qid)) {
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

                narrow_phase_dfs_body<is_vf, conservative, N, T, I>(overlap0,
                                                      overlap1,
                                                      sp,
                                                      ep,
                                                      element_stride,
                                                      elements,
                                                      tol,
                                                      max_depth,
                                                      toi,
                                                      toi_stride,
                                                      g_out,
                                                      TC(0),
                                                      b_qid,
                                                      b_cur,
                                                      b_level,
                                                      /*do_hard_defer=*/false);

                // b_cur is about to be overwritten by the next pop.
                __syncthreads();
            }
        }

        template <bool is_vf, bool conservative, typename T, typename I>
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

            // SCCD_READ_ENV stringifies the variable name, so this used to read a
            // bare lowercase `alpha` from the environment -- an unprefixed name in
            // a library's process, which is a collision waiting to happen. Read
            // the prefixed name instead.
            double SCCD_NP_ALPHA = 0.5;
            SCCD_READ_ENV(SCCD_NP_ALPHA, atof);
            const T np_alpha = (T)SCCD_NP_ALPHA;


            // ----------------------------------------------------------------
            // Auto-sized hyperparameters.
            //
            //   SCCD_BLOCKS_PER_SM    -> from CUDA occupancy API
            //   SCCD_GSTACK_CAP_MAX   -> soft cap on a single grow step
            //   SCCD_BATCH_SIZE       -> candidates per outer iteration
            //                            (default: noverlaps)
            //
            // The global stack starts empty and is grown on demand from the
            // deficit the kernel reports. There is no SCCD_GSTACK_CAP_INIT: it
            // was documented here for a while but never read anywhere.
            //
            // Any of these can be overridden via the matching env var.
            // ----------------------------------------------------------------
            int dev = 0;
            cudaGetDevice(&dev);
            cudaDeviceProp prop{};
            cudaGetDeviceProperties(&prop, dev);

            constexpr int N = SCCD_NP_THREADS_PER_BLOCK;
#ifdef SCCD_NP_COUNT_BOXES
            const int S_CAP_DIAG = SharedStackCap<T>::value;
            int drain_rounds_total = 0;
#endif

            int SCCD_BLOCKS_PER_SM = 4;
            {
                int occ = 0;
                if (cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                        &occ, (const void*)narrow_phase_dfs_zero_stride_from_stack_kernel<is_vf, conservative, N, T, I>, N, 0) ==
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

            // Keyed on the compute type, not the interface type: the stack carries
            // box bounds, which the search holds in double whatever the caller's
            // geometry is stored as. A float and a double caller therefore share
            // one persistent stack rather than each growing their own.
            PersistentDfsStack<TC>& gstack = persistent_dfs_stack<TC>();

            if (!gstack.counters) {
                SCCD_CHECK_CUDA(cudaMalloc(&gstack.counters, 4 * sizeof(int)));
                // Pinned staging: these two words are read back once per drain
                // iteration, and a pageable destination forces the driver
                // through an extra bounce buffer each time.
                SCCD_CHECK_CUDA(cudaHostAlloc(&gstack.host_counters, 4 * sizeof(int), cudaHostAllocDefault));
            }
            // Two write cursors, one per buffer, and one deficit counter.
            int* const g_cursor[2] = {gstack.counters, gstack.counters + 1};
            int* const g_request = gstack.counters + 2;

            TC* g_tlower = gstack.tlower;
            TC* g_tupper = gstack.tupper;
            TC* g_ulower = gstack.ulower;
            TC* g_uupper = gstack.uupper;
            TC* g_vlower = gstack.vlower;
            TC* g_vupper = gstack.vupper;
            int* g_level = gstack.level;
            int* g_qid = gstack.qid;
            gstack_cap = gstack.cap;

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
                cudaMalloc(&g_tlower, 2 * (size_t)new_cap * sizeof(TC));
                cudaMalloc(&g_tupper, 2 * (size_t)new_cap * sizeof(TC));
                cudaMalloc(&g_ulower, 2 * (size_t)new_cap * sizeof(TC));
                cudaMalloc(&g_uupper, 2 * (size_t)new_cap * sizeof(TC));
                cudaMalloc(&g_vlower, 2 * (size_t)new_cap * sizeof(TC));
                cudaMalloc(&g_vupper, 2 * (size_t)new_cap * sizeof(TC));
                cudaMalloc(&g_level, 2 * (size_t)new_cap * sizeof(int));
                cudaMalloc(&g_qid, 2 * (size_t)new_cap * sizeof(int));
                gstack_cap = new_cap;
            };

#ifdef SCCD_NP_COUNT_BOXES
            unsigned long long* d_perq = nullptr;
            {
                const unsigned long long zero = 0;
                SCCD_CHECK_CUDA(cudaMemcpyToSymbol(g_np_evals, &zero, sizeof(zero)));
                SCCD_CHECK_CUDA(cudaMemcpyToSymbol(g_np_depth_accept, &zero, sizeof(zero)));
                SCCD_CHECK_CUDA(cudaMemcpyToSymbol(g_np_push_shared, &zero, sizeof(zero)));
                SCCD_CHECK_CUDA(cudaMemcpyToSymbol(g_np_push_global, &zero, sizeof(zero)));
                SCCD_CHECK_CUDA(cudaMemcpyToSymbol(g_np_push_lost, &zero, sizeof(zero)));
                unsigned long long lvl_zero[80] = {0};
                SCCD_CHECK_CUDA(cudaMemcpyToSymbol(g_np_level, lvl_zero, sizeof(lvl_zero)));
                SCCD_CHECK_CUDA(cudaMalloc(&d_perq, noverlaps * sizeof(unsigned long long)));
                SCCD_CHECK_CUDA(cudaMemset(d_perq, 0, noverlaps * sizeof(unsigned long long)));
                SCCD_CHECK_CUDA(cudaMemcpyToSymbol(g_np_perq, &d_perq, sizeof(d_perq)));
            }
#endif

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
            // Pass 2 drain loop on whatever spilled to g_out.  If any
            // push overflowed (g_request > 0) the entire batch is
            // retried after growing the stack -- TOIs are preserved
            // across retries, so subsequent attempts prune more
            // aggressively.
            for (size_t begin = 0; begin < noverlaps; begin += batch_size) {
                const size_t end = (begin + batch_size < noverlaps) ? (begin + batch_size) : noverlaps;
                const int this_batch = (int)(end - begin);

                int retry_rounds = 0;
                while (true) {
                    // g_top and g_request are adjacent ints, so one memset
                    // replaces what was a full kernel launch to zero two words.
                    SCCD_CHECK_CUDA(cudaMemsetAsync(gstack.counters, 0, 3 * sizeof(int)));
#ifdef SCCD_NP_COUNT_BOXES
                    fprintf(stderr,
                            "sccd-np-retry batch=[%zu,%zu) attempt=%d gstack_cap=%d\n",
                            begin,
                            end,
                            retry_rounds,
                            gstack_cap);
#endif

                    // The two halves of the queue. A launch reads one and writes
                    // the other, so `make_buf(w)` is the write buffer for round w
                    // and `make_buf(w ^ 1)` the read buffer, with the roles
                    // swapping each drain round.
                    const auto make_buf = [&](const int half, const int count) {
                        const size_t off = (size_t)half * (size_t)gstack_cap;
                        Stack<TC> b = {g_tlower + off,
                                       g_tupper + off,
                                       g_ulower + off,
                                       g_uupper + off,
                                       g_vlower + off,
                                       g_vupper + off,
                                       g_level + off,
                                       g_qid + off,
                                       g_cursor[half],
                                       g_request,
                                       gstack_cap,
                                       count};
                        return b;
                    };

                    int write_half = 0;
                    Stack<TC> g_out = make_buf(write_half, 0);

                    // Pass 1: seed-driven.
                    if (toi_stride == 0) {
                        const int grid_blocks_zs = (this_batch + N - 1) / N;
                        dim3 grid_pass1_zs(grid_blocks_zs, 1, 1);
                        narrow_phase_dfs_zero_stride_kernel<is_vf, conservative, N, T, I>
                            <<<grid_pass1_zs, block_pass1>>>(overlap0,
                                                             overlap1,
                                                             v0,
                                                             v1,
                                                             element_stride,
                                                             elements,
                                                             tol,
                                                             max_depth,
                                                             d_toi,
                                                             g_out,
                                                             (int)begin,
                                                             (int)end);
                    } else {
                        dim3 grid_pass1(this_batch, 1, 1);
                        narrow_phase_dfs_kernel<is_vf, conservative, SCCD_NP_THREADS_PER_BLOCK, T, I>
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
                                                          g_out,
                                                          np_alpha,
                                                          (int)begin,
                                                          (int)end);
                    }
                    SCCD_CUDA_LAST_ERROR();

                    auto read_cursor = [&](const int half) {
                        SCCD_CHECK_CUDA(cudaMemcpy(
                            gstack.host_counters, g_cursor[half], sizeof(int), cudaMemcpyDeviceToHost));
                        const int n = gstack.host_counters[0];
                        return n > gstack_cap ? gstack_cap : n;
                    };

                    int h_g_top = read_cursor(write_half);

                    // Pass 2: drain what Pass 1 wrote. Each round reads the buffer
                    // the previous round filled and writes the other one, so a
                    // launch never touches a buffer anyone else is touching and
                    // neither side has to wait for the other. Relaunch until the
                    // round produces nothing.
                    while (h_g_top > 0) {
#ifdef SCCD_NP_COUNT_BOXES
                        ++drain_rounds_total;
#endif
                        const Stack<TC> g_in = make_buf(write_half, h_g_top);
                        write_half ^= 1;
                        // The write cursor for the buffer we are about to fill has
                        // to start at zero, and the read cursor for g_in is the
                        // same word we just read the count from -- reset it so the
                        // claims in pop_global start from the first entry.
                        SCCD_CHECK_CUDA(cudaMemsetAsync(g_cursor[write_half], 0, sizeof(int)));
                        SCCD_CHECK_CUDA(cudaMemsetAsync(g_cursor[write_half ^ 1], 0, sizeof(int)));
                        g_out = make_buf(write_half, 0);

                        // The grid must cover every entry in the read buffer. It
                        // used to be capped at base_grid_blocks and the loop
                        // relaunched against the same stack until it emptied; with
                        // two buffers the unread remainder would be discarded at
                        // the swap instead, so the cap has to go. One thread per
                        // entry for the zero-stride kernel; the block-per-query
                        // kernel claims up to SCCD_NP_DRAIN_PER_BLOCK each.
                        long long need = (toi_stride == 0)
                                             ? ((long long)h_g_top + N - 1) / N
                                             : ((long long)h_g_top + SCCD_NP_DRAIN_PER_BLOCK - 1) /
                                                   SCCD_NP_DRAIN_PER_BLOCK;
                        if (need < 1) need = 1;
                        if (need > 2147483647LL) need = 2147483647LL;
                        dim3 grid_pass2((unsigned)need, 1, 1);

                        if (toi_stride == 0) {
                            narrow_phase_dfs_zero_stride_from_stack_kernel<is_vf, conservative, N, T, I><<<grid_pass2, block_pass1>>>(
                                overlap0, overlap1, v0, v1, element_stride, elements, tol, max_depth, d_toi, g_in, g_out);
                        } else {
                            narrow_phase_dfs_from_stack_kernel<is_vf, conservative, N, T, I>
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
                                                              g_in,
                                                              g_out);
                        }
                        SCCD_CUDA_LAST_ERROR();

                        h_g_top = read_cursor(write_half);
                    }

                    // Did anything overflow during Pass 1 or Pass 2?
                    // Grow and retry the entire batch from seeds; the
                    // tighter TOIs from this attempt cut the work tree
                    // on the next pass.
                    SCCD_CHECK_CUDA(
                        cudaMemcpy(gstack.host_counters + 1, g_request, sizeof(int), cudaMemcpyDeviceToHost));
                    const int h_g_request = gstack.host_counters[1];
                    if (h_g_request <= 0) break;

                    // Bound the grow-and-retry. Each round reruns the whole batch,
                    // and the growth step is only as large as the deficit the last
                    // round happened to report, so a search that keeps overflowing
                    // can spin here for an unbounded number of full-batch reruns --
                    // which is how a kernel that produced too many boxes showed up
                    // as a hang rather than as a slow run. Giving up is safe: boxes
                    // dropped on overflow can only leave the time of impact too
                    // large, never too small, and the toi carried over from the
                    // rounds already done is kept.
                    if (++retry_rounds >= SCCD_NP_MAX_RETRY_ROUNDS) {
                        fprintf(stderr,
                                "sccd: narrow phase gave up growing the global stack after %d rounds "
                                "(deficit %d, capacity %d). The result is still conservative, but this "
                                "batch did more work than the stack can hold.\n",
                                retry_rounds,
                                h_g_request,
                                gstack_cap);
                        break;
                    }

                    // printf("Overflowed: %d\n", h_g_request);

                    // Grow by the deficit, deliberately, not geometrically.
                    //
                    // Doubling looks like the textbook fix for repeated reallocs,
                    // and it is a regression here: every call memsets the whole
                    // stack capacity to restore the empty slot marker, so an
                    // oversized stack is paid for on every subsequent call, not
                    // once. Measured on armadillo-rollers edge-edge, doubling cost
                    // 803 -> ~950 ms. The retry count is bounded by
                    // SCCD_NP_MAX_RETRY_ROUNDS instead, which is what actually
                    // needed fixing.
                    int grow_by = h_g_request;
                    if (grow_by > SCCD_GSTACK_CAP_MAX) grow_by = SCCD_GSTACK_CAP_MAX;
                    const long long target_ll = (long long)gstack_cap + (long long)grow_by;
                    const int new_cap = (target_ll > (long long)INT_MAX) ? INT_MAX : (int)target_ll;
                    grow_stack(new_cap);
                }
            }

            // Hand the (possibly grown) buffers back for the next call.
            gstack.tlower = g_tlower;
            gstack.tupper = g_tupper;
            gstack.ulower = g_ulower;
            gstack.uupper = g_uupper;
            gstack.vlower = g_vlower;
            gstack.vupper = g_vupper;
            gstack.level = g_level;
            gstack.qid = g_qid;
            gstack.cap = gstack_cap;

            SCCD_CUDA_LAST_ERROR();

#ifdef SCCD_NP_COUNT_BOXES
            {
                unsigned long long evals = 0;
                SCCD_CHECK_CUDA(cudaDeviceSynchronize());
                SCCD_CHECK_CUDA(cudaMemcpyFromSymbol(&evals, g_np_evals, sizeof(evals)));
                fprintf(stderr,
                        "sccd-np-count %s %s stride=%d queries=%zu corner_evals=%llu per_query=%.1f\n",
                        is_vf ? "vf" : "ee",
                        conservative ? "conservative" : "mode0",
                        toi_stride,
                        noverlaps,
                        evals,
                        noverlaps ? (double)evals / (double)noverlaps : 0.0);

                std::vector<unsigned long long> perq(noverlaps, 0ull);
                SCCD_CHECK_CUDA(cudaMemcpy(perq.data(),
                                           d_perq,
                                           noverlaps * sizeof(unsigned long long),
                                           cudaMemcpyDeviceToHost));
                unsigned long long hist[24] = {0}, worst = 0;
                size_t worst_q = 0;
                for (size_t q = 0; q < noverlaps; ++q) {
                    unsigned long long n = perq[q];
                    if (n > worst) {
                        worst = n;
                        worst_q = q;
                    }
                    int b = 0;
                    while (n > 1 && b < 23) {
                        n >>= 1;
                        ++b;
                    }
                    ++hist[b];
                }
                fprintf(stderr,
                        "sccd-np-hist device stride=%d queries=%zu worst=%llu at=%zu hist=",
                        toi_stride,
                        noverlaps,
                        worst,
                        worst_q);
                for (int b = 0; b < 24; ++b) fprintf(stderr, "%llu%s", hist[b], b == 23 ? "\n" : ",");

                unsigned long long lvl[80] = {0}, depth_acc = 0;
                SCCD_CHECK_CUDA(cudaMemcpyFromSymbol(lvl, g_np_level, sizeof(lvl)));
                SCCD_CHECK_CUDA(cudaMemcpyFromSymbol(&depth_acc, g_np_depth_accept, sizeof(depth_acc)));
                unsigned long long ps = 0, pg = 0, pl = 0;
                SCCD_CHECK_CUDA(cudaMemcpyFromSymbol(&ps, g_np_push_shared, sizeof(ps)));
                SCCD_CHECK_CUDA(cudaMemcpyFromSymbol(&pg, g_np_push_global, sizeof(pg)));
                SCCD_CHECK_CUDA(cudaMemcpyFromSymbol(&pl, g_np_push_lost, sizeof(pl)));
                fprintf(stderr,
                        "sccd-np-push stride=%d s_cap=%d shared=%llu global=%llu lost=%llu drains=%d\n",
                        toi_stride,
                        S_CAP_DIAG,
                        ps,
                        pg,
                        pl,
                        drain_rounds_total);
                fprintf(stderr, "sccd-np-level device depth_accept=%llu levels=", depth_acc);
                for (int b = 0; b < 80; ++b) fprintf(stderr, "%llu%s", lvl[b], b == 79 ? "\n" : ",");

                const unsigned long long* np = nullptr;
                SCCD_CHECK_CUDA(cudaMemcpyToSymbol(g_np_perq, &np, sizeof(np)));
                SCCD_CHECK_CUDA(cudaFree(d_perq));
            }
#endif
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
            if (narrow_phase_mode_is_tight(narrow_phase_mode())) {
                return narrow_phase_generic<false, true, T, I>(
                    noverlaps, overlap0, overlap1, v0, v1, edge_stride, edges, max_toi, toi, max_depth, tol, toi_stride);
            }
            return narrow_phase_generic<false, false, T, I>(
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
            if (narrow_phase_mode_is_tight(narrow_phase_mode())) {
                return narrow_phase_generic<true, true, T, I>(
                    noverlaps, voveralp, foveralp, v0, v1, face_stride, faces, max_toi, toi, max_depth, tol, toi_stride);
            }
            return narrow_phase_generic<true, false, T, I>(
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

#define SCCD_NP_NARROW_PHASE_EE(T, I)                                                   \
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

#define SCCD_NP_NARROW_PHASE_VF(NXE, T, I)                                                   \
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

SCCD_NP_NARROW_PHASE_EE(float, int32_t);
SCCD_NP_NARROW_PHASE_EE(float, int64_t);
SCCD_NP_NARROW_PHASE_EE(double, int32_t);
SCCD_NP_NARROW_PHASE_EE(double, int64_t);

SCCD_NP_NARROW_PHASE_VF(3, float, int32_t);
SCCD_NP_NARROW_PHASE_VF(3, float, int64_t);
SCCD_NP_NARROW_PHASE_VF(3, double, int32_t);
SCCD_NP_NARROW_PHASE_VF(3, double, int64_t);

template int sccd::device::minmax<float>(const float* const SCCD_RESTRICT data,
                                         const size_t n,
                                         float* const h_min,
                                         float* const h_max);
template int sccd::device::minmax<double>(const double* const SCCD_RESTRICT data,
                                          const size_t n,
                                          double* const h_min,
                                          double* const h_max);

#undef SCCD_NP_NARROW_PHASE_EE
#undef SCCD_NP_NARROW_PHASE_VF
