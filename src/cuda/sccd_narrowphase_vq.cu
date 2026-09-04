#include "sccd_narrowphase_vq.cuh"

#include "sccd_cuda_base.cuh"
#include "sccd_tolerance.hpp"

#include <cfloat>
#include <cstdint>
#include <cstdio>

namespace sccd {
    namespace device {

        namespace {

            // The search computes in double whatever the geometry is stored as.
            // See the note in the header: in single precision the certified error
            // bound and the tolerances that terminate the search are close enough
            // together that the guarantee does not survive.
            using TC = double;

            // Splits per axis, matching SCCD_ADAPTIVE_NUM_SPLITS on the host so the two
            // explore the same tree.
            constexpr int kSplits = 2;
            constexpr int kSamples = kSplits + 2;

            // Per-thread stack, sized so that overflow cannot happen.
            //
            // Each pop pushes at most kSplits + 1 children, so the stack grows by
            // at most kSplits per level and a depth-D search needs 1 + kSplits * D
            // entries. The cap and the depth limit are therefore tied together:
            // kMaxDeviceDepth below is derived from this number, not chosen.
            //
            // Sizing this by intuition instead cost accuracy, and it is worth
            // being exact about what it did and did not cost.
            //
            // The host averages about 1.2 boxes live per query when a shared time
            // of impact is pruning and about 11 when it is not, so 32 looked
            // ample. It was not. Overflow accepts, and accepting is safe under
            // any circumstances: the reported time is the box's t lower bound,
            // which is at or before any root inside it. So a small cap cannot
            // miss a collision -- conservativeness here is structural and a cap
            // does not touch it.
            //
            // What it does is accept boxes that hold no root, which is a false
            // positive, and report a time earlier than the true one. Both are on
            // the acceptable side of the asymmetry, so this was an accuracy
            // regression rather than a correctness one: the answer was 1.13e-3
            // before the true value where the triangle path manages 3.0e-7, which
            // costs a solver step size and never a contact. Raising the cap walked
            // it back exactly:
            //
            //     cap  32 -> 0.6655369599
            //     cap 128 -> 0.6666666238   (the host's answer)
            //     cap 512 -> 0.6666666238
            //
            // ## Why this array is large, and why that is free
            //
            // 257 entries is a 14.7 KB per-thread stack frame, which looks
            // alarming and measures as nothing. Only the first few entries are
            // ever touched -- the rest is address space, not traffic -- and on
            // GH200, 400k queries at depth 69:
            //
            //     cap   frame     registers  spill      stride 1   stride 0
            //     140   8128 B    255        112/216    7.589 ms   0.904 ms
            //     257  14688 B    255        112/216    7.441 ms   0.893 ms
            //     513  29024 B    255        112/216    7.355 ms   0.922 ms
            //
            // The register count and the spill are the inclusion function's
            // working set -- 30 coordinates, two Frames, eight corner triples --
            // and do not move with this array at all: they are 255 and 112/216
            // from cap 4 to cap 513.
            //
            // What the cap does change is the depth limit derived from it, and
            // that is worth a great deal: at a matched depth the array size is
            // invisible (cap 140 and cap 32 both run 3.7 ms at depth 15, to the
            // same answer), while cap 32's implied depth of 15 costs 5.5e-5 of
            // mean accuracy and cap 8's implied depth of 3 collapses the shared
            // minimum to zero.
            //
            // So the cap buys search depth and headroom costs nothing. It is
            // sized for depth 128 rather than for the default 69, which is what
            // lets a caller raise SCCD_MAX_DEPTH and actually get it.
#ifndef SCCD_VQ_MAX_DEPTH
#define SCCD_VQ_MAX_DEPTH 128
#endif
            constexpr int kStackCap = 1 + kSplits * SCCD_VQ_MAX_DEPTH;

            // The deepest search this stack can hold without overflowing.
            constexpr int kMaxDeviceDepth = (kStackCap - 1) / kSplits;

            struct Domain {
                TC tlo, thi, ulo, uhi, vlo, vhi;
                int depth;
            };

            struct Frame {
                TC V[3];
                TC P[4][3];
            };

            // ----------------------------------------------------------------
            // The inclusion function, factored over the quad's parameter domain.
            //
            //   F(t,u,v) = V(t) - SUM_k w_k(u,v) P_k(t)
            //
            // The weights do not depend on t, so P_k and V can be computed once
            // per t and reused across the (u, v) corners. Same regrouping as the
            // host, and the host test checks it against the reference form.
            // ----------------------------------------------------------------
            static inline __device__ void frame_at(const TC sv[3],
                                                   const TC s[4][3],
                                                   const TC ev[3],
                                                   const TC e[4][3],
                                                   const TC t,
                                                   Frame& f) {
                const TC omt = TC(1) - t;
#pragma unroll
                for (int d = 0; d < 3; ++d) f.V[d] = omt * sv[d] + t * ev[d];
#pragma unroll
                for (int k = 0; k < 4; ++k) {
#pragma unroll
                    for (int d = 0; d < 3; ++d) f.P[k][d] = omt * s[k][d] + t * e[k][d];
                }
            }

            static inline __device__ void eval_frame(const Frame& f, const TC u, const TC v, TC F[3]) {
                const TC omu = TC(1) - u;
                const TC omv = TC(1) - v;
                const TC w0 = omu * omv;
                const TC w1 = u * omv;
                const TC w2 = omu * v;
                const TC w3 = u * v;
#pragma unroll
                for (int d = 0; d < 3; ++d) {
                    F[d] = f.V[d] - (w0 * f.P[0][d] + w1 * f.P[1][d] + w2 * f.P[2][d] + w3 * f.P[3][d]);
                }
            }

            static inline __device__ void diff_vq(const TC sv[3],
                                                  const TC s[4][3],
                                                  const TC ev[3],
                                                  const TC e[4][3],
                                                  const TC t,
                                                  const TC u,
                                                  const TC v,
                                                  TC F[3]) {
                Frame f;
                frame_at(sv, s, ev, e, t, f);
                eval_frame(f, u, v, F);
            }

            // ----------------------------------------------------------------
            // Per-query reductions: the Lipschitz constants (which are also the
            // codomain widths), the per-axis tolerances derived from them, and the
            // certified numerical error bound. One pass, as on the host.
            // ----------------------------------------------------------------
            static inline __device__ void prepare(const TC codomain_tol,
                                                  const TC sv[3],
                                                  const TC s[4][3],
                                                  const TC ev[3],
                                                  const TC e[4][3],
                                                  TC tols[3],
                                                  TC widths[3],
                                                  TC err[3]) {
                TC lip[3] = {TC(0), TC(0), TC(0)};
                TC maxc[3];

#pragma unroll
                for (int d = 0; d < 3; ++d) {
                    const TC vt = ev[d] - sv[d];

                    TC l0 = TC(0), l1 = TC(0), l2 = TC(0);
#pragma unroll
                    for (int k = 0; k < 4; ++k) {
                        l0 = fmax(l0, fabs(vt - (e[k][d] - s[k][d])));
                    }
                    l1 = fmax(fmax(fabs(s[1][d] - s[0][d]), fabs(s[3][d] - s[2][d])),
                              fmax(fabs(e[1][d] - e[0][d]), fabs(e[3][d] - e[2][d])));
                    l2 = fmax(fmax(fabs(s[2][d] - s[0][d]), fabs(s[3][d] - s[1][d])),
                              fmax(fabs(e[2][d] - e[0][d]), fabs(e[3][d] - e[1][d])));

                    lip[0] = fmax(lip[0], l0);
                    lip[1] = fmax(lip[1], l1);
                    lip[2] = fmax(lip[2], l2);

                    TC m = fmax(fabs(sv[d]), fabs(ev[d]));
#pragma unroll
                    for (int k = 0; k < 4; ++k) {
                        m = fmax(m, fmax(fabs(s[k][d]), fabs(e[k][d])));
                    }
                    // fmin, not fmax: TightInclusion's bound is
                    // `filter * min(max_coord, 1)^3`, so the cube is 1 on any
                    // scene at unit scale or larger. Clamping the other way
                    // grows the pad as the cube of the scene size.
                    maxc[d] = fmin(TC(1), m);
                }

                const TC axis_tol = codomain_tol / TC(3);
                const TC eps = DBL_EPSILON;
                // 38 eps max_coord^3 for the five-point vertex-quad form, matching
                // sccd_get_numerical_error_vq_soa on the host. This is the bound
                // the rejection pads by, and substituting the user's tolerance for
                // it is exactly the unsound rejection this project has already
                // found once.
                const TC filter = TC(38) * eps;

                // The same caps as the host's compute_vertex_quad_tolerance and
                // vq_prepare: axis 0 is time, axes 1 and 2 are the quad's
                // parameters. Without them the division grows without bound as
                // the motion slows, and below epsilon it fell back to 1.0 -- a
                // tolerance as wide as the whole domain.
                // Spelled out rather than calling sccd::clamp_domain_tol, which
                // is host-only; this is the same expression, and the same
                // clamp_tol lambda the triangle kernel uses.
                const TC caps[3] = {TC(SCCD_MAX_TIME_TOL), TC(SCCD_MAX_COORD_TOL),
                                    TC(SCCD_MAX_COORD_TOL)};
                auto clamp_tol = [](const TC v, const TC cap) { return (v < cap) ? v : cap; };
#pragma unroll
                for (int d = 0; d < 3; ++d) {
                    const TC raw = lip[d] > eps ? axis_tol / lip[d] : caps[d];
                    tols[d] = clamp_tol(raw, caps[d]);
                    widths[d] = lip[d];
                    err[d] = maxc[d] * maxc[d] * maxc[d] * filter;
                }

                // Normalise exactly as normalize_vertex_quad_codomain_widths does:
                // share of the total, floored at 1e-4 so a degenerate axis still
                // gets split eventually rather than never being chosen.
                const TC total = widths[0] + widths[1] + widths[2] + eps;
#pragma unroll
                for (int d = 0; d < 3; ++d) widths[d] = fmax(TC(1e-4), widths[d] / total);
            }

            static inline __device__ bool codomain_acceptance(const TC fmin[3],
                                                              const TC fmax_[3],
                                                              const TC tol,
                                                              const TC tols[3],
                                                              const TC err[3],
                                                              bool& accept) {
                bool contains_zero = true;
                bool smaller_than_axis_tol = true;
                bool inside_epsilon_box = true;
                bool smaller_than_scalar_tol = true;
                bool degenerate_interval = true;

#pragma unroll
                for (int d = 0; d < 3; ++d) {
                    const TC width = fmax_[d] - fmin[d];
                    contains_zero = contains_zero && (fmin[d] <= err[d]) && (fmax_[d] >= -err[d]);
                    smaller_than_axis_tol = smaller_than_axis_tol && (width <= tols[d]);
                    inside_epsilon_box =
                        inside_epsilon_box && (fmin[d] >= -err[d]) && (fmax_[d] <= err[d]);
                    smaller_than_scalar_tol = smaller_than_scalar_tol && (width < tol);
                    degenerate_interval = degenerate_interval && (fmin[d] >= fmax_[d]);
                }

                accept = contains_zero && (smaller_than_axis_tol || inside_epsilon_box ||
                                           smaller_than_scalar_tol || degenerate_interval);
                return contains_zero;
            }

            // The Gauss-Newton splitters, ported from the host so the two explore
            // the same tree rather than merely the same answer.
            template <int SplitDim>
            static inline __device__ void axis_splitters(const Domain& box,
                                                         const TC sv[3],
                                                         const TC s[4][3],
                                                         const TC ev[3],
                                                         const TC e[4][3],
                                                         TC splitters[kSplits]) {
                const TC lo = (SplitDim == 0) ? box.tlo : (SplitDim == 1 ? box.ulo : box.vlo);
                const TC hi = (SplitDim == 0) ? box.thi : (SplitDim == 1 ? box.uhi : box.vhi);
                const TC h = (hi - lo) / TC(kSplits + 1);
                const TC radius = h * TC(0.45);

                const TC mid_t = (box.tlo + box.thi) * TC(0.5);
                const TC mid_u = (box.ulo + box.uhi) * TC(0.5);
                const TC mid_v = (box.vlo + box.vhi) * TC(0.5);

                TC F_base[3];
                const TC base_t = (SplitDim == 0) ? TC(0) : mid_t;
                const TC base_u = (SplitDim == 1) ? TC(0) : mid_u;
                const TC base_v = (SplitDim == 2) ? TC(0) : mid_v;
                diff_vq(sv, s, ev, e, base_t, base_u, base_v, F_base);

                const TC omt = TC(1) - mid_t;
                const TC omu = TC(1) - mid_u;
                const TC omv = TC(1) - mid_v;

                TC J[3];
                TC H = TC(0);
#pragma unroll
                for (int d = 0; d < 3; ++d) {
                    if (SplitDim == 0) {
                        const TC sq = omu * omv * s[0][d] + mid_u * omv * s[1][d] +
                                      omu * mid_v * s[2][d] + mid_u * mid_v * s[3][d];
                        const TC eq = omu * omv * e[0][d] + mid_u * omv * e[1][d] +
                                      omu * mid_v * e[2][d] + mid_u * mid_v * e[3][d];
                        J[d] = (ev[d] - sv[d]) - (eq - sq);
                    } else if (SplitDim == 1) {
                        const TC du_s = omv * (s[1][d] - s[0][d]) + mid_v * (s[3][d] - s[2][d]);
                        const TC du_e = omv * (e[1][d] - e[0][d]) + mid_v * (e[3][d] - e[2][d]);
                        J[d] = -(omt * du_s + mid_t * du_e);
                    } else {
                        const TC dv_s = omu * (s[2][d] - s[0][d]) + mid_u * (s[3][d] - s[1][d]);
                        const TC dv_e = omu * (e[2][d] - e[0][d]) + mid_u * (e[3][d] - e[1][d]);
                        J[d] = -(omt * dv_s + mid_t * dv_e);
                    }
                    H += J[d] * J[d];
                }

                const TC step_scale = H > DBL_EPSILON ? TC(1) / H : TC(0.00001);
#pragma unroll
                for (int i = 0; i < kSplits; ++i) {
                    const TC x0 = lo + h * TC(i + 1);
                    const TC xmin = fmax(lo, x0 - radius);
                    const TC xmax = fmin(hi, x0 + radius);
                    TC g = TC(0);
#pragma unroll
                    for (int d = 0; d < 3; ++d) g += (F_base[d] + x0 * J[d]) * J[d];
                    splitters[i] = fmin(xmax, fmax(xmin, x0 - g * step_scale));
                }
            }

            // Conservative atomic minimum on the shared time of impact. CUDA has no
            // atomicMin for double, so this is the usual compare-and-swap loop.
            static inline __device__ double atomic_min_double(double* addr, const double value) {
                unsigned long long int* as_ull = (unsigned long long int*)addr;
                unsigned long long int old = *as_ull;
                unsigned long long int assumed;
                do {
                    const double current = __longlong_as_double(old);
                    if (current <= value) return current;
                    assumed = old;
                    old = atomicCAS(as_ull, assumed, __double_as_longlong(value));
                } while (assumed != old);
                return __longlong_as_double(old);
            }

            // ----------------------------------------------------------------
            // One thread per query.
            // ----------------------------------------------------------------
            template <typename T, typename I>
            __global__ void narrow_phase_vq_kernel(const size_t noverlaps,
                                                   const I* const SCCD_RESTRICT voverlap,
                                                   const I* const SCCD_RESTRICT qoverlap,
                                                   T* const SCCD_RESTRICT v0x,
                                                   T* const SCCD_RESTRICT v0y,
                                                   T* const SCCD_RESTRICT v0z,
                                                   T* const SCCD_RESTRICT v1x,
                                                   T* const SCCD_RESTRICT v1y,
                                                   T* const SCCD_RESTRICT v1z,
                                                   const size_t quad_stride,
                                                   const I* const SCCD_RESTRICT q0,
                                                   const I* const SCCD_RESTRICT q1,
                                                   const I* const SCCD_RESTRICT q2,
                                                   const I* const SCCD_RESTRICT q3,
                                                   const TC max_toi,
                                                   double* const SCCD_RESTRICT shared_toi,
                                                   T* const SCCD_RESTRICT toi,
                                                   const int max_depth,
                                                   const TC tol,
                                                   const int toi_stride) {
                const size_t i = blockIdx.x * (size_t)blockDim.x + threadIdx.x;
                if (i >= noverlaps) return;

                const I vi = voverlap[i];
                const I qi = qoverlap[i];
                const I n[4] = {q0[qi * quad_stride], q1[qi * quad_stride],
                                q2[qi * quad_stride], q3[qi * quad_stride]};

                // Widening float storage to double on load is exact.
                const TC sv[3] = {(TC)v0x[vi], (TC)v0y[vi], (TC)v0z[vi]};
                const TC ev[3] = {(TC)v1x[vi], (TC)v1y[vi], (TC)v1z[vi]};
                TC s[4][3], e[4][3];
#pragma unroll
                for (int k = 0; k < 4; ++k) {
                    s[k][0] = (TC)v0x[n[k]];
                    s[k][1] = (TC)v0y[n[k]];
                    s[k][2] = (TC)v0z[n[k]];
                    e[k][0] = (TC)v1x[n[k]];
                    e[k][1] = (TC)v1y[n[k]];
                    e[k][2] = (TC)v1z[n[k]];
                }

                TC tols[3], widths[3], err[3];
                prepare(tol, sv, s, ev, e, tols, widths, err);

                // Clamp to what the stack can hold. Going deeper than this could
                // only end in an overflow, and an overflow accepts, so the clamp
                // costs nothing a deeper search would have gained -- it just makes
                // the limit explicit and keeps the reported time of impact a
                // property of the depth rather than of an array size.
                const int depth_limit = max_depth < kMaxDeviceDepth ? max_depth : kMaxDeviceDepth;

                TC t = (toi_stride == 0) ? *shared_toi : (TC)max_toi;
                if (t > (TC)max_toi) t = (TC)max_toi;
                const TC t_upper = fmin(t, TC(1));

                Domain stack[kStackCap];
                int top = 0;
                stack[top++] = Domain{TC(0), t_upper, TC(0), TC(1), TC(0), TC(1), 0};

                bool found = false;

                while (top > 0) {
                    Domain box = stack[--top];

                    if (box.tlo >= t) continue;
                    box.thi = fmin(box.thi, t);

                    // Widest axis in the codomain metric. The comparison order is
                    // Box::widest_dimension's, ties included: t wins a tie with u,
                    // and u wins a tie with v. Getting this subtly different would
                    // make the device explore a different tree from the host for no
                    // reason, which is exactly what makes a divergence hard to read.
                    const TC dt = (box.thi - box.tlo) * widths[0];
                    const TC du = (box.uhi - box.ulo) * widths[1];
                    const TC dv = (box.vhi - box.vlo) * widths[2];
                    const int dim = (du > dt && du >= dv) ? 1 : ((dv > dt && dv > du) ? 2 : 0);

                    TC splitters[kSplits];
                    if (dim == 0) {
                        axis_splitters<0>(box, sv, s, ev, e, splitters);
                    } else if (dim == 1) {
                        axis_splitters<1>(box, sv, s, ev, e, splitters);
                    } else {
                        axis_splitters<2>(box, sv, s, ev, e, splitters);
                    }

                    const TC lo = (dim == 0) ? box.tlo : (dim == 1 ? box.ulo : box.vlo);
                    const TC hi = (dim == 0) ? box.thi : (dim == 1 ? box.uhi : box.vhi);
                    TC samples[kSamples];
                    samples[0] = lo;
                    samples[kSamples - 1] = hi;
#pragma unroll
                    for (int k = 0; k < kSplits; ++k) samples[k + 1] = splitters[k];

                    // Sub-boxes in order; the caller's stack is LIFO, so pushing in
                    // reverse would visit the earliest t last. The t-cutoff below
                    // is re-read each iteration, exactly as on the host, because it
                    // tightens whenever a sub-box accepts.
                    TC plane_a[4][3];
                    bool have_a = false;

                    for (int k = 0; k < kSplits + 1; ++k) {
                        const TC smin = samples[k];
                        const TC smax = samples[k + 1];
                        const TC tt_min = (dim == 0) ? smin : box.tlo;
                        if (tt_min >= t) { have_a = false; continue; }

                        // The four corners on each of the two bounding planes.
                        TC plane_b[4][3];
                        const TC a_lo = (dim == 0) ? box.ulo : box.tlo;
                        const TC a_hi = (dim == 0) ? box.uhi : box.thi;
                        const TC b_lo = (dim == 2) ? box.ulo : box.vlo;
                        const TC b_hi = (dim == 2) ? box.uhi : box.vhi;

                        for (int pass = 0; pass < 2; ++pass) {
                            if (pass == 0 && have_a) continue;
                            const TC sample = (pass == 0) ? smin : smax;
                            TC (*dst)[3] = (pass == 0) ? plane_a : plane_b;

                            if (dim == 0) {
                                Frame f;
                                frame_at(sv, s, ev, e, sample, f);
#pragma unroll
                                for (int c = 0; c < 4; ++c) {
                                    eval_frame(f, (c & 1) ? a_hi : a_lo, (c & 2) ? b_hi : b_lo, dst[c]);
                                }
                            } else {
                                Frame flo, fhi;
                                frame_at(sv, s, ev, e, box.tlo, flo);
                                frame_at(sv, s, ev, e, box.thi, fhi);
#pragma unroll
                                for (int c = 0; c < 4; ++c) {
                                    const Frame& f = (c & 1) ? fhi : flo;
                                    const TC other = (c & 2) ? b_hi : b_lo;
                                    const TC uu = (dim == 1) ? sample : other;
                                    const TC vv = (dim == 1) ? other : sample;
                                    eval_frame(f, uu, vv, dst[c]);
                                }
                            }
                        }
                        have_a = true;

                        TC fmin_[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
                        TC fmax_[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
#pragma unroll
                        for (int c = 0; c < 4; ++c) {
#pragma unroll
                            for (int d = 0; d < 3; ++d) {
                                fmin_[d] = fmin(fmin_[d], fmin(plane_a[c][d], plane_b[c][d]));
                                fmax_[d] = fmax(fmax_[d], fmax(plane_a[c][d], plane_b[c][d]));
                            }
                        }

                        bool accepted = false;
                        const bool may_contain =
                            codomain_acceptance(fmin_, fmax_, tol, tols, err, accepted);

                        // Carry this sub-box's upper face forward as the next
                        // sub-box's lower face.
#pragma unroll
                        for (int c = 0; c < 4; ++c) {
#pragma unroll
                            for (int d = 0; d < 3; ++d) plane_a[c][d] = plane_b[c][d];
                        }

                        if (!may_contain) continue;

                        // TightInclusion's no_zero_toi policy, as on the host: a
                        // contact reported at exactly t == 0 stalls a solver.
                        accepted = accepted && (tt_min > TC(0));

                        Domain sub = box;
                        if (dim == 0) { sub.tlo = smin; sub.thi = smax; }
                        else if (dim == 1) { sub.ulo = smin; sub.uhi = smax; }
                        else { sub.vlo = smin; sub.vhi = smax; }
                        sub.depth = box.depth + 1;

                        // Exhaustion accepts. Dropping a box that may contain a
                        // root is the one way this algorithm loses a collision.
                        //
                        // The stack test is unreachable given the depth clamp above
                        // -- that is the point of tying the two together -- and it
                        // stays as the guard that makes the invariant hold even if
                        // someone changes kSplits without redoing the arithmetic.
                        const bool exhausted = (sub.depth > depth_limit) || (top >= kStackCap);
                        if (accepted || exhausted) {
                            if (sub.tlo < t) {
                                t = sub.tlo;
                                found = true;
                            }
                            continue;
                        }

                        stack[top++] = sub;
                    }
                }

                if (found) {
                    if (toi_stride == 0) {
                        const double previous = atomic_min_double(shared_toi, (double)t);
                        (void)previous;
                    } else {
                        // Narrow toward negative infinity: an earlier time of impact
                        // is safe, a later one is the failure this exists to prevent.
                        toi[i] = sizeof(T) == sizeof(float) ? (T)__double2float_rd(t) : (T)t;
                    }
                } else if (toi_stride == 1) {
                    toi[i] = sizeof(T) == sizeof(float) ? (T)__double2float_rd((double)max_toi)
                                                        : (T)max_toi;
                }
            }

            __global__ void write_shared_toi_kernel(const double* const shared_toi,
                                                    float* const out) {
                if (threadIdx.x == 0 && blockIdx.x == 0) *out = __double2float_rd(*shared_toi);
            }

            __global__ void write_shared_toi_kernel_d(const double* const shared_toi,
                                                      double* const out) {
                if (threadIdx.x == 0 && blockIdx.x == 0) *out = *shared_toi;
            }

        }  // namespace

        template <typename T, typename I>
        int narrow_phase_vq(const size_t noverlaps,
                            const I* const SCCD_RESTRICT voverlap,
                            const I* const SCCD_RESTRICT qoverlap,
                            T** const SCCD_RESTRICT v0,
                            T** const SCCD_RESTRICT v1,
                            const size_t quad_stride,
                            I** const SCCD_RESTRICT quads,
                            const T max_toi,
                            T* const SCCD_RESTRICT toi,
                            const int max_depth,
                            const T tol,
                            const int toi_stride) {
            if (noverlaps == 0) {
                if (toi != nullptr && toi_stride == 0) {
                    SCCD_CHECK_CUDA(cudaMemcpy(toi, &max_toi, sizeof(T), cudaMemcpyHostToDevice));
                }
                return 0;
            }

            double* d_shared = nullptr;
            SCCD_CHECK_CUDA(cudaMalloc(&d_shared, sizeof(double)));
            const double h_shared = (double)max_toi;
            SCCD_CHECK_CUDA(cudaMemcpy(d_shared, &h_shared, sizeof(double), cudaMemcpyHostToDevice));

            // A depth the device cannot reach must not pass unremarked. The
            // clamp itself is safe -- exhaustion accepts at the box's t lower
            // bound, so a shallower search reports an earlier time of impact and
            // never a later one -- but it is an accuracy divergence from the
            // host, which honours whatever depth it is given, and a silent one
            // is the kind nobody finds. SCCD_VQ_MAX_DEPTH sizes the stack, and
            // the measurements above say raising it is free.
            if (max_depth > kMaxDeviceDepth) {
                static bool warned = false;
                if (!warned) {
                    warned = true;
                    std::fprintf(stderr,
                                 "SCCD: device vertex-quad narrow phase clamps max_depth %d to %d "
                                 "(the compiled stack holds %d boxes). The reported time of impact "
                                 "will be earlier than the host's, never later. Rebuild with "
                                 "-DSCCD_VQ_MAX_DEPTH=%d to match.\n",
                                 max_depth, kMaxDeviceDepth, kStackCap, max_depth);
                }
            }

            const int block = 128;
            const int grid = (int)((noverlaps + block - 1) / block);

            // v0, v1 and quads are *device* arrays of device pointers -- the same
            // convention narrow_phase_vf and narrow_phase_ee take, and what the
            // CCD integration hands all three. The kernel wants the rows as
            // separate arguments, so the pointers have to be brought back to the
            // host first: indexing v0[0] here would dereference device memory on
            // the host, which faults. It had gone unnoticed because no test
            // reached this entry point and the quad device path is new.
            T* h_v0[3];
            T* h_v1[3];
            I* h_quads[4];
            SCCD_CHECK_CUDA(cudaMemcpy(h_v0, v0, sizeof(T*) * 3, cudaMemcpyDeviceToHost));
            SCCD_CHECK_CUDA(cudaMemcpy(h_v1, v1, sizeof(T*) * 3, cudaMemcpyDeviceToHost));
            SCCD_CHECK_CUDA(cudaMemcpy(h_quads, quads, sizeof(I*) * 4, cudaMemcpyDeviceToHost));

            narrow_phase_vq_kernel<T, I><<<grid, block>>>(noverlaps,
                                                          voverlap,
                                                          qoverlap,
                                                          h_v0[0], h_v0[1], h_v0[2],
                                                          h_v1[0], h_v1[1], h_v1[2],
                                                          quad_stride,
                                                          h_quads[0], h_quads[1], h_quads[2], h_quads[3],
                                                          (TC)max_toi,
                                                          d_shared,
                                                          toi,
                                                          max_depth,
                                                          (TC)tol,
                                                          toi_stride);
            SCCD_CUDA_LAST_ERROR();

            if (toi_stride == 0) {
                if constexpr (sizeof(T) == sizeof(float)) {
                    write_shared_toi_kernel<<<1, 1>>>(d_shared, (float*)toi);
                } else {
                    write_shared_toi_kernel_d<<<1, 1>>>(d_shared, (double*)toi);
                }
                SCCD_CUDA_LAST_ERROR();
            }

            SCCD_CHECK_CUDA(cudaDeviceSynchronize());
            SCCD_CHECK_CUDA(cudaFree(d_shared));
            return 0;
        }

    }  // namespace device
}  // namespace sccd

#define SCCD_VQ_INSTANTIATE(T, I)                                                        \
    template int sccd::device::narrow_phase_vq<T, I>(const size_t,                       \
                                                     const I* const,                     \
                                                     const I* const,                     \
                                                     T** const,                          \
                                                     T** const,                          \
                                                     const size_t,                       \
                                                     I** const,                          \
                                                     const T,                            \
                                                     T* const,                           \
                                                     const int,                           \
                                                     const T,                            \
                                                     const int);

SCCD_VQ_INSTANTIATE(float, int32_t)
SCCD_VQ_INSTANTIATE(double, int32_t)

#undef SCCD_VQ_INSTANTIATE
