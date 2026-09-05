#ifndef SCCD_NARROWPHASE_QUAD_HPP
#define SCCD_NARROWPHASE_QUAD_HPP

/**
 * \file
 * \brief Vertex-quad narrow phase.
 *
 * ## Node order: lexicographic, not cyclic
 *
 * The four quad nodes handed to every entry point here -- `quads[0..3][qi]` --
 * are in **lexicographic** order over the parameter square, not wound around the
 * quad. The bilinear form weights node k by
 *
 *     w1 = (1-u)(1-v)   w2 = u(1-v)   w3 = (1-u)v   w4 = uv
 *
 * so the nodes are the corners `(0,0)`, `(1,0)`, `(0,1)`, `(1,1)` in that order.
 * A cyclic winding swaps the last two.
 *
 * This is worth stating loudly because getting it wrong is silent. The solver
 * does not validate the quad; it simply searches a bilinear patch through the
 * nodes in the order given, which for a cyclic winding is a saddle stretched
 * across the diagonal rather than the quad the caller meant. Contacts against
 * the real quad are then not found, and nothing reports an error -- the answer
 * is a correct time of impact for the wrong surface.
 */

#include "sccd_base.hpp"
#include "sccd_parallel.hpp"
#include "sccd_narrowphase_mode.hpp"
#include "sccd_rootfinder_quad.hpp"

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <vector>

namespace sccd {

    template <typename T, typename I>
    int narrow_phase_vq(const size_t noverlaps,
                        const I *const SCCD_RESTRICT voverlap,
                        const I *const SCCD_RESTRICT qoverlap,
                        T **const SCCD_RESTRICT v0,
                        T **const SCCD_RESTRICT v1,
                        const size_t element_stride,
                        I **const SCCD_RESTRICT quads,
                        const T max_toi,
                        T *const SCCD_RESTRICT toi,
                        const int max_depth,
                        const T tol,
                        const ToiOutput toi_output = ToiOutput::Earliest) {
        using T_HP = double;

        narrow_phase_mode_note_quads_ignore();
        if (noverlaps == 0) {
            if (toi != nullptr && toi_output == ToiOutput::Earliest) {
                toi[0] = max_toi;
            }
            return 0;
        }
        assert(toi != nullptr);


        std::atomic<T> min_t = max_toi;

        if (toi_output == ToiOutput::PerPair) {
            sccd::parallel_for_br_dynamic(0, (ptrdiff_t)noverlaps,
                                          [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
                                              for (ptrdiff_t i = rbegin; i < rend; ++i) toi[i] = max_toi;
                                          });
        }

        if (toi_output == ToiOutput::Earliest && min_t.load(std::memory_order_relaxed) == T(0)) {
            toi[0] = T(0);
            return 0;
        }

        if (toi_output == ToiOutput::Earliest) {
            toi[0] = max_toi;
        }

        // Per-candidate cost varies by orders of magnitude, hence the dynamic
        // schedule. Going through parallel_for_br_dynamic rather than a raw
        // OpenMP pragma is what lets a TBB build use TBB: with SCCD_ENABLE_TBB
        // set and no -fopenmp the pragmas below were ignored and the whole quad
        // narrow phase ran on one thread, and with both it oversubscribed. The
        // stack is thread_local rather than block-local so it keeps its capacity
        // across blocks instead of reallocating per chunk.
        sccd::parallel_for_br_dynamic(0, (ptrdiff_t)noverlaps, [&](const ptrdiff_t rbegin,
                                                                   const ptrdiff_t rend) {
            static thread_local std::vector<Box<T_HP>> stack;
            stack.reserve(64);

            for (ptrdiff_t i = rbegin; i < rend; ++i) {
                if (toi_output == ToiOutput::PerPair) {
                    toi[i] = max_toi;
                }

                const I vi = voverlap[i];
                const I qi = qoverlap[i];

                const I nodes[4] = {quads[0][qi * element_stride],
                                    quads[1][qi * element_stride],
                                    quads[2][qi * element_stride],
                                    quads[3][qi * element_stride]};

                const T_HP sv[3] = {v0[0][vi], v0[1][vi], v0[2][vi]};
                const T_HP ev[3] = {v1[0][vi], v1[1][vi], v1[2][vi]};

                const T_HP s1[3] = {v0[0][nodes[0]], v0[1][nodes[0]], v0[2][nodes[0]]};
                const T_HP s2[3] = {v0[0][nodes[1]], v0[1][nodes[1]], v0[2][nodes[1]]};
                const T_HP s3[3] = {v0[0][nodes[2]], v0[1][nodes[2]], v0[2][nodes[2]]};
                const T_HP s4[3] = {v0[0][nodes[3]], v0[1][nodes[3]], v0[2][nodes[3]]};

                const T_HP e1[3] = {v1[0][nodes[0]], v1[1][nodes[0]], v1[2][nodes[0]]};
                const T_HP e2[3] = {v1[0][nodes[1]], v1[1][nodes[1]], v1[2][nodes[1]]};
                const T_HP e3[3] = {v1[0][nodes[2]], v1[1][nodes[2]], v1[2][nodes[2]]};
                const T_HP e4[3] = {v1[0][nodes[3]], v1[1][nodes[3]], v1[2][nodes[3]]};

                T_HP t = toi_output == ToiOutput::Earliest ? T_HP(min_t.load(std::memory_order_relaxed)) : T_HP(toi[i]);
                T_HP u = T_HP(0);
                T_HP v = T_HP(0);

                const T_HP t_upper = sccd::min<T_HP>(t, T_HP(1));
                const Box<T_HP> initial_domain(Interval<T_HP>{T_HP(0), t_upper},
                                               Interval<T_HP>{T_HP(0), T_HP(1)},
                                               Interval<T_HP>{T_HP(0), T_HP(1)},
                                               0);

#if !SCCD_ENABLE_CODOMAIN_SCALING
                T_HP codomain_widths[3] = {T_HP(1), T_HP(1), T_HP(1)};
                T_HP tols[3];
                compute_vertex_quad_tolerance<T_HP>(tol, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, tols);
                T_HP numerical_error[3];
                vq_numerical_error<T_HP>(sv, s1, s2, s3, s4, ev, e1, e2, e3, e4, numerical_error);
#else
                // One pass for the tolerances, the codomain widths and the
                // certified error bound. With a shared time of impact the search
                // explores about 1.2 boxes per query, so this preamble is most of
                // the query rather than something amortised over a deep search.
                T_HP codomain_widths[3];
                T_HP tols[3];
                T_HP numerical_error[3];
                // t_upper, not 1: the split scale should follow the window the
                // search is actually confined to, which the shared time of
                // impact has already pruned. Same as the triangle path.
                vq_prepare<T_HP>(tol, t_upper, sv, s1, s2, s3, s4, ev, e1, e2, e3, e4,
                                 tols, codomain_widths, numerical_error);
#endif

                stack.clear();
                stack.push_back(initial_domain);

                while (!stack.empty()) {
                    Box<T_HP> box = stack.back();
                    stack.pop_back();

                    if (box.tuv[0].lower >= t) {
                        continue;
                    }

                    box.tuv[0].upper = sccd::min<T_HP>(box.tuv[0].upper, t);

                    const bool found = find_root_grid_adaptive_split_vq<T_HP>(max_depth,
                                                                              tol,
                                                                              tols,
                                                                              numerical_error,
                                                                              codomain_widths,
                                                                              sv,
                                                                              s1,
                                                                              s2,
                                                                              s3,
                                                                              s4,
                                                                              ev,
                                                                              e1,
                                                                              e2,
                                                                              e3,
                                                                              e4,
                                                                              box,
                                                                              t,
                                                                              u,
                                                                              v,
                                                                              stack,
                                                                              /*refine=*/false);

                    if (found) {
                        if (toi_output == ToiOutput::Earliest) {
                            const T previous = sccd::atomic_min<T>(min_t, T(t));
                            if (previous < T(t)) {
                                t = T_HP(previous);
                            }
                        } else {
                            toi[i] = T(t);
                        }
                    } else if (!stack.empty() && toi_output == ToiOutput::Earliest) {
                        t = sccd::min<T_HP>(t, T_HP(min_t.load(std::memory_order_relaxed)));
                    }
                }
            }
        });

        if (toi_output == ToiOutput::Earliest) {
            toi[0] = min_t.load(std::memory_order_relaxed);
        }
        return 0;
    }


} // namespace sccd

#endif
