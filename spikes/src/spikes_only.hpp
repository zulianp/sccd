#ifndef SCCD_SPIKES_ONLY_HPP
#define SCCD_SPIKES_ONLY_HPP

// Shipped-header code whose only callers live in spikes/.
//
// Distinct from dead.hpp: these are called, just never by anything that ships.
// Leaving them in src/ meant the installed public headers carried code the
// library never reaches and a consumer compiles for nothing.
//
// spikes/README.md applies: not built by default, not installed, deletable
// alongside the spike that uses them.

#include "sccd_aabb.hpp"
#include "sccd_base.hpp"
#include "sccd_broadphase_sweep.hpp"
#include "sccd_parallel.hpp"

#include <cstddef>
#include <cstdint>

namespace sccd {
    namespace detail {

        // ---- build_disjoint_mask_for_block
        // from sccd_broadphase_sweep.hpp; callers: spikes/src/broadphase_lb.hpp, cell_broadphase.hpp

        /**
         * \brief Build per-lane disjoint mask for a block of B against a single A.
         * \param second_aabbs SoA arrays for the second list.
         * \param start First B index to test.
         * \param chunk_len Number of lanes to process (<=
         * SCCD_AABB_DISJOINT_CHUNK_SIZE).
         * \param aminx Scalar A min x.
         * \param aminy Scalar A min y.
         * \param aminz Scalar A min z.
         * \param amaxx Scalar A max x.
         * \param amaxy Scalar A max y.
         * \param amaxz Scalar A max z (also used for tail fill).
         * \param mask_out Output mask: 1=disjoint, 0=potential overlap.
         */
        template <typename T>
        static inline void build_disjoint_mask_for_block(T **const SCCD_RESTRICT second_aabbs,
                                                         const ptrdiff_t start,
                                                         const ptrdiff_t chunk_len,
                                                         const T aminx,
                                                         const T aminy,
                                                         const T aminz,
                                                         const T amaxx,
                                                         const T amaxy,
                                                         const T amaxz,
                                                         uint32_t *const SCCD_RESTRICT mask_out) {
            if (chunk_len != SCCD_AABB_DISJOINT_CHUNK_SIZE) {
                alignas(64) T B_minx[SCCD_AABB_DISJOINT_CHUNK_SIZE];
                alignas(64) T B_miny[SCCD_AABB_DISJOINT_CHUNK_SIZE];
                alignas(64) T B_minz[SCCD_AABB_DISJOINT_CHUNK_SIZE];
                alignas(64) T B_maxx[SCCD_AABB_DISJOINT_CHUNK_SIZE];
                alignas(64) T B_maxy[SCCD_AABB_DISJOINT_CHUNK_SIZE];
                alignas(64) T B_maxz[SCCD_AABB_DISJOINT_CHUNK_SIZE];

                prepare_B_block(second_aabbs, start, chunk_len, B_minx, B_miny, B_minz, B_maxx, B_maxy, B_maxz);
                tail_fill_B(amaxx, amaxy, amaxz, chunk_len, B_minx, B_miny, B_minz, B_maxx, B_maxy, B_maxz);
                vaabb_disjoint_one_to_many(
                    aminx, aminy, aminz, amaxx, amaxy, amaxz, B_minx, B_miny, B_minz, B_maxx, B_maxy, B_maxz, mask_out);
            } else {
                vaabb_disjoint_one_to_many(aminx,
                                           aminy,
                                           aminz,
                                           amaxx,
                                           amaxy,
                                           amaxz,
                                           second_aabbs[0] + start,
                                           second_aabbs[1] + start,
                                           second_aabbs[2] + start,
                                           second_aabbs[3] + start,
                                           second_aabbs[4] + start,
                                           second_aabbs[5] + start,
                                           mask_out);
            }
        }

        // ---- mask_out_shared_two_lists
        // from sccd_broadphase_sweep.hpp; callers: spikes/src/broadphase_lb.hpp, cell_broadphase.hpp

        /**
         * \brief Mark lanes where A and B share a vertex (invalid pairs).
         * \tparam F Number of vertices in first element.
         * \tparam S Number of vertices in second element.
         * \param dmask In/out lane mask; set to 1 when a vertex is shared.
         * \param chunk_len Number of valid lanes.
         * \param noffset Starting index in B for this chunk.
         * \param ev Vertex indices of the current A element.
         * \param second_idx Index mapping for second list elements.
         * \param second_elements SoA vertex arrays for second list.
         * \param second_element_stride Stride between elements in second arrays.
         */
        template <int F, int S, typename I>
        static inline void mask_out_shared_two_lists(uint32_t *const SCCD_RESTRICT dmask,
                                                     const ptrdiff_t chunk_len,
                                                     const ptrdiff_t noffset,
                                                     const I (&ev)[F],
                                                     const I *const SCCD_RESTRICT second_idx,
                                                     I **const SCCD_RESTRICT second_elements,
                                                     const ptrdiff_t second_element_stride) {
            for (ptrdiff_t lane = 0; lane < chunk_len; ++lane) {
                if (dmask[lane]) continue;

                const ptrdiff_t j = noffset + lane;
                const I jidx = second_idx[j];
                int match = 0;
                if constexpr (S > 1) {
                    I sev[S];
                    for (int v = 0; v < S; ++v) {
                        sev[v] = second_elements[v][jidx * second_element_stride];
                    }
                    for (int a = 0; a < F; ++a) {
                        for (int b = 0; b < S; ++b) {
                            match |= (ev[a] == sev[b]);
                        }
                    }
                } else {
                    for (int a = 0; a < F; ++a) {
                        match |= (ev[a] == jidx);
                    }
                }
                dmask[lane] |= (uint32_t)match;
            }
        }

        // ---- parallel_cum_max_br
        // from sccd_parallel.hpp; caller: spikes/hacks/conversions.hpp

    template <typename T>
    void parallel_cum_max_br(T* const begin, T* const end) {
        const ptrdiff_t len = end - begin;
        if (len <= 0) {
            return;
        }

#ifdef SCCD_ENABLE_TBB
        tbb::parallel_scan(
            tbb::blocked_range<ptrdiff_t>(0, len),
            begin[0],
            [=](const tbb::blocked_range<ptrdiff_t>& r, T acc, bool is_final_scan) -> T {
                if (!is_final_scan) {
                    T temp = acc;
                    for (ptrdiff_t i = r.begin(); i < r.end(); ++i) {
                        temp = sccd::max(temp, begin[i]);
                    }
                    return temp;
                } else {
                    begin[r.begin()] = sccd::max(begin[r.begin()], acc);
                    for (ptrdiff_t i = r.begin() + 1; i < r.end(); ++i) {
                        begin[i] = sccd::max(begin[i], begin[i - 1]);
                    }

                    return begin[r.end() - 1];
                }
            },
            [](T left, T right) { return sccd::max(left, right); });
#else
        detail::parallel_inclusive_scan<T>(begin, len, [](const T a, const T b) { return sccd::max(a, b); });
#endif  // SCCD_ENABLE_TBB
    }

    }  // namespace detail

    namespace detail {
        // ---- from src/broadphase/sccd_broadphase_sweep.hpp ----
        //
        // Four helpers of the sorted-interval sweep whose only callers are the
        // spikes in this directory: broadphase_lb.hpp, cell_broadphase.hpp and
        // spikes_only.hpp itself. The shipped sweep does not use any of them.
        // They stayed in an installed public header, so a consumer compiled
        // them for nothing.

        template <typename T>
        static inline void prepare_B_block(T **const SCCD_RESTRICT aabbs,
                                           const ptrdiff_t start,
                                           const ptrdiff_t len,
                                           T *const SCCD_RESTRICT B_minx,
                                           T *const SCCD_RESTRICT B_miny,
                                           T *const SCCD_RESTRICT B_minz,
                                           T *const SCCD_RESTRICT B_maxx,
                                           T *const SCCD_RESTRICT B_maxy,
                                           T *const SCCD_RESTRICT B_maxz) {
            for (ptrdiff_t lane = 0; lane < len; ++lane) {
                const ptrdiff_t j = start + lane;
                B_minx[lane] = aabbs[0][j];
                B_miny[lane] = aabbs[1][j];
                B_minz[lane] = aabbs[2][j];
                B_maxx[lane] = aabbs[3][j];
                B_maxy[lane] = aabbs[4][j];
                B_maxz[lane] = aabbs[5][j];
            }
        }

        template <typename T>
        static inline void tail_fill_B(const T amaxx0,
                                       const T amaxy0,
                                       const T amaxz0,
                                       const ptrdiff_t len,
                                       T *const SCCD_RESTRICT B_minx,
                                       T *const SCCD_RESTRICT B_miny,
                                       T *const SCCD_RESTRICT B_minz,
                                       T *const SCCD_RESTRICT B_maxx,
                                       T *const SCCD_RESTRICT B_maxy,
                                       T *const SCCD_RESTRICT B_maxz) {
            for (ptrdiff_t lane = len; lane < SCCD_AABB_DISJOINT_CHUNK_SIZE; ++lane) {
                B_minx[lane] = amaxx0 + 1;
                B_miny[lane] = amaxy0 + 1;
                B_minz[lane] = amaxz0 + 1;
                B_maxx[lane] = amaxx0;
                B_maxy[lane] = amaxy0;
                B_maxz[lane] = amaxz0;
            }
        }

        template <int F, int S, typename T, typename I>
        static inline ptrdiff_t scalar_count_range_two_lists(T **const SCCD_RESTRICT first_aabbs,
                                                             const ptrdiff_t fi,
                                                             T **const SCCD_RESTRICT second_aabbs,
                                                             const I *const SCCD_RESTRICT second_idx,
                                                             I **const SCCD_RESTRICT second_elements,
                                                             const ptrdiff_t second_element_stride,
                                                             const I (&ev)[F],
                                                             const ptrdiff_t begin,
                                                             const ptrdiff_t end) {
            ptrdiff_t count = 0;
            const T aminx = first_aabbs[0][fi];
            const T aminy = first_aabbs[1][fi];
            const T aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi];
            const T amaxy = first_aabbs[4][fi];
            const T amaxz = first_aabbs[5][fi];
            for (ptrdiff_t j = begin; j < end; ++j) {
                if (disjoint(aminx,
                             aminy,
                             aminz,
                             amaxx,
                             amaxy,
                             amaxz,
                             second_aabbs[0][j],
                             second_aabbs[1][j],
                             second_aabbs[2][j],
                             second_aabbs[3][j],
                             second_aabbs[4][j],
                             second_aabbs[5][j])) {
                    continue;
                }
                bool share = false;
                if constexpr (S > 1) {
                    const I jidx = second_idx[j];
                    I sev[S];
                    for (int v = 0; v < S; ++v) {
                        sev[v] = second_elements[v][jidx * second_element_stride];
                    }
                    share = shares_vertex<F, S>(ev, sev);
                } else {
                    for (int a = 0; a < F; ++a) {
                        if (ev[a] == second_idx[j]) {
                            share = true;
                            break;
                        }
                    }
                }
                count += share ? 0 : 1;
            }
            return count;
        }

        template <int F, int S, typename T, typename I>
        static inline ptrdiff_t scalar_collect_range_two_lists(T **const SCCD_RESTRICT first_aabbs,
                                                               const ptrdiff_t fi,
                                                               const I first_idxi,
                                                               T **const SCCD_RESTRICT second_aabbs,
                                                               const I *const SCCD_RESTRICT second_idx,
                                                               I **const SCCD_RESTRICT second_elements,
                                                               const ptrdiff_t second_element_stride,
                                                               const I (&ev)[F],
                                                               const ptrdiff_t begin,
                                                               const ptrdiff_t end,
                                                               I *const SCCD_RESTRICT first_out,
                                                               I *const SCCD_RESTRICT second_out) {
            ptrdiff_t count = 0;
            const T aminx = first_aabbs[0][fi];
            const T aminy = first_aabbs[1][fi];
            const T aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi];
            const T amaxy = first_aabbs[4][fi];
            const T amaxz = first_aabbs[5][fi];
            for (ptrdiff_t j = begin; j < end; ++j) {
                if (disjoint(aminx,
                             aminy,
                             aminz,
                             amaxx,
                             amaxy,
                             amaxz,
                             second_aabbs[0][j],
                             second_aabbs[1][j],
                             second_aabbs[2][j],
                             second_aabbs[3][j],
                             second_aabbs[4][j],
                             second_aabbs[5][j])) {
                    continue;
                }
                const I jidx = second_idx[j];
                int match = 0;
                if constexpr (S > 1) {
                    I sev[S];
                    for (int v = 0; v < S; ++v) {
                        sev[v] = second_elements[v][jidx * second_element_stride];
                    }
                    for (int a = 0; a < F; ++a) {
                        for (int b = 0; b < S; ++b) {
                            match |= (ev[a] == sev[b]);
                        }
                    }
                } else {
                    for (int a = 0; a < F; ++a) {
                        match |= (ev[a] == jidx);
                    }
                }
                if (!match) {
                    first_out[count] = first_idxi;
                    second_out[count] = jidx;
                    count += 1;
                }
            }
            return count;
        }

    }  // namespace detail

}  // namespace sccd

#endif  // SCCD_SPIKES_ONLY_HPP
