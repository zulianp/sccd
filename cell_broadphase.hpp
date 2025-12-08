#ifndef CELL_BROADPHASE_HPP
#define CELL_BROADPHASE_HPP

#include "broadphase.hpp"
#include "sparallel.hpp"

namespace sccd {
    namespace sccd_detail {

        template <int F, int S, typename T, typename I>
        static inline void cell_mask_out_shared_two_lists(uint32_t *const SCCD_RESTRICT dmask,
                                                          const ptrdiff_t chunk_len,
                                                          const ptrdiff_t noffset,
                                                          const I (&ev)[F],
                                                          const I *const SCCD_RESTRICT second_idx,
                                                          I **const SCCD_RESTRICT second_elements,
                                                          const ptrdiff_t second_stride,
                                                          const I *const SCCD_RESTRICT cellidx) {
            for (ptrdiff_t lane = 0; lane < chunk_len; ++lane) {
                if (dmask[lane]) continue;

                const ptrdiff_t j = cellidx[noffset + lane];
                const I jidx = second_idx[j];
                int match = 0;
                if constexpr (S > 1) {
                    I sev[S];
                    for (int v = 0; v < S; ++v) {
                        sev[v] = second_elements[v][jidx * second_stride];
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

        template <int F, int S, typename T, typename I>
        static inline ptrdiff_t cell_scalar_count_range_two_lists(T **const SCCD_RESTRICT first_aabbs,
                                                                  const ptrdiff_t fi,
                                                                  T **const SCCD_RESTRICT second_aabbs,
                                                                  const I *const SCCD_RESTRICT second_idx,
                                                                  I **const SCCD_RESTRICT second_elements,
                                                                  const ptrdiff_t second_stride,
                                                                  const I (&ev)[F],
                                                                  const ptrdiff_t begin_k,
                                                                  const ptrdiff_t end_k,
                                                                  const I *const SCCD_RESTRICT cellidx) {
            ptrdiff_t count = 0;
            const T aminx = first_aabbs[0][fi];
            const T aminy = first_aabbs[1][fi];
            const T aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi];
            const T amaxy = first_aabbs[4][fi];
            const T amaxz = first_aabbs[5][fi];
            for (ptrdiff_t k = begin_k; k < end_k; ++k) {
                const ptrdiff_t j = cellidx[k];
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
                        sev[v] = second_elements[v][jidx * second_stride];
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

        /**
         * \brief Scalar reference: collect candidate overlaps in [begin,end) for
         * two lists. \return Number of pairs written to \p first_out and \p
         * second_out.
         */
        template <int F, int S, typename T, typename I>
        static inline ptrdiff_t cell_scalar_collect_range_two_lists(T **const SCCD_RESTRICT first_aabbs,      // 0
                                                                    const ptrdiff_t fi,                       // 1
                                                                    const I first_idxi,                       // 2
                                                                    T **const SCCD_RESTRICT second_aabbs,     // 3
                                                                    const I *const SCCD_RESTRICT second_idx,  // 4
                                                                    I **const SCCD_RESTRICT second_elements,  // 5
                                                                    const ptrdiff_t second_stride,            // 6
                                                                    const I (&ev)[F],                         // 7
                                                                    const ptrdiff_t begin_k,                  // 8
                                                                    const ptrdiff_t end_k,                    // 9
                                                                    const I *const SCCD_RESTRICT cellidx,     // 10
                                                                    I *const SCCD_RESTRICT first_out,         // 11
                                                                    I *const SCCD_RESTRICT second_out) {      // 12
            ptrdiff_t count = 0;
            const T aminx = first_aabbs[0][fi];
            const T aminy = first_aabbs[1][fi];
            const T aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi];
            const T amaxy = first_aabbs[4][fi];
            const T amaxz = first_aabbs[5][fi];
            for (ptrdiff_t k = begin_k; k < end_k; ++k) {
                const ptrdiff_t j = cellidx[k];
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
                        sev[v] = second_elements[v][jidx * second_stride];
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

    }  // namespace sccd_detail

    template <typename T>
    static void cell_setup(const ptrdiff_t n,
                           const T *const SCCD_RESTRICT xmin,
                           const T *const SCCD_RESTRICT xmax,
                           ptrdiff_t *const SCCD_RESTRICT inout_ncells,
                           T *const SCCD_RESTRICT cell_min,
                           T *const SCCD_RESTRICT cell_size) {
        T gmin = xmin[0];
        T gmax = xmax[0];
        T max_cell_size = xmax[0] - xmin[0];
        for (ptrdiff_t i = 0; i < n; i++) {
            T xmini = xmin[i] - 1e-6;
            T xmaxi = xmax[i] + 1e-6;
            gmin = std::min(gmin, xmini);
            gmax = std::max(gmax, xmaxi);
            max_cell_size = std::max(max_cell_size, xmaxi - xmini);
        }

        *inout_ncells =
            std::min(*inout_ncells, std::max((ptrdiff_t)1, (ptrdiff_t)floor((gmax - gmin) / max_cell_size)));
        *cell_min = gmin;
        *cell_size = (gmax - gmin) / *inout_ncells;  // Recompute due to floating point perturbation
    }

    template <typename T>
    static void cell_starts_setup(const ptrdiff_t n,
                                  const T *const SCCD_RESTRICT xmin,
                                  const T *const SCCD_RESTRICT xmax,
                                  const ptrdiff_t ncells,
                                  T *const SCCD_RESTRICT cell_min,
                                  T *const SCCD_RESTRICT cell_size) {
        T gmin = xmin[0];
        T gmax = xmax[0];
        for (ptrdiff_t i = 0; i < n; i++) {
            T xmini = xmin[i] - T(1e-6);
            T xmaxi = xmax[i] + T(1e-6);
            gmin = std::min(gmin, xmini);
            gmax = std::max(gmax, xmaxi);
        }
        *cell_min = gmin;
        *cell_size = (gmax - gmin) / ncells;
    }

    template <typename T>
    static void cell_starts(const ptrdiff_t ncells,
                            const T cell_min,
                            const T cell_size,
                            const ptrdiff_t n,
                            const T *const SCCD_RESTRICT xmin,
                            ptrdiff_t *const SCCD_RESTRICT starts) {
        starts[0] = 0;
        T current_xmin = xmin[0];
        ptrdiff_t current_start = 0;

        for (ptrdiff_t i = 1; i < n; i++) {
            const T x_i = xmin[i];
            const int cell_idx = (x_i - cell_min) / cell_size;
            assert(cell_idx >= 0 && cell_idx < ncells);
            assert(cell_idx >= current_start);
            if (cell_idx != current_start) {
                starts[cell_idx] = i;
                current_start = cell_idx;
            }
        }
        // Backfill empty cells to make starts non-decreasing
        for (ptrdiff_t c = 1; c < ncells; ++c) {
            if (starts[c] == 0) starts[c] = starts[c - 1];
        }
    }
    template <typename T>
    static void cell_count(
        // Cell list
        const ptrdiff_t ncells,
        const T cell_min,
        const T cell_size,
        // Points
        const ptrdiff_t n,
        const T *const SCCD_RESTRICT xmin,
        // Output
        int *const SCCD_RESTRICT cellptr) {
        memset(cellptr, 0, sizeof(int) * (ncells + 1));

        for (ptrdiff_t i = 0; i < n; i++) {
            const T x_i = xmin[i];
            const int cell_idx = nextafter_down((x_i - cell_min) / cell_size);
            assert(cell_idx >= 0 && cell_idx < ncells);
            cellptr[cell_idx + 1]++;
        }

        for (ptrdiff_t i = 0; i < ncells; i++) {
            cellptr[i + 1] += cellptr[i];
        }
    }

    template <typename T, typename I>
    static void cell_populate(
        // Cell list
        const ptrdiff_t ncells,
        const T cell_min,
        const T cell_size,
        // Points
        const ptrdiff_t n,
        const T *const SCCD_RESTRICT xmin,
        const int *const SCCD_RESTRICT cellptr,
        I *const SCCD_RESTRICT idx,
        I *const SCCD_RESTRICT bookkeeping) {
        memset(bookkeeping, 0, sizeof(I) * ncells);
        for (ptrdiff_t i = 0; i < n; i++) {
            const T x_i = xmin[i];
            const int cell_idx = nextafter_down((x_i - cell_min) / cell_size);
            assert(cell_idx >= 0 && cell_idx < ncells);
            idx[cellptr[cell_idx] + bookkeeping[cell_idx]++] = i;
        }

#ifndef NDEBUG
        for (ptrdiff_t i = 0; i < ncells; i++) {
            assert(cellptr[i + 1] - cellptr[i] == bookkeeping[i]);

            I prev = idx[cellptr[i]];
            for (ptrdiff_t j = cellptr[i] + 1; j < cellptr[i + 1]; j++) {
                assert(idx[j] > prev);
                prev = idx[j];
            }
        }
#endif
    }

    template <int first_nxe, int second_nxe, typename T, typename I>
    bool cell_count_overlaps(const int sort_axis,
                             const count_t first_count,
                             T **const SCCD_RESTRICT first_aabbs,
                             I *const SCCD_RESTRICT first_idx,
                             const ptrdiff_t first_stride,
                             I **const SCCD_RESTRICT first_elements,
                             const count_t second_count,
                             T **const SCCD_RESTRICT second_aabbs,
                             I *const SCCD_RESTRICT second_idx,
                             const ptrdiff_t second_stride,
                             I **const SCCD_RESTRICT second_elements,
                             // Cell list
                             const int cell_list_axis,
                             const ptrdiff_t ncells,
                             const T cell_min,
                             const T cell_size,
                             const I *const SCCD_RESTRICT cellptr,
                             const I *const SCCD_RESTRICT cellidx,
                             ptrdiff_t *const SCCD_RESTRICT ccdptr) {
        const T *const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
        const T *const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
        const T *const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
        const T *const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

        const T *const SCCD_RESTRICT cell_xmin = first_aabbs[cell_list_axis];
        const T *const SCCD_RESTRICT cell_xmax = first_aabbs[3 + cell_list_axis];

        ccdptr[0] = 0;
        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; fi++) {
                const T fimin = first_xmin[fi];
                const T fimax = first_xmax[fi];
                const I first_idxi = first_idx[fi];

                ptrdiff_t cell_start =
                    std::max((int)0, (int)(floor(nextafter_down((cell_xmin[fi] - cell_min) / cell_size) - 1)));
                ptrdiff_t cell_end =
                    std::min((int)ncells, (int)(floor(nextafter_up((cell_xmax[fi] - cell_min) / cell_size))) + 1);

                I ev[first_nxe];
                for (int v = 0; v < first_nxe; v++) {
                    ev[v] = first_elements[v][first_idxi * first_stride];
                }

                ptrdiff_t count = 0;
                for (ptrdiff_t cell_i = cell_start; cell_i < cell_end; cell_i++) {
                    ptrdiff_t cell_i_begin = cellptr[cell_i];
                    ptrdiff_t cell_i_end = cellptr[cell_i + 1];

                    ptrdiff_t begin_k = cell_i_begin;
                    for (; begin_k < cell_i_end; begin_k++) {
                        if (fimin <= second_xmax[cellidx[begin_k]]) {
                            break;
                        }
                    }

                    ptrdiff_t end_k = begin_k;
                    for (; end_k < cell_i_end; end_k++) {
                        if (fimax < second_xmin[cellidx[end_k]]) {
                            break;
                        }
                    }

                    if (begin_k >= end_k) {
                        continue;
                    }

                    if (end_k - begin_k < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
                        count += sccd_detail::cell_scalar_count_range_two_lists<first_nxe, second_nxe>(first_aabbs,
                                                                                                       fi,
                                                                                                       second_aabbs,
                                                                                                       second_idx,
                                                                                                       second_elements,
                                                                                                       second_stride,
                                                                                                       ev,
                                                                                                       begin_k,
                                                                                                       end_k,
                                                                                                       cellidx);
                        continue;
                    }

                    for (; begin_k < end_k;) {
                        const ptrdiff_t chunk_len = std::min((ptrdiff_t)AABB_DISJOINT_CHUNK_SIZE, end_k - begin_k);

                        uint32_t dmask[AABB_DISJOINT_CHUNK_SIZE] = {0};

                        const T aminx = first_aabbs[0][fi];
                        const T aminy = first_aabbs[1][fi];
                        const T aminz = first_aabbs[2][fi];
                        const T amaxx = first_aabbs[3][fi];
                        const T amaxy = first_aabbs[4][fi];
                        const T amaxz = first_aabbs[5][fi];

                        alignas(64) T B_minx[AABB_DISJOINT_CHUNK_SIZE];
                        alignas(64) T B_miny[AABB_DISJOINT_CHUNK_SIZE];
                        alignas(64) T B_minz[AABB_DISJOINT_CHUNK_SIZE];
                        alignas(64) T B_maxx[AABB_DISJOINT_CHUNK_SIZE];
                        alignas(64) T B_maxy[AABB_DISJOINT_CHUNK_SIZE];
                        alignas(64) T B_maxz[AABB_DISJOINT_CHUNK_SIZE];

                        for (ptrdiff_t lane = 0; lane < chunk_len; ++lane) {
                            const ptrdiff_t j = cellidx[begin_k + lane];
                            B_minx[lane] = second_aabbs[0][j];
                            B_miny[lane] = second_aabbs[1][j];
                            B_minz[lane] = second_aabbs[2][j];
                            B_maxx[lane] = second_aabbs[3][j];
                            B_maxy[lane] = second_aabbs[4][j];
                            B_maxz[lane] = second_aabbs[5][j];
                        }

                        sccd_detail::tail_fill_B(
                            amaxx, amaxy, amaxz, chunk_len, B_minx, B_miny, B_minz, B_maxx, B_maxy, B_maxz);

                        vaabb_disjoint_one_to_many(aminx,
                                                   aminy,
                                                   aminz,
                                                   amaxx,
                                                   amaxy,
                                                   amaxz,
                                                   B_minx,
                                                   B_miny,
                                                   B_minz,
                                                   B_maxx,
                                                   B_maxy,
                                                   B_maxz,
                                                   dmask);

                        sccd_detail::cell_mask_out_shared_two_lists<first_nxe, second_nxe, T, I>(
                            dmask, chunk_len, begin_k, ev, second_idx, second_elements, second_stride, cellidx);

                        for (ptrdiff_t lane = 0; lane < chunk_len; ++lane) {
                            count += dmask[lane] ? 0 : 1;
                        }

                        begin_k += chunk_len;
                    }
                }

                ccdptr[fi + 1] = count;
            }
        });

        for (ptrdiff_t fi = 0; fi < first_count; fi++) {
            ccdptr[fi + 1] += ccdptr[fi];
        }

        return ccdptr[first_count] > 0;
    }

    template <int first_nxe, int second_nxe, typename T, typename I>
    void cell_collect_overlaps(const int sort_axis,
                               const count_t first_count,
                               T **const SCCD_RESTRICT first_aabbs,
                               I *const SCCD_RESTRICT first_idx,
                               const ptrdiff_t first_stride,
                               I **const SCCD_RESTRICT first_elements,
                               const count_t second_count,
                               T **const SCCD_RESTRICT second_aabbs,
                               I *const SCCD_RESTRICT second_idx,
                               const ptrdiff_t second_stride,
                               I **const SCCD_RESTRICT second_elements,
                               // Cell list
                               const int cell_list_axis,
                               const ptrdiff_t ncells,
                               const T cell_min,
                               const T cell_size,
                               const I *const SCCD_RESTRICT cellptr,
                               const I *const SCCD_RESTRICT cellidx,
                               const ptrdiff_t *const SCCD_RESTRICT ccdptr,
                               I *SCCD_RESTRICT foverlap,
                               I *SCCD_RESTRICT noverlap) {
        const T *const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
        const T *const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
        const T *const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
        const T *const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

        const T *const SCCD_RESTRICT cell_xmin = first_aabbs[cell_list_axis];
        const T *const SCCD_RESTRICT cell_xmax = first_aabbs[3 + cell_list_axis];

        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; fi++) {
                const ptrdiff_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
                if (expected_count == 0) {
                    continue;
                }

                const T fimin = first_xmin[fi];
                const T fimax = first_xmax[fi];
                const I first_idxi = first_idx[fi];

                I *SCCD_RESTRICT const first_local_elements = &foverlap[ccdptr[fi]];
                I *SCCD_RESTRICT const second_local_elements = &noverlap[ccdptr[fi]];

                ptrdiff_t cell_start =
                    std::max((int)0, (int)(floor(nextafter_down((cell_xmin[fi] - cell_min) / cell_size) - 1)));
                ptrdiff_t cell_end =
                    std::min((int)ncells, (int)(floor(nextafter_up((cell_xmax[fi] - cell_min) / cell_size))) + 1);

                I ev[first_nxe];
                for (int v = 0; v < first_nxe; v++) {
                    ev[v] = first_elements[v][first_idxi * first_stride];
                }

                ptrdiff_t count = 0;
                for (ptrdiff_t cell_i = cell_start; cell_i < cell_end; cell_i++) {
                    ptrdiff_t cell_i_begin = cellptr[cell_i];
                    ptrdiff_t cell_i_end = cellptr[cell_i + 1];

                    ptrdiff_t begin_k = cell_i_begin;
                    for (; begin_k < cell_i_end; begin_k++) {
                        if (fimin <= second_xmax[cellidx[begin_k]]) {
                            break;
                        }
                    }

                    ptrdiff_t end_k = begin_k;
                    for (; end_k < cell_i_end; end_k++) {
                        if (fimax < second_xmin[cellidx[end_k]]) {
                            break;
                        }
                    }

                    if (begin_k >= end_k) {
                        continue;
                    }

                    count += sccd_detail::cell_scalar_collect_range_two_lists<first_nxe,
                                                                              second_nxe>(
                        first_aabbs,                     // 0
                        fi,                              // 1
                        first_idxi,                      // 2
                        second_aabbs,                    // 3
                        second_idx,                      // 4
                        second_elements,                 // 5
                        second_stride,                   // 6
                        ev,                              // 7
                        begin_k,                         // 8
                        end_k,                           // 9
                        cellidx,                         // 10
                        &first_local_elements[count],    // 11
                        &second_local_elements[count]);  // 12
                }
                assert(count == expected_count);
            }
        });
    }

    template <int first_nxe, int second_nxe, typename T, typename I>
    bool count_overlaps_with_starts(const int sort_axis,
                                    const count_t first_count,
                                    T **const SCCD_RESTRICT first_aabbs,
                                    I *const SCCD_RESTRICT first_idx,
                                    const ptrdiff_t first_stride,
                                    I **const SCCD_RESTRICT first_elements,
                                    const count_t second_count,
                                    T **const SCCD_RESTRICT second_aabbs,
                                    I *const SCCD_RESTRICT second_idx,
                                    const ptrdiff_t second_stride,
                                    I **const SCCD_RESTRICT second_elements,
                                    const ptrdiff_t ncells,
                                    const T cell_min,
                                    const T cell_size,
                                    const ptrdiff_t *const SCCD_RESTRICT starts,
                                    ptrdiff_t *const SCCD_RESTRICT ccdptr) {
        const T *const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
        const T *const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
        const T *const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
        const T *const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

        ccdptr[0] = 0;

        if (first_xmax[first_count - 1] < second_xmin[0]) return false;

        if (second_xmax[second_count - 1] < first_xmin[0]) return false;

        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            ptrdiff_t cell = std::max(0, (int)floor(nextafter_down((first_xmin[rbegin] - cell_min) / cell_size)) - 1);
            ptrdiff_t ni = starts[cell];

            for (; ni >= 0; ni--) {
                if (first_xmin[0] < second_xmax[ni]) {
                    ni++;
                    break;
                }
            }

#ifndef NDEBUG
            {
                ptrdiff_t test = 0;
                for (; test < second_count; test++) {
                    if (second_xmax[test] > first_xmin[0]) {
                        break;
                    }
                }

                if (test != ni) {
                    printf("start: %lu, test: %lu, ni: %lu\n", starts[cell], test, ni);
                }
                assert(test == ni);
            }
#endif

            for (ptrdiff_t fi = rbegin; fi < rend; fi++) {
                const T fimin = first_xmin[fi];
                const T fimax = first_xmax[fi];
                const I first_idxi = first_idx[fi];

                I ev[first_nxe];
                for (int v = 0; v < first_nxe; v++) {
                    ev[v] = first_elements[v][first_idxi * first_stride];
                }

                ptrdiff_t end = ni;
                sccd_detail::compute_candidate_window_progressive(
                    fimin, fimax, second_xmax, second_xmin, second_count, ni, end);

                if (ni >= end) {
                    continue;
                }

                if (end - ni < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
                    ptrdiff_t count = sccd_detail::scalar_count_range_two_lists<first_nxe, second_nxe>(
                        first_aabbs, fi, second_aabbs, second_idx, second_elements, second_stride, ev, ni, end);
                    ccdptr[fi + 1] = count;
                    continue;
                }

                ptrdiff_t count = 0;
                ptrdiff_t noffset = ni;

                for (; noffset < end;) {
                    const ptrdiff_t chunk_len = std::min((ptrdiff_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    uint32_t dmask[AABB_DISJOINT_CHUNK_SIZE] = {0};

                    sccd_detail::build_disjoint_mask_for_block(second_aabbs,
                                                               noffset,
                                                               chunk_len,
                                                               first_aabbs[0][fi],
                                                               first_aabbs[1][fi],
                                                               first_aabbs[2][fi],
                                                               first_aabbs[3][fi],
                                                               first_aabbs[4][fi],
                                                               first_aabbs[5][fi],
                                                               dmask);

                    sccd_detail::mask_out_shared_two_lists<first_nxe, second_nxe>(
                        dmask, chunk_len, noffset, ev, second_idx, second_elements, second_stride);

                    for (ptrdiff_t lane = 0; lane < chunk_len; ++lane) {
                        count += dmask[lane] ? 0 : 1;
                    }

                    noffset += chunk_len;
                }

                ccdptr[fi + 1] = count;
            }
        });

        for (ptrdiff_t fi = 0; fi < first_count; fi++) {
            ccdptr[fi + 1] += ccdptr[fi];
        }

        return ccdptr[first_count] > 0;
    }

}  // namespace sccd

#endif  // CELL_BROADPHASE_HPP