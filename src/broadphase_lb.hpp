#ifndef SCCC_BROADPHASE_LB_HPP
#define SCCC_BROADPHASE_LB_HPP

#include "broadphase.hpp"

namespace sccd {
    template <int first_nxe, int second_nxe, typename T, typename I>
    bool count_overlaps_lb(const int sort_axis,
                           const ptrdiff_t first_count,
                           T **const SCCD_RESTRICT first_aabbs,
                           I *const SCCD_RESTRICT first_idx,
                           const ptrdiff_t first_stride,
                           I **const SCCD_RESTRICT first_elements,
                           const ptrdiff_t second_count,
                           T **const SCCD_RESTRICT second_aabbs,
                           I *const SCCD_RESTRICT second_idx,
                           const ptrdiff_t second_stride,
                           I **const SCCD_RESTRICT second_elements,
                           const T *const SCCD_RESTRICT lb,
                           ptrdiff_t *const SCCD_RESTRICT ccdptr) {
        const T *const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
        const T *const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
        const T *const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
        const T *const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

        ccdptr[0] = 0;

        if (first_xmax[first_count - 1] < second_xmin[0]) return false;

        if (second_xmax[second_count - 1] < first_xmin[0]) return false;

        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            ptrdiff_t ni = std::distance(lb, std::lower_bound(lb, lb + second_count, first_xmin[rbegin]));

#ifndef NDEBUG
            ptrdiff_t test_ni = 0;
            for (; test_ni < second_count; test_ni++) {
                if (second_xmax[test_ni] > first_xmin[rbegin]) {
                    break;
                }
            }

            assert(ni == test_ni);
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
                    const ptrdiff_t chunk_len = sccd::min((ptrdiff_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

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

        sccd::parallel_cum_sum_br(ccdptr, ccdptr + first_count + 1);
        return ccdptr[first_count] > 0;
    }

    template <int first_nxe, int second_nxe, typename T, typename I>
    void collect_overlaps_lb(const int sort_axis,
                             const ptrdiff_t first_count,
                             T **const SCCD_RESTRICT first_aabbs,
                             I *const SCCD_RESTRICT first_idx,
                             const ptrdiff_t first_stride,
                             I **SCCD_RESTRICT const first_elements,
                             const ptrdiff_t second_count,
                             T **const SCCD_RESTRICT second_aabbs,
                             I *const SCCD_RESTRICT second_idx,
                             const ptrdiff_t second_stride,
                             I **SCCD_RESTRICT const second_elements,
                             const T *const SCCD_RESTRICT lb,
                             const ptrdiff_t *const SCCD_RESTRICT ccdptr,
                             I *SCCD_RESTRICT foverlap,
                             I *SCCD_RESTRICT noverlap) {
        const T *const SCCD_RESTRICT first_xmin = first_aabbs[sort_axis];
        const T *const SCCD_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
        const T *const SCCD_RESTRICT second_xmin = second_aabbs[sort_axis];
        const T *const SCCD_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

        if (first_xmax[first_count - 1] < second_xmin[0]) return;

        if (second_xmax[second_count - 1] < first_xmin[0]) return;

        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            ptrdiff_t ni = std::distance(lb, std::lower_bound(lb, lb + second_count, first_xmin[rbegin]));

#ifndef NDEBUG
            ptrdiff_t test_ni = 0;
            for (; test_ni < second_count; test_ni++) {
                if (second_xmax[test_ni] > first_xmin[rbegin]) {
                    break;
                }
            }

            assert(ni == test_ni);
#endif

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
                    ptrdiff_t count =
                        sccd_detail::scalar_collect_range_two_lists<first_nxe, second_nxe>(first_aabbs,
                                                                                           fi,
                                                                                           first_idxi,
                                                                                           second_aabbs,
                                                                                           second_idx,
                                                                                           second_elements,
                                                                                           second_stride,
                                                                                           ev,
                                                                                           ni,
                                                                                           end,
                                                                                           first_local_elements,
                                                                                           second_local_elements);
                    assert(expected_count == count);
                    continue;
                }

                ptrdiff_t count = 0;
                ptrdiff_t noffset = ni;

                for (; noffset < end;) {
                    const ptrdiff_t chunk_len = sccd::min((ptrdiff_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

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
                        if (dmask[lane]) continue;
                        const ptrdiff_t j = noffset + lane;
                        first_local_elements[count] = first_idxi;
                        second_local_elements[count] = second_idx[j];
                        count += 1;
                    }

                    noffset += chunk_len;
                }

                assert(expected_count == count);
            }
        });
    }

}  // namespace sccd
#endif  // SCCC_BROADPHASE_LB_HPP
