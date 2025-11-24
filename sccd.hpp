#ifndef LEAN_CCD_HPP
#define LEAN_CCD_HPP

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <vector>
#include <algorithm>
#include <filesystem>
#include <cfloat>

#include "vaabb.h"

using geom_t = float;
using idx_t = int;
using count_t = int;

#ifndef _WIN32
#define SFEM_RESTRICT __restrict__
#else
#define SFEM_RESTRICT __restrict
#endif

inline geom_t lean_nextafter_up(const geom_t x)
{
    return nextafterf(x, FLT_MAX);
}

inline geom_t lean_nextafter_down(const geom_t x)
{
    return nextafterf(x, -FLT_MAX);
}

// inline geom_t lean_nextafter_up(const geom_t x)
// {
//     return nextafter(x,DBL_MAX);
// }

// inline geom_t lean_nextafter_down(const geom_t x)
// {
//     return nextafter(x, DBL_MAX);
// }

inline static bool lean_disjoint(
    const geom_t aminx,
    const geom_t aminy,
    const geom_t aminz,
    const geom_t amaxx,
    const geom_t amaxy,
    const geom_t amaxz,
    const geom_t bminx,
    const geom_t bminy,
    const geom_t bminz,
    const geom_t bmaxx,
    const geom_t bmaxy,
    const geom_t bmaxz)
{
    return 
        aminx > bmaxx |
        aminy > bmaxy |
        aminz > bmaxz |
        bminx > amaxx |
        bminy > amaxy |
        bminz > amaxz;
}

int lean_choose_axis(const size_t n, geom_t** const SFEM_RESTRICT aabb)
{
    geom_t mean[3] = { 0 };
    geom_t var[3] = { 0 };
    for (int d = 0; d < 3; d++) {
        for (size_t i = 0; i < n; i++) {
            const geom_t c = (aabb[d + 3][i] + aabb[d][i]) / 2;
            mean[d] += c;
        }

        mean[d] /= n;
        for (size_t i = 0; i < n; i++) {
            const geom_t c = (aabb[d + 3][i] + aabb[d][i]) / 2;
            var[d] += (c - mean[d]) * (c - mean[d]);
        }
    }

    int fargmax = 0;
    geom_t fmax = var[0];

    for (int d = 1; d < 3; d++) {
        if (fmax < var[d]) {
            fmax = var[d];
            fargmax = d;
        }
    }

    return fargmax;
}

void lean_sort_along_axis(
    const size_t n,
    const int sort_axis,
    geom_t** const SFEM_RESTRICT arrays,
    idx_t* const SFEM_RESTRICT idx,
    geom_t* const SFEM_RESTRICT scratch)
{
    for (size_t i = 0; i < n; i++) {
        idx[i] = i;
    }

    const geom_t* const SFEM_RESTRICT x = arrays[sort_axis];
    // std::sort
    tbb::parallel_sort(idx, idx + n, [x](const idx_t l, const idx_t r) {
        return x[l] < x[r];
    });

    for (int d = 0; d < 6; d++) {
        memcpy(scratch, arrays[d], sizeof(geom_t) * n);
        for (size_t i = 0; i < n; i++) {
            arrays[d][i] = scratch[idx[i]];
        }
    }
}

static void remap_idx(
    const size_t n,
    const idx_t* const SFEM_RESTRICT idx,
    idx_t* const SFEM_RESTRICT remapped)
{
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, n),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t i = r.begin(); i < r.end(); i++) {
                remapped[i] = idx[remapped[i]];
            }
        });
}

template <int first_nxe, int second_nxe>
bool lean_count_overlaps(
    const int sort_axis,
    const count_t first_count,
    geom_t** const SFEM_RESTRICT first_aabbs,
    idx_t* const SFEM_RESTRICT first_idx,
    const size_t first_stride,
    idx_t** const first_elements,
    const count_t second_count,
    geom_t** const SFEM_RESTRICT second_aabbs,
    idx_t* const SFEM_RESTRICT second_idx,
    const size_t second_stride,
    idx_t** const second_elements,
    size_t* const ccdptr)
{
    const geom_t* const SFEM_RESTRICT first_xmin = first_aabbs[sort_axis];
    const geom_t* const SFEM_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
    const geom_t* const SFEM_RESTRICT second_xmin = second_aabbs[sort_axis];
    const geom_t* const SFEM_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

    if (first_xmax[first_count - 1] < second_xmin[0])
        return false;

    if (second_xmax[second_count - 1] < first_xmin[0])
        return false;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, first_count),
        [&](const tbb::blocked_range<size_t>& r) {
            size_t ni = 0;
            for (; ni < second_count; ni++) {
                if (second_xmax[ni] > first_xmin[0]) {
                    break;
                }
            }

            for (size_t fi = r.begin(); fi < r.end(); fi++) {
                const geom_t fimin = first_xmin[fi];
                const geom_t fimax = first_xmax[fi];
                const idx_t first_idxi = first_idx[fi];

                const geom_t fi_min[3] = { first_aabbs[0][fi],
                                           first_aabbs[1][fi],
                                           first_aabbs[2][fi] };

                const geom_t fi_max[3] = { first_aabbs[3][fi],
                                           first_aabbs[4][fi],
                                           first_aabbs[5][fi] };

                idx_t ev[first_nxe];
                for (int v = 0; v < first_nxe; v++) {
                    ev[v] = first_elements[v][first_idxi * first_stride];
                }

                for (; ni < second_count; ni++) {
                    if (fimin < second_xmax[ni]) {
                        break;
                    }
                }

                // Find count potential overlaps
                size_t count = 0;
                size_t noffset = ni;
                for (; noffset < second_count; noffset++) {
                    if (fimax < second_xmin[noffset]) {
                        break;
                    }

                    const geom_t si_min[3] = { second_aabbs[0][noffset],
                                               second_aabbs[1][noffset],
                                               second_aabbs[2][noffset] };

                    const geom_t si_max[3] = { second_aabbs[3][noffset],
                                               second_aabbs[4][noffset],
                                               second_aabbs[5][noffset] };

                    bool skip = false;
                    for (int d = 0; d < 3; d++) {
                        skip |= si_min[d] > fi_max[d] || si_max[d] < fi_min[d];
                    }

                    if (!skip) {
                        if (second_nxe > 1) {
                            idx_t second_ev[second_nxe];
                            for (int v = 0; v < second_nxe; v++) {
                                second_ev[v] =
                                    second_elements[v][noffset * second_stride];
                            }

                            for (int first_v = 0; first_v < first_nxe;
                                 first_v++) {
                                for (int second_v = 0; second_v < second_nxe;
                                     second_v++) {
                                    skip = ev[first_v] == second_ev[second_v]
                                        ? 1
                                        : skip;
                                }
                            }

                        } else {
                            for (int first_v = 0; first_v < first_nxe;
                                 first_v++) {
                                skip = ev[first_v] == second_idx[noffset]
                                    ? 1
                                    : skip;
                            }
                        }
                    }

                    count += skip ? 0 : 1;
                }

                ccdptr[fi + 1] = count;
            }
        });

    for (size_t fi = 0; fi < first_count; fi++) {
        ccdptr[fi + 1] += ccdptr[fi];
    }

    return ccdptr[first_count] > 0;
}

template <int first_nxe, int second_nxe>
void lean_collect_overlaps(
    const int sort_axis,
    const count_t first_count,
    geom_t** const SFEM_RESTRICT first_aabbs,
    idx_t* const SFEM_RESTRICT first_idx,
    const size_t first_stride,
    idx_t** const
        first_elements, // first_elements are accessed through first_idx
    const count_t second_count,
    geom_t** const SFEM_RESTRICT second_aabbs,
    idx_t* const SFEM_RESTRICT second_idx,
    const size_t second_stride,
    idx_t** const second_elements,
    const size_t* const SFEM_RESTRICT ccdptr,
    idx_t* SFEM_RESTRICT foverlap,
    idx_t* SFEM_RESTRICT noverlap)
{
    const geom_t* const SFEM_RESTRICT first_xmin = first_aabbs[sort_axis];
    const geom_t* const SFEM_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
    const geom_t* const SFEM_RESTRICT second_xmin = second_aabbs[sort_axis];
    const geom_t* const SFEM_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

    if (first_xmax[first_count - 1] < second_xmin[0])
        return;

    if (second_xmax[second_count - 1] < first_xmin[0])
        return;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, first_count),
        [&](const tbb::blocked_range<size_t>& r) {
            size_t ni = 0; // This can be improved
            for (; ni < second_count; ni++) {
                if (second_xmax[ni] > first_xmin[0]) {
                    break;
                }
            }

            for (size_t fi = r.begin(); fi < r.end(); fi++) {
                const geom_t fimin = first_xmin[fi];
                const geom_t fimax = first_xmax[fi];
                const idx_t first_idxi = first_idx[fi];

                const geom_t fi_min[3] = { first_aabbs[0][fi],
                                           first_aabbs[1][fi],
                                           first_aabbs[2][fi] };

                const geom_t fi_max[3] = { first_aabbs[3][fi],
                                           first_aabbs[4][fi],
                                           first_aabbs[5][fi] };

                idx_t* SFEM_RESTRICT const first_local_elements =
                    &foverlap[ccdptr[fi]];
                idx_t* SFEM_RESTRICT const second_local_elements =
                    &noverlap[ccdptr[fi]];

                const size_t expected_count = ccdptr[fi + 1] - ccdptr[fi];

                idx_t ev[first_nxe];
                for (int v = 0; v < first_nxe; v++) {
                    ev[v] = first_elements[v][first_idxi * first_stride];
                }

                for (; ni < second_count; ni++) {
                    if (fimin < second_xmax[ni]) {
                        break;
                    }
                }

                // Find count potential overlaps
                size_t count = 0;
                size_t noffset = ni;
                for (; noffset < second_count; noffset++) {
                    if (fimax < second_xmin[noffset]) {
                        break;
                    }

                    const geom_t si_min[3] = { second_aabbs[0][noffset],
                                               second_aabbs[1][noffset],
                                               second_aabbs[2][noffset] };

                    const geom_t si_max[3] = { second_aabbs[3][noffset],
                                               second_aabbs[4][noffset],
                                               second_aabbs[5][noffset] };

                    bool skip = false;
                    for (int d = 0; d < 3; d++) {
                        skip |= si_min[d] > fi_max[d] || si_max[d] < fi_min[d];
                    }

                    const idx_t second_idxi = second_idx[noffset];
                    if (!skip) {
                        if (second_nxe > 1) {
                            idx_t second_ev[second_nxe];
                            for (int v = 0; v < second_nxe; v++) {
                                second_ev[v] =
                                    second_elements[v][second_idxi * second_stride];
                            }

                            for (int first_v = 0; first_v < first_nxe;
                                 first_v++) {
                                for (int second_v = 0; second_v < second_nxe;
                                     second_v++) {
                                    skip = ev[first_v] == second_ev[second_v]
                                        ? 1
                                        : skip;
                                }
                            }

                        } else {
                            for (int first_v = 0; first_v < first_nxe;
                                 first_v++) {
                                skip = ev[first_v] == second_idx[noffset]
                                    ? 1
                                    : skip;
                            }
                        }
                    }

                    if (!skip) {
                        first_local_elements[count] = first_idxi;
                        second_local_elements[count] = second_idxi;
                    }

                    count += skip ? 0 : 1;
                }

                assert(expected_count == count);
            }
        });
}

// --------------------------------------

template <int nxe>
bool lean_count_self_overlaps(
    const int sort_axis,
    const count_t element_count,
    geom_t** const SFEM_RESTRICT aabbs,
    idx_t* const SFEM_RESTRICT idx,
    const size_t stride,
    idx_t** const elements,
    size_t* const ccdptr)
{
    const geom_t* const SFEM_RESTRICT xmin = aabbs[sort_axis];
    const geom_t* const SFEM_RESTRICT xmax = aabbs[3 + sort_axis];

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, element_count),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t fi = r.begin(); fi < r.end(); fi++) {
                const geom_t fimin = xmin[fi];
                const geom_t fimax = xmax[fi];
                const idx_t idxi = idx[fi];
                const geom_t fi_min[3] = { aabbs[0][fi], aabbs[1][fi],
                                           aabbs[2][fi] };
                const geom_t fi_max[3] = { aabbs[3][fi], aabbs[4][fi],
                                           aabbs[5][fi] };

                idx_t ev[nxe];
                for (int v = 0; v < nxe; v++) {
                    ev[v] = elements[v][idxi * stride];
                }

                size_t noffset = fi + 1;
                for (; noffset < element_count; noffset++) {
                    if (fimin < xmax[noffset]) {
                        break;
                    }
                }

                // Count potential overlaps using vectorized disjoint test
                size_t count = 0;

                // Prepare repeated A values for SIMD chunk
                geom_t A_minx[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_miny[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_minz[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxx[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxy[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxz[AABB_DISJOINT_CHUNK_SIZE];
                for (int k = 0; k < AABB_DISJOINT_CHUNK_SIZE; ++k) {
                    A_minx[k] = fi_min[0];
                    A_miny[k] = fi_min[1];
                    A_minz[k] = fi_min[2];
                    A_maxx[k] = fi_max[0];
                    A_maxy[k] = fi_max[1];
                    A_maxz[k] = fi_max[2];
                }

                // Determine end of candidate range using original loop semantics
                size_t end = noffset;
                for (; end < element_count; end++) {
                    if (fimax < xmin[end]) {
                        break;
                    }
                }

                for (; noffset < end;) {
                    const size_t chunk_len =
                        std::min((size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    geom_t B_minx[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_miny[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_minz[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_maxx[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_maxy[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_maxz[AABB_DISJOINT_CHUNK_SIZE];

                    for (size_t lane = 0; lane < chunk_len; ++lane) {
                        const size_t j = noffset + lane;
                        B_minx[lane] = aabbs[0][j];
                        B_miny[lane] = aabbs[1][j];
                        B_minz[lane] = aabbs[2][j];
                        B_maxx[lane] = aabbs[3][j];
                        B_maxy[lane] = aabbs[4][j];
                        B_maxz[lane] = aabbs[5][j];
                    }
                    // Tail-safe fill: force disjoint
                    for (size_t lane = chunk_len; lane < AABB_DISJOINT_CHUNK_SIZE; ++lane) {
                        B_minx[lane] = A_maxx[0] + 1;
                        B_miny[lane] = A_maxy[0] + 1;
                        B_minz[lane] = A_maxz[0] + 1;
                        B_maxx[lane] = A_maxx[0];
                        B_maxy[lane] = A_maxy[0];
                        B_maxz[lane] = A_maxz[0];
                    }

                    // Candidate mask = not disjoint
                    const unsigned disj_bits = vdisjoint_mask(
                        A_minx, A_miny, A_minz, A_maxx, A_maxy, A_maxz, B_minx,
                        B_miny, B_minz, B_maxx, B_maxy, B_maxz, chunk_len);
                    unsigned cand_bits =
                        ((chunk_len >= 32 ? 0xFFFFFFFFu : ((1u << chunk_len) - 1u)))
                        & ~disj_bits;
                    // Shared-vertex mask using original indices
                    const unsigned share_bits = vshare_edge_shares_vertex_mask(
                        idx, noffset, chunk_len, (int)stride,
                        (const int*)elements[0], (const int*)elements[1],
                        ev[0], ev[1]);
                    // Valid = candidate and not shared
                    unsigned valid_bits = cand_bits & ~share_bits;
                    count += __builtin_popcount(valid_bits);

                    noffset += chunk_len;
                }

                ccdptr[fi + 1] = count;
            }
        });

    for (size_t fi = 0; fi < element_count; fi++) {
        ccdptr[fi + 1] += ccdptr[fi];
    }

    return ccdptr[element_count] > 0;
}

template <int nxe>
void lean_collect_self_overlaps(
    const int sort_axis,
    const count_t element_count,
    geom_t** const SFEM_RESTRICT aabbs,
    idx_t* const SFEM_RESTRICT idx,
    const size_t stride,
    idx_t** const elements,
    const size_t* const SFEM_RESTRICT ccdptr,
    idx_t* SFEM_RESTRICT foverlap,
    idx_t* SFEM_RESTRICT noverlap)
{
    const geom_t* const SFEM_RESTRICT xmin = aabbs[sort_axis];
    const geom_t* const SFEM_RESTRICT xmax = aabbs[3 + sort_axis];

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, element_count),
        [&](const tbb::blocked_range<size_t>& r) {
            for (size_t fi = r.begin(); fi < r.end(); fi++) {
                const geom_t fimin = xmin[fi];
                const geom_t fimax = xmax[fi];
                const idx_t idxi = idx[fi];

                const geom_t fi_min[3] = { aabbs[0][fi], aabbs[1][fi],
                                           aabbs[2][fi] };

                const geom_t fi_max[3] = { aabbs[3][fi], aabbs[4][fi],
                                           aabbs[5][fi] };

                idx_t* SFEM_RESTRICT const first_local_elements =
                    &foverlap[ccdptr[fi]];

                idx_t* SFEM_RESTRICT const second_local_elements =
                    &noverlap[ccdptr[fi]];

                const size_t expected_count = ccdptr[fi + 1] - ccdptr[fi];

                idx_t ev[nxe];
                for (int v = 0; v < nxe; v++) {
                    ev[v] = elements[v][idxi * stride];
                }

                size_t noffset = fi + 1;
                for (; noffset < element_count; noffset++) {
                    if (fimin < xmax[noffset]) {
                        break;
                    }
                }

                // Collect potential overlaps using vectorized disjoint test
                size_t count = 0;

                // Prepare repeated A values for SIMD chunk
                geom_t A_minx[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_miny[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_minz[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxx[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxy[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxz[AABB_DISJOINT_CHUNK_SIZE];
                for (int k = 0; k < AABB_DISJOINT_CHUNK_SIZE; ++k) {
                    A_minx[k] = fi_min[0];
                    A_miny[k] = fi_min[1];
                    A_minz[k] = fi_min[2];
                    A_maxx[k] = fi_max[0];
                    A_maxy[k] = fi_max[1];
                    A_maxz[k] = fi_max[2];
                }

                // Determine end of candidate range using original loop semantics
                size_t end = noffset;
                for (; end < element_count; end++) {
                    if (fimax < xmin[end]) {
                        break;
                    }
                }

                for (; noffset < end;) {
                    const size_t chunk_len =
                        std::min((size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    geom_t B_minx[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_miny[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_minz[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_maxx[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_maxy[AABB_DISJOINT_CHUNK_SIZE];
                    geom_t B_maxz[AABB_DISJOINT_CHUNK_SIZE];

                    for (size_t lane = 0; lane < chunk_len; ++lane) {
                        const size_t j = noffset + lane;
                        B_minx[lane] = aabbs[0][j];
                        B_miny[lane] = aabbs[1][j];
                        B_minz[lane] = aabbs[2][j];
                        B_maxx[lane] = aabbs[3][j];
                        B_maxy[lane] = aabbs[4][j];
                        B_maxz[lane] = aabbs[5][j];
                    }
                    // Tail-safe fill: force disjoint
                    for (size_t lane = chunk_len; lane < AABB_DISJOINT_CHUNK_SIZE; ++lane) {
                        B_minx[lane] = A_maxx[0] + 1;
                        B_miny[lane] = A_maxy[0] + 1;
                        B_minz[lane] = A_maxz[0] + 1;
                        B_maxx[lane] = A_maxx[0];
                        B_maxy[lane] = A_maxy[0];
                        B_maxz[lane] = A_maxz[0];
                    }

                    // Candidate mask = not disjoint
                    const unsigned disj_bits = vdisjoint_mask(
                        A_minx, A_miny, A_minz, A_maxx, A_maxy, A_maxz, B_minx,
                        B_miny, B_minz, B_maxx, B_maxy, B_maxz, chunk_len);
                    unsigned cand_bits =
                        ((chunk_len >= 32 ? 0xFFFFFFFFu : ((1u << chunk_len) - 1u)))
                        & ~disj_bits;

                    // Compute share mask via SIMD helper and write valid pairs
                    const unsigned share_bits = vshare_edge_shares_vertex_mask(
                        idx, noffset, chunk_len, (int)stride,
                        (const int*)elements[0], (const int*)elements[1],
                        ev[0], ev[1]);
                        
                    unsigned valid_bits = cand_bits & ~share_bits;
                    while (valid_bits) {
                        unsigned lane = __builtin_ctz(valid_bits);
                        valid_bits &= valid_bits - 1;
                        const size_t j = noffset + lane;
                        const idx_t jidx = idx[j];
                        first_local_elements[count] = std::min(idxi, jidx);
                        second_local_elements[count] = std::max(idxi, jidx);
                        count += 1;
                    }

                    noffset += chunk_len;
                }

                assert(expected_count == count);
            }
        });
}

#endif
