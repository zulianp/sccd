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

namespace sccd_detail {

// Compute begin/end candidate window indices for the current AABB
static inline void compute_candidate_window(
    const geom_t fimin,
    const geom_t fimax,
    const geom_t* const SFEM_RESTRICT second_xmax,
    const geom_t* const SFEM_RESTRICT second_xmin,
    const size_t second_count,
    const size_t start_guess,
    size_t& out_begin,
    size_t& out_end)
{
    size_t begin = start_guess;
    for (; begin < second_count; ++begin) {
        if (fimin < second_xmax[begin]) {
            break;
        }
    }
    size_t end = begin;
    for (; end < second_count; ++end) {
        if (fimax < second_xmin[end]) {
            break;
        }
    }
    out_begin = begin;
    out_end = end;
}

// Replicate A (fi) box values into SoA scratch for SIMD chunk processing
static inline void prepare_A_replicated(
    geom_t** const SFEM_RESTRICT aabbs,
    const size_t fi,
    geom_t* const SFEM_RESTRICT A_minx,
    geom_t* const SFEM_RESTRICT A_miny,
    geom_t* const SFEM_RESTRICT A_minz,
    geom_t* const SFEM_RESTRICT A_maxx,
    geom_t* const SFEM_RESTRICT A_maxy,
    geom_t* const SFEM_RESTRICT A_maxz)
{
    const geom_t aminx = aabbs[0][fi];
    const geom_t aminy = aabbs[1][fi];
    const geom_t aminz = aabbs[2][fi];
    const geom_t amaxx = aabbs[3][fi];
    const geom_t amaxy = aabbs[4][fi];
    const geom_t amaxz = aabbs[5][fi];
    for (int k = 0; k < AABB_DISJOINT_CHUNK_SIZE; ++k) {
        A_minx[k] = aminx;
        A_miny[k] = aminy;
        A_minz[k] = aminz;
        A_maxx[k] = amaxx;
        A_maxy[k] = amaxy;
        A_maxz[k] = amaxz;
    }
}

// Load a block of B boxes [start, start+len) into SoA buffers
static inline void prepare_B_block(
    geom_t** const SFEM_RESTRICT aabbs,
    const size_t start,
    const size_t len,
    geom_t* const SFEM_RESTRICT B_minx,
    geom_t* const SFEM_RESTRICT B_miny,
    geom_t* const SFEM_RESTRICT B_minz,
    geom_t* const SFEM_RESTRICT B_maxx,
    geom_t* const SFEM_RESTRICT B_maxy,
    geom_t* const SFEM_RESTRICT B_maxz)
{
    for (size_t lane = 0; lane < len; ++lane) {
        const size_t j = start + lane;
        B_minx[lane] = aabbs[0][j];
        B_miny[lane] = aabbs[1][j];
        B_minz[lane] = aabbs[2][j];
        B_maxx[lane] = aabbs[3][j];
        B_maxy[lane] = aabbs[4][j];
        B_maxz[lane] = aabbs[5][j];
    }
}

// Force remaining lanes to be disjoint by setting B outside A
static inline void tail_fill_B(
    const geom_t amaxx0,
    const geom_t amaxy0,
    const geom_t amaxz0,
    const size_t len,
    geom_t* const SFEM_RESTRICT B_minx,
    geom_t* const SFEM_RESTRICT B_miny,
    geom_t* const SFEM_RESTRICT B_minz,
    geom_t* const SFEM_RESTRICT B_maxx,
    geom_t* const SFEM_RESTRICT B_maxy,
    geom_t* const SFEM_RESTRICT B_maxz)
{
    for (size_t lane = len; lane < AABB_DISJOINT_CHUNK_SIZE; ++lane) {
        B_minx[lane] = amaxx0 + 1;
        B_miny[lane] = amaxy0 + 1;
        B_minz[lane] = amaxz0 + 1;
        B_maxx[lane] = amaxx0;
        B_maxy[lane] = amaxy0;
        B_maxz[lane] = amaxz0;
    }
}

// Prepare B block (with tail fill) and compute disjoint mask in one call
static inline void build_disjoint_mask_for_block(
    geom_t** const SFEM_RESTRICT second_aabbs,
    const size_t start,
    const size_t chunk_len,
    const geom_t* const SFEM_RESTRICT A_minx,
    const geom_t* const SFEM_RESTRICT A_miny,
    const geom_t* const SFEM_RESTRICT A_minz,
    const geom_t* const SFEM_RESTRICT A_maxx,
    const geom_t* const SFEM_RESTRICT A_maxy,
    const geom_t* const SFEM_RESTRICT A_maxz,
    const geom_t amaxx0,
    const geom_t amaxy0,
    const geom_t amaxz0,
    uint32_t* const SFEM_RESTRICT mask_out)
{
    geom_t B_minx[AABB_DISJOINT_CHUNK_SIZE];
    geom_t B_miny[AABB_DISJOINT_CHUNK_SIZE];
    geom_t B_minz[AABB_DISJOINT_CHUNK_SIZE];
    geom_t B_maxx[AABB_DISJOINT_CHUNK_SIZE];
    geom_t B_maxy[AABB_DISJOINT_CHUNK_SIZE];
    geom_t B_maxz[AABB_DISJOINT_CHUNK_SIZE];

    prepare_B_block(
        second_aabbs, start, chunk_len, B_minx, B_miny, B_minz, B_maxx, B_maxy,
        B_maxz);
    tail_fill_B(
        amaxx0, amaxy0, amaxz0, chunk_len, B_minx, B_miny, B_minz, B_maxx,
        B_maxy, B_maxz);
    vdisjoint(
        A_minx, A_miny, A_minz, A_maxx, A_maxy, A_maxz, B_minx, B_miny, B_minz,
        B_maxx, B_maxy, B_maxz, mask_out);
}

template <int nxe>
static inline void load_ev(
    idx_t** const SFEM_RESTRICT elements,
    const idx_t elem_idx,
    const size_t stride,
    idx_t (&out)[nxe])
{
    for (int v = 0; v < nxe; ++v) {
        out[v] = elements[v][elem_idx * stride];
    }
}

template <int n1, int n2>
static inline bool shares_vertex(const idx_t (&a)[n1], const idx_t (&b)[n2])
{
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            if (a[i] == b[j]) {
                return true;
            }
        }
    }
    return false;
}

// -----------------------------
// Scalar fallbacks (two lists)
template <int F, int S>
static inline size_t scalar_count_range_two_lists(
    geom_t** const SFEM_RESTRICT first_aabbs,
    const size_t fi,
    geom_t** const SFEM_RESTRICT second_aabbs,
    const idx_t* const SFEM_RESTRICT second_idx,
    idx_t** const SFEM_RESTRICT second_elements,
    const size_t second_stride,
    const idx_t (&ev)[F],
    const size_t begin,
    const size_t end)
{
    size_t count = 0;
    const geom_t aminx = first_aabbs[0][fi];
    const geom_t aminy = first_aabbs[1][fi];
    const geom_t aminz = first_aabbs[2][fi];
    const geom_t amaxx = first_aabbs[3][fi];
    const geom_t amaxy = first_aabbs[4][fi];
    const geom_t amaxz = first_aabbs[5][fi];
    for (size_t j = begin; j < end; ++j) {
        if (disjoint(
                aminx, aminy, aminz, amaxx, amaxy, amaxz, second_aabbs[0][j],
                second_aabbs[1][j], second_aabbs[2][j], second_aabbs[3][j],
                second_aabbs[4][j], second_aabbs[5][j])) {
            continue;
        }
        bool share = false;
        if constexpr (S > 1) {
            const idx_t jidx = second_idx[j];
            idx_t sev[S];
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

template <int F, int S>
static inline size_t scalar_collect_range_two_lists(
    geom_t** const SFEM_RESTRICT first_aabbs,
    const size_t fi,
    const idx_t first_idxi,
    geom_t** const SFEM_RESTRICT second_aabbs,
    const idx_t* const SFEM_RESTRICT second_idx,
    idx_t** const SFEM_RESTRICT second_elements,
    const size_t second_stride,
    const idx_t (&ev)[F],
    const size_t begin,
    const size_t end,
    idx_t* const SFEM_RESTRICT first_out,
    idx_t* const SFEM_RESTRICT second_out)
{
    size_t count = 0;
    const geom_t aminx = first_aabbs[0][fi];
    const geom_t aminy = first_aabbs[1][fi];
    const geom_t aminz = first_aabbs[2][fi];
    const geom_t amaxx = first_aabbs[3][fi];
    const geom_t amaxy = first_aabbs[4][fi];
    const geom_t amaxz = first_aabbs[5][fi];
    for (size_t j = begin; j < end; ++j) {
        if (disjoint(
                aminx, aminy, aminz, amaxx, amaxy, amaxz, second_aabbs[0][j],
                second_aabbs[1][j], second_aabbs[2][j], second_aabbs[3][j],
                second_aabbs[4][j], second_aabbs[5][j])) {
            continue;
        }
        bool share = false;
        const idx_t jidx = second_idx[j];
        if constexpr (S > 1) {
            idx_t sev[S];
            for (int v = 0; v < S; ++v) {
                sev[v] = second_elements[v][jidx * second_stride];
            }
            share = shares_vertex<F, S>(ev, sev);
        } else {
            for (int a = 0; a < F; ++a) {
                if (ev[a] == jidx) {
                    share = true;
                    break;
                }
            }
        }
        if (!share) {
            first_out[count] = first_idxi;
            second_out[count] = jidx;
            count += 1;
        }
    }
    return count;
}

// -----------------------------
// Scalar fallbacks (self)
template <int N>
static inline size_t scalar_count_range_self(
    geom_t** const SFEM_RESTRICT aabbs,
    const size_t fi,
    idx_t** const SFEM_RESTRICT elements,
    const idx_t* const SFEM_RESTRICT idx,
    const size_t stride,
    const idx_t (&ev)[N],
    const size_t begin,
    const size_t end)
{
    size_t count = 0;
    const geom_t aminx = aabbs[0][fi];
    const geom_t aminy = aabbs[1][fi];
    const geom_t aminz = aabbs[2][fi];
    const geom_t amaxx = aabbs[3][fi];
    const geom_t amaxy = aabbs[4][fi];
    const geom_t amaxz = aabbs[5][fi];
    for (size_t j = begin; j < end; ++j) {
        if (disjoint(
                aminx, aminy, aminz, amaxx, amaxy, amaxz, aabbs[0][j], aabbs[1][j],
                aabbs[2][j], aabbs[3][j], aabbs[4][j], aabbs[5][j])) {
            continue;
        }
        const idx_t jidx = idx[j];
        idx_t sev[N];
        load_ev<N>(elements, jidx, stride, sev);
        const bool share = shares_vertex<N, N>(ev, sev);
        count += share ? 0 : 1;
    }
    return count;
}

template <int N>
static inline size_t scalar_collect_range_self(
    geom_t** const SFEM_RESTRICT aabbs,
    const size_t fi,
    const idx_t idxi,
    idx_t** const SFEM_RESTRICT elements,
    const idx_t* const SFEM_RESTRICT idx,
    const size_t stride,
    const idx_t (&ev)[N],
    const size_t begin,
    const size_t end,
    idx_t* const SFEM_RESTRICT first_out,
    idx_t* const SFEM_RESTRICT second_out)
{
    size_t count = 0;
    const geom_t aminx = aabbs[0][fi];
    const geom_t aminy = aabbs[1][fi];
    const geom_t aminz = aabbs[2][fi];
    const geom_t amaxx = aabbs[3][fi];
    const geom_t amaxy = aabbs[4][fi];
    const geom_t amaxz = aabbs[5][fi];
    for (size_t j = begin; j < end; ++j) {
        if (disjoint(
                aminx, aminy, aminz, amaxx, amaxy, amaxz, aabbs[0][j], aabbs[1][j],
                aabbs[2][j], aabbs[3][j], aabbs[4][j], aabbs[5][j])) {
            continue;
        }
        const idx_t jidx = idx[j];
        idx_t sev[N];
        load_ev<N>(elements, jidx, stride, sev);
        if (!shares_vertex<N, N>(ev, sev)) {
            first_out[count] = std::min(idxi, jidx);
            second_out[count] = std::max(idxi, jidx);
            count += 1;
        }
    }
    return count;
}

} // namespace sccd_detail
template <int first_nxe, int second_nxe>
bool lean_count_overlaps(
    const int sort_axis,
    const count_t first_count,
    geom_t** const SFEM_RESTRICT first_aabbs,
    idx_t* const SFEM_RESTRICT first_idx,
    const size_t first_stride,
    idx_t** const SFEM_RESTRICT first_elements,
    const count_t second_count,
    geom_t** const SFEM_RESTRICT second_aabbs,
    idx_t* const SFEM_RESTRICT second_idx,
    const size_t second_stride,
    idx_t** const SFEM_RESTRICT second_elements,
    size_t* const SFEM_RESTRICT ccdptr)
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

                // Determine candidate window [ni, end)
                size_t end = ni;
                sccd_detail::compute_candidate_window(
                    fimin, fimax, second_xmax, second_xmin, second_count, ni, ni, end);

                if (ni >= end) {
                    continue;
                }
                // Scalar fallback for small candidate ranges
                if (end - ni < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
                    size_t count = sccd_detail::scalar_count_range_two_lists<
                        first_nxe, second_nxe>(first_aabbs, fi, second_aabbs,
                                              second_idx, second_elements,
                                              second_stride, ev, ni, end);
                    ccdptr[fi + 1] = count;
                    continue;
                }

                // Vectorized count of potential overlaps using vdisjoint
                size_t count = 0;
                size_t noffset = ni;

                // Prepare repeated A values for SIMD chunk
                geom_t A_minx[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_miny[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_minz[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxx[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxy[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxz[AABB_DISJOINT_CHUNK_SIZE];
                sccd_detail::prepare_A_replicated(
                    first_aabbs, fi, A_minx, A_miny, A_minz, A_maxx, A_maxy, A_maxz);

                for (; noffset < end;) {
                    const size_t chunk_len = std::min(
                        (size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    uint32_t dmask[AABB_DISJOINT_CHUNK_SIZE];
                    sccd_detail::build_disjoint_mask_for_block(
                        second_aabbs, noffset, chunk_len, A_minx, A_miny, A_minz,
                        A_maxx, A_maxy, A_maxz, A_maxx[0], A_maxy[0], A_maxz[0],
                        dmask);

                    for (size_t lane = 0; lane < chunk_len; ++lane) {
                        if (dmask[lane]) {
                            continue; // disjoint
                        }
                        bool skip = false; // share-a-vertex
                        if (second_nxe > 1) {
                            const idx_t jidx = second_idx[noffset + lane];
                            idx_t second_ev[second_nxe];
                            for (int v = 0; v < second_nxe; v++) {
                                second_ev[v] =
                                    second_elements[v][jidx * second_stride];
                            }
                            for (int first_v = 0; first_v < first_nxe && !skip;
                                 first_v++) {
                                for (int second_v = 0; second_v < second_nxe;
                                     second_v++) {
                                    if (ev[first_v] == second_ev[second_v]) {
                                        skip = true;
                                        break;
                                    }
                                }
                            }
                        } else {
                            for (int first_v = 0; first_v < first_nxe;
                                 first_v++) {
                                if (ev[first_v] == second_idx[noffset + lane]) {
                                    skip = true;
                                    break;
                                }
                            }
                        }
                        count += skip ? 0 : 1;
                    }

                    noffset += chunk_len;
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
    idx_t** SFEM_RESTRICT const
        first_elements, // first_elements are accessed through first_idx
    const count_t second_count,
    geom_t** const SFEM_RESTRICT second_aabbs,
    idx_t* const SFEM_RESTRICT second_idx,
    const size_t second_stride,
    idx_t** SFEM_RESTRICT const second_elements,
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
                const size_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
                if(expected_count == 0) {
                    continue;
                }

                const geom_t fimin = first_xmin[fi];
                const geom_t fimax = first_xmax[fi];
                const idx_t first_idxi = first_idx[fi];

                idx_t* SFEM_RESTRICT const first_local_elements =
                    &foverlap[ccdptr[fi]];
                idx_t* SFEM_RESTRICT const second_local_elements =
                    &noverlap[ccdptr[fi]];

                idx_t ev[first_nxe];
                for (int v = 0; v < first_nxe; v++) {
                    ev[v] = first_elements[v][first_idxi * first_stride];
                }

                for (; ni < second_count; ni++) {
                    if (fimin < second_xmax[ni]) {
                        break;
                    }
                }

                // Determine end of candidate range
                size_t end = ni;
                for (; end < second_count; end++) {
                    if (fimax < second_xmin[end]) {
                        break;
                    }
                }

                if (ni >= end) {
                    continue;
                }

                // Scalar fallback for small candidate ranges
                if (end - ni < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
                    size_t count = sccd_detail::scalar_collect_range_two_lists<
                        first_nxe, second_nxe>(
                        first_aabbs, fi, first_idxi, second_aabbs, second_idx,
                        second_elements, second_stride, ev, ni, end,
                        first_local_elements, second_local_elements);
                    assert(expected_count == count);
                    continue;
                }

                // Vectorized collect of potential overlaps using vdisjoint
                size_t count = 0;
                size_t noffset = ni;

                // Prepare repeated A values for SIMD chunk
                geom_t A_minx[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_miny[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_minz[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxx[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxy[AABB_DISJOINT_CHUNK_SIZE];
                geom_t A_maxz[AABB_DISJOINT_CHUNK_SIZE];
                sccd_detail::prepare_A_replicated(
                    first_aabbs, fi, A_minx, A_miny, A_minz, A_maxx, A_maxy, A_maxz);

                for (; noffset < end;) {
                    const size_t chunk_len = std::min(
                        (size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    uint32_t dmask[AABB_DISJOINT_CHUNK_SIZE];
                    sccd_detail::build_disjoint_mask_for_block(
                        second_aabbs, noffset, chunk_len, A_minx, A_miny, A_minz,
                        A_maxx, A_maxy, A_maxz, A_maxx[0], A_maxy[0], A_maxz[0],
                        dmask);

                    for (size_t lane = 0; lane < chunk_len; ++lane) {
                        if (dmask[lane]) {
                            continue; // disjoint
                        }
                        bool skip = false; // share-a-vertex
                        const size_t j = noffset + lane;
                        const idx_t second_idxi = second_idx[j];
                        if (second_nxe > 1) {
                            idx_t second_ev[second_nxe];
                            for (int v = 0; v < second_nxe; v++) {
                                second_ev[v] = second_elements[v]
                                                              [second_idxi
                                                               * second_stride];
                            }
                            for (int first_v = 0; first_v < first_nxe && !skip;
                                 first_v++) {
                                for (int second_v = 0; second_v < second_nxe;
                                     second_v++) {
                                    if (ev[first_v] == second_ev[second_v]) {
                                        skip = true;
                                        break;
                                    }
                                }
                            }
                        } else {
                            for (int first_v = 0; first_v < first_nxe;
                                 first_v++) {
                                if (ev[first_v] == second_idxi) {
                                    skip = true;
                                    break;
                                }
                            }
                        }

                        if (!skip) {
                            first_local_elements[count] = first_idxi;
                            second_local_elements[count] = second_idxi;
                            count += 1;
                        }
                    }

                    noffset += chunk_len;
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
    idx_t** const SFEM_RESTRICT elements,
    size_t* const SFEM_RESTRICT ccdptr)
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

                // Determine end of candidate range using original loop
                // semantics
                size_t end = noffset;
                for (; end < element_count; end++) {
                    if (fimax < xmin[end]) {
                        break;
                    }
                }

                if (noffset >= end) {
                    ccdptr[fi + 1] = 0;
                    continue;
                }

                // Scalar fallback for small candidate ranges
                if (end - noffset < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
                    size_t count =
                        sccd_detail::scalar_count_range_self<nxe>(
                            aabbs, fi, elements, idx, stride, ev, noffset, end);
                    ccdptr[fi + 1] = count;
                    continue;
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
                sccd_detail::prepare_A_replicated(
                    aabbs, fi, A_minx, A_miny, A_minz, A_maxx, A_maxy, A_maxz);

                for (; noffset < end;) {
                    const size_t chunk_len = std::min(
                        (size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    // Disjoint array mask -> per-lane skip logic
                    uint32_t mask[AABB_DISJOINT_CHUNK_SIZE];
                    sccd_detail::build_disjoint_mask_for_block(
                        aabbs, noffset, chunk_len, A_minx, A_miny, A_minz, A_maxx,
                        A_maxy, A_maxz, A_maxx[0], A_maxy[0], A_maxz[0], mask);

                    for (size_t lane = 0; lane < chunk_len; ++lane) {
                        if (mask[lane]) {
                            continue; // disjoint
                        }
                        const size_t j = noffset + lane;
                        const idx_t jidx = idx[j];
                        idx_t second_ev[nxe];
                        for (int v = 0; v < nxe; v++) {
                            second_ev[v] = elements[v][jidx * stride];
                        }
                        bool share = false;
                        for (int a = 0; a < nxe && !share; ++a) {
                            for (int b = 0; b < nxe; ++b) {
                                if (ev[a] == second_ev[b]) {
                                    share = true;
                                    break;
                                }
                            }
                        }
                        count += share ? 0 : 1;
                    }

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
                const size_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
                if(expected_count == 0) {
                    continue;
                }

                const geom_t fimin = xmin[fi];
                const geom_t fimax = xmax[fi];
                const idx_t idxi = idx[fi];

                idx_t* SFEM_RESTRICT const first_local_elements =
                    &foverlap[ccdptr[fi]];

                idx_t* SFEM_RESTRICT const second_local_elements =
                    &noverlap[ccdptr[fi]];

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

                size_t end = noffset;
                for (; end < element_count; end++) {
                    if (fimax < xmin[end]) {
                        break;
                    }
                }

                if(noffset >= end) {
                    continue;
                }

                // Scalar fallback for small candidate ranges
                if (end - noffset < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
                    size_t count = sccd_detail::scalar_collect_range_self<nxe>(
                        aabbs, fi, idxi, elements, idx, stride, ev, noffset, end,
                        first_local_elements, second_local_elements);
                    assert(expected_count == count);
                    continue;
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
                sccd_detail::prepare_A_replicated(
                    aabbs, fi, A_minx, A_miny, A_minz, A_maxx, A_maxy, A_maxz);

                // Determine end of candidate range using original loop
                // semantics
              
                for (; noffset < end;) {
                    const size_t chunk_len = std::min(
                        (size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

                    // Disjoint array mask -> per-lane skip logic and write pairs
                    uint32_t mask[AABB_DISJOINT_CHUNK_SIZE];
                    sccd_detail::build_disjoint_mask_for_block(
                        aabbs, noffset, chunk_len, A_minx, A_miny, A_minz, A_maxx,
                        A_maxy, A_maxz, A_maxx[0], A_maxy[0], A_maxz[0], mask);

                    for (size_t lane = 0; lane < chunk_len; ++lane) {
                        if (mask[lane]) {
                            continue; // disjoint
                        }
                        const size_t j = noffset + lane;
                        const idx_t jidx = idx[j];
                        idx_t second_ev[nxe];
                        for (int v = 0; v < nxe; v++) {
                            second_ev[v] = elements[v][jidx * stride];
                        }
                        bool share = false;
                        for (int a = 0; a < nxe && !share; ++a) {
                            for (int b = 0; b < nxe; ++b) {
                                if (ev[a] == second_ev[b]) {
                                    share = true;
                                    break;
                                }
                            }
                        }
                        if (!share) {
                            first_local_elements[count] = std::min(idxi, jidx);
                            second_local_elements[count] = std::max(idxi, jidx);
                            count += 1;
                        }
                    }

                    noffset += chunk_len;
                }

                assert(expected_count == count);
            }
        });
}

#endif
