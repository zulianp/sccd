#ifndef LEAN_CCD_HPP
#define LEAN_CCD_HPP

/**
 * \file sccd.hpp
 * \brief SIMD-friendly AABB broad-phase helpers and lean CCD utilities (test
 * scope).
 *
 * Provides counting and collection of candidate overlaps for two lists and
 * self-overlaps, with vectorized disjoint tests and scalar fallbacks.
 * The AABB layout is structure-of-arrays: [minx, miny, minz, maxx, maxy, maxz].
 */

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <algorithm>
#include <cfloat>
#include <cstdio>
#include <filesystem>
#include <vector>

#include "vaabb.h"

namespace sccd {
/// Geometry scalar type used for coordinates and AABB values.
using geom_t = float;
/// Integer type used for element/vertex indices.
using idx_t = int;
/// Integer type used for counts and prefix sums.
using count_t = int;

#ifndef _WIN32
#define SFEM_RESTRICT __restrict__
#else
#define SFEM_RESTRICT __restrict
#endif

/**
 * \brief Return the next representable value toward +infinity.
 * \param x Input value.
 * \return Next representable float greater than or equal to x.
 */
inline geom_t nextafter_up(const geom_t x) { return nextafterf(x, FLT_MAX); }

/**
 * \brief Return the next representable value toward -infinity.
 * \param x Input value.
 * \return Next representable float less than or equal to x.
 */
inline geom_t nextafter_down(const geom_t x) { return nextafterf(x, -FLT_MAX); }

// inline geom_t nextafter_up(const geom_t x)
// {
//     return nextafter(x,DBL_MAX);
// }

// inline geom_t nextafter_down(const geom_t x)
// {
//     return nextafter(x, DBL_MAX);
// }

/**
 * \brief Choose the axis (0=x,1=y,2=z) with largest variance of AABB centers.
 * \param n Number of AABBs.
 * \param aabb SoA arrays of size 6: minx,miny,minz,maxx,maxy,maxz; each of
 * length n. \return Axis index in {0,1,2}.
 */
int choose_axis(const size_t n, geom_t **const SFEM_RESTRICT aabb) {
  geom_t mean[3] = {0};
  geom_t var[3] = {0};
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

/**
 * \brief Sort AABBs along \p sort_axis and permute all six SoA arrays
 * coherently. \param n Number of AABBs. \param sort_axis Axis to sort by
 * (0=x,1=y,2=z). \param arrays SoA arrays [6][n]:
 * minx,miny,minz,maxx,maxy,maxz. \param idx Output permutation array of size n
 * (initialized to 0..n-1 then sorted). \param scratch Scratch buffer of length
 * n used to permute each component.
 *
 * Arrays must be valid and sufficiently sized. Uses parallel sort.
 */
void sort_along_axis(const size_t n, const int sort_axis,
                     geom_t **const SFEM_RESTRICT arrays,
                     idx_t *const SFEM_RESTRICT idx,
                     geom_t *const SFEM_RESTRICT scratch) {
  for (size_t i = 0; i < n; i++) {
    idx[i] = i;
  }

  const geom_t *const SFEM_RESTRICT x = arrays[sort_axis];
  tbb::parallel_sort(idx, idx + n,
                     [x](const idx_t l, const idx_t r) { return x[l] < x[r]; });

  for (int d = 0; d < 6; d++) {
    memcpy(scratch, arrays[d], sizeof(geom_t) * n);
    for (size_t i = 0; i < n; i++) {
      arrays[d][i] = scratch[idx[i]];
    }
  }
}

/**
 * \brief Remap indices in-place through a permutation table.
 * \param n Number of entries.
 * \param idx Permutation table mapping old index -> new index.
 * \param remapped Array of indices to update; each entry is replaced by
 * idx[entry].
 */
static void remap_idx(const size_t n, const idx_t *const SFEM_RESTRICT idx,
                      idx_t *const SFEM_RESTRICT remapped) {
  tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
                    [&](const tbb::blocked_range<size_t> &r) {
                      for (size_t i = r.begin(); i < r.end(); i++) {
                        remapped[i] = idx[remapped[i]];
                      }
                    });
}

/**
 * \namespace sccd_detail
 * \brief Internal helpers for vectorized AABB disjoint tests and overlap
 * filtering.
 */
namespace sccd_detail {

/**
 * \brief Load the \p nxe vertex indices of element \p elem_idx.
 * \tparam nxe Number of indices per element.
 * \param elements SoA arrays for element indices.
 * \param elem_idx Logical element index.
 * \param stride Stride between consecutive elements in the arrays.
 * \param out Output array of size nxe.
 */
template <int nxe>
static inline void load_ev(idx_t **const SFEM_RESTRICT elements,
                           const idx_t elem_idx, const size_t stride,
                           idx_t (&out)[nxe]) {
  for (int v = 0; v < nxe; ++v) {
    out[v] = elements[v][elem_idx * stride];
  }
}

/**
 * \brief Test if two index tuples share any vertex id.
 * \tparam n1 Length of first tuple.
 * \tparam n2 Length of second tuple.
 * \param a First tuple of indices.
 * \param b Second tuple of indices.
 * \return True if any index matches.
 */
template <int n1, int n2>
static inline bool shares_vertex(const idx_t (&a)[n1], const idx_t (&b)[n2]) {
  for (int i = 0; i < n1; ++i) {
    for (int j = 0; j < n2; ++j) {
      if (a[i] == b[j]) {
        return true;
      }
    }
  }
  return false;
}

/**
 * \brief Advance begin and compute end of candidate window along a sorted axis.
 * \param fimin Min coordinate of first AABB on the sort axis.
 * \param fimax Max coordinate of first AABB on the sort axis.
 * \param second_xmax Sorted xmax array of the second list along the sort axis.
 * \param second_xmin Sorted xmin array of the second list along the sort axis.
 * \param second_count Number of AABBs in the second list.
 * \param begin In/out: advanced to first index where second_xmax[begin] >
 * fimin. \param end Out: first index where second_xmin[end] > fimax.
 */
static inline void compute_candidate_window_progressive(
    const geom_t fimin, const geom_t fimax,
    const geom_t *const SFEM_RESTRICT second_xmax,
    const geom_t *const SFEM_RESTRICT second_xmin, const size_t second_count,
    size_t &begin, size_t &end) {
  for (; begin < second_count; ++begin) {
    if (fimin < second_xmax[begin]) {
      break;
    }
  }
  end = begin;
  for (; end < second_count; ++end) {
    if (fimax < second_xmin[end]) {
      break;
    }
  }
}

/**
 * \brief Load a contiguous block of B AABBs into SoA buffers.
 * \param aabbs SoA arrays [6][...].
 * \param start Starting index in B.
 * \param len Number of AABBs to load.
 * \param B_minx..B_maxz Output arrays of length at least \p len.
 */
static inline void prepare_B_block(
    geom_t **const SFEM_RESTRICT aabbs, const size_t start, const size_t len,
    geom_t *const SFEM_RESTRICT B_minx, geom_t *const SFEM_RESTRICT B_miny,
    geom_t *const SFEM_RESTRICT B_minz, geom_t *const SFEM_RESTRICT B_maxx,
    geom_t *const SFEM_RESTRICT B_maxy, geom_t *const SFEM_RESTRICT B_maxz) {
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

/**
 * \brief Force remaining SIMD lanes (len..chunk_size) to be disjoint.
 * \param amaxx0,amaxy0,amaxz0 Max components of A used to place B outside.
 * \param len Number of valid lanes already filled [0..len).
 * \param B_minx..B_maxz In/out buffers to tail-fill.
 */
static inline void tail_fill_B(const geom_t amaxx0, const geom_t amaxy0,
                               const geom_t amaxz0, const size_t len,
                               geom_t *const SFEM_RESTRICT B_minx,
                               geom_t *const SFEM_RESTRICT B_miny,
                               geom_t *const SFEM_RESTRICT B_minz,
                               geom_t *const SFEM_RESTRICT B_maxx,
                               geom_t *const SFEM_RESTRICT B_maxy,
                               geom_t *const SFEM_RESTRICT B_maxz) {
  for (size_t lane = len; lane < AABB_DISJOINT_CHUNK_SIZE; ++lane) {
    B_minx[lane] = amaxx0 + 1;
    B_miny[lane] = amaxy0 + 1;
    B_minz[lane] = amaxz0 + 1;
    B_maxx[lane] = amaxx0;
    B_maxy[lane] = amaxy0;
    B_maxz[lane] = amaxz0;
  }
}

/**
 * \brief Build per-lane disjoint mask for a block of B against broadcast A.
 * \param second_aabbs SoA arrays for the second list.
 * \param start First B index to test.
 * \param chunk_len Number of lanes to process (<= AABB_DISJOINT_CHUNK_SIZE).
 * \param A_minx..A_maxz Broadcast A components.
 * \param amaxx0,amaxy0,amaxz0 Max components of A used for tail fill.
 * \param mask_out Output mask: 1=disjoint, 0=potential overlap.
 */
static inline void build_disjoint_mask_for_block(
    geom_t **const SFEM_RESTRICT second_aabbs, const size_t start,
    const size_t chunk_len, const geom_t *const SFEM_RESTRICT A_minx,
    const geom_t *const SFEM_RESTRICT A_miny,
    const geom_t *const SFEM_RESTRICT A_minz,
    const geom_t *const SFEM_RESTRICT A_maxx,
    const geom_t *const SFEM_RESTRICT A_maxy,
    const geom_t *const SFEM_RESTRICT A_maxz, const geom_t amaxx0,
    const geom_t amaxy0, const geom_t amaxz0,
    uint32_t *const SFEM_RESTRICT mask_out) {
  alignas(64) geom_t B_minx[AABB_DISJOINT_CHUNK_SIZE];
  alignas(64) geom_t B_miny[AABB_DISJOINT_CHUNK_SIZE];
  alignas(64) geom_t B_minz[AABB_DISJOINT_CHUNK_SIZE];
  alignas(64) geom_t B_maxx[AABB_DISJOINT_CHUNK_SIZE];
  alignas(64) geom_t B_maxy[AABB_DISJOINT_CHUNK_SIZE];
  alignas(64) geom_t B_maxz[AABB_DISJOINT_CHUNK_SIZE];

  prepare_B_block(second_aabbs, start, chunk_len, B_minx, B_miny, B_minz,
                  B_maxx, B_maxy, B_maxz);
  tail_fill_B(amaxx0, amaxy0, amaxz0, chunk_len, B_minx, B_miny, B_minz, B_maxx,
              B_maxy, B_maxz);
  vdisjoint(A_minx, A_miny, A_minz, A_maxx, A_maxy, A_maxz, B_minx, B_miny,
            B_minz, B_maxx, B_maxy, B_maxz, mask_out);
}

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
 * \param second_stride Stride between elements in second arrays.
 */
template <int F, int S>
static inline void mask_out_shared_two_lists(
    uint32_t *const SFEM_RESTRICT dmask, const size_t chunk_len,
    const size_t noffset, const idx_t (&ev)[F],
    const idx_t *const SFEM_RESTRICT second_idx,
    idx_t **const SFEM_RESTRICT second_elements, const size_t second_stride) {
  for (size_t lane = 0; lane < chunk_len; ++lane) {
    if (dmask[lane])
      continue;

    const size_t j = noffset + lane;
    const idx_t jidx = second_idx[j];
    int match = 0;
    if constexpr (S > 1) {
      idx_t sev[S];
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

/**
 * \brief Mark lanes where elements share a vertex in self-overlap path.
 * \tparam N Number of vertices per element.
 * \param dmask In/out lane mask; set to 1 when a vertex is shared.
 * \param chunk_len Number of valid lanes.
 * \param noffset Starting j index (j > i).
 * \param ev Vertex indices of element i.
 * \param idx Mapping from sorted position to element id.
 * \param elements SoA vertex arrays.
 * \param stride Stride between elements in the arrays.
 */
template <int N>
static inline void
mask_out_shared_self(uint32_t *const SFEM_RESTRICT dmask,
                     const size_t chunk_len, const size_t noffset,
                     const idx_t (&ev)[N], const idx_t *const SFEM_RESTRICT idx,
                     idx_t **const SFEM_RESTRICT elements,
                     const size_t stride) {
  for (size_t lane = 0; lane < chunk_len; ++lane) {
    if (dmask[lane])
      continue;

    const size_t j = noffset + lane;
    const idx_t jidx = idx[j];
    idx_t sev[N];
    load_ev<N>(elements, jidx, stride, sev);
    if (shares_vertex<N, N>(ev, sev)) {
      dmask[lane] = 1;
    }
  }
}

// -----------------------------

/**
 * \brief Scalar reference: count candidate overlaps in [begin,end) for two
 * lists. \return Number of non-disjoint, non-shared-vertex candidates.
 */
template <int F, int S>
static inline size_t scalar_count_range_two_lists(
    geom_t **const SFEM_RESTRICT first_aabbs, const size_t fi,
    geom_t **const SFEM_RESTRICT second_aabbs,
    const idx_t *const SFEM_RESTRICT second_idx,
    idx_t **const SFEM_RESTRICT second_elements, const size_t second_stride,
    const idx_t (&ev)[F], const size_t begin, const size_t end) {
  size_t count = 0;
  const geom_t aminx = first_aabbs[0][fi];
  const geom_t aminy = first_aabbs[1][fi];
  const geom_t aminz = first_aabbs[2][fi];
  const geom_t amaxx = first_aabbs[3][fi];
  const geom_t amaxy = first_aabbs[4][fi];
  const geom_t amaxz = first_aabbs[5][fi];
  for (size_t j = begin; j < end; ++j) {
    if (disjoint(aminx, aminy, aminz, amaxx, amaxy, amaxz, second_aabbs[0][j],
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

/**
 * \brief Scalar reference: collect candidate overlaps in [begin,end) for two
 * lists. \return Number of pairs written to \p first_out and \p second_out.
 */
template <int F, int S>
static inline size_t scalar_collect_range_two_lists(
    geom_t **const SFEM_RESTRICT first_aabbs, const size_t fi,
    const idx_t first_idxi, geom_t **const SFEM_RESTRICT second_aabbs,
    const idx_t *const SFEM_RESTRICT second_idx,
    idx_t **const SFEM_RESTRICT second_elements, const size_t second_stride,
    const idx_t (&ev)[F], const size_t begin, const size_t end,
    idx_t *const SFEM_RESTRICT first_out,
    idx_t *const SFEM_RESTRICT second_out) {
  size_t count = 0;
  const geom_t aminx = first_aabbs[0][fi];
  const geom_t aminy = first_aabbs[1][fi];
  const geom_t aminz = first_aabbs[2][fi];
  const geom_t amaxx = first_aabbs[3][fi];
  const geom_t amaxy = first_aabbs[4][fi];
  const geom_t amaxz = first_aabbs[5][fi];
  for (size_t j = begin; j < end; ++j) {
    if (disjoint(aminx, aminy, aminz, amaxx, amaxy, amaxz, second_aabbs[0][j],
                 second_aabbs[1][j], second_aabbs[2][j], second_aabbs[3][j],
                 second_aabbs[4][j], second_aabbs[5][j])) {
      continue;
    }
    const idx_t jidx = second_idx[j];
    int match = 0;
    if constexpr (S > 1) {
      idx_t sev[S];
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

// -----------------------------

/**
 * \brief Scalar reference: count self-overlaps in [begin,end) for element i.
 * \return Number of non-disjoint, non-shared-vertex candidates with j>i.
 */
template <int N>
static inline size_t
scalar_count_range_self(geom_t **const SFEM_RESTRICT aabbs, const size_t fi,
                        idx_t **const SFEM_RESTRICT elements,
                        const idx_t *const SFEM_RESTRICT idx,
                        const size_t stride, const idx_t (&ev)[N],
                        const size_t begin, const size_t end) {
  size_t count = 0;
  const geom_t aminx = aabbs[0][fi];
  const geom_t aminy = aabbs[1][fi];
  const geom_t aminz = aabbs[2][fi];
  const geom_t amaxx = aabbs[3][fi];
  const geom_t amaxy = aabbs[4][fi];
  const geom_t amaxz = aabbs[5][fi];
  for (size_t j = begin; j < end; ++j) {
    if (disjoint(aminx, aminy, aminz, amaxx, amaxy, amaxz, aabbs[0][j],
                 aabbs[1][j], aabbs[2][j], aabbs[3][j], aabbs[4][j],
                 aabbs[5][j])) {
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

/**
 * \brief Scalar reference: collect self-overlaps in [begin,end) for element i.
 * \return Number of pairs written to outputs, with (min(idxi,jidx), max(...)).
 */
template <int N>
static inline size_t scalar_collect_range_self(
    geom_t **const SFEM_RESTRICT aabbs, const size_t fi, const idx_t idxi,
    idx_t **const SFEM_RESTRICT elements, const idx_t *const SFEM_RESTRICT idx,
    const size_t stride, const idx_t (&ev)[N], const size_t begin,
    const size_t end, idx_t *const SFEM_RESTRICT first_out,
    idx_t *const SFEM_RESTRICT second_out) {
  size_t count = 0;
  const geom_t aminx = aabbs[0][fi];
  const geom_t aminy = aabbs[1][fi];
  const geom_t aminz = aabbs[2][fi];
  const geom_t amaxx = aabbs[3][fi];
  const geom_t amaxy = aabbs[4][fi];
  const geom_t amaxz = aabbs[5][fi];
  for (size_t j = begin; j < end; ++j) {
    if (disjoint(aminx, aminy, aminz, amaxx, amaxy, amaxz, aabbs[0][j],
                 aabbs[1][j], aabbs[2][j], aabbs[3][j], aabbs[4][j],
                 aabbs[5][j])) {
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

/**
 * \brief Count candidate overlaps between two sorted AABB lists.
 *
 * Excludes pairs that are axis-disjoint or share a vertex. Both AABB lists
 * must be sorted by \p sort_axis and their index arrays aligned accordingly.
 *
 * \tparam first_nxe Number of vertices per element in the first list.
 * \tparam second_nxe Number of vertices per element in the second list.
 * \param sort_axis Sort axis (0=x,1=y,2=z).
 * \param first_count Number of AABBs in the first list.
 * \param first_aabbs SoA arrays [6][first_count].
 * \param first_idx Mapping from sorted position to element id for the first
 * list. \param first_stride Stride between elements in the first element
 * arrays. \param first_elements SoA element-vertex arrays for the first list.
 * \param second_count Number of AABBs in the second list.
 * \param second_aabbs SoA arrays [6][second_count].
 * \param second_idx Mapping from sorted position to element id for the second
 * list. \param second_stride Stride between elements in the second element
 * arrays. \param second_elements SoA element-vertex arrays for the second list.
 * \param ccdptr Prefix sum array of size first_count+1. On return:
 *               ccdptr[i+1]-ccdptr[i] = candidates for first i, and
 *               ccdptr[first_count] = total candidates.
 * \return True if any candidates exist.
 */
template <int first_nxe, int second_nxe>
bool count_overlaps(const int sort_axis, const count_t first_count,
                    geom_t **const SFEM_RESTRICT first_aabbs,
                    idx_t *const SFEM_RESTRICT first_idx,
                    const size_t first_stride,
                    idx_t **const SFEM_RESTRICT first_elements,
                    const count_t second_count,
                    geom_t **const SFEM_RESTRICT second_aabbs,
                    idx_t *const SFEM_RESTRICT second_idx,
                    const size_t second_stride,
                    idx_t **const SFEM_RESTRICT second_elements,
                    size_t *const SFEM_RESTRICT ccdptr) {
  const geom_t *const SFEM_RESTRICT first_xmin = first_aabbs[sort_axis];
  const geom_t *const SFEM_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
  const geom_t *const SFEM_RESTRICT second_xmin = second_aabbs[sort_axis];
  const geom_t *const SFEM_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

  ccdptr[0] = 0;

  if (first_xmax[first_count - 1] < second_xmin[0])
    return false;

  if (second_xmax[second_count - 1] < first_xmin[0])
    return false;

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, first_count),
      [&](const tbb::blocked_range<size_t> &r) {
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

          idx_t ev[first_nxe];
          for (int v = 0; v < first_nxe; v++) {
            ev[v] = first_elements[v][first_idxi * first_stride];
          }

          size_t end = ni;
          sccd_detail::compute_candidate_window_progressive(
              fimin, fimax, second_xmax, second_xmin, second_count, ni, end);

          if (ni >= end) {
            continue;
          }

          if (end - ni < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
            size_t count =
                sccd_detail::scalar_count_range_two_lists<first_nxe,
                                                          second_nxe>(
                    first_aabbs, fi, second_aabbs, second_idx, second_elements,
                    second_stride, ev, ni, end);
            ccdptr[fi + 1] = count;
            continue;
          }

          size_t count = 0;
          size_t noffset = ni;

          alignas(64) geom_t A_minx[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_miny[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_minz[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxx[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxy[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxz[AABB_DISJOINT_CHUNK_SIZE];

          vaabb_broadcast(first_aabbs, fi, A_minx, A_miny, A_minz, A_maxx,
                          A_maxy, A_maxz);

          for (; noffset < end;) {
            const size_t chunk_len =
                std::min((size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

            uint32_t dmask[AABB_DISJOINT_CHUNK_SIZE] = {0};

            sccd_detail::build_disjoint_mask_for_block(
                second_aabbs, noffset, chunk_len, A_minx, A_miny, A_minz,
                A_maxx, A_maxy, A_maxz, A_maxx[0], A_maxy[0], A_maxz[0], dmask);

            sccd_detail::mask_out_shared_two_lists<first_nxe, second_nxe>(
                dmask, chunk_len, noffset, ev, second_idx, second_elements,
                second_stride);

            for (size_t lane = 0; lane < chunk_len; ++lane) {
              count += dmask[lane] ? 0 : 1;
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

/**
 * \brief Collect candidate overlaps between two sorted AABB lists.
 *
 * Uses the prefix offsets computed by count_overlaps to write pairs.
 *
 * \param sort_axis Sort axis (0=x,1=y,2=z).
 * \param first_count Number of AABBs in the first list.
 * \param first_aabbs SoA arrays [6][first_count].
 * \param first_idx Mapping from sorted position to element id for the first
 * list. \param first_stride Stride between elements in the first element
 * arrays. \param first_elements SoA element-vertex arrays for the first list.
 * \param second_count Number of AABBs in the second list.
 * \param second_aabbs SoA arrays [6][second_count].
 * \param second_idx Mapping from sorted position to element id for the second
 * list. \param second_stride Stride between elements in the second element
 * arrays. \param second_elements SoA element-vertex arrays for the second list.
 * \param ccdptr Prefix offsets from the count pass (size first_count+1).
 * \param foverlap Output array (size ccdptr[first_count]) for first indices.
 * \param noverlap Output array (size ccdptr[first_count]) for second indices.
 */
template <int first_nxe, int second_nxe>
void collect_overlaps(
    const int sort_axis, const count_t first_count,
    geom_t **const SFEM_RESTRICT first_aabbs,
    idx_t *const SFEM_RESTRICT first_idx, const size_t first_stride,
    idx_t **SFEM_RESTRICT const first_elements, const count_t second_count,
    geom_t **const SFEM_RESTRICT second_aabbs,
    idx_t *const SFEM_RESTRICT second_idx, const size_t second_stride,
    idx_t **SFEM_RESTRICT const second_elements,
    const size_t *const SFEM_RESTRICT ccdptr, idx_t *SFEM_RESTRICT foverlap,
    idx_t *SFEM_RESTRICT noverlap) {
  const geom_t *const SFEM_RESTRICT first_xmin = first_aabbs[sort_axis];
  const geom_t *const SFEM_RESTRICT first_xmax = first_aabbs[3 + sort_axis];
  const geom_t *const SFEM_RESTRICT second_xmin = second_aabbs[sort_axis];
  const geom_t *const SFEM_RESTRICT second_xmax = second_aabbs[3 + sort_axis];

  if (first_xmax[first_count - 1] < second_xmin[0])
    return;

  if (second_xmax[second_count - 1] < first_xmin[0])
    return;

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, first_count),
      [&](const tbb::blocked_range<size_t> &r) {
        size_t ni = 0;
        for (; ni < second_count; ni++) {
          if (second_xmax[ni] > first_xmin[0]) {
            break;
          }
        }

        for (size_t fi = r.begin(); fi < r.end(); fi++) {
          const size_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
          if (expected_count == 0) {
            continue;
          }

          const geom_t fimin = first_xmin[fi];
          const geom_t fimax = first_xmax[fi];
          const idx_t first_idxi = first_idx[fi];

          idx_t *SFEM_RESTRICT const first_local_elements =
              &foverlap[ccdptr[fi]];
          idx_t *SFEM_RESTRICT const second_local_elements =
              &noverlap[ccdptr[fi]];

          idx_t ev[first_nxe];
          for (int v = 0; v < first_nxe; v++) {
            ev[v] = first_elements[v][first_idxi * first_stride];
          }

          size_t end = ni;
          sccd_detail::compute_candidate_window_progressive(
              fimin, fimax, second_xmax, second_xmin, second_count, ni, end);

          if (ni >= end) {
            continue;
          }

          if (end - ni < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
            size_t count =
                sccd_detail::scalar_collect_range_two_lists<first_nxe,
                                                            second_nxe>(
                    first_aabbs, fi, first_idxi, second_aabbs, second_idx,
                    second_elements, second_stride, ev, ni, end,
                    first_local_elements, second_local_elements);
            assert(expected_count == count);
            continue;
          }

          size_t count = 0;
          size_t noffset = ni;

          alignas(64) geom_t A_minx[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_miny[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_minz[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxx[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxy[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxz[AABB_DISJOINT_CHUNK_SIZE];

          vaabb_broadcast(first_aabbs, fi, A_minx, A_miny, A_minz, A_maxx,
                          A_maxy, A_maxz);

          for (; noffset < end;) {
            const size_t chunk_len =
                std::min((size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

            uint32_t dmask[AABB_DISJOINT_CHUNK_SIZE] = {0};
            sccd_detail::build_disjoint_mask_for_block(
                second_aabbs, noffset, chunk_len, A_minx, A_miny, A_minz,
                A_maxx, A_maxy, A_maxz, A_maxx[0], A_maxy[0], A_maxz[0], dmask);

            sccd_detail::mask_out_shared_two_lists<first_nxe, second_nxe>(
                dmask, chunk_len, noffset, ev, second_idx, second_elements,
                second_stride);

            for (size_t lane = 0; lane < chunk_len; ++lane) {
              if (dmask[lane])
                continue;
              const size_t j = noffset + lane;
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

// --------------------------------------

/**
 * \brief Count candidate self-overlaps (upper triangle) within one sorted AABB
 * list.
 *
 * Excludes pairs that are axis-disjoint or share a vertex. The AABBs must be
 * sorted by \p sort_axis and \p idx maps sorted positions to element ids.
 *
 * \tparam nxe Number of vertices per element.
 * \param sort_axis Sort axis (0=x,1=y,2=z).
 * \param element_count Number of AABBs/elements.
 * \param aabbs SoA arrays [6][element_count].
 * \param idx Mapping from sorted position to element id.
 * \param stride Stride between elements in the vertex arrays.
 * \param elements SoA element-vertex arrays.
 * \param ccdptr Prefix sum array size element_count+1; filled as in the
 * two-lists case. \return True if any candidates exist.
 */
template <int nxe>
bool count_self_overlaps(const int sort_axis, const count_t element_count,
                         geom_t **const SFEM_RESTRICT aabbs,
                         idx_t *const SFEM_RESTRICT idx, const size_t stride,
                         idx_t **const SFEM_RESTRICT elements,
                         size_t *const SFEM_RESTRICT ccdptr) {
  const geom_t *const SFEM_RESTRICT xmin = aabbs[sort_axis];
  const geom_t *const SFEM_RESTRICT xmax = aabbs[3 + sort_axis];

  ccdptr[0] = 0;

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, element_count),
      [&](const tbb::blocked_range<size_t> &r) {
        for (size_t fi = r.begin(); fi < r.end(); fi++) {
          const geom_t fimin = xmin[fi];
          const geom_t fimax = xmax[fi];
          const idx_t idxi = idx[fi];

          assert(idxi >= 0);
          assert(idxi < element_count);

          idx_t ev[nxe];
          for (int v = 0; v < nxe; v++) {
            ev[v] = elements[v][idxi * stride];
          }

          size_t noffset = fi + 1;
          size_t end = noffset;
          sccd_detail::compute_candidate_window_progressive(
              fimin, fimax, xmax, xmin, element_count, noffset, end);

          if (noffset >= end) {
            ccdptr[fi + 1] = 0;
            continue;
          }

          if (end - noffset < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
            size_t count = sccd_detail::scalar_count_range_self<nxe>(
                aabbs, fi, elements, idx, stride, ev, noffset, end);
            ccdptr[fi + 1] = count;
            continue;
          }

          size_t count = 0;

          alignas(64) geom_t A_minx[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_miny[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_minz[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxx[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxy[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxz[AABB_DISJOINT_CHUNK_SIZE];
          vaabb_broadcast(aabbs, fi, A_minx, A_miny, A_minz, A_maxx, A_maxy,
                          A_maxz);

          for (; noffset < end;) {
            const size_t chunk_len =
                std::min((size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

            uint32_t mask[AABB_DISJOINT_CHUNK_SIZE] = {0};
            sccd_detail::build_disjoint_mask_for_block(
                aabbs, noffset, chunk_len, A_minx, A_miny, A_minz, A_maxx,
                A_maxy, A_maxz, A_maxx[0], A_maxy[0], A_maxz[0], mask);

            sccd_detail::mask_out_shared_self<nxe>(mask, chunk_len, noffset, ev,
                                                   idx, elements, stride);

            for (size_t lane = 0; lane < chunk_len; ++lane) {
              count += mask[lane] ? 0 : 1;
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

/**
 * \brief Collect candidate self-overlap pairs using prefix offsets.
 *
 * Writes pairs (min(idxi,jidx), max(idxi,jidx)) for j>i into output arrays.
 *
 * \param sort_axis Sort axis (0=x,1=y,2=z).
 * \param element_count Number of elements/AABBs.
 * \param aabbs SoA arrays [6][element_count].
 * \param idx Mapping from sorted position to element id.
 * \param stride Stride between elements in the vertex arrays.
 * \param elements SoA element-vertex arrays.
 * \param ccdptr Prefix offsets from the self count pass (size element_count+1).
 * \param foverlap Output array of first indices.
 * \param noverlap Output array of second indices.
 */
template <int nxe>
void collect_self_overlaps(const int sort_axis, const count_t element_count,
                           geom_t **const SFEM_RESTRICT aabbs,
                           idx_t *const SFEM_RESTRICT idx, const size_t stride,
                           idx_t **const elements,
                           const size_t *const SFEM_RESTRICT ccdptr,
                           idx_t *SFEM_RESTRICT foverlap,
                           idx_t *SFEM_RESTRICT noverlap) {
  const geom_t *const SFEM_RESTRICT xmin = aabbs[sort_axis];
  const geom_t *const SFEM_RESTRICT xmax = aabbs[3 + sort_axis];

  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, element_count),
      [&](const tbb::blocked_range<size_t> &r) {
        for (size_t fi = r.begin(); fi < r.end(); fi++) {
          const size_t expected_count = ccdptr[fi + 1] - ccdptr[fi];
          if (!expected_count)
            continue;

          const geom_t fimin = xmin[fi];
          const geom_t fimax = xmax[fi];
          const idx_t idxi = idx[fi];

          idx_t *SFEM_RESTRICT const first_local_elements =
              &foverlap[ccdptr[fi]];

          idx_t *SFEM_RESTRICT const second_local_elements =
              &noverlap[ccdptr[fi]];

          idx_t ev[nxe];
          for (int v = 0; v < nxe; v++) {
            ev[v] = elements[v][idxi * stride];
          }

          size_t noffset = fi + 1;
          size_t end = noffset;
          sccd_detail::compute_candidate_window_progressive(
              fimin, fimax, xmax, xmin, element_count, noffset, end);

          if (noffset >= end) {
            continue;
          }

          if (end - noffset < AABB_DISJOINT_NOVECTORIZE_THRESHOLD) {
            size_t count = sccd_detail::scalar_collect_range_self<nxe>(
                aabbs, fi, idxi, elements, idx, stride, ev, noffset, end,
                first_local_elements, second_local_elements);
            assert(expected_count == count);
            continue;
          }

          size_t count = 0;

          alignas(64) geom_t A_minx[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_miny[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_minz[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxx[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxy[AABB_DISJOINT_CHUNK_SIZE];
          alignas(64) geom_t A_maxz[AABB_DISJOINT_CHUNK_SIZE];

          vaabb_broadcast(aabbs, fi, A_minx, A_miny, A_minz, A_maxx, A_maxy,
                          A_maxz);

          for (; noffset < end;) {
            const size_t chunk_len =
                std::min((size_t)AABB_DISJOINT_CHUNK_SIZE, end - noffset);

            uint32_t mask[AABB_DISJOINT_CHUNK_SIZE] = {0};

            sccd_detail::build_disjoint_mask_for_block(
                aabbs, noffset, chunk_len, A_minx, A_miny, A_minz, A_maxx,
                A_maxy, A_maxz, A_maxx[0], A_maxy[0], A_maxz[0], mask);

            sccd_detail::mask_out_shared_self<nxe>(mask, chunk_len, noffset, ev,
                                                   idx, elements, stride);

            for (size_t lane = 0; lane < chunk_len; ++lane) {
              if (mask[lane])
                continue;
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

} // namespace sccd
#endif
