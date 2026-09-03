#ifndef SCCD_CELL2D_BROADPHASE_HPP
#define SCCD_CELL2D_BROADPHASE_HPP

#include "broadphase.hpp"
#include "sparallel.hpp"
#include "vaabb.hpp"

#include <cmath>
#include <cstring>
#include <vector>

/**
 * \file
 * \brief Broad phase over a two-dimensional cell list, with no sorting.
 *
 * The sweep-and-prune this replaces sorts three AABB lists (O(n log n)) and then
 * prunes on one axis. Measured on a refined cloth-ball frame, that walks about
 * 2,100 candidates for every pair it finds, and no choice of axis helps -- the
 * two good axes give 2,882 and 2,838, the bad one 126,063. One-dimensional
 * pruning cannot separate a draped surface.
 *
 * This bins instead. Two points of design worth stating because both were
 * deliberate:
 *
 * **Two axes, not three.** The geometry is surfaces, so with N elements of size h
 * in a box of side L, N ~ (L/h)^2. A 2D grid at that resolution has ~N cells and
 * is mostly occupied; a 3D grid has (L/h)^3 = N^1.5 -- at N = 1.5M that is 1.8
 * billion cells at 0.08% occupancy, memory and empty-cell iteration bought for
 * pruning along an axis the surface barely occupies. The third axis is rejected
 * by the ordinary AABB test inside the cell, so pairs are still culled in 3D.
 *
 * **No ordering assumption anywhere.** Cells are filled by counting sort and the
 * query tests every entry of a touched cell. There is no break-on-first-miss, so
 * nothing needs the lists sorted, which is the point: it removes the sort rather
 * than moving it.
 *
 * Both the binning and the pair collection are count-then-fill, matching the rest
 * of the broad phase: the allocation is exact and known before anything is
 * written, and the CRS output is deterministic.
 */

namespace sccd {

    /** \brief A uniform 2D cell list over two of the three coordinate axes. */
    template <typename T>
    struct Cell2DGrid {
        int axis0 = 0;
        int axis1 = 1;
        T min0 = 0, min1 = 0;
        T inv0 = 1, inv1 = 1;
        int n0 = 1, n1 = 1;

        ptrdiff_t ncells() const { return (ptrdiff_t)n0 * (ptrdiff_t)n1; }

        int clamp0(const T v) const {
            const T f = (v - min0) * inv0;
            if (!(f > T(0))) return 0;  // also catches NaN
            const int c = (int)f;
            return c >= n0 ? n0 - 1 : c;
        }

        int clamp1(const T v) const {
            const T f = (v - min1) * inv1;
            if (!(f > T(0))) return 0;
            const int c = (int)f;
            return c >= n1 ? n1 - 1 : c;
        }

        ptrdiff_t cell_of(const int c0, const int c1) const { return (ptrdiff_t)c1 * n0 + c0; }
    };

    namespace cell2d_detail {

        /** \brief The two axes with the largest centre spread, widest first. */
        template <typename T>
        static void choose_two_axes(const ptrdiff_t n, T** const SCCD_RESTRICT aabb, int& a0, int& a1) {
            T var[3];
            sccd::center_variance(n, aabb, var);

            int order[3] = {0, 1, 2};
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 3; ++j) {
                    if (var[order[j]] > var[order[i]]) {
                        const int t = order[i];
                        order[i] = order[j];
                        order[j] = t;
                    }
                }
            }
            a0 = order[0];
            a1 = order[1];
        }

    }  // namespace cell2d_detail

    /**
     * \brief Size a grid for \p n boxes so a box spans O(1) cells.
     *
     * The cell is the mean AABB extent, which keeps occupancy near one box per
     * cell and the number of cells near n. Sweep AABBs are much larger than the
     * elements they came from, so sizing on the element size instead would make
     * every box touch a large block of cells and give the counting pass more work
     * than the query saves.
     */
    template <typename T>
    static void cell2d_setup(const ptrdiff_t n, T** const SCCD_RESTRICT aabb, Cell2DGrid<T>& grid) {
        cell2d_detail::choose_two_axes<T>(n, aabb, grid.axis0, grid.axis1);

        const int d0 = grid.axis0;
        const int d1 = grid.axis1;

        T lo0 = aabb[d0][0], hi0 = aabb[3 + d0][0];
        T lo1 = aabb[d1][0], hi1 = aabb[3 + d1][0];
        T sum0 = 0, sum1 = 0;
        for (ptrdiff_t i = 0; i < n; ++i) {
            lo0 = sccd::min<T>(lo0, aabb[d0][i]);
            hi0 = sccd::max<T>(hi0, aabb[3 + d0][i]);
            lo1 = sccd::min<T>(lo1, aabb[d1][i]);
            hi1 = sccd::max<T>(hi1, aabb[3 + d1][i]);
            sum0 += aabb[3 + d0][i] - aabb[d0][i];
            sum1 += aabb[3 + d1][i] - aabb[d1][i];
        }

        const T span0 = sccd::max<T>(hi0 - lo0, std::numeric_limits<T>::min());
        const T span1 = sccd::max<T>(hi1 - lo1, std::numeric_limits<T>::min());
        const T mean0 = sccd::max<T>(sum0 / (T)n, span0 / (T)(1 << 20));
        const T mean1 = sccd::max<T>(sum1 / (T)n, span1 / (T)(1 << 20));

        // Cap the cell count so a degenerate extent cannot allocate unboundedly.
        // 4n cells is already past the point where more resolution pays.
        const double cap = 4.0 * (double)sccd::max<ptrdiff_t>(n, 1);
        double want0 = (double)(span0 / mean0);
        double want1 = (double)(span1 / mean1);
        if (want0 < 1) want0 = 1;
        if (want1 < 1) want1 = 1;
        const double total = want0 * want1;
        if (total > cap) {
            const double s = std::sqrt(cap / total);
            want0 *= s;
            want1 *= s;
            if (want0 < 1) want0 = 1;
            if (want1 < 1) want1 = 1;
        }

        grid.n0 = (int)sccd::max<double>(1.0, std::floor(want0));
        grid.n1 = (int)sccd::max<double>(1.0, std::floor(want1));
        grid.min0 = lo0;
        grid.min1 = lo1;
        // Nudge the span so the largest coordinate lands inside the last cell.
        grid.inv0 = (T)grid.n0 / (span0 * (T)(1 + std::numeric_limits<T>::epsilon()));
        grid.inv1 = (T)grid.n1 / (span1 * (T)(1 + std::numeric_limits<T>::epsilon()));
    }

    /**
     * \brief Bin \p n boxes into \p grid: count, prefix sum, then scatter.
     *
     * \p cellptr must hold ncells + 1 entries, \p cellidx the total span count
     * which is cellptr[ncells] after the prefix sum, so this is called in two
     * steps by the caller the same way the pair collection is.
     */
    template <typename T, typename I>
    static void cell2d_count(const ptrdiff_t n,
                             T** const SCCD_RESTRICT aabb,
                             const Cell2DGrid<T>& grid,
                             ptrdiff_t* const SCCD_RESTRICT cellptr) {
        const ptrdiff_t ncells = grid.ncells();
        std::memset(cellptr, 0, sizeof(ptrdiff_t) * (size_t)(ncells + 1));

        const T* const SCCD_RESTRICT lo0 = aabb[grid.axis0];
        const T* const SCCD_RESTRICT hi0 = aabb[3 + grid.axis0];
        const T* const SCCD_RESTRICT lo1 = aabb[grid.axis1];
        const T* const SCCD_RESTRICT hi1 = aabb[3 + grid.axis1];

        // Serial: the scatter is a histogram, and at ~n cells the contention of a
        // parallel version costs more than the pass itself.
        for (ptrdiff_t i = 0; i < n; ++i) {
            const int a = grid.clamp0(lo0[i]);
            const int b = grid.clamp0(hi0[i]);
            const int c = grid.clamp1(lo1[i]);
            const int d = grid.clamp1(hi1[i]);
            for (int j = c; j <= d; ++j) {
                for (int k = a; k <= b; ++k) {
                    cellptr[grid.cell_of(k, j) + 1] += 1;
                }
            }
        }

        for (ptrdiff_t i = 0; i < ncells; ++i) {
            cellptr[i + 1] += cellptr[i];
        }
    }

    /** \brief Scatter box indices into the cells counted by cell2d_count. */
    template <typename T, typename I>
    static void cell2d_fill(const ptrdiff_t n,
                            T** const SCCD_RESTRICT aabb,
                            const Cell2DGrid<T>& grid,
                            const ptrdiff_t* const SCCD_RESTRICT cellptr,
                            I* const SCCD_RESTRICT cellidx,
                            ptrdiff_t* const SCCD_RESTRICT cursor) {
        const ptrdiff_t ncells = grid.ncells();
        std::memcpy(cursor, cellptr, sizeof(ptrdiff_t) * (size_t)ncells);

        const T* const SCCD_RESTRICT lo0 = aabb[grid.axis0];
        const T* const SCCD_RESTRICT hi0 = aabb[3 + grid.axis0];
        const T* const SCCD_RESTRICT lo1 = aabb[grid.axis1];
        const T* const SCCD_RESTRICT hi1 = aabb[3 + grid.axis1];

        for (ptrdiff_t i = 0; i < n; ++i) {
            const int a = grid.clamp0(lo0[i]);
            const int b = grid.clamp0(hi0[i]);
            const int c = grid.clamp1(lo1[i]);
            const int d = grid.clamp1(hi1[i]);
            for (int j = c; j <= d; ++j) {
                for (int k = a; k <= b; ++k) {
                    cellidx[cursor[grid.cell_of(k, j)]++] = (I)i;
                }
            }
        }
    }

    namespace cell2d_detail {

        /**
         * \brief Walk the cells a box touches, reporting each overlapping partner once.
         *
         * A box spans several cells, so a pair can be met in several of them. The
         * pair is attributed to the cell holding the minimum corner of the two
         * boxes' overlap: that corner lies inside both boxes, so both are binned
         * there and the cell is unique. This is exact and costs two clamps, where
         * a hash or a mark array would cost memory proportional to the pair count.
         */
        template <int F, int S, typename T, typename I, typename Visit>
        static inline void for_each_unique_partner(T** const SCCD_RESTRICT first_aabbs,
                                                   const ptrdiff_t fi,
                                                   T** const SCCD_RESTRICT second_aabbs,
                                                   const I* const SCCD_RESTRICT second_idx,
                                                   I** const SCCD_RESTRICT second_elements,
                                                   const ptrdiff_t second_stride,
                                                   const I (&ev)[F],
                                                   const Cell2DGrid<T>& grid,
                                                   const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                   const I* const SCCD_RESTRICT cellidx,
                                                   Visit&& visit) {
            const T aminx = first_aabbs[0][fi], aminy = first_aabbs[1][fi], aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi], amaxy = first_aabbs[4][fi], amaxz = first_aabbs[5][fi];

            const T amin0 = first_aabbs[grid.axis0][fi];
            const T amax0 = first_aabbs[3 + grid.axis0][fi];
            const T amin1 = first_aabbs[grid.axis1][fi];
            const T amax1 = first_aabbs[3 + grid.axis1][fi];

            const int c0b = grid.clamp0(amin0), c0e = grid.clamp0(amax0);
            const int c1b = grid.clamp1(amin1), c1e = grid.clamp1(amax1);

            for (int c1 = c1b; c1 <= c1e; ++c1) {
                for (int c0 = c0b; c0 <= c0e; ++c0) {
                    const ptrdiff_t cell = grid.cell_of(c0, c1);
                    const ptrdiff_t begin = cellptr[cell];
                    const ptrdiff_t end = cellptr[cell + 1];

                    for (ptrdiff_t k = begin; k < end; ++k) {
                        const ptrdiff_t j = (ptrdiff_t)cellidx[k];

                        if (sccd::disjoint<T>(aminx,
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

                        // Report only from the cell owning the overlap's min corner.
                        const T o0 = sccd::max<T>(amin0, second_aabbs[grid.axis0][j]);
                        const T o1 = sccd::max<T>(amin1, second_aabbs[grid.axis1][j]);
                        if (grid.clamp0(o0) != c0 || grid.clamp1(o1) != c1) {
                            continue;
                        }

                        bool share = false;
                        if constexpr (S > 1) {
                            const I jidx = second_idx[j];
                            I sev[S];
                            for (int v = 0; v < S; ++v) {
                                sev[v] = second_elements[v][jidx * second_stride];
                            }
                            share = sccd::sccd_detail::shares_vertex<F, S>(ev, sev);
                        } else {
                            for (int a = 0; a < F; ++a) {
                                if (ev[a] == second_idx[j]) {
                                    share = true;
                                    break;
                                }
                            }
                        }
                        if (share) {
                            continue;
                        }

                        visit(j);
                    }
                }
            }
        }

        /**
         * \brief The self-overlap form: one list against itself.
         *
         * Two filters, and both are needed. `j > fi` makes each unordered pair
         * appear once rather than twice; the canonical cell makes it appear once
         * rather than once per cell the two boxes share. The sweep gets the first
         * for free by starting its window at fi + 1 in sorted order, which is not
         * available here because nothing is sorted.
         */
        template <int NXE, typename T, typename I, typename Visit>
        static inline void for_each_unique_self_partner(T** const SCCD_RESTRICT aabbs,
                                                        const ptrdiff_t fi,
                                                        I* const SCCD_RESTRICT idx,
                                                        I** const SCCD_RESTRICT elements,
                                                        const ptrdiff_t stride,
                                                        const I (&ev)[NXE],
                                                        const Cell2DGrid<T>& grid,
                                                        const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                        const I* const SCCD_RESTRICT cellidx,
                                                        Visit&& visit) {
            const T aminx = aabbs[0][fi], aminy = aabbs[1][fi], aminz = aabbs[2][fi];
            const T amaxx = aabbs[3][fi], amaxy = aabbs[4][fi], amaxz = aabbs[5][fi];

            const T amin0 = aabbs[grid.axis0][fi];
            const T amax0 = aabbs[3 + grid.axis0][fi];
            const T amin1 = aabbs[grid.axis1][fi];
            const T amax1 = aabbs[3 + grid.axis1][fi];

            const int c0b = grid.clamp0(amin0), c0e = grid.clamp0(amax0);
            const int c1b = grid.clamp1(amin1), c1e = grid.clamp1(amax1);

            for (int c1 = c1b; c1 <= c1e; ++c1) {
                for (int c0 = c0b; c0 <= c0e; ++c0) {
                    const ptrdiff_t cell = grid.cell_of(c0, c1);
                    const ptrdiff_t begin = cellptr[cell];
                    const ptrdiff_t end = cellptr[cell + 1];

                    for (ptrdiff_t k = begin; k < end; ++k) {
                        const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                        if (j <= fi) {
                            continue;
                        }

                        if (sccd::disjoint<T>(aminx,
                                              aminy,
                                              aminz,
                                              amaxx,
                                              amaxy,
                                              amaxz,
                                              aabbs[0][j],
                                              aabbs[1][j],
                                              aabbs[2][j],
                                              aabbs[3][j],
                                              aabbs[4][j],
                                              aabbs[5][j])) {
                            continue;
                        }

                        const T o0 = sccd::max<T>(amin0, aabbs[grid.axis0][j]);
                        const T o1 = sccd::max<T>(amin1, aabbs[grid.axis1][j]);
                        if (grid.clamp0(o0) != c0 || grid.clamp1(o1) != c1) {
                            continue;
                        }

                        const I jidx = idx[j];
                        I sev[NXE];
                        for (int v = 0; v < NXE; ++v) {
                            sev[v] = elements[v][jidx * stride];
                        }
                        if (sccd::sccd_detail::shares_vertex<NXE, NXE>(ev, sev)) {
                            continue;
                        }

                        visit(j, jidx);
                    }
                }
            }
        }

    }  // namespace cell2d_detail

    /** \brief Self-overlap count: one list against itself, CRS offsets out. */
    template <int nxe, typename T, typename I>
    bool cell2d_count_self_overlaps(const ptrdiff_t element_count,
                                    T** const SCCD_RESTRICT aabbs,
                                    I* const SCCD_RESTRICT idx,
                                    const ptrdiff_t stride,
                                    I** const SCCD_RESTRICT elements,
                                    const Cell2DGrid<T>& grid,
                                    const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                    const I* const SCCD_RESTRICT cellidx,
                                    ptrdiff_t* const SCCD_RESTRICT ccdptr) {
        ccdptr[0] = 0;

        sccd::parallel_for_br(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I idxi = idx[fi];
                I ev[nxe];
                for (int v = 0; v < nxe; ++v) {
                    ev[v] = elements[v][idxi * stride];
                }

                ptrdiff_t count = 0;
                cell2d_detail::for_each_unique_self_partner<nxe, T, I>(
                    aabbs, fi, idx, elements, stride, ev, grid, cellptr, cellidx,
                    [&](const ptrdiff_t, const I) { ++count; });
                ccdptr[fi + 1] = count;
            }
        });

        sccd::parallel_cum_sum_br(ccdptr, ccdptr + element_count + 1);
        return ccdptr[element_count] > 0;
    }

    /** \brief Write the self-overlap pairs, as (min, max) to match the sweep. */
    template <int nxe, typename T, typename I>
    void cell2d_fill_self_overlaps(const ptrdiff_t element_count,
                                   T** const SCCD_RESTRICT aabbs,
                                   I* const SCCD_RESTRICT idx,
                                   const ptrdiff_t stride,
                                   I** const SCCD_RESTRICT elements,
                                   const Cell2DGrid<T>& grid,
                                   const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                   const I* const SCCD_RESTRICT cellidx,
                                   const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                   I* const SCCD_RESTRICT out0,
                                   I* const SCCD_RESTRICT out1) {
        sccd::parallel_for_br(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I idxi = idx[fi];
                I ev[nxe];
                for (int v = 0; v < nxe; ++v) {
                    ev[v] = elements[v][idxi * stride];
                }

                ptrdiff_t at = ccdptr[fi];
                cell2d_detail::for_each_unique_self_partner<nxe, T, I>(
                    aabbs, fi, idx, elements, stride, ev, grid, cellptr, cellidx,
                    [&](const ptrdiff_t, const I jidx) {
                        out0[at] = sccd::min<I>(idxi, jidx);
                        out1[at] = sccd::max<I>(idxi, jidx);
                        ++at;
                    });
            }
        });
    }

    /** \brief Count candidate pairs per first-list box, then prefix-sum into CRS offsets. */
    template <int first_nxe, int second_nxe, typename T, typename I>
    bool cell2d_count_overlaps(const ptrdiff_t first_count,
                               T** const SCCD_RESTRICT first_aabbs,
                               I* const SCCD_RESTRICT first_idx,
                               const ptrdiff_t first_stride,
                               I** const SCCD_RESTRICT first_elements,
                               T** const SCCD_RESTRICT second_aabbs,
                               I* const SCCD_RESTRICT second_idx,
                               const ptrdiff_t second_stride,
                               I** const SCCD_RESTRICT second_elements,
                               const Cell2DGrid<T>& grid,
                               const ptrdiff_t* const SCCD_RESTRICT cellptr,
                               const I* const SCCD_RESTRICT cellidx,
                               ptrdiff_t* const SCCD_RESTRICT ccdptr) {
        ccdptr[0] = 0;

        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I first_idxi = first_idx[fi];
                I ev[first_nxe];
                for (int v = 0; v < first_nxe; ++v) {
                    ev[v] = first_elements[v][first_idxi * first_stride];
                }

                ptrdiff_t count = 0;
                cell2d_detail::for_each_unique_partner<first_nxe, second_nxe, T, I>(first_aabbs,
                                                                                    fi,
                                                                                    second_aabbs,
                                                                                    second_idx,
                                                                                    second_elements,
                                                                                    second_stride,
                                                                                    ev,
                                                                                    grid,
                                                                                    cellptr,
                                                                                    cellidx,
                                                                                    [&](const ptrdiff_t) { ++count; });
                ccdptr[fi + 1] = count;
            }
        });

        sccd::parallel_cum_sum_br(ccdptr, ccdptr + first_count + 1);
        return ccdptr[first_count] > 0;
    }

    /** \brief Write the pairs counted by cell2d_count_overlaps into the CRS arrays. */
    template <int first_nxe, int second_nxe, typename T, typename I>
    void cell2d_fill_overlaps(const ptrdiff_t first_count,
                              T** const SCCD_RESTRICT first_aabbs,
                              I* const SCCD_RESTRICT first_idx,
                              const ptrdiff_t first_stride,
                              I** const SCCD_RESTRICT first_elements,
                              T** const SCCD_RESTRICT second_aabbs,
                              I* const SCCD_RESTRICT second_idx,
                              const ptrdiff_t second_stride,
                              I** const SCCD_RESTRICT second_elements,
                              const Cell2DGrid<T>& grid,
                              const ptrdiff_t* const SCCD_RESTRICT cellptr,
                              const I* const SCCD_RESTRICT cellidx,
                              const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                              I* const SCCD_RESTRICT first_out,
                              I* const SCCD_RESTRICT second_out) {
        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I first_idxi = first_idx[fi];
                I ev[first_nxe];
                for (int v = 0; v < first_nxe; ++v) {
                    ev[v] = first_elements[v][first_idxi * first_stride];
                }

                ptrdiff_t at = ccdptr[fi];
                cell2d_detail::for_each_unique_partner<first_nxe, second_nxe, T, I>(
                    first_aabbs,
                    fi,
                    second_aabbs,
                    second_idx,
                    second_elements,
                    second_stride,
                    ev,
                    grid,
                    cellptr,
                    cellidx,
                    [&](const ptrdiff_t j) {
                        first_out[at] = first_idxi;
                        second_out[at] = second_idx[j];
                        ++at;
                    });
            }
        });
    }

}  // namespace sccd

#endif  // SCCD_CELL2D_BROADPHASE_HPP
