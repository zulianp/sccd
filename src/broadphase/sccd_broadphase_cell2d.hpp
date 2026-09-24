#ifndef SCCD_BROADPHASE_CELL2D_HPP
#define SCCD_BROADPHASE_CELL2D_HPP

#include "sccd_broadphase_sweep.hpp"
#include "sccd_parallel.hpp"
#include "sccd_aabb.hpp"

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
 * This bins instead. Three points of design worth stating because all three were
 * deliberate:
 *
 * **The centroid decides the cell, and the cell is at least a box wide.** A box
 * goes in one cell, the one holding its centroid, and a cell is as wide as the
 * widest box the grid will see. Two overlapping boxes then have centroids at most
 * one maximum extent apart on each axis, so their cells differ by at most one and
 * a fixed 3x3 stencil finds every partner. The search radius is a property of the
 * grid rather than of the data, which is what the cell size buys.
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

    namespace detail {

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

        /** \brief The midpoint of a box on one axis, written so it cannot overflow. */
        template <typename T>
        static inline T midpoint(const T lo, const T hi) {
            return lo + (hi - lo) * T(0.5);
        }

        /** \brief The column of the cell holding box \p i's centroid. */
        template <typename T>
        static inline int centroid_col(T** const SCCD_RESTRICT aabb, const Cell2DGrid<T>& grid, const ptrdiff_t i) {
            return grid.clamp0(midpoint<T>(aabb[grid.axis0][i], aabb[3 + grid.axis0][i]));
        }

        /** \brief The row of the cell holding box \p i's centroid. */
        template <typename T>
        static inline int centroid_row(T** const SCCD_RESTRICT aabb, const Cell2DGrid<T>& grid, const ptrdiff_t i) {
            return grid.clamp1(midpoint<T>(aabb[grid.axis1][i], aabb[3 + grid.axis1][i]));
        }

    }  // namespace detail

    /** \brief The largest AABB extent on each axis, for sizing a grid. */
    template <typename T>
    static void max_box_extent(const ptrdiff_t n, T** const SCCD_RESTRICT aabb, T out[3]) {
        struct Ext3 {
            T e[3];
        };

        if (n <= 0) {
            out[0] = out[1] = out[2] = T(0);
            return;
        }

        const Ext3 m = sccd::parallel_tiled_reduce<Ext3>(
            0,
            n,
            [&](const ptrdiff_t begin, const ptrdiff_t end) {
                Ext3 e{{T(0), T(0), T(0)}};
                for (ptrdiff_t i = begin; i < end; ++i) {
                    for (int d = 0; d < 3; ++d) {
                        e.e[d] = sccd::max<T>(e.e[d], aabb[3 + d][i] - aabb[d][i]);
                    }
                }
                return e;
            },
            [](const Ext3 a, const Ext3 b) {
                return Ext3{{sccd::max<T>(a.e[0], b.e[0]),
                             sccd::max<T>(a.e[1], b.e[1]),
                             sccd::max<T>(a.e[2], b.e[2])}};
            });

        for (int d = 0; d < 3; ++d) out[d] = m.e[d];
    }

    /**
     * \brief Shrink a grid until it holds at most `4 * n` cells.
     *
     * The caller sizes its cell array from this bound, so it has to hold
     * exactly, not approximately. Both sides come in at least one, and the
     * larger is reduced to whatever the smaller leaves room for.
     */
    inline void cap_cells(const ptrdiff_t n, int& n0, int& n1) {
        const ptrdiff_t cap = sccd::max<ptrdiff_t>(4 * n, 1);
        if (n0 < 1) n0 = 1;
        if (n1 < 1) n1 = 1;
        if ((ptrdiff_t)n0 > cap) n0 = (int)cap;
        if ((ptrdiff_t)n1 > cap) n1 = (int)cap;
        if ((ptrdiff_t)n0 * (ptrdiff_t)n1 <= cap) return;
        if (n0 >= n1) {
            n0 = (int)sccd::max<ptrdiff_t>(cap / (ptrdiff_t)n1, 1);
        } else {
            n1 = (int)sccd::max<ptrdiff_t>(cap / (ptrdiff_t)n0, 1);
        }
    }

    /**
     * \brief Size a grid over \p n boxes so that a cell is at least a box wide.
     *
     * The cell width is the largest AABB extent on the axis, which is what makes
     * the query a fixed 3x3 stencil: two overlapping boxes have centroids no more
     * than one maximum extent apart on each axis, so their centroid cells differ
     * by at most one. Every subsequent shrink of the grid widens the cells, so the
     * bound survives the cell cap.
     *
     * \p query_extent is the per-axis largest extent of the list that will query
     * this grid, for the two-list form where the querying boxes are not the boxes
     * in the cells. A face is wider than a vertex, so the face-vertex query needs
     * the face width here for its stencil to be complete. Null means the list
     * queries itself.
     */
    template <typename T>
    static void cell2d_setup(const ptrdiff_t n,
                             T** const SCCD_RESTRICT aabb,
                             Cell2DGrid<T>& grid,
                             const T* const query_extent = nullptr) {
        detail::choose_two_axes<T>(n, aabb, grid.axis0, grid.axis1);

        const int d0 = grid.axis0;
        const int d1 = grid.axis1;

        struct Extent {
            T lo0, hi0, lo1, hi1, ext0, ext1;
        };

        const Extent ext = sccd::parallel_tiled_reduce<Extent>(
            0,
            n,
            [&](const ptrdiff_t begin, const ptrdiff_t end) {
                Extent e{aabb[d0][begin], aabb[3 + d0][begin], aabb[d1][begin], aabb[3 + d1][begin], 0, 0};
                for (ptrdiff_t i = begin; i < end; ++i) {
                    e.lo0 = sccd::min<T>(e.lo0, aabb[d0][i]);
                    e.hi0 = sccd::max<T>(e.hi0, aabb[3 + d0][i]);
                    e.lo1 = sccd::min<T>(e.lo1, aabb[d1][i]);
                    e.hi1 = sccd::max<T>(e.hi1, aabb[3 + d1][i]);
                    e.ext0 = sccd::max<T>(e.ext0, aabb[3 + d0][i] - aabb[d0][i]);
                    e.ext1 = sccd::max<T>(e.ext1, aabb[3 + d1][i] - aabb[d1][i]);
                }
                return e;
            },
            [](const Extent a, const Extent b) {
                return Extent{sccd::min<T>(a.lo0, b.lo0),
                              sccd::max<T>(a.hi0, b.hi0),
                              sccd::min<T>(a.lo1, b.lo1),
                              sccd::max<T>(a.hi1, b.hi1),
                              sccd::max<T>(a.ext0, b.ext0),
                              sccd::max<T>(a.ext1, b.ext1)};
            });

        const T lo0 = ext.lo0, hi0 = ext.hi0;
        const T lo1 = ext.lo1, hi1 = ext.hi1;

        const T span0 = sccd::max<T>(hi0 - lo0, std::numeric_limits<T>::min());
        const T span1 = sccd::max<T>(hi1 - lo1, std::numeric_limits<T>::min());

        // The widest box either list can put in a cell, and never finer than the
        // span divided by a million so a scene of points still gets a grid.
        const T q0 = query_extent ? query_extent[d0] : T(0);
        const T q1 = query_extent ? query_extent[d1] : T(0);
        const T cell0 = sccd::max<T>(sccd::max<T>(ext.ext0, q0), span0 / (T)(1 << 20));
        const T cell1 = sccd::max<T>(sccd::max<T>(ext.ext1, q1), span1 / (T)(1 << 20));

        // Cap the cell count so a degenerate extent cannot allocate unboundedly.
        // 4n cells is already past the point where more resolution pays.
        const double cap = 4.0 * (double)sccd::max<ptrdiff_t>(n, 1);
        double want0 = (double)(span0 / cell0);
        double want1 = (double)(span1 / cell1);
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
        // Scaling both sides by the same factor does not enforce the cap on a
        // scene whose two axes differ wildly: a side scaled below one is raised
        // back to one, and the product grows past the cap again by exactly that
        // much. The bound is what the caller sizes its cell array from, so it is
        // enforced here on the integers that array is sized by.
        cap_cells(sccd::max<ptrdiff_t>(n, 1), grid.n0, grid.n1);
        grid.min0 = lo0;
        grid.min1 = lo1;
        // Nudge the span so the largest coordinate lands inside the last cell.
        grid.inv0 = (T)grid.n0 / (span0 * (T)(1 + std::numeric_limits<T>::epsilon()));
        grid.inv1 = (T)grid.n1 / (span1 * (T)(1 + std::numeric_limits<T>::epsilon()));
    }

    /**
     * \brief Which boxes a block of cell rows has to look at.
     *
     * The histogram and the scatter both write a counter per box, so neither is
     * safe to run box-parallel without either atomics or one private copy of the
     * whole cell array per worker. Neither is affordable: the cell count is of the
     * order of the box count, so private copies cost workers times that.
     *
     * Partitioning by cell row avoids both. A block owns a contiguous
     * range of rows on the grid's second axis, and no two blocks own a cell, so
     * each can write its own slice of the cell array with no synchronisation at
     * all. A box is binned by its centroid, so it lands in exactly one block's
     * list and the lists partition the boxes.
     *
     * Boxes stay in index order within a block's list, and blocks cover the rows
     * in order, so both passes visit each cell's boxes in exactly the order the
     * serial form does. The cell array and the index array that come out are
     * therefore identical to the serial ones, not merely equivalent.
     */
    struct Cell2DPartition {
        int nblocks = 1;
        int rows_per_block = 1;
        std::vector<ptrdiff_t> blockptr;  ///< nblocks + 1 offsets into blockbox
        std::vector<int> blockbox;        ///< box indices, grouped by block
        std::vector<ptrdiff_t> counts;    ///< per-chunk, per-block counters, reused across steps

        /** \brief True when the two binning passes should just run serially. */
        bool serial() const { return nblocks <= 1; }
    };

    /**
     * \brief Below this many boxes the partition costs more than it saves.
     *
     * Measured on a draped sheet at ten workers, the partitioned binning is
     * $0.42\times$ the speed of the serial one at $8{,}000$ boxes, level at
     * $16{,}000$, $1.61\times$ at $32{,}000$ and $3.10\times$ at $256{,}000$. The
     * threshold sits just past the crossover: a mesh small enough to be below it
     * has a binning pass measured in tenths of a millisecond either way.
     *
     * A total, not a count per block. Deriving the block count from the boxes
     * instead -- so that a wider machine asks for fewer, fatter blocks -- was
     * four times slower on armadillo-rollers at seventy-two threads, because it
     * capped the partition at a handful of blocks on exactly the mesh sizes
     * preparation spends its time on.
     */
    static const ptrdiff_t SCCD_CELL2D_MIN_PARALLEL = 32768;

    /**
     * \brief Group \p n boxes by the block of cell rows they touch.
     *
     * Cheap relative to what it enables: one pass over the boxes reading two
     * coordinates each to count, and a second to write. Both are box-parallel
     * and neither touches the cell array.
     */
    template <typename T>
    static void cell2d_partition(const ptrdiff_t n,
                                 T** const SCCD_RESTRICT aabb,
                                 const Cell2DGrid<T>& grid,
                                 Cell2DPartition& part) {
        part.blockptr.clear();
        part.blockbox.clear();

        const int max_workers = sccd::max_concurrency();
        if (n < SCCD_CELL2D_MIN_PARALLEL || max_workers <= 1 || grid.n1 <= 1) {
            part.nblocks = 1;
            part.rows_per_block = grid.n1;
            return;
        }

        // More blocks than workers so that a row band holding an unusually dense
        // patch of the surface does not become the whole critical path.
        const int want = sccd::min<int>(grid.n1, max_workers * 4);
        part.nblocks = sccd::max<int>(1, want);
        part.rows_per_block = (grid.n1 + part.nblocks - 1) / part.nblocks;
        part.nblocks = (grid.n1 + part.rows_per_block - 1) / part.rows_per_block;

        const int nblocks = part.nblocks;
        const int rpb = part.rows_per_block;

        const T* const SCCD_RESTRICT lo1 = aabb[grid.axis1];
        const T* const SCCD_RESTRICT hi1 = aabb[3 + grid.axis1];

        // Chunks of boxes, one private counter row each. The counters are over
        // blocks and not over cells, so this is nchunks * nblocks entries and
        // not nchunks * ncells.
        const int nchunks = nblocks;
        const ptrdiff_t chunk = (n + nchunks - 1) / nchunks;
        std::vector<ptrdiff_t> &counts = part.counts;
        counts.assign((size_t)nchunks * (size_t)nblocks, 0);

        sccd::parallel_for_chunks(0, nchunks, [&](const ptrdiff_t c) {
            ptrdiff_t* const row = counts.data() + c * (ptrdiff_t)nblocks;
            const ptrdiff_t begin = c * chunk;
            const ptrdiff_t end = sccd::min<ptrdiff_t>(begin + chunk, n);
            for (ptrdiff_t i = begin; i < end; ++i) {
                row[grid.clamp1(detail::midpoint<T>(lo1[i], hi1[i])) / rpb] += 1;
            }
        });

        // Offsets, block-major so a block's list is contiguous, chunk-minor so
        // the boxes inside it stay in index order.
        part.blockptr.assign((size_t)nblocks + 1, 0);
        ptrdiff_t running = 0;
        for (int b = 0; b < nblocks; ++b) {
            part.blockptr[(size_t)b] = running;
            for (int c = 0; c < nchunks; ++c) {
                ptrdiff_t& slot = counts[(size_t)c * (size_t)nblocks + (size_t)b];
                const ptrdiff_t take = slot;
                slot = running;
                running += take;
            }
        }
        part.blockptr[(size_t)nblocks] = running;

        part.blockbox.resize((size_t)running);
        sccd::parallel_for_chunks(0, nchunks, [&](const ptrdiff_t c) {
            ptrdiff_t* const row = counts.data() + c * (ptrdiff_t)nblocks;
            const ptrdiff_t begin = c * chunk;
            const ptrdiff_t end = sccd::min<ptrdiff_t>(begin + chunk, n);
            for (ptrdiff_t i = begin; i < end; ++i) {
                const int b = grid.clamp1(detail::midpoint<T>(lo1[i], hi1[i])) / rpb;
                part.blockbox[(size_t)row[b]++] = (int)i;
            }
        });
    }

    namespace detail {

        /**
         * \brief Visit the one cell box \p i belongs to, if a block of rows owns it.
         *
         * A box is binned by its centroid and so belongs to exactly one cell. The
         * query relies on that: the cell is at least as wide as the widest box, so
         * a partner's centroid is at most one cell away and the stencil is fixed.
         *
         * Shared by the counting and the scatter pass so that the two cannot
         * disagree about which cell a box belongs to, which is the failure that
         * would corrupt the cell array silently.
         */
        template <typename T, typename Visit>
        static inline void for_each_incidence(T** const SCCD_RESTRICT aabb,
                                              const Cell2DGrid<T>& grid,
                                              const ptrdiff_t i,
                                              const int row_begin,
                                              const int row_end,
                                              Visit&& visit) {
            const int c1 = centroid_row<T>(aabb, grid, i);
            if (c1 < row_begin || c1 >= row_end) return;
            visit(grid.cell_of(centroid_col<T>(aabb, grid, i), c1));
        }

    }  // namespace detail

    /**
     * \brief Bin \p n boxes into \p grid: count, prefix sum, then scatter.
     *
     * \p cellptr must hold ncells + 1 entries, \p cellidx the total span count
     * which is cellptr[ncells] after the prefix sum, so this is called in two
     * steps by the caller the same way the pair collection is.
     *
     * \p part comes from cell2d_partition on the same boxes and grid. It decides
     * only how the work is divided; the result does not depend on it.
     */
    template <typename T>
    static void cell2d_count(const ptrdiff_t n,
                             T** const SCCD_RESTRICT aabb,
                             const Cell2DGrid<T>& grid,
                             const Cell2DPartition& part,
                             ptrdiff_t* const SCCD_RESTRICT cellptr) {
        const ptrdiff_t ncells = grid.ncells();

        if (part.serial()) {
            std::memset(cellptr, 0, sizeof(ptrdiff_t) * (size_t)(ncells + 1));
            for (ptrdiff_t i = 0; i < n; ++i) {
                detail::for_each_incidence<T>(
                    aabb, grid, i, 0, grid.n1, [&](const ptrdiff_t cell) { cellptr[cell + 1] += 1; });
            }
            for (ptrdiff_t i = 0; i < ncells; ++i) {
                cellptr[i + 1] += cellptr[i];
            }
            return;
        }

        sccd::parallel_for_br(0, ncells + 1, [&](const ptrdiff_t begin, const ptrdiff_t end) {
            std::memset(cellptr + begin, 0, sizeof(ptrdiff_t) * (size_t)(end - begin));
        });

        const int rpb = part.rows_per_block;
        sccd::parallel_for_chunks(0, part.nblocks, [&](const ptrdiff_t b) {
            const int row_begin = (int)b * rpb;
            const int row_end = sccd::min<int>(row_begin + rpb, grid.n1);
            const ptrdiff_t from = part.blockptr[(size_t)b];
            const ptrdiff_t to = part.blockptr[(size_t)b + 1];
            for (ptrdiff_t e = from; e < to; ++e) {
                detail::for_each_incidence<T>(aabb,
                                              grid,
                                              (ptrdiff_t)part.blockbox[(size_t)e],
                                              row_begin,
                                              row_end,
                                              [&](const ptrdiff_t cell) { cellptr[cell + 1] += 1; });
            }
        });

        sccd::parallel_cum_sum_br(cellptr, cellptr + ncells + 1);
    }

    /** \brief Scatter box indices into the cells counted by cell2d_count. */
    template <typename T, typename I>
    static void cell2d_fill(const ptrdiff_t n,
                            T** const SCCD_RESTRICT aabb,
                            const Cell2DGrid<T>& grid,
                            const Cell2DPartition& part,
                            const ptrdiff_t* const SCCD_RESTRICT cellptr,
                            I* const SCCD_RESTRICT cellidx,
                            ptrdiff_t* const SCCD_RESTRICT cursor) {
        const ptrdiff_t ncells = grid.ncells();

        if (part.serial()) {
            std::memcpy(cursor, cellptr, sizeof(ptrdiff_t) * (size_t)ncells);
            for (ptrdiff_t i = 0; i < n; ++i) {
                detail::for_each_incidence<T>(aabb, grid, i, 0, grid.n1, [&](const ptrdiff_t cell) {
                    cellidx[cursor[cell]++] = (I)i;
                });
            }
            return;
        }

        sccd::parallel_for_br(0, ncells, [&](const ptrdiff_t begin, const ptrdiff_t end) {
            std::memcpy(cursor + begin, cellptr + begin, sizeof(ptrdiff_t) * (size_t)(end - begin));
        });

        const int rpb = part.rows_per_block;
        sccd::parallel_for_chunks(0, part.nblocks, [&](const ptrdiff_t b) {
            const int row_begin = (int)b * rpb;
            const int row_end = sccd::min<int>(row_begin + rpb, grid.n1);
            const ptrdiff_t from = part.blockptr[(size_t)b];
            const ptrdiff_t to = part.blockptr[(size_t)b + 1];
            for (ptrdiff_t e = from; e < to; ++e) {
                const ptrdiff_t i = (ptrdiff_t)part.blockbox[(size_t)e];
                detail::for_each_incidence<T>(aabb, grid, i, row_begin, row_end, [&](const ptrdiff_t cell) {
                    cellidx[cursor[cell]++] = (I)i;
                });
            }
        });
    }

    namespace detail {

        /**
         * \brief Walk the nine cells around a box, reporting every overlapping partner.
         *
         * The two lists are distinct, so every pair is reported by its first-list
         * box and once only; nothing has to be deduplicated. The nine cells are
         * enough because the grid was sized on the wider of the two lists, which
         * puts a partner's centroid at most one cell away on each axis.
         */
        template <int F, int S, typename T, typename I, typename Visit>
        static inline void for_each_unique_partner(T** const SCCD_RESTRICT first_aabbs,
                                                   const ptrdiff_t fi,
                                                   T** const SCCD_RESTRICT second_aabbs,
                                                   const I* const SCCD_RESTRICT second_idx,
                                                   I** const SCCD_RESTRICT second_elements,
                                                   const ptrdiff_t second_element_stride,
                                                   const I (&ev)[F],
                                                   const Cell2DGrid<T>& grid,
                                                   const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                   const I* const SCCD_RESTRICT cellidx,
                                                   Visit&& visit) {
            const T aminx = first_aabbs[0][fi], aminy = first_aabbs[1][fi], aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi], amaxy = first_aabbs[4][fi], amaxz = first_aabbs[5][fi];

            const int col = centroid_col<T>(first_aabbs, grid, fi);
            const int row = centroid_row<T>(first_aabbs, grid, fi);

            for (int dr = -1; dr <= 1; ++dr) {
                const int c1 = row + dr;
                if (c1 < 0 || c1 >= grid.n1) continue;

                for (int dc = -1; dc <= 1; ++dc) {
                    const int c0 = col + dc;
                    if (c0 < 0 || c0 >= grid.n0) continue;

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

                        bool share = false;
                        if constexpr (S > 1) {
                            const I jidx = second_idx[j];
                            I sev[S];
                            for (int v = 0; v < S; ++v) {
                                sev[v] = second_elements[v][jidx * second_element_stride];
                            }
                            share = sccd::detail::shares_vertex<F, S>(ev, sev);
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
         * Half the stencil, and one index test. The walk takes the row above and
         * its own row from its own column on, which is the half of the nine cells
         * whose linear index is at least its own. A pair in two different cells is
         * therefore reported by the box in the lower cell, and a pair inside one
         * cell by the lower index, so each unordered pair is reported once by
         * construction and the query reads five cells rather than nine.
         */
        template <int NXE, typename T, typename I, typename Visit>
        static inline void for_each_unique_self_partner(T** const SCCD_RESTRICT aabbs,
                                                        const ptrdiff_t fi,
                                                        I* const SCCD_RESTRICT idx,
                                                        I** const SCCD_RESTRICT elements,
                                                        const ptrdiff_t element_stride,
                                                        const I (&ev)[NXE],
                                                        const Cell2DGrid<T>& grid,
                                                        const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                        const I* const SCCD_RESTRICT cellidx,
                                                        Visit&& visit) {
            const T aminx = aabbs[0][fi], aminy = aabbs[1][fi], aminz = aabbs[2][fi];
            const T amaxx = aabbs[3][fi], amaxy = aabbs[4][fi], amaxz = aabbs[5][fi];

            const int col = centroid_col<T>(aabbs, grid, fi);
            const int row = centroid_row<T>(aabbs, grid, fi);

            for (int dr = 0; dr <= 1; ++dr) {
                const int c1 = row + dr;
                if (c1 >= grid.n1) continue;

                for (int dc = (dr == 0 ? 0 : -1); dc <= 1; ++dc) {
                    const int c0 = col + dc;
                    if (c0 < 0 || c0 >= grid.n0) continue;

                    const bool own = (dr == 0 && dc == 0);
                    const ptrdiff_t cell = grid.cell_of(c0, c1);
                    const ptrdiff_t begin = cellptr[cell];
                    const ptrdiff_t end = cellptr[cell + 1];

                    for (ptrdiff_t k = begin; k < end; ++k) {
                        const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                        if (own && j <= fi) {
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

                        const I jidx = idx[j];
                        I sev[NXE];
                        for (int v = 0; v < NXE; ++v) {
                            sev[v] = elements[v][jidx * element_stride];
                        }
                        if (sccd::detail::shares_vertex<NXE, NXE>(ev, sev)) {
                            continue;
                        }

                        visit(j, jidx);
                    }
                }
            }
        }

    }  // namespace detail

    /** \brief Self-overlap count: one list against itself, CRS offsets out. */
    template <int nxe, typename T, typename I>
    bool cell2d_count_self_overlaps(const ptrdiff_t element_count,
                                    T** const SCCD_RESTRICT aabbs,
                                    I* const SCCD_RESTRICT idx,
                                    const ptrdiff_t element_stride,
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
                    ev[v] = elements[v][idxi * element_stride];
                }

                ptrdiff_t count = 0;
                detail::for_each_unique_self_partner<nxe, T, I>(
                    aabbs, fi, idx, elements, element_stride, ev, grid, cellptr, cellidx,
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
                                   const ptrdiff_t element_stride,
                                   I** const SCCD_RESTRICT elements,
                                   const Cell2DGrid<T>& grid,
                                   const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                   const I* const SCCD_RESTRICT cellidx,
                                   const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                   I* const SCCD_RESTRICT first_out,
                                   I* const SCCD_RESTRICT second_out) {
        sccd::parallel_for_br(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I idxi = idx[fi];
                I ev[nxe];
                for (int v = 0; v < nxe; ++v) {
                    ev[v] = elements[v][idxi * element_stride];
                }

                ptrdiff_t at = ccdptr[fi];
                detail::for_each_unique_self_partner<nxe, T, I>(
                    aabbs, fi, idx, elements, element_stride, ev, grid, cellptr, cellidx,
                    [&](const ptrdiff_t, const I jidx) {
                        first_out[at] = sccd::min<I>(idxi, jidx);
                        second_out[at] = sccd::max<I>(idxi, jidx);
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
                               const ptrdiff_t first_element_stride,
                               I** const SCCD_RESTRICT first_elements,
                               T** const SCCD_RESTRICT second_aabbs,
                               I* const SCCD_RESTRICT second_idx,
                               const ptrdiff_t second_element_stride,
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
                    ev[v] = first_elements[v][first_idxi * first_element_stride];
                }

                ptrdiff_t count = 0;
                detail::for_each_unique_partner<first_nxe, second_nxe, T, I>(first_aabbs,
                                                                                    fi,
                                                                                    second_aabbs,
                                                                                    second_idx,
                                                                                    second_elements,
                                                                                    second_element_stride,
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
                              const ptrdiff_t first_element_stride,
                              I** const SCCD_RESTRICT first_elements,
                              T** const SCCD_RESTRICT second_aabbs,
                              I* const SCCD_RESTRICT second_idx,
                              const ptrdiff_t second_element_stride,
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
                    ev[v] = first_elements[v][first_idxi * first_element_stride];
                }

                ptrdiff_t at = ccdptr[fi];
                detail::for_each_unique_partner<first_nxe, second_nxe, T, I>(
                    first_aabbs,
                    fi,
                    second_aabbs,
                    second_idx,
                    second_elements,
                    second_element_stride,
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

#endif  // SCCD_BROADPHASE_CELL2D_HPP
