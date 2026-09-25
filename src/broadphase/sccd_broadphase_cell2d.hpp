#ifndef SCCD_BROADPHASE_CELL2D_HPP
#define SCCD_BROADPHASE_CELL2D_HPP

#include "sccd_broadphase_sweep.hpp"
#include "sccd_parallel.hpp"
#include "sccd_aabb.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
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
 *
 * ## A second binning, for the self query
 *
 * The one list against itself -- edge against edge -- can be done more cheaply,
 * and the `cell2dmin_` functions below do it. A box enters the **one** cell
 * holding its minimum corner instead of every cell it touches, so the cell array
 * holds one entry per box rather than one per covered cell, and a partner is met
 * at most once. The duplicate rule disappears with it.
 *
 * **The walk has to be over a total order.** Reading the box's own footprint --
 * its columns crossed with its rows -- loses every pair whose minimum-corner
 * cells are *incomparable*, one box ahead on the first axis and behind on the
 * second, because then neither footprint holds the other's cell and neither ever
 * reads the other. Componentwise order on cells is partial. Row-major linear
 * index is total, and over it the walk is complete: for overlapping boxes the
 * partner's cell never passes the cell of this box's maximum corner, so whichever
 * of the two comes first in that order holds the other in its range.
 *
 * **The bounds are what make it cheap.** Read literally, "forward in linear
 * order" means whole rows. Each cell therefore carries the largest upper bound of
 * the boxes binned in it, and each row a prefix maximum of those, so a binary
 * search skips the columns holding nothing that reaches back to the querying box
 * -- the sweep's running maximum, applied per row. What is left is close to the
 * footprint again, and the measured candidate count is below what the binning
 * above walks. A wide box is the case that separates them: it costs the extent
 * binning a place in every cell it crosses and costs this one a single entry.
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

        /** \brief The axis the grid does not use, which the cell test culls on. */
        int axis2() const { return 3 - axis0 - axis1; }
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

    }  // namespace detail

    /**
     * \brief Size a grid for \p n boxes so a box spans O(1) cells.
     *
     * The cell is the mean AABB extent, which keeps occupancy near one box per
     * cell and the number of cells near n. Sweep AABBs are much larger than the
     * elements they came from, so sizing on the element size instead would make
     * every box touch a large block of cells and give the counting pass more work
     * than the query saves.
     */
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

    template <typename T>
    static void cell2d_setup(const ptrdiff_t n, T** const SCCD_RESTRICT aabb, Cell2DGrid<T>& grid) {
        detail::choose_two_axes<T>(n, aabb, grid.axis0, grid.axis1);

        const int d0 = grid.axis0;
        const int d1 = grid.axis1;

        struct Extent {
            T lo0, hi0, lo1, hi1, sum0, sum1;
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
                    e.sum0 += aabb[3 + d0][i] - aabb[d0][i];
                    e.sum1 += aabb[3 + d1][i] - aabb[d1][i];
                }
                return e;
            },
            [](const Extent a, const Extent b) {
                return Extent{sccd::min<T>(a.lo0, b.lo0),
                              sccd::max<T>(a.hi0, b.hi0),
                              sccd::min<T>(a.lo1, b.lo1),
                              sccd::max<T>(a.hi1, b.hi1),
                              a.sum0 + b.sum0,
                              a.sum1 + b.sum1};
            });

        const T lo0 = ext.lo0, hi0 = ext.hi0;
        const T lo1 = ext.lo1, hi1 = ext.hi1;
        const T sum0 = ext.sum0, sum1 = ext.sum1;

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
     * The histogram and the scatter both write one counter per (box, cell)
     * incidence, so neither is safe to run box-parallel without either atomics
     * or one private copy of the whole cell array per worker. Neither is
     * affordable: the cell count is of the order of the box count, so private
     * copies cost workers times that.
     *
     * Partitioning by cell row avoids both. A block owns a contiguous
     * range of rows on the grid's second axis, and no two blocks own a cell, so
     * each can write its own slice of the cell array with no synchronisation at
     * all. A box that spans several row blocks appears in each of their lists,
     * which is why the lists are built by the same count-then-fill pass as
     * everything else here.
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
                const int b0 = grid.clamp1(lo1[i]) / rpb;
                const int b1 = grid.clamp1(hi1[i]) / rpb;
                for (int b = b0; b <= b1; ++b) {
                    row[b] += 1;
                }
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
                const int b0 = grid.clamp1(lo1[i]) / rpb;
                const int b1 = grid.clamp1(hi1[i]) / rpb;
                for (int b = b0; b <= b1; ++b) {
                    part.blockbox[(size_t)row[b]++] = (int)i;
                }
            }
        });
    }

    namespace detail {

        /**
         * \brief Visit every (cell, box) incidence a block of rows owns.
         *
         * Shared by the counting and the scatter pass so that the two cannot
         * disagree about which cells a box belongs to, which is the failure that
         * would corrupt the cell array silently.
         */
        template <typename T, typename Visit>
        static inline void for_each_incidence(T** const SCCD_RESTRICT aabb,
                                              const Cell2DGrid<T>& grid,
                                              const ptrdiff_t i,
                                              const int row_begin,
                                              const int row_end,
                                              Visit&& visit) {
            const int a = grid.clamp0(aabb[grid.axis0][i]);
            const int b = grid.clamp0(aabb[3 + grid.axis0][i]);
            const int c = sccd::max<int>(grid.clamp1(aabb[grid.axis1][i]), row_begin);
            const int d = sccd::min<int>(grid.clamp1(aabb[3 + grid.axis1][i]), row_end - 1);
            for (int j = c; j <= d; ++j) {
                for (int k = a; k <= b; ++k) {
                    visit(grid.cell_of(k, j));
                }
            }
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
                                                   const ptrdiff_t second_element_stride,
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
                                                        const ptrdiff_t element_stride,
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

    // ---------------------------------------------------------------------
    // The self query over a minimum-corner binning. See the \file block.
    // ---------------------------------------------------------------------

    namespace detail {

        /** \brief The column of the cell holding box \p i's minimum corner. */
        template <typename T>
        static inline int min_col(T** const SCCD_RESTRICT aabb, const Cell2DGrid<T>& grid, const ptrdiff_t i) {
            return grid.clamp0(aabb[grid.axis0][i]);
        }

        /** \brief The row of the cell holding box \p i's minimum corner. */
        template <typename T>
        static inline int min_row(T** const SCCD_RESTRICT aabb, const Cell2DGrid<T>& grid, const ptrdiff_t i) {
            return grid.clamp1(aabb[grid.axis1][i]);
        }

    }  // namespace detail

    /**
     * \brief Group \p n boxes by the block of cell rows their minimum corner is in.
     *
     * The same device as cell2d_partition and for the same reason -- the counting
     * and the scatter write per-cell counters, so they are split by cell row and
     * no two blocks touch a cell. The difference is that a box has one cell here,
     * so it appears in exactly one block's list and the lists *partition* the
     * boxes rather than covering them. That also makes the per-cell bounds below
     * safe to write from the same blocks.
     */
    template <typename T>
    static void cell2dmin_partition(const ptrdiff_t n,
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

        const int want = sccd::min<int>(grid.n1, max_workers * 4);
        part.nblocks = sccd::max<int>(1, want);
        part.rows_per_block = (grid.n1 + part.nblocks - 1) / part.nblocks;
        part.nblocks = (grid.n1 + part.rows_per_block - 1) / part.rows_per_block;

        const int nblocks = part.nblocks;
        const int rpb = part.rows_per_block;

        const int nchunks = nblocks;
        const ptrdiff_t chunk = (n + nchunks - 1) / nchunks;
        std::vector<ptrdiff_t>& counts = part.counts;
        counts.assign((size_t)nchunks * (size_t)nblocks, 0);

        sccd::parallel_for_chunks(0, nchunks, [&](const ptrdiff_t c) {
            ptrdiff_t* const row = counts.data() + c * (ptrdiff_t)nblocks;
            const ptrdiff_t begin = c * chunk;
            const ptrdiff_t end = sccd::min<ptrdiff_t>(begin + chunk, n);
            for (ptrdiff_t i = begin; i < end; ++i) {
                row[detail::min_row<T>(aabb, grid, i) / rpb] += 1;
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
                part.blockbox[(size_t)(row[detail::min_row<T>(aabb, grid, i) / rpb]++)] = (int)i;
            }
        });
    }

    /** \brief Count the boxes per cell by minimum corner, then prefix-sum into CRS offsets. */
    template <typename T>
    static void cell2dmin_count(const ptrdiff_t n,
                                T** const SCCD_RESTRICT aabb,
                                const Cell2DGrid<T>& grid,
                                const Cell2DPartition& part,
                                ptrdiff_t* const SCCD_RESTRICT cellptr) {
        const ptrdiff_t ncells = grid.ncells();

        if (part.serial()) {
            std::memset(cellptr, 0, sizeof(ptrdiff_t) * (size_t)(ncells + 1));
            for (ptrdiff_t i = 0; i < n; ++i) {
                cellptr[grid.cell_of(detail::min_col<T>(aabb, grid, i), detail::min_row<T>(aabb, grid, i)) + 1] += 1;
            }
            for (ptrdiff_t i = 0; i < ncells; ++i) {
                cellptr[i + 1] += cellptr[i];
            }
            return;
        }

        sccd::parallel_for_br(0, ncells + 1, [&](const ptrdiff_t begin, const ptrdiff_t end) {
            std::memset(cellptr + begin, 0, sizeof(ptrdiff_t) * (size_t)(end - begin));
        });

        sccd::parallel_for_chunks(0, part.nblocks, [&](const ptrdiff_t b) {
            const ptrdiff_t from = part.blockptr[(size_t)b];
            const ptrdiff_t to = part.blockptr[(size_t)b + 1];
            for (ptrdiff_t e = from; e < to; ++e) {
                const ptrdiff_t i = (ptrdiff_t)part.blockbox[(size_t)e];
                cellptr[grid.cell_of(detail::min_col<T>(aabb, grid, i), detail::min_row<T>(aabb, grid, i)) + 1] += 1;
            }
        });

        sccd::parallel_cum_sum_br(cellptr, cellptr + ncells + 1);
    }

    /** \brief Scatter each box index into the cell its minimum corner is in. */
    template <typename T, typename I>
    static void cell2dmin_fill(const ptrdiff_t n,
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
                const ptrdiff_t cell =
                    grid.cell_of(detail::min_col<T>(aabb, grid, i), detail::min_row<T>(aabb, grid, i));
                cellidx[cursor[cell]++] = (I)i;
            }
            return;
        }

        sccd::parallel_for_br(0, ncells, [&](const ptrdiff_t begin, const ptrdiff_t end) {
            std::memcpy(cursor + begin, cellptr + begin, sizeof(ptrdiff_t) * (size_t)(end - begin));
        });

        sccd::parallel_for_chunks(0, part.nblocks, [&](const ptrdiff_t b) {
            const ptrdiff_t from = part.blockptr[(size_t)b];
            const ptrdiff_t to = part.blockptr[(size_t)b + 1];
            for (ptrdiff_t e = from; e < to; ++e) {
                const ptrdiff_t i = (ptrdiff_t)part.blockbox[(size_t)e];
                const ptrdiff_t cell =
                    grid.cell_of(detail::min_col<T>(aabb, grid, i), detail::min_row<T>(aabb, grid, i));
                cellidx[cursor[cell]++] = (I)i;
            }
        });
    }

    /**
     * \brief The bounds that let the query skip a cell without reading it.
     *
     * \p row_prefix comes out as the running maximum, left to right along each
     * row, of the largest upper bound on the first axis of the boxes binned in
     * each cell. A binary search on it gives the first column of a row that holds
     * anything reaching back to a given coordinate, so every column before it is
     * skipped at once rather than one at a time. \p cell_hi1 is the same largest
     * upper bound on the second axis, per cell, and is tested directly.
     *
     * A cell with no boxes carries the lowest representable value, which fails
     * both tests, so an empty grid costs the query nothing.
     *
     * The maximum is written before it is scanned, and `sccd::cummax` reads
     * `in[i]` before writing `out[i]`, so the scan runs in place over
     * \p row_prefix and no third array is needed. Rows are independent.
     */
    template <typename T>
    static void cell2dmin_bounds(const ptrdiff_t n,
                                 T** const SCCD_RESTRICT aabb,
                                 const Cell2DGrid<T>& grid,
                                 const Cell2DPartition& part,
                                 T* const SCCD_RESTRICT row_prefix,
                                 T* const SCCD_RESTRICT cell_hi1) {
        const ptrdiff_t ncells = grid.ncells();
        const T floor = std::numeric_limits<T>::lowest();
        const int d0 = grid.axis0;
        const int d1 = grid.axis1;

        sccd::parallel_for_br(0, ncells, [&](const ptrdiff_t begin, const ptrdiff_t end) {
            for (ptrdiff_t c = begin; c < end; ++c) {
                row_prefix[c] = floor;
                cell_hi1[c] = floor;
            }
        });

        // Each box writes only the cell its minimum corner is in, and the blocks
        // own disjoint rows, so these are ordinary writes and not atomics.
        const auto accumulate = [&](const ptrdiff_t i) {
            const ptrdiff_t cell =
                grid.cell_of(detail::min_col<T>(aabb, grid, i), detail::min_row<T>(aabb, grid, i));
            row_prefix[cell] = sccd::max<T>(row_prefix[cell], aabb[3 + d0][i]);
            cell_hi1[cell] = sccd::max<T>(cell_hi1[cell], aabb[3 + d1][i]);
        };

        if (part.serial()) {
            for (ptrdiff_t i = 0; i < n; ++i) accumulate(i);
        } else {
            sccd::parallel_for_chunks(0, part.nblocks, [&](const ptrdiff_t b) {
                const ptrdiff_t from = part.blockptr[(size_t)b];
                const ptrdiff_t to = part.blockptr[(size_t)b + 1];
                for (ptrdiff_t e = from; e < to; ++e) accumulate((ptrdiff_t)part.blockbox[(size_t)e]);
            });
        }

        sccd::parallel_for_chunks(0, grid.n1, [&](const ptrdiff_t r) {
            T* const row = row_prefix + r * (ptrdiff_t)grid.n0;
            sccd::cummax<T>(grid.n0, row, row);
        });
    }

    /**
     * \brief Order each cell's entries on the axis the grid does not use.
     *
     * The grid culls on two axes and the cell test on the third, so inside a cell
     * the entries are in no useful order and every one of them is tested. Sorting
     * them by their minimum on that third axis makes the scan monotone: once an
     * entry begins past the querying box's maximum, so does every entry after it,
     * and the scan stops. \p cell_key holds those minima in the sorted order, one
     * per entry, so the test reads a contiguous array rather than chasing the
     * index into the box arrays.
     *
     * \p cell_hi2 is the largest maximum on the same axis, per cell, which rules
     * out a whole cell before its entries are touched at all.
     *
     * Cells are independent, so this is one parallel pass over them. The
     * comparison breaks ties on the box index, which makes the order total and
     * the result the same on every run whatever the sort does with equal keys.
     */
    template <typename T, typename I>
    static void cell2dmin_sort_cells(T** const SCCD_RESTRICT aabb,
                                     const Cell2DGrid<T>& grid,
                                     const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                     I* const SCCD_RESTRICT cellidx,
                                     T* const SCCD_RESTRICT cell_key,
                                     T* const SCCD_RESTRICT cell_hi2) {
        const int d2 = grid.axis2();
        const T* const SCCD_RESTRICT lo2 = aabb[d2];
        const T* const SCCD_RESTRICT hi2 = aabb[3 + d2];
        const T floor = std::numeric_limits<T>::lowest();

        sccd::parallel_for_br(0, grid.ncells(), [&](const ptrdiff_t begin, const ptrdiff_t end) {
            for (ptrdiff_t c = begin; c < end; ++c) {
                const ptrdiff_t from = cellptr[c];
                const ptrdiff_t to = cellptr[c + 1];

                std::sort(cellidx + from, cellidx + to, [&](const I a, const I b) {
                    return lo2[a] != lo2[b] ? lo2[a] < lo2[b] : a < b;
                });

                T top = floor;
                for (ptrdiff_t k = from; k < to; ++k) {
                    const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                    cell_key[k] = lo2[j];
                    top = sccd::max<T>(top, hi2[j]);
                }
                cell_hi2[c] = top;
            }
        });
    }

    namespace detail {

        /**
         * \brief Walk forward in linear cell order, reporting each partner once.
         *
         * A box sits in the one cell holding its minimum corner, so a partner is
         * met at most once and nothing has to be deduplicated across cells. The
         * pair is emitted by whichever of the two boxes comes first in row-major
         * order, with the index deciding inside one cell -- and that is the only
         * index test here.
         *
         * The walk spans the rows the box covers, and within a row the columns
         * from the first one holding anything that reaches back to the box, found
         * by binary search on the row's prefix maximum, up to the column of the
         * box's own maximum corner. A partner cannot begin past that column
         * without beginning past the box itself.
         *
         * \tparam SORTED The cells are ordered on the axis the grid does not use,
         * by cell2dmin_sort_cells. Then a whole cell can be ruled out by its
         * largest maximum on that axis, and the scan inside one stops at the first
         * entry that begins past this box -- every entry after it begins later
         * still. Without it the cell is scanned whole and the ordinary box test
         * does the culling. Both are compiled, so neither carries the other's
         * branches.
         */
        template <int NXE, bool SORTED, typename T, typename I, typename Visit>
        static inline void for_each_forward_self_partner(T** const SCCD_RESTRICT aabbs,
                                                         const ptrdiff_t fi,
                                                         I* const SCCD_RESTRICT idx,
                                                         I** const SCCD_RESTRICT elements,
                                                         const ptrdiff_t element_stride,
                                                         const I (&ev)[NXE],
                                                         const Cell2DGrid<T>& grid,
                                                         const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                                         const I* const SCCD_RESTRICT cellidx,
                                                         const T* const SCCD_RESTRICT row_prefix,
                                                         const T* const SCCD_RESTRICT cell_hi1,
                                                         const T* const SCCD_RESTRICT cell_key,
                                                         const T* const SCCD_RESTRICT cell_hi2,
                                                         Visit&& visit) {
            const T aminx = aabbs[0][fi], aminy = aabbs[1][fi], aminz = aabbs[2][fi];
            const T amaxx = aabbs[3][fi], amaxy = aabbs[4][fi], amaxz = aabbs[5][fi];

            const T amin0 = aabbs[grid.axis0][fi];
            const T amin1 = aabbs[grid.axis1][fi];
            T amin2 = 0, amax2 = 0;
            if constexpr (SORTED) {
                amin2 = aabbs[grid.axis2()][fi];
                amax2 = aabbs[3 + grid.axis2()][fi];
            }

            const int c0b = grid.clamp0(amin0), c0e = grid.clamp0(aabbs[3 + grid.axis0][fi]);
            const int c1b = grid.clamp1(amin1), c1e = grid.clamp1(aabbs[3 + grid.axis1][fi]);
            const ptrdiff_t own = grid.cell_of(c0b, c1b);

            for (int c1 = c1b; c1 <= c1e; ++c1) {
                const ptrdiff_t row = (ptrdiff_t)c1 * grid.n0;
                const T* const pre = row_prefix + row;
                int c0 = (int)(std::lower_bound(pre, pre + c0e + 1, amin0) - pre);
                if (c1 == c1b) {
                    c0 = sccd::max<int>(c0, c0b);
                }

                for (; c0 <= c0e; ++c0) {
                    const ptrdiff_t cell = row + c0;
                    if (cell_hi1[cell] < amin1) {
                        continue;
                    }
                    if constexpr (SORTED) {
                        if (cell_hi2[cell] < amin2) {
                            continue;
                        }
                    }

                    const ptrdiff_t begin = cellptr[cell];
                    const ptrdiff_t end = cellptr[cell + 1];

                    for (ptrdiff_t k = begin; k < end; ++k) {
                        if constexpr (SORTED) {
                            // Sorted ascending, so nothing after this one begins
                            // any earlier either.
                            if (cell_key[k] > amax2) {
                                break;
                            }
                        }

                        const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                        if (cell == own && j <= fi) {
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

    /**
     * \brief Self-overlap count over a minimum-corner binning, CRS offsets out.
     *
     * \tparam sorted The cells were ordered by cell2dmin_sort_cells, and
     * \p cell_key and \p cell_hi2 hold what it produced. False leaves them
     * unread, and null is the right thing to pass. It is a template parameter and
     * not an argument so that neither variant carries a branch the other needs.
     */
    template <int nxe, bool sorted, typename T, typename I>
    bool cell2dmin_count_self_overlaps(const ptrdiff_t element_count,
                                       T** const SCCD_RESTRICT aabbs,
                                       I* const SCCD_RESTRICT idx,
                                       const ptrdiff_t element_stride,
                                       I** const SCCD_RESTRICT elements,
                                       const Cell2DGrid<T>& grid,
                                       const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                       const I* const SCCD_RESTRICT cellidx,
                                       const T* const SCCD_RESTRICT row_prefix,
                                       const T* const SCCD_RESTRICT cell_hi1,
                                       const T* const SCCD_RESTRICT cell_key,
                                       const T* const SCCD_RESTRICT cell_hi2,
                                       ptrdiff_t* const SCCD_RESTRICT ccdptr) {
        ccdptr[0] = 0;

        // Dynamic, not the usual block schedule: a box spanning many rows walks
        // far more cells than one inside a single row, and a few of those in a
        // block would otherwise hold the whole block up.
        sccd::parallel_for_br_dynamic(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I idxi = idx[fi];
                I ev[nxe];
                for (int v = 0; v < nxe; ++v) {
                    ev[v] = elements[v][idxi * element_stride];
                }

                ptrdiff_t count = 0;
                detail::for_each_forward_self_partner<nxe, sorted, T, I>(
                    aabbs, fi, idx, elements, element_stride, ev, grid, cellptr, cellidx, row_prefix,
                    cell_hi1, cell_key, cell_hi2, [&](const ptrdiff_t, const I) { ++count; });
                ccdptr[fi + 1] = count;
            }
        });

        sccd::parallel_cum_sum_br(ccdptr, ccdptr + element_count + 1);
        return ccdptr[element_count] > 0;
    }

    /** \brief Write those pairs, as (min, max) to match the sweep. */
    template <int nxe, bool sorted, typename T, typename I>
    void cell2dmin_fill_self_overlaps(const ptrdiff_t element_count,
                                      T** const SCCD_RESTRICT aabbs,
                                      I* const SCCD_RESTRICT idx,
                                      const ptrdiff_t element_stride,
                                      I** const SCCD_RESTRICT elements,
                                      const Cell2DGrid<T>& grid,
                                      const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                      const I* const SCCD_RESTRICT cellidx,
                                      const T* const SCCD_RESTRICT row_prefix,
                                      const T* const SCCD_RESTRICT cell_hi1,
                                      const T* const SCCD_RESTRICT cell_key,
                                      const T* const SCCD_RESTRICT cell_hi2,
                                      const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                      I* const SCCD_RESTRICT first_out,
                                      I* const SCCD_RESTRICT second_out) {
        sccd::parallel_for_br_dynamic(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I idxi = idx[fi];
                I ev[nxe];
                for (int v = 0; v < nxe; ++v) {
                    ev[v] = elements[v][idxi * element_stride];
                }

                ptrdiff_t at = ccdptr[fi];
                detail::for_each_forward_self_partner<nxe, sorted, T, I>(
                    aabbs, fi, idx, elements, element_stride, ev, grid, cellptr, cellidx, row_prefix,
                    cell_hi1, cell_key, cell_hi2, [&](const ptrdiff_t, const I jidx) {
                        first_out[at] = sccd::min<I>(idxi, jidx);
                        second_out[at] = sccd::max<I>(idxi, jidx);
                        ++at;
                    });
            }
        });
    }

}  // namespace sccd

#endif  // SCCD_BROADPHASE_CELL2D_HPP
