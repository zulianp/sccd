#ifndef SCCD_BROADPHASE_HGRID2_HPP
#define SCCD_BROADPHASE_HGRID2_HPP

#include "sccd_broadphase_cell2d.hpp"

#include <cstdio>
#include <cstdlib>

/**
 * \file
 * \brief Broad phase over two nested cell lists, binned by centroid.
 *
 * A box goes into the one cell holding its centroid, and a cell is at least as
 * wide as any box binned into it. Two overlapping boxes then have centroids at
 * most one cell apart on each axis, so a fixed `3x3` stencil around a box's own
 * cell reaches every partner: the search radius is a property of the grid rather
 * than of the data, and the pair collection needs no duplicate test.
 *
 * One grid cannot have that and stay fine. Sizing a single cell to the widest
 * box hands every box the worst case: on armadillo-rollers one swept edge spans
 * 11% of the scene, which caps the grid at `8x9` and puts 36,363 edges into 72
 * cells. So the boxes are split in two by size.
 *
 * **The fine level** takes every box no wider than its cell, and that width is
 * chosen by costing every candidate against a histogram of box sizes, so the
 * level holds whatever share of the boxes pays for itself.
 *
 * **The coarse level** takes the rest, with a cell as wide as the widest box in
 * the scene. It is the grid a single level would have been forced to use, but
 * only the tail is in it, so its occupancy is the tail count and not \p n.
 *
 * Three pair classes follow, and each is emitted exactly once:
 *
 * - fine-fine and coarse-coarse, by a box walking the five cells of its own
 *   level whose linear index is at least its own, taking `j > i` inside its own
 *   cell;
 * - fine-coarse, by the **fine** box walking the coarse level's full `3x3`. It
 *   is complete there because a fine box is no wider than a coarse cell, so the
 *   two centroids are at most one coarse cell apart, and it is unique because a
 *   coarse box never looks at the fine level.
 *
 * A query from a separate list -- a face against the vertex grid -- is not
 * binned anywhere, so it walks both levels, with a radius per level taken from
 * its own extent. A box narrower than the level's cell gets the same `3x3`
 * everything else does, and a wider one pays in proportion to its own size
 * rather than making every other box pay.
 *
 * Both levels are binned and collected count-then-fill, like the rest of the
 * broad phase.
 */

namespace sccd {

    /**
     * \brief Two uniform grids over the same axes and origin, coarse over fine.
     *
     * \p w0 and \p w1 are the fine level's realised cell widths, which is what
     * the split is taken at: a box is fine when it fits in a fine cell, so the
     * test uses the width the grid ended up with and not the width it asked for.
     */
    template <typename T>
    struct HGrid2 {
        Cell2DGrid<T> fine;
        Cell2DGrid<T> coarse;
        T w0 = 0, w1 = 0;

        /** \brief True when box \p i fits in a fine cell on both axes. */
        bool is_fine(T** const SCCD_RESTRICT aabb, const ptrdiff_t i) const {
            return (aabb[3 + fine.axis0][i] - aabb[fine.axis0][i]) <= w0 &&
                   (aabb[3 + fine.axis1][i] - aabb[fine.axis1][i]) <= w1;
        }
    };

    namespace detail {

        /**
         * \brief The grid of cells at least \p w0 by \p w1, inside the cell cap.
         *
         * The cap has to be met on the integers the caller sizes its cell array
         * from, and meeting it by clamping each side in turn would collapse one
         * of them: a grid asked for at ten thousand columns comes back one column
         * wide, which prunes on one axis only and reads as cheap to anything
         * costing it by cell count. So both sides are scaled by the same factor
         * first, which keeps the cell's shape, and the clamp is left as the last
         * resort it was written to be.
         */
        static inline void fit_cells(const double span0,
                                     const double span1,
                                     const double w0,
                                     const double w1,
                                     const ptrdiff_t n,
                                     int& n0,
                                     int& n1) {
            const double cap = 4.0 * (double)sccd::max<ptrdiff_t>(n, 1);
            double want0 = span0 / w0;
            double want1 = span1 / w1;
            if (!(want0 >= 1)) want0 = 1;  // the negated form also catches NaN
            if (!(want1 >= 1)) want1 = 1;

            const double total = want0 * want1;
            if (total > cap) {
                const double s = std::sqrt(cap / total);
                want0 *= s;
                want1 *= s;
                if (want0 < 1) want0 = 1;
                if (want1 < 1) want1 = 1;
            }
            if (want0 > cap) want0 = cap;
            if (want1 > cap) want1 = cap;

            n0 = (int)sccd::max<double>(1.0, std::floor(want0));
            n1 = (int)sccd::max<double>(1.0, std::floor(want1));
            cap_cells(sccd::max<ptrdiff_t>(n, 1), n0, n1);
        }

        /** \brief Octaves of the size histogram, and which one holds a box of mean size. */
        static const int SCCD_HGRID2_BUCKETS = 32;
        static const int SCCD_HGRID2_UNIT = 16;

        /**
         * \brief The octave of a box's size relative to the mean box.
         *
         * One number for the two axes, taken as the larger of the two ratios, so
         * that a single threshold on it is exactly the test "fits in a cell whose
         * sides are that multiple of the mean extents". Keeping the cell's aspect
         * at the mean box's aspect means a cell is square in the units the boxes
         * are square in, which the domain's own aspect would not give.
         */
        template <typename T>
        static inline int size_bucket(const T e0, const T e1, const T mean0, const T mean1) {
            const double s = sccd::max<double>((double)(e0 / mean0), (double)(e1 / mean1));
            if (!(s > 0)) return 0;  // also catches NaN
            int exponent = 0;
            std::frexp(s, &exponent);
            // frexp gives s = m * 2^exponent with m in [0.5, 1), so the floor of
            // the base-two logarithm is exponent - 1.
            const int b = SCCD_HGRID2_UNIT + exponent - 1;
            if (b < 0) return 0;
            return b >= SCCD_HGRID2_BUCKETS ? SCCD_HGRID2_BUCKETS - 1 : b;
        }

        /** \brief How many boxes fall in each octave of size. */
        struct SizeHistogram {
            ptrdiff_t count[SCCD_HGRID2_BUCKETS];
        };

    }  // namespace detail

    /**
     * \brief Size both grids over \p n boxes and pick the split between them.
     *
     * Where to cut is the only free parameter of the structure, and a fixed
     * fraction of boxes is the wrong way to fix it: how much mass the tail
     * carries varies by scene and by step, and a cut above the tail puts every
     * box back on one coarse grid. So the cut is chosen by cost.
     *
     * Each candidate cut is one octave of box size. A histogram of sizes gives
     * exactly how many boxes fall either side of it, the cell counts follow from
     * the widths, and the estimate is the work the three pair classes would do:
     * the fine level's own query, the coarse level's own query, and the fine
     * boxes' walk over the coarse level, plus the cost of touching the two cell
     * arrays at all, which is what stops the search from asking for a grid finer
     * than it can pay for. The cut with the smallest estimate wins.
     *
     * The estimate uses mean occupancy, so it understates a clustered level. It
     * is a means of ordering the candidates rather than a prediction.
     */
    template <typename T>
    static void hgrid2_setup(const ptrdiff_t n, T** const SCCD_RESTRICT aabb, HGrid2<T>& grid) {
        detail::choose_two_axes<T>(n, aabb, grid.fine.axis0, grid.fine.axis1);
        grid.coarse.axis0 = grid.fine.axis0;
        grid.coarse.axis1 = grid.fine.axis1;

        const int d0 = grid.fine.axis0;
        const int d1 = grid.fine.axis1;

        struct Bounds {
            T lo0, hi0, lo1, hi1, ext0, ext1, sum0, sum1;
        };

        const Bounds b = sccd::parallel_tiled_reduce<Bounds>(
            0,
            n,
            [&](const ptrdiff_t begin, const ptrdiff_t end) {
                Bounds e{aabb[d0][begin], aabb[3 + d0][begin], aabb[d1][begin], aabb[3 + d1][begin], 0, 0, 0, 0};
                for (ptrdiff_t i = begin; i < end; ++i) {
                    e.lo0 = sccd::min<T>(e.lo0, aabb[d0][i]);
                    e.hi0 = sccd::max<T>(e.hi0, aabb[3 + d0][i]);
                    e.lo1 = sccd::min<T>(e.lo1, aabb[d1][i]);
                    e.hi1 = sccd::max<T>(e.hi1, aabb[3 + d1][i]);
                    e.ext0 = sccd::max<T>(e.ext0, aabb[3 + d0][i] - aabb[d0][i]);
                    e.ext1 = sccd::max<T>(e.ext1, aabb[3 + d1][i] - aabb[d1][i]);
                    e.sum0 += aabb[3 + d0][i] - aabb[d0][i];
                    e.sum1 += aabb[3 + d1][i] - aabb[d1][i];
                }
                return e;
            },
            [](const Bounds x, const Bounds y) {
                return Bounds{sccd::min<T>(x.lo0, y.lo0),
                              sccd::max<T>(x.hi0, y.hi0),
                              sccd::min<T>(x.lo1, y.lo1),
                              sccd::max<T>(x.hi1, y.hi1),
                              sccd::max<T>(x.ext0, y.ext0),
                              sccd::max<T>(x.ext1, y.ext1),
                              x.sum0 + y.sum0,
                              x.sum1 + y.sum1};
            });

        const ptrdiff_t safe_n = sccd::max<ptrdiff_t>(n, 1);
        const T span0 = sccd::max<T>(b.hi0 - b.lo0, std::numeric_limits<T>::min());
        const T span1 = sccd::max<T>(b.hi1 - b.lo1, std::numeric_limits<T>::min());
        const T mean0 = sccd::max<T>(b.sum0 / (T)safe_n, span0 / (T)(1 << 20));
        const T mean1 = sccd::max<T>(b.sum1 / (T)safe_n, span1 / (T)(1 << 20));

        const detail::SizeHistogram hist = sccd::parallel_tiled_reduce<detail::SizeHistogram>(
            0,
            n,
            [&](const ptrdiff_t begin, const ptrdiff_t end) {
                detail::SizeHistogram h{};
                for (ptrdiff_t i = begin; i < end; ++i) {
                    h.count[detail::size_bucket<T>(aabb[3 + d0][i] - aabb[d0][i],
                                                   aabb[3 + d1][i] - aabb[d1][i],
                                                   mean0,
                                                   mean1)] += 1;
                }
                return h;
            },
            [](detail::SizeHistogram x, const detail::SizeHistogram y) {
                for (int k = 0; k < detail::SCCD_HGRID2_BUCKETS; ++k) x.count[k] += y.count[k];
                return x;
            });

        // The coarse level's shape does not depend on the cut: its cell holds the
        // widest box there is either way, so it is the same for every candidate
        // and can be costed once.
        int cn0 = 1, cn1 = 1;
        detail::fit_cells((double)span0,
                          (double)span1,
                          (double)sccd::max<T>(b.ext0, mean0),
                          (double)sccd::max<T>(b.ext1, mean1),
                          safe_n,
                          cn0,
                          cn1);
        const double coarse_cells = (double)cn0 * (double)cn1;

        int best_bucket = detail::SCCD_HGRID2_BUCKETS - 1;
        double best_cost = -1;
        ptrdiff_t below = 0;
        for (int k = 0; k < detail::SCCD_HGRID2_BUCKETS; ++k) {
            below += hist.count[k];
            if (below == 0) continue;  // no box this small; the cut would be empty

            const double m = std::ldexp(1.0, k + 1 - detail::SCCD_HGRID2_UNIT);
            int fn0 = 1, fn1 = 1;
            detail::fit_cells(
                (double)span0, (double)span1, m * (double)mean0, m * (double)mean1, safe_n, fn0, fn1);

            const double fine_cells = (double)fn0 * (double)fn1;
            const double nf = (double)below;
            const double nc = (double)(n - below);
            const double of = nf / fine_cells;
            const double oc = nc / coarse_cells;

            // Five cells for a level's own forward-half walk, nine for a fine
            // box's walk over the coarse level, and one touch per cell of each
            // array to bin at all.
            const double cost = 5.0 * nf * of + 5.0 * nc * oc + 9.0 * nf * oc + fine_cells + coarse_cells;
            if (best_cost < 0 || cost < best_cost) {
                best_cost = cost;
                best_bucket = k;
            }
            if (below == n) break;  // every box is fine from here on
        }

        const double cut = std::ldexp(1.0, best_bucket + 1 - detail::SCCD_HGRID2_UNIT);
        const T eps = (T)1 + std::numeric_limits<T>::epsilon();

        detail::fit_cells((double)span0,
                          (double)span1,
                          cut * (double)mean0,
                          cut * (double)mean1,
                          safe_n,
                          grid.fine.n0,
                          grid.fine.n1);
        grid.fine.min0 = b.lo0;
        grid.fine.min1 = b.lo1;
        grid.fine.inv0 = (T)grid.fine.n0 / (span0 * eps);
        grid.fine.inv1 = (T)grid.fine.n1 / (span1 * eps);

        // The realised widths, which is what membership is tested against: the
        // floor and the cell cap both make a cell wider than it asked to be, and
        // a box that fits the wider cell belongs on the fine level.
        grid.w0 = span0 * eps / (T)grid.fine.n0;
        grid.w1 = span1 * eps / (T)grid.fine.n1;

        // The coarse cell holds the widest box there is, and is never finer than
        // the fine cell, so a box rejected by the fine level always fits here.
        detail::fit_cells((double)span0,
                          (double)span1,
                          (double)sccd::max<T>(b.ext0, grid.w0),
                          (double)sccd::max<T>(b.ext1, grid.w1),
                          safe_n,
                          grid.coarse.n0,
                          grid.coarse.n1);
        grid.coarse.min0 = b.lo0;
        grid.coarse.min1 = b.lo1;
        grid.coarse.inv0 = (T)grid.coarse.n0 / (span0 * eps);
        grid.coarse.inv1 = (T)grid.coarse.n1 / (span1 * eps);

        // `SCCD_HGRID2_VERBOSE` reports where the cut fell and what each level
        // came out as. How much of the scene is on the coarse level is the
        // number that decides whether two levels were worth having, and it is a
        // property of the step's geometry rather than of the code.
        if (getenv("SCCD_HGRID2_VERBOSE")) {
            ptrdiff_t coarse_n = 0;
            for (int k = best_bucket + 1; k < detail::SCCD_HGRID2_BUCKETS; ++k) coarse_n += hist.count[k];
            fprintf(stderr,
                    "sccd hgrid2: n %ld  axes %d,%d  mean box %.4g/%.4g  widest %.4g/%.4g  "
                    "cut %.4g means  fine %dx%d  coarse %dx%d  on coarse %ld (%.2f%%)\n",
                    (long)n,
                    d0,
                    d1,
                    (double)mean0,
                    (double)mean1,
                    (double)b.ext0,
                    (double)b.ext1,
                    cut,
                    grid.fine.n0,
                    grid.fine.n1,
                    grid.coarse.n0,
                    grid.coarse.n1,
                    (long)coarse_n,
                    100.0 * (double)coarse_n / (double)safe_n);
        }
    }

    namespace detail {

        /** \brief The midpoint of a box on one axis, written so it cannot overflow. */
        template <typename T>
        static inline T midpoint(const T lo, const T hi) {
            return lo + (hi - lo) * T(0.5);
        }

        /** \brief The column of the cell holding box \p i's centroid. */
        template <typename T>
        static inline int centroid_col(T** const SCCD_RESTRICT aabb, const Cell2DGrid<T>& g, const ptrdiff_t i) {
            return g.clamp0(midpoint<T>(aabb[g.axis0][i], aabb[3 + g.axis0][i]));
        }

        /** \brief The row of the cell holding box \p i's centroid. */
        template <typename T>
        static inline int centroid_row(T** const SCCD_RESTRICT aabb, const Cell2DGrid<T>& g, const ptrdiff_t i) {
            return g.clamp1(midpoint<T>(aabb[g.axis1][i], aabb[3 + g.axis1][i]));
        }

    }  // namespace detail

    /**
     * \brief Group one level's boxes by the block of cell rows they land in.
     *
     * The same device as the cell list's partition and for the same reason: the
     * histogram and the scatter write per-cell counters, so they are split by
     * cell row and no two blocks touch a cell. A box is binned by its centroid,
     * so it lands in exactly one block's list; boxes of the other level are left
     * out of the lists entirely.
     */
    template <typename T>
    static void hgrid2_partition(const ptrdiff_t n,
                                 T** const SCCD_RESTRICT aabb,
                                 const HGrid2<T>& grid,
                                 const bool fine_level,
                                 Cell2DPartition& part) {
        const Cell2DGrid<T>& g = fine_level ? grid.fine : grid.coarse;

        part.blockptr.clear();
        part.blockbox.clear();

        const int max_workers = sccd::max_concurrency();
        if (n < SCCD_CELL2D_MIN_PARALLEL || max_workers <= 1 || g.n1 <= 1) {
            part.nblocks = 1;
            part.rows_per_block = g.n1;
            return;
        }

        const int want = sccd::min<int>(g.n1, max_workers * 4);
        part.nblocks = sccd::max<int>(1, want);
        part.rows_per_block = (g.n1 + part.nblocks - 1) / part.nblocks;
        part.nblocks = (g.n1 + part.rows_per_block - 1) / part.rows_per_block;

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
                if (grid.is_fine(aabb, i) != fine_level) continue;
                row[detail::centroid_row<T>(aabb, g, i) / rpb] += 1;
            }
        });

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
                if (grid.is_fine(aabb, i) != fine_level) continue;
                part.blockbox[(size_t)(row[detail::centroid_row<T>(aabb, g, i) / rpb]++)] = (int)i;
            }
        });
    }

    /** \brief Count one level's boxes per cell, then prefix-sum into CRS offsets. */
    template <typename T>
    static void hgrid2_count(const ptrdiff_t n,
                             T** const SCCD_RESTRICT aabb,
                             const HGrid2<T>& grid,
                             const bool fine_level,
                             const Cell2DPartition& part,
                             ptrdiff_t* const SCCD_RESTRICT cellptr) {
        const Cell2DGrid<T>& g = fine_level ? grid.fine : grid.coarse;
        const ptrdiff_t ncells = g.ncells();

        if (part.serial()) {
            std::memset(cellptr, 0, sizeof(ptrdiff_t) * (size_t)(ncells + 1));
            for (ptrdiff_t i = 0; i < n; ++i) {
                if (grid.is_fine(aabb, i) != fine_level) continue;
                cellptr[g.cell_of(detail::centroid_col<T>(aabb, g, i), detail::centroid_row<T>(aabb, g, i)) + 1] += 1;
            }
            for (ptrdiff_t i = 0; i < ncells; ++i) cellptr[i + 1] += cellptr[i];
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
                cellptr[g.cell_of(detail::centroid_col<T>(aabb, g, i), detail::centroid_row<T>(aabb, g, i)) + 1] += 1;
            }
        });

        sccd::parallel_cum_sum_br(cellptr, cellptr + ncells + 1);
    }

    /** \brief Scatter one level's box indices into the cells counted above. */
    template <typename T, typename I>
    static void hgrid2_fill(const ptrdiff_t n,
                            T** const SCCD_RESTRICT aabb,
                            const HGrid2<T>& grid,
                            const bool fine_level,
                            const Cell2DPartition& part,
                            const ptrdiff_t* const SCCD_RESTRICT cellptr,
                            I* const SCCD_RESTRICT cellidx,
                            ptrdiff_t* const SCCD_RESTRICT cursor) {
        const Cell2DGrid<T>& g = fine_level ? grid.fine : grid.coarse;
        const ptrdiff_t ncells = g.ncells();

        if (part.serial()) {
            std::memcpy(cursor, cellptr, sizeof(ptrdiff_t) * (size_t)ncells);
            for (ptrdiff_t i = 0; i < n; ++i) {
                if (grid.is_fine(aabb, i) != fine_level) continue;
                const ptrdiff_t cell =
                    g.cell_of(detail::centroid_col<T>(aabb, g, i), detail::centroid_row<T>(aabb, g, i));
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
                    g.cell_of(detail::centroid_col<T>(aabb, g, i), detail::centroid_row<T>(aabb, g, i));
                cellidx[cursor[cell]++] = (I)i;
            }
        });
    }

    namespace detail {

        /**
         * \brief Visit a rectangle of cells around (\p col, \p row), clipped to the grid.
         *
         * \p half restricts the walk to the cells whose linear index is at least
         * the middle one's: the row above from the column before, and its own row
         * from its own column. That is the half of the stencil that makes a self
         * query report each unordered pair once.
         */
        template <typename T, typename Visit>
        static inline void for_each_stencil_cell(const Cell2DGrid<T>& g,
                                                 const int col,
                                                 const int row,
                                                 const int r0,
                                                 const int r1,
                                                 const bool half,
                                                 Visit&& visit) {
            const int dr_begin = half ? 0 : -r1;
            for (int dr = dr_begin; dr <= r1; ++dr) {
                const int c1 = row + dr;
                if (c1 < 0 || c1 >= g.n1) continue;

                const int dc_begin = (half && dr == 0) ? 0 : -r0;
                for (int dc = dc_begin; dc <= r0; ++dc) {
                    const int c0 = col + dc;
                    if (c0 < 0 || c0 >= g.n0) continue;
                    visit(g.cell_of(c0, c1), dr == 0 && dc == 0);
                }
            }
        }

        /**
         * \brief How many cells of \p g a box of extent \p e has to reach.
         *
         * A partner in \p g is no wider than one of its cells, so two overlapping
         * centroids are at most `(e + w) / 2` apart and the walk spans that many
         * cells either side. A box narrower than the cell gets the ordinary
         * radius of one.
         */
        template <typename T>
        static inline int stencil_radius(const T e, const T w) {
            if (!(e > w)) return 1;
            const double r = std::ceil((double)((e + w) / (w + w)));
            return r > 1e6 ? 1000000 : (int)r;
        }

        /**
         * \brief Report the partners of box \p fi found on one level.
         *
         * \p own_level says the query box is itself binned in this grid, which
         * turns the walk into its forward half and adds the `j > fi` test in its
         * own cell. Those two together are what make a self query report each
         * unordered pair once.
         */
        template <int F, int S, typename T, typename I, typename Visit>
        static inline void hgrid2_walk_level(T** const SCCD_RESTRICT first_aabbs,
                                             const ptrdiff_t fi,
                                             T** const SCCD_RESTRICT second_aabbs,
                                             const I* const SCCD_RESTRICT second_idx,
                                             I** const SCCD_RESTRICT second_elements,
                                             const ptrdiff_t second_element_stride,
                                             const I (&ev)[F],
                                             const Cell2DGrid<T>& g,
                                             const ptrdiff_t* const SCCD_RESTRICT cellptr,
                                             const I* const SCCD_RESTRICT cellidx,
                                             const int r0,
                                             const int r1,
                                             const bool own_level,
                                             Visit&& visit) {
            const T aminx = first_aabbs[0][fi], aminy = first_aabbs[1][fi], aminz = first_aabbs[2][fi];
            const T amaxx = first_aabbs[3][fi], amaxy = first_aabbs[4][fi], amaxz = first_aabbs[5][fi];

            const int col = centroid_col<T>(first_aabbs, g, fi);
            const int row = centroid_row<T>(first_aabbs, g, fi);

            for_each_stencil_cell<T>(
                g, col, row, r0, r1, own_level, [&](const ptrdiff_t cell, const bool centre) {
                    const ptrdiff_t begin = cellptr[cell];
                    const ptrdiff_t end = cellptr[cell + 1];

                    for (ptrdiff_t k = begin; k < end; ++k) {
                        const ptrdiff_t j = (ptrdiff_t)cellidx[k];
                        if (own_level && centre && j <= fi) {
                            continue;
                        }

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

                        const I jidx = second_idx[j];
                        bool share = false;
                        if constexpr (S > 1) {
                            I sev[S];
                            for (int v = 0; v < S; ++v) {
                                sev[v] = second_elements[v][jidx * second_element_stride];
                            }
                            share = sccd::detail::shares_vertex<F, S>(ev, sev);
                        } else {
                            for (int a = 0; a < F; ++a) {
                                if (ev[a] == jidx) {
                                    share = true;
                                    break;
                                }
                            }
                        }
                        if (share) {
                            continue;
                        }

                        visit(j, jidx);
                    }
                });
        }

        /**
         * \brief Every partner of box \p fi in a self query, across both levels.
         *
         * A fine box takes its own level's forward half and then the whole of the
         * coarse level, because a coarse box never looks down. A coarse box takes
         * only its own level's forward half.
         */
        template <int NXE, typename T, typename I, typename Visit>
        static inline void hgrid2_self_partners(T** const SCCD_RESTRICT aabbs,
                                                const ptrdiff_t fi,
                                                I* const SCCD_RESTRICT idx,
                                                I** const SCCD_RESTRICT elements,
                                                const ptrdiff_t element_stride,
                                                const I (&ev)[NXE],
                                                const HGrid2<T>& grid,
                                                const ptrdiff_t* const SCCD_RESTRICT fine_cellptr,
                                                const I* const SCCD_RESTRICT fine_cellidx,
                                                const ptrdiff_t* const SCCD_RESTRICT coarse_cellptr,
                                                const I* const SCCD_RESTRICT coarse_cellidx,
                                                Visit&& visit) {
            const bool fine = grid.is_fine(aabbs, fi);

            hgrid2_walk_level<NXE, NXE, T, I>(aabbs,
                                              fi,
                                              aabbs,
                                              idx,
                                              elements,
                                              element_stride,
                                              ev,
                                              fine ? grid.fine : grid.coarse,
                                              fine ? fine_cellptr : coarse_cellptr,
                                              fine ? fine_cellidx : coarse_cellidx,
                                              1,
                                              1,
                                              true,
                                              visit);

            // Nine cell reads per box is worth a branch when there is nothing
            // in them: a scene whose boxes are all of a size leaves the coarse
            // level empty, and every fine box would otherwise walk it.
            if (fine && coarse_cellptr[grid.coarse.ncells()] > 0) {
                // A fine box is no wider than a coarse cell, so the two centroids
                // are at most one coarse cell apart and the ordinary radius holds.
                hgrid2_walk_level<NXE, NXE, T, I>(aabbs,
                                                  fi,
                                                  aabbs,
                                                  idx,
                                                  elements,
                                                  element_stride,
                                                  ev,
                                                  grid.coarse,
                                                  coarse_cellptr,
                                                  coarse_cellidx,
                                                  1,
                                                  1,
                                                  false,
                                                  visit);
            }
        }

        /**
         * \brief Every partner of a box from another list, across both levels.
         *
         * The query box is binned nowhere, so it walks both levels whole, with a
         * radius per level taken from its own extent. Nothing has to be
         * deduplicated: a box of the second list sits in one cell of one level.
         */
        template <int F, int S, typename T, typename I, typename Visit>
        static inline void hgrid2_partners(T** const SCCD_RESTRICT first_aabbs,
                                           const ptrdiff_t fi,
                                           T** const SCCD_RESTRICT second_aabbs,
                                           const I* const SCCD_RESTRICT second_idx,
                                           I** const SCCD_RESTRICT second_elements,
                                           const ptrdiff_t second_element_stride,
                                           const I (&ev)[F],
                                           const HGrid2<T>& grid,
                                           const ptrdiff_t* const SCCD_RESTRICT fine_cellptr,
                                           const I* const SCCD_RESTRICT fine_cellidx,
                                           const ptrdiff_t* const SCCD_RESTRICT coarse_cellptr,
                                           const I* const SCCD_RESTRICT coarse_cellidx,
                                           Visit&& visit) {
            const T e0 = first_aabbs[3 + grid.fine.axis0][fi] - first_aabbs[grid.fine.axis0][fi];
            const T e1 = first_aabbs[3 + grid.fine.axis1][fi] - first_aabbs[grid.fine.axis1][fi];

            const T cw0 = (T)1 / grid.coarse.inv0;
            const T cw1 = (T)1 / grid.coarse.inv1;
            const bool coarse_occupied = coarse_cellptr[grid.coarse.ncells()] > 0;

            hgrid2_walk_level<F, S, T, I>(first_aabbs,
                                          fi,
                                          second_aabbs,
                                          second_idx,
                                          second_elements,
                                          second_element_stride,
                                          ev,
                                          grid.fine,
                                          fine_cellptr,
                                          fine_cellidx,
                                          stencil_radius<T>(e0, grid.w0),
                                          stencil_radius<T>(e1, grid.w1),
                                          false,
                                          visit);

            if (coarse_occupied) {
                hgrid2_walk_level<F, S, T, I>(first_aabbs,
                                              fi,
                                              second_aabbs,
                                              second_idx,
                                              second_elements,
                                              second_element_stride,
                                              ev,
                                              grid.coarse,
                                              coarse_cellptr,
                                              coarse_cellidx,
                                              stencil_radius<T>(e0, cw0),
                                              stencil_radius<T>(e1, cw1),
                                              false,
                                              visit);
            }
        }

    }  // namespace detail

    /** \brief Self-overlap count: one list against itself, CRS offsets out. */
    template <int nxe, typename T, typename I>
    bool hgrid2_count_self_overlaps(const ptrdiff_t element_count,
                                    T** const SCCD_RESTRICT aabbs,
                                    I* const SCCD_RESTRICT idx,
                                    const ptrdiff_t element_stride,
                                    I** const SCCD_RESTRICT elements,
                                    const HGrid2<T>& grid,
                                    const ptrdiff_t* const SCCD_RESTRICT fine_cellptr,
                                    const I* const SCCD_RESTRICT fine_cellidx,
                                    const ptrdiff_t* const SCCD_RESTRICT coarse_cellptr,
                                    const I* const SCCD_RESTRICT coarse_cellidx,
                                    ptrdiff_t* const SCCD_RESTRICT ccdptr) {
        ccdptr[0] = 0;

        sccd::parallel_for_br(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I idxi = idx[fi];
                I ev[nxe];
                for (int v = 0; v < nxe; ++v) ev[v] = elements[v][idxi * element_stride];

                ptrdiff_t count = 0;
                detail::hgrid2_self_partners<nxe, T, I>(aabbs,
                                                        fi,
                                                        idx,
                                                        elements,
                                                        element_stride,
                                                        ev,
                                                        grid,
                                                        fine_cellptr,
                                                        fine_cellidx,
                                                        coarse_cellptr,
                                                        coarse_cellidx,
                                                        [&](const ptrdiff_t, const I) { ++count; });
                ccdptr[fi + 1] = count;
            }
        });

        sccd::parallel_cum_sum_br(ccdptr, ccdptr + element_count + 1);
        return ccdptr[element_count] > 0;
    }

    /** \brief Write the self-overlap pairs, as (min, max) to match the sweep. */
    template <int nxe, typename T, typename I>
    void hgrid2_fill_self_overlaps(const ptrdiff_t element_count,
                                   T** const SCCD_RESTRICT aabbs,
                                   I* const SCCD_RESTRICT idx,
                                   const ptrdiff_t element_stride,
                                   I** const SCCD_RESTRICT elements,
                                   const HGrid2<T>& grid,
                                   const ptrdiff_t* const SCCD_RESTRICT fine_cellptr,
                                   const I* const SCCD_RESTRICT fine_cellidx,
                                   const ptrdiff_t* const SCCD_RESTRICT coarse_cellptr,
                                   const I* const SCCD_RESTRICT coarse_cellidx,
                                   const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                                   I* const SCCD_RESTRICT first_out,
                                   I* const SCCD_RESTRICT second_out) {
        sccd::parallel_for_br(0, element_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I idxi = idx[fi];
                I ev[nxe];
                for (int v = 0; v < nxe; ++v) ev[v] = elements[v][idxi * element_stride];

                ptrdiff_t at = ccdptr[fi];
                detail::hgrid2_self_partners<nxe, T, I>(aabbs,
                                                        fi,
                                                        idx,
                                                        elements,
                                                        element_stride,
                                                        ev,
                                                        grid,
                                                        fine_cellptr,
                                                        fine_cellidx,
                                                        coarse_cellptr,
                                                        coarse_cellidx,
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
    bool hgrid2_count_overlaps(const ptrdiff_t first_count,
                               T** const SCCD_RESTRICT first_aabbs,
                               I* const SCCD_RESTRICT first_idx,
                               const ptrdiff_t first_element_stride,
                               I** const SCCD_RESTRICT first_elements,
                               T** const SCCD_RESTRICT second_aabbs,
                               I* const SCCD_RESTRICT second_idx,
                               const ptrdiff_t second_element_stride,
                               I** const SCCD_RESTRICT second_elements,
                               const HGrid2<T>& grid,
                               const ptrdiff_t* const SCCD_RESTRICT fine_cellptr,
                               const I* const SCCD_RESTRICT fine_cellidx,
                               const ptrdiff_t* const SCCD_RESTRICT coarse_cellptr,
                               const I* const SCCD_RESTRICT coarse_cellidx,
                               ptrdiff_t* const SCCD_RESTRICT ccdptr) {
        ccdptr[0] = 0;

        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I first_idxi = first_idx[fi];
                I ev[first_nxe];
                for (int v = 0; v < first_nxe; ++v) ev[v] = first_elements[v][first_idxi * first_element_stride];

                ptrdiff_t count = 0;
                detail::hgrid2_partners<first_nxe, second_nxe, T, I>(first_aabbs,
                                                                     fi,
                                                                     second_aabbs,
                                                                     second_idx,
                                                                     second_elements,
                                                                     second_element_stride,
                                                                     ev,
                                                                     grid,
                                                                     fine_cellptr,
                                                                     fine_cellidx,
                                                                     coarse_cellptr,
                                                                     coarse_cellidx,
                                                                     [&](const ptrdiff_t, const I) { ++count; });
                ccdptr[fi + 1] = count;
            }
        });

        sccd::parallel_cum_sum_br(ccdptr, ccdptr + first_count + 1);
        return ccdptr[first_count] > 0;
    }

    /** \brief Write the pairs counted by hgrid2_count_overlaps into the CRS arrays. */
    template <int first_nxe, int second_nxe, typename T, typename I>
    void hgrid2_fill_overlaps(const ptrdiff_t first_count,
                              T** const SCCD_RESTRICT first_aabbs,
                              I* const SCCD_RESTRICT first_idx,
                              const ptrdiff_t first_element_stride,
                              I** const SCCD_RESTRICT first_elements,
                              T** const SCCD_RESTRICT second_aabbs,
                              I* const SCCD_RESTRICT second_idx,
                              const ptrdiff_t second_element_stride,
                              I** const SCCD_RESTRICT second_elements,
                              const HGrid2<T>& grid,
                              const ptrdiff_t* const SCCD_RESTRICT fine_cellptr,
                              const I* const SCCD_RESTRICT fine_cellidx,
                              const ptrdiff_t* const SCCD_RESTRICT coarse_cellptr,
                              const I* const SCCD_RESTRICT coarse_cellidx,
                              const ptrdiff_t* const SCCD_RESTRICT ccdptr,
                              I* const SCCD_RESTRICT first_out,
                              I* const SCCD_RESTRICT second_out) {
        sccd::parallel_for_br(0, first_count, [&](const ptrdiff_t rbegin, const ptrdiff_t rend) {
            for (ptrdiff_t fi = rbegin; fi < rend; ++fi) {
                const I first_idxi = first_idx[fi];
                I ev[first_nxe];
                for (int v = 0; v < first_nxe; ++v) ev[v] = first_elements[v][first_idxi * first_element_stride];

                ptrdiff_t at = ccdptr[fi];
                detail::hgrid2_partners<first_nxe, second_nxe, T, I>(first_aabbs,
                                                                     fi,
                                                                     second_aabbs,
                                                                     second_idx,
                                                                     second_elements,
                                                                     second_element_stride,
                                                                     ev,
                                                                     grid,
                                                                     fine_cellptr,
                                                                     fine_cellidx,
                                                                     coarse_cellptr,
                                                                     coarse_cellidx,
                                                                     [&](const ptrdiff_t, const I jidx) {
                                                                         first_out[at] = first_idxi;
                                                                         second_out[at] = jidx;
                                                                         ++at;
                                                                     });
            }
        });
    }

}  // namespace sccd

#endif  // SCCD_BROADPHASE_HGRID2_HPP
