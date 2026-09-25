// The 2D cell broad phase must find exactly the pair set the sweep finds.
//
// It is an optimisation, so the only thing that makes it acceptable is producing
// the same answer. Pairs are compared as sets, because the two produce them in
// different orders by construction: the sweep walks a sorted axis, the cell list
// walks cells.
//
// The cases are chosen to exercise what the two disagree about most easily --
// boxes spanning many cells (where the cell list must not report a pair once per
// shared cell), boxes outside the grid (where clamping must stay consistent
// between the binning and the duplicate check), and degenerate extents.

#include "sccd_broadphase_sweep.hpp"
#include "sccd_broadphase_cell2d.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <random>
#include <set>
#include <vector>

using scalar_t = double;
using idx_t = std::int32_t;

namespace {

    struct Boxes {
        std::vector<scalar_t> data[6];
        std::vector<idx_t> idx;
        std::vector<idx_t> elem[4];
        scalar_t* ptr[6];
        idx_t* elem_ptr[4];
        ptrdiff_t n = 0;

        void bind() {
            for (int d = 0; d < 6; ++d) ptr[d] = data[d].data();
            for (int d = 0; d < 4; ++d) elem_ptr[d] = elem[d].data();
        }
    };

    // nxe == 1 means the "element" is a single node (the vertex list); nxe == 3 a
    // triangle and nxe == 4 a quad. All three shapes appear in the real broad
    // phase, and the vertex count is what the pair collector uses to skip a
    // face's own vertices -- so it has to be exercised at each width.
    Boxes make_boxes(std::mt19937& rng, const ptrdiff_t n, const int nxe, const double spread, const double size) {
        std::uniform_real_distribution<double> pos(0.0, spread);
        std::uniform_real_distribution<double> ext(0.0, size);
        std::uniform_int_distribution<int> node(0, 999);

        Boxes b;
        b.n = n;
        for (int d = 0; d < 6; ++d) b.data[d].resize(n);
        b.idx.resize(n);
        for (int v = 0; v < nxe; ++v) b.elem[v].resize(n);

        for (ptrdiff_t i = 0; i < n; ++i) {
            for (int d = 0; d < 3; ++d) {
                const double lo = pos(rng);
                b.data[d][i] = lo;
                b.data[3 + d][i] = lo + ext(rng);
            }
            b.idx[i] = (idx_t)i;
            for (int v = 0; v < nxe; ++v) b.elem[v][i] = (idx_t)node(rng);
        }
        b.bind();
        return b;
    }

    using PairSet = std::set<std::pair<idx_t, idx_t>>;

    // An independent reference: every pair whose boxes overlap, minus the pairs
    // that share a node. Written out longhand on purpose -- the sweep and the
    // cell list are cross-checked against each other elsewhere in this file, and
    // two implementations agreeing says nothing if they share a bug in the
    // shared-node masking, which is exactly the code second_nxe > 1 turns on.
    template <int first_nxe, int second_nxe>
    PairSet brute_pairs(const Boxes& first, const Boxes& second, ptrdiff_t* masked = nullptr) {
        PairSet out;
        if (masked) *masked = 0;
        for (ptrdiff_t i = 0; i < first.n; ++i) {
            for (ptrdiff_t j = 0; j < second.n; ++j) {
                bool disjoint = false;
                for (int d = 0; d < 3; ++d) {
                    if (first.data[3 + d][i] < second.data[d][j] ||
                        second.data[3 + d][j] < first.data[d][i]) {
                        disjoint = true;
                        break;
                    }
                }
                if (disjoint) continue;

                bool shares = false;
                for (int a = 0; a < first_nxe && !shares; ++a) {
                    for (int b = 0; b < second_nxe; ++b) {
                        if (first.elem[a][i] == second.elem[b][j]) { shares = true; break; }
                    }
                }
                if (shares) {
                    if (masked) ++(*masked);
                } else {
                    out.insert({(idx_t)i, (idx_t)j});
                }
            }
        }
        return out;
    }

    template <int first_nxe, int second_nxe = 1>
    PairSet sweep_pairs(Boxes& first, Boxes& second) {
        // The sweep needs both lists sorted along a common axis; the cell list
        // does not, which is the point of it.
        std::vector<scalar_t> scratch(std::max(first.n, second.n) * 2);
        const int axis = sccd::choose_axis<scalar_t>(second.n, second.ptr);
        sccd::sort_along_axis(first.n, axis, first.ptr, first.idx.data(), scratch.data());
        sccd::sort_along_axis(second.n, axis, second.ptr, second.idx.data(), scratch.data());

        std::vector<scalar_t> cummax(second.n);
        sccd::cummax(second.n, second.ptr[3 + axis], cummax.data());

        std::vector<ptrdiff_t> ccdptr(first.n + 1, 0);
        const bool any = sccd::count_overlaps<first_nxe, second_nxe, scalar_t, idx_t>(axis,
                                                                             first.n,
                                                                             first.ptr,
                                                                             first.idx.data(),
                                                                             1,
                                                                             first.elem_ptr,
                                                                             second.n,
                                                                             second.ptr,
                                                                             second.idx.data(),
                                                                             second_nxe > 1 ? 1 : 0,
                                                                             second_nxe > 1 ? second.elem_ptr : nullptr,
                                                                             ccdptr.data(),
                                                                             cummax.data());
        PairSet out;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[first.n]), b(ccdptr[first.n]);
        sccd::collect_overlaps<first_nxe, second_nxe, scalar_t, idx_t>(axis,
                                                      first.n,
                                                      first.ptr,
                                                      first.idx.data(),
                                                      1,
                                                      first.elem_ptr,
                                                      second.n,
                                                      second.ptr,
                                                      second.idx.data(),
                                                      second_nxe > 1 ? 1 : 0,
                                                      second_nxe > 1 ? second.elem_ptr : nullptr,
                                                      ccdptr.data(),
                                                      cummax.data(),
                                                      a.data(),
                                                      b.data());
        for (size_t i = 0; i < a.size(); ++i) out.insert({a[i], b[i]});
        return out;
    }

    template <int first_nxe, int second_nxe = 1>
    PairSet cell2d_pairs(Boxes& first, Boxes& second) {
        sccd::Cell2DGrid<scalar_t> grid;
        sccd::cell2d_setup<scalar_t>(second.n, second.ptr, grid);

        sccd::Cell2DPartition part;
        sccd::cell2d_partition<scalar_t>(second.n, second.ptr, grid, part);

        std::vector<ptrdiff_t> cellptr(grid.ncells() + 1);
        sccd::cell2d_count<scalar_t>(second.n, second.ptr, grid, part, cellptr.data());

        std::vector<idx_t> cellidx(cellptr[grid.ncells()]);
        std::vector<ptrdiff_t> cursor(grid.ncells());
        std::vector<scalar_t> boxdata[6];
        scalar_t* cellbox[6];
        for (int d = 0; d < 6; ++d) {
            boxdata[d].resize((size_t)cellptr[grid.ncells()]);
            cellbox[d] = boxdata[d].data();
        }
        sccd::cell2d_fill<scalar_t, idx_t>(
            second.n, second.ptr, grid, part, cellptr.data(), cellidx.data(), cursor.data());
        sccd::cell2d_pack_boxes<scalar_t, idx_t>(grid, second.ptr, cellptr.data(), cellidx.data(), cellbox);

        std::vector<ptrdiff_t> ccdptr(first.n + 1, 0);
        const bool any = sccd::cell2d_count_overlaps<first_nxe, second_nxe, scalar_t, idx_t>(first.n,
                                                                            first.ptr,
                                                                            first.idx.data(),
                                                                            1,
                                                                            first.elem_ptr,
                                                                            second.ptr,
                                                                            second.idx.data(),
                                                                            second_nxe > 1 ? 1 : 0,
                                                                            second_nxe > 1 ? second.elem_ptr : nullptr,
                                                                            grid,
                                                                            cellptr.data(),
                                                                            cellidx.data(),
                                                            cellbox,
                                                                            ccdptr.data());
        PairSet out;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[first.n]), b(ccdptr[first.n]);
        sccd::cell2d_fill_overlaps<first_nxe, second_nxe, scalar_t, idx_t>(first.n,
                                                          first.ptr,
                                                          first.idx.data(),
                                                          1,
                                                          first.elem_ptr,
                                                          second.ptr,
                                                          second.idx.data(),
                                                          second_nxe > 1 ? 1 : 0,
                                                          second_nxe > 1 ? second.elem_ptr : nullptr,
                                                          grid,
                                                          cellptr.data(),
                                                          cellidx.data(),
                                                            cellbox,
                                                          ccdptr.data(),
                                                          a.data(),
                                                          b.data());
        for (size_t i = 0; i < a.size(); ++i) out.insert({a[i], b[i]});
        return out;
    }

    // force_axis pins the sort axis instead of letting choose_axis pick. The
    // sweep must return the same pair set whichever axis it sorts on, and that
    // invariant is what exposes a candidate window inconsistent with the overlap
    // predicate: the disagreement only appears when the degenerate axis is the
    // one being swept, which choose_axis will normally avoid.
    PairSet sweep_self_pairs(Boxes& e, const int force_axis = -1, ptrdiff_t* emitted = nullptr) {
        std::vector<scalar_t> scratch(e.n * 2);
        const int axis = force_axis >= 0 ? force_axis : sccd::choose_axis<scalar_t>(e.n, e.ptr);
        sccd::sort_along_axis(e.n, axis, e.ptr, e.idx.data(), scratch.data());

        std::vector<ptrdiff_t> ccdptr(e.n + 1, 0);
        const bool any = sccd::count_self_overlaps<2, scalar_t, idx_t>(
            axis, e.n, e.ptr, e.idx.data(), 1, e.elem_ptr, ccdptr.data());
        PairSet out;
        if (emitted) *emitted = 0;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[e.n]), b(ccdptr[e.n]);
        sccd::collect_self_overlaps<2, scalar_t, idx_t>(
            axis, e.n, e.ptr, e.idx.data(), 1, e.elem_ptr, ccdptr.data(), a.data(), b.data());
        for (size_t i = 0; i < a.size(); ++i) out.insert({a[i], b[i]});
        if (emitted) *emitted = (ptrdiff_t)a.size();
        return out;
    }

    PairSet cell2d_self_pairs(Boxes& e, ptrdiff_t* emitted = nullptr) {
        sccd::Cell2DGrid<scalar_t> grid;
        sccd::cell2d_setup<scalar_t>(e.n, e.ptr, grid);

        sccd::Cell2DPartition part;
        sccd::cell2d_partition<scalar_t>(e.n, e.ptr, grid, part);

        std::vector<ptrdiff_t> cellptr(grid.ncells() + 1);
        sccd::cell2d_count<scalar_t>(e.n, e.ptr, grid, part, cellptr.data());
        std::vector<idx_t> cellidx(cellptr[grid.ncells()]);
        std::vector<ptrdiff_t> cursor(grid.ncells());
        std::vector<scalar_t> boxdata[6];
        scalar_t* cellbox[6];
        for (int d = 0; d < 6; ++d) {
            boxdata[d].resize((size_t)cellptr[grid.ncells()]);
            cellbox[d] = boxdata[d].data();
        }
        sccd::cell2d_fill<scalar_t, idx_t>(
            e.n, e.ptr, grid, part, cellptr.data(), cellidx.data(), cursor.data());
        sccd::cell2d_pack_boxes<scalar_t, idx_t>(grid, e.ptr, cellptr.data(), cellidx.data(), cellbox);

        std::vector<ptrdiff_t> ccdptr(e.n + 1, 0);
        const bool any = sccd::cell2d_count_self_overlaps<2, scalar_t, idx_t>(
            e.n, e.ptr, e.idx.data(), 1, e.elem_ptr, grid, cellptr.data(), cellidx.data(),
                                                            cellbox, ccdptr.data());
        PairSet out;
        if (emitted) *emitted = 0;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[e.n]), b(ccdptr[e.n]);
        sccd::cell2d_fill_self_overlaps<2, scalar_t, idx_t>(e.n,
                                                            e.ptr,
                                                            e.idx.data(),
                                                            1,
                                                            e.elem_ptr,
                                                            grid,
                                                            cellptr.data(),
                                                            cellidx.data(),
                                                            cellbox,
                                                            ccdptr.data(),
                                                            a.data(),
                                                            b.data());
        for (size_t i = 0; i < a.size(); ++i) out.insert({a[i], b[i]});
        if (emitted) *emitted = (ptrdiff_t)a.size();
        return out;
    }

    // The same list against itself over a minimum-corner binning: one cell per
    // box, a forward walk in linear cell order, and the per-cell bounds pruning
    // it. It must return the sweep's pair set exactly, like the cell list does.
    template <bool sorted = false>
    PairSet cell2dmin_self_pairs(Boxes& e, ptrdiff_t* emitted = nullptr) {
        sccd::Cell2DGrid<scalar_t> grid;
        sccd::cell2d_setup<scalar_t>(e.n, e.ptr, grid);

        sccd::Cell2DPartition part;
        sccd::cell2dmin_partition<scalar_t>(e.n, e.ptr, grid, part);

        std::vector<ptrdiff_t> cellptr(grid.ncells() + 1);
        sccd::cell2dmin_count<scalar_t>(e.n, e.ptr, grid, part, cellptr.data());
        std::vector<idx_t> cellidx(cellptr[grid.ncells()]);
        std::vector<ptrdiff_t> cursor(grid.ncells());
        sccd::cell2dmin_fill<scalar_t, idx_t>(
            e.n, e.ptr, grid, part, cellptr.data(), cellidx.data(), cursor.data());

        std::vector<scalar_t> row_prefix(grid.ncells()), cell_hi1(grid.ncells());
        sccd::cell2dmin_bounds<scalar_t>(e.n, e.ptr, grid, part, row_prefix.data(), cell_hi1.data());

        std::vector<scalar_t> cell_key, cell_hi2;
        if constexpr (sorted) {
            cell_key.assign((size_t)cellptr[grid.ncells()], 0);
            cell_hi2.assign((size_t)grid.ncells(), 0);
            sccd::cell2dmin_sort_cells<scalar_t, idx_t>(
                e.ptr, grid, cellptr.data(), cellidx.data(), cell_key.data(), cell_hi2.data());
        }
        std::vector<scalar_t> boxdata[6];
        scalar_t* cellbox[6];
        for (int d = 0; d < 6; ++d) {
            boxdata[d].resize((size_t)cellptr[grid.ncells()]);
            cellbox[d] = boxdata[d].data();
        }
        sccd::cell2d_pack_boxes<scalar_t, idx_t>(grid, e.ptr, cellptr.data(), cellidx.data(), cellbox);

        const scalar_t* const key = sorted ? cell_key.data() : nullptr;
        const scalar_t* const hi2 = sorted ? cell_hi2.data() : nullptr;

        std::vector<ptrdiff_t> ccdptr(e.n + 1, 0);
        const bool any = sccd::cell2dmin_count_self_overlaps<2, sorted, scalar_t, idx_t>(e.n,
                                                                                 e.ptr,
                                                                                 e.idx.data(),
                                                                                 1,
                                                                                 e.elem_ptr,
                                                                                 grid,
                                                                                 cellptr.data(),
                                                                                 cellidx.data(),
                                                                                 cellbox,
                                                                                 row_prefix.data(),
                                                                                 cell_hi1.data(),
                                                                                 key,
                                                                                 hi2,
                                                                                 ccdptr.data());
        PairSet out;
        if (emitted) *emitted = 0;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[e.n]), b(ccdptr[e.n]);
        sccd::cell2dmin_fill_self_overlaps<2, sorted, scalar_t, idx_t>(e.n,
                                                               e.ptr,
                                                               e.idx.data(),
                                                               1,
                                                               e.elem_ptr,
                                                               grid,
                                                               cellptr.data(),
                                                               cellidx.data(),
                                                               cellbox,
                                                               row_prefix.data(),
                                                               cell_hi1.data(),
                                                               key,
                                                               hi2,
                                                               ccdptr.data(),
                                                               a.data(),
                                                               b.data());
        for (size_t i = 0; i < a.size(); ++i) out.insert({a[i], b[i]});
        if (emitted) *emitted = (ptrdiff_t)a.size();
        return out;
    }

    // Boxes that are flat on one axis and sit at a handful of repeated
    // coordinates, so that one box's xmax lands exactly on another's xmin.
    //
    // This is the case the random generator above will essentially never
    // produce and that real geometry produces constantly: an axis-aligned face
    // sweeps to a zero-extent AABB. The overlap predicate counts touching boxes
    // as overlapping, and the sweep's candidate window used a strict comparison
    // that skipped them, so it silently dropped real pairs -- 20 of 2220 on a
    // refined cube. A missed pair is a collision the narrow phase never sees,
    // which the conservativeness invariant does not allow.
    Boxes make_flat_boxes(std::mt19937& rng, const ptrdiff_t n, const int nxe, const int planes) {
        std::uniform_int_distribution<int> plane(0, planes - 1);
        std::uniform_real_distribution<double> pos(0.0, 10.0);
        std::uniform_int_distribution<int> node(0, (int)(n * nxe));

        Boxes b;
        for (int d = 0; d < 6; ++d) b.data[d].resize(n);
        b.idx.resize(n);
        for (int v = 0; v < nxe; ++v) b.elem[v].resize(n);
        b.n = n;

        for (ptrdiff_t i = 0; i < n; ++i) {
            // Flat on x, at one of a few shared planes: coincidence by design.
            const scalar_t x = (scalar_t)plane(rng);
            b.data[0][i] = x;
            b.data[3][i] = x;
            for (int d = 1; d < 3; ++d) {
                const scalar_t lo = (scalar_t)pos(rng);
                b.data[d][i] = lo;
                b.data[3 + d][i] = lo + (scalar_t)2.0;
            }
            b.idx[i] = (idx_t)i;
            for (int v = 0; v < nxe; ++v) b.elem[v][i] = (idx_t)node(rng);
        }
        b.bind();
        return b;
    }

    /**
     * \brief The partitioned binning must reproduce the serial binning exactly.
     *
     * Not "produce an equivalent cell list": the same bytes. A block owns a
     * contiguous range of cell rows and the boxes inside a block's list stay in
     * index order, so every cell sees its boxes in the order the serial scatter
     * would have written them. Comparing the arrays element by element is what
     * makes that claim testable, and it catches a partition that drops or
     * duplicates an incidence, which a pair-set comparison could still absorb.
     */
    int run_binning_case(const char* name, const ptrdiff_t n, const double spread, const double size) {
        std::mt19937 rng(20250913);
        Boxes b = make_boxes(rng, n, 2, spread, size);

        sccd::Cell2DGrid<scalar_t> grid;
        sccd::cell2d_setup<scalar_t>(b.n, b.ptr, grid);

        sccd::Cell2DPartition parallel_part;
        sccd::cell2d_partition<scalar_t>(b.n, b.ptr, grid, parallel_part);
        const sccd::Cell2DPartition serial_part;  // nblocks == 1: the serial path

        auto bin = [&](const sccd::Cell2DPartition& part,
                       std::vector<ptrdiff_t>& cellptr,
                       std::vector<idx_t>& cellidx) {
            cellptr.assign((size_t)grid.ncells() + 1, -1);
            sccd::cell2d_count<scalar_t>(b.n, b.ptr, grid, part, cellptr.data());
            cellidx.assign((size_t)cellptr[grid.ncells()], -1);
            std::vector<ptrdiff_t> cursor((size_t)grid.ncells());
            sccd::cell2d_fill<scalar_t, idx_t>(
                b.n, b.ptr, grid, part, cellptr.data(), cellidx.data(), cursor.data());
        };

        std::vector<ptrdiff_t> sptr, pptr;
        std::vector<idx_t> sidx, pidx;
        bin(serial_part, sptr, sidx);
        bin(parallel_part, pptr, pidx);

        // The minimum-corner binning and its bounds, through the same two paths.
        // This is the only case large enough to cross the parallel threshold, so
        // without it that path is never run.
        sccd::Cell2DPartition min_parallel_part;
        sccd::cell2dmin_partition<scalar_t>(b.n, b.ptr, grid, min_parallel_part);
        const sccd::Cell2DPartition min_serial_part;

        auto bin_min = [&](const sccd::Cell2DPartition& part,
                           std::vector<ptrdiff_t>& cellptr,
                           std::vector<idx_t>& cellidx,
                           std::vector<scalar_t>& row_prefix,
                           std::vector<scalar_t>& cell_hi1) {
            cellptr.assign((size_t)grid.ncells() + 1, -1);
            sccd::cell2dmin_count<scalar_t>(b.n, b.ptr, grid, part, cellptr.data());
            cellidx.assign((size_t)cellptr[grid.ncells()], -1);
            std::vector<ptrdiff_t> cursor((size_t)grid.ncells());
            sccd::cell2dmin_fill<scalar_t, idx_t>(
                b.n, b.ptr, grid, part, cellptr.data(), cellidx.data(), cursor.data());
            row_prefix.assign((size_t)grid.ncells(), 0);
            cell_hi1.assign((size_t)grid.ncells(), 0);
            sccd::cell2dmin_bounds<scalar_t>(b.n, b.ptr, grid, part, row_prefix.data(), cell_hi1.data());
        };

        std::vector<ptrdiff_t> msptr, mpptr;
        std::vector<idx_t> msidx, mpidx;
        std::vector<scalar_t> mspre, mppre, mshi, mphi;
        bin_min(min_serial_part, msptr, msidx, mspre, mshi);
        bin_min(min_parallel_part, mpptr, mpidx, mppre, mphi);

        const bool min_ok = (msptr == mpptr) && (msidx == mpidx) && (mspre == mppre) && (mshi == mphi);
        // One entry per box, which is the whole point of that binning.
        const bool min_entries_ok = msptr[grid.ncells()] == n;

        const bool ok = (sptr == pptr) && (sidx == pidx) && min_ok && min_entries_ok;
        std::printf("%-24s boxes=%-7ld cells=%-8ld spans=%-9ld mincorner=%-9ld blocks=%-4d  %s\n",
                    name, (long)n, (long)grid.ncells(), (long)sptr[grid.ncells()],
                    (long)msptr[grid.ncells()], parallel_part.nblocks, ok ? "ok" : "MISMATCH");
        if (!ok) {
            if (sptr != pptr) std::printf("    cell offsets differ\n");
            if (sidx != pidx) std::printf("    cell contents differ\n");
            if (msptr != mpptr) std::printf("    mincorner cell offsets differ\n");
            if (msidx != mpidx) std::printf("    mincorner cell contents differ\n");
            if (mspre != mppre) std::printf("    mincorner row prefix maxima differ\n");
            if (mshi != mphi) std::printf("    mincorner per-cell bounds differ\n");
            if (!min_entries_ok) {
                std::printf("    mincorner binned %ld entries for %ld boxes, expected one each\n",
                            (long)msptr[grid.ncells()], (long)n);
            }
        }
        // A case that never crosses into the partitioned path proves nothing
        // about it, so say so rather than passing quietly.
        if (parallel_part.serial() && n >= sccd::SCCD_CELL2D_MIN_PARALLEL) {
            std::printf("    note: ran serial anyway (one worker, or a grid one row deep)\n");
        }
        return ok ? 0 : 1;
    }

    /**
     * \brief sort_along_axis must produce exactly the order std::sort produces.
     *
     * The sweep's window walk assumes a total order on (key, index): equal keys
     * broken by index. The radix sort supplies the tie-break through stability
     * instead of through a comparator, so the two have to be checked against
     * each other directly. A pair-set comparison would not do it, because a
     * wrong order can still emit the right pairs on geometry that never puts two
     * boxes at the same coordinate.
     */
    int run_sort_case(const char* name, const ptrdiff_t n, const double spread, const double size) {
        std::mt19937 rng(31337);
        Boxes b = make_boxes(rng, n, 2, spread, size);

        struct KI { scalar_t key; idx_t idx; };
        int bad = 0;
        for (int axis = 0; axis < 3; ++axis) {
            std::vector<KI> want((size_t)n);
            for (ptrdiff_t i = 0; i < n; ++i) want[(size_t)i] = {b.data[axis][(size_t)i], (idx_t)i};
            std::sort(want.begin(), want.end(), [](const KI& l, const KI& r) {
                if (l.key < r.key) return true;
                if (r.key < l.key) return false;
                return l.idx < r.idx;
            });

            Boxes got = b;
            got.bind();
            std::vector<idx_t> idx((size_t)n);
            std::vector<scalar_t> scratch((size_t)n);
            sccd::sort_along_axis<scalar_t, idx_t>(n, axis, got.ptr, idx.data(), scratch.data());

            bool perm_ok = true, sorted_ok = true, rows_ok = true;
            for (ptrdiff_t i = 0; i < n; ++i) {
                if (idx[(size_t)i] != want[(size_t)i].idx) perm_ok = false;
                if (i && got.data[axis][(size_t)i] < got.data[axis][(size_t)(i - 1)]) sorted_ok = false;
                // Every one of the six rows must carry the same permutation.
                for (int d = 0; d < 6; ++d) {
                    if (got.data[d][(size_t)i] != b.data[d][(size_t)idx[(size_t)i]]) rows_ok = false;
                }
            }
            const bool ok = perm_ok && sorted_ok && rows_ok;
            std::printf("%-24s axis=%d boxes=%-7ld  %s%s%s\n", name, axis, (long)n,
                        ok ? "ok" : "MISMATCH",
                        perm_ok ? "" : " (permutation differs from std::sort)",
                        rows_ok ? "" : " (a row is permuted differently)");
            bad |= ok ? 0 : 1;
        }
        return bad;
    }

    // Emitted entries against distinct pairs. PairSet is a std::set, so equality
    // with the sweep proves no pair is missing and says nothing about a pair
    // arriving twice; the count pass is what makes a duplicate visible. The
    // edge-edge walk reads each cell once and each partner from one cell, so
    // the two numbers must agree exactly.
    int report_duplicates(const char* name, const char* which, const ptrdiff_t emitted,
                          const PairSet& pairs) {
        if (emitted == (ptrdiff_t)pairs.size()) return 0;
        std::printf("    %s: %s emitted %ld entries for %zu distinct pairs -- %ld DUPLICATE\n",
                    name, which, (long)emitted, pairs.size(), (long)(emitted - (ptrdiff_t)pairs.size()));
        return 1;
    }


    // ---------------------------------------------------------------------
    // Face-vertex with the vertex as the segment it really is.
    // ---------------------------------------------------------------------

    // Vertex trajectories, and the boxes that are their hulls. The kernel reads
    // the endpoints and the reference reads the same ones, so the two describe
    // one geometry.
    struct Moving {
        std::vector<scalar_t> p0[3], p1[3];
        scalar_t* p0_ptr[3];
        scalar_t* p1_ptr[3];

        void bind() {
            for (int d = 0; d < 3; ++d) {
                p0_ptr[d] = p0[d].data();
                p1_ptr[d] = p1[d].data();
            }
        }
    };

    Moving make_moving(std::mt19937& rng, Boxes& v, const double spread, const double travel) {
        std::uniform_real_distribution<double> pos(0.0, spread);
        std::uniform_real_distribution<double> step(-travel, travel);

        Moving m;
        for (int d = 0; d < 3; ++d) {
            m.p0[d].resize(v.n);
            m.p1[d].resize(v.n);
        }
        for (ptrdiff_t i = 0; i < v.n; ++i) {
            for (int d = 0; d < 3; ++d) {
                const scalar_t a = (scalar_t)pos(rng);
                const scalar_t b = a + (scalar_t)step(rng);
                m.p0[d][i] = a;
                m.p1[d][i] = b;
                v.data[d][i] = std::min(a, b);
                v.data[3 + d][i] = std::max(a, b);
            }
        }
        m.bind();
        v.bind();
        return m;
    }

    // Every (face, vertex) whose segment meets the face box, minus the faces the
    // vertex belongs to. The predicate the kernel must reproduce exactly.
    template <int nxe>
    PairSet brute_segment_pairs(const Boxes& f, const Boxes& v, const Moving& m) {
        PairSet out;
        for (ptrdiff_t i = 0; i < v.n; ++i) {
            scalar_t p0[3], d[3];
            for (int a = 0; a < 3; ++a) {
                p0[a] = m.p0[a][i];
                d[a] = m.p1[a][i] - p0[a];
            }
            for (ptrdiff_t j = 0; j < f.n; ++j) {
                const scalar_t bmin[3] = {f.data[0][j], f.data[1][j], f.data[2][j]};
                const scalar_t bmax[3] = {f.data[3][j], f.data[4][j], f.data[5][j]};
                scalar_t t = 0;
                scalar_t inv[3];
                for (int a = 0; a < 3; ++a) inv[a] = d[a] == 0 ? 0 : scalar_t(1) / d[a];
                if (!sccd::detail::segment_box_entry<scalar_t>(p0, d, inv, bmin, bmax, t)) continue;
                bool share = false;
                for (int a = 0; a < nxe; ++a) {
                    if (f.elem[a][j] == v.idx[i]) { share = true; break; }
                }
                if (share) continue;
                out.insert({f.idx[j], v.idx[i]});
            }
        }
        return out;
    }

    template <int nxe>
    PairSet cell2dseg_pairs(Boxes& f, Boxes& v, Moving& m, ptrdiff_t* emitted = nullptr) {
        sccd::Cell2DGrid<scalar_t> grid;
        sccd::cell2d_setup<scalar_t>(f.n, f.ptr, grid);

        sccd::Cell2DPartition part;
        sccd::cell2d_partition<scalar_t>(f.n, f.ptr, grid, part);

        std::vector<ptrdiff_t> cellptr(grid.ncells() + 1);
        sccd::cell2d_count<scalar_t>(f.n, f.ptr, grid, part, cellptr.data());
        std::vector<idx_t> cellidx(cellptr[grid.ncells()]);
        std::vector<ptrdiff_t> cursor(grid.ncells());
        std::vector<scalar_t> boxdata[6];
        scalar_t* cellbox[6];
        for (int d = 0; d < 6; ++d) {
            boxdata[d].resize((size_t)cellptr[grid.ncells()]);
            cellbox[d] = boxdata[d].data();
        }
        sccd::cell2d_fill<scalar_t, idx_t>(
            f.n, f.ptr, grid, part, cellptr.data(), cellidx.data(), cursor.data());
        sccd::cell2d_pack_boxes<scalar_t, idx_t>(grid, f.ptr, cellptr.data(), cellidx.data(), cellbox);

        std::vector<ptrdiff_t> ccdptr(v.n + 1, 0);
        const bool any = sccd::cell2dseg_count_vf_overlaps<nxe, scalar_t, idx_t>(
            v.n, m.p0_ptr, m.p1_ptr, v.idx.data(), f.idx.data(), 1, f.elem_ptr,
            grid, cellptr.data(), cellidx.data(), cellbox, ccdptr.data());

        PairSet out;
        if (emitted) *emitted = 0;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[v.n]), b(ccdptr[v.n]);
        sccd::cell2dseg_fill_vf_overlaps<nxe, scalar_t, idx_t>(
            v.n, m.p0_ptr, m.p1_ptr, v.idx.data(), f.idx.data(), 1, f.elem_ptr,
            grid, cellptr.data(), cellidx.data(), cellbox, ccdptr.data(), a.data(), b.data());
        for (size_t k = 0; k < a.size(); ++k) out.insert({a[k], b[k]});
        if (emitted) *emitted = (ptrdiff_t)a.size();
        return out;
    }

    // The segment query against brute force, and against the box query it would
    // replace: it must be a subset of the boxes' answer, since a segment meeting
    // a box implies the segment's own box meeting it.
    template <int nxe>
    int run_segment_case(const char* name, const ptrdiff_t n_faces, const ptrdiff_t n_verts,
                         const double spread, const double face_size, const double travel) {
        std::mt19937 rng(9001);
        Boxes f = make_boxes(rng, n_faces, nxe, spread, face_size);
        Boxes v = make_boxes(rng, n_verts, 1, spread, 0.0);
        Moving m = make_moving(rng, v, spread, travel);

        Boxes f_q = f, v_q = v;
        f_q.bind();
        v_q.bind();

        ptrdiff_t emitted = 0;
        const PairSet got = cell2dseg_pairs<nxe>(f_q, v_q, m, &emitted);
        const PairSet want = brute_segment_pairs<nxe>(f, v, m);

        Boxes f_b = f, v_b = v;
        f_b.bind();
        v_b.bind();
        const PairSet boxes = cell2d_pairs<nxe, 1>(f_b, v_b);

        std::vector<std::pair<idx_t, idx_t>> extra;
        std::set_difference(got.begin(), got.end(), boxes.begin(), boxes.end(),
                            std::back_inserter(extra));

        const bool dup = emitted != (ptrdiff_t)got.size();
        const bool ok = (got == want) && extra.empty() && !dup;
        std::printf("%-26s faces=%-6ld verts=%-6ld travel=%-5.1f boxes=%-8zu segment=%-8zu "
                    "tighter=%5.1f%%  %s\n",
                    name, (long)n_faces, (long)n_verts, travel, boxes.size(), got.size(),
                    boxes.size() ? 100.0 * (1.0 - (double)got.size() / boxes.size()) : 0.0,
                    ok ? "ok" : "MISMATCH");
        if (dup) {
            std::printf("    emitted %ld entries for %zu distinct pairs -- %ld DUPLICATE\n",
                        (long)emitted, got.size(), (long)(emitted - (ptrdiff_t)got.size()));
        }
        if (got != want) {
            std::vector<std::pair<idx_t, idx_t>> missed;
            std::set_difference(want.begin(), want.end(), got.begin(), got.end(),
                                std::back_inserter(missed));
            std::printf("    MISSED %zu pairs brute force found\n", missed.size());
        }
        if (!extra.empty()) {
            std::printf("    emitted %zu pairs the box query does not -- not a subset\n", extra.size());
        }
        return ok ? 0 : 1;
    }

    int run_flat_self_case(const char* name, const ptrdiff_t n, const int planes) {
        std::mt19937 rng(4242);
        Boxes e = make_flat_boxes(rng, n, 2, planes);
        Boxes e_c = e, e_m = e;
        e_c.bind();
        e_m.bind();

        Boxes e_s = e;
        e_s.bind();
        ptrdiff_t n_cell = 0, n_min = 0, n_sort = 0;
        const PairSet cell = cell2d_self_pairs(e_c, &n_cell);
        const PairSet mincorner = cell2dmin_self_pairs(e_m, &n_min);
        const PairSet minsorted = cell2dmin_self_pairs<true>(e_s, &n_sort);

        int bad = 0;
        bad |= report_duplicates(name, "cell list", n_cell, cell);
        bad |= report_duplicates(name, "mincorner", n_min, mincorner);
        bad |= report_duplicates(name, "sorted", n_sort, minsorted);

        for (int axis = 0; axis < 3; ++axis) {
            Boxes e_axis = e;
            e_axis.bind();
            ptrdiff_t n_sweep = 0;
            const PairSet sweep = sweep_self_pairs(e_axis, axis, &n_sweep);
            bad |= report_duplicates(name, "sweep", n_sweep, sweep);

            std::vector<std::pair<idx_t, idx_t>> only_cell;
            std::set_difference(cell.begin(), cell.end(), sweep.begin(), sweep.end(),
                                std::back_inserter(only_cell));

            const bool ok = (cell == sweep) && (mincorner == sweep) && (minsorted == sweep);
            std::printf("%-24s axis=%d boxes=%-6ld planes=%-3d sweep=%-8zu cell=%-8zu mincorner=%-8zu sorted=%-8zu  %s\n",
                        name, axis, (long)n, planes, sweep.size(), cell.size(), mincorner.size(),
                        minsorted.size(), ok ? "ok" : "MISMATCH");
            if (!only_cell.empty()) {
                std::printf("    the sweep MISSED %zu pairs the cell list found -- "
                            "a missed pair is a collision the narrow phase never sees\n",
                            only_cell.size());
            }
            bad |= ok ? 0 : 1;
        }
        return bad;
    }

    int run_self_boxes(const char* name, Boxes& e, const ptrdiff_t n) {
        Boxes e_c = e, e_m = e, e_s = e;
        e_c.bind();
        e_m.bind();
        e_s.bind();

        ptrdiff_t n_cell = 0, n_min = 0, n_sort = 0, n_sweep = 0;
        const PairSet cell = cell2d_self_pairs(e_c, &n_cell);
        const PairSet mincorner = cell2dmin_self_pairs(e_m, &n_min);
        const PairSet minsorted = cell2dmin_self_pairs<true>(e_s, &n_sort);
        const PairSet sweep = sweep_self_pairs(e, -1, &n_sweep);

        int dup = 0;
        dup |= report_duplicates(name, "cell list", n_cell, cell);
        dup |= report_duplicates(name, "mincorner", n_min, mincorner);
        dup |= report_duplicates(name, "sorted", n_sort, minsorted);
        dup |= report_duplicates(name, "sweep", n_sweep, sweep);

        const bool ok = (cell == sweep) && (mincorner == sweep) && (minsorted == sweep) && !dup;
        std::printf("%-26s edges=%-6ld sweep=%-8zu cell=%-8zu mincorner=%-8zu sorted=%-8zu  %s\n",
                    name,
                    (long)n,
                    sweep.size(),
                    cell.size(),
                    mincorner.size(),
                    minsorted.size(),
                    ok ? "ok" : "MISMATCH");
        if (!ok) {
            const PairSet& bad =
                (cell != sweep) ? cell : ((mincorner != sweep) ? mincorner : minsorted);
            const char* which =
                (cell != sweep) ? "cell list" : ((mincorner != sweep) ? "mincorner" : "sorted");
            std::vector<std::pair<idx_t, idx_t>> only_sweep, only_bad;
            std::set_difference(
                sweep.begin(), sweep.end(), bad.begin(), bad.end(), std::back_inserter(only_sweep));
            std::set_difference(
                bad.begin(), bad.end(), sweep.begin(), sweep.end(), std::back_inserter(only_bad));
            std::printf("    missed by %s: %zu   extra: %zu\n", which, only_sweep.size(), only_bad.size());
            for (size_t i = 0; i < only_sweep.size() && i < 5; ++i) {
                std::printf("    missing (%d,%d)\n", only_sweep[i].first, only_sweep[i].second);
            }
        }
        return ok ? 0 : 1;
    }

    int run_self_case(const char* name, const ptrdiff_t n, const double spread, const double size) {
        std::mt19937 rng(999);
        Boxes e = make_boxes(rng, n, 2, spread, size);
        return run_self_boxes(name, e, n);
    }

    // Most boxes small, a handful enormous.
    //
    // This is what a swept AABB set looks like when a few elements move fast,
    // and it is where the duplicate rule earns its keep: a box a hundred times
    // the mean covers hundreds of cells, meets each partner in as many of them
    // as the two share, and must still emit each pair once. The uniform
    // generator above never produces that, and the real scenes do -- on
    // armadillo-rollers the widest swept edge is 7.5 times the mean.
    Boxes make_spread_boxes(std::mt19937& rng,
                            const ptrdiff_t n,
                            const int nxe,
                            const double spread,
                            const double small,
                            const double big,
                            const ptrdiff_t n_big) {
        Boxes b = make_boxes(rng, n, nxe, spread, small);
        std::uniform_real_distribution<double> pos(0.0, spread);
        for (ptrdiff_t k = 0; k < n_big && k < n; ++k) {
            const ptrdiff_t i = (k * n) / (n_big > 0 ? n_big : 1);
            for (int d = 0; d < 3; ++d) {
                const double lo = pos(rng);
                b.data[d][i] = lo;
                b.data[3 + d][i] = lo + big;
            }
        }
        b.bind();
        return b;
    }

    int run_spread_self_case(const char* name,
                             const ptrdiff_t n,
                             const double spread,
                             const double small,
                             const double big,
                             const ptrdiff_t n_big) {
        std::mt19937 rng(777);
        Boxes e = make_spread_boxes(rng, n, 2, spread, small, big, n_big);
        return run_self_boxes(name, e, n);
    }

    int run_spread_case(const char* name,
                        const ptrdiff_t nf,
                        const ptrdiff_t nv,
                        const double spread,
                        const double small,
                        const double big,
                        const ptrdiff_t n_big) {
        std::mt19937 rng(31337);
        // The spread is on both lists, so wide faces walk many cells and wide
        // vertices sit in many.
        Boxes faces = make_spread_boxes(rng, nf, 3, spread, small, big, n_big);
        Boxes verts = make_spread_boxes(rng, nv, 1, spread, small * 0.1, big * 0.5, n_big);

        Boxes faces_c = faces, verts_c = verts;
        faces_c.bind();
        verts_c.bind();

        const PairSet cell = cell2d_pairs<3>(faces_c, verts_c);
        const PairSet sweep = sweep_pairs<3>(faces, verts);

        const bool ok = (cell == sweep);
        std::printf("%-28s faces=%-7ld verts=%-7ld sweep=%-8zu cell=%-8zu  %s\n",
                    name, (long)nf, (long)nv, sweep.size(), cell.size(), ok ? "ok" : "MISMATCH");
        if (!ok) {
            std::vector<std::pair<idx_t, idx_t>> only_sweep;
            std::set_difference(
                sweep.begin(), sweep.end(), cell.begin(), cell.end(), std::back_inserter(only_sweep));
            std::printf("    missed by cell list: %zu\n", only_sweep.size());
            for (size_t i = 0; i < only_sweep.size() && i < 5; ++i) {
                std::printf("    missing (%d,%d)\n", only_sweep[i].first, only_sweep[i].second);
            }
        }
        return ok ? 0 : 1;
    }

    template <int nxe>
    int run_case(const char* name, const ptrdiff_t nf, const ptrdiff_t nv, const double spread, const double size) {
        std::mt19937 rng(12345);
        Boxes faces = make_boxes(rng, nf, nxe, spread, size);
        Boxes verts = make_boxes(rng, nv, 1, spread, size * 0.1);

        // Both take sorted input for the sweep, so copy before it reorders them.
        Boxes faces_c = faces, verts_c = verts;
        faces_c.bind();
        verts_c.bind();

        const PairSet cell = cell2d_pairs<nxe>(faces_c, verts_c);
        const PairSet sweep = sweep_pairs<nxe>(faces, verts);

        std::vector<std::pair<idx_t, idx_t>> only_sweep, only_cell;
        std::set_difference(sweep.begin(), sweep.end(), cell.begin(), cell.end(), std::back_inserter(only_sweep));
        std::set_difference(cell.begin(), cell.end(), sweep.begin(), sweep.end(), std::back_inserter(only_cell));

        const bool ok = only_sweep.empty() && only_cell.empty();
        std::printf("%-20s nxe=%d faces=%-7ld verts=%-7ld sweep=%-8zu cell=%-8zu  %s\n",
                    name,
                    nxe,
                    (long)nf,
                    (long)nv,
                    sweep.size(),
                    cell.size(),
                    ok ? "ok" : "MISMATCH");
        if (!ok) {
            std::printf("    missed by cell list: %zu   extra in cell list: %zu\n",
                        only_sweep.size(),
                        only_cell.size());
            for (size_t i = 0; i < only_sweep.size() && i < 5; ++i) {
                std::printf("    missing (%d,%d)\n", only_sweep[i].first, only_sweep[i].second);
            }
            for (size_t i = 0; i < only_cell.size() && i < 5; ++i) {
                std::printf("    extra   (%d,%d)\n", only_cell[i].first, only_cell[i].second);
            }
        }
        return ok ? 0 : 1;
    }

    // Element against element.
    //
    // Every caller in the repository passes second_nxe == 1, which leaves the
    // `if constexpr (second_nxe > 1)` branch of the shared-node masking
    // unreachable -- so the broad phase's generality over the second list's
    // element width was carried by the signature and by nothing else. This runs
    // it, against a reference that does the masking independently.
    template <int first_nxe, int second_nxe>
    int run_two_element_case(const char* name, const ptrdiff_t na, const ptrdiff_t nb) {
        std::mt19937 rng(9001);
        Boxes a = make_boxes(rng, na, first_nxe, 20.0, 2.0);
        Boxes b = make_boxes(rng, nb, second_nxe, 20.0, 2.0);
        // A small node pool, so shared nodes between the two lists are common
        // and the masking actually has work to do.
        for (ptrdiff_t i = 0; i < a.n; ++i)
            for (int v = 0; v < first_nxe; ++v) a.elem[v][i] %= 64;
        for (ptrdiff_t i = 0; i < b.n; ++i)
            for (int v = 0; v < second_nxe; ++v) b.elem[v][i] %= 64;
        a.bind();
        b.bind();

        ptrdiff_t masked = 0;
        const PairSet expected = brute_pairs<first_nxe, second_nxe>(a, b, &masked);

        Boxes a_c = a, b_c = b;
        a_c.bind();
        b_c.bind();
        const PairSet cell = cell2d_pairs<first_nxe, second_nxe>(a_c, b_c);
        const PairSet sweep = sweep_pairs<first_nxe, second_nxe>(a, b);

        // If the masking removed nothing, the branch ran but proved nothing --
        // the scene has to contain shared nodes for this to be a test.
        const bool ok = (sweep == expected) && (cell == expected) && (masked > 0);
        std::printf("%-22s <%d,%d> a=%-6ld b=%-6ld brute=%-7zu masked=%-6ld sweep=%-7zu cell=%-7zu  %s\n",
                    name, first_nxe, second_nxe, (long)na, (long)nb,
                    expected.size(), (long)masked, sweep.size(), cell.size(),
                    ok ? "ok" : "MISMATCH");
        return ok ? 0 : 1;
    }

}  // namespace

int main() {
    int bad = 0;
    // Every case runs at both face widths. Triangles and quads differ only in how
    // many of a face's own vertices the collector must skip, but that is exactly
    // the step where the cell list and the sweep could disagree, and quads are a
    // supported topology -- so they are checked, not assumed.
    for (int pass = 0; pass < 2; ++pass) {
        const bool quad = (pass == 1);
        // Ordinary: small boxes, most spanning one cell.
        bad |= quad ? run_case<4>("small boxes", 2000, 4000, 100.0, 1.0)
                    : run_case<3>("small boxes", 2000, 4000, 100.0, 1.0);
        // Boxes spanning many cells: the duplicate-suppression case.
        bad |= quad ? run_case<4>("spanning many cells", 1500, 3000, 100.0, 25.0)
                    : run_case<3>("spanning many cells", 1500, 3000, 100.0, 25.0);
        // One box covering the whole domain.
        bad |= quad ? run_case<4>("very large boxes", 800, 1500, 10.0, 10.0)
                    : run_case<3>("very large boxes", 800, 1500, 10.0, 10.0);
        // Dense: everything overlaps everything.
        bad |= quad ? run_case<4>("dense overlap", 600, 1200, 1.0, 1.0)
                    : run_case<3>("dense overlap", 600, 1200, 1.0, 1.0);
        // Degenerate: zero-extent boxes, all coincident on two axes.
        bad |= quad ? run_case<4>("tiny spread", 500, 1000, 0.0001, 0.00001)
                    : run_case<3>("tiny spread", 500, 1000, 0.0001, 0.00001);
    }

    // Element against element, both lists carrying connectivity. Nothing in the
    // library calls the broad phase this way today; the capability is real and
    // was untested, which is how it came to look removable.
    bad |= run_two_element_case<3, 2>("tri vs edge", 1200, 2000);
    bad |= run_two_element_case<3, 3>("tri vs tri", 1000, 1000);
    bad |= run_two_element_case<4, 2>("quad vs edge", 900, 1500);

    // sort_along_axis at sizes on both sides of the radix cutoff, so both the
    // radix path and the comparison path are checked. "ties" puts many boxes at
    // the same coordinate, which is where the tie-break rule shows.
    bad |= run_sort_case("sort: random", 50000, 100.0, 1.0);
    bad |= run_sort_case("sort: ties", 50000, 0.001, 0.0001);
    bad |= run_sort_case("sort: below cutoff", 1000, 100.0, 1.0);

    // The partitioned binning against the serial binning, at sizes that cross
    // the per-block minimum so the partitioned path actually runs.
    bad |= run_binning_case("binning: one cell each", 60000, 100.0, 0.5);
    bad |= run_binning_case("binning: many cells", 60000, 100.0, 6.0);
    bad |= run_binning_case("binning: dense", 40000, 1.0, 1.0);
    bad |= run_binning_case("binning: degenerate", 40000, 0.0001, 0.00001);

    // Self-overlap: edge-edge, where each unordered pair must appear once.
    bad |= run_self_case("self: small boxes", 3000, 100.0, 1.0);
    bad |= run_self_case("self: many cells", 2000, 100.0, 25.0);
    bad |= run_self_case("self: dense", 900, 1.0, 1.0);

    // Degenerate: flat boxes sharing exact coordinates, where one box's xmax is
    // another's xmin. Regression for a sweep that dropped touching pairs.
    bad |= run_flat_self_case("self: flat, coincident", 2000, 4);
    bad |= run_flat_self_case("self: flat, one plane", 500, 1);

    // A size spread, where the duplicate rule does real work: a box far wider
    // than a cell covers many of them and shares several with each partner.
    bad |= run_spread_self_case("spread: one outlier", 4000, 100.0, 1.0, 60.0, 1);
    bad |= run_spread_self_case("spread: 1 in 200", 4000, 100.0, 0.5, 40.0, 20);
    bad |= run_spread_self_case("spread: heavy tail", 3000, 100.0, 0.5, 20.0, 300);
    bad |= run_spread_self_case("spread: dense and mixed", 900, 4.0, 0.5, 4.0, 40);
    bad |= run_spread_case("spread: faces vs verts", 900, 2000, 60.0, 1.0, 30.0, 12);

    // Face-vertex queried with the vertex trajectory. The travel is what decides
    // whether the segment is worth anything: at rest it is its own box, and at
    // several cells it is a diagonal through a rectangle it barely touches.
    bad |= run_segment_case<3>("segment: short travel", 3000, 4000, 100.0, 2.0, 0.5);
    bad |= run_segment_case<3>("segment: one cell", 3000, 4000, 100.0, 2.0, 2.0);
    bad |= run_segment_case<3>("segment: far travel", 2000, 3000, 100.0, 2.0, 20.0);
    bad |= run_segment_case<3>("segment: crossing", 1500, 2000, 100.0, 2.0, 100.0);
    bad |= run_segment_case<4>("segment: quads", 2000, 3000, 100.0, 2.0, 10.0);
    bad |= run_segment_case<3>("segment: dense", 800, 1200, 4.0, 1.0, 4.0);
    bad |= run_spread_case("spread: wide queries", 700, 1800, 40.0, 0.5, 25.0, 120);

    std::printf("%s\n", bad ? "FAIL" : "OK: cell list and sweep agree on every case");
    return bad;
}
