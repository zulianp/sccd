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

    template <int first_nxe>
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
        const bool any = sccd::count_overlaps<first_nxe, 1, scalar_t, idx_t>(axis,
                                                                             first.n,
                                                                             first.ptr,
                                                                             first.idx.data(),
                                                                             1,
                                                                             first.elem_ptr,
                                                                             second.n,
                                                                             second.ptr,
                                                                             second.idx.data(),
                                                                             0,
                                                                             nullptr,
                                                                             ccdptr.data(),
                                                                             cummax.data());
        PairSet out;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[first.n]), b(ccdptr[first.n]);
        sccd::collect_overlaps<first_nxe, 1, scalar_t, idx_t>(axis,
                                                      first.n,
                                                      first.ptr,
                                                      first.idx.data(),
                                                      1,
                                                      first.elem_ptr,
                                                      second.n,
                                                      second.ptr,
                                                      second.idx.data(),
                                                      0,
                                                      nullptr,
                                                      ccdptr.data(),
                                                      cummax.data(),
                                                      a.data(),
                                                      b.data());
        for (size_t i = 0; i < a.size(); ++i) out.insert({a[i], b[i]});
        return out;
    }

    template <int first_nxe>
    PairSet cell2d_pairs(Boxes& first, Boxes& second) {
        sccd::Cell2DGrid<scalar_t> grid;
        sccd::cell2d_setup<scalar_t>(second.n, second.ptr, grid);

        std::vector<ptrdiff_t> cellptr(grid.ncells() + 1);
        sccd::cell2d_count<scalar_t>(second.n, second.ptr, grid, cellptr.data());

        std::vector<idx_t> cellidx(cellptr[grid.ncells()]);
        std::vector<ptrdiff_t> cursor(grid.ncells());
        sccd::cell2d_fill<scalar_t, idx_t>(
            second.n, second.ptr, grid, cellptr.data(), cellidx.data(), cursor.data());

        std::vector<ptrdiff_t> ccdptr(first.n + 1, 0);
        const bool any = sccd::cell2d_count_overlaps<first_nxe, 1, scalar_t, idx_t>(first.n,
                                                                            first.ptr,
                                                                            first.idx.data(),
                                                                            1,
                                                                            first.elem_ptr,
                                                                            second.ptr,
                                                                            second.idx.data(),
                                                                            0,
                                                                            nullptr,
                                                                            grid,
                                                                            cellptr.data(),
                                                                            cellidx.data(),
                                                                            ccdptr.data());
        PairSet out;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[first.n]), b(ccdptr[first.n]);
        sccd::cell2d_fill_overlaps<first_nxe, 1, scalar_t, idx_t>(first.n,
                                                          first.ptr,
                                                          first.idx.data(),
                                                          1,
                                                          first.elem_ptr,
                                                          second.ptr,
                                                          second.idx.data(),
                                                          0,
                                                          nullptr,
                                                          grid,
                                                          cellptr.data(),
                                                          cellidx.data(),
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
    PairSet sweep_self_pairs(Boxes& e, const int force_axis = -1) {
        std::vector<scalar_t> scratch(e.n * 2);
        const int axis = force_axis >= 0 ? force_axis : sccd::choose_axis<scalar_t>(e.n, e.ptr);
        sccd::sort_along_axis(e.n, axis, e.ptr, e.idx.data(), scratch.data());

        std::vector<ptrdiff_t> ccdptr(e.n + 1, 0);
        const bool any = sccd::count_self_overlaps<2, scalar_t, idx_t>(
            axis, e.n, e.ptr, e.idx.data(), 1, e.elem_ptr, ccdptr.data());
        PairSet out;
        if (!any) return out;

        std::vector<idx_t> a(ccdptr[e.n]), b(ccdptr[e.n]);
        sccd::collect_self_overlaps<2, scalar_t, idx_t>(
            axis, e.n, e.ptr, e.idx.data(), 1, e.elem_ptr, ccdptr.data(), a.data(), b.data());
        for (size_t i = 0; i < a.size(); ++i) out.insert({a[i], b[i]});
        return out;
    }

    PairSet cell2d_self_pairs(Boxes& e) {
        sccd::Cell2DGrid<scalar_t> grid;
        sccd::cell2d_setup<scalar_t>(e.n, e.ptr, grid);

        std::vector<ptrdiff_t> cellptr(grid.ncells() + 1);
        sccd::cell2d_count<scalar_t>(e.n, e.ptr, grid, cellptr.data());
        std::vector<idx_t> cellidx(cellptr[grid.ncells()]);
        std::vector<ptrdiff_t> cursor(grid.ncells());
        sccd::cell2d_fill<scalar_t, idx_t>(e.n, e.ptr, grid, cellptr.data(), cellidx.data(), cursor.data());

        std::vector<ptrdiff_t> ccdptr(e.n + 1, 0);
        const bool any = sccd::cell2d_count_self_overlaps<2, scalar_t, idx_t>(
            e.n, e.ptr, e.idx.data(), 1, e.elem_ptr, grid, cellptr.data(), cellidx.data(), ccdptr.data());
        PairSet out;
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
                                                            ccdptr.data(),
                                                            a.data(),
                                                            b.data());
        for (size_t i = 0; i < a.size(); ++i) out.insert({a[i], b[i]});
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

    int run_flat_self_case(const char* name, const ptrdiff_t n, const int planes) {
        std::mt19937 rng(4242);
        Boxes e = make_flat_boxes(rng, n, 2, planes);
        Boxes e_c = e;
        e_c.bind();

        const PairSet cell = cell2d_self_pairs(e_c);

        int bad = 0;
        for (int axis = 0; axis < 3; ++axis) {
            Boxes e_axis = e;
            e_axis.bind();
            const PairSet sweep = sweep_self_pairs(e_axis, axis);

            std::vector<std::pair<idx_t, idx_t>> only_cell;
            std::set_difference(cell.begin(), cell.end(), sweep.begin(), sweep.end(),
                                std::back_inserter(only_cell));

            const bool ok = (cell == sweep);
            std::printf("%-24s axis=%d boxes=%-6ld planes=%-3d sweep=%-8zu cell=%-8zu  %s\n",
                        name, axis, (long)n, planes, sweep.size(), cell.size(),
                        ok ? "ok" : "MISMATCH");
            if (!only_cell.empty()) {
                std::printf("    the sweep MISSED %zu pairs the cell list found -- "
                            "a missed pair is a collision the narrow phase never sees\n",
                            only_cell.size());
            }
            bad |= ok ? 0 : 1;
        }
        return bad;
    }

    int run_self_case(const char* name, const ptrdiff_t n, const double spread, const double size) {
        std::mt19937 rng(999);
        Boxes e = make_boxes(rng, n, 2, spread, size);
        Boxes e_c = e;
        e_c.bind();

        const PairSet cell = cell2d_self_pairs(e_c);
        const PairSet sweep = sweep_self_pairs(e);

        const bool ok = (cell == sweep);
        std::printf("%-28s edges=%-7ld            sweep=%-8zu cell=%-8zu  %s\n",
                    name,
                    (long)n,
                    sweep.size(),
                    cell.size(),
                    ok ? "ok" : "MISMATCH");
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

    // Self-overlap: edge-edge, where each unordered pair must appear once.
    bad |= run_self_case("self: small boxes", 3000, 100.0, 1.0);
    bad |= run_self_case("self: many cells", 2000, 100.0, 25.0);
    bad |= run_self_case("self: dense", 900, 1.0, 1.0);

    // Degenerate: flat boxes sharing exact coordinates, where one box's xmax is
    // another's xmin. Regression for a sweep that dropped touching pairs.
    bad |= run_flat_self_case("self: flat, coincident", 2000, 4);
    bad |= run_flat_self_case("self: flat, one plane", 500, 1);

    std::printf("%s\n", bad ? "FAIL" : "OK: cell list and sweep agree on every case");
    return bad;
}
