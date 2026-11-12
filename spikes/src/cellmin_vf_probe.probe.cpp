// Can the vertex-face query bin at the minimum corner the way the edge-edge one
// does, and what does it cost?
//
// The shipped vertex-face query bins the VERTICES by extent and iterates over
// the FACES. `cell2d_setup` sizes the cells from the list it bins, so the cells
// are vertex-sized and a face -- the larger of the two -- walks a rectangle of
// them, paying the duplicate rule (the cell holding the minimum corner of the
// overlap) on every surviving candidate.
//
// The proposal here is the mirror image, and it is the one `Cell2DSeg` already
// showed has a cheaper query: bin the FACES and iterate over the VERTICES. What
// sank `Cell2DSeg` was not the query but the structure -- extent-binning faces
// costs about 4.3 entries each over twice as many elements, and the extra
// preparation ate the saving on three scenes of five. Binning the faces at their
// MINIMUM CORNER costs one entry each, which is the whole point of measuring
// this again.
//
// Completeness, which is the part that differs from the edge-edge walk. There is
// no self-symmetry to exploit: vertices never walk, so a vertex must find every
// face wherever it is binned. Two bounds make that finite.
//
//   A face f can overlap a vertex v only if f.min <= v.max componentwise, so
//   cell(f.min) <= cell(v.max) in BOTH axes. Nothing right of col(v.max) or
//   above row(v.max) is ever read. That is exact, not a heuristic.
//
//   Downwards and leftwards the walk is bounded by what the cells hold. Per row,
//   a prefix maximum of f.max on the first axis gives, by one binary search, the
//   first column holding anything that reaches back to v -- the same array and
//   the same search the edge-edge walk uses. Per row range, K1 is the largest
//   number of rows any face spans, so a face binned below row(v.max) - K1 cannot
//   reach v.
//
// No duplicate rule is needed on either side: a face has exactly one cell, so a
// pair is met in exactly one place.
//
// Usage: cellmin_vf_probe [faces] [verts] [face extent] [vertex extent]
// Exit status is non-zero if either scheme misses a pair brute force found.
// Written up in wip/DECISIONS.md.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <random>
#include <set>
#include <vector>

struct Box { double lo[2], hi[2]; };

static bool overlap(const Box& a, const Box& b) {
    for (int d = 0; d < 2; ++d)
        if (a.lo[d] > b.hi[d] || b.lo[d] > a.hi[d]) return false;
    return true;
}

// A grid sized the way cell2d_setup sizes one: the cell side is the mean extent
// of the list being binned, so the resolution follows whichever list is indexed.
struct Grid {
    double cell[2];
    int n[2];
    Grid(const std::vector<Box>& b, double span) {
        double mean[2] = {0, 0};
        for (const Box& x : b) for (int d = 0; d < 2; ++d) mean[d] += x.hi[d] - x.lo[d];
        for (int d = 0; d < 2; ++d) {
            mean[d] = mean[d] / (double)b.size();
            if (mean[d] <= 0) mean[d] = span;
            cell[d] = mean[d];
            n[d] = (int)(span / cell[d]) + 1;
        }
        // cell2d_setup caps the cell count at 4n and coarsens uniformly to get
        // there. Without this the grid is far finer than the shipped one and the
        // occupancies come out about twice what the real scenes show.
        const double cap = 4.0 * (double)b.size();
        double have = (double)n[0] * (double)n[1];
        if (have > cap) {
            const double s = std::sqrt(cap / have);
            for (int d = 0; d < 2; ++d) {
                n[d] = (int)((double)n[d] * s);
                if (n[d] < 1) n[d] = 1;
                cell[d] = span / (double)n[d];
            }
        }
    }
    int c0(double v) const { int c = (int)(v / cell[0]); return c < 0 ? 0 : (c >= n[0] ? n[0] - 1 : c); }
    int c1(double v) const { int c = (int)(v / cell[1]); return c < 0 ? 0 : (c >= n[1] ? n[1] - 1 : c); }
    size_t of(int a, int b) const { return (size_t)b * n[0] + a; }
    size_t ncells() const { return (size_t)n[0] * n[1]; }
};

int main(int argc, char** argv) {
    const int nf = argc > 1 ? atoi(argv[1]) : 8000;
    const int nv = argc > 2 ? atoi(argv[2]) : 4000;
    const double span = 100.0;
    // Calibrated against the occupancies recorded in wip/DECISIONS.md: on the
    // vertex-sized grid a vertex touches about 2.0 cells and a face about 4.0.
    const double fext = argc > 3 ? atof(argv[3]) : 0.90;
    const double vext = argc > 4 ? atof(argv[4]) : 0.45;

    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos(0.0, span);
    std::uniform_real_distribution<double> fe(0.0, 2.0 * fext), ve(0.0, 2.0 * vext);

    // A heavy tail, as the real scenes have: DECISIONS records the widest swept
    // element at 7.5x the mean on armadillo-rollers and 24x on cloth-funnel.
    // This is the case that decides the row bound, since K1 is a maximum and one
    // oversized face raises it for every vertex.
    const int n_big = argc > 5 ? atoi(argv[5]) : 0;
    const double big = argc > 6 ? atof(argv[6]) : 20.0 * fext;

    std::vector<Box> f(nf), v(nv);
    for (int i = 0; i < nf; ++i)
        for (int d = 0; d < 2; ++d) { f[i].lo[d] = pos(rng); f[i].hi[d] = f[i].lo[d] + fe(rng); }
    for (int k = 0; k < n_big && k < nf; ++k) {
        const int i = (int)((long)k * nf / (n_big ? n_big : 1));
        for (int d = 0; d < 2; ++d) { f[i].lo[d] = pos(rng); f[i].hi[d] = f[i].lo[d] + big; }
    }
    for (int i = 0; i < nv; ++i)
        for (int d = 0; d < 2; ++d) { v[i].lo[d] = pos(rng); v[i].hi[d] = v[i].lo[d] + ve(rng); }

    std::set<std::pair<int, int>> brute;   // (face, vertex)
    for (int i = 0; i < nf; ++i)
        for (int j = 0; j < nv; ++j)
            if (overlap(f[i], v[j])) brute.insert({i, j});

    // ---------------- A: what ships -- vertices binned by extent, faces walk --
    const Grid ga(v, span);
    long long a_entries = 0, a_cells = 0, a_tested = 0;
    std::set<std::pair<int, int>> a_got;
    double a_vcells = 0, a_fcells = 0;
    {
        std::vector<std::vector<int>> bin(ga.ncells());
        for (int j = 0; j < nv; ++j) {
            const int r0 = ga.c1(v[j].lo[1]), r1 = ga.c1(v[j].hi[1]);
            const int k0 = ga.c0(v[j].lo[0]), k1 = ga.c0(v[j].hi[0]);
            a_vcells += (double)(r1 - r0 + 1) * (k1 - k0 + 1);
            for (int r = r0; r <= r1; ++r)
                for (int c = k0; c <= k1; ++c) { bin[ga.of(c, r)].push_back(j); ++a_entries; }
        }
        for (int i = 0; i < nf; ++i) {
            const int r0 = ga.c1(f[i].lo[1]), r1 = ga.c1(f[i].hi[1]);
            const int k0 = ga.c0(f[i].lo[0]), k1 = ga.c0(f[i].hi[0]);
            a_fcells += (double)(r1 - r0 + 1) * (k1 - k0 + 1);
            for (int r = r0; r <= r1; ++r)
                for (int c = k0; c <= k1; ++c) {
                    ++a_cells;
                    for (int j : bin[ga.of(c, r)]) {
                        ++a_tested;
                        if (!overlap(f[i], v[j])) continue;
                        // the shipped duplicate rule: the cell holding the
                        // minimum corner of the overlap
                        const double o0 = std::max(f[i].lo[0], v[j].lo[0]);
                        const double o1 = std::max(f[i].lo[1], v[j].lo[1]);
                        if (ga.c0(o0) != c || ga.c1(o1) != r) continue;
                        a_got.insert({i, j});
                    }
                }
        }
        a_vcells /= nv; a_fcells /= nf;
    }

    // ---------- B: proposed -- faces binned at the minimum corner, verts walk --
    const Grid gb(f, span);
    long long b_entries = nf, b_cells = 0, b_tested = 0, b_cols_skipped = 0,
              b_rows_skipped = 0, b_searches = 0;
    std::set<std::pair<int, int>> b_got;
    int K1 = 0, K0 = 0;
    {
        const double NEG = -1e300;
        std::vector<std::vector<int>> bin(gb.ncells());
        std::vector<double> cell_hi0(gb.ncells(), NEG), cell_hi1(gb.ncells(), NEG),
                            row_prefix(gb.ncells(), NEG);
        for (int i = 0; i < nf; ++i) {
            const int c = gb.c0(f[i].lo[0]), r = gb.c1(f[i].lo[1]);
            bin[gb.of(c, r)].push_back(i);
            cell_hi0[gb.of(c, r)] = std::max(cell_hi0[gb.of(c, r)], f[i].hi[0]);
            cell_hi1[gb.of(c, r)] = std::max(cell_hi1[gb.of(c, r)], f[i].hi[1]);
            K0 = std::max(K0, gb.c0(f[i].hi[0]) - c);
            K1 = std::max(K1, gb.c1(f[i].hi[1]) - r);
        }
        for (int r = 0; r < gb.n[1]; ++r) {
            double run = NEG;
            for (int c = 0; c < gb.n[0]; ++c) {
                run = std::max(run, cell_hi0[gb.of(c, r)]);
                row_prefix[gb.of(c, r)] = run;
            }
        }

        for (int j = 0; j < nv; ++j) {
            const int cmax = gb.c0(v[j].hi[0]), rmax = gb.c1(v[j].hi[1]);
            // A face reaching v on axis 1 has row(f.max) >= row(v.min), and it
            // spans at most K1 rows, so it is binned no lower than
            // row(v.min) - K1. Anchoring on row(v.max) instead would cut off
            // valid rows whenever the vertex itself straddles a boundary.
            const int rlo = std::max(0, gb.c1(v[j].lo[1]) - K1);
            b_rows_skipped += rlo;
            for (int r = rlo; r <= rmax; ++r) {
                // One binary search per row scanned, whether or not it yields a
                // cell. This is the cost a large K1 adds, so it is counted.
                ++b_searches;
                const double* pre = row_prefix.data() + (size_t)r * gb.n[0];
                const int start =
                    (int)(std::lower_bound(pre, pre + cmax + 1, v[j].lo[0]) - pre);
                b_cols_skipped += start;
                for (int c = start; c <= cmax; ++c) {
                    const size_t cell = gb.of(c, r);
                    if (cell_hi1[cell] < v[j].lo[1]) { ++b_cols_skipped; continue; }
                    ++b_cells;
                    for (int i : bin[cell]) {
                        ++b_tested;
                        if (overlap(f[i], v[j])) b_got.insert({i, j});   // no dedup rule
                    }
                }
            }
        }
    }

    size_t a_missed = 0, b_missed = 0;
    for (const auto& p : brute) { if (!a_got.count(p)) ++a_missed; if (!b_got.count(p)) ++b_missed; }

    printf("faces=%d verts=%d  extents f=%.2f v=%.2f  pairs=%zu\n",
           nf, nv, fext, vext, brute.size());
    printf("  grid A (vertex-sized) %dx%d   vertex touches %.2f cells, face touches %.2f\n",
           ga.n[0], ga.n[1], a_vcells, a_fcells);
    printf("  grid B (face-sized)   %dx%d   K0=%d K1=%d rows per vertex <= %d\n",
           gb.n[0], gb.n[1], K0, K1, K1 + 1);
    printf("  A shipped : MISSED %-5zu entries %-9lld cells %-9lld candidates %lld\n",
           a_missed, a_entries, a_cells, a_tested);
    printf("  B proposed: MISSED %-5zu entries %-9lld cells %-9lld candidates %lld"
           "  (%lld searches, skipped %lld cols, %lld rows)\n",
           b_missed, b_entries, b_cells, b_tested, b_searches, b_cols_skipped, b_rows_skipped);
    if (a_tested && a_cells && a_entries)
        printf("  B/A: entries %.2fx  cells %.2fx  candidates %.2fx\n",
               (double)b_entries / (double)a_entries,
               (double)b_cells / (double)a_cells,
               (double)b_tested / (double)a_tested);
    return (a_missed || b_missed) ? 1 : 0;
}
