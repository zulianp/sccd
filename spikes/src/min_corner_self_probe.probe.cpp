// Can the self query bin at the minimum corner and read only forward?
//
// The proposal, for the edge-edge broad phase where one list is queried against
// itself: bin each box into the one cell holding its *minimum* corner rather
// than into every cell its extent touches, and read only cells that come after
// its own. Each box is then in one cell, so a partner is met at most once and
// the shipped duplicate rule -- attribute the pair to the cell holding the
// minimum corner of the overlap -- is not needed. It also shrinks the cell array
// from one entry per covered cell to one per box.
//
// It depends entirely on what "after" means, and this program runs both
// readings against a brute-force reference.
//
//   RECTANGLE. Read the cells of the box's own footprint: columns from its
//   minimum to its maximum, crossed with rows likewise. This is NOT complete.
//   It loses exactly the pairs whose minimum-corner cells are incomparable --
//   one box ahead on the first axis and behind on the second -- because then
//   neither box's footprint contains the other's bin cell and neither ever
//   reads the other. Componentwise ordering on cells is partial, and a forward
//   walk needs a total one.
//
//   BAND. Read the cells whose row-major linear index lies between the box's
//   own cell and the cell of its maximum corner. Linear index IS a total order,
//   and this is complete: for overlapping i and j, L(cell(j.min)) is never past
//   L(cell(i.max)), so whichever of the two has the smaller L holds the other
//   inside its band. Equal L means the same cell, where the index breaks the
//   tie.
//
// The band walk is the one that works, and its cost is why the shipped code
// does not use it: "forward" in a total order on cells means whole rows, so a
// box spanning R rows reads about R times the grid width instead of its own
// footprint. The two cost counters below say by how much on a given shape.
//
// Usage: min_corner_self_probe [boxes] [max box extent]
// Exit status is non-zero if the band walk ever misses a pair.
// Written up in wip/DECISIONS.md.
#include <cstdio>
#include <cstdlib>
#include <random>
#include <set>
#include <vector>

struct Box { double lo[2], hi[2]; };

// The rectangle walk's failure in two boxes, with the arithmetic written out,
// so it depends on no random draw and on no index test. At a cell width of 1:
//
//   A = [0.5, 1.5] x [1.5, 2.5]   minimum corner (0.5, 1.5) -> cell (0, 1)
//   B = [1.2, 2.2] x [0.8, 1.8]   minimum corner (1.2, 0.8) -> cell (1, 0)
//
// They overlap: x on [1.2, 1.5], y on [1.5, 1.8]. B sits one cell to the LEFT
// of A and one cell ABOVE it, so neither footprint contains the other's bin
// cell. The band walk does find them, because B's linear index is the smaller
// and A's cell falls inside B's band.
static void two_box_counterexample(int n0) {
    const double w = 1.0;
    const double alo[2] = {0.5, 1.5}, ahi[2] = {1.5, 2.5};
    const double blo[2] = {1.2, 0.8}, bhi[2] = {2.2, 1.8};
    auto c = [&](double v) { return (int)(v / w); };
    auto lin = [&](double x, double y) { return c(y) * n0 + c(x); };

    const bool hit = !(alo[0] > bhi[0] || blo[0] > ahi[0] ||
                       alo[1] > bhi[1] || blo[1] > ahi[1]);
    const bool a_rect_reads_b = c(blo[0]) >= c(alo[0]) && c(blo[0]) <= c(ahi[0]) &&
                                c(blo[1]) >= c(alo[1]) && c(blo[1]) <= c(ahi[1]);
    const bool b_rect_reads_a = c(alo[0]) >= c(blo[0]) && c(alo[0]) <= c(bhi[0]) &&
                                c(alo[1]) >= c(blo[1]) && c(alo[1]) <= c(bhi[1]);
    const bool a_band_reads_b = lin(blo[0], blo[1]) >= lin(alo[0], alo[1]) &&
                                lin(blo[0], blo[1]) <= lin(ahi[0], ahi[1]);
    const bool b_band_reads_a = lin(alo[0], alo[1]) >= lin(blo[0], blo[1]) &&
                                lin(alo[0], alo[1]) <= lin(bhi[0], bhi[1]);

    printf("two boxes at cell width 1, grid %d wide: overlap=%s\n", n0, hit ? "yes" : "no");
    printf("  A bin cell (%d,%d) linear %d, footprint cols %d-%d rows %d-%d\n",
           c(alo[0]), c(alo[1]), lin(alo[0], alo[1]), c(alo[0]), c(ahi[0]), c(alo[1]), c(ahi[1]));
    printf("  B bin cell (%d,%d) linear %d, footprint cols %d-%d rows %d-%d\n",
           c(blo[0]), c(blo[1]), lin(blo[0], blo[1]), c(blo[0]), c(bhi[0]), c(blo[1]), c(bhi[1]));
    printf("  rectangle walk: A reads B=%s, B reads A=%s  -> pair %s\n",
           a_rect_reads_b ? "yes" : "no", b_rect_reads_a ? "yes" : "no",
           (a_rect_reads_b || b_rect_reads_a) ? "found" : "LOST");
    printf("  band walk:      A reads B=%s, B reads A=%s  -> pair %s\n\n",
           a_band_reads_b ? "yes" : "no", b_band_reads_a ? "yes" : "no",
           (a_band_reads_b || b_band_reads_a) ? "found" : "LOST");
}

int main(int argc, char** argv) {
    two_box_counterexample(3);

    const int n = argc > 1 ? atoi(argv[1]) : 4000;
    const double span = 100.0;
    const double size = argc > 2 ? atof(argv[2]) : 6.0;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos(0.0, span), ext(0.0, size);

    std::vector<Box> b(n);
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < 2; ++d) { b[i].lo[d] = pos(rng); b[i].hi[d] = b[i].lo[d] + ext(rng); }

    // Cell width: the mean extent, as the shipped grid uses.
    double mean[2] = {0, 0};
    for (int i = 0; i < n; ++i) for (int d = 0; d < 2; ++d) mean[d] += b[i].hi[d] - b[i].lo[d];
    for (int d = 0; d < 2; ++d) mean[d] /= n;
    const int nc[2] = {(int)(span / mean[0]) + 1, (int)(span / mean[1]) + 1};
    auto cell = [&](double v, int d) {
        int c = (int)(v / mean[d]);
        return c < 0 ? 0 : (c >= nc[d] ? nc[d] - 1 : c);
    };
    auto overlap = [&](int i, int j) {
        for (int d = 0; d < 2; ++d)
            if (b[i].lo[d] > b[j].hi[d] || b[j].lo[d] > b[i].hi[d]) return false;
        return true;
    };

    std::set<std::pair<int, int>> brute;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (overlap(i, j)) brute.insert({i, j});

    // One cell per box, at its minimum corner. Both walks read this.
    std::vector<std::vector<int>> bin((size_t)nc[0] * nc[1]);
    for (int i = 0; i < n; ++i)
        bin[(size_t)cell(b[i].lo[1], 1) * nc[0] + cell(b[i].lo[0], 0)].push_back(i);

    std::set<std::pair<int, int>> rect, band;
    long long rect_cells = 0, rect_tested = 0, band_cells = 0, band_tested = 0;

    for (int i = 0; i < n; ++i) {
        const int c0b = cell(b[i].lo[0], 0), c0e = cell(b[i].hi[0], 0);
        const int c1b = cell(b[i].lo[1], 1), c1e = cell(b[i].hi[1], 1);

        for (int c1 = c1b; c1 <= c1e; ++c1)
            for (int c0 = c0b; c0 <= c0e; ++c0) {
                ++rect_cells;
                for (int j : bin[(size_t)c1 * nc[0] + c0]) {
                    if (c0 == c0b && c1 == c1b && j <= i) continue;
                    ++rect_tested;
                    if (overlap(i, j)) rect.insert({i < j ? i : j, i < j ? j : i});
                }
            }

        const long lo = (long)c1b * nc[0] + c0b, hi = (long)c1e * nc[0] + c0e;
        for (long c = lo; c <= hi; ++c) {
            ++band_cells;
            for (int j : bin[(size_t)c]) {
                if (c == lo && j <= i) continue;
                ++band_tested;
                if (overlap(i, j)) band.insert({i < j ? i : j, i < j ? j : i});
            }
        }
    }

    // What the shipped scheme costs on the same boxes: bin over the whole
    // extent, read the same footprint, settle duplicates with the corner rule.
    long long ship_cells = 0, ship_tested = 0, ship_entries = 0;
    {
        std::vector<std::vector<int>> ext_bin((size_t)nc[0] * nc[1]);
        for (int i = 0; i < n; ++i)
            for (int c1 = cell(b[i].lo[1], 1); c1 <= cell(b[i].hi[1], 1); ++c1)
                for (int c0 = cell(b[i].lo[0], 0); c0 <= cell(b[i].hi[0], 0); ++c0)
                    { ext_bin[(size_t)c1 * nc[0] + c0].push_back(i); ++ship_entries; }
        for (int i = 0; i < n; ++i)
            for (int c1 = cell(b[i].lo[1], 1); c1 <= cell(b[i].hi[1], 1); ++c1)
                for (int c0 = cell(b[i].lo[0], 0); c0 <= cell(b[i].hi[0], 0); ++c0) {
                    ++ship_cells;
                    for (int j : ext_bin[(size_t)c1 * nc[0] + c0]) if (j > i) ++ship_tested;
                }
    }

    size_t rect_missed = 0, band_missed = 0;
    for (const auto& p : brute) {
        if (!rect.count(p)) ++rect_missed;
        if (!band.count(p)) ++band_missed;
    }

    printf("boxes=%d  grid=%dx%d  overlapping pairs=%zu\n", n, nc[0], nc[1], brute.size());
    printf("  %-22s found %-8zu MISSED %-8zu cells read %-10lld candidates %lld\n",
           "min corner, rectangle", rect.size(), rect_missed, rect_cells, rect_tested);
    printf("  %-22s found %-8zu MISSED %-8zu cells read %-10lld candidates %lld\n",
           "min corner, band", band.size(), band_missed, band_cells, band_tested);
    printf("  %-22s %-23s cells read %-10lld candidates %lld  (%lld cell entries)\n",
           "shipped, full extent", "", ship_cells, ship_tested, ship_entries);
    return band_missed ? 1 : 0;
}
