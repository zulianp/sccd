// Can the edge-edge self query bin at the minimum corner and read only forward?
// Yes, and on these inputs it is cheaper than what ships.
//
// The proposal: bin each box into the one cell holding its *minimum* corner
// instead of into every cell its extent touches. One entry per box, a partner
// met at most once, and the shipped duplicate rule -- attribute the pair to the
// cell holding the minimum corner of the overlap -- is not needed.
//
// Everything turns on what "forward" means, and two readings are not the same.
//
//   The FOOTPRINT rectangle -- this box's columns crossed with its rows -- is
//   NOT complete. It loses exactly the pairs whose minimum-corner cells are
//   incomparable, one box ahead on the first axis and behind on the second,
//   because then neither footprint holds the other's bin cell and neither ever
//   reads the other. Componentwise order on cells is partial; a forward walk
//   needs a total one. two_box_counterexample() below is that failure in two
//   boxes, with no randomness and no index test involved.
//
//   ROW-MAJOR LINEAR index is total, and over it the walk is complete. For
//   overlapping i and j, L(cell(j.min)) never passes L(cell(i.max)), so the one
//   with the smaller L holds the other in its band. Equal L means one cell,
//   where the index breaks the tie -- the only index test still needed.
//
// Read naively the band is whole rows and costs a grid width per row spanned.
// The point of this probe is that it need not be read naively. Each cell carries
// the largest upper bound of what it holds; a prefix maximum along each row then
// says which columns hold nothing reaching back to this box, and a binary search
// skips them at once. That is the sweep's cummax, per row. What is left is the
// columns that can hold a partner, and the measured cost falls below the shipped
// scheme's on every shape tried -- by more, not less, as the box sizes grow a
// heavy tail, because a wide box costs the shipped binning hundreds of cells and
// costs this one a single entry.
//
// Usage: min_corner_self_probe [boxes] [extent] [big boxes] [big extent]
// Exit status is non-zero if the pruned walk ever misses a pair.
// Written up in wip/DECISIONS.md.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <set>
#include <vector>

struct Box { double lo[2], hi[2]; };


// The footprint reading's failure in two boxes, at a cell width of 1:
//
//   A = [0.5, 1.5] x [1.5, 2.5]   min corner -> cell (0, 1), linear 3 on a grid 3 wide
//   B = [1.2, 2.2] x [0.8, 1.8]   min corner -> cell (1, 0), linear 1
//
// They overlap -- x on [1.2, 1.5], y on [1.5, 1.8] -- and B sits one cell LEFT
// of A and one cell ABOVE it, so neither footprint holds the other's bin cell.
// In linear order B comes first and A's cell falls inside B's band, so the walk
// this probe measures does find them.
static void two_box_counterexample() {
    const int n0 = 3;
    const double alo[2] = {0.5, 1.5}, ahi[2] = {1.5, 2.5};
    const double blo[2] = {1.2, 0.8}, bhi[2] = {2.2, 1.8};
    auto c = [](double v) { return (int)v; };
    auto lin = [&](double x, double y) { return c(y) * n0 + c(x); };

    const bool a_rect = c(blo[0]) >= c(alo[0]) && c(blo[0]) <= c(ahi[0]) &&
                        c(blo[1]) >= c(alo[1]) && c(blo[1]) <= c(ahi[1]);
    const bool b_rect = c(alo[0]) >= c(blo[0]) && c(alo[0]) <= c(bhi[0]) &&
                        c(alo[1]) >= c(blo[1]) && c(alo[1]) <= c(bhi[1]);
    const bool a_band = lin(blo[0], blo[1]) >= lin(alo[0], alo[1]) &&
                        lin(blo[0], blo[1]) <= lin(ahi[0], ahi[1]);
    const bool b_band = lin(alo[0], alo[1]) >= lin(blo[0], blo[1]) &&
                        lin(alo[0], alo[1]) <= lin(bhi[0], bhi[1]);

    printf("two overlapping boxes, cell width 1, grid %d wide\n", n0);
    printf("  A cell (%d,%d) linear %d   B cell (%d,%d) linear %d\n",
           c(alo[0]), c(alo[1]), lin(alo[0], alo[1]),
           c(blo[0]), c(blo[1]), lin(blo[0], blo[1]));
    printf("  footprint walk: A reads B=%s B reads A=%s -> %s\n",
           a_rect ? "yes" : "no", b_rect ? "yes" : "no", (a_rect || b_rect) ? "found" : "LOST");
    printf("  linear walk:    A reads B=%s B reads A=%s -> %s\n\n",
           a_band ? "yes" : "no", b_band ? "yes" : "no", (a_band || b_band) ? "found" : "LOST");
}

int main(int argc, char** argv) {
    two_box_counterexample();

    const int n = argc > 1 ? atoi(argv[1]) : 4000;
    const double span = 100.0;
    const double size = argc > 2 ? atof(argv[2]) : 6.0;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos(0.0, span), ext(0.0, size);

    // A heavy tail, as the real scenes have: a few elements move fast and sweep
    // a box many times the mean. On armadillo-rollers the widest swept edge is
    // 7.5x the mean and on cloth-funnel's edges 24x.
    const int n_big = argc > 3 ? atoi(argv[3]) : 0;
    const double big = argc > 4 ? atof(argv[4]) : 20.0 * size;
    std::vector<Box> b(n);
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < 2; ++d) { b[i].lo[d] = pos(rng); b[i].hi[d] = b[i].lo[d] + ext(rng); }
    for (int k = 0; k < n_big && k < n; ++k) {
        const int i = (int)((long)k * n / (n_big ? n_big : 1));
        for (int d = 0; d < 2; ++d) { b[i].lo[d] = pos(rng); b[i].hi[d] = b[i].lo[d] + big; }
    }

    double mean[2] = {0, 0};
    for (int i = 0; i < n; ++i) for (int d = 0; d < 2; ++d) mean[d] += b[i].hi[d] - b[i].lo[d];
    for (int d = 0; d < 2; ++d) mean[d] /= n;
    const int n0 = (int)(span / mean[0]) + 1, n1 = (int)(span / mean[1]) + 1;
    auto cx = [&](double v) { int c = (int)(v / mean[0]); return c < 0 ? 0 : (c >= n0 ? n0 - 1 : c); };
    auto cy = [&](double v) { int c = (int)(v / mean[1]); return c < 0 ? 0 : (c >= n1 ? n1 - 1 : c); };
    auto overlap = [&](int i, int j) {
        for (int d = 0; d < 2; ++d)
            if (b[i].lo[d] > b[j].hi[d] || b[j].lo[d] > b[i].hi[d]) return false;
        return true;
    };

    std::set<std::pair<int, int>> brute;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (overlap(i, j)) brute.insert({i, j});

    // --- bin at the minimum corner: one entry per box ---
    const size_t ncells = (size_t)n0 * n1;
    std::vector<std::vector<int>> bin(ncells);
    for (int i = 0; i < n; ++i) bin[(size_t)cy(b[i].lo[1]) * n0 + cx(b[i].lo[0])].push_back(i);

    // Per cell, the largest upper bound of what it holds; then a prefix maximum
    // along each row on the first axis.
    const double NEG = -1e300;
    std::vector<double> cell_hi0(ncells, NEG), cell_hi1(ncells, NEG), row_prefix(ncells, NEG);
    for (int i = 0; i < n; ++i) {
        const size_t c = (size_t)cy(b[i].lo[1]) * n0 + cx(b[i].lo[0]);
        cell_hi0[c] = std::max(cell_hi0[c], b[i].hi[0]);
        cell_hi1[c] = std::max(cell_hi1[c], b[i].hi[1]);
    }
    for (int r = 0; r < n1; ++r) {
        double run = NEG;
        for (int c = 0; c < n0; ++c) {
            run = std::max(run, cell_hi0[(size_t)r * n0 + c]);
            row_prefix[(size_t)r * n0 + c] = run;
        }
    }

    std::set<std::pair<int, int>> got;
    long long cells_read = 0, tested = 0, cells_skipped = 0;
    for (int i = 0; i < n; ++i) {
        const int c0b = cx(b[i].lo[0]), c0e = cx(b[i].hi[0]);
        const int c1b = cy(b[i].lo[1]), c1e = cy(b[i].hi[1]);
        const long lin_self = (long)c1b * n0 + c0b;

        for (int r = c1b; r <= c1e; ++r) {
            // First column whose row prefix reaches back to this box's minimum.
            // Everything to its left holds only boxes ending before it.
            const double* pre = row_prefix.data() + (size_t)r * n0;
            int start = (int)(std::lower_bound(pre, pre + c0e + 1, b[i].lo[0]) - pre);
            if (r == c1b) start = std::max(start, c0b);
            cells_skipped += start;

            for (int c = start; c <= c0e; ++c) {
                const size_t cell = (size_t)r * n0 + c;
                const long lin = (long)r * n0 + c;
                if (lin < lin_self) continue;
                if (cell_hi1[cell] < b[i].lo[1]) { ++cells_skipped; continue; }
                ++cells_read;
                for (int j : bin[cell]) {
                    if (lin == lin_self && j <= i) continue;
                    ++tested;
                    if (overlap(i, j)) got.insert({i < j ? i : j, i < j ? j : i});
                }
            }
        }
    }

    // The shipped scheme on the same boxes, for cost.
    long long ship_cells = 0, ship_tested = 0, ship_entries = 0;
    {
        std::vector<std::vector<int>> eb(ncells);
        for (int i = 0; i < n; ++i)
            for (int r = cy(b[i].lo[1]); r <= cy(b[i].hi[1]); ++r)
                for (int c = cx(b[i].lo[0]); c <= cx(b[i].hi[0]); ++c)
                    { eb[(size_t)r * n0 + c].push_back(i); ++ship_entries; }
        for (int i = 0; i < n; ++i)
            for (int r = cy(b[i].lo[1]); r <= cy(b[i].hi[1]); ++r)
                for (int c = cx(b[i].lo[0]); c <= cx(b[i].hi[0]); ++c) {
                    ++ship_cells;
                    for (int j : eb[(size_t)r * n0 + c]) if (j > i) ++ship_tested;
                }
    }

    size_t missed = 0;
    for (const auto& p : brute) if (!got.count(p)) ++missed;
    printf("boxes=%d (%d big, extent %.0f vs %.1f mean) grid=%dx%d pairs=%zu\n",
           n, n_big, big, 0.5 * size, n0, n1, brute.size());
    printf("  pruned band: MISSED %-6zu cells read %-9lld candidates %-10lld entries %d\n",
           missed, cells_read, tested, n);
    printf("  shipped:              cells read %-9lld candidates %-10lld entries %lld\n",
           ship_cells, ship_tested, ship_entries);
    if (tested) printf("  candidates: pruned/shipped = %.2fx   cell entries = %.2fx\n",
                       (double)tested / (double)ship_tested, (double)n / (double)ship_entries);
    return missed ? 1 : 0;
}
