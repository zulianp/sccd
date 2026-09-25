// Does min-corner binning with a forward-only walk find every self pair? No.
//
// The proposal this answers, for the edge-edge broad phase where one list is
// queried against itself: bin each box into the one cell holding its *minimum*
// corner rather than into every cell its extent touches, and have a box read
// only its own footprint, from its min cell up to its max cell. Each box is
// then in one cell, so a partner is met at most once and the shipped duplicate
// rule -- attribute the pair to the cell holding the minimum corner of the
// overlap -- is not needed at all. It also halves the cell array.
//
// It is not conservative, and this program is the evidence. Against brute force
// it misses the pairs whose minimum-corner cells are *incomparable*: one box
// ahead on the first axis and behind on the second, so "forward" from either
// one points away from the other. It reports how many of each it lost, and the
// two counts come out equal -- every comparable pair is found and every
// incomparable one is lost.
//
// Reading forward is sound when the keys are totally ordered, which is what
// sweep and prune exploits on its single axis: the minima sort, so starting the
// window at i+1 yields each unordered pair once. Minimum corners in two
// dimensions are only partially ordered.
//
// Usage: min_corner_self_probe [boxes] [max box extent]
// Exit status is non-zero when a pair was missed, which is always.
// The argument is written up in wip/DECISIONS.md.
#include <cstdio>
#include <cstdlib>
#include <random>
#include <set>
#include <vector>

struct Box { double lo[2], hi[2]; };

// The failure in two boxes, with the arithmetic written out, so it does not
// depend on a random draw or on any index test. At a cell width of 1:
//
//   A = [0.5, 1.5] x [1.5, 2.5]   min corner (0.5, 1.5) -> cell (0, 1)
//   B = [1.2, 2.2] x [0.8, 1.8]   min corner (1.2, 0.8) -> cell (1, 0)
//
// They overlap: x on [1.2, 1.5], y on [1.5, 1.8]. But B sits one cell to the
// LEFT of A and one cell ABOVE it, so neither box's footprint contains the
// other's bin cell and neither ever reads the other. No index test is reached,
// because the inner loop never yields the other index in the first place.
static int two_box_counterexample() {
    const double w = 1.0;
    const double alo[2] = {0.5, 1.5}, ahi[2] = {1.5, 2.5};
    const double blo[2] = {1.2, 0.8}, bhi[2] = {2.2, 1.8};
    auto c = [&](double v) { return (int)(v / w); };

    const bool hit = !(alo[0] > bhi[0] || blo[0] > ahi[0] ||
                       alo[1] > bhi[1] || blo[1] > ahi[1]);
    // Does A's footprint contain B's bin cell, or B's contain A's?
    const bool a_reads_b = c(blo[0]) >= c(alo[0]) && c(blo[0]) <= c(ahi[0]) &&
                           c(blo[1]) >= c(alo[1]) && c(blo[1]) <= c(ahi[1]);
    const bool b_reads_a = c(alo[0]) >= c(blo[0]) && c(alo[0]) <= c(bhi[0]) &&
                           c(alo[1]) >= c(blo[1]) && c(alo[1]) <= c(bhi[1]);

    printf("two-box counterexample: overlap=%s  A reads B=%s  B reads A=%s\n",
           hit ? "yes" : "no", a_reads_b ? "yes" : "no", b_reads_a ? "yes" : "no");
    printf("  A bin cell (%d,%d), footprint cols %d-%d rows %d-%d\n",
           c(alo[0]), c(alo[1]), c(alo[0]), c(ahi[0]), c(alo[1]), c(ahi[1]));
    printf("  B bin cell (%d,%d), footprint cols %d-%d rows %d-%d\n",
           c(blo[0]), c(blo[1]), c(blo[0]), c(bhi[0]), c(blo[1]), c(bhi[1]));
    return (hit && !a_reads_b && !b_reads_a) ? 1 : 0;
}

int main(int argc, char** argv) {
    const int lost_deterministically = two_box_counterexample();
    printf("\n");
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

    // Bin each box at its minimum corner: one cell per box.
    std::vector<std::vector<int>> bin((size_t)nc[0] * nc[1]);
    for (int i = 0; i < n; ++i)
        bin[(size_t)cell(b[i].lo[1], 1) * nc[0] + cell(b[i].lo[0], 0)].push_back(i);

    // Walk the footprint: own min cell up to own max cell, both axes.
    std::set<std::pair<int, int>> got;
    long long emitted = 0;
    for (int i = 0; i < n; ++i) {
        for (int c1 = cell(b[i].lo[1], 1); c1 <= cell(b[i].hi[1], 1); ++c1)
            for (int c0 = cell(b[i].lo[0], 0); c0 <= cell(b[i].hi[0], 0); ++c0)
                for (int j : bin[(size_t)c1 * nc[0] + c0]) {
                    if (j == i) continue;
                    // Same cell: both ends fire, so break the tie on the index.
                    const bool same = cell(b[j].lo[0], 0) == cell(b[i].lo[0], 0) &&
                                      cell(b[j].lo[1], 1) == cell(b[i].lo[1], 1);
                    if (same && j < i) continue;
                    if (!overlap(i, j)) continue;
                    ++emitted;
                    got.insert({i < j ? i : j, i < j ? j : i});
                }
    }

    std::vector<std::pair<int, int>> missed;
    long long comparable_missed = 0, incomparable_total = 0, incomparable_missed = 0;
    for (const auto& p : brute) {
        const int i = p.first, j = p.second;
        const int di0 = cell(b[i].lo[0], 0) - cell(b[j].lo[0], 0);
        const int di1 = cell(b[i].lo[1], 1) - cell(b[j].lo[1], 1);
        // Incomparable: one box's min cell is ahead on one axis and behind on
        // the other, so neither min corner dominates the other.
        const bool incomparable = (di0 > 0 && di1 < 0) || (di0 < 0 && di1 > 0);
        if (incomparable) ++incomparable_total;
        if (!got.count(p)) {
            missed.push_back(p);
            if (incomparable) ++incomparable_missed; else ++comparable_missed;
        }
    }
    printf("  of the missed: %lld have incomparable min cells, %lld comparable\n",
           incomparable_missed, comparable_missed);
    printf("  incomparable pairs in total: %lld\n", incomparable_total);

    printf("boxes=%d  grid=%dx%d  brute=%zu  found=%zu  emitted=%lld  MISSED=%zu\n",
           n, nc[0], nc[1], brute.size(), got.size(), emitted, missed.size());
    for (size_t k = 0; k < missed.size() && k < 3; ++k) {
        const int i = missed[k].first, j = missed[k].second;
        printf("  missed (%d,%d): A=[%.2f,%.2f]x[%.2f,%.2f] cell(%d,%d)  "
               "B=[%.2f,%.2f]x[%.2f,%.2f] cell(%d,%d)\n",
               i, j, b[i].lo[0], b[i].hi[0], b[i].lo[1], b[i].hi[1],
               cell(b[i].lo[0], 0), cell(b[i].lo[1], 1),
               b[j].lo[0], b[j].hi[0], b[j].lo[1], b[j].hi[1],
               cell(b[j].lo[0], 0), cell(b[j].lo[1], 1));
    }
    return (missed.empty() && !lost_deterministically) ? 0 : 1;
}
