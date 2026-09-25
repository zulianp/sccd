// Must the walk look left of the querying box's own column on later rows? Yes.
//
// The question comes up because the cells to the left are the ones a drawing of
// the algorithm cannot justify: nothing about the querying box suggests they
// would ever be read. They are read because a box binned in a lower column can
// still be wide enough to reach back, and the box holding the lower cell is the
// one that must report that pair.
//
// Two variants, both binning at the minimum corner and both deduplicating by
// row-major linear cell index -- the pair is emitted by the box in the lower
// cell, ties on the box index. They differ only in the left bound of the column
// scan on rows after the querying box's own:
//
//   WIDE : columns 0 .. c0e          (what ships)
//   NARROW: columns c0b .. c0e       (never left of the box's own column)
//
// Both are checked against brute force.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <set>
#include <vector>

struct Box { double lo[2], hi[2]; };

int main(int argc, char** argv) {
    const int n = argc > 1 ? atoi(argv[1]) : 4000;
    const double span = 100.0;
    const double size = argc > 2 ? atof(argv[2]) : 6.0;
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> pos(0.0, span), ext(0.0, size);

    std::vector<Box> b(n);
    for (int i = 0; i < n; ++i)
        for (int d = 0; d < 2; ++d) { b[i].lo[d] = pos(rng); b[i].hi[d] = b[i].lo[d] + ext(rng); }

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

    const size_t ncells = (size_t)n0 * n1;
    std::vector<std::vector<int>> bin(ncells);
    for (int i = 0; i < n; ++i) bin[(size_t)cy(b[i].lo[1]) * n0 + cx(b[i].lo[0])].push_back(i);

    const double NEG = -1e300;
    std::vector<double> hi0(ncells, NEG), hi1(ncells, NEG), pre(ncells, NEG);
    for (int i = 0; i < n; ++i) {
        const size_t c = (size_t)cy(b[i].lo[1]) * n0 + cx(b[i].lo[0]);
        hi0[c] = std::max(hi0[c], b[i].hi[0]);
        hi1[c] = std::max(hi1[c], b[i].hi[1]);
    }
    for (int r = 0; r < n1; ++r) {
        double run = NEG;
        for (int c = 0; c < n0; ++c) { run = std::max(run, hi0[(size_t)r * n0 + c]); pre[(size_t)r * n0 + c] = run; }
    }

    auto walk = [&](const bool narrow) {
        std::set<std::pair<int, int>> got;
        for (int i = 0; i < n; ++i) {
            const int c0b = cx(b[i].lo[0]), c0e = cx(b[i].hi[0]);
            const int c1b = cy(b[i].lo[1]), c1e = cy(b[i].hi[1]);
            const long own = (long)c1b * n0 + c0b;
            for (int r = c1b; r <= c1e; ++r) {
                const double* p = pre.data() + (size_t)r * n0;
                int c = (int)(std::lower_bound(p, p + c0e + 1, b[i].lo[0]) - p);
                if (r == c1b || narrow) c = std::max(c, c0b);
                for (; c <= c0e; ++c) {
                    const size_t cell = (size_t)r * n0 + c;
                    if (hi1[cell] < b[i].lo[1]) continue;
                    for (int j : bin[cell]) {
                        if ((long)cell == own && j <= i) continue;
                        if (overlap(i, j)) got.insert({i < j ? i : j, i < j ? j : i});
                    }
                }
            }
        }
        return got;
    };

    for (const bool narrow : {false, true}) {
        const auto got = walk(narrow);
        size_t missed = 0;
        std::pair<int, int> first{-1, -1};
        for (const auto& p : brute)
            if (!got.count(p)) { if (missed == 0) first = p; ++missed; }
        printf("%-7s found %-8zu of %-8zu  MISSED %zu", narrow ? "NARROW" : "WIDE",
               got.size(), brute.size(), missed);
        if (missed) {
            const int i = first.first, j = first.second;
            printf("   e.g. (%d,%d): cells %d and %d",
                   i, j, cy(b[i].lo[1]) * n0 + cx(b[i].lo[0]), cy(b[j].lo[1]) * n0 + cx(b[j].lo[0]));
        }
        printf("\n");
    }
    return 0;
}
