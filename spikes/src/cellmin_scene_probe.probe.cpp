// The scene drawn in the paper's edge-edge figure, run through the shipped
// bounds, printing which columns each row of each walk actually scans.
//
// The figure claims a cell index, a walk range and a set of skipped columns for
// each of seven boxes. Those are easy to draw wrong and hard to check by eye --
// a first version marked the columns the running maximum skips on one row of
// one panel and left the same columns unmarked on the row below it. This prints
// them, so the drawing can be checked against the code rather than against
// somebody's reading of the code.
//
// The grid is set by hand to the one the figure draws rather than derived from
// the boxes, so the two agree by construction.
//
// Usage: cellmin_scene_probe
#include "sccd_broadphase_cell2d.hpp"

#include <algorithm>
#include <cstdio>
#include <vector>

using T = double;
using I = int;

int main() {
    const double S = 0.62;
    const char* name = "abcdefg";
    const double lo0[7] = {0.35, 1.05, 0.80, 1.95, 2.40, 0.30, 0.75};
    const double lo1[7] = {0.30, 0.20, 0.85, 0.45, 1.10, 1.75, 1.95};
    const double hi0[7] = {1.10, 1.65, 1.50, 2.60, 3.10, 0.90, 1.45};
    const double hi1[7] = {1.00, 0.75, 1.50, 1.55, 1.70, 2.30, 2.40};
    const int n = 7;

    std::vector<T> rows[6];
    for (int d = 0; d < 6; ++d) rows[d].resize(n);
    for (int i = 0; i < n; ++i) {
        rows[0][i] = lo0[i]; rows[1][i] = lo1[i]; rows[2][i] = 0.0;
        rows[3][i] = hi0[i]; rows[4][i] = hi1[i]; rows[5][i] = 0.0;
    }
    T* aabb[6];
    for (int d = 0; d < 6; ++d) aabb[d] = rows[d].data();

    // The figure's grid, set by hand rather than derived, so the drawing and
    // this agree by construction.
    sccd::Cell2DGrid<T> g;
    g.axis0 = 0; g.axis1 = 1;
    g.min0 = 0; g.min1 = 0;
    g.n0 = 6; g.n1 = 4;
    g.inv0 = 1.0 / S; g.inv1 = 1.0 / S;

    sccd::Cell2DPartition part;
    sccd::cell2dmin_partition<T>(n, aabb, g, part);
    std::vector<ptrdiff_t> cellptr(g.ncells() + 1);
    sccd::cell2dmin_count<T>(n, aabb, g, part, cellptr.data());
    std::vector<I> cellidx(cellptr[g.ncells()]);
    std::vector<ptrdiff_t> cursor(g.ncells());
    sccd::cell2dmin_fill<T, I>(n, aabb, g, part, cellptr.data(), cellidx.data(), cursor.data());
    std::vector<T> row_prefix(g.ncells()), cell_hi1(g.ncells());
    sccd::cell2dmin_bounds<T>(n, aabb, g, part, row_prefix.data(), cell_hi1.data());

    printf("binning:");
    for (int i = 0; i < n; ++i)
        printf(" %c->%ld", name[i], (long)g.cell_of(g.clamp0(lo0[i]), g.clamp1(lo1[i])));
    printf("\n\n");

    for (int i = 0; i < n; ++i) {
        const int c0b = g.clamp0(lo0[i]), c0e = g.clamp0(hi0[i]);
        const int c1b = g.clamp1(lo1[i]), c1e = g.clamp1(hi1[i]);
        printf("%c: cells %ld..%ld, rows %d..%d, cols %d..%d\n", name[i],
               (long)g.cell_of(c0b, c1b), (long)g.cell_of(c0e, c1e), c1b, c1e, c0b, c0e);
        for (int c1 = c1b; c1 <= c1e; ++c1) {
            const T* pre = row_prefix.data() + (ptrdiff_t)c1 * g.n0;
            int c0 = (int)(std::lower_bound(pre, pre + c0e + 1, lo0[i]) - pre);
            if (c1 == c1b) c0 = std::max(c0, c0b);
            printf("   row %d: considered cols %d..%d, start %d -> ",
                   c1, (c1 == c1b ? c0b : 0), c0e, c0);
            if (c0 > c0e) { printf("row skipped\n"); continue; }
            printf("scans");
            for (int c = c0; c <= c0e; ++c) {
                const ptrdiff_t cell = g.cell_of(c, c1);
                printf(" %ld%s", (long)cell, cell_hi1[cell] < lo1[i] ? "(dropped)" : "");
            }
            printf("\n");
        }
    }
    return 0;
}
