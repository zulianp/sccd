// sccd::numerical_error_bound must reproduce TightInclusion's
// get_numerical_error exactly. The predicate that decides whether a box can
// contain a root pads the corner bounds by this value, so if ours is ever
// smaller than TightInclusion's we can reject a box that really does contain a
// root -- a missed collision. Bitwise equality is the cheap way to keep that
// guarantee, since it is the same closed-form filter * min(max|coord|, 1)^3.
#ifndef SCCD_ENABLE_TIGHT_INCLUSION
int main() { return 0; }
#else
#include "snumerical_error.hpp"
#include "tight_inclusion/interval_root_finder.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdint>

static uint64_t st = 0x243F6A8885A308D3ull;
static double urand(double lo, double hi) {
    st ^= st << 13; st ^= st >> 7; st ^= st << 17;
    return lo + (double(st >> 11) * (1.0 / 9007199254740992.0)) * (hi - lo);
}

template <bool VF>
static int check(const char* name, int trials, double lo, double hi) {
    int mism = 0;
    double worst = 0;
    for (int k = 0; k < trials; ++k) {
        double P[8][3];
        for (int i = 0; i < 8; ++i)
            for (int d = 0; d < 3; ++d) P[i][d] = urand(lo, hi);

        double ours[3];
        sccd::numerical_error_bound<VF, double>(P[0],P[1],P[2],P[3],P[4],P[5],P[6],P[7], ours);

        std::vector<ticcd::Vector3> vs;
        for (int i = 0; i < 8; ++i) vs.emplace_back(P[i][0], P[i][1], P[i][2]);
        const ticcd::Array3 ti = ticcd::get_numerical_error(vs, VF, false);

        for (int d = 0; d < 3; ++d) {
            const double diff = std::abs(ours[d] - ti[d]);
            if (diff != 0.0) {
                ++mism;
                worst = std::max(worst, diff / std::max(std::abs(ti[d]), 1e-300));
                break;
            }
        }
    }
    printf("%-26s trials=%d  bitwise-mismatches=%d  worst_rel=%.3e\n", name, trials, mism, worst);
    return mism;
}

int main() {
    int bad = 0;
    bad += check<true >("VF, coords in [-1,1]",   20000, -1.0, 1.0);
    bad += check<true >("VF, coords in [-100,100]",20000, -100.0, 100.0);
    bad += check<true >("VF, coords in [-1e-3,1e-3]",20000, -1e-3, 1e-3);
    bad += check<false>("EE, coords in [-1,1]",   20000, -1.0, 1.0);
    bad += check<false>("EE, coords in [-100,100]",20000, -100.0, 100.0);
    printf("\n%s\n", bad == 0 ? "MATCH: identical to TightInclusion on every trial"
                              : "DIFFERS from TightInclusion");
    return bad != 0;
}
#endif
