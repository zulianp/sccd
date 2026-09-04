// The per-axis domain tolerances must match TightInclusion's exactly.
//
// They decide when a box is tight enough to accept, so a tolerance that is too
// large accepts a box before its t lower bound has converged and the reported
// time of impact lands after the true one -- a conservativeness violation. This
// caught a real defect: the generated edge-edge form returned the same value for
// tol0 and tol1, and was off by up to 2.5x and 6.9x on the u and v axes.
#ifndef SCCD_ENABLE_TIGHT_INCLUSION
int main() { return 0; }
#else
#include "sccd_tolerance.hpp"
#include "tight_inclusion/ccd.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>

static uint64_t rng = 0x9E3779B97F4A7C15ull;
static double urand(const double lo, const double hi) {
    rng ^= rng << 13;
    rng ^= rng >> 7;
    rng ^= rng << 17;
    return lo + (double(rng >> 11) * (1.0 / 9007199254740992.0)) * (hi - lo);
}

// Tolerances are a division, so exact bit equality is not required; anything at
// rounding level means the same formula.
static constexpr double MAX_REL_DIFF = 1e-14;

static int check(const char* name, const int trials, const double lo, const double hi, const double tol_in) {
    double worst[3] = {0, 0, 0};
    int failures = 0;

    for (int k = 0; k < trials; ++k) {
        double P[8][3];
        for (int i = 0; i < 8; ++i) {
            for (int d = 0; d < 3; ++d) {
                P[i][d] = urand(lo, hi);
            }
        }
        ticcd::Vector3 V[8];
        for (int i = 0; i < 8; ++i) {
            V[i] = ticcd::Vector3(P[i][0], P[i][1], P[i][2]);
        }

        double ours_vf[3], ours_ee[3];
        compute_face_vertex_tolerance<double>(tol_in, P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], ours_vf);
        compute_edge_edge_tolerance<double>(tol_in, P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], ours_ee);

        const ticcd::Array3 ti_vf =
            ticcd::compute_vertex_face_tolerances(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], tol_in);
        const ticcd::Array3 ti_ee =
            ticcd::compute_edge_edge_tolerances(V[0], V[1], V[2], V[3], V[4], V[5], V[6], V[7], tol_in);

        for (int d = 0; d < 3; ++d) {
            for (int which = 0; which < 2; ++which) {
                const double ours = which ? ours_ee[d] : ours_vf[d];
                const double ti = which ? ti_ee[d] : ti_vf[d];
                const double rel = std::abs(ours - ti) / std::max(std::abs(ti), 1e-300);
                if (rel > worst[d]) {
                    worst[d] = rel;
                }
                if (!(rel <= MAX_REL_DIFF)) {
                    ++failures;
                }
            }
        }
    }

    std::printf("%-30s worst rel diff per axis: %.3e %.3e %.3e  failures=%d\n",
                name, worst[0], worst[1], worst[2], failures);
    return failures;
}

int main() {
    int failures = 0;
    failures += check("coords in [-1, 1]", 20000, -1.0, 1.0, 3e-8);
    failures += check("coords in [-100, 100]", 20000, -100.0, 100.0, 3e-8);
    failures += check("coords in [-1e-3, 1e-3]", 20000, -1e-3, 1e-3, 3e-8);
    failures += check("looser codomain tolerance", 20000, -1.0, 1.0, 1e-6);
    std::printf("\n%s\n", failures == 0 ? "MATCH: tolerances agree with TightInclusion"
                                        : "DIFFERS from TightInclusion");
    return failures != 0;
}
#endif
