// The C ABI, which is the library's only compiled translation unit and the basis
// for python/sccd.py, had no test at all.
//
// That mattered more than an untested corner usually does. The C entry points do
// not share a dispatcher with the C++ API -- each one assembles its own search --
// so nothing that covers narrowphase.hpp covers them. They drifted, and the drift
// was invisible: until recently every C entry point defaulted SCCD_ADAPTIVE_SPLIT
// to 0, so the installed ABI ran the uniform splitter that the assessment
// measured as never winning, while the C++ path ran the adaptive one.
//
// The queries here have times of impact that can be written down, so the exports
// are checked against numbers rather than against each other. Conservativeness is
// the property under test: a reported time of impact must never be later than the
// true one, misses are illegal, and earlier is fine.
//
// The declarations come from sccd.h rather than being re-typed here. They used
// to be re-typed, which meant this file could compile and link against a
// signature the library no longer had -- C linkage does not mangle, so the
// mismatch would have surfaced as wrong answers rather than as a build error.

#include "sccd.h"

#include <cmath>
#include <cstdio>

namespace {

    int failures = 0;

    // A vertex dropping through the triangle (0,0,0), (1,0,0), (0,1,0), which lies
    // in z = 0. The triangle does not move, so the crossing time is where the
    // vertex's z reaches zero, and the barycentric position is just (x, y).
    struct VfCase {
        double x, y, z0, z1;
        double expect;   // true time of impact, when there is one
        bool hit;
        const char* name;
    };

    const VfCase kVfCases[] = {
        {0.25, 0.25,  1.0, -1.0, 0.500, true,  "vf: crosses at t=1/2, inside"},
        {0.20, 0.30,  1.0, -3.0, 0.250, true,  "vf: crosses at t=1/4, inside"},
        {0.90, 0.90,  1.0, -1.0, 0.000, false, "vf: crosses outside the triangle"},
        {0.25, 0.25,  1.0,  0.5, 0.000, false, "vf: never reaches the plane"},
    };

    void report(const char* name, const char* fn, const bool hit, const double toi,
                const VfCase& c) {
        if (hit != c.hit) {
            std::printf("  FAIL %-34s [%s] expected %s, got %s\n",
                        name, fn, c.hit ? "hit" : "miss", hit ? "hit" : "miss");
            ++failures;
            return;
        }
        if (!c.hit) {
            std::printf("  %-34s [%-28s] miss, as expected\n", name, fn);
            return;
        }
        // Late is illegal. Early is safe, but wildly early would mean the search
        // is not converging, so that is flagged too.
        const double err = toi - c.expect;
        if (err > 1e-5) {
            std::printf("  FAIL %-34s [%s] LATE: %.9f vs true %.6f\n", name, fn, toi, c.expect);
            ++failures;
        } else if (err < -0.05) {
            std::printf("  FAIL %-34s [%s] far too early: %.9f vs true %.6f\n",
                        name, fn, toi, c.expect);
            ++failures;
        } else {
            std::printf("  %-34s [%-28s] toi=%.9f  err=%+.2e\n", name, fn, toi, err);
        }
    }

    void run_vf() {
        const double t3d[3][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};

        for (const VfCase& c : kVfCases) {
            const double sv[3] = {c.x, c.y, c.z0};
            const double ev[3] = {c.x, c.y, c.z1};

            {
                double t = 1, u = 0, v = 0;
                const bool hit = sccd_find_root_vf_d(69, 1e-8, sv, t3d[0], t3d[1], t3d[2],
                                                     ev, t3d[0], t3d[1], t3d[2], &t, &u, &v) != 0;
                report(c.name, "sccd_find_root_vf_d", hit, t, c);
            }
            {
                double t = 1, u = 0, v = 0;
                const bool hit = sccd_find_root_bisection_vf_d(69, 1e-8, sv, t3d[0], t3d[1], t3d[2],
                                                               ev, t3d[0], t3d[1], t3d[2],
                                                               &t, &u, &v) != 0;
                report(c.name, "sccd_find_root_bisection_vf_d", hit, t, c);
            }
            {
                const float svf[3] = {(float)sv[0], (float)sv[1], (float)sv[2]};
                const float evf[3] = {(float)ev[0], (float)ev[1], (float)ev[2]};
                float f0[3] = {0, 0, 0}, f1[3] = {1, 0, 0}, f2[3] = {0, 1, 0};
                float t = 1, u = 0, v = 0;
                const bool hit = sccd_find_root_vf_f(69, 1e-6f, svf, f0, f1, f2,
                                                     evf, f0, f1, f2, &t, &u, &v) != 0;
                report(c.name, "sccd_find_root_vf_f", hit, (double)t, c);
            }
            {
                const float svf[3] = {(float)sv[0], (float)sv[1], (float)sv[2]};
                const float evf[3] = {(float)ev[0], (float)ev[1], (float)ev[2]};
                float f0[3] = {0, 0, 0}, f1[3] = {1, 0, 0}, f2[3] = {0, 1, 0};
                float t = 1, u = 0, v = 0;
                const bool hit = sccd_find_root_bisection_vf_f(69, 1e-6f, svf, f0, f1, f2,
                                                               evf, f0, f1, f2, &t, &u, &v) != 0;
                report(c.name, "sccd_find_root_bisection_vf_f", hit, (double)t, c);
            }
        }
    }

    // Two segments, one along x at z = 0 and one along y that starts above and
    // drops through it. They cross at the time the moving segment's z reaches
    // zero, provided they overlap in the other two axes when it does.
    void run_ee() {
        struct EeCase { double z0, z1, yoff, expect; bool hit; const char* name; };
        const EeCase cases[] = {
            { 1.0, -1.0, 0.0, 0.500, true,  "ee: crosses at t=1/2"},
            { 1.0, -3.0, 0.0, 0.250, true,  "ee: crosses at t=1/4"},
            { 1.0, -1.0, 5.0, 0.000, false, "ee: passes far to the side"},
            { 1.0,  0.5, 0.0, 0.000, false, "ee: never reaches the plane"},
        };

        // Static segment: (-1,0,0) -> (1,0,0).
        const double a0[3] = {-1, 0, 0};
        const double a1[3] = { 1, 0, 0};

        for (const EeCase& c : cases) {
            // Moving segment: (0,-1,z) -> (0,1,z), swept in z, offset in y.
            const double b0s[3] = {0, -1 + c.yoff, c.z0};
            const double b1s[3] = {0,  1 + c.yoff, c.z0};
            const double b0e[3] = {0, -1 + c.yoff, c.z1};
            const double b1e[3] = {0,  1 + c.yoff, c.z1};

            double t = 1, u = 0, v = 0;
            const bool hit = sccd_find_root_ee_d(69, 1e-8, a0, a1, b0s, b1s,
                                                 a0, a1, b0e, b1e, &t, &u, &v) != 0;
            const VfCase as_vf{0, 0, 0, 0, c.expect, c.hit, c.name};
            report(c.name, "sccd_find_root_ee_d", hit, t, as_vf);
        }
    }

}  // namespace

int main() {
    std::printf("C ABI, vertex-face:\n");
    run_vf();
    std::printf("C ABI, edge-edge:\n");
    run_ee();
    std::printf("%s\n", failures ? "FAIL" : "OK: every C export is conservative on every case");
    return failures ? 1 : 0;
}
