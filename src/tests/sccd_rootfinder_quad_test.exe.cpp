// Accuracy gate for the vertex-quad root finder.
//
// The quad path had no such gate. `ti_oracle` reads triangle query sets only, so
// the one code path in this repository where an unsound rejection was actually
// found -- `codomain_acceptance_vq` padding with machine epsilon instead of the
// certified bound -- is the one path nothing was checking. That is also the path
// the optimised kernel changes, so it needs checking before the optimisation can
// be trusted.
//
// Three checks, in increasing strength:
//
//  1. The factored evaluation agrees with the reference `diff_vq` to well inside
//     the certified numerical error bound. This is the soundness link for the
//     regrouping: the rejection test pads by that bound, so a difference safely
//     inside it cannot turn an accept into a reject.
//
//  2. Constructed queries with an analytic time of impact. A vertex crossing a
//     stationary planar quad has a root that can be written down, so the kernel
//     can be checked against a number rather than against itself.
//
//  3. Random queries with a root planted at a known parameter point. A root
//     provably exists there, so the true time of impact is at or before it, and
//     the kernel must report a hit no later. Earlier is allowed -- that is the
//     safe direction and the whole design leans on it, which is also why the
//     obvious reference does not work; see the note above planted_query.

#include "sccd_narrowphase_quad.hpp"
#include "sccd_rootfinder_quad.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

using S = double;

namespace {

    int failures = 0;

    // How far before the true time of impact a report is allowed to be.
    //
    // These are accuracy budgets, not safety ones: reporting early is the safe
    // direction and the design leans on it. They exist because "earlier is safe"
    // is not a licence to be arbitrarily early -- a time of impact half a step
    // before the contact costs a solver its step size just as surely as a missed
    // collision costs it correctness. Nothing asserted earliness on this path
    // before, so an accuracy regression of any size passed silently.
    //
    // The values are the tolerance caps the search is supposed to honour:
    // SCCD_MAX_TIME_TOL is 1e-3, and the searches below resolve t to about that.
    //
    // The scale check instead asserts a ratio, because the property it tests is
    // that the answer does not change with the size of the scene at all: a
    // similar scene has the same time of impact, so the earliness at every scale
    // should match the earliness at scale 1 rather than merely stay under a cap.
    constexpr double kScaleDriftFactor = 4.0;
    constexpr double kSlowEarlyBudget  = 1e-3;
    constexpr double kEntryEarlyBudget = 1e-3;

    void fail(const char* what, const double a, const double b) {
        std::printf("  FAIL %s: %.17g vs %.17g\n", what, a, b);
        ++failures;
    }

    struct Query {
        S sv[3], ev[3];
        S s[4][3], e[4][3];
    };

    Query random_query(std::mt19937& rng, const double spread) {
        std::uniform_real_distribution<double> d(-spread, spread);
        Query q;
        for (int c = 0; c < 3; ++c) {
            q.sv[c] = (S)d(rng);
            q.ev[c] = (S)d(rng);
        }
        for (int k = 0; k < 4; ++k) {
            for (int c = 0; c < 3; ++c) {
                q.s[k][c] = (S)d(rng);
                q.e[k][c] = (S)d(rng);
            }
        }
        return q;
    }

    // ---------------------------------------------------------------------
    // 1. The factored evaluation against the reference.
    // ---------------------------------------------------------------------
    int check_evaluation_agrees(const int n_queries, const int n_samples) {
        std::mt19937 rng(20260903);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        double worst_ratio = 0.0;

        for (int i = 0; i < n_queries; ++i) {
            const Query q = random_query(rng, 3.0);

            S err[3];
            sccd::vq_numerical_error<S>(
                q.sv, q.s[0], q.s[1], q.s[2], q.s[3], q.ev, q.e[0], q.e[1], q.e[2], q.e[3], err);

            for (int j = 0; j < n_samples; ++j) {
                const S t = (S)unit(rng), u = (S)unit(rng), v = (S)unit(rng);

                S ref[3];
                sccd::diff_vq<S>(q.sv, q.s[0], q.s[1], q.s[2], q.s[3],
                                 q.ev, q.e[0], q.e[1], q.e[2], q.e[3], t, u, v, ref);

                sccd::VQFrame<S> frame;
                sccd::vq_frame_at<S>(q.sv, q.s[0], q.s[1], q.s[2], q.s[3],
                                     q.ev, q.e[0], q.e[1], q.e[2], q.e[3], t, frame);
                S fact[3];
                sccd::vq_eval_frame<S>(frame, u, v, fact);

                for (int c = 0; c < 3; ++c) {
                    const double diff = std::fabs((double)ref[c] - (double)fact[c]);
                    const double ratio = diff / (double)err[c];
                    if (ratio > worst_ratio) worst_ratio = ratio;
                    if (!(diff <= (double)err[c])) {
                        fail("factored evaluation outside the certified bound", diff, (double)err[c]);
                        return 1;
                    }
                }
            }
        }
        std::printf("  evaluation agreement: worst difference is %.3g of the certified bound\n",
                    worst_ratio);
        return 0;
    }

    // ---------------------------------------------------------------------
    // Run the shipped kernel on one query.
    // ---------------------------------------------------------------------
    bool run_kernel(const Query& q, const int max_depth, const S tol, S& toi_out) {
        S tols[3], widths[3], err[3];
        sccd::compute_vertex_quad_tolerance<S>(
            tol, q.sv, q.s[0], q.s[1], q.s[2], q.s[3], q.ev, q.e[0], q.e[1], q.e[2], q.e[3], tols);
        sccd::compute_vertex_quad_codomain_widths<S>(
            q.sv, q.s[0], q.s[1], q.s[2], q.s[3], q.ev, q.e[0], q.e[1], q.e[2], q.e[3], widths);
        sccd::normalize_vertex_quad_codomain_widths<S>(widths);
        sccd::vq_numerical_error<S>(
            q.sv, q.s[0], q.s[1], q.s[2], q.s[3], q.ev, q.e[0], q.e[1], q.e[2], q.e[3], err);

        S t = S(1), u = S(0), v = S(0);
        std::vector<sccd::Box<S>> stack;
        stack.push_back(sccd::unit_domain_box<S>());

        bool hit = false;
        while (!stack.empty()) {
            sccd::Box<S> box = stack.back();
            stack.pop_back();
            if (box.tuv[0].lower >= t) continue;
            box.tuv[0].upper = sccd::min<S>(box.tuv[0].upper, t);
            if (sccd::find_root_grid_adaptive_split_vq<S>(
                    max_depth, tol, tols, err, widths,
                    q.sv, q.s[0], q.s[1], q.s[2], q.s[3],
                    q.ev, q.e[0], q.e[1], q.e[2], q.e[3],
                    box, t, u, v, stack)) {
                hit = true;
            }
        }
        toi_out = t;
        return hit;
    }

    // ---------------------------------------------------------------------
    // A query with a root planted at a known (t, u, v).
    //
    // The obvious reference -- run the same corner-hull test on a fine uniform
    // grid and compare -- does not work, and it is worth saying why. That test is
    // a conservative OVER-approximation: it accepts cells containing no root, and
    // it reports a cell's t lower bound, which is at or before any root inside.
    // So it gives a LOWER bound on the true time of impact, and the kernel is
    // also required to be at or below the true value. Two lower bounds cannot be
    // ordered against each other, and comparing them reports "late" and "missed"
    // for queries where nothing is wrong.
    //
    // Planting the root instead gives an upper bound to check against. Pick the
    // quad freely, pick a parameter point, evaluate the quad there, and choose
    // the vertex trajectory to pass exactly through that point at that time. A
    // root then provably exists at t*, so the true time of impact is at most t*,
    // and the kernel must report a hit at or before it.
    // ---------------------------------------------------------------------
    struct Planted { Query q; S t_star; bool usable; };

    Planted planted_query(std::mt19937& rng) {
        std::uniform_real_distribution<double> coord(-2.0, 2.0);
        std::uniform_real_distribution<double> param(0.05, 0.95);
        std::uniform_real_distribution<double> dir(-1.0, 1.0);

        Planted p;
        for (int k = 0; k < 4; ++k) {
            for (int d = 0; d < 3; ++d) {
                p.q.s[k][d] = (S)coord(rng);
                // A moving quad, but not so fast that the swept cell is enormous.
                p.q.e[k][d] = p.q.s[k][d] + (S)(0.5 * dir(rng));
            }
        }
        const S t_star = (S)param(rng);
        const S u_star = (S)param(rng);
        const S v_star = (S)param(rng);

        sccd::VQFrame<S> frame;
        sccd::vq_frame_at<S>(p.q.sv, p.q.s[0], p.q.s[1], p.q.s[2], p.q.s[3],
                             p.q.ev, p.q.e[0], p.q.e[1], p.q.e[2], p.q.e[3], t_star, frame);
        // vq_frame_at only reads the quad for P; V needs the vertex, which is what
        // we are about to choose, so recompute the quad point directly.
        const S omu = S(1) - u_star, omv = S(1) - v_star;
        const S w[4] = {omu * omv, u_star * omv, omu * v_star, u_star * v_star};
        S Q[3];
        for (int d = 0; d < 3; ++d) {
            Q[d] = w[0] * frame.P[0][d] + w[1] * frame.P[1][d] + w[2] * frame.P[2][d] + w[3] * frame.P[3][d];
        }

        S dvec[3];
        for (int d = 0; d < 3; ++d) dvec[d] = (S)dir(rng);

        // V(t) = sv + t * dvec, so that V(t_star) == Q.
        for (int d = 0; d < 3; ++d) {
            p.q.sv[d] = Q[d] - t_star * dvec[d];
            p.q.ev[d] = p.q.sv[d] + dvec[d];
        }
        p.t_star = t_star;

        // Confirm the plant survived rounding before relying on it.
        S F[3];
        sccd::diff_vq<S>(p.q.sv, p.q.s[0], p.q.s[1], p.q.s[2], p.q.s[3],
                         p.q.ev, p.q.e[0], p.q.e[1], p.q.e[2], p.q.e[3],
                         t_star, u_star, v_star, F);
        p.usable = (std::fabs((double)F[0]) < 1e-12 &&
                    std::fabs((double)F[1]) < 1e-12 &&
                    std::fabs((double)F[2]) < 1e-12);
        return p;
    }

    // ---------------------------------------------------------------------
    // 2. Analytic cases.
    // ---------------------------------------------------------------------
    int check_analytic() {
        // Stationary unit square in z = 0, parameterised so that (u, v) are the
        // square's own coordinates: s1 at the origin, s2 along u, s3 along v.
        const S sq[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};

        struct Case { double x, y, z0, z1, expect; bool hit; const char* name; };
        const Case cases[] = {
            {0.30, 0.40,  1.0, -1.0, 0.50, true,  "crosses at t=1/2, inside"},
            {0.50, 0.50,  1.0, -3.0, 0.25, true,  "crosses at t=1/4, centre"},
            {0.10, 0.90,  0.5, -0.5, 0.50, true,  "crosses at t=1/2, near a corner"},
            {1.70, 0.50,  1.0, -1.0, 0.00, false, "crosses the plane outside the quad"},
            {0.50, 0.50,  1.0,  0.5, 0.00, false, "never reaches the plane"},
        };

        int bad = 0;
        for (const Case& c : cases) {
            Query q;
            for (int k = 0; k < 4; ++k) {
                for (int d = 0; d < 3; ++d) { q.s[k][d] = sq[k][d]; q.e[k][d] = sq[k][d]; }
            }
            q.sv[0] = (S)c.x; q.sv[1] = (S)c.y; q.sv[2] = (S)c.z0;
            q.ev[0] = (S)c.x; q.ev[1] = (S)c.y; q.ev[2] = (S)c.z1;

            S toi = S(1);
            const bool hit = run_kernel(q, 69, S(1e-8), toi);

            if (hit != c.hit) {
                std::printf("  FAIL %s: expected %s, got %s\n",
                            c.name, c.hit ? "hit" : "miss", hit ? "hit" : "miss");
                ++failures; ++bad;
                continue;
            }
            if (!c.hit) { std::printf("  %-42s miss, as expected\n", c.name); continue; }

            // Conservative: never later than the true root, and not absurdly early.
            const double err = (double)toi - c.expect;
            if (err > 1e-6) { fail(c.name, (double)toi, c.expect); ++bad; }
            else if (err < -0.05) { fail(std::string(std::string(c.name) + " (far too early)").c_str(),
                                         (double)toi, c.expect); ++bad; }
            else {
                std::printf("  %-42s toi=%.9f  true=%.3f  err=%+.2e\n", c.name, (double)toi, c.expect, err);
            }
        }
        return bad;
    }

    // ---------------------------------------------------------------------
    // 3. Random queries with a planted root.
    // ---------------------------------------------------------------------
    int check_planted_roots(const int n_queries) {
        std::mt19937 rng(7);
        int usable = 0, missed = 0, late = 0;
        double worst_late = 0.0, worst_early = 0.0;

        for (int i = 0; i < n_queries; ++i) {
            const Planted p = planted_query(rng);
            if (!p.usable) continue;
            ++usable;

            S toi = S(1);
            const bool hit = run_kernel(p.q, 69, S(1e-8), toi);

            if (!hit) {
                ++missed;
                std::printf("  MISSED query %d: a root exists at t=%.9f\n", i, (double)p.t_star);
                ++failures;
                continue;
            }
            const double delta = (double)toi - (double)p.t_star;
            if (delta > 1e-9) {
                ++late;
                worst_late = std::max(worst_late, delta);
                std::printf("  LATE query %d: kernel %.9f, planted root %.9f\n",
                            i, (double)toi, (double)p.t_star);
                ++failures;
            } else {
                worst_early = std::min(worst_early, delta);
            }
        }
        std::printf("  planted roots: %d usable queries, %d missed, %d late; "
                    "earliest report was %+.2e before the planted root\n",
                    usable, missed, late, worst_early);
        return missed + late;
    }


    // ---------------------------------------------------------------------
    // 4. Scale invariance.
    //
    // Multiplying every coordinate of a scene by the same factor produces a
    // similar scene: the vertex meets the quad at the same parameter, so the
    // time of impact is unchanged. Scaling the distance tolerance by the same
    // factor keeps the two acceptance conditions that measure a distance --
    // `smaller_than_scalar_tol` and the per-axis tolerances, which divide by
    // Lipschitz constants that scale the same way -- invariant too.
    //
    // So everything in the search is scale-free except the certified numerical
    // error bound, and any drift in the reported time of impact across the
    // scales below is that bound and nothing else.
    // ---------------------------------------------------------------------
    int check_scale_invariance() {
        const S sq[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
        const double scales[] = {1.0, 1e1, 1e2, 1e3, 1e4};
        const double expect = 0.5;

        int bad = 0;
        double reference_error = 0.0;  // the scale-1 result, which the rest must match
        for (const double scale : scales) {
            Query q;
            for (int k = 0; k < 4; ++k) {
                for (int d = 0; d < 3; ++d) {
                    q.s[k][d] = (S)(sq[k][d] * scale);
                    q.e[k][d] = q.s[k][d];
                }
            }
            q.sv[0] = (S)(0.30 * scale); q.sv[1] = (S)(0.40 * scale); q.sv[2] = (S)( 1.0 * scale);
            q.ev[0] = q.sv[0];           q.ev[1] = q.sv[1];           q.ev[2] = (S)(-1.0 * scale);

            S toi = S(1);
            const bool hit = run_kernel(q, 69, (S)(1e-8 * scale), toi);
            if (!hit) {
                std::printf("  FAIL scale %.0e: missed a crossing at t=1/2\n", scale);
                ++failures; ++bad;
                continue;
            }
            const double err = (double)toi - expect;
            if (err > 1e-6) { fail("scale invariance (late)", (double)toi, expect); ++bad; continue; }

            if (scale == 1.0) {
                reference_error = -err;
                std::printf("  scale %.0e: toi=%.9f  err=%+.2e  (the reference)\n",
                            scale, (double)toi, err);
                continue;
            }
            // The scene is similar, so the answer should be too. Allow a small
            // factor for the search taking a different path, not a growth trend.
            const double ratio = -err / reference_error;
            if (ratio > kScaleDriftFactor) {
                std::printf("  FAIL scale %.0e: toi=%.9f is %.2e early, %.0fx the scale-1 result\n",
                            scale, (double)toi, -err, ratio);
                ++failures; ++bad;
            } else {
                std::printf("  scale %.0e: toi=%.9f  err=%+.2e  (%.1fx the reference)\n",
                            scale, (double)toi, err, ratio);
            }
        }
        return bad;
    }

    // 5. Slow motion.
    //
    // A vertex crossing a stationary quad at t = 1/2 crosses at t = 1/2 however
    // slowly it travels. The geometry below shrinks the vertex's displacement
    // over the step through eighteen orders of magnitude while leaving the quad,
    // and therefore the u and v Lipschitz constants, at unit scale.
    //
    // The per-axis tolerance is `codomain_tol / 3 / lipschitz[d]`, so a small
    // Lipschitz constant in t makes the time tolerance enormous -- and below
    // machine epsilon the quad path abandons the division and returns 1.0, a
    // tolerance as wide as the whole domain. `SCCD_MAX_TIME_TOL` and
    // `SCCD_MAX_COORD_TOL` exist to bound exactly that, and no quad file
    // references them.
    // ---------------------------------------------------------------------
    int check_slow_motion() {
        const S sq[4][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}};
        const double travels[] = {2.0, 2e-3, 2e-6, 2e-9, 2e-12, 2e-18};
        const double expect = 0.5;

        int bad = 0;
        for (const double d : travels) {
            Query q;
            for (int k = 0; k < 4; ++k) {
                for (int c = 0; c < 3; ++c) { q.s[k][c] = sq[k][c]; q.e[k][c] = sq[k][c]; }
            }
            // Straddles z = 0 symmetrically, so the crossing is at t = 1/2 for
            // every displacement.
            q.sv[0] = S(0.30); q.sv[1] = S(0.40); q.sv[2] = (S)( 0.5 * d);
            q.ev[0] = S(0.30); q.ev[1] = S(0.40); q.ev[2] = (S)(-0.5 * d);

            S toi = S(1);
            const bool hit = run_kernel(q, 69, S(1e-8), toi);
            if (!hit) {
                std::printf("  FAIL travel %.0e: missed a crossing at t=1/2\n", d);
                ++failures; ++bad;
                continue;
            }
            const double err = (double)toi - expect;
            if (err > 1e-6) { fail("slow motion (late)", (double)toi, expect); ++bad; }
            else if (err < -kSlowEarlyBudget) {
                std::printf("  FAIL travel %.0e: toi=%.9f is %.2e before the true 0.5\n",
                            d, (double)toi, -err);
                ++failures; ++bad;
            } else {
                std::printf("  travel %.0e: toi=%.9f  err=%+.2e\n", d, (double)toi, err);
            }
        }
        return bad;
    }

    // ---------------------------------------------------------------------
    // 6. The shipped entry point, not just the free root finder.
    //
    // Checks 1-5 call `find_root_grid_adaptive_split_vq` directly and build the
    // tolerances with the three separate `compute_vertex_quad_tolerance` /
    // `compute_vertex_quad_codomain_widths` / `vq_numerical_error` calls.
    // `narrow_phase_vq` instead goes through the one-pass `vq_prepare`, which
    // recomputes all three. Nothing else in the repository runs both forms on the
    // same query, so a divergence between them would be invisible.
    //
    // This also covers ToiOutput::Earliest, the shared-minimum path, which is the
    // only thing that exercises the quad atomics.
    // ---------------------------------------------------------------------
    int check_entry_point() {
        // Two vertices against one stationary unit square: one crosses at t=1/2,
        // one at t=1/4. The shared minimum is 1/4; the per-query answers are both.
        S px0[6] = {0, 1, 0, 1, S(0.30), S(0.50)};
        S py0[6] = {0, 0, 1, 1, S(0.40), S(0.50)};
        S pz0[6] = {0, 0, 0, 0, S(1.0),  S(1.0)};
        S px1[6], py1[6], pz1[6];
        for (int i = 0; i < 6; ++i) { px1[i] = px0[i]; py1[i] = py0[i]; pz1[i] = pz0[i]; }
        pz1[4] = S(-1.0);   // crosses z=0 at t=1/2
        pz1[5] = S(-3.0);   // crosses z=0 at t=1/4

        S* v0[3] = {px0, py0, pz0};
        S* v1[3] = {px1, py1, pz1};

        int q0[1] = {0}, q1[1] = {1}, q2[1] = {2}, q3[1] = {3};
        int* quads[4] = {q0, q1, q2, q3};

        const int vov[2] = {4, 5};
        const int qov[2] = {0, 0};

        int bad = 0;

        // ToiOutput::PerPair: one answer per candidate.
        {
            S toi[2] = {S(1), S(1)};
            sccd::narrow_phase_vq<S, int>(2, vov, qov, v0, v1, 1, quads,
                                             S(1), toi, 69, S(1e-8), sccd::ToiOutput::PerPair);
            const double expect[2] = {0.5, 0.25};
            for (int i = 0; i < 2; ++i) {
                const double err = (double)toi[i] - expect[i];
                if (err > 1e-6 || err < -kEntryEarlyBudget) {
                    std::printf("  FAIL narrow_phase_vq stride 1, query %d: toi=%.9f, true=%.3f\n",
                                i, (double)toi[i], expect[i]);
                    ++failures; ++bad;
                } else {
                    std::printf("  narrow_phase_vq stride 1, query %d: toi=%.9f  err=%+.2e\n",
                                i, (double)toi[i], err);
                }
            }
        }

        // ToiOutput::Earliest: the shared minimum over both candidates.
        {
            S toi = S(1);
            sccd::narrow_phase_vq<S, int>(2, vov, qov, v0, v1, 1, quads,
                                             S(1), &toi, 69, S(1e-8), sccd::ToiOutput::Earliest);
            const double err = (double)toi - 0.25;
            if (err > 1e-6 || err < -kEntryEarlyBudget) {
                std::printf("  FAIL narrow_phase_vq stride 0: toi=%.9f, true=%.3f\n",
                            (double)toi, 0.25);
                ++failures; ++bad;
            } else {
                std::printf("  narrow_phase_vq stride 0: toi=%.9f  err=%+.2e\n", (double)toi, err);
            }
        }
        return bad;
    }

}  // namespace

int main() {
    std::printf("1. factored evaluation vs the reference\n");
    check_evaluation_agrees(400, 40);

    std::printf("2. analytic times of impact\n");
    check_analytic();

    std::printf("3. random queries with a planted root\n");
    check_planted_roots(2000);

    std::printf("4. scale invariance\n");
    check_scale_invariance();

    std::printf("5. slow motion\n");
    check_slow_motion();

    std::printf("6. the shipped entry point\n");
    check_entry_point();

    std::printf("%s\n", failures ? "FAIL" : "OK: the vertex-quad root finder is conservative on every check");
    return failures ? 1 : 0;
}
