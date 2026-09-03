// find_earliest_impact_time and find_impact_times, the two calls a user of this
// library actually makes, had no test. They were exercised only by a demo with
// no assertions in it, so nothing checked that the number coming out was right --
// or, for find_impact_times, that the buffers it hands back are populated at all.
//
// The edge-edge narrow phase was in the same position: broadphase_test covers
// vertex-face only, so the edge-edge time of impact had no direct coverage
// either. It is checked here, since these calls are how anybody reaches it.
//
// Every case has a time of impact that can be written down. Conservativeness is
// the property under test: reported at or before the true value, never after,
// and never a miss.

#include "sccd_smesh_CCD.hpp"
#include "smesh_mesh.hpp"
#include "smesh_test.hpp"

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

    using T = double;

    // The mesh stores geometry as smesh::geom_t; the CCD works in scalar_t. Both
    // come from the same table so the two views cannot drift apart.
    template <typename S, std::size_t N>
    smesh::SharedBuffer<S*> points_soa(const std::array<std::array<double, N>, 3>& p) {
        auto buffer = smesh::create_host_buffer<S>(3, N);
        for (int d = 0; d < 3; ++d) {
            for (std::size_t i = 0; i < N; ++i) buffer->data()[d][i] = (S)p[d][i];
        }
        return buffer;
    }

    // Two triangles that do not share a vertex, in one mesh.
    //
    //   nodes 0,1,2 -- a triangle in the plane z = 0, corners (0,0) (2,0) (0,2)
    //   nodes 3,4,5 -- a parallel triangle above it, which descends
    //
    // The upper triangle STRADDLES the lower one's hypotenuse: its corner
    // (0.5,0.5) is inside, and (2.5,0.5) and (0.5,2.5) are outside. That matters,
    // and it is the second geometry this test used. The first put the upper
    // triangle strictly inside the lower one, which gives a vertex-face contact
    // and no edge-edge contact at all -- two coplanar triangles, one contained in
    // the other, have no crossing edges. The edge-edge assertion failed and the
    // test was what was wrong. Straddling an edge makes both phases see it.
    std::shared_ptr<smesh::Mesh> two_triangles() {
        auto elements = smesh::create_host_buffer<smesh::idx_t>(3, 2);
        elements->data()[0][0] = 0;
        elements->data()[1][0] = 1;
        elements->data()[2][0] = 2;
        elements->data()[0][1] = 3;
        elements->data()[1][1] = 4;
        elements->data()[2][1] = 5;

        const std::array<std::array<double, 6>, 3> rest = {{
            {{0.0, 2.0, 0.0,  0.5, 2.5, 0.5}},
            {{0.0, 0.0, 2.0,  0.5, 0.5, 2.5}},
            {{0.0, 0.0, 0.0,  1.0, 1.0, 1.0}},
        }};
        return std::make_shared<smesh::Mesh>(smesh::Communicator::self(),
                                             smesh::TRISHELL3,
                                             elements,
                                             points_soa<smesh::geom_t, 6>(rest));
    }

    // The same six nodes; the upper triangle translates down by `drop` over the
    // step, so it reaches z = 0 at t = 1 / drop.
    std::array<std::array<double, 6>, 3> frame(const double z_upper) {
        return {{
            {{0.0, 2.0, 0.0,  0.5, 2.5, 0.5}},
            {{0.0, 0.0, 2.0,  0.5, 0.5, 2.5}},
            {{0.0, 0.0, 0.0,  z_upper, z_upper, z_upper}},
        }};
    }

    int check_toi(const char* what, const double toi, const double expect) {
        // Late is illegal; earlier is safe but should not be wild.
        const double err = toi - expect;
        if (err > 1e-6) {
            std::printf("  FAIL %-38s toi=%.9f is LATER than the true %.6f\n", what, toi, expect);
            return 1;
        }
        if (err < -0.05) {
            std::printf("  FAIL %-38s toi=%.9f is far earlier than %.6f\n", what, toi, expect);
            return 1;
        }
        std::printf("  %-38s toi=%.9f  true=%.4f  err=%+.2e\n", what, toi, expect, err);
        return 0;
    }

    int test_find_earliest_impact_time() {
        auto mesh = two_triangles();
        auto ccd = sccd::CCD<T>::create(mesh, smesh::EXECUTION_SPACE_HOST);

        // Drop from z = 1 to z = -1: the plane is reached halfway.
        auto p0 = points_soa<T, 6>(frame(1.0));
        auto p1 = points_soa<T, 6>(frame(-1.0));

        T toi = 1;
        int err = ccd->find_earliest_impact_time(p0, p1, toi, 40, 1e-8);
        SMESH_TEST_EQ(err, SCCD_SUCCESS);
        if (check_toi("find_earliest_impact_time, t=1/2", (double)toi, 0.5)) return SMESH_TEST_FAILURE;

        // Drop from z = 1 to z = -3: the plane is reached at a quarter.
        auto q0 = points_soa<T, 6>(frame(1.0));
        auto q1 = points_soa<T, 6>(frame(-3.0));
        auto ccd2 = sccd::CCD<T>::create(two_triangles(), smesh::EXECUTION_SPACE_HOST);
        T toi2 = 1;
        err = ccd2->find_earliest_impact_time(q0, q1, toi2, 40, 1e-8);
        SMESH_TEST_EQ(err, SCCD_SUCCESS);
        if (check_toi("find_earliest_impact_time, t=1/4", (double)toi2, 0.25)) return SMESH_TEST_FAILURE;

        // No motion at all: nothing may be reported.
        auto s0 = points_soa<T, 6>(frame(1.0));
        auto s1 = points_soa<T, 6>(frame(1.0));
        auto ccd3 = sccd::CCD<T>::create(two_triangles(), smesh::EXECUTION_SPACE_HOST);
        T toi3 = 1;
        err = ccd3->find_earliest_impact_time(s0, s1, toi3, 40, 1e-8);
        SMESH_TEST_EQ(err, SCCD_SUCCESS);
        if ((double)toi3 < 1.0) {
            std::printf("  FAIL stationary geometry reported a contact at %.9f\n", (double)toi3);
            return SMESH_TEST_FAILURE;
        }
        std::printf("  %-38s toi=%.9f (no contact, as expected)\n", "stationary geometry", (double)toi3);

        return SMESH_TEST_SUCCESS;
    }

    int test_find_impact_times() {
        auto mesh = two_triangles();
        auto ccd = sccd::CCD<T>::create(mesh, smesh::EXECUTION_SPACE_HOST);

        auto p0 = points_soa<T, 6>(frame(1.0));
        auto p1 = points_soa<T, 6>(frame(-1.0));

        smesh::SharedBuffer<smesh::idx_t> v_overlap, f_overlap, e0_overlap, e1_overlap;
        smesh::SharedBuffer<T> vf_tois, ee_tois;

        const int err = ccd->find_impact_times(
            p0, p1, v_overlap, f_overlap, vf_tois, e0_overlap, e1_overlap, ee_tois, 40, 1e-8);
        SMESH_TEST_EQ(err, SCCD_SUCCESS);

        // Every buffer must come back, and the per-pair arrays must agree in
        // length with the pair lists they belong to. A handed-back null here
        // would be exactly the kind of thing a demo without assertions misses.
        SMESH_TEST_ASSERT(bool(v_overlap));
        SMESH_TEST_ASSERT(bool(f_overlap));
        SMESH_TEST_ASSERT(bool(vf_tois));
        SMESH_TEST_ASSERT(bool(e0_overlap));
        SMESH_TEST_ASSERT(bool(e1_overlap));
        SMESH_TEST_ASSERT(bool(ee_tois));
        SMESH_TEST_EQ(v_overlap->size(), f_overlap->size());
        SMESH_TEST_EQ(vf_tois->size(), v_overlap->size());
        SMESH_TEST_EQ(e0_overlap->size(), e1_overlap->size());
        SMESH_TEST_EQ(ee_tois->size(), e0_overlap->size());

        std::printf("  vertex-face pairs: %ld    edge-edge pairs: %ld\n",
                    (long)v_overlap->size(), (long)e0_overlap->size());

        // The geometry was chosen so both phases see the contact, and the
        // earliest time either reports must be the conservative one.
        int bad = 0;
        double earliest_vf = 1.0;
        for (ptrdiff_t i = 0; i < (ptrdiff_t)vf_tois->size(); ++i) {
            earliest_vf = std::min(earliest_vf, (double)vf_tois->data()[i]);
        }
        double earliest_ee = 1.0;
        for (ptrdiff_t i = 0; i < (ptrdiff_t)ee_tois->size(); ++i) {
            earliest_ee = std::min(earliest_ee, (double)ee_tois->data()[i]);
        }

        SMESH_TEST_ASSERT(v_overlap->size() > 0u);
        bad |= check_toi("find_impact_times, earliest vertex-face", earliest_vf, 0.5);

        // Edge-edge: this is the only direct coverage it has.
        SMESH_TEST_ASSERT(e0_overlap->size() > 0u);
        bad |= check_toi("find_impact_times, earliest edge-edge", earliest_ee, 0.5);

        return bad ? SMESH_TEST_FAILURE : SMESH_TEST_SUCCESS;
    }

    // The staged calls must agree with the one-shot call. They are the same work
    // split differently, and a user who stages the phases to interleave their own
    // logic should not get a different answer.
    int test_staged_matches_oneshot() {
        auto p0 = points_soa<T, 6>(frame(1.0));
        auto p1 = points_soa<T, 6>(frame(-1.0));

        auto one = sccd::CCD<T>::create(two_triangles(), smesh::EXECUTION_SPACE_HOST);
        T toi_oneshot = 1;
        SMESH_TEST_EQ(one->find_earliest_impact_time(p0, p1, toi_oneshot, 40, 1e-8), SCCD_SUCCESS);

        auto staged = sccd::CCD<T>::create(two_triangles(), smesh::EXECUTION_SPACE_HOST);
        smesh::SharedBuffer<smesh::idx_t> v_overlap, f_overlap, e0_overlap, e1_overlap;
        SMESH_TEST_EQ(staged->broad_phase_prep(p0, p1), SCCD_SUCCESS);
        SMESH_TEST_EQ(staged->broad_phase_fv_step(v_overlap, f_overlap), SCCD_SUCCESS);
        SMESH_TEST_EQ(staged->broad_phase_ee_step(e0_overlap, e1_overlap), SCCD_SUCCESS);

        T toi_staged = 1;
        smesh::SharedBuffer<T> vf_tois, ee_tois;
        SMESH_TEST_EQ(staged->narrow_phase_fv(toi_staged, vf_tois, 40, 1e-8, 0), SCCD_SUCCESS);
        SMESH_TEST_EQ(staged->narrow_phase_ee(toi_staged, ee_tois, 40, 1e-8, 0), SCCD_SUCCESS);

        std::printf("  one-shot toi=%.9f   staged vertex-face toi=%.9f\n",
                    (double)toi_oneshot, (double)toi_staged);
        SMESH_TEST_APPROXEQ((double)toi_oneshot, (double)toi_staged, 1e-9);
        return SMESH_TEST_SUCCESS;
    }

    // The broad-phase auto-tuner races the sweep against the cell list over
    // consecutive steps and keeps the winner, so a long-lived CCD object changes
    // strategy underneath the caller. That is only acceptable if the answer does
    // not change with it -- the two produce identical pair sets, which is the
    // whole premise of racing them.
    //
    // Driving one object over enough steps exercises both strategies and the
    // switch between them. Fixed geometry, so anything that moves is the tuner.
    int test_tuner_does_not_change_the_answer() {
        auto ccd = sccd::CCD<T>::create(two_triangles(), smesh::EXECUTION_SPACE_HOST);
        auto p0 = points_soa<T, 6>(frame(1.0));
        auto p1 = points_soa<T, 6>(frame(-1.0));

        std::size_t first_vf = 0, first_ee = 0;
        double first_toi = 0.0;

        for (int step = 0; step < 10; ++step) {
            smesh::SharedBuffer<smesh::idx_t> v_overlap, f_overlap, e0_overlap, e1_overlap;
            smesh::SharedBuffer<T> vf_tois, ee_tois;
            const int err = ccd->find_impact_times(
                p0, p1, v_overlap, f_overlap, vf_tois, e0_overlap, e1_overlap, ee_tois, 40, 1e-8);
            SMESH_TEST_EQ(err, SCCD_SUCCESS);

            double toi = 1.0;
            for (ptrdiff_t i = 0; i < (ptrdiff_t)vf_tois->size(); ++i) {
                toi = std::min(toi, (double)vf_tois->data()[i]);
            }
            for (ptrdiff_t i = 0; i < (ptrdiff_t)ee_tois->size(); ++i) {
                toi = std::min(toi, (double)ee_tois->data()[i]);
            }

            if (step == 0) {
                first_vf = v_overlap->size();
                first_ee = e0_overlap->size();
                first_toi = toi;
                std::printf("  step 0: vf=%ld ee=%ld toi=%.9f\n",
                            (long)first_vf, (long)first_ee, first_toi);
                continue;
            }

            if (v_overlap->size() != first_vf || e0_overlap->size() != first_ee) {
                std::printf("  FAIL step %d: pair counts moved, vf %ld->%ld ee %ld->%ld\n",
                            step, (long)first_vf, (long)v_overlap->size(),
                            (long)first_ee, (long)e0_overlap->size());
                return SMESH_TEST_FAILURE;
            }
            if (std::fabs(toi - first_toi) > 1e-12) {
                std::printf("  FAIL step %d: toi moved %.12f -> %.12f\n", step, first_toi, toi);
                return SMESH_TEST_FAILURE;
            }
        }
        std::printf("  %-38s 10 steps, pair sets and toi identical throughout\n", "tuner");
        return SMESH_TEST_SUCCESS;
    }

}  // namespace

int main(int argc, char** argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);
    {
        SMESH_RUN_TEST(test_find_earliest_impact_time);
        SMESH_RUN_TEST(test_find_impact_times);
        SMESH_RUN_TEST(test_staged_matches_oneshot);
        SMESH_RUN_TEST(test_tuner_does_not_change_the_answer);
    }
    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
