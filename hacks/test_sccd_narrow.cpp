#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/catch_approx.hpp>

#include "io.hpp"
#include "ground_truth.hpp"

#include <scalable_ccd/config.hpp>

#include <scalable_ccd/utils/logger.hpp>
#include <scalable_ccd/utils/profiler.hpp>

#include "sccd.hpp"
#include "conversions.hpp"

#include <filesystem>
namespace fs = std::filesystem;

TEST_CASE("Test Lean Narrow phase", "[narrow_phase]")
{
    using namespace scalable_ccd;

    const fs::path data(SCALABLE_CCD_DATA_DIR);

    fs::path file_t0, file_t1, vf_ground_truth, ee_ground_truth;

    int CASE = 0;
    SCCD_READ_ENV(CASE, atoi);

    switch (CASE) {
    case 0: {
        file_t0 = data / "cloth-ball" / "frames" / "cloth_ball92.ply";
        file_t1 = data / "cloth-ball" / "frames" / "cloth_ball93.ply";
        vf_ground_truth = data / "cloth-ball" / "boxes" / "92vf.json";
        ee_ground_truth = data / "cloth-ball" / "boxes" / "92ee.json";
        break;
    }
    case 1: {
        file_t0 = data / "rod-twist" / "frames" / "3036.ply";
        file_t1 = data / "rod-twist" / "frames" / "3037.ply";
        vf_ground_truth = data / "rod-twist" / "boxes" / "3036vf.json";
        ee_ground_truth = data / "rod-twist" / "boxes" / "3036ee.json";
        break;
    }
    case 2: {
        file_t0 = data / "armadillo-rollers" / "frames" / "326.ply";
        file_t1 = data / "armadillo-rollers" / "frames" / "327.ply";
        vf_ground_truth = data / "armadillo-rollers" / "boxes" / "326vf.json";
        ee_ground_truth = data / "armadillo-rollers" / "boxes" / "326ee.json";
        break;
    }
    case 3: {
        file_t0 = data / "cloth-funnel" / "frames" / "227.ply";
        file_t1 = data / "cloth-funnel" / "frames" / "228.ply";
        vf_ground_truth = data / "cloth-funnel" / "boxes" / "227vf.json";
        ee_ground_truth = data / "cloth-funnel" / "boxes" / "227ee.json";
        break;
    }
    case 4: {
        file_t0 = data / "n-body-simulation" / "frames" / "balls16_18.ply";
        file_t1 = data / "n-body-simulation" / "frames" / "balls16_19.ply";
        vf_ground_truth = data / "n-body-simulation" / "boxes" / "18vf.json";
        ee_ground_truth = data / "n-body-simulation" / "boxes" / "18ee.json";
        break;
    }
    case 5: {
        file_t0 = data / "puffer-ball" / "frames" / "20.ply";
        file_t1 = data / "puffer-ball" / "frames" / "21.ply";
        vf_ground_truth = "";
        ee_ground_truth = "";
        break;
    }
    default: {
        throw std::runtime_error("Invalid case");
        break;
    }
    }

    // ------------------------------------------------------------------------
    // Load meshes

    Eigen::MatrixXd vertices_t0, vertices_t1;
    Eigen::MatrixXi edges, faces;
    parse_mesh(file_t0, file_t1, vertices_t0, vertices_t1, faces, edges);

    std::vector<AABB> vertex_boxes, edge_boxes, face_boxes;
    build_vertex_boxes(vertices_t0, vertices_t1, vertex_boxes);
    build_edge_boxes(vertex_boxes, edges, edge_boxes);
    build_face_boxes(vertex_boxes, faces, face_boxes);

    constexpr bool allow_zero_toi = true;
    constexpr Scalar min_distance = 0;
    constexpr int max_iterations = -1;
    constexpr Scalar tolerance = 1e-6;
    constexpr Scalar memory_limit_GB = 0;

    std::vector<std::tuple<int, int, Scalar>> collisions;

    Timer timer_ccd;

    sccd::SCCD_t ccd_algo;
    ccd_algo.init(vertex_boxes, edge_boxes, face_boxes);
    ccd_algo.init_point_data(vertices_t0, vertices_t1);

    timer_ccd.start();
    ccd_algo.broad_phase();
    timer_ccd.stop();

    double tts_broad_phase = timer_ccd.getElapsedTimeInMilliSec();

    printf(
        "Broad time: %g [ms]\n", tts_broad_phase);

    timer_ccd.start();
    Scalar toi = ccd_algo.narrow_phase();
    timer_ccd.stop();

    double tts_narrow_phase = timer_ccd.getElapsedTimeInMilliSec();
    printf(
        "Narrow time: %g [ms]\n", tts_narrow_phase);

    printf(
        "Total time: %g [ms]\n", tts_narrow_phase + tts_broad_phase);

    printf("Toi: %g\n", toi);

    ccd_algo.export_narrowphase_results(collisions);

    //     Scalar toi =
    //         ccd(vertices_t0, vertices_t1, edges, faces, min_distance,
    //             max_iterations, tolerance, allow_zero_toi,
    // #ifdef SCALABLE_CCD_TOI_PER_QUERY
    //             collisions,
    // #endif
    //             memory_limit_GB);

    for (const auto& [i, j, _toi] : collisions) {
        CHECK(toi <= _toi);
    }

    // if(!CASE)
    //     CHECK(toi == Catch::Approx(3.814697265625e-06));

#ifdef SCALABLE_CCD_WITH_PROFILER
    profiler().data()["memory_limit_GB"] = memory_limit_GB;
    profiler().data()["toi"] = toi;
    profiler().print();
#endif
}