#include <catch2/catch_test_macros.hpp>

#include "io.hpp"
#include "ground_truth.hpp"

#include <scalable_ccd/broad_phase/sort_and_sweep.hpp>
#include <scalable_ccd/utils/timer.hpp>
#include <scalable_ccd/utils/logger.hpp>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <iostream>

#include "sccd.hpp"
#include "conversions.hpp"

namespace fs = std::filesystem;

TEST_CASE("Test CPU broad phase", "[cpu][broad_phase]")
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

    std::vector<AABB> vertex_boxes, edge_boxes, face_boxes;
    parse_mesh(file_t0, file_t1, vertex_boxes, edge_boxes, face_boxes);

    // ------------------------------------------------------------------------
    // Run SCCD

    const size_t nfaces = face_boxes.size();
    const size_t nedges = edge_boxes.size();
    const size_t nnodes = vertex_boxes.size();

    printf(
        "faces %ld, vertices %ld, edges: %ld\n", face_boxes.size(),
        vertex_boxes.size(), edge_boxes.size());


    Timer timer;
    timer.start();

    Timer timer_ccd;

    sccd::SCCD_t ccd_algo;
    ccd_algo.init(vertex_boxes, edge_boxes, face_boxes);

    timer_ccd.start();
    ccd_algo.broad_phase();
    // ccd_algo.find_with_cell_list();
    timer_ccd.stop();

    std::vector<std::pair<int, int>> lean_vf_overlaps, lean_ee_overlaps;
    ccd_algo.export_broadphase_results(lean_vf_overlaps, lean_ee_overlaps);

    timer.stop();

    double lean_time = timer.getElapsedTimeInMilliSec();


    printf(
        "Lean time (%lu, %lu): %g [ms], ccd only %g [ms] \n",
        lean_vf_overlaps.size(), lean_ee_overlaps.size(), lean_time,
        (double)timer_ccd.getElapsedTimeInMilliSec());

    // ------------------------------------------------------------------------
    // Run original

    timer.start();
    int sort_axis = ccd_algo.sort_axis;
    std::vector<std::pair<int, int>> vf_overlaps;
    sort_and_sweep(vertex_boxes, face_boxes, sort_axis, vf_overlaps);
    CHECK(sort_axis == ccd_algo.sort_axis); // check output sort axis

    sort_axis = ccd_algo.sort_axis; // Reset sort axis
    std::vector<std::pair<int, int>> ee_overlaps;
    sort_and_sweep(edge_boxes, sort_axis, ee_overlaps);
    CHECK(sort_axis == ccd_algo.sort_axis); // check output sort axis

    timer.stop();

    double orignal_time = (double)timer.getElapsedTimeInMilliSec();
    printf(
        "Elapsed time (%lu, %lu): %g [ms]\n", vf_overlaps.size(),
        ee_overlaps.size(), orignal_time);

    logger().trace("Elapsed time: {:.6f} ms", orignal_time);
    printf("Lean/Original speed up %gx\n", orignal_time / lean_time);

    // ------------------------------------------------------------------------
    // Compare

    // Sort to check if they are the same
    std::sort(vf_overlaps.begin(), vf_overlaps.end());
    std::sort(lean_vf_overlaps.begin(), lean_vf_overlaps.end());
    std::sort(lean_ee_overlaps.begin(), lean_ee_overlaps.end());
    std::sort(ee_overlaps.begin(), ee_overlaps.end());

    for (size_t i = 0; i < lean_vf_overlaps.size(); i++) {
        if (vf_overlaps[i] != lean_vf_overlaps[i]) {
            printf(
                "[%lu] (%d, %d) == (%d, %d)\n", i, vf_overlaps[i].first,
                vf_overlaps[i].second, lean_vf_overlaps[i].first,
                lean_vf_overlaps[i].second);
            assert(false);
        }
    }

    CHECK(vf_overlaps == lean_vf_overlaps);
    CHECK(ee_overlaps == lean_ee_overlaps);

    // Offset the boxes to match the way ground truth was originally generated.
    int offset = vertex_boxes.size();
    for (auto& [a, b] : ee_overlaps) {
        a += offset;
        b += offset;
    }
    offset += edge_boxes.size();
    for (auto& [v, f] : vf_overlaps) {
        f += offset;
    }

    if (!vf_ground_truth.empty()) {
        compare_mathematica(vf_overlaps, vf_ground_truth);
    }
    if (!ee_ground_truth.empty()) {
        compare_mathematica(ee_overlaps, ee_ground_truth);
    }
}