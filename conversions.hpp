#ifndef CONVERSION_HPP
#define CONVERSION_HPP

#include <scalable_ccd/broad_phase/sort_and_sweep.hpp>
#include <scalable_ccd/utils/logger.hpp>
#include <scalable_ccd/utils/timer.hpp>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

using namespace scalable_ccd;

#include "cell_broadphase.hpp"

namespace sccd {

// (Ugly) Conversion between Scalable-CCD and sccd and back

typedef struct SCCD {
  int sort_axis;
  size_t nfaces;
  size_t nedges;
  size_t nnodes;

  geom_t *faabb[6];
  geom_t *eaabb[6];
  geom_t *vaabb[6];

  idx_t *soafaces[3];
  idx_t *soaedges[2];

  // FV overalps
  std::vector<idx_t> foverlap;
  std::vector<idx_t> voverlap;

  // Edge overlaps
  std::vector<idx_t> e0_overlap, e1_overlap;

  // Indices
  std::vector<idx_t> fidx;
  std::vector<idx_t> vidx;
  std::vector<idx_t> eidx;

  // CCD poter
  std::vector<size_t> ccdptr;

  // Scratch space for sorting
  std::vector<geom_t> scratch;

  // Workspace allocations
  std::vector<std::vector<geom_t>> face_aabb;
  std::vector<std::vector<geom_t>> edge_aabb;
  std::vector<std::vector<geom_t>> vertex_aabb;

  std::vector<idx_t> faces;
  std::vector<idx_t> edges;
  bool verbose{true};

  void init(const std::vector<AABB> &vertex_boxes,
            const std::vector<AABB> &edge_boxes,
            const std::vector<AABB> &face_boxes) {
    nfaces = face_boxes.size();
    nedges = edge_boxes.size();
    nnodes = vertex_boxes.size();

    fidx.resize(nfaces);
    vidx.resize(nnodes);
    eidx.resize(nedges);

    scratch.resize(std::max(nfaces, std::max(nedges, nnodes)));

    size_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
    ccdptr.resize(max_ccdptr_size);

    face_aabb.resize(6);
    for (auto &v : face_aabb) {
      v.resize(nfaces);
    }

    edge_aabb.resize(6);
    for (auto &v : edge_aabb) {
      v.resize(nedges);
    }

    vertex_aabb.resize(6);
    for (auto &v : vertex_aabb) {
      v.resize(nnodes);
    }

    // C++ to C conversion
    for (int d = 0; d < 6; d++) {
      faabb[d] = face_aabb[d].data();
      eaabb[d] = edge_aabb[d].data();
      vaabb[d] = vertex_aabb[d].data();
    }

    for (size_t i = 0; i < nfaces; i++) {
      for (int d = 0; d < 3; d++) {
        faabb[d][i] = face_boxes[i].min[d];
        faabb[3 + d][i] = face_boxes[i].max[d];
      }
    }

    for (size_t i = 0; i < nnodes; i++) {
      for (int d = 0; d < 3; d++) {
        vaabb[d][i] = vertex_boxes[i].min[d];
        vaabb[3 + d][i] = vertex_boxes[i].max[d];
      }
    }

    for (size_t i = 0; i < nedges; i++) {
      for (int d = 0; d < 3; d++) {
        eaabb[d][i] = edge_boxes[i].min[d];
        eaabb[3 + d][i] = edge_boxes[i].max[d];
      }
    }

    faces.resize(nfaces * 3);
    for (size_t i = 0; i < nfaces; i++) {
      for (int d = 0; d < 3; d++) {
        faces[i * 3 + d] = face_boxes[i].vertex_ids[d];
      }
    }

    edges.resize(nedges * 2);
    for (size_t i = 0; i < nedges; i++) {
      for (int d = 0; d < 2; d++) {
        edges[i * 2 + d] = edge_boxes[i].vertex_ids[d];
      }
    }

    soafaces[0] = &faces[0];
    soafaces[1] = &faces[1];
    soafaces[2] = &faces[2];

    soaedges[0] = &edges[0];
    soaedges[1] = &edges[1];

    sort_axis = -1;
  }

  // void find_with_starts()
  // {
  //     Timer timer;

  //     timer.start();

  //     size_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
  //     ccdptr.resize(max_ccdptr_size);

  //     int axes[3];
  //     largest_variance_axes_sort(nnodes, vaabb, axes);

  //     sort_axis = axes[0];
  //     sort_along_axis(nfaces, sort_axis, faabb, fidx.data(), scratch.data());
  //     sort_along_axis(nnodes, sort_axis, vaabb, vidx.data(), scratch.data());
  //     sort_along_axis(nedges, sort_axis, eaabb, eidx.data(), scratch.data());

  //     timer.stop();
  //     if (verbose)
  //         printf(
  //             "SCCD, Sorting: %g [ms]\n", timer.getElapsedTimeInMilliSec());
  //     timer.start();

  //     int cell_list_axis = sort_axis;
  //     size_t ncells = 1024; // Max amount
  //     geom_t cell_min;
  //     geom_t cell_size;
  //     cell_starts_setup(
  //         nnodes, vaabb[cell_list_axis], vaabb[cell_list_axis + 3], ncells,
  //         &cell_min, &cell_size);
  //     std::vector<size_t> starts(ncells);
  //     cell_starts(
  //         ncells, cell_min, cell_size, nnodes, vaabb[cell_list_axis],
  //         starts.data());

  //     timer.stop();
  //     if (verbose)
  //         printf(
  //             "SCCD, Cell Start(%lu): %g [ms]\n", ncells,
  //             timer.getElapsedTimeInMilliSec());
  //     timer.start();

  //     // F2V
  //     count_overlaps_with_starts<3, 1>(
  //         sort_axis, nfaces, faabb, fidx.data(), 3, soafaces, nnodes, vaabb,
  //         vidx.data(), 0, nullptr, ncells, cell_min, cell_size,
  //         starts.data(), ccdptr.data());

  //     timer.stop();
  //     if (verbose)
  //         printf(
  //             "SCCD count cell starts(%lu), F2V: %g [ms]\n", ccdptr[nfaces],
  //             timer.getElapsedTimeInMilliSec());
  // }

  void find_with_cell_list() {
    Timer timer;

    timer.start();

    size_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
    ccdptr.resize(max_ccdptr_size);

    int axes[3];
    largest_variance_axes_sort(nnodes, vaabb, axes);

    sort_axis = axes[0];
    sort_along_axis(nfaces, sort_axis, faabb, fidx.data(), scratch.data());
    sort_along_axis(nnodes, sort_axis, vaabb, vidx.data(), scratch.data());
    sort_along_axis(nedges, sort_axis, eaabb, eidx.data(), scratch.data());

    int cell_list_axis = axes[1];
    size_t ncells = 2048; // Max amount
    geom_t cell_min;
    geom_t cell_size;
    cell_setup(nnodes, vaabb[cell_list_axis], vaabb[cell_list_axis + 3],
               &ncells, &cell_min, &cell_size);
    std::vector<idx_t> cellptr(ncells + 1), bookkeeping(ncells);
    cell_count(ncells, cell_min, cell_size, nnodes, vaabb[cell_list_axis],
               cellptr.data());

    std::vector<idx_t> cellidx(cellptr[ncells]);
    cell_populate(ncells, cell_min, cell_size, nnodes, vaabb[cell_list_axis],
                  cellptr.data(), cellidx.data(), bookkeeping.data());

    timer.stop();
    if (verbose)
      printf("SCCD(CELL), Sorting: %g [ms]\n",
             timer.getElapsedTimeInMilliSec());
    timer.start();

    // E2E
    std::fill(ccdptr.begin(), ccdptr.end(), 0);

    count_self_overlaps<2>(sort_axis, nedges, eaabb, eidx.data(), 2, soaedges,
                           ccdptr.data());

    const size_t ee_n_intersections = ccdptr[nedges];
    e0_overlap.resize(ee_n_intersections);
    e1_overlap.resize(ee_n_intersections);

    collect_self_overlaps<2>(sort_axis, nedges, eaabb, eidx.data(), 2, soaedges,
                             ccdptr.data(), e0_overlap.data(),
                             e1_overlap.data());

    timer.stop();
    if (verbose)
      printf("SCCD, E2E: %g [ms]\n", timer.getElapsedTimeInMilliSec());
    timer.start();

    Timer co;
    co.start();
    // F2V
    cell_count_overlaps<3, 1>(sort_axis, nfaces, faabb, fidx.data(), 3,
                              soafaces, nnodes, vaabb, vidx.data(), 0, nullptr,
                              cell_list_axis, ncells, cell_min, cell_size,
                              cellptr.data(), cellidx.data(), ccdptr.data());

    co.stop();
    printf("SCCD(CELL): F2V Count %g [ms]\n", co.getElapsedTimeInMilliSec());

    // Allocation (expensive, could be expanded dynamically in CCD)
    const size_t fv_nintersections = ccdptr[nfaces];
    foverlap.resize(fv_nintersections);
    voverlap.resize(fv_nintersections);

    cell_collect_overlaps<3, 1>(
        sort_axis, nfaces, faabb, fidx.data(), 3, soafaces, nnodes, vaabb,
        vidx.data(), 0, nullptr, cell_list_axis, ncells, cell_min, cell_size,
        cellptr.data(), cellidx.data(), ccdptr.data(), foverlap.data(),
        voverlap.data());

    timer.stop();
    if (verbose)
      printf("SCCD(CELL), F2V: %g [ms]\n", timer.getElapsedTimeInMilliSec());
  }

  void broad_phase() {
    Timer timer;

    timer.start();

    size_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
    ccdptr.resize(max_ccdptr_size);

    sort_axis = choose_axis(nnodes, vaabb);
    sort_along_axis(nfaces, sort_axis, faabb, fidx.data(), scratch.data());
    sort_along_axis(nnodes, sort_axis, vaabb, vidx.data(), scratch.data());
    sort_along_axis(nedges, sort_axis, eaabb, eidx.data(), scratch.data());

    timer.stop();
    if (verbose)
      printf("SCCD, Sorting: %g [ms]\n", timer.getElapsedTimeInMilliSec());
    timer.start();

    // E2E
    std::fill(ccdptr.begin(), ccdptr.end(), 0);

    count_self_overlaps<2>(sort_axis, nedges, eaabb, eidx.data(), 2, soaedges,
                           ccdptr.data());

    const size_t ee_n_intersections = ccdptr[nedges];
    e0_overlap.resize(ee_n_intersections);
    e1_overlap.resize(ee_n_intersections);

    collect_self_overlaps<2>(sort_axis, nedges, eaabb, eidx.data(), 2, soaedges,
                             ccdptr.data(), e0_overlap.data(),
                             e1_overlap.data());

    timer.stop();
    if (verbose)
      printf("SCCD, E2E: %g [ms]\n", timer.getElapsedTimeInMilliSec());
    timer.start();

    Timer co;
    co.start();
    // F2V
    count_overlaps<3, 1>(sort_axis, nfaces, faabb, fidx.data(), 3, soafaces,
                         nnodes, vaabb, vidx.data(), 0, nullptr, ccdptr.data());

    co.stop();
    printf("SCCD: F2V Count %g [ms]\n", co.getElapsedTimeInMilliSec());

    // Allocation (expensive, could be expanded dynamically in CCD)
    const size_t fv_nintersections = ccdptr[nfaces];
    foverlap.resize(fv_nintersections);
    voverlap.resize(fv_nintersections);

    collect_overlaps<3, 1>(sort_axis, nfaces, faabb, fidx.data(), 3, soafaces,
                           nnodes, vaabb, vidx.data(), 0, nullptr,
                           ccdptr.data(), foverlap.data(), voverlap.data());

    timer.stop();
    if (verbose)
      printf("SCCD, F2V: %g [ms]\n", timer.getElapsedTimeInMilliSec());
  }

  void
  export_broadphase_results(std::vector<std::pair<int, int>> &vf_overlaps,
                            std::vector<std::pair<int, int>> &ee_overlaps) {
    const size_t fv_nintersections = foverlap.size();
    vf_overlaps.resize(fv_nintersections);
    for (ptrdiff_t i = 0; i < fv_nintersections; i++) {
      vf_overlaps[i].first = voverlap[i];
      vf_overlaps[i].second = foverlap[i];
    }

    const size_t ee_n_intersections = e0_overlap.size();
    ee_overlaps.resize(ee_n_intersections);
    for (ptrdiff_t i = 0; i < ee_n_intersections; i++) {
      ee_overlaps[i].first = e0_overlap[i];
      ee_overlaps[i].second = e1_overlap[i];
    }
  }

  void narrow_phase() {
    // TODO
  }

  void export_narrowphase_results(
      std::vector<std::tuple<int, int, Scalar>> &collisions) {
    // TODO
  }

} SCCD_t;
} // namespace sccd

#endif // CONVERSION_HPP
