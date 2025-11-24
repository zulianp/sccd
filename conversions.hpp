#ifndef CONVERSION_HPP
#define CONVERSION_HPP

#include <scalable_ccd/broad_phase/sort_and_sweep.hpp>
#include <scalable_ccd/utils/timer.hpp>
#include <scalable_ccd/utils/logger.hpp>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

using namespace scalable_ccd;

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

    std::vector<idx_t> faces(nfaces * 3);
    for (size_t i = 0; i < nfaces; i++) {
      for (int d = 0; d < 3; d++) {
        faces[i * 3 + d] = face_boxes[i].vertex_ids[d];
      }
    }

    std::vector<idx_t> edges(nedges * 2);
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

  void find() {
    sort_axis = lean_choose_axis(nnodes, vaabb);
    lean_sort_along_axis(nfaces, sort_axis, faabb, fidx.data(), scratch.data());
    lean_sort_along_axis(nnodes, sort_axis, vaabb, vidx.data(), scratch.data());
    lean_sort_along_axis(nedges, sort_axis, eaabb, eidx.data(), scratch.data());

    // F2V
    size_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
    ccdptr.resize(max_ccdptr_size);
    lean_count_overlaps<3, 1>(sort_axis, nfaces, faabb, fidx.data(), 3,
                              soafaces, nnodes, vaabb, vidx.data(), 0, nullptr,
                              ccdptr.data());

    // Allocation (expensive, could be expanded dynamically in CCD)
    const size_t fv_nintersections = ccdptr[nfaces];
    foverlap.resize(fv_nintersections);
    voverlap.resize(fv_nintersections);

    lean_collect_overlaps<3, 1>(sort_axis, nfaces, faabb, fidx.data(), 3,
                                soafaces, nnodes, vaabb, vidx.data(), 0,
                                nullptr, ccdptr.data(), foverlap.data(),
                                voverlap.data());

    // E2E
    std::fill(ccdptr.begin(), ccdptr.end(), 0);

    lean_count_self_overlaps<2>(sort_axis, nedges, eaabb, eidx.data(), 2,
                                soaedges, ccdptr.data());

    const size_t ee_n_intersections = ccdptr[nedges];
    e0_overlap.resize(ee_n_intersections);
    e1_overlap.resize(ee_n_intersections);

    lean_collect_self_overlaps<2>(sort_axis, nedges, eaabb, eidx.data(), 2,
                                  soaedges, ccdptr.data(), e0_overlap.data(),
                                  e1_overlap.data());
  }

  void finalize(std::vector<std::pair<int, int>> &vf_overlaps,
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
} SCCD_t;

#endif // CONVERSION_HPP
