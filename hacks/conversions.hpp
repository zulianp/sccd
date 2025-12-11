#ifndef CONVERSION_HPP
#define CONVERSION_HPP

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <scalable_ccd/broad_phase/sort_and_sweep.hpp>
#include <scalable_ccd/utils/logger.hpp>
#include <scalable_ccd/utils/timer.hpp>

using namespace scalable_ccd;

#include "broadphase_lb.hpp"
#include "cell_broadphase.hpp"
#include "narrowphase.hpp"

/// Geometry scalar type used for coordinates and AABB values.
using geom_t = float;
/// Integer type used for element/vertex indices.
using idx_t = int;
/// Integer type used for counts and prefix sums.
using count_t = int;

namespace sccd {

    // (Ugly) Conversion between Scalable-CCD and sccd and back

    typedef struct SCCD {
        int sort_axis;
        ptrdiff_t nfaces;
        ptrdiff_t nedges;
        ptrdiff_t nnodes;

        geom_t* faabb[6];
        geom_t* eaabb[6];
        geom_t* vaabb[6];

        idx_t* soafaces[3];
        idx_t* soaedges[2];

        // Mesh data
        std::vector<geom_t> x0;
        std::vector<geom_t> y0;
        std::vector<geom_t> z0;

        std::vector<geom_t> x1;
        std::vector<geom_t> y1;
        std::vector<geom_t> z1;

        // Time of Impact
        std::vector<geom_t> vf_toi;
        std::vector<geom_t> ee_toi;

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
        std::vector<ptrdiff_t> ccdptr;

        // Scratch space for sorting
        std::vector<geom_t> scratch;

        // Workspace allocations
        std::vector<std::vector<geom_t>> face_aabb;
        std::vector<std::vector<geom_t>> edge_aabb;
        std::vector<std::vector<geom_t>> vertex_aabb;

        std::vector<idx_t> faces;
        std::vector<idx_t> edges;
        bool verbose{true};

        void init(const std::vector<AABB>& vertex_boxes,
                  const std::vector<AABB>& edge_boxes,
                  const std::vector<AABB>& face_boxes) {
            nfaces = face_boxes.size();
            nedges = edge_boxes.size();
            nnodes = vertex_boxes.size();

            fidx.resize(nfaces);
            vidx.resize(nnodes);
            eidx.resize(nedges);

            scratch.resize(std::max(nfaces, std::max(nedges, nnodes)));

            ptrdiff_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
            ccdptr.resize(max_ccdptr_size);

            face_aabb.resize(6);
            for (auto& v : face_aabb) {
                v.resize(nfaces);
            }

            edge_aabb.resize(6);
            for (auto& v : edge_aabb) {
                v.resize(nedges);
            }

            vertex_aabb.resize(6);
            for (auto& v : vertex_aabb) {
                v.resize(nnodes);
            }

            // C++ to C conversion
            for (int d = 0; d < 6; d++) {
                faabb[d] = face_aabb[d].data();
                eaabb[d] = edge_aabb[d].data();
                vaabb[d] = vertex_aabb[d].data();
            }

            for (ptrdiff_t i = 0; i < nfaces; i++) {
                for (int d = 0; d < 3; d++) {
                    faabb[d][i] = face_boxes[i].min[d];
                    faabb[3 + d][i] = face_boxes[i].max[d];
                }
            }

            for (ptrdiff_t i = 0; i < nnodes; i++) {
                for (int d = 0; d < 3; d++) {
                    vaabb[d][i] = vertex_boxes[i].min[d];
                    vaabb[3 + d][i] = vertex_boxes[i].max[d];
                }
            }

            for (ptrdiff_t i = 0; i < nedges; i++) {
                for (int d = 0; d < 3; d++) {
                    eaabb[d][i] = edge_boxes[i].min[d];
                    eaabb[3 + d][i] = edge_boxes[i].max[d];
                }
            }

            faces.resize(nfaces * 3);
            for (ptrdiff_t i = 0; i < nfaces; i++) {
                for (int d = 0; d < 3; d++) {
                    faces[i * 3 + d] = face_boxes[i].vertex_ids[d];
                }
            }

            edges.resize(nedges * 2);
            for (ptrdiff_t i = 0; i < nedges; i++) {
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

        //     ptrdiff_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
        //     ccdptr.resize(max_ccdptr_size);

        //     int axes[3];
        //     largest_variance_axes_sort(nnodes, vaabb, axes);

        //     sort_axis = axes[0];
        //     sort_along_axis(nfaces, sort_axis, faabb, fidx.data(),
        //     scratch.data()); sort_along_axis(nnodes, sort_axis, vaabb,
        //     vidx.data(), scratch.data()); sort_along_axis(nedges, sort_axis,
        //     eaabb, eidx.data(), scratch.data());

        //     timer.stop();
        //     if (verbose)
        //         printf(
        //             "SCCD, Sorting: %g [ms]\n",
        //             timer.getElapsedTimeInMilliSec());
        //     timer.start();

        //     int cell_list_axis = sort_axis;
        //     ptrdiff_t ncells = 1024; // Max amount
        //     geom_t cell_min;
        //     geom_t cell_size;
        //     cell_starts_setup(
        //         nnodes, vaabb[cell_list_axis], vaabb[cell_list_axis + 3], ncells,
        //         &cell_min, &cell_size);
        //     std::vector<ptrdiff_t> starts(ncells);
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
        //         sort_axis, nfaces, faabb, fidx.data(), 3, soafaces, nnodes,
        //         vaabb, vidx.data(), 0, nullptr, ncells, cell_min, cell_size,
        //         starts.data(), ccdptr.data());

        //     timer.stop();
        //     if (verbose)
        //         printf(
        //             "SCCD count cell starts(%lu), F2V: %g [ms]\n",
        //             ccdptr[nfaces], timer.getElapsedTimeInMilliSec());
        // }

        void find_with_cell_list() {
            Timer timer;

            timer.start();

            ptrdiff_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
            ccdptr.resize(max_ccdptr_size);

            int axes[3];
            largest_variance_axes_sort(nnodes, vaabb, axes);

            sort_axis = axes[0];
            sort_along_axis(nfaces, sort_axis, faabb, fidx.data(), scratch.data());
            sort_along_axis(nnodes, sort_axis, vaabb, vidx.data(), scratch.data());
            sort_along_axis(nedges, sort_axis, eaabb, eidx.data(), scratch.data());

            int cell_list_axis = axes[1];
            ptrdiff_t ncells = 2048;  // Max amount
            geom_t cell_min;
            geom_t cell_size;
            cell_setup(nnodes, vaabb[cell_list_axis], vaabb[cell_list_axis + 3], &ncells, &cell_min, &cell_size);
            std::vector<idx_t> cellptr(ncells + 1), bookkeeping(ncells);
            cell_count(ncells, cell_min, cell_size, nnodes, vaabb[cell_list_axis], cellptr.data());

            std::vector<idx_t> cellidx(cellptr[ncells]);
            cell_populate(ncells,
                          cell_min,
                          cell_size,
                          nnodes,
                          vaabb[cell_list_axis],
                          cellptr.data(),
                          cellidx.data(),
                          bookkeeping.data());

            timer.stop();
            if (verbose) printf("SCCD(CELL), Sorting: %g [ms]\n", timer.getElapsedTimeInMilliSec());
            timer.start();

            // E2E
            std::fill(ccdptr.begin(), ccdptr.end(), 0);

            count_self_overlaps<2>(sort_axis, nedges, eaabb, eidx.data(), 2, soaedges, ccdptr.data());

            const ptrdiff_t ee_n_intersections = ccdptr[nedges];
            e0_overlap.resize(ee_n_intersections);
            e1_overlap.resize(ee_n_intersections);

            collect_self_overlaps<2>(sort_axis,
                                     nedges,
                                     eaabb,
                                     eidx.data(),
                                     2,
                                     soaedges,
                                     ccdptr.data(),
                                     e0_overlap.data(),
                                     e1_overlap.data());

            timer.stop();
            if (verbose) printf("SCCD, E2E: %g [ms]\n", timer.getElapsedTimeInMilliSec());
            timer.start();

            Timer co;
            co.start();
            // F2V
            cell_count_overlaps<3, 1, geom_t, idx_t>(sort_axis,
                                                     nfaces,
                                                     faabb,
                                                     fidx.data(),
                                                     3,
                                                     soafaces,
                                                     nnodes,
                                                     vaabb,
                                                     vidx.data(),
                                                     0,
                                                     nullptr,
                                                     cell_list_axis,
                                                     ncells,
                                                     cell_min,
                                                     cell_size,
                                                     cellptr.data(),
                                                     cellidx.data(),
                                                     ccdptr.data());

            co.stop();
            printf("SCCD(CELL): F2V Count %g [ms]\n", co.getElapsedTimeInMilliSec());

            // Allocation (expensive, could be expanded dynamically in CCD)
            const ptrdiff_t fv_nintersections = ccdptr[nfaces];
            foverlap.resize(fv_nintersections);
            voverlap.resize(fv_nintersections);

            cell_collect_overlaps<3, 1, geom_t, idx_t>(sort_axis,
                                                       nfaces,
                                                       faabb,
                                                       fidx.data(),
                                                       3,
                                                       soafaces,
                                                       nnodes,
                                                       vaabb,
                                                       vidx.data(),
                                                       0,
                                                       nullptr,
                                                       cell_list_axis,
                                                       ncells,
                                                       cell_min,
                                                       cell_size,
                                                       cellptr.data(),
                                                       cellidx.data(),
                                                       ccdptr.data(),
                                                       foverlap.data(),
                                                       voverlap.data());

            timer.stop();
            if (verbose) printf("SCCD(CELL), F2V: %g [ms]\n", timer.getElapsedTimeInMilliSec());
        }

        void broad_phase() {
            Timer timer;

            timer.start();

            ptrdiff_t max_ccdptr_size = std::max(nfaces, nedges) + 1;
            ccdptr.resize(max_ccdptr_size);

            sort_axis = choose_axis(nnodes, vaabb);
            sort_along_axis(nfaces, sort_axis, faabb, fidx.data(), scratch.data());
            sort_along_axis(nnodes, sort_axis, vaabb, vidx.data(), scratch.data());
            sort_along_axis(nedges, sort_axis, eaabb, eidx.data(), scratch.data());

            timer.stop();
            if (verbose) printf("SCCD, Sorting: %g [ms]\n", timer.getElapsedTimeInMilliSec());
            timer.start();

            // E2E
            std::fill(ccdptr.begin(), ccdptr.end(), 0);

            count_self_overlaps<2>(sort_axis, nedges, eaabb, eidx.data(), 2, soaedges, ccdptr.data());

            const ptrdiff_t ee_n_intersections = ccdptr[nedges];
            e0_overlap.resize(ee_n_intersections);
            e1_overlap.resize(ee_n_intersections);

            collect_self_overlaps<2>(sort_axis,
                                     nedges,
                                     eaabb,
                                     eidx.data(),
                                     2,
                                     soaedges,
                                     ccdptr.data(),
                                     e0_overlap.data(),
                                     e1_overlap.data());

            timer.stop();
            if (verbose) printf("SCCD, E2E: %g [ms]\n", timer.getElapsedTimeInMilliSec());
            timer.start();

            Timer co;
            co.start();

#define USE_LB
#ifndef USE_LB
            // F2V
            count_overlaps<3, 1, geom_t, idx_t>(sort_axis,
                                                nfaces,
                                                faabb,
                                                fidx.data(),
                                                3,
                                                soafaces,
                                                nnodes,
                                                vaabb,
                                                vidx.data(),
                                                0,
                                                nullptr,
                                                ccdptr.data());

#else

            std::vector<geom_t> lb = vertex_aabb[3+sort_axis];
            sccd::parallel_cum_max_br(lb.data(), lb.data() + lb.size());

            count_overlaps_lb<3, 1, geom_t, idx_t>(sort_axis,
                                                   nfaces,
                                                   faabb,
                                                   fidx.data(),
                                                   3,
                                                   soafaces,
                                                   nnodes,
                                                   vaabb,
                                                   vidx.data(),
                                                   0,
                                                   nullptr,
                                                   lb.data(),
                                                   ccdptr.data());
#endif

            co.stop();
            printf("SCCD: F2V Count %g [ms]\n", co.getElapsedTimeInMilliSec());

            // Allocation (expensive, could be expanded dynamically in CCD)
            const ptrdiff_t fv_nintersections = ccdptr[nfaces];
            foverlap.resize(fv_nintersections);
            voverlap.resize(fv_nintersections);

#ifndef USE_LB
            collect_overlaps<3, 1, geom_t, idx_t>(sort_axis,
                                                  nfaces,
                                                  faabb,
                                                  fidx.data(),
                                                  3,
                                                  soafaces,
                                                  nnodes,
                                                  vaabb,
                                                  vidx.data(),
                                                  0,
                                                  nullptr,
                                                  ccdptr.data(),
                                                  foverlap.data(),
                                                  voverlap.data());
#else

            collect_overlaps_lb<3, 1, geom_t, idx_t>(sort_axis,
                                                     nfaces,
                                                     faabb,
                                                     fidx.data(),
                                                     3,
                                                     soafaces,
                                                     nnodes,
                                                     vaabb,
                                                     vidx.data(),
                                                     0,
                                                     nullptr,
                                                     lb.data(),
                                                     ccdptr.data(),
                                                     foverlap.data(),
                                                     voverlap.data());

#endif

            timer.stop();
            if (verbose) printf("SCCD, F2V: %g [ms]\n", timer.getElapsedTimeInMilliSec());
        }

        void export_broadphase_results(std::vector<std::pair<int, int>>& vf_overlaps,
                                       std::vector<std::pair<int, int>>& ee_overlaps) {
            const ptrdiff_t fv_nintersections = foverlap.size();
            vf_overlaps.resize(fv_nintersections);
            for (ptrdiff_t i = 0; i < fv_nintersections; i++) {
                vf_overlaps[i].first = voverlap[i];
                vf_overlaps[i].second = foverlap[i];
            }

            const ptrdiff_t ee_n_intersections = e0_overlap.size();
            ee_overlaps.resize(ee_n_intersections);
            for (ptrdiff_t i = 0; i < ee_n_intersections; i++) {
                ee_overlaps[i].first = e0_overlap[i];
                ee_overlaps[i].second = e1_overlap[i];
            }
        }

        void init_point_data(const Eigen::MatrixXd& v0, const Eigen::MatrixXd& v1) {
            x0.resize(v0.rows());
            y0.resize(v0.rows());
            z0.resize(v0.rows());

            x1.resize(v1.rows());
            y1.resize(v1.rows());
            z1.resize(v1.rows());

            tbb::parallel_for(tbb::blocked_range<long>(0, v0.rows()), [&](const tbb::blocked_range<long>& r) {
                for (long i = r.begin(); i < r.end(); i++) {
                    auto row0 = v0.row(i);
                    auto row1 = v1.row(i);

                    x0[i] = row0[0];
                    y0[i] = row0[1];
                    z0[i] = row0[2];

                    x1[i] = row1[0];
                    y1[i] = row1[1];
                    z1[i] = row1[2];
                }
            });
        }

        geom_t narrow_phase() {
            geom_t* v0[3] = {x0.data(), y0.data(), z0.data()};
            geom_t* v1[3] = {x1.data(), y1.data(), z1.data()};

            Timer timer;
            timer.start();

            vf_toi.resize(voverlap.size());

            geom_t toi_vf = sccd::narrow_phase_vf<3, geom_t>(voverlap.size(),
                                                             voverlap.data(),
                                                             foverlap.data(),
                                                             // Geometric data
                                                             v0,
                                                             v1,
                                                             3,
                                                             soafaces,
                                                             // Output
                                                             vf_toi.data());

            timer.stop();
            printf("NP VF(%lu): %g [ms]\n", vf_toi.size(), timer.getElapsedTimeInMilliSec());
            timer.start();

            ee_toi.resize(e0_overlap.size());

            geom_t toi_ee = sccd::narrow_phase_ee<geom_t>(e0_overlap.size(),
                                                          e0_overlap.data(),
                                                          e1_overlap.data(),
                                                          // Geometric data
                                                          v0,
                                                          v1,
                                                          2,
                                                          soaedges,
                                                          // Output
                                                          ee_toi.data());

            timer.stop();
            printf("NP EE(%lu): %g [ms]\n", ee_toi.size(), timer.getElapsedTimeInMilliSec());

            return sccd::min(toi_ee, toi_vf);
        }

        void export_narrowphase_results(std::vector<std::tuple<int, int, Scalar>>& vf_collisions,
                                        std::vector<std::tuple<int, int, Scalar>>& ee_collisions) {
            const ptrdiff_t vf_size = vf_toi.size();
            const ptrdiff_t ee_size = ee_toi.size();

            vf_collisions.reserve(vf_size);
            ee_collisions.reserve(ee_size);

            for (ptrdiff_t i = 0; i < vf_size; i++) {
                if (vf_toi[i] < 1.1) {
                    vf_collisions.push_back({voverlap[i], foverlap[i], vf_toi[i]});
                }
            }

            for (ptrdiff_t i = 0; i < ee_size; i++) {
                if (ee_toi[i] < 1.1) {
                    ee_collisions.push_back({e0_overlap[i], e1_overlap[i], ee_toi[i]});
                }
            }
        }

    } SCCD_t;
}  // namespace sccd

#endif  // CONVERSION_HPP
