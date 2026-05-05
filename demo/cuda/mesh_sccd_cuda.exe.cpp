#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "broadphase.hpp"
#include "narrowphase.hpp"
#include "sccd_base.hpp"
#include "sccd_config.hpp"

#include "sccd_broadphase.cuh"
#include "sccd_narrowphase.cuh"
#include "sccd_vaabb.cuh"

#include <cuda_runtime_api.h>

#include <algorithm>
#include <limits>
#include <type_traits>

// using scalar_t = smesh::geom_t;
using scalar_t = double;

template <typename T, typename MeshT>
smesh::SharedBuffer<T*> device_points_as(const MeshT& mesh) {
    if constexpr (std::is_same<T, smesh::geom_t>::value) {
        return mesh->device_points_SoA();
    } else {
        auto points = smesh::astype<T>(mesh->points());
        return smesh::to_device(points);
    }
}

int main(int argc, char** argv) {
    auto ctx = smesh::initialize(argc, argv);
    SMESH_TRACE_SCOPE("mesh_sccd_cuda.exe");

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <mesh_t0> <mesh_t1>\n", argv[0]);
        return 1;
    } else {
        double tick = smesh::time_seconds();

        auto comm = ctx->communicator();

        // Read surface meshes
        auto t0 = smesh::Mesh::create_from_file(comm, smesh::Path(argv[1]));
        auto t1 = smesh::Mesh::create_from_file(comm, smesh::Path(argv[2]));

        const int dim = t0->spatial_dimension();

        if (dim != t1->spatial_dimension()) {
            SMESH_ERROR("Mesh t0 and t1 must have the same spatial dimension\n");
            return 1;
        }

        if (t0->n_elements() != t1->n_elements()) {
            SMESH_ERROR("Mesh t0 and t1 must have the same number of elements\n");
            return 1;
        }

        auto faces = t0->block(0)->device_elements_SoA();

        auto points0 = device_points_as<scalar_t>(t0);
        auto points1 = device_points_as<scalar_t>(t1);

        // AABB faces
        const ptrdiff_t n_faces = t0->block(0)->n_elements();

        // AABB nodes
        const ptrdiff_t n_nodes = t0->n_nodes();

        // AABB edges
        auto n2n_crs = t0->edge_graph();
        auto row_idx_temp = smesh::create_host_buffer<smesh::idx_t>(n2n_crs->nnz());
        smesh::crs_to_coo(t0->n_nodes(), n2n_crs->rowptr()->data(), row_idx_temp->data());

        const ptrdiff_t n_edges = n2n_crs->nnz();

        auto row_idx = smesh::to_device(row_idx_temp);
        auto col_idx = smesh::to_device(n2n_crs->colidx());
        auto edges = smesh::create_2d(std::vector{row_idx, col_idx});

        // CCD: Broadphase
        smesh::SharedBuffer<smesh::idx_t> e0_overlap, e1_overlap, f_overlap, v_overlap;
        {
            SMESH_TRACE_SCOPE("Broadphase");

            auto vaabb = smesh::create_device_buffer<scalar_t>(2 * dim, n_nodes);
            sccd::device::compute_aabbs(dim, t0->n_nodes(), points0->data(), points1->data(), vaabb->data());

            auto faabb = smesh::create_device_buffer<scalar_t>(2 * dim, t0->block(0)->n_elements());
            sccd::device::compute_aabbs(t0->block(0)->n_nodes_per_element(),
                                        n_faces,
                                        faces->data(),
                                        dim,
                                        points0->data(),
                                        points1->data(),
                                        faabb->data());

            auto eaabb = smesh::create_device_buffer<scalar_t>(2 * dim, n_edges);
            sccd::device::compute_aabbs(
                2, n_edges, edges->data(), dim, points0->data(), points1->data(), eaabb->data());

            int sort_axis = sccd::device::choose_axis(dim, n_nodes, vaabb->data());

            auto vidx = smesh::create_device_buffer<smesh::idx_t>(n_nodes);
            auto fidx = smesh::create_device_buffer<smesh::idx_t>(n_faces);
            auto eidx = smesh::create_device_buffer<smesh::idx_t>(n_edges);

            auto scratch = smesh::create_device_buffer<scalar_t>(std::max(n_nodes, std::max(n_faces, n_edges)));

            {
                SMESH_TRACE_SCOPE("Sorting AABBs");

                sccd::device::sort_along_axis(dim, n_nodes, sort_axis, vaabb->data(), vidx->data(), scratch->data());
                sccd::device::sort_along_axis(dim, n_faces, sort_axis, faabb->data(), fidx->data(), scratch->data());
                sccd::device::sort_along_axis(dim, n_edges, sort_axis, eaabb->data(), eidx->data(), scratch->data());
            }

            ptrdiff_t max_ccdptr_size = std::max(n_faces, n_edges) + 1;
            auto ccdptr = smesh::create_device_buffer<ptrdiff_t>(max_ccdptr_size);

            // Results of the broadphase
            {
                SMESH_TRACE_SCOPE("Broadphase: E2E");
                sccd::device::count_self_overlaps<2>(
                    sort_axis, n_edges, eaabb->data(), eidx->data(), 1, edges->data(), ccdptr->data());

                ptrdiff_t n_edge_overlaps = 0;
                cudaMemcpy(&n_edge_overlaps, ccdptr->data() + n_edges, sizeof(n_edge_overlaps), cudaMemcpyDeviceToHost);

                e0_overlap = smesh::create_device_buffer<smesh::idx_t>(n_edge_overlaps);
                e1_overlap = smesh::create_device_buffer<smesh::idx_t>(n_edge_overlaps);

                sccd::device::collect_self_overlaps<2>(sort_axis,
                                                       n_edges,
                                                       eaabb->data(),
                                                       eidx->data(),
                                                       1,
                                                       edges->data(),
                                                       ccdptr->data(),
                                                       e0_overlap->data(),
                                                       e1_overlap->data());
            }

            {
                SMESH_TRACE_SCOPE("Broadphase: F2V");
                auto cm = smesh::create_device_buffer<scalar_t>(n_nodes);
                scalar_t* vaabb_max_axis = nullptr;
                cudaMemcpy(
                    &vaabb_max_axis, vaabb->data() + dim + sort_axis, sizeof(vaabb_max_axis), cudaMemcpyDeviceToHost);

                sccd::device::cummax(n_nodes, vaabb_max_axis, cm->data());

                sccd::device::count_overlaps<3, 1>(sort_axis,
                                                   n_faces,
                                                   faabb->data(),
                                                   fidx->data(),
                                                   1,
                                                   faces->data(),
                                                   n_nodes,
                                                   vaabb->data(),
                                                   vidx->data(),
                                                   0,
                                                   (smesh::idx_t**)nullptr,
                                                   ccdptr->data(),
                                                   cm->data());

                ptrdiff_t n_face_overlaps = 0;
                cudaMemcpy(&n_face_overlaps, ccdptr->data() + n_faces, sizeof(n_face_overlaps), cudaMemcpyDeviceToHost);

                f_overlap = smesh::create_device_buffer<smesh::idx_t>(n_face_overlaps);
                v_overlap = smesh::create_device_buffer<smesh::idx_t>(n_face_overlaps);

                sccd::device::collect_overlaps<3, 1>(sort_axis,
                                                     n_faces,
                                                     faabb->data(),
                                                     fidx->data(),
                                                     1,
                                                     faces->data(),
                                                     n_nodes,
                                                     vaabb->data(),
                                                     vidx->data(),
                                                     0,
                                                     (smesh::idx_t**)nullptr,
                                                     ccdptr->data(),
                                                     cm->data(),
                                                     f_overlap->data(),
                                                     v_overlap->data());
            }
        }

        // Narrow phase
        scalar_t toi = 1;
        scalar_t toi_vf = toi;
        scalar_t toi_ee = toi;

        int SCCD_NP_TOI_STRIDE = 0;
        SCCD_READ_ENV(SCCD_NP_TOI_STRIDE, atoi);
        const auto toi_storage_size = [SCCD_NP_TOI_STRIDE](const size_t noverlaps) {
            return SCCD_NP_TOI_STRIDE == 0 ? size_t(1) : noverlaps;
        };

        const auto scalar_device_to_host = [](scalar_t* const d_toi) {
            scalar_t h_toi = std::numeric_limits<scalar_t>::max();
            cudaMemcpy(&h_toi, d_toi, sizeof(h_toi), cudaMemcpyDeviceToHost);
            return h_toi;
        };

        {
            SMESH_TRACE_SCOPE("Narrow phase");
            smesh::SharedBuffer<scalar_t> vf_toi;
            smesh::SharedBuffer<scalar_t> ee_toi;

            {
                SMESH_TRACE_SCOPE("Narrow phase: F2V");
                vf_toi = smesh::create_device_buffer<scalar_t>(toi_storage_size(v_overlap->size()));
                sccd::device::narrow_phase_vf<3>(v_overlap->size(),
                                                 v_overlap->data(),
                                                 f_overlap->data(),
                                                 points0->data(),
                                                 points1->data(),
                                                 1,
                                                 faces->data(),
                                                 toi,
                                                 vf_toi->data(),
                                                 SCCD_NP_TOI_STRIDE);

                if (SCCD_NP_TOI_STRIDE == 0 && v_overlap->size() != 0) {
                    toi_vf = scalar_device_to_host(vf_toi->data());
                    toi = toi_vf;
                }
            }

            {
                SMESH_TRACE_SCOPE("Narrow phase: E2E");
                ee_toi = smesh::create_device_buffer<scalar_t>(toi_storage_size(e0_overlap->size()));
                sccd::device::narrow_phase_ee(e0_overlap->size(),
                                              e0_overlap->data(),
                                              e1_overlap->data(),
                                              points0->data(),
                                              points1->data(),
                                              1,
                                              edges->data(),
                                              toi,
                                              ee_toi->data(),
                                              SCCD_NP_TOI_STRIDE);
                if (SCCD_NP_TOI_STRIDE == 0 && e0_overlap->size() != 0) {
                    toi_ee = scalar_device_to_host(ee_toi->data());
                    toi = toi_ee;
                }
            }

            if (SCCD_NP_TOI_STRIDE == 1) {
                const auto device_minmax_to_host =
                    [](scalar_t* const d_values, const size_t nvalues, scalar_t* const h_min, scalar_t* const h_max) {
                        if (nvalues == 0) {
                            *h_min = std::numeric_limits<scalar_t>::max();
                            *h_max = std::numeric_limits<scalar_t>::lowest();
                            return;
                        }

                        sccd::device::minmax(d_values, nvalues, h_min, h_max);
                    };

                scalar_t vf_toi_max = std::numeric_limits<scalar_t>::lowest();
                scalar_t ee_toi_max = std::numeric_limits<scalar_t>::lowest();
                device_minmax_to_host(vf_toi->data(), v_overlap->size(), &toi_vf, &vf_toi_max);
                device_minmax_to_host(ee_toi->data(), e0_overlap->size(), &toi_ee, &ee_toi_max);
                toi = std::min(toi_vf, toi_ee);

                printf("vt_toi_min: %g, ee_toi_min: %g\n", toi_vf, toi_ee);
                printf("vf_toi_max: %g, ee_toi_max: %g\n", vf_toi_max, ee_toi_max);
            }
        }

        double tock = smesh::time_seconds();
        printf("#faces %ld #edges %ld #nodes %ld, #e2e %ld #f2v %ld, %g [s], toi %g\n",
               n_faces,
               n_edges,
               n_nodes,
               e0_overlap->size(),
               f_overlap->size(),
               tock - tick,
               (double)toi);
    }

    return 0;
}
