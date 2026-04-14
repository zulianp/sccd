#include "smesh_context.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "broadphase.hpp"
#include "narrowphase.hpp"
#include "sccd_config.hpp"

#include "sccd_broadphase.cuh"
#include "sccd_vaabb.cuh"

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
        auto points0 = t0->device_points_SoA();
        auto points1 = t1->device_points_SoA();

        auto faabb = smesh::create_device_buffer<smesh::geom_t>(2 * dim, t0->block(0)->n_elements());

        // AABB faces
        const ptrdiff_t n_faces = t0->block(0)->n_elements();
        sccd::device::compute_aabbs<smesh::idx_t, smesh::geom_t, smesh::geom_t>(t0->block(0)->n_nodes_per_element(),
                                                                                n_faces,
                                                                                faces->data(),
                                                                                dim,
                                                                                points0->data(),
                                                                                points1->data(),
                                                                                faabb->data());

        // AABB nodes
        const ptrdiff_t n_nodes = t0->n_nodes();
        auto vaabb = smesh::create_device_buffer<smesh::geom_t>(2 * dim, n_nodes);
        sccd::device::compute_aabbs<smesh::geom_t, smesh::geom_t>(
            dim, t0->n_nodes(), points0->data(), points1->data(), vaabb->data());

        // AABB edges
        auto n2n_crs = t0->edge_graph();
        auto row_idx_temp = smesh::create_host_buffer<smesh::idx_t>(n2n_crs->nnz());
        smesh::crs_to_coo(t0->n_nodes(), n2n_crs->rowptr()->data(), row_idx_temp->data());

        const ptrdiff_t n_edges = n2n_crs->nnz();

        auto row_idx = smesh::to_device(row_idx_temp);
        auto col_idx = smesh::to_device(n2n_crs->colidx());
        auto edges = smesh::create_2d(std::vector{row_idx, col_idx});

        auto eaabb = smesh::create_device_buffer<smesh::geom_t>(2 * dim, n_edges);
        sccd::device::compute_aabbs<smesh::idx_t, smesh::geom_t, smesh::geom_t>(
            2, n_edges, edges->data(), dim, points0->data(), points1->data(), eaabb->data());

        // CCD: Broadphase
        int sort_axis = sccd::device::choose_axis(dim, n_nodes, vaabb->data());

        auto vidx = smesh::create_device_buffer<smesh::idx_t>(n_nodes);
        auto fidx = smesh::create_device_buffer<smesh::idx_t>(n_faces);
        auto eidx = smesh::create_device_buffer<smesh::idx_t>(n_edges);

        auto scratch = smesh::create_device_buffer<smesh::geom_t>(std::max(n_nodes, std::max(n_faces, n_edges)));

        sccd::device::sort_along_axis(dim, n_nodes, sort_axis, vaabb->data(), vidx->data(), scratch->data());
        sccd::device::sort_along_axis(dim, n_faces, sort_axis, faabb->data(), fidx->data(), scratch->data());
        sccd::device::sort_along_axis(dim, n_edges, sort_axis, eaabb->data(), eidx->data(), scratch->data());

        //     ptrdiff_t max_ccdptr_size = std::max(n_faces, n_edges) + 1;
        //     auto ccdptr = smesh::create_host_buffer<ptrdiff_t>(max_ccdptr_size);

        //     // Results of the broadphase
        //     std::shared_ptr<smesh::Buffer<smesh::idx_t>> e0_overlap, e1_overlap, f_overlap, v_overlap;

        //     {
        //         SMESH_TRACE_SCOPE("Broadphase: E2E");
        //         sccd::count_self_overlaps<2>(sort_axis, n_edges, eaabb, eidx->data(), 1, edges, ccdptr->data());

        //         e0_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_edges]);
        //         e1_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_edges]);

        //         sccd::collect_self_overlaps<2>(sort_axis,
        //                                        n_edges,
        //                                        eaabb,
        //                                        eidx->data(),
        //                                        1,
        //                                        edges,
        //                                        ccdptr->data(),
        //                                        e0_overlap->data(),
        //                                        e1_overlap->data());
        //     }

        //     {
        //         SMESH_TRACE_SCOPE("Broadphase: F2V");
        //         sccd::count_overlaps<3, 1, smesh::geom_t, smesh::idx_t>(sort_axis,
        //                                                                 n_faces,
        //                                                                 faabb,
        //                                                                 fidx->data(),
        //                                                                 1,
        //                                                                 t0->elements(0)->data(),
        //                                                                 n_nodes,
        //                                                                 vaabb,
        //                                                                 vidx->data(),
        //                                                                 0,
        //                                                                 nullptr,
        //                                                                 ccdptr->data());

        //         f_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_faces]);
        //         v_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_faces]);

        //         sccd::collect_overlaps<3, 1, smesh::geom_t, smesh::idx_t>(sort_axis,
        //                                                                   n_faces,
        //                                                                   faabb,
        //                                                                   fidx->data(),
        //                                                                   1,
        //                                                                   t0->elements(0)->data(),
        //                                                                   n_nodes,
        //                                                                   vaabb,
        //                                                                   vidx->data(),
        //                                                                   0,
        //                                                                   nullptr,
        //                                                                   ccdptr->data(),
        //                                                                   f_overlap->data(),
        //                                                                   v_overlap->data());
        //     }

        //     // Narrow phase
        //     smesh::geom_t toi = std::numeric_limits<smesh::geom_t>::max();
        //     smesh::geom_t toi_vf, toi_ee;

        //     {
        //         SMESH_TRACE_SCOPE("Narrow phase: F2V");
        //         toi_vf = sccd::narrow_phase_vf<3, smesh::geom_t>(
        //             v_overlap->size(), v_overlap->data(), f_overlap->data(), p0, p1, 1, t0->elements(0)->data(),
        //             toi);
        //         toi = toi_vf;
        //     }

        //     {
        //         SMESH_TRACE_SCOPE("Narrow phase: E2E");
        //         toi_ee = sccd::narrow_phase_ee<smesh::geom_t>(
        //             e0_overlap->size(), e0_overlap->data(), e1_overlap->data(), p0, p1, 1, edges, toi);
        //         toi = toi_ee;
        //     }

        //     // toi = sccd::min(toi_vf, toi_ee);

        //     double tock = smesh::time_seconds();
        //     printf("#faces %ld #edges %ld $nodes %ld, #e2e %ld #f2v %ld, %g [s], toi %g\n",
        //            n_faces,
        //            n_edges,
        //            n_nodes,
        //            e0_overlap->size(),
        //            f_overlap->size(),
        //            tock - tick,
        //            (double)toi);
    }

    return 0;
}
