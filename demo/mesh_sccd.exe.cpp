#include "smesh_context.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "broadphase.hpp"
#include "narrowphase.hpp"
#include "sccd_config.hpp"

int main(int argc, char** argv) {
    auto ctx = smesh::initialize(argc, argv);
    SMESH_TRACE_SCOPE("mesh_sccd.exe");

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

        auto elements = t0->elements(0)->data();
        auto p0 = t0->points()->data();
        auto p1 = t1->points()->data();

        // AABB nodes
        const ptrdiff_t n_nodes = t0->n_nodes();
        auto aabb_min_nodes = smesh::create_host_buffer<smesh::geom_t>(dim, n_nodes);
        auto aabb_max_nodes = smesh::create_host_buffer<smesh::geom_t>(dim, n_nodes);
        sccd::compute_aabbs(dim, t0->n_nodes(), p0, p1, aabb_min_nodes->data(), aabb_max_nodes->data());

        // AABB faces
        const ptrdiff_t n_faces = t0->n_elements();
        auto aabb_min_faces = smesh::create_host_buffer<smesh::geom_t>(dim, n_faces);
        auto aabb_max_faces = smesh::create_host_buffer<smesh::geom_t>(dim, n_faces);

        ptrdiff_t element_offset = 0;
        for (auto b : t0->blocks()) {
            int nxe = b->n_nodes_per_element();
            const ptrdiff_t n_elements = b->n_elements();

            auto amin = smesh::view(aabb_min_faces, 0, dim, element_offset, element_offset + n_elements);
            auto amax = smesh::view(aabb_max_faces, 0, dim, element_offset, element_offset + n_elements);

            sccd::compute_aabbs(nxe, n_elements, b->elements()->data(), dim, p0, p1, amin->data(), amax->data());
            element_offset += n_elements;
        }

        // AABB edges
        auto n2n_crs = t0->edge_graph();
        auto row_idx = smesh::create_host_buffer<smesh::idx_t>(n2n_crs->nnz());
        smesh::crs_to_coo(t0->n_nodes(), n2n_crs->rowptr()->data(), row_idx->data());

        const ptrdiff_t n_edges = n2n_crs->nnz();
        auto aabb_min_edges = smesh::create_host_buffer<smesh::geom_t>(dim, n_edges);
        auto aabb_max_edges = smesh::create_host_buffer<smesh::geom_t>(dim, n_edges);

        smesh::idx_t* edges[2] = {row_idx->data(), n2n_crs->colidx()->data()};
        sccd::compute_aabbs(2, row_idx->size(), edges, dim, p0, p1, aabb_min_edges->data(), aabb_max_edges->data());

        // CCD: Broadphase
        auto scratch = smesh::create_host_buffer<smesh::geom_t>(std::max(n_nodes, std::max(n_faces, n_edges)));

        smesh::geom_t* vaabb[6] = {aabb_min_nodes->data()[0],
                                   aabb_min_nodes->data()[1],
                                   aabb_min_nodes->data()[2],
                                   aabb_max_nodes->data()[0],
                                   aabb_max_nodes->data()[1],
                                   aabb_max_nodes->data()[2]};

        smesh::geom_t* faabb[6] = {aabb_min_faces->data()[0],
                                   aabb_min_faces->data()[1],
                                   aabb_min_faces->data()[2],
                                   aabb_max_faces->data()[0],
                                   aabb_max_faces->data()[1],
                                   aabb_max_faces->data()[2]};

        smesh::geom_t* eaabb[6] = {aabb_min_edges->data()[0],
                                   aabb_min_edges->data()[1],
                                   aabb_min_edges->data()[2],
                                   aabb_max_edges->data()[0],
                                   aabb_max_edges->data()[1],
                                   aabb_max_edges->data()[2]};

        int sort_axis = sccd::choose_axis(n_nodes, vaabb);

        auto vidx = smesh::create_host_buffer<smesh::idx_t>(n_nodes);
        for (int i = 0; i < n_nodes; i++) {
            vidx->data()[i] = i;
        }

        auto fidx = smesh::create_host_buffer<smesh::idx_t>(n_faces);
        for (int i = 0; i < n_faces; i++) {
            fidx->data()[i] = i;
        }

        auto eidx = smesh::create_host_buffer<smesh::idx_t>(n_edges);
        for (int i = 0; i < n_edges; i++) {
            eidx->data()[i] = i;
        }

        sccd::sort_along_axis(n_nodes, sort_axis, vaabb, vidx->data(), scratch->data());
        sccd::sort_along_axis(n_faces, sort_axis, faabb, fidx->data(), scratch->data());
        sccd::sort_along_axis(n_edges, sort_axis, eaabb, eidx->data(), scratch->data());

        ptrdiff_t max_ccdptr_size = std::max(n_faces, n_edges) + 1;
        auto ccdptr = smesh::create_host_buffer<ptrdiff_t>(max_ccdptr_size);

        // Results of the broadphase
        std::shared_ptr<smesh::Buffer<smesh::idx_t>> e0_overlap, e1_overlap, f_overlap, v_overlap;

        {
            SMESH_TRACE_SCOPE("Broadphase: E2E");

            sccd::count_self_overlaps<2>(sort_axis, n_edges, eaabb, eidx->data(), 1, edges, ccdptr->data());

            e0_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_edges]);
            e1_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_edges]);

            sccd::collect_self_overlaps<2>(sort_axis,
                                           n_edges,
                                           eaabb,
                                           eidx->data(),
                                           1,
                                           edges,
                                           ccdptr->data(),
                                           e0_overlap->data(),
                                           e1_overlap->data());
        }

        {
            SMESH_TRACE_SCOPE("Broadphase: F2V");

            auto cm = smesh::create_host_buffer<smesh::geom_t>(n_nodes);
            sccd::cummax(n_nodes, vaabb[dim + sort_axis], cm->data());

            sccd::count_overlaps<3, 1, smesh::geom_t, smesh::idx_t>(sort_axis,
                                                                    n_faces,
                                                                    faabb,
                                                                    fidx->data(),
                                                                    1,
                                                                    t0->elements(0)->data(),
                                                                    n_nodes,
                                                                    vaabb,
                                                                    vidx->data(),
                                                                    0,
                                                                    nullptr,
                                                                    ccdptr->data(),
                                                                    cm->data());

            f_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_faces]);
            v_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_faces]);

            sccd::collect_overlaps<3, 1, smesh::geom_t, smesh::idx_t>(sort_axis,
                                                                      n_faces,
                                                                      faabb,
                                                                      fidx->data(),
                                                                      1,
                                                                      t0->elements(0)->data(),
                                                                      n_nodes,
                                                                      vaabb,
                                                                      vidx->data(),
                                                                      0,
                                                                      nullptr,
                                                                      ccdptr->data(),
                                                                      cm->data(),
                                                                      f_overlap->data(),
                                                                      v_overlap->data());
        }

        // Narrow phase
        smesh::geom_t toi = std::numeric_limits<smesh::geom_t>::max();
        smesh::geom_t toi_vf, toi_ee;

        {
            SMESH_TRACE_SCOPE("Narrow phase: F2V");
            toi_vf = sccd::narrow_phase_vf<3, smesh::geom_t>(
                v_overlap->size(), v_overlap->data(), f_overlap->data(), p0, p1, 1, t0->elements(0)->data(), toi);
            toi = toi_vf;
        }

        {
            SMESH_TRACE_SCOPE("Narrow phase: E2E");
            toi_ee = sccd::narrow_phase_ee<smesh::geom_t>(
                e0_overlap->size(), e0_overlap->data(), e1_overlap->data(), p0, p1, 1, edges, toi);
            toi = toi_ee;
        }

        // toi = sccd::min(toi_vf, toi_ee);

        double tock = smesh::time_seconds();
        printf("#faces %ld #edges %ld $nodes %ld, #e2e %ld #f2v %ld, %g [s], toi %g\n",
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
