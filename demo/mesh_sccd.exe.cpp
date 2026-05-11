#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "sccd_config.hpp"
#include "sccd_smesh_CCD.hpp"

#include <algorithm>
#include <limits>

// using scalar_t = float;
using scalar_t = double;

static scalar_t buffer_min_or(const smesh::SharedBuffer<scalar_t>& values, const scalar_t fallback) {
    if (!values || values->size() == 0) {
        return fallback;
    }

    auto host_values = smesh::to_host(values);
    const auto begin = host_values->data();
    return *std::min_element(begin, begin + host_values->size());
}

int main(int argc, char** argv) {
    auto ctx = smesh::initialize(argc, argv);
    SMESH_TRACE_SCOPE("mesh_sccd.exe");

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <mesh_t0> <mesh_t1>\n", argv[0]);
        return 1;
    }

    double tick = smesh::time_seconds();
    auto comm = ctx->communicator();

    auto t0 = smesh::Mesh::create_from_file(comm, smesh::Path(argv[1]));
    auto t1 = smesh::Mesh::create_from_file(comm, smesh::Path(argv[2]));

    if (t0->spatial_dimension() != t1->spatial_dimension()) {
        SMESH_ERROR("Mesh t0 and t1 must have the same spatial dimension\n");
        return 1;
    }

    if (t0->n_elements() != t1->n_elements()) {
        SMESH_ERROR("Mesh t0 and t1 must have the same number of elements\n");
        return 1;
    }

    auto points0 = smesh::astype<scalar_t>(t0->points());
    auto points1 = smesh::astype<scalar_t>(t1->points());

    smesh::SharedBuffer<smesh::idx_t> v_overlap;
    smesh::SharedBuffer<smesh::idx_t> f_overlap;
    smesh::SharedBuffer<smesh::idx_t> e0_overlap;
    smesh::SharedBuffer<smesh::idx_t> e1_overlap;
    smesh::SharedBuffer<scalar_t> vf_toi;
    smesh::SharedBuffer<scalar_t> ee_toi;

    auto ccd = sccd::CCD<scalar_t>::create(t0);
    const int err = ccd->find_impact_times(points0, points1, v_overlap, f_overlap, vf_toi, e0_overlap, e1_overlap, ee_toi);
    if (err != SCCD_SUCCESS) {
        return err;
    }

    const scalar_t toi_vf = buffer_min_or(vf_toi, scalar_t(1));
    const scalar_t toi_ee = buffer_min_or(ee_toi, scalar_t(1));
    const scalar_t toi = std::min(toi_vf, toi_ee);
    const ptrdiff_t n_edges = t0->edge_graph()->nnz();

    double tock = smesh::time_seconds();
    printf("#faces %ld #edges %ld $nodes %ld, #e2e %ld #f2v %ld, %g [s], toi %g, toi_vf %g\n",
           t0->block(0)->n_elements(),
           n_edges,
           t0->n_nodes(),
           e0_overlap->size(),
           f_overlap->size(),
           tock - tick,
           (double)toi,
           (double)toi_vf);

    return 0;
}
