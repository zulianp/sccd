#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "sccd_base.hpp"
#include "sccd_config.hpp"
#include "sccd_smesh_CCD.hpp"

#include <algorithm>
#include <filesystem>
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
    SMESH_TRACE_SCOPE("mesh_sccd_cuda.exe");

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

    auto points0 = device_points_as<scalar_t>(t0);
    auto points1 = device_points_as<scalar_t>(t1);

    auto ccd = sccd::CCD<scalar_t>::create(t0, smesh::EXECUTION_SPACE_DEVICE);

    int SCCD_USE_FIND_EARLIEST_IMPACT_TIME = 1;
    SCCD_READ_ENV(SCCD_USE_FIND_EARLIEST_IMPACT_TIME, atoi);

    int err = SCCD_SUCCESS;
    if (SCCD_USE_FIND_EARLIEST_IMPACT_TIME) {
        scalar_t toi = 1;

        err = ccd->find_earliest_impact_time(points0, points1, toi);

        double tock = smesh::time_seconds();
        printf("#faces %ld #edges %ld #nodes %ld, %g [s], toi %g\n",
               t0->block(0)->n_elements(),
               t0->edge_graph()->nnz(),
               t0->n_nodes(),
               (tock - tick),
               (double)toi);

    } else {
        smesh::SharedBuffer<smesh::idx_t> v_overlap;
        smesh::SharedBuffer<smesh::idx_t> f_overlap;
        smesh::SharedBuffer<smesh::idx_t> e0_overlap;
        smesh::SharedBuffer<smesh::idx_t> e1_overlap;
        smesh::SharedBuffer<scalar_t> vf_toi;
        smesh::SharedBuffer<scalar_t> ee_toi;
        scalar_t toi = 1;
        scalar_t toi_vf = 1;
        scalar_t toi_ee = 1;
        ptrdiff_t n_e2e = -1;
        ptrdiff_t n_f2v = -1;

        err = ccd->find_impact_times(points0, points1, v_overlap, f_overlap, vf_toi, e0_overlap, e1_overlap, ee_toi);

        double tock = smesh::time_seconds();

        toi_vf = buffer_min_or(vf_toi, scalar_t(1));
        toi_ee = buffer_min_or(ee_toi, scalar_t(1));
        toi = std::min(toi_vf, toi_ee);
        n_e2e = e0_overlap->size();
        n_f2v = f_overlap->size();

        printf("#faces %ld #edges %ld #nodes %ld, #e2e %ld #f2v %ld, %g [s], toi %g, toi_vf %g, toi_ee %g\n",
               t0->block(0)->n_elements(),
               t0->edge_graph()->nnz(),
               t0->n_nodes(),
               n_e2e,
               n_f2v,
               (tock - tick),
               (double)toi,
               (double)toi_vf,
               (double)toi_ee);

        int SCCD_EXPORT_COLLISIONS = 0;
        SCCD_READ_ENV(SCCD_EXPORT_COLLISIONS, atoi);
        if (SCCD_EXPORT_COLLISIONS && !SCCD_USE_FIND_EARLIEST_IMPACT_TIME) {
            auto folder = smesh::Path("collisions");
            std::filesystem::create_directories("collisions");
            smesh::to_host(v_overlap)->to_file(folder / "v_overlap.int32");
            smesh::to_host(f_overlap)->to_file(folder / "f_overlap.int32");
            smesh::to_host(vf_toi)->to_file(folder / "vf_toi.float64");
            smesh::to_host(e0_overlap)->to_file(folder / "e0_overlap.int32");
            smesh::to_host(e1_overlap)->to_file(folder / "e1_overlap.int32");
            smesh::to_host(ee_toi)->to_file(folder / "ee_toi.float64");
        }
    }

    return 0;
}
