#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include "sccd_base.hpp"
#include "sccd_config.hpp"

#include "sccd_lower_bound_all_to_all.cuh"

#include "sccd_broadphase.cuh"
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
    auto ctx = smesh::initialize_serial(argc, argv);
    SMESH_TRACE_SCOPE("mesh_lower_bound_cuda.exe");

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <mesh_t0> <mesh_t1>\n", argv[0]);
        return 1;
    } else {
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

        int sort_axis = sccd::device::choose_axis(dim, n_nodes, vaabb->data());

        auto vidx = smesh::create_device_buffer<smesh::idx_t>(n_nodes);
        auto fidx = smesh::create_device_buffer<smesh::idx_t>(n_faces);

        auto scratch = smesh::create_device_buffer<scalar_t>(std::max(n_nodes, n_faces));

        {
            SMESH_TRACE_SCOPE("sort");

            sccd::device::sort_along_axis(dim, n_nodes, sort_axis, vaabb->data(), vidx->data(), scratch->data());
            sccd::device::sort_along_axis(dim, n_faces, sort_axis, faabb->data(), fidx->data(), scratch->data());
        }

        auto lb = smesh::create_device_buffer<ptrdiff_t>(n_faces);

        scalar_t* vaabb_max_axis = nullptr;
        auto cm = smesh::create_device_buffer<scalar_t>(n_nodes);
        {
            SMESH_TRACE_SCOPE("cummax");

            cudaMemcpy(
                &vaabb_max_axis, vaabb->data() + dim + sort_axis, sizeof(vaabb_max_axis), cudaMemcpyDeviceToHost);

            sccd::device::cummax(n_nodes, vaabb_max_axis, cm->data());
        }

        auto h_vaabb = smesh::to_host(vaabb);
        auto h_faabb = smesh::to_host(faabb);
        auto h_cm = smesh::to_host(cm);
        auto h_lb = smesh::create_host_buffer<ptrdiff_t>(n_faces);

        int SCCD_LB_REPEAT = 1;
        SCCD_READ_ENV(SCCD_LB_REPEAT, atoi);
        for (int i = 0; i < SCCD_LB_REPEAT; i++) {
            SMESH_TRACE_SCOPE("host");

            sccd::host::lower_bound_all_to_all(
                n_faces, h_faabb->data()[sort_axis], n_nodes, h_cm->data(), h_lb->data());
        }

        for (int i = 0; i < SCCD_LB_REPEAT; i++) {
            SMESH_TRACE_SCOPE("device");

            scalar_t* faabb_min_axis = nullptr;
            cudaMemcpy(&faabb_min_axis, h_faabb->data() + sort_axis, sizeof(faabb_min_axis), cudaMemcpyDeviceToHost);
            sccd::device::lower_bound_all_to_all(n_faces, faabb_min_axis, n_nodes, cm->data(), lb->data());
        }

        auto h_acutal = smesh::to_host(lb);
        for (int i = 0; i < n_faces; i++) {
            if (h_acutal->data()[i] != h_lb->data()[i]) {
                SMESH_ERROR("Lower bound mismatch at face %d: %d != %d\n", i, h_acutal->data()[i], h_lb->data()[i]);
                return 1;
            }
        }
    }

    return 0;
}
