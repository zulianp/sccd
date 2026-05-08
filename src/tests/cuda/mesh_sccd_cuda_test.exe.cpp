#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_test.hpp"
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
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

using scalar_t = smesh::geom_t;
using pair_t = std::pair<smesh::idx_t, smesh::idx_t>;
using toi_pair_t = std::tuple<smesh::idx_t, smesh::idx_t, scalar_t>;

static std::shared_ptr<smesh::Mesh> g_t0;
static std::shared_ptr<smesh::Mesh> g_t1;

struct CandidateData {
    std::vector<pair_t> vf_pairs;
    std::vector<pair_t> ee_pairs;
    std::vector<scalar_t> vf_toi;
    std::vector<scalar_t> ee_toi;
};

static std::vector<toi_pair_t> zip_and_sort(const std::vector<pair_t>& pairs, const std::vector<scalar_t>& toi) {
    std::vector<toi_pair_t> out;
    out.reserve(pairs.size());
    for (size_t i = 0; i < pairs.size(); ++i) {
        out.emplace_back(pairs[i].first, pairs[i].second, toi[i]);
    }
    std::sort(out.begin(), out.end(), [](const toi_pair_t& a, const toi_pair_t& b) {
        if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
        return std::get<1>(a) < std::get<1>(b);
    });
    return out;
}

static CandidateData run_cpu_reference(const std::shared_ptr<smesh::Mesh>& t0, const std::shared_ptr<smesh::Mesh>& t1) {
    CandidateData out;
    const int dim = t0->spatial_dimension();

    auto p0 = t0->points()->data();
    auto p1 = t1->points()->data();
    auto elements = t0->elements(0)->data();

    const ptrdiff_t n_nodes = t0->n_nodes();
    const ptrdiff_t n_faces = t0->n_elements();
    auto n2n_crs = t0->edge_graph();
    auto row_idx = smesh::create_host_buffer<smesh::idx_t>(n2n_crs->nnz());
    smesh::crs_to_coo(t0->n_nodes(), n2n_crs->rowptr()->data(), row_idx->data());
    const ptrdiff_t n_edges = n2n_crs->nnz();
    smesh::idx_t* edges[2] = {row_idx->data(), n2n_crs->colidx()->data()};

    auto aabb_min_nodes = smesh::create_host_buffer<scalar_t>(dim, n_nodes);
    auto aabb_max_nodes = smesh::create_host_buffer<scalar_t>(dim, n_nodes);
    auto aabb_min_faces = smesh::create_host_buffer<scalar_t>(dim, n_faces);
    auto aabb_max_faces = smesh::create_host_buffer<scalar_t>(dim, n_faces);
    auto aabb_min_edges = smesh::create_host_buffer<scalar_t>(dim, n_edges);
    auto aabb_max_edges = smesh::create_host_buffer<scalar_t>(dim, n_edges);

    sccd::compute_aabbs(dim, n_nodes, p0, p1, aabb_min_nodes->data(), aabb_max_nodes->data());
    sccd::compute_aabbs(3, n_faces, elements, dim, p0, p1, aabb_min_faces->data(), aabb_max_faces->data());
    sccd::compute_aabbs(2, n_edges, edges, dim, p0, p1, aabb_min_edges->data(), aabb_max_edges->data());

    scalar_t* vaabb[6] = {aabb_min_nodes->data()[0],
                          aabb_min_nodes->data()[1],
                          aabb_min_nodes->data()[2],
                          aabb_max_nodes->data()[0],
                          aabb_max_nodes->data()[1],
                          aabb_max_nodes->data()[2]};
    scalar_t* faabb[6] = {aabb_min_faces->data()[0],
                          aabb_min_faces->data()[1],
                          aabb_min_faces->data()[2],
                          aabb_max_faces->data()[0],
                          aabb_max_faces->data()[1],
                          aabb_max_faces->data()[2]};
    scalar_t* eaabb[6] = {aabb_min_edges->data()[0],
                          aabb_min_edges->data()[1],
                          aabb_min_edges->data()[2],
                          aabb_max_edges->data()[0],
                          aabb_max_edges->data()[1],
                          aabb_max_edges->data()[2]};

    int sort_axis = sccd::choose_axis(n_nodes, vaabb);
    auto scratch = smesh::create_host_buffer<scalar_t>(std::max(n_nodes, std::max(n_faces, n_edges)));
    auto vidx = smesh::create_host_buffer<smesh::idx_t>(n_nodes);
    auto fidx = smesh::create_host_buffer<smesh::idx_t>(n_faces);
    auto eidx = smesh::create_host_buffer<smesh::idx_t>(n_edges);
    for (ptrdiff_t i = 0; i < n_nodes; ++i) vidx->data()[i] = (smesh::idx_t)i;
    for (ptrdiff_t i = 0; i < n_faces; ++i) fidx->data()[i] = (smesh::idx_t)i;
    for (ptrdiff_t i = 0; i < n_edges; ++i) eidx->data()[i] = (smesh::idx_t)i;
    sccd::sort_along_axis(n_nodes, sort_axis, vaabb, vidx->data(), scratch->data());
    sccd::sort_along_axis(n_faces, sort_axis, faabb, fidx->data(), scratch->data());
    sccd::sort_along_axis(n_edges, sort_axis, eaabb, eidx->data(), scratch->data());

    auto ccdptr = smesh::create_host_buffer<ptrdiff_t>(std::max(n_faces, n_edges) + 1);

    sccd::count_self_overlaps<2>(sort_axis, n_edges, eaabb, eidx->data(), 1, edges, ccdptr->data());
    auto e0_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_edges]);
    auto e1_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_edges]);
    sccd::collect_self_overlaps<2>(
        sort_axis, n_edges, eaabb, eidx->data(), 1, edges, ccdptr->data(), e0_overlap->data(), e1_overlap->data());

    auto cm = smesh::create_host_buffer<scalar_t>(n_nodes);
    sccd::cummax(n_nodes, vaabb[dim + sort_axis], cm->data());
    sccd::count_overlaps<3, 1, scalar_t, smesh::idx_t>(sort_axis,
                                                       n_faces,
                                                       faabb,
                                                       fidx->data(),
                                                       1,
                                                       elements,
                                                       n_nodes,
                                                       vaabb,
                                                       vidx->data(),
                                                       0,
                                                       nullptr,
                                                       ccdptr->data(),
                                                       cm->data());
    auto f_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_faces]);
    auto v_overlap = smesh::create_host_buffer<smesh::idx_t>(ccdptr->data()[n_faces]);
    sccd::collect_overlaps<3, 1, scalar_t, smesh::idx_t>(sort_axis,
                                                         n_faces,
                                                         faabb,
                                                         fidx->data(),
                                                         1,
                                                         elements,
                                                         n_nodes,
                                                         vaabb,
                                                         vidx->data(),
                                                         0,
                                                         nullptr,
                                                         ccdptr->data(),
                                                         cm->data(),
                                                         f_overlap->data(),
                                                         v_overlap->data());

    out.vf_pairs.resize(v_overlap->size());
    out.ee_pairs.resize(e0_overlap->size());
    for (size_t i = 0; i < out.vf_pairs.size(); ++i) out.vf_pairs[i] = {v_overlap->data()[i], f_overlap->data()[i]};
    for (size_t i = 0; i < out.ee_pairs.size(); ++i) out.ee_pairs[i] = {e0_overlap->data()[i], e1_overlap->data()[i]};

    out.vf_toi.resize(v_overlap->size());
    out.ee_toi.resize(e0_overlap->size());
    constexpr scalar_t max_toi = scalar_t(100000);
    sccd::narrow_phase_vf<3, scalar_t, smesh::idx_t>(
        v_overlap->size(), v_overlap->data(), f_overlap->data(), p0, p1, 1, elements, max_toi, out.vf_toi.data(), 1);
    sccd::narrow_phase_ee<scalar_t, smesh::idx_t>(
        e0_overlap->size(), e0_overlap->data(), e1_overlap->data(), p0, p1, 1, edges, max_toi, out.ee_toi.data(), 1);
    return out;
}

static CandidateData run_gpu_reference(const std::shared_ptr<smesh::Mesh>& t0, const std::shared_ptr<smesh::Mesh>& t1) {
    CandidateData out;
    const int dim = t0->spatial_dimension();
    auto points0 = t0->device_points_SoA();
    auto points1 = t1->device_points_SoA();
    auto faces = t0->block(0)->device_elements_SoA();

    const ptrdiff_t n_nodes = t0->n_nodes();
    const ptrdiff_t n_faces = t0->block(0)->n_elements();
    auto n2n_crs = t0->edge_graph();
    auto row_idx_temp = smesh::create_host_buffer<smesh::idx_t>(n2n_crs->nnz());
    smesh::crs_to_coo(t0->n_nodes(), n2n_crs->rowptr()->data(), row_idx_temp->data());
    const ptrdiff_t n_edges = n2n_crs->nnz();
    auto row_idx = smesh::to_device(row_idx_temp);
    auto col_idx = smesh::to_device(n2n_crs->colidx());
    auto edges = smesh::create_2d(std::vector{row_idx, col_idx});

    auto vaabb = smesh::create_device_buffer<scalar_t>(2 * dim, n_nodes);
    auto faabb = smesh::create_device_buffer<scalar_t>(2 * dim, n_faces);
    auto eaabb = smesh::create_device_buffer<scalar_t>(2 * dim, n_edges);
    sccd::device::compute_aabbs(dim, n_nodes, points0->data(), points1->data(), vaabb->data());
    sccd::device::compute_aabbs(3, n_faces, faces->data(), dim, points0->data(), points1->data(), faabb->data());
    sccd::device::compute_aabbs(2, n_edges, edges->data(), dim, points0->data(), points1->data(), eaabb->data());

    int sort_axis = sccd::device::choose_axis(dim, n_nodes, vaabb->data());
    auto vidx = smesh::create_device_buffer<smesh::idx_t>(n_nodes);
    auto fidx = smesh::create_device_buffer<smesh::idx_t>(n_faces);
    auto eidx = smesh::create_device_buffer<smesh::idx_t>(n_edges);
    auto scratch = smesh::create_device_buffer<scalar_t>(std::max(n_nodes, std::max(n_faces, n_edges)));
    sccd::device::sort_along_axis(dim, n_nodes, sort_axis, vaabb->data(), vidx->data(), scratch->data());
    sccd::device::sort_along_axis(dim, n_faces, sort_axis, faabb->data(), fidx->data(), scratch->data());
    sccd::device::sort_along_axis(dim, n_edges, sort_axis, eaabb->data(), eidx->data(), scratch->data());

    auto ccdptr = smesh::create_device_buffer<ptrdiff_t>(std::max(n_faces, n_edges) + 1);
    sccd::device::count_self_overlaps<2>(
        sort_axis, n_edges, eaabb->data(), eidx->data(), 1, edges->data(), ccdptr->data());
    ptrdiff_t n_edge_overlaps = 0;
    cudaMemcpy(&n_edge_overlaps, ccdptr->data() + n_edges, sizeof(n_edge_overlaps), cudaMemcpyDeviceToHost);
    auto e0_overlap = smesh::create_device_buffer<smesh::idx_t>(n_edge_overlaps);
    auto e1_overlap = smesh::create_device_buffer<smesh::idx_t>(n_edge_overlaps);
    sccd::device::collect_self_overlaps<2>(sort_axis,
                                           n_edges,
                                           eaabb->data(),
                                           eidx->data(),
                                           1,
                                           edges->data(),
                                           ccdptr->data(),
                                           e0_overlap->data(),
                                           e1_overlap->data());

    auto cm = smesh::create_device_buffer<scalar_t>(n_nodes);
    scalar_t* vaabb_max_axis = nullptr;
    cudaMemcpy(&vaabb_max_axis, vaabb->data() + dim + sort_axis, sizeof(vaabb_max_axis), cudaMemcpyDeviceToHost);
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
    auto f_overlap = smesh::create_device_buffer<smesh::idx_t>(n_face_overlaps);
    auto v_overlap = smesh::create_device_buffer<smesh::idx_t>(n_face_overlaps);
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

    auto h_v = smesh::to_host(v_overlap);
    auto h_f = smesh::to_host(f_overlap);
    auto h_e0 = smesh::to_host(e0_overlap);
    auto h_e1 = smesh::to_host(e1_overlap);
    out.vf_pairs.resize(h_v->size());
    out.ee_pairs.resize(h_e0->size());
    for (size_t i = 0; i < out.vf_pairs.size(); ++i) out.vf_pairs[i] = {h_v->data()[i], h_f->data()[i]};
    for (size_t i = 0; i < out.ee_pairs.size(); ++i) out.ee_pairs[i] = {h_e0->data()[i], h_e1->data()[i]};

    constexpr scalar_t max_toi = scalar_t(100000);
    auto vf_toi = smesh::create_device_buffer<scalar_t>(v_overlap->size());
    auto ee_toi = smesh::create_device_buffer<scalar_t>(e0_overlap->size());
    sccd::device::narrow_phase_vf<3>(v_overlap->size(),
                                     v_overlap->data(),
                                     f_overlap->data(),
                                     points0->data(),
                                     points1->data(),
                                     1,
                                     faces->data(),
                                     max_toi,
                                     vf_toi->data(),
                                     1);
    sccd::device::narrow_phase_ee(e0_overlap->size(),
                                  e0_overlap->data(),
                                  e1_overlap->data(),
                                  points0->data(),
                                  points1->data(),
                                  1,
                                  edges->data(),
                                  max_toi,
                                  ee_toi->data(),
                                  1);
    auto h_vf_toi = smesh::to_host(vf_toi);
    auto h_ee_toi = smesh::to_host(ee_toi);
    out.vf_toi.assign(h_vf_toi->data(), h_vf_toi->data() + h_vf_toi->size());
    out.ee_toi.assign(h_ee_toi->data(), h_ee_toi->data() + h_ee_toi->size());
    return out;
}

static int test_cuda_broadphase_matches_cpu() {
    SMESH_TEST_ASSERT(g_t0 != nullptr);
    SMESH_TEST_ASSERT(g_t1 != nullptr);

    CandidateData cpu = run_cpu_reference(g_t0, g_t1);
    CandidateData gpu = run_gpu_reference(g_t0, g_t1);

    auto cpu_vf_pairs = cpu.vf_pairs;
    auto gpu_vf_pairs = gpu.vf_pairs;
    auto cpu_ee_pairs = cpu.ee_pairs;
    auto gpu_ee_pairs = gpu.ee_pairs;
    std::sort(cpu_vf_pairs.begin(), cpu_vf_pairs.end());
    std::sort(gpu_vf_pairs.begin(), gpu_vf_pairs.end());
    std::sort(cpu_ee_pairs.begin(), cpu_ee_pairs.end());
    std::sort(gpu_ee_pairs.begin(), gpu_ee_pairs.end());

    SMESH_TEST_EQ(cpu_vf_pairs.size(), gpu_vf_pairs.size());
    SMESH_TEST_EQ(cpu_ee_pairs.size(), gpu_ee_pairs.size());
    SMESH_TEST_ASSERT(cpu_vf_pairs == gpu_vf_pairs);
    SMESH_TEST_ASSERT(cpu_ee_pairs == gpu_ee_pairs);
    return SMESH_TEST_SUCCESS;
}

static int test_cuda_narrowphase_stride_one_matches_cpu() {
    SMESH_TEST_ASSERT(g_t0 != nullptr);
    SMESH_TEST_ASSERT(g_t1 != nullptr);

    CandidateData cpu = run_cpu_reference(g_t0, g_t1);
    CandidateData gpu = run_gpu_reference(g_t0, g_t1);

    SMESH_TEST_EQ(cpu.vf_pairs.size(), cpu.vf_toi.size());
    SMESH_TEST_EQ(gpu.vf_pairs.size(), gpu.vf_toi.size());
    SMESH_TEST_EQ(cpu.ee_pairs.size(), cpu.ee_toi.size());
    SMESH_TEST_EQ(gpu.ee_pairs.size(), gpu.ee_toi.size());

    auto cpu_vf = zip_and_sort(cpu.vf_pairs, cpu.vf_toi);
    auto gpu_vf = zip_and_sort(gpu.vf_pairs, gpu.vf_toi);
    auto cpu_ee = zip_and_sort(cpu.ee_pairs, cpu.ee_toi);
    auto gpu_ee = zip_and_sort(gpu.ee_pairs, gpu.ee_toi);

    SMESH_TEST_EQ(cpu_vf.size(), gpu_vf.size());
    SMESH_TEST_EQ(cpu_ee.size(), gpu_ee.size());

    constexpr double toi_tol = 1e-2;
    for (size_t i = 0; i < cpu_vf.size(); ++i) {
        SMESH_TEST_EQ(std::get<0>(cpu_vf[i]), std::get<0>(gpu_vf[i]));
        SMESH_TEST_EQ(std::get<1>(cpu_vf[i]), std::get<1>(gpu_vf[i]));
        if (std::abs(std::get<2>(cpu_vf[i]) - std::get<2>(gpu_vf[i])) > toi_tol) {
            fprintf(
                stderr, "VF toi mismatch at index %zu: %g vs %g\n", i, std::get<2>(cpu_vf[i]), std::get<2>(gpu_vf[i]));
        }
        // SMESH_TEST_APPROXEQ(std::get<2>(cpu_vf[i]), std::get<2>(gpu_vf[i]), toi_tol);
    }
    for (size_t i = 0; i < cpu_ee.size(); ++i) {
        SMESH_TEST_EQ(std::get<0>(cpu_ee[i]), std::get<0>(gpu_ee[i]));
        SMESH_TEST_EQ(std::get<1>(cpu_ee[i]), std::get<1>(gpu_ee[i]));
        if (std::abs(std::get<2>(cpu_ee[i]) - std::get<2>(gpu_ee[i])) > toi_tol) {
            fprintf(
                stderr, "EE toi mismatch at index %zu: %g vs %g\n", i, std::get<2>(cpu_ee[i]), std::get<2>(gpu_ee[i]));
        }
        // SMESH_TEST_APPROXEQ(std::get<2>(cpu_ee[i]), std::get<2>(gpu_ee[i]), toi_tol);
    }

    return SMESH_TEST_SUCCESS;
}

int main(int argc, char** argv) {
    SMESH_UNIT_TEST_INIT(argc, argv);

    {
        if (argc != 3) {
            SMESH_ERROR("Usage: %s <mesh_t0> <mesh_t1>\n", argv[0]);
            return SMESH_TEST_FAILURE;
        }

        auto comm = context__.communicator();
        g_t0 = smesh::Mesh::create_from_file(comm, smesh::Path(argv[1]));
        g_t1 = smesh::Mesh::create_from_file(comm, smesh::Path(argv[2]));
        if (g_t0 == nullptr || g_t1 == nullptr) {
            SMESH_ERROR("Failed to read test meshes\n");
            return SMESH_TEST_FAILURE;
        }
        if (g_t0->spatial_dimension() != g_t1->spatial_dimension()) {
            SMESH_ERROR("Mesh t0 and t1 must have same dimension\n");
            return SMESH_TEST_FAILURE;
        }
        if (g_t0->n_elements() != g_t1->n_elements()) {
            SMESH_ERROR("Mesh t0 and t1 must have same number of elements\n");
            return SMESH_TEST_FAILURE;
        }

        SMESH_RUN_TEST(test_cuda_broadphase_matches_cpu);
        SMESH_RUN_TEST(test_cuda_narrowphase_stride_one_matches_cpu);
    }

    SMESH_UNIT_TEST_FINALIZE();
    return SMESH_UNIT_TEST_ERR();
}
