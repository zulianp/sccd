#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_test.hpp"

#include "sccd_base.hpp"
#include "sccd_config.hpp"
#include "sccd_smesh_ccd.hpp"

#include <algorithm>
#include <cmath>
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
    int err{SCCD_SUCCESS};
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

static void copy_pairs(const smesh::SharedBuffer<smesh::idx_t>& first,
                       const smesh::SharedBuffer<smesh::idx_t>& second,
                       std::vector<pair_t>& out) {
    auto h_first = smesh::to_host(first);
    auto h_second = smesh::to_host(second);

    out.resize(h_first->size());
    for (size_t i = 0; i < out.size(); ++i) {
        out[i] = {h_first->data()[i], h_second->data()[i]};
    }
}

static void copy_values(const smesh::SharedBuffer<scalar_t>& values, std::vector<scalar_t>& out) {
    auto h_values = smesh::to_host(values);
    out.assign(h_values->data(), h_values->data() + h_values->size());
}

static CandidateData run_ccd(const std::shared_ptr<smesh::Mesh>& t0,
                             const std::shared_ptr<smesh::Mesh>& t1,
                             const smesh::ExecutionSpace execution_space) {
    CandidateData out;

    smesh::SharedBuffer<scalar_t*> points0;
    smesh::SharedBuffer<scalar_t*> points1;

    if (execution_space == smesh::EXECUTION_SPACE_DEVICE) {
        points0 = t0->device_points_SoA();
        points1 = t1->device_points_SoA();
    } else {
        points0 = t0->points();
        points1 = t1->points();
    }

    smesh::SharedBuffer<smesh::idx_t> v_overlap;
    smesh::SharedBuffer<smesh::idx_t> f_overlap;
    smesh::SharedBuffer<smesh::idx_t> e0_overlap;
    smesh::SharedBuffer<smesh::idx_t> e1_overlap;
    smesh::SharedBuffer<scalar_t> vf_toi;
    smesh::SharedBuffer<scalar_t> ee_toi;

    auto ccd = sccd::CCD<scalar_t>::create(t0, execution_space);
    out.err = ccd->find_impact_times(
        points0, points1, v_overlap, f_overlap, vf_toi, e0_overlap, e1_overlap, ee_toi, 69, scalar_t(3e-8));
    if (out.err != SCCD_SUCCESS) {
        return out;
    }

    copy_pairs(v_overlap, f_overlap, out.vf_pairs);
    copy_pairs(e0_overlap, e1_overlap, out.ee_pairs);
    copy_values(vf_toi, out.vf_toi);
    copy_values(ee_toi, out.ee_toi);
    return out;
}

static CandidateData run_ccd_cpu(const std::shared_ptr<smesh::Mesh>& t0, const std::shared_ptr<smesh::Mesh>& t1) {
    return run_ccd(t0, t1, smesh::EXECUTION_SPACE_HOST);
}

static CandidateData run_ccd_gpu(const std::shared_ptr<smesh::Mesh>& t0, const std::shared_ptr<smesh::Mesh>& t1) {
    return run_ccd(t0, t1, smesh::EXECUTION_SPACE_DEVICE);
}

static int test_cuda_broadphase_matches_cpu() {
    SMESH_TEST_ASSERT(g_t0 != nullptr);
    SMESH_TEST_ASSERT(g_t1 != nullptr);

    CandidateData cpu = run_ccd_cpu(g_t0, g_t1);
    CandidateData gpu = run_ccd_gpu(g_t0, g_t1);
    SMESH_TEST_EQ(cpu.err, SCCD_SUCCESS);
    SMESH_TEST_EQ(gpu.err, SCCD_SUCCESS);

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

    CandidateData cpu = run_ccd_cpu(g_t0, g_t1);
    CandidateData gpu = run_ccd_gpu(g_t0, g_t1);
    SMESH_TEST_EQ(cpu.err, SCCD_SUCCESS);
    SMESH_TEST_EQ(gpu.err, SCCD_SUCCESS);

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
    }
    for (size_t i = 0; i < cpu_ee.size(); ++i) {
        SMESH_TEST_EQ(std::get<0>(cpu_ee[i]), std::get<0>(gpu_ee[i]));
        SMESH_TEST_EQ(std::get<1>(cpu_ee[i]), std::get<1>(gpu_ee[i]));
        if (std::abs(std::get<2>(cpu_ee[i]) - std::get<2>(gpu_ee[i])) > toi_tol) {
            fprintf(
                stderr, "EE toi mismatch at index %zu: %g vs %g\n", i, std::get<2>(cpu_ee[i]), std::get<2>(gpu_ee[i]));
        }
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
