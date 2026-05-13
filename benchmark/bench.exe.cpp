#include "sccd_smesh_CCD.hpp"
#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {

using scalar_t = double;
using idx_t = smesh::idx_t;

struct CaseFile {
    fs::path dir;
    std::string key;
    bool is_vf = false;
};

struct MeshPair {
    std::shared_ptr<smesh::Mesh> t0;
    std::shared_ptr<smesh::Mesh> t1;
};

std::uint64_t pair_key(const std::int64_t a, const std::int64_t b)
{
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32)
           | static_cast<std::uint32_t>(b);
}

template <typename T>
bool read_raw(const fs::path& path, std::vector<T>& values)
{
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) {
        std::cerr << "error: failed to open " << path << "\n";
        return false;
    }

    const auto bytes = in.tellg();
    if (bytes < 0 || (static_cast<std::uintmax_t>(bytes) % sizeof(T)) != 0) {
        std::cerr << "error: invalid raw size for " << path << "\n";
        return false;
    }

    values.resize(static_cast<std::size_t>(bytes) / sizeof(T));
    in.seekg(0);
    in.read(reinterpret_cast<char*>(values.data()), bytes);
    if (!in) {
        std::cerr << "error: failed to read " << path << "\n";
        return false;
    }
    return true;
}

template <typename T>
bool write_raw(const fs::path& path, const std::vector<T>& values)
{
    std::error_code ec;
    fs::create_directories(path.parent_path(), ec);
    if (ec) {
        std::cerr << "error: failed to create " << path.parent_path() << ": " << ec.message() << "\n";
        return false;
    }

    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        std::cerr << "error: failed to open " << path << "\n";
        return false;
    }
    out.write(reinterpret_cast<const char*>(values.data()),
              static_cast<std::streamsize>(values.size() * sizeof(T)));
    return static_cast<bool>(out);
}

std::vector<CaseFile> scan_cases(const fs::path& boxes_dir)
{
    std::vector<CaseFile> cases;
    if (!fs::is_directory(boxes_dir)) {
        return cases;
    }

    for (const auto& entry : fs::directory_iterator(boxes_dir)) {
        if (!entry.is_directory()) {
            continue;
        }
        const fs::path dir = entry.path();
        const std::string key = dir.filename().string();
        if ((key.size() < 3) || !fs::exists(dir / "c0.int32") || !fs::exists(dir / "c1.int32")) {
            continue;
        }
        if (key.rfind("vf") == key.size() - 2) {
            cases.push_back({dir, key, true});
        } else if (key.rfind("ee") == key.size() - 2) {
            cases.push_back({dir, key, false});
        }
    }

    std::sort(cases.begin(), cases.end(), [](const CaseFile& a, const CaseFile& b) { return a.key < b.key; });
    return cases;
}

int case_step(const std::string& key)
{
    const std::size_t suffix = (key.size() >= 2 && (key.rfind("vf") == key.size() - 2 || key.rfind("ee") == key.size() - 2)) ? 2 : 0;
    return std::stoi(key.substr(0, key.size() - suffix));
}

fs::path frame_path(const fs::path& dataset_dir, const std::string& dataset, const int step)
{
    const fs::path frames = dataset_dir / "frames";
    if (dataset == "cloth-ball") {
        return frames / ("cloth_ball" + std::to_string(step) + ".ply");
    }
    if (dataset == "n-body-simulation") {
        return frames / ("balls16_" + std::to_string(step) + ".ply");
    }
    return frames / (std::to_string(step) + ".ply");
}

bool load_mesh_pair(const std::shared_ptr<smesh::Communicator>& comm,
                    const fs::path& dataset_dir,
                    const std::string& dataset,
                    const std::string& key,
                    MeshPair& meshes)
{
    const int step = case_step(key);
    const fs::path t0_path = frame_path(dataset_dir, dataset, step);
    const fs::path t1_path = frame_path(dataset_dir, dataset, step + 1);
    if (!fs::exists(t0_path) || !fs::exists(t1_path)) {
        std::cerr << "error: missing frame pair " << t0_path << " / " << t1_path << "\n";
        return false;
    }

    meshes.t0 = smesh::Mesh::create_from_file(comm, smesh::Path(t0_path.string()));
    meshes.t1 = smesh::Mesh::create_from_file(comm, smesh::Path(t1_path.string()));
    return meshes.t0 && meshes.t1;
}

bool normalize_pairs(const std::vector<std::int32_t>& c0,
                     const std::vector<std::int32_t>& c1,
                     const bool is_vf,
                     const ptrdiff_t n_nodes,
                     const ptrdiff_t n_edges,
                     std::vector<idx_t>& a,
                     std::vector<idx_t>& b)
{
    if (c0.size() != c1.size()) {
        std::cerr << "error: c0/c1 size mismatch\n";
        return false;
    }

    const ptrdiff_t edge_offset = n_nodes;
    const ptrdiff_t face_offset = n_nodes + n_edges;
    a.resize(c0.size());
    b.resize(c1.size());

    for (std::size_t i = 0; i < c0.size(); ++i) {
        const ptrdiff_t x = c0[i];
        const ptrdiff_t y = c1[i];
        if (is_vf) {
            if (x < n_nodes && y >= face_offset) {
                a[i] = static_cast<idx_t>(x);
                b[i] = static_cast<idx_t>(y - face_offset);
            } else if (y < n_nodes && x >= face_offset) {
                a[i] = static_cast<idx_t>(y);
                b[i] = static_cast<idx_t>(x - face_offset);
            } else {
                std::cerr << "error: invalid vf pair (" << x << ", " << y << ")\n";
                return false;
            }
        } else {
            a[i] = static_cast<idx_t>(x - edge_offset);
            b[i] = static_cast<idx_t>(y - edge_offset);
        }
    }
    return true;
}

template <typename Buffer>
std::vector<typename Buffer::element_type> to_vector(const Buffer& buffer)
{
    auto host = smesh::to_host(buffer);
    return std::vector<typename Buffer::element_type>(host->data(), host->data() + host->size());
}

std::unordered_set<std::uint64_t> expected_pairs(const std::vector<std::int32_t>& c0,
                                                 const std::vector<std::int32_t>& c1)
{
    std::unordered_set<std::uint64_t> set;
    set.reserve(c0.size() * 2 + 1);
    for (std::size_t i = 0; i < c0.size(); ++i) {
        set.insert(pair_key(c0[i], c1[i]));
    }
    return set;
}

std::unordered_set<std::uint64_t> broad_pairs_for(const MeshPair& meshes,
                                                  const std::vector<std::int32_t>& c0,
                                                  const std::vector<std::int32_t>& c1,
                                                  const bool is_vf,
                                                  std::uint64_t& false_positives)
{
    SMESH_TRACE_SCOPE("benchmark broadphase");

    auto points0 = smesh::astype<scalar_t>(meshes.t0->points());
    auto points1 = smesh::astype<scalar_t>(meshes.t1->points());
    const int dim = meshes.t0->spatial_dimension();
    const ptrdiff_t n_nodes = meshes.t0->n_nodes();
    const ptrdiff_t n_faces = meshes.t0->block(0)->n_elements();
    auto graph = meshes.t0->edge_graph();
    const ptrdiff_t n_edges = graph->nnz();
    const ptrdiff_t face_offset = n_nodes + n_edges;
    const ptrdiff_t edge_offset = n_nodes;
    auto faces = meshes.t0->block(0)->elements();
    auto e0 = smesh::create_host_buffer<idx_t>(n_edges);
    smesh::crs_to_coo(n_nodes, graph->rowptr()->data(), e0->data());
    auto edges = smesh::create_2d(std::vector{e0, graph->colidx()});

    auto vaabb = smesh::create_buffer<scalar_t>(2 * dim, n_nodes, smesh::EXECUTION_SPACE_HOST);
    auto faabb = smesh::create_buffer<scalar_t>(2 * dim, n_faces, smesh::EXECUTION_SPACE_HOST);
    auto eaabb = smesh::create_buffer<scalar_t>(2 * dim, n_edges, smesh::EXECUTION_SPACE_HOST);
    auto vidx = smesh::create_buffer<idx_t>(n_nodes, smesh::EXECUTION_SPACE_HOST);
    auto fidx = smesh::create_buffer<idx_t>(n_faces, smesh::EXECUTION_SPACE_HOST);
    auto eidx = smesh::create_buffer<idx_t>(n_edges, smesh::EXECUTION_SPACE_HOST);
    auto scratch = smesh::create_buffer<scalar_t>(std::max(n_nodes, std::max(n_faces, n_edges)), smesh::EXECUTION_SPACE_HOST);
    auto ccdptr = smesh::create_buffer<ptrdiff_t>(std::max(n_faces, n_edges) + 1, smesh::EXECUTION_SPACE_HOST);
    auto cumulative_max = smesh::create_buffer<scalar_t>(n_nodes, smesh::EXECUTION_SPACE_HOST);

    sccd::compute_aabbs(dim, n_nodes, points0->data(), points1->data(), vaabb->data(), vaabb->data() + dim);
    sccd::compute_aabbs(
        meshes.t0->block(0)->n_nodes_per_element(), n_faces, faces->data(), dim, points0->data(), points1->data(), faabb->data(), faabb->data() + dim);
    sccd::compute_aabbs(2, n_edges, edges->data(), dim, points0->data(), points1->data(), eaabb->data(), eaabb->data() + dim);

    const int sort_axis = sccd::choose_axis(n_nodes, vaabb->data());
    sccd::sort_along_axis(n_nodes, sort_axis, vaabb->data(), vidx->data(), scratch->data());

    std::unordered_set<std::uint64_t> actual;
    if (is_vf) {
        sccd::sort_along_axis(n_faces, sort_axis, faabb->data(), fidx->data(), scratch->data());
        sccd::cummax(n_nodes, vaabb->data()[dim + sort_axis], cumulative_max->data());
        sccd::count_overlaps<3, 1, scalar_t, idx_t>(sort_axis,
                                                    n_faces,
                                                    faabb->data(),
                                                    fidx->data(),
                                                    1,
                                                    faces->data(),
                                                    n_nodes,
                                                    vaabb->data(),
                                                    vidx->data(),
                                                    0,
                                                    nullptr,
                                                    ccdptr->data(),
                                                    cumulative_max->data());
        const ptrdiff_t n_overlaps = ccdptr->data()[n_faces];
        auto f_overlap = smesh::create_buffer<idx_t>(n_overlaps, smesh::EXECUTION_SPACE_HOST);
        auto v_overlap = smesh::create_buffer<idx_t>(n_overlaps, smesh::EXECUTION_SPACE_HOST);
        sccd::collect_overlaps<3, 1, scalar_t, idx_t>(sort_axis,
                                                      n_faces,
                                                      faabb->data(),
                                                      fidx->data(),
                                                      1,
                                                      faces->data(),
                                                      n_nodes,
                                                      vaabb->data(),
                                                      vidx->data(),
                                                      0,
                                                      nullptr,
                                                      ccdptr->data(),
                                                      cumulative_max->data(),
                                                      f_overlap->data(),
                                                      v_overlap->data());
        actual.reserve(static_cast<std::size_t>(n_overlaps) * 2 + 1);
        for (ptrdiff_t i = 0; i < n_overlaps; ++i) {
            actual.insert(pair_key(v_overlap->data()[i], f_overlap->data()[i] + face_offset));
        }
    } else {
        sccd::sort_along_axis(n_edges, sort_axis, eaabb->data(), eidx->data(), scratch->data());
        sccd::count_self_overlaps<2>(sort_axis, n_edges, eaabb->data(), eidx->data(), 1, edges->data(), ccdptr->data());
        const ptrdiff_t n_overlaps = ccdptr->data()[n_edges];
        auto e0_overlap = smesh::create_buffer<idx_t>(n_overlaps, smesh::EXECUTION_SPACE_HOST);
        auto e1_overlap = smesh::create_buffer<idx_t>(n_overlaps, smesh::EXECUTION_SPACE_HOST);
        sccd::collect_self_overlaps<2>(sort_axis,
                                       n_edges,
                                       eaabb->data(),
                                       eidx->data(),
                                       1,
                                       edges->data(),
                                       ccdptr->data(),
                                       e0_overlap->data(),
                                       e1_overlap->data());
        actual.reserve(static_cast<std::size_t>(n_overlaps) * 2 + 1);
        for (ptrdiff_t i = 0; i < n_overlaps; ++i) {
            actual.insert(pair_key(e0_overlap->data()[i] + edge_offset, e1_overlap->data()[i] + edge_offset));
        }
    }

    const auto expected = expected_pairs(c0, c1);
    false_positives = 0;
    for (const auto pair : actual) {
        if (expected.find(pair) == expected.end()) {
            ++false_positives;
        }
    }
    return actual;
}

bool run_case(const std::shared_ptr<smesh::Communicator>& comm,
              const fs::path& dataset_dir,
              const std::string& dataset,
              const CaseFile& case_file,
              std::vector<double>& timings_ms)
{
    std::vector<std::int32_t> c0;
    std::vector<std::int32_t> c1;
    if (!read_raw(case_file.dir / "c0.int32", c0) || !read_raw(case_file.dir / "c1.int32", c1)) {
        return false;
    }

    MeshPair meshes;
    if (!load_mesh_pair(comm, dataset_dir, dataset, case_file.key, meshes)) {
        return false;
    }

    const std::string trace_file = dataset + "/" + case_file.key;
    setenv("SMESH_TRACE_FILE", trace_file.c_str(), 1);

    const ptrdiff_t n_nodes = meshes.t0->n_nodes();
    const ptrdiff_t n_edges = meshes.t0->edge_graph()->nnz();
    std::vector<idx_t> q0;
    std::vector<idx_t> q1;
    if (!normalize_pairs(c0, c1, case_file.is_vf, n_nodes, n_edges, q0, q1)) {
        return false;
    }

    std::uint64_t broad_fp_count = 0;
    const auto broad_actual = broad_pairs_for(meshes, c0, c1, case_file.is_vf, broad_fp_count);

    std::vector<std::uint8_t> mma;
    const fs::path mma_path = dataset_dir / "mma_bool" / case_file.key / "mma_bool.uint8";
    if (fs::exists(mma_path) && !read_raw(mma_path, mma)) {
        return false;
    }

    std::vector<scalar_t> root_toi;
    const fs::path roots_dir = dataset_dir / "roots" / case_file.key;
    const fs::path root_path = roots_dir / "toi.float64";
    if (fs::exists(root_path) && !read_raw(root_path, root_toi)) {
        return false;
    }

    auto points0 = smesh::astype<scalar_t>(meshes.t0->points());
    auto points1 = smesh::astype<scalar_t>(meshes.t1->points());
    std::vector<scalar_t> sccd_toi(q0.size(), scalar_t(1));

    const auto start = std::chrono::steady_clock::now();
    if (case_file.is_vf) {
        auto faces = meshes.t0->block(0)->elements();
        sccd::narrow_phase_vf<3, scalar_t, idx_t>(
            q0.size(), q0.data(), q1.data(), points0->data(), points1->data(), 1, faces->data(), 1.0, sccd_toi.data(), 1);
    } else {
        auto graph = meshes.t0->edge_graph();
        auto e0 = smesh::create_host_buffer<idx_t>(graph->nnz());
        smesh::crs_to_coo(meshes.t0->n_nodes(), graph->rowptr()->data(), e0->data());
        auto edges = smesh::create_2d(std::vector{e0, graph->colidx()});
        sccd::narrow_phase_ee<scalar_t, idx_t>(
            q0.size(), q0.data(), q1.data(), points0->data(), points1->data(), 1, edges->data(), 1.0, sccd_toi.data(), 1);
    }
    const auto stop = std::chrono::steady_clock::now();
    timings_ms.push_back(std::chrono::duration<double, std::milli>(stop - start).count());

    std::vector<std::uint8_t> fp(q0.size(), 0);
    std::vector<std::uint8_t> fn(q0.size(), 0);
    std::vector<std::uint8_t> fp_broad(q0.size(), 0);
    std::vector<std::uint8_t> fn_broad(q0.size(), 0);

    for (std::size_t i = 0; i < q0.size(); ++i) {
        bool expected = false;
        if (i < mma.size()) {
            expected = mma[i] != 0;
        } else if (i < root_toi.size()) {
            expected = std::isfinite(root_toi[i]);
        }

        const bool found = sccd_toi[i] < scalar_t(1);
        fp[i] = static_cast<std::uint8_t>(found && !expected);
        fn[i] = static_cast<std::uint8_t>(!found && expected);

        const bool broad_found = broad_actual.find(pair_key(c0[i], c1[i])) != broad_actual.end()
                                 || broad_actual.find(pair_key(c1[i], c0[i])) != broad_actual.end();
        fn_broad[i] = static_cast<std::uint8_t>(!broad_found);
    }

    if (broad_fp_count != 0) {
        std::cerr << "warning: " << dataset << "/" << case_file.key << " has " << broad_fp_count
                  << " broad-phase false positives outside the query order\n";
    }

    return write_raw(roots_dir / "sccd_toi.float64", sccd_toi)
           && write_raw(roots_dir / "sccd_fp.uint8", fp)
           && write_raw(roots_dir / "sccd_fn.uint8", fn)
           && write_raw(roots_dir / "sccd_fp_broad.uint8", fp_broad)
           && write_raw(roots_dir / "sccd_fn_broad.uint8", fn_broad);
}

} // namespace

int main(int argc, char** argv)
{
    auto ctx = smesh::initialize(argc, argv);
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <data-dir> <dataset> [<dataset> ...]\n";
        return EXIT_FAILURE;
    }

    const fs::path data_dir = argv[1];
    bool ok = true;
    for (int i = 2; i < argc; ++i) {
        const std::string dataset = argv[i];
        const fs::path dataset_dir = data_dir / dataset;
        std::vector<double> vf_timings;
        std::vector<double> ee_timings;

        for (const CaseFile& case_file : scan_cases(dataset_dir / "boxes")) {
            ok = run_case(ctx->communicator(), dataset_dir, dataset, case_file, case_file.is_vf ? vf_timings : ee_timings) && ok;
        }

        const fs::path out_dir = dataset_dir / "benchmark";
        ok = write_raw(out_dir / (dataset + "-vf.float64"), vf_timings) && ok;
        ok = write_raw(out_dir / (dataset + "-ee.float64"), ee_timings) && ok;
    }

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
