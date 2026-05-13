#include "sccd_smesh_CCD.hpp"
#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <algorithm>
#include <array>
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

    struct QueryGeometry {
        std::array<std::vector<scalar_t>, 3> points0;
        std::array<std::vector<scalar_t>, 3> points1;
        std::vector<idx_t> q0;
        std::vector<idx_t> q1;
        std::array<std::vector<idx_t>, 3> faces;
        std::array<std::vector<idx_t>, 2> edges;
    };

    std::uint64_t pair_key(const std::int64_t a, const std::int64_t b) {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) | static_cast<std::uint32_t>(b);
    }

    template <typename T>
    bool read_raw(const fs::path& path, std::vector<T>& values) {
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
    bool write_raw(const fs::path& path, const std::vector<T>& values) {
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

    std::vector<CaseFile> scan_cases(const fs::path& boxes_dir) {
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

    int case_step(const std::string& key) {
        const std::size_t suffix =
            (key.size() >= 2 && (key.rfind("vf") == key.size() - 2 || key.rfind("ee") == key.size() - 2)) ? 2 : 0;
        return std::stoi(key.substr(0, key.size() - suffix));
    }

    fs::path frame_path(const fs::path& dataset_dir, const std::string& dataset, const int step) {
        const fs::path frames = dataset_dir / "frames";
        if (dataset == "cloth-ball") {
            return frames / ("cloth_ball" + std::to_string(step) + ".ply");
        }
        if (dataset == "n-body-simulation") {
            return frames / ("balls16_" + std::to_string(step) + ".ply");
        }
        return frames / (std::to_string(step) + ".ply");
    }

    std::string shell_quote(const fs::path& path) {
        std::string quoted = "'";
        for (const char c : path.string()) {
            if (c == '\'') {
                quoted += "'\\''";
            } else {
                quoted += c;
            }
        }
        quoted += "'";
        return quoted;
    }

    bool raw_mesh_is_current(const fs::path& input_mesh, const fs::path& output_dir) {
        std::error_code ec;
        if (!fs::is_directory(output_dir, ec)) {
            return false;
        }

        const fs::file_time_type input_time = fs::last_write_time(input_mesh, ec);
        if (ec) {
            return false;
        }

        for (const char* required : {"x.float32", "y.float32", "z.float32", "i0.int32", "i1.int32", "i2.int32"}) {
            const fs::path raw_file = output_dir / required;
            if (!fs::exists(raw_file, ec)) {
                return false;
            }
            const fs::file_time_type raw_time = fs::last_write_time(raw_file, ec);
            if (ec || raw_time < input_time) {
                return false;
            }
        }
        return true;
    }

    bool convert_mesh_to_raw(const fs::path& input_mesh, const fs::path& output_dir) {
        if (raw_mesh_is_current(input_mesh, output_dir)) {
            return true;
        }

        const char* db_to_raw_env = std::getenv("SCCD_DB_TO_RAW");
        const fs::path db_to_raw =
            (db_to_raw_env && db_to_raw_env[0] != '\0') ? fs::path(db_to_raw_env) : fs::path("db_to_raw");

        std::error_code ec;
        fs::create_directories(output_dir, ec);
        if (ec) {
            std::cerr << "error: failed to create " << output_dir << ": " << ec.message() << "\n";
            return false;
        }

        const std::string command =
            shell_quote(db_to_raw) + " " + shell_quote(input_mesh) + " " + shell_quote(output_dir);
        const int status = std::system(command.c_str());
        if (status != 0) {
            std::cerr << "error: db_to_raw failed for " << input_mesh << "\n";
            return false;
        }
        return true;
    }

    fs::path raw_frame_path(const fs::path& dataset_dir, const fs::path& frame) {
        return dataset_dir / "frames_raw" / frame.stem();
    }

    fs::path query_path(const fs::path& dataset_dir, const std::string& key) {
        return dataset_dir / "queries" / (key + ".csv");
    }

    bool parse_query_row(const std::string& line, scalar_t coords[3]) {
        long double values[6];
        const char* cursor = line.c_str();
        for (int i = 0; i < 6; ++i) {
            char* end = nullptr;
            values[i] = std::strtold(cursor, &end);
            if (end == cursor || (i < 5 && *end != ',')) {
                return false;
            }
            cursor = (i < 5) ? end + 1 : end;
        }
        coords[0] = static_cast<scalar_t>(values[0] / values[1]);
        coords[1] = static_cast<scalar_t>(values[2] / values[3]);
        coords[2] = static_cast<scalar_t>(values[4] / values[5]);
        return true;
    }

    bool read_query_geometry(const fs::path& path, const bool is_vf, QueryGeometry& query) {
        std::ifstream in(path);
        if (!in) {
            std::cerr << "error: failed to open " << path << "\n";
            return false;
        }

        std::string line;
        std::size_t row = 0;
        while (std::getline(in, line)) {
            if (line.empty()) {
                continue;
            }

            scalar_t coords[3];
            if (!parse_query_row(line, coords)) {
                std::cerr << "error: invalid query row in " << path << "\n";
                return false;
            }

            auto& points = (row % 8 < 4) ? query.points0 : query.points1;
            for (int d = 0; d < 3; ++d) {
                points[d].push_back(coords[d]);
            }
            ++row;
        }

        if (row % 8 != 0 || query.points0[0].size() != query.points1[0].size()) {
            std::cerr << "error: invalid query row count in " << path << "\n";
            return false;
        }

        const std::size_t n_queries = row / 8;
        query.q0.resize(n_queries);
        query.q1.resize(n_queries);
        if (is_vf) {
            for (std::size_t i = 0; i < n_queries; ++i) {
                const idx_t base = static_cast<idx_t>(4 * i);
                query.q0[i] = base;
                query.q1[i] = static_cast<idx_t>(i);
                query.faces[0].push_back(base + 1);
                query.faces[1].push_back(base + 2);
                query.faces[2].push_back(base + 3);
            }
        } else {
            for (std::size_t i = 0; i < n_queries; ++i) {
                const idx_t base = static_cast<idx_t>(4 * i);
                query.q0[i] = static_cast<idx_t>(2 * i);
                query.q1[i] = static_cast<idx_t>(2 * i + 1);
                query.edges[0].push_back(base);
                query.edges[1].push_back(base + 1);
                query.edges[0].push_back(base + 2);
                query.edges[1].push_back(base + 3);
            }
        }

        return true;
    }

    bool load_mesh_pair(const std::shared_ptr<smesh::Communicator>& comm,
                        const fs::path& dataset_dir,
                        const std::string& dataset,
                        const std::string& key,
                        MeshPair& meshes) {
        const int step = case_step(key);
        const fs::path t0_path = frame_path(dataset_dir, dataset, step);
        const fs::path t1_path = frame_path(dataset_dir, dataset, step + 1);
        if (!fs::exists(t0_path) || !fs::exists(t1_path)) {
            std::cerr << "error: missing frame pair " << t0_path << " / " << t1_path << "\n";
            return false;
        }

        const fs::path t0_raw = raw_frame_path(dataset_dir, t0_path);
        const fs::path t1_raw = raw_frame_path(dataset_dir, t1_path);
        if (!convert_mesh_to_raw(t0_path, t0_raw) || !convert_mesh_to_raw(t1_path, t1_raw)) {
            return false;
        }

        meshes.t0 = smesh::Mesh::create_from_file(comm, smesh::Path(t0_raw.string()));
        meshes.t1 = smesh::Mesh::create_from_file(comm, smesh::Path(t1_raw.string()));
        return meshes.t0 && meshes.t1;
    }

    bool normalize_pairs(const std::vector<std::int32_t>& c0,
                         const std::vector<std::int32_t>& c1,
                         const bool is_vf,
                         const ptrdiff_t n_nodes,
                         const ptrdiff_t n_edges,
                         std::vector<idx_t>& a,
                         std::vector<idx_t>& b) {
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
    std::vector<typename Buffer::element_type> to_vector(const Buffer& buffer) {
        auto host = smesh::to_host(buffer);
        return std::vector<typename Buffer::element_type>(host->data(), host->data() + host->size());
    }

    std::unordered_set<std::uint64_t> expected_pairs(const std::vector<std::int32_t>& c0,
                                                     const std::vector<std::int32_t>& c1) {
        std::unordered_set<std::uint64_t> set;
        set.reserve(c0.size() * 2 + 1);
        for (std::size_t i = 0; i < c0.size(); ++i) {
            set.insert(pair_key(c0[i], c1[i]));
        }
        return set;
    }

    bool contains_pair(const std::unordered_set<std::uint64_t>& set, const std::uint64_t pair) {
        const std::uint64_t reversed = (pair << 32) | (pair >> 32);
        return set.find(pair) != set.end() || set.find(reversed) != set.end();
    }

    std::array<std::vector<idx_t>, 2> dataset_ordered_edges(const std::shared_ptr<smesh::Mesh>& mesh) {
        auto faces = mesh->block(0)->elements();
        const ptrdiff_t n_faces = mesh->block(0)->n_elements();
        std::unordered_set<std::uint64_t> seen;
        seen.reserve(static_cast<std::size_t>(3 * n_faces));
        std::vector<std::uint64_t> keys;
        keys.reserve(static_cast<std::size_t>(3 * n_faces));

        for (ptrdiff_t i = 0; i < n_faces; ++i) {
            const idx_t tri[3] = {faces->data()[0][i], faces->data()[1][i], faces->data()[2][i]};
            for (int e = 0; e < 3; ++e) {
                const idx_t a = std::min(tri[e], tri[(e + 1) % 3]);
                const idx_t b = std::max(tri[e], tri[(e + 1) % 3]);
                const std::uint64_t key = pair_key(a, b);
                if (seen.insert(key).second) {
                    keys.push_back(key);
                }
            }
        }

        std::sort(keys.begin(), keys.end(), [](const std::uint64_t a, const std::uint64_t b) {
            const std::uint32_t a0 = static_cast<std::uint32_t>(a >> 32);
            const std::uint32_t a1 = static_cast<std::uint32_t>(a);
            const std::uint32_t b0 = static_cast<std::uint32_t>(b >> 32);
            const std::uint32_t b1 = static_cast<std::uint32_t>(b);
            return a1 < b1 || (a1 == b1 && a0 < b0);
        });

        std::array<std::vector<idx_t>, 2> edges;
        edges[0].reserve(keys.size());
        edges[1].reserve(keys.size());
        for (const std::uint64_t key : keys) {
            edges[0].push_back(static_cast<idx_t>(key >> 32));
            edges[1].push_back(static_cast<idx_t>(key & 0xffffffffu));
        }
        return edges;
    }

    std::unordered_set<std::uint64_t> broad_pairs_for(const MeshPair& meshes,
                                                      const bool is_vf,
                                                      const std::unordered_set<std::uint64_t>& expected,
                                                      std::uint64_t& false_positives) {
        SMESH_TRACE_SCOPE("benchmark broadphase");

        auto points0 = smesh::astype<scalar_t>(meshes.t0->points());
        auto points1 = smesh::astype<scalar_t>(meshes.t1->points());
        const int dim = meshes.t0->spatial_dimension();
        const ptrdiff_t n_nodes = meshes.t0->n_nodes();
        const ptrdiff_t n_faces = meshes.t0->block(0)->n_elements();
        auto ordered_edges = dataset_ordered_edges(meshes.t0);
        const ptrdiff_t n_edges = static_cast<ptrdiff_t>(ordered_edges[0].size());
        const ptrdiff_t face_offset = n_nodes + n_edges;
        const ptrdiff_t edge_offset = n_nodes;
        auto faces = meshes.t0->block(0)->elements();
        idx_t* edges[2] = {ordered_edges[0].data(), ordered_edges[1].data()};

        auto vaabb = smesh::create_buffer<scalar_t>(2 * dim, n_nodes, smesh::EXECUTION_SPACE_HOST);
        auto faabb = smesh::create_buffer<scalar_t>(2 * dim, n_faces, smesh::EXECUTION_SPACE_HOST);
        auto eaabb = smesh::create_buffer<scalar_t>(2 * dim, n_edges, smesh::EXECUTION_SPACE_HOST);
        auto vidx = smesh::create_buffer<idx_t>(n_nodes, smesh::EXECUTION_SPACE_HOST);
        auto fidx = smesh::create_buffer<idx_t>(n_faces, smesh::EXECUTION_SPACE_HOST);
        auto eidx = smesh::create_buffer<idx_t>(n_edges, smesh::EXECUTION_SPACE_HOST);
        auto scratch =
            smesh::create_buffer<scalar_t>(std::max(n_nodes, std::max(n_faces, n_edges)), smesh::EXECUTION_SPACE_HOST);
        auto ccdptr = smesh::create_buffer<ptrdiff_t>(std::max(n_faces, n_edges) + 1, smesh::EXECUTION_SPACE_HOST);
        auto cumulative_max = smesh::create_buffer<scalar_t>(n_nodes, smesh::EXECUTION_SPACE_HOST);

        sccd::compute_aabbs(dim, n_nodes, points0->data(), points1->data(), vaabb->data(), vaabb->data() + dim);
        sccd::compute_aabbs(meshes.t0->block(0)->n_nodes_per_element(),
                            n_faces,
                            faces->data(),
                            dim,
                            points0->data(),
                            points1->data(),
                            faabb->data(),
                            faabb->data() + dim);
        sccd::compute_aabbs(
            2, n_edges, edges, dim, points0->data(), points1->data(), eaabb->data(), eaabb->data() + dim);

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
            sccd::count_self_overlaps<2>(
                sort_axis, n_edges, eaabb->data(), eidx->data(), 1, edges, ccdptr->data());
            const ptrdiff_t n_overlaps = ccdptr->data()[n_edges];
            auto e0_overlap = smesh::create_buffer<idx_t>(n_overlaps, smesh::EXECUTION_SPACE_HOST);
            auto e1_overlap = smesh::create_buffer<idx_t>(n_overlaps, smesh::EXECUTION_SPACE_HOST);
            sccd::collect_self_overlaps<2>(sort_axis,
                                           n_edges,
                                           eaabb->data(),
                                           eidx->data(),
                                           1,
                                           edges,
                                           ccdptr->data(),
                                           e0_overlap->data(),
                                           e1_overlap->data());
            actual.reserve(static_cast<std::size_t>(n_overlaps) * 2 + 1);
            for (ptrdiff_t i = 0; i < n_overlaps; ++i) {
                actual.insert(pair_key(e0_overlap->data()[i] + edge_offset, e1_overlap->data()[i] + edge_offset));
            }
        }

        false_positives = 0;
        for (const auto pair : actual) {
            if (!contains_pair(expected, pair)) {
                ++false_positives;
            }
        }
        return actual;
    }

    bool run_case(const std::shared_ptr<smesh::Communicator>& comm,
                  const fs::path& dataset_dir,
                  const std::string& dataset,
                  const CaseFile& case_file,
                  std::vector<double>& timings_ms) {
        std::vector<std::int32_t> c0;
        std::vector<std::int32_t> c1;
        if (!read_raw(case_file.dir / "c0.int32", c0) || !read_raw(case_file.dir / "c1.int32", c1)) {
            return false;
        }

        MeshPair meshes;
        if (!load_mesh_pair(comm, dataset_dir, dataset, case_file.key, meshes)) {
            return false;
        }

        const fs::path trace_file = dataset_dir / "benchmark" / "traces" / (case_file.key + ".csv");
        std::error_code trace_ec;
        fs::create_directories(trace_file.parent_path(), trace_ec);
        if (trace_ec) {
            std::cerr << "error: failed to create " << trace_file.parent_path() << ": " << trace_ec.message() << "\n";
            return false;
        }
        setenv("SMESH_TRACE_FILE", trace_file.string().c_str(), 1);

        const ptrdiff_t n_nodes = meshes.t0->n_nodes();
        const ptrdiff_t n_edges = meshes.t0->edge_graph()->nnz();
        std::vector<idx_t> q0;
        std::vector<idx_t> q1;
        if (!normalize_pairs(c0, c1, case_file.is_vf, n_nodes, n_edges, q0, q1)) {
            return false;
        }

        QueryGeometry query_geometry;
        const fs::path query_file = query_path(dataset_dir, case_file.key);
        if (!read_query_geometry(query_file, case_file.is_vf, query_geometry)) {
            return false;
        }
        if (query_geometry.q0.size() != q0.size()) {
            std::cerr << "error: query count mismatch for " << query_file << "\n";
            return false;
        }

        const auto broad_expected = expected_pairs(c0, c1);
        std::uint64_t broad_fp_count = 0;
        const auto broad_start = std::chrono::steady_clock::now();
        const auto broad_actual = broad_pairs_for(meshes, case_file.is_vf, broad_expected, broad_fp_count);
        const auto broad_stop = std::chrono::steady_clock::now();
        const double broad_ms = std::chrono::duration<double, std::milli>(broad_stop - broad_start).count();

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

        std::vector<scalar_t> sccd_toi(q0.size(), scalar_t(1));

        const auto start = std::chrono::steady_clock::now();
        scalar_t* points0[3] = {
            query_geometry.points0[0].data(), query_geometry.points0[1].data(), query_geometry.points0[2].data()};
        scalar_t* points1[3] = {
            query_geometry.points1[0].data(), query_geometry.points1[1].data(), query_geometry.points1[2].data()};
        if (case_file.is_vf) {
            idx_t* faces[3] = {
                query_geometry.faces[0].data(), query_geometry.faces[1].data(), query_geometry.faces[2].data()};
            sccd::narrow_phase_vf<3, scalar_t, idx_t>(query_geometry.q0.size(),
                                                      query_geometry.q0.data(),
                                                      query_geometry.q1.data(),
                                                      points0,
                                                      points1,
                                                      1,
                                                      faces,
                                                      1.0,
                                                      sccd_toi.data(),
                                                      1);
        } else {
            idx_t* edges[2] = {query_geometry.edges[0].data(), query_geometry.edges[1].data()};
            sccd::narrow_phase_ee<scalar_t, idx_t>(query_geometry.q0.size(),
                                                   query_geometry.q0.data(),
                                                   query_geometry.q1.data(),
                                                   points0,
                                                   points1,
                                                   1,
                                                   edges,
                                                   1.0,
                                                   sccd_toi.data(),
                                                   1);
        }
        const auto stop = std::chrono::steady_clock::now();
        const double narrow_ms = std::chrono::duration<double, std::milli>(stop - start).count();
        timings_ms.push_back(narrow_ms);

        std::vector<std::uint8_t> fp(q0.size(), 0);
        std::vector<std::uint8_t> fn(q0.size(), 0);
        std::vector<std::uint8_t> fp_broad;
        std::vector<std::uint8_t> fn_broad(q0.size(), 0);
        std::vector<std::int32_t> fp_broad_c0;
        std::vector<std::int32_t> fp_broad_c1;
        std::vector<std::int32_t> fn_broad_c0;
        std::vector<std::int32_t> fn_broad_c1;
        fp_broad.reserve(static_cast<std::size_t>(broad_fp_count));
        fp_broad_c0.reserve(static_cast<std::size_t>(broad_fp_count));
        fp_broad_c1.reserve(static_cast<std::size_t>(broad_fp_count));

        std::uint64_t fp_count = 0;
        std::uint64_t fn_count = 0;
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
            fp_count += fp[i];
            fn_count += fn[i];

            const bool broad_found = contains_pair(broad_actual, pair_key(c0[i], c1[i]));
            if (!broad_found) {
                fn_broad[i] = 1;
                fn_broad_c0.push_back(c0[i]);
                fn_broad_c1.push_back(c1[i]);
            }
        }

        if (broad_fp_count != 0) {
            for (const auto pair : broad_actual) {
                if (!contains_pair(broad_expected, pair)) {
                    fp_broad.push_back(1);
                    fp_broad_c0.push_back(static_cast<std::int32_t>(pair >> 32));
                    fp_broad_c1.push_back(static_cast<std::int32_t>(pair & 0xffffffffu));
                }
            }
        }

        const std::uint64_t broad_fn_count = static_cast<std::uint64_t>(fn_broad_c0.size());
        std::cout << dataset << ',' << case_file.key << ',' << (case_file.is_vf ? "vf" : "ee") << ',' << q0.size()
                  << ',' << broad_ms << ',' << narrow_ms << ',' << fp_count << ',' << fn_count << ','
                  << broad_fp_count << ',' << broad_fn_count << '\n';

        const bool wrote =
            write_raw(roots_dir / "sccd_toi.float64", sccd_toi) && write_raw(roots_dir / "sccd_fp.uint8", fp) &&
            write_raw(roots_dir / "sccd_fn.uint8", fn) && write_raw(roots_dir / "sccd_fp_broad.uint8", fp_broad) &&
            write_raw(roots_dir / "sccd_fn_broad.uint8", fn_broad) &&
            write_raw(roots_dir / "sccd_fp_broad_c0.int32", fp_broad_c0) &&
            write_raw(roots_dir / "sccd_fp_broad_c1.int32", fp_broad_c1) &&
            write_raw(roots_dir / "sccd_fn_broad_c0.int32", fn_broad_c0) &&
            write_raw(roots_dir / "sccd_fn_broad_c1.int32", fn_broad_c1);
        if (fn_count != 0 || broad_fn_count != 0) {
            std::cerr << "error: false negatives in " << dataset << "/" << case_file.key << ": fn=" << fn_count
                      << " broad_fn=" << broad_fn_count << "\n";
            return false;
        }
        return wrote;
    }

}  // namespace

int main(int argc, char** argv) {
    auto ctx = smesh::initialize(argc, argv);
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <data-dir> <dataset> [<dataset> ...]\n";
        return EXIT_FAILURE;
    }

    const fs::path data_dir = argv[1];
    bool ok = true;
    std::cout << "dataset,case,type,queries,broad_ms,narrow_ms,fp,fn,broad_fp,broad_fn\n";
    for (int i = 2; i < argc; ++i) {
        const std::string dataset = argv[i];
        const fs::path dataset_dir = data_dir / dataset;
        std::vector<double> vf_timings;
        std::vector<double> ee_timings;

        for (const CaseFile& case_file : scan_cases(dataset_dir / "boxes")) {
            ok = run_case(
                     ctx->communicator(), dataset_dir, dataset, case_file, case_file.is_vf ? vf_timings : ee_timings) &&
                 ok;
        }

        const fs::path out_dir = dataset_dir / "benchmark";
        ok = write_raw(out_dir / (dataset + "-vf.float64"), vf_timings) && ok;
        ok = write_raw(out_dir / (dataset + "-ee.float64"), ee_timings) && ok;
    }

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
