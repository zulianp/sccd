#include "sccd_smesh_CCD.hpp"
#include "smesh_buffer.hpp"
#include "smesh_context.hpp"
#include "smesh_device_buffer.hpp"
#include "smesh_graph.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"
#include "smesh_tracer.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fs = std::filesystem;

namespace {

    using scalar_t = double;
    using idx_t = smesh::idx_t;

    int narrowphase_max_depth = 69;
    scalar_t narrowphase_tol = scalar_t(3e-8);

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

    struct BroadphaseResult {
        int err = SCCD_SUCCESS;
        double prep_elapsed_ms = 0.0;
        double elapsed_ms = 0.0;
        std::unordered_set<std::uint64_t> pairs;
        std::uint64_t false_positives = 0;
        smesh::SharedBuffer<idx_t> v_overlap;
        smesh::SharedBuffer<idx_t> f_overlap;
        smesh::SharedBuffer<idx_t> e0_overlap;
        smesh::SharedBuffer<idx_t> e1_overlap;
    };

    struct CCDRun {
        std::shared_ptr<sccd::CCD<scalar_t>> ccd;
        smesh::SharedBuffer<scalar_t*> points0;
        smesh::SharedBuffer<scalar_t*> points1;
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

        const fs::path queries_dir = boxes_dir.parent_path() / "queries";
        for (const auto& entry : fs::directory_iterator(boxes_dir)) {
            if (!entry.is_directory()) {
                continue;
            }
            const fs::path dir = entry.path();
            const std::string key = dir.filename().string();
            if ((key.size() < 3) || !fs::exists(dir / "c0.int32") || !fs::exists(dir / "c1.int32")) {
                continue;
            }
            if (!fs::exists(queries_dir / (key + ".csv"))) {
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

    std::array<std::vector<idx_t>, 2> benchmark_ordered_edges(const std::shared_ptr<smesh::Mesh>& mesh) {
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

    bool make_benchmark_edge_id_map(const std::shared_ptr<smesh::Mesh>& mesh, std::vector<idx_t>& edge_id_map) {
        const auto edges = benchmark_ordered_edges(mesh);
        std::unordered_map<std::uint64_t, idx_t> benchmark_ids;
        benchmark_ids.reserve(edges[0].size() * 2 + 1);
        for (std::size_t i = 0; i < edges[0].size(); ++i) {
            benchmark_ids.emplace(pair_key(edges[0][i], edges[1][i]), static_cast<idx_t>(i));
        }

        auto graph = mesh->edge_graph();
        auto row_idx = smesh::create_host_buffer<idx_t>(graph->nnz());
        smesh::crs_to_coo(mesh->n_nodes(), graph->rowptr()->data(), row_idx->data());

        edge_id_map.resize(static_cast<std::size_t>(graph->nnz()));
        for (ptrdiff_t i = 0; i < graph->nnz(); ++i) {
            const idx_t a = std::min(row_idx->data()[i], graph->colidx()->data()[i]);
            const idx_t b = std::max(row_idx->data()[i], graph->colidx()->data()[i]);
            const auto it = benchmark_ids.find(pair_key(a, b));
            if (it == benchmark_ids.end()) {
                std::cerr << "error: failed to map CCD edge " << i << " (" << a << ", " << b
                          << ") to benchmark edge numbering\n";
                return false;
            }
            edge_id_map[static_cast<std::size_t>(i)] = it->second;
        }
        return true;
    }

    smesh::ExecutionSpace benchmark_execution_space() {
        const char* env = std::getenv("SCCD_BENCH_EXECUTION_SPACE");
        if (!env || env[0] == '\0') {
            env = std::getenv("SCCD_EXECUTION_SPACE");
        }
        if (!env || env[0] == '\0') {
            return smesh::EXECUTION_SPACE_HOST;
        }

        std::string value(env);
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        });

        if (value == "cuda" || value == "gpu" || value == "device") {
#if defined(SCCD_ENABLE_CUDA)
            return smesh::EXECUTION_SPACE_DEVICE;
#else
            static bool warned = false;
            if (!warned) {
                std::cerr << "warning: SCCD_BENCH_EXECUTION_SPACE=" << env
                          << " requested but SCCD was built without CUDA; using host\n";
                warned = true;
            }
#endif
        }
        return smesh::EXECUTION_SPACE_HOST;
    }

    CCDRun make_ccd_run(const MeshPair& meshes, const smesh::ExecutionSpace execution_space) {
        CCDRun run;
        run.ccd = sccd::CCD<scalar_t>::create(meshes.t0, execution_space);
        run.points0 = smesh::astype<scalar_t>(meshes.t0->points());
        run.points1 = smesh::astype<scalar_t>(meshes.t1->points());
        if (execution_space == smesh::EXECUTION_SPACE_DEVICE) {
            run.points0 = smesh::to_device(run.points0);
            run.points1 = smesh::to_device(run.points1);
        }
        return run;
    }

    BroadphaseResult run_broadphase(const bool is_vf, CCDRun& ccd_run) {
        SMESH_TRACE_SCOPE("benchmark broadphase");

        BroadphaseResult result;
        auto start = std::chrono::steady_clock::now();
        result.err = ccd_run.ccd->broad_phase_prep(ccd_run.points0, ccd_run.points1);
        auto stop = std::chrono::steady_clock::now();
        result.prep_elapsed_ms = std::chrono::duration<double, std::milli>(stop - start).count();
        if (result.err != SCCD_SUCCESS) {
            return result;
        }

        start = std::chrono::steady_clock::now();
        if (is_vf) {
            result.err = ccd_run.ccd->broad_phase_fv_step(result.v_overlap, result.f_overlap);
        } else {
            result.err = ccd_run.ccd->broad_phase_ee_step(result.e0_overlap, result.e1_overlap);
        }
        stop = std::chrono::steady_clock::now();
        result.elapsed_ms = std::chrono::duration<double, std::milli>(stop - start).count();

        return result;
    }

    void populate_broadphase_pairs(const MeshPair& meshes,
                                   const bool is_vf,
                                   const std::unordered_set<std::uint64_t>& expected,
                                   const std::vector<idx_t>& edge_id_map,
                                   BroadphaseResult& result) {
        const ptrdiff_t n_nodes = meshes.t0->n_nodes();
        const ptrdiff_t n_edges =
            edge_id_map.empty() ? meshes.t0->edge_graph()->nnz() : static_cast<ptrdiff_t>(edge_id_map.size());
        const ptrdiff_t face_offset = n_nodes + n_edges;
        const ptrdiff_t edge_offset = n_nodes;

        if (is_vf) {
            auto v_overlap = smesh::to_host(result.v_overlap);
            auto f_overlap = smesh::to_host(result.f_overlap);
            result.pairs.reserve(v_overlap->size() * 2 + 1);
            for (std::size_t i = 0; i < v_overlap->size(); ++i) {
                result.pairs.insert(pair_key(v_overlap->data()[i], f_overlap->data()[i] + face_offset));
            }
        } else {
            auto e0_overlap = smesh::to_host(result.e0_overlap);
            auto e1_overlap = smesh::to_host(result.e1_overlap);
            result.pairs.reserve(e0_overlap->size() * 2 + 1);
            for (std::size_t i = 0; i < e0_overlap->size(); ++i) {
                const idx_t e0 = edge_id_map.empty() ? e0_overlap->data()[i] : edge_id_map[e0_overlap->data()[i]];
                const idx_t e1 = edge_id_map.empty() ? e1_overlap->data()[i] : edge_id_map[e1_overlap->data()[i]];
                result.pairs.insert(pair_key(e0 + edge_offset, e1 + edge_offset));
            }
        }

        result.false_positives = 0;
        for (const auto pair : result.pairs) {
            if (!contains_pair(expected, pair)) {
                ++result.false_positives;
            }
        }
    }

    double time_narrowphase_zero_stride(const bool is_vf, CCDRun& ccd_run, int& err) {
        scalar_t toi = scalar_t(1);
        smesh::SharedBuffer<scalar_t> vf_tois;
        smesh::SharedBuffer<scalar_t> ee_tois;

        const auto start = std::chrono::steady_clock::now();
        if (is_vf) {
            err = ccd_run.ccd->narrow_phase_fv(toi, vf_tois, narrowphase_max_depth, narrowphase_tol, 0);
        } else {
            err = ccd_run.ccd->narrow_phase_ee(toi, ee_tois, narrowphase_max_depth, narrowphase_tol, 0);
        }
        const auto stop = std::chrono::steady_clock::now();
        static volatile scalar_t toi_sink;
        toi_sink = toi;
        return std::chrono::duration<double, std::milli>(stop - start).count();
    }

    template <typename T>
    smesh::SharedBuffer<T> vector_buffer(std::vector<T>& values) {
        return smesh::Buffer<T>::wrap(values.size(), values.data(), smesh::MEMORY_SPACE_HOST);
    }

    template <typename T, std::size_t N>
    smesh::SharedBuffer<T*> make_2d_buffer(std::array<std::vector<T>, N>& values,
                                           const smesh::ExecutionSpace execution_space) {
        std::vector<smesh::SharedBuffer<T>> buffers;
        buffers.reserve(N);
        for (std::vector<T>& value : values) {
            buffers.push_back(vector_buffer(value));
        }

        if (execution_space == smesh::EXECUTION_SPACE_DEVICE) {
            buffers = smesh::to_device(buffers);
        }
        return smesh::create_2d(buffers);
    }

    bool compute_query_toi(const bool is_vf,
                           QueryGeometry& query_geometry,
                           const smesh::ExecutionSpace execution_space,
                           std::vector<scalar_t>& sccd_toi) {
        sccd_toi.assign(query_geometry.q0.size(), scalar_t(1));

        if (execution_space == smesh::EXECUTION_SPACE_DEVICE) {
#if defined(SCCD_ENABLE_CUDA)
            smesh::SharedBuffer<idx_t> q0 = smesh::to_device(vector_buffer(query_geometry.q0));
            smesh::SharedBuffer<idx_t> q1 = smesh::to_device(vector_buffer(query_geometry.q1));
            smesh::SharedBuffer<scalar_t*> points0 = make_2d_buffer(query_geometry.points0, execution_space);
            smesh::SharedBuffer<scalar_t*> points1 = make_2d_buffer(query_geometry.points1, execution_space);
            smesh::SharedBuffer<scalar_t> toi = smesh::to_device(vector_buffer(sccd_toi));

            if (is_vf) {
                smesh::SharedBuffer<idx_t*> faces = make_2d_buffer(query_geometry.faces, execution_space);
                sccd::device::narrow_phase_vf<3>(query_geometry.q0.size(),
                                                 q0->data(),
                                                 q1->data(),
                                                 points0->data(),
                                                 points1->data(),
                                                 1,
                                                 faces->data(),
                                                 scalar_t(1),
                                                 toi->data(),
                                                 narrowphase_max_depth,
                                                 narrowphase_tol,
                                                 1);
            } else {
                smesh::SharedBuffer<idx_t*> edges = make_2d_buffer(query_geometry.edges, execution_space);
                sccd::device::narrow_phase_ee(query_geometry.q0.size(),
                                              q0->data(),
                                              q1->data(),
                                              points0->data(),
                                              points1->data(),
                                              1,
                                              edges->data(),
                                              scalar_t(1),
                                              toi->data(),
                                              narrowphase_max_depth,
                                              narrowphase_tol,
                                              1);
            }

            auto host_toi = smesh::to_host(toi);
            sccd_toi.assign(host_toi->data(), host_toi->data() + host_toi->size());
            return true;
#else
            std::cerr << "error: CUDA query narrowphase requested but SCCD was built without CUDA\n";
            return false;
#endif
        }

        scalar_t* points0[3] = {
            query_geometry.points0[0].data(), query_geometry.points0[1].data(), query_geometry.points0[2].data()};
        scalar_t* points1[3] = {
            query_geometry.points1[0].data(), query_geometry.points1[1].data(), query_geometry.points1[2].data()};
        if (is_vf) {
            idx_t* faces[3] = {
                query_geometry.faces[0].data(), query_geometry.faces[1].data(), query_geometry.faces[2].data()};
            sccd::narrow_phase_vf<3, scalar_t, idx_t>(query_geometry.q0.size(),
                                                      query_geometry.q0.data(),
                                                      query_geometry.q1.data(),
                                                      points0,
                                                      points1,
                                                      1,
                                                      faces,
                                                      scalar_t(1),
                                                      sccd_toi.data(),
                                                      narrowphase_max_depth,
                                                      narrowphase_tol,
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
                                                   scalar_t(1),
                                                   sccd_toi.data(),
                                                   narrowphase_max_depth,
                                                   narrowphase_tol,
                                                   1);
        }
        return true;
    }

    bool write_case_outputs(const fs::path& dataset_dir,
                            const std::string& dataset,
                            const CaseFile& case_file,
                            const MeshPair& meshes,
                            const std::vector<std::int32_t>& c0,
                            const std::vector<std::int32_t>& c1,
                            const std::vector<idx_t>& q0,
                            QueryGeometry& query_geometry,
                            const std::vector<idx_t>& edge_id_map,
                            const std::unordered_set<std::uint64_t>& broad_expected,
                            const double prep_ms,
                            const double broad_ms,
                            const double narrow_ms,
                            const smesh::ExecutionSpace execution_space,
                            BroadphaseResult& broadphase,
                            std::vector<double>& timings_ms) {
        const std::size_t narrow_queries =
            case_file.is_vf ? broadphase.v_overlap->size() : broadphase.e0_overlap->size();
        timings_ms.push_back(narrow_ms);
        populate_broadphase_pairs(meshes, case_file.is_vf, broad_expected, edge_id_map, broadphase);

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

        std::vector<scalar_t> sccd_toi;
        if (!compute_query_toi(case_file.is_vf, query_geometry, execution_space, sccd_toi)) {
            return false;
        }

        std::vector<std::uint8_t> fp(q0.size(), 0);
        std::vector<std::uint8_t> fn(q0.size(), 0);
        std::vector<std::uint8_t> fp_broad;
        std::vector<std::uint8_t> fn_broad(q0.size(), 0);
        std::vector<std::int32_t> fp_broad_c0;
        std::vector<std::int32_t> fp_broad_c1;
        std::vector<std::int32_t> fn_broad_c0;
        std::vector<std::int32_t> fn_broad_c1;
        fp_broad.reserve(static_cast<std::size_t>(broadphase.false_positives));
        fp_broad_c0.reserve(static_cast<std::size_t>(broadphase.false_positives));
        fp_broad_c1.reserve(static_cast<std::size_t>(broadphase.false_positives));

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

            const bool broad_found = contains_pair(broadphase.pairs, pair_key(c0[i], c1[i]));
            if (!broad_found) {
                fn_broad[i] = 1;
                fn_broad_c0.push_back(c0[i]);
                fn_broad_c1.push_back(c1[i]);
            }
        }

        if (broadphase.false_positives != 0) {
            for (const auto pair : broadphase.pairs) {
                if (!contains_pair(broad_expected, pair)) {
                    fp_broad.push_back(1);
                    fp_broad_c0.push_back(static_cast<std::int32_t>(pair >> 32));
                    fp_broad_c1.push_back(static_cast<std::int32_t>(pair & 0xffffffffu));
                }
            }
        }

        const std::uint64_t broad_fn_count = static_cast<std::uint64_t>(fn_broad_c0.size());
        std::cout << dataset << ',' << case_file.key << ',' << (case_file.is_vf ? "vf" : "ee") << ',' << narrow_queries
                  << ',' << prep_ms << ',' << broad_ms << ',' << narrow_ms << ',' << fp_count << ',' << fn_count << ','
                  << broadphase.false_positives << ',' << broad_fn_count << '\n';

        const bool wrote = write_raw(roots_dir / "sccd_toi.float64", sccd_toi) &&
                           write_raw(roots_dir / "sccd_fp.uint8", fp) && write_raw(roots_dir / "sccd_fn.uint8", fn) &&
                           write_raw(roots_dir / "sccd_fp_broad.uint8", fp_broad) &&
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

        std::vector<idx_t> edge_id_map;
        if (!case_file.is_vf && !make_benchmark_edge_id_map(meshes.t0, edge_id_map)) {
            return false;
        }

        const smesh::ExecutionSpace execution_space = benchmark_execution_space();
        CCDRun ccd_run = make_ccd_run(meshes, execution_space);
        const auto broad_expected = expected_pairs(c0, c1);
        const bool warmup_first_case = timings_ms.empty();
        if (warmup_first_case) {
            BroadphaseResult warmup = run_broadphase(case_file.is_vf, ccd_run);
            if (warmup.err != SCCD_SUCCESS) {
                std::cerr << "error: CCD broadphase warmup failed for " << dataset << "/" << case_file.key << "\n";
                return false;
            }
        }
        BroadphaseResult broadphase = run_broadphase(case_file.is_vf, ccd_run);
        if (broadphase.err != SCCD_SUCCESS) {
            std::cerr << "error: CCD broadphase failed for " << dataset << "/" << case_file.key << "\n";
            return false;
        }
        const double prep_ms = broadphase.prep_elapsed_ms;
        const double broad_ms = broadphase.elapsed_ms;
        int narrow_err = SCCD_SUCCESS;
        if (warmup_first_case) {
            time_narrowphase_zero_stride(case_file.is_vf, ccd_run, narrow_err);
            if (narrow_err != SCCD_SUCCESS) {
                std::cerr << "error: CCD narrowphase warmup failed for " << dataset << "/" << case_file.key << "\n";
                return false;
            }
        }
        narrow_err = SCCD_SUCCESS;
        const double narrow_ms = time_narrowphase_zero_stride(case_file.is_vf, ccd_run, narrow_err);
        if (narrow_err != SCCD_SUCCESS) {
            std::cerr << "error: CCD narrowphase failed for " << dataset << "/" << case_file.key << "\n";
            return false;
        }

        return write_case_outputs(dataset_dir,
                                  dataset,
                                  case_file,
                                  meshes,
                                  c0,
                                  c1,
                                  q0,
                                  query_geometry,
                                  edge_id_map,
                                  broad_expected,
                                  prep_ms,
                                  broad_ms,
                                  narrow_ms,
                                  execution_space,
                                  broadphase,
                                  timings_ms);
    }

}  // namespace

int main(int argc, char** argv) {
    auto ctx = smesh::initialize(argc, argv);
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <data-dir> <dataset> [<dataset> ...]\n";
        return EXIT_FAILURE;
    }

    int SCCD_MAX_DEPTH = narrowphase_max_depth;
    SCCD_READ_ENV(SCCD_MAX_DEPTH, atoi);
    narrowphase_max_depth = SCCD_MAX_DEPTH;

    scalar_t SCCD_TOL = narrowphase_tol;
    SCCD_READ_ENV(SCCD_TOL, atof);
    narrowphase_tol = SCCD_TOL;

    const fs::path data_dir = argv[1];
    bool ok = true;
    std::cout << "dataset,case,type,queries,prep_ms,broad_ms,narrow_ms,fp,fn,broad_fp,broad_fn\n";
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
