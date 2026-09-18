// Additive CCD over the SCCD benchmark.
//
// Same scenes, same cases, same CSV columns as benchmark/bench.exe.cpp, so a run
// of this binary concatenates with a run of that one.
//
//     SCCD_BENCH_CASE_BEGIN=0 SCCD_BENCH_CASE_END=200 \
//         accd_bench <data-dir> <scene> [<scene> ...]
//
//     --header   print the CSV schema and exit
//
// **Additive CCD is a narrow phase, compared with SCCD per collision pair.**
// Given a pair of primitives it advances conservatively until they touch; it
// produces no candidates of its own. Two measurements come out of each case:
//
// * **Cost**, over SCCD's broad-phase candidates on the mesh. `broad_ms` is SCCD's
//   broad phase and `narrow_ms` is additive CCD over what it produced, so `ns/q`
//   sits beside SCCD's per-pair narrow phase (`narrow_ms_s1`) over the identical
//   list. Additive CCD answers one pair at a time and carries no parallelism of
//   its own, so the loop over candidates is a `tbb::parallel_for`, which is how
//   `Candidates::compute_collision_free_stepsize` drives it; against SCCD's
//   threaded narrow phase a serial loop would be timing one thread against 72.
// * **Accuracy**, per curated query at the query's own coordinates -- the ones
//   the exact roots were computed for -- which is where `sccd_bench` scores SCCD.
//   `query_narrow_ms` is the cost of that pass.
//
// **Format conversion is outside both measurements.** ACCD takes eight points per
// query and SCCD's broad phase returns pairs of element indices, so the gather
// from indices to coordinates sits between the two timed regions and is charged
// to neither. It is an artefact of bolting two libraries together, not a cost
// either of them would pay in use.
//
// **It is conservative by construction and deliberately not tight.** The method
// stops short of contact by design -- `conservative_rescaling` is the knob, 0.9
// by default -- so `toi_late` must be zero and the interesting columns are
// `toi_max_early` and `toi_med_early`: how much earlier than the true root it
// stops. A non-zero `toi_late` means the harness is wrong, because the method
// cannot overshoot.

#include "competitor_bench.hpp"

#include "sccd_base.hpp"
#include "sccd_broadphase_strategy.hpp"
#include "sccd_smesh_ccd.hpp"

#include "smesh_buffer.hpp"
#include "smesh_graph.hpp"
#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"

#include <Eigen/Core>

#include <ipc/ccd/additive_ccd.hpp>
#include <ipc/utils/logger.hpp>

#include <spdlog/sinks/stdout_sinks.h>

#include <tbb/parallel_for.h>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

    using namespace sccd_competitor;

    // The type SCCD's broad phase computes in, and it has to be the one
    // `bench.exe.cpp` uses: the mesh stores float32, both harnesses widen it,
    // and widening to float instead of double leaves the outward-ULP box
    // rounding coarser. That produced 64 more candidates on armadillo 100ee than
    // the SCCD rows saw, so the two narrow phases were not scoring the same list.
    using scalar_t = double;
    using Clock = std::chrono::steady_clock;

    double ms_since(const Clock::time_point t) {
        return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
    }

    // The method's own parameters, as the IPC toolkit ships them.
    //
    // `conservative_rescaling` is what makes additive CCD stop short of contact,
    // and 0.9 is the value the toolkit and [Li et al. 2021] use. It is the knob
    // whose effect this comparison reports, so it stays at the default rather
    // than being tuned to flatter either side. Iterations are unlimited, for the
    // same reason the Scalable CCD harness gives: a budget that truncates the
    // search measures the budget.
    constexpr long kMaxIterations = ipc::AdditiveCCD::UNLIMITTED_ITERATIONS;
    constexpr double kConservativeRescaling = ipc::AdditiveCCD::DEFAULT_CONSERVATIVE_RESCALING;
    constexpr double kTMax = 1.0;

    /// Minimum separation, the ξ of the method. Zero by default.
    ///
    /// Zero asks for exact continuous collision detection, which is the query
    /// SCCD and Scalable CCD answer and so the only setting that compares. A
    /// positive ξ is what IPC uses and answers a different question, so a number
    /// from it must be labelled with the ξ it used.
    ///
    /// Zero is well behaved here. It is worth recording that it briefly looked
    /// otherwise: fed the wrong edge endpoints, additive CCD saw coincident
    /// points, reported `toi = 0`, warned that its gap was below machine epsilon
    /// and "can lead to missed collisions", and spent ten million iterations on a
    /// single pair -- some candidates taking minutes. All three warnings were the
    /// toolkit correctly reporting that its precondition `d > min_distance` had
    /// been violated by its caller. With the numbering fixed they do not appear.
    double min_distance() {
        if (const char* v = std::getenv("SCCD_ACCD_MIN_DISTANCE")) {
            const double d = std::atof(v);
            if (d > 0) return d;
        }
        return 0.0;
    }

    /// Cap on candidates scored per case, 0 for all and the default.
    ///
    /// Kept as a diagnostic. It was what localised a pathology to a particular
    /// candidate index by bisection, and it is the tool to reach for again if one
    /// case ever runs long. `narrow_ms` and `queries` stay an honest pair under a
    /// cap -- the measured cost of the measured number of candidates, so `ns/q`
    /// remains a true rate -- but the per-case total is then not the cost of the
    /// case, and any report drawn from a capped run has to say so.
    std::size_t candidate_cap() {
        if (const char* v = std::getenv("SCCD_ACCD_MAX_CANDIDATES")) {
            const long long n = std::atoll(v);
            if (n > 0) return static_cast<std::size_t>(n);
        }
        return 0;
    }

    struct MeshPair {
        std::shared_ptr<smesh::Mesh> t0;
        std::shared_ptr<smesh::Mesh> t1;
    };

    bool load_meshes(const std::shared_ptr<smesh::Communicator>& comm, const fs::path& dataset_dir,
                     const std::string& dataset, const int step, MeshPair& out) {
        const fs::path d0 = raw_frame_path(dataset_dir, dataset, step);
        const fs::path d1 = raw_frame_path(dataset_dir, dataset, step + 1);
        if (!fs::is_directory(d0) || !fs::is_directory(d1)) {
            std::cerr << "error: no prepared mesh for " << dataset << " step " << step
                      << "; expected " << d0 << " -- run benchmark/scripts/prepare_data.sh\n";
            return false;
        }
        out.t0 = smesh::Mesh::create_from_file(comm, smesh::Path(d0.string()));
        out.t1 = smesh::Mesh::create_from_file(comm, smesh::Path(d1.string()));
        return out.t0 && out.t1;
    }

    /// Eight points per candidate: the four primitives at t=0, then at t=1.
    ///
    /// Materialised between the two timed regions. Gathering indices into
    /// coordinates is what the seam between the two libraries costs, and charging
    /// it to either one would misreport that library.
    struct Candidates {
        std::vector<Eigen::Vector3d> p;  // 8 per candidate
        std::size_t n = 0;

        void resize(const std::size_t count) {
            n = count;
            p.assign(count * 8, Eigen::Vector3d::Zero());
        }
        Eigen::Vector3d* at(const std::size_t i) { return p.data() + i * 8; }
        const Eigen::Vector3d* at(const std::size_t i) const { return p.data() + i * 8; }
    };

    /// The edges as SCCD numbers them, and the map onto the benchmark's numbering.
    ///
    /// SCCD's broad phase returns indices into its own edge list, which comes
    /// from the mesh's edge graph, while the curated queries are in the
    /// benchmark's ordering. `bench.exe.cpp` reconciles the two with
    /// `make_benchmark_edge_id_map`; the same has to happen here or every curated
    /// pair reads as missed and, worse, the coordinate gather picks up the wrong
    /// edges entirely.
    struct EdgeTables {
        std::vector<smesh::idx_t> a, b;      ///< endpoints, in SCCD's order
        std::vector<idx_t> to_benchmark;     ///< SCCD's id -> the benchmark's
    };

    bool build_edge_tables(const std::shared_ptr<smesh::Mesh>& mesh, EdgeTables& out) {
        const auto ordered = benchmark_ordered_edges(mesh);
        std::unordered_map<std::uint64_t, idx_t> benchmark_ids;
        benchmark_ids.reserve(ordered[0].size() * 2 + 1);
        for (std::size_t i = 0; i < ordered[0].size(); ++i) {
            benchmark_ids.emplace(pair_key(ordered[0][i], ordered[1][i]), static_cast<idx_t>(i));
        }

        auto graph = mesh->edge_graph();
        auto row_idx = smesh::create_host_buffer<smesh::idx_t>(graph->nnz());
        smesh::crs_to_coo(mesh->n_nodes(), graph->rowptr()->data(), row_idx->data());

        const ptrdiff_t nnz = graph->nnz();
        out.a.resize(static_cast<std::size_t>(nnz));
        out.b.resize(static_cast<std::size_t>(nnz));
        out.to_benchmark.resize(static_cast<std::size_t>(nnz));
        for (ptrdiff_t i = 0; i < nnz; ++i) {
            const smesh::idx_t x = row_idx->data()[i];
            const smesh::idx_t y = graph->colidx()->data()[i];
            out.a[static_cast<std::size_t>(i)] = x;
            out.b[static_cast<std::size_t>(i)] = y;
            const auto it = benchmark_ids.find(pair_key(std::min(x, y), std::max(x, y)));
            if (it == benchmark_ids.end()) {
                std::cerr << "error: CCD edge " << i << " (" << x << ", " << y
                          << ") has no benchmark numbering\n";
                return false;
            }
            out.to_benchmark[static_cast<std::size_t>(i)] = it->second;
        }
        return true;
    }

    Eigen::Vector3d point_of(scalar_t* const* pts, const smesh::idx_t v) {
        return Eigen::Vector3d(static_cast<double>(pts[0][v]), static_cast<double>(pts[1][v]),
                               static_cast<double>(pts[2][v]));
    }

}  // namespace

int main(int argc, char** argv) {
    auto ctx = smesh::initialize(argc, argv);

    // The toolkit's logger writes to stdout by default, which is where this
    // harness writes its CSV: a warning about one query then becomes a row that
    // parses as a scene named "[timestamp] [ipctk] [warning] ...". The warnings
    // are worth keeping -- they say the method was handed a pair that already
    // touches -- so they are moved to stderr rather than silenced.
    ipc::set_logger(std::make_shared<spdlog::logger>(
        "ipctk", std::make_shared<spdlog::sinks::stderr_sink_mt>()));

    if (argc == 2 && std::string(argv[1]) == "--header") {
        std::cout << kCsvHeader << "\n";
        return EXIT_SUCCESS;
    }
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <data-dir> <scene> [<scene> ...]\n";
        return EXIT_FAILURE;
    }

    const fs::path data_dir(argv[1]);
    const ipc::AdditiveCCD accd(kMaxIterations, kConservativeRescaling);

    std::cout << kCsvHeader << "\n";
    std::cout.precision(6);

    int failures = 0;
    for (int a = 2; a < argc; ++a) {
        const std::string scene = argv[a];
        const fs::path dataset_dir = data_dir / scene;

        auto cases = scan_cases(dataset_dir / "boxes");
        if (cases.empty()) {
            std::cerr << "error: no prepared cases under " << (dataset_dir / "boxes") << "\n";
            ++failures;
            continue;
        }
        apply_case_range(cases);

        for (const CaseFile& c : cases) {
            const int step = case_step(c.key);
            if (step < 0) continue;

            CaseTruth truth;
            if (!load_case_truth(dataset_dir, c, truth)) {
                ++failures;
                continue;
            }

            MeshPair meshes;
            if (!load_meshes(ctx->communicator(), dataset_dir, scene, step, meshes)) {
                ++failures;
                continue;
            }

            // --- SCCD's broad phase, set up exactly as bench.exe.cpp does -----
            auto ccd = sccd::CCD<scalar_t>::create(meshes.t0, smesh::EXECUTION_SPACE_HOST);
            ccd->set_box_rounding(sccd::BoxRounding::OutwardUlp);
            auto points0 = smesh::astype<scalar_t>(meshes.t0->points());
            auto points1 = smesh::astype<scalar_t>(meshes.t1->points());
            // The edge graph and the working buffers are a function of the mesh,
            // so a simulation pays for them once. Without this they would land
            // inside the first timed region of every case.
            ccd->initialize();

            Row row;
            row.dataset = scene;
            row.mode = "accd";
            row.broadphase = sccd::broadphase_strategy_name(sccd::broadphase_strategy_setting());
            row.key = c.key;
            row.type = c.is_vf ? "vf" : "ee";
            row.gt_earliest = truth.gt_earliest;
            row.root_n = truth.root_toi.size();

            EdgeTables et;
            if (!build_edge_tables(meshes.t0, et)) {
                ++failures;
                continue;
            }

            if (!normalize_pairs(truth, c.is_vf, meshes.t0->n_nodes(),
                                 static_cast<ptrdiff_t>(
                                     benchmark_ordered_edges(meshes.t0)[0].size()))) {
                ++failures;
                continue;
            }

            smesh::SharedBuffer<smesh::idx_t> a0, a1;
            {
                const auto t = Clock::now();
                int err = ccd->broad_phase_prep(points0, points1);
                row.prep_ms = ms_since(t);
                const auto t2 = Clock::now();
                err |= c.is_vf ? ccd->broad_phase_fv_step(a0, a1)
                               : ccd->broad_phase_ee_step(a0, a1);
                row.broad_ms = ms_since(t2);
                if (err != SCCD_SUCCESS) {
                    std::cerr << "error: broad phase failed on " << scene << "/" << c.key << "\n";
                    ++failures;
                    continue;
                }
            }

            // --- the seam: indices to coordinates, charged to neither phase ---
            std::size_t n = a0 ? static_cast<std::size_t>(a0->size()) : 0;
            const std::size_t cap = candidate_cap();
            const std::size_t produced = n;
            if (cap && n > cap) n = cap;
            BroadAccounting acct(truth);
            Candidates cand;
            cand.resize(n);
            std::size_t skipped_adjacent = 0;
            std::size_t kept = 0;
            std::vector<std::size_t> keep;
            keep.reserve(n);
            {
                auto faces = meshes.t0->block(0)->elements();
                scalar_t* const* p0 = points0->data();
                scalar_t* const* p1 = points1->data();

                for (std::size_t i = 0; i < n; ++i) {
                    const smesh::idx_t x = a0->data()[i];
                    const smesh::idx_t y = a1->data()[i];
                    smesh::idx_t v[4];
                    if (c.is_vf) {
                        // (vertex, face): the point, then the triangle. Vertex and
                        // face numbering is the mesh's, so it needs no mapping.
                        acct.add(static_cast<int>(x), static_cast<int>(y), row.broad_fp);
                        v[0] = x;
                        for (int k = 0; k < 3; ++k) v[k + 1] = faces->data()[k][y];
                    } else {
                        // (edge, edge): endpoints from SCCD's own numbering, and
                        // the pair recorded in the benchmark's.
                        acct.add(static_cast<int>(et.to_benchmark[static_cast<std::size_t>(x)]),
                                 static_cast<int>(et.to_benchmark[static_cast<std::size_t>(y)]),
                                 row.broad_fp);
                        v[0] = et.a[static_cast<std::size_t>(x)];
                        v[1] = et.b[static_cast<std::size_t>(x)];
                        v[2] = et.a[static_cast<std::size_t>(y)];
                        v[3] = et.b[static_cast<std::size_t>(y)];
                    }
                    // Additive CCD opens with `assert(d > min_distance)`, so a
                    // pair of primitives sharing a vertex violates its
                    // precondition. SCCD's broad phase already excludes adjacency,
                    // so this never fires on its output and the count is zero on
                    // every case measured -- it is a guard against handing the
                    // method something it documents as out of contract, not a
                    // filter the comparison depends on. Outside the timed region,
                    // and counted so that if it ever does fire it is visible.
                    bool adjacent = false;
                    for (int j = 0; j < 2 && !adjacent; ++j) {
                        for (int k = 2; k < 4; ++k) {
                            if (v[j] == v[k]) { adjacent = true; break; }
                        }
                    }
                    if (c.is_vf) {
                        adjacent = (v[0] == v[1] || v[0] == v[2] || v[0] == v[3]);
                    }
                    if (adjacent) {
                        ++skipped_adjacent;
                        continue;
                    }

                    keep.push_back(i);
                    Eigen::Vector3d* q = cand.at(kept);
                    for (int k = 0; k < 4; ++k) {
                        q[k] = point_of(p0, v[k]);
                        q[k + 4] = point_of(p1, v[k]);
                    }
                    ++kept;
                }
            }
            row.queries = kept;
            // Accounting is over everything the broad phase produced, because
            // that is a property of the broad phase and not of how much of its
            // output was scored.
            for (std::size_t i = n; i < produced; ++i) {
                const smesh::idx_t x = a0->data()[i];
                const smesh::idx_t y = a1->data()[i];
                if (c.is_vf) {
                    acct.add(static_cast<int>(x), static_cast<int>(y), row.broad_fp);
                } else {
                    acct.add(static_cast<int>(et.to_benchmark[static_cast<std::size_t>(x)]),
                             static_cast<int>(et.to_benchmark[static_cast<std::size_t>(y)]),
                             row.broad_fp);
                }
            }
            row.broad_fn = acct.missing();

            // --- additive CCD over exactly those candidates -------------------
            std::vector<double> cand_toi(kept, 1.0);
            {
                const double ms = min_distance();
                const auto t = Clock::now();
                // Over the candidates in parallel, which is how the toolkit
                // drives this narrow phase: `Candidates::compute_collision_free_stepsize`
                // wraps the per-candidate call in exactly this `tbb::parallel_for`.
                // A serial loop here would time one thread against SCCD's 72.
                tbb::parallel_for(std::size_t(0), kept, [&](const std::size_t i) {
                    const Eigen::Vector3d* q = cand.at(i);
                    double toi = 1.0;
                    const bool hit =
                        c.is_vf ? accd.point_triangle_ccd(q[0], q[1], q[2], q[3], q[4], q[5], q[6],
                                                          q[7], toi, ms, kTMax)
                                : accd.edge_edge_ccd(q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7],
                                                     toi, ms, kTMax);
                    cand_toi[i] = hit ? toi : 1.0;
                });
                row.narrow_ms = ms_since(t);
            }

            if (skipped_adjacent) {
                std::cerr << "note: " << scene << "/" << c.key << ": skipped "
                          << skipped_adjacent << " of " << produced
                          << " candidates whose primitives share a vertex; additive CCD "
                             "requires a positive initial distance\n";
            }
            double earliest = 1.0;
            for (const double t : cand_toi) earliest = std::min(earliest, t);
            row.earliest_toi = earliest;

            // --- per pair: the curated queries at their own coordinates -------
            //
            // Scored exactly where `sccd_bench` scores SCCD: each curated query
            // at the coordinates its exact root was computed for, one time of
            // impact per query. `query_narrow_ms` is that pass, as it is there.
            std::vector<double> xyz;
            if (!read_query_points(dataset_dir, c.key, xyz)) {
                ++failures;
                continue;
            }
            const std::size_t nq = xyz.size() / 24;
            if (nq != truth.size()) {
                std::cerr << "error: " << scene << "/" << c.key << ": " << nq
                          << " query geometries against " << truth.size() << " curated pairs\n";
                ++failures;
                continue;
            }
            std::vector<double> toi(nq, 1.0);
            {
                const double ms = min_distance();
                const auto t = Clock::now();
                tbb::parallel_for(std::size_t(0), nq, [&](const std::size_t i) {
                    Eigen::Vector3d q[8];
                    for (int k = 0; k < 8; ++k) {
                        q[k] = Eigen::Vector3d(xyz[24 * i + 3 * k], xyz[24 * i + 3 * k + 1],
                                               xyz[24 * i + 3 * k + 2]);
                    }
                    double t_i = 1.0;
                    const bool hit =
                        c.is_vf ? accd.point_triangle_ccd(q[0], q[1], q[2], q[3], q[4], q[5], q[6],
                                                          q[7], t_i, ms, kTMax)
                                : accd.edge_edge_ccd(q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7],
                                                     t_i, ms, kTMax);
                    toi[i] = hit ? t_i : 1.0;
                });
                row.query_narrow_ms = ms_since(t);
            }
            row.acc = score(truth, toi);

            // The method cannot overshoot, so a late time of impact is a harness
            // fault, not a measurement. Say so rather than publishing it.
            if (row.acc.toi_late != 0) {
                std::cerr << "error: " << scene << "/" << c.key << ": additive CCD reported "
                          << row.acc.toi_late << " time(s) of impact after the true root, by up to "
                          << row.acc.toi_max_late
                          << ". The method advances conservatively and cannot do this, so the "
                             "harness is reading the candidates wrongly.\n";
                ++failures;
            }

            write_row(std::cout, row);
        }
    }

    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
