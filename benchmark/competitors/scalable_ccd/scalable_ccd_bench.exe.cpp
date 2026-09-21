// Scalable CCD over the SCCD benchmark.
//
// Same scenes, same cases, same phase split and same CSV columns as
// benchmark/bench.exe.cpp, so a run of this binary concatenates with a run of
// that one and benchmark/report/ reads both. What differs is only the library
// doing the work.
//
// Usage, and the environment it reads, match sccd_bench:
//
//     SCCD_BENCH_CASE_BEGIN=0 SCCD_BENCH_CASE_END=200 \
//         scalable_ccd_bench <data-dir> <scene> [<scene> ...]
//
//     --header                                   print the CSV schema and exit
//     --their-ccd    <data-dir> <scene> <step>   their `cuda::ccd` on one step
//     --buffer-probe <data-dir> <scene> <step> [<pad>]
//                                                the curated queries, bare and padded
//
// Three things about the library shape the harness.
//
// **Its narrow phase is CUDA only.** `scalable_ccd/broad_phase/` is host code;
// the root finding lives entirely under `scalable_ccd/cuda/narrow_phase/`. So a
// host build of this harness measures the broad phase and leaves the narrow
// columns empty, and only a CUDA build produces a full pipeline row.
//
// **It is compared with SCCD on the earliest time of impact.** Built as its
// authors ship it, the narrow phase is a branch and bound for the earliest
// contact: every box whose time lower bound is at or after the running earliest
// time of impact is pruned. `sccd_bench` times SCCD's earliest-time-of-impact
// path per case -- one query kind of one step, starting from a bound of 1 -- so
// this harness runs Scalable CCD's the same way: each case starts from a fresh
// bound, and `s0_toi` is that case's earliest time of impact. A step's answer is
// the minimum over its two cases, for both libraries alike.
//
// `SCALABLE_CCD_TOI_PER_QUERY` is a compile-time option of the library, off by
// default, that removes the cross-query pruning and reports every colliding
// candidate. It is not part of the comparison. CMake exposes it as
// `SCCD_SCALABLE_CCD_TOI_PER_QUERY` for diagnosis only; rows from such a build
// carry the mode `scalable-ccd-perquery`.
//
// **It interleaves the two phases.** `cuda::ccd` loops
// `while (!broad_phase.is_complete())`, alternating `detect_overlaps_partial`
// with `narrow_phase`, so that a scene whose overlap list does not fit in device
// memory is processed in batches. Timing therefore accumulates per phase across
// the loop.

#include "competitor_bench.hpp"

#include "sccd_base.hpp"

#include "smesh_context.hpp"
#include "smesh_mesh.hpp"
#include "smesh_path.hpp"

#include <Eigen/Core>

#include <scalable_ccd/config.hpp>
#include <scalable_ccd/broad_phase/aabb.hpp>
#include <scalable_ccd/broad_phase/sort_and_sweep.hpp>

#if defined(SCCD_COMPETITOR_WITH_CUDA)
#include <scalable_ccd/cuda/broad_phase/aabb.cuh>
#include <scalable_ccd/cuda/broad_phase/broad_phase.cuh>
#include <scalable_ccd/cuda/memory_handler.hpp>
#include <scalable_ccd/cuda/narrow_phase/narrow_phase.cuh>
#include <scalable_ccd/cuda/utils/device_matrix.cuh>
#include <scalable_ccd/cuda/ccd.cuh>
#include <cuda_runtime.h>
#include <thrust/copy.h>
#endif

#include <chrono>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

    using namespace sccd_competitor;

    using Clock = std::chrono::steady_clock;

#if defined(SCALABLE_CCD_TOI_PER_QUERY)
    constexpr bool kPerQueryBuild = true;
#else
    constexpr bool kPerQueryBuild = false;
#endif

    using Collisions = std::vector<std::tuple<int, int, scalable_ccd::Scalar>>;

    // Scalable CCD's own parameters, not SCCD's.
    //
    // `max_iterations` bounds the number of box checks a single query may make,
    // and exceeding it returns without recording anything: the box is dropped.
    // `-1` means unlimited and is what their own tests use; any finite value
    // makes the search non-conservative and measures a starved search.
    constexpr int kMaxIterations = -1;

    // Their tolerance is a co-domain tolerance, on the inclusion function's
    // value; SCCD's 3e-8 is a domain tolerance on (t, u, v). Different
    // quantities with no exact conversion, so each library runs at the accuracy
    // target its authors chose. 1e-6 is the value their tests use.
    constexpr double kDefaultTolerance = 1e-6;

    /// Overridable, so a dependence of the result on the tolerance can be tested
    /// rather than assumed.
    double tolerance() {
        if (const char* v = std::getenv("SCCD_COMPETITOR_TOLERANCE")) {
            const double t = std::atof(v);
            if (t > 0) return t;
        }
        return kDefaultTolerance;
    }

    // A plain continuous collision query: no minimum separation, and a contact
    // already present at t=0 is reported rather than skipped. Both are what
    // their narrow-phase test uses for the same question.
    constexpr double kMinDistance = 0.0;
    constexpr bool kAllowZeroToi = true;

    // The thread counts `cuda::ccd` uses.
    constexpr int kBroadPhaseThreads = 32;
    constexpr int kNarrowPhaseThreads = 1024;

    double ms_since(const Clock::time_point t) {
        return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
    }

    /// The mesh as Scalable CCD wants it: vertices rowwise, and the index arrays
    /// in the benchmark's own numbering so the pairs it returns need no mapping.
    struct MeshArrays {
        Eigen::MatrixXd V0;
        Eigen::MatrixXd V1;
        Eigen::MatrixXi E;
        Eigen::MatrixXi F;
    };

    bool load_mesh_arrays(const std::shared_ptr<smesh::Communicator>& comm, const fs::path& dataset_dir,
                          const std::string& dataset, const int step, MeshArrays& out) {
        const fs::path d0 = raw_frame_path(dataset_dir, dataset, step);
        const fs::path d1 = raw_frame_path(dataset_dir, dataset, step + 1);
        if (!fs::is_directory(d0) || !fs::is_directory(d1)) {
            std::cerr << "error: no prepared mesh for " << dataset << " step " << step
                      << "; expected " << d0 << " -- run benchmark/scripts/prepare_data.sh\n";
            return false;
        }

        auto m0 = smesh::Mesh::create_from_file(comm, smesh::Path(d0.string()));
        auto m1 = smesh::Mesh::create_from_file(comm, smesh::Path(d1.string()));
        if (!m0 || !m1) {
            std::cerr << "error: failed to read meshes for step " << step << "\n";
            return false;
        }
        if (m0->n_nodes() != m1->n_nodes()) {
            std::cerr << "error: frame " << step << " and " << step + 1 << " differ in vertex count\n";
            return false;
        }

        const ptrdiff_t n = m0->n_nodes();
        const ptrdiff_t nf = m0->block(0)->n_elements();

        // smesh keeps coordinates as structure-of-arrays and may store them as
        // float; Eigen wants rows of double. Widening is exact, so this only
        // changes the layout.
        out.V0.resize(n, 3);
        out.V1.resize(n, 3);
        auto p0 = m0->points();
        auto p1 = m1->points();
        for (int d = 0; d < 3; ++d) {
            for (ptrdiff_t i = 0; i < n; ++i) {
                out.V0(i, d) = static_cast<double>(p0->data()[d][i]);
                out.V1(i, d) = static_cast<double>(p1->data()[d][i]);
            }
        }

        auto faces = m0->block(0)->elements();
        out.F.resize(nf, 3);
        for (ptrdiff_t i = 0; i < nf; ++i) {
            for (int v = 0; v < 3; ++v) {
                out.F(i, v) = static_cast<int>(faces->data()[v][i]);
            }
        }

        const auto edges = benchmark_ordered_edges(m0);
        out.E.resize(static_cast<Eigen::Index>(edges[0].size()), 2);
        for (std::size_t i = 0; i < edges[0].size(); ++i) {
            out.E(static_cast<Eigen::Index>(i), 0) = static_cast<int>(edges[0][i]);
            out.E(static_cast<Eigen::Index>(i), 1) = static_cast<int>(edges[1][i]);
        }
        return true;
    }

#if defined(SCCD_COMPETITOR_WITH_CUDA)
    /// The mesh on the device. A format conversion, not work either library does
    /// in use: `sccd_bench` uploads its points before anything is timed, and so
    /// does this harness.
    struct DeviceMesh {
        scalable_ccd::cuda::DeviceMatrix<scalable_ccd::Scalar> V0;
        scalable_ccd::cuda::DeviceMatrix<scalable_ccd::Scalar> V1;
        scalable_ccd::cuda::DeviceMatrix<int> E;
        scalable_ccd::cuda::DeviceMatrix<int> F;

        explicit DeviceMesh(const MeshArrays& m) : V0(m.V0), V1(m.V1), E(m.E), F(m.F) {
            cudaDeviceSynchronize();
        }
    };

    /// One call of their narrow phase. The collision list only exists in the
    /// per-query build, which is the one difference between the two signatures.
    void call_narrow_phase(const DeviceMesh& d, const bool is_vf,
                           const thrust::device_vector<int2>& pairs,
                           const std::shared_ptr<scalable_ccd::cuda::MemoryHandler>& memory_handler,
                           Collisions& collisions, scalable_ccd::Scalar& toi) {
        namespace scc = scalable_ccd::cuda;
#if defined(SCALABLE_CCD_TOI_PER_QUERY)
        if (is_vf) {
            scc::narrow_phase<true>(d.V0, d.V1, d.E, d.F, pairs, kNarrowPhaseThreads, kMaxIterations,
                                    tolerance(), kMinDistance, kAllowZeroToi, memory_handler, collisions, toi);
        } else {
            scc::narrow_phase<false>(d.V0, d.V1, d.E, d.F, pairs, kNarrowPhaseThreads, kMaxIterations,
                                     tolerance(), kMinDistance, kAllowZeroToi, memory_handler, collisions, toi);
        }
#else
        (void)collisions;
        if (is_vf) {
            scc::narrow_phase<true>(d.V0, d.V1, d.E, d.F, pairs, kNarrowPhaseThreads, kMaxIterations,
                                    tolerance(), kMinDistance, kAllowZeroToi, memory_handler, toi);
        } else {
            scc::narrow_phase<false>(d.V0, d.V1, d.E, d.F, pairs, kNarrowPhaseThreads, kMaxIterations,
                                     tolerance(), kMinDistance, kAllowZeroToi, memory_handler, toi);
        }
#endif
    }

    /// One case on the device, composed the way `cuda::ccd` composes a pass,
    /// with each phase timed separately.
    ///
    /// `toi` is the running bound, in and out: it starts at 1 and comes back as
    /// the case's earliest time of impact.
    bool run_case_device(const MeshArrays& mesh, const DeviceMesh& d, const bool is_vf, Row& row,
                         BroadAccounting& acct, Collisions& collisions, scalable_ccd::Scalar& toi) {
        // Qualified rather than pulled in wholesale: `scalable_ccd::AABB` is the
        // host box and `scalable_ccd::cuda::AABB` the device one.
        namespace scc = scalable_ccd::cuda;

        auto memory_handler = std::make_shared<scc::MemoryHandler>();

        // --- prep: the acceleration structure ------------------------------
        //
        // Boxes, their upload and `BroadPhase::build`, which sorts, are the
        // structure and line up with `sccd_bench`'s `prep_ms`; the sweep below
        // is the query and lines up with `broad_ms`.
        const auto prep_begin = Clock::now();
        std::vector<scc::AABB> vertex_boxes, edge_boxes, face_boxes;
        scc::build_vertex_boxes(mesh.V0, mesh.V1, vertex_boxes, kMinDistance);
        scc::BroadPhase broad_phase(memory_handler);
        broad_phase.threads_per_block = kBroadPhaseThreads;
        if (is_vf) {
            scc::build_face_boxes(vertex_boxes, mesh.F, face_boxes);
            broad_phase.build(std::make_shared<scc::DeviceAABBs>(vertex_boxes),
                              std::make_shared<scc::DeviceAABBs>(face_boxes));
        } else {
            scc::build_edge_boxes(vertex_boxes, mesh.E, edge_boxes);
            broad_phase.build(std::make_shared<scc::DeviceAABBs>(edge_boxes));
        }
        cudaDeviceSynchronize();
        row.prep_ms += ms_since(prep_begin);

        // --- alternating batches -------------------------------------------
        while (!broad_phase.is_complete()) {
            const auto bp = Clock::now();
            broad_phase.detect_overlaps_partial();
            cudaDeviceSynchronize();
            row.broad_ms += ms_since(bp);

            // The candidates come back to the host to be counted against the
            // curated pairs, outside both timed regions.
            const thrust::device_vector<int2>& d_overlaps = broad_phase.overlaps();
            row.queries += d_overlaps.size();
            std::vector<int2> h_overlaps(d_overlaps.size());
            thrust::copy(d_overlaps.begin(), d_overlaps.end(), h_overlaps.begin());
            for (const int2& o : h_overlaps) acct.add(o.x, o.y, row.broad_fp);

            const auto np = Clock::now();
            call_narrow_phase(d, is_vf, d_overlaps, memory_handler, collisions, toi);
            cudaDeviceSynchronize();
            row.narrow_ms += ms_since(np);
        }
        return true;
    }
#endif

    /// One case on the host. Scalable CCD has no host narrow phase, so this
    /// measures the broad phase alone.
    ///
    /// Unreachable: the comparison is on the device only, and the caller refuses
    /// the host space. Kept because it is the only record of how their host
    /// broad phase would be timed, and deleting it would invite someone to
    /// rebuild it differently.
    [[maybe_unused]] bool run_case_host(const MeshArrays& mesh, const bool is_vf, Row& row,
                                        BroadAccounting& acct) {
        namespace sc = scalable_ccd;

        // The boxes are the structure, so they are prep; `sort_and_sweep` sorts
        // and sweeps in one call and cannot be split further, so it is the query.
        const auto prep_begin = Clock::now();
        std::vector<sc::AABB> vertex_boxes, edge_boxes, face_boxes;
        sc::build_vertex_boxes(mesh.V0, mesh.V1, vertex_boxes, kMinDistance);
        if (is_vf) {
            sc::build_face_boxes(vertex_boxes, mesh.F, face_boxes);
        } else {
            sc::build_edge_boxes(vertex_boxes, mesh.E, edge_boxes);
        }
        row.prep_ms += ms_since(prep_begin);

        std::vector<std::pair<int, int>> overlaps;
        int sort_axis = 0;
        const auto broad_begin = Clock::now();
        if (is_vf) {
            sc::sort_and_sweep(vertex_boxes, face_boxes, sort_axis, overlaps);
        } else {
            sc::sort_and_sweep(edge_boxes, sort_axis, overlaps);
        }
        row.broad_ms += ms_since(broad_begin);
        row.queries = overlaps.size();
        for (const auto& [a, b] : overlaps) acct.add(a, b, row.broad_fp);
        return true;
    }

#if defined(SCCD_COMPETITOR_WITH_CUDA)
    /// A pair of primitives whose swept extents are disjoint along x, so the
    /// narrow phase rejects it on its first inclusion test and never splits it.
    /// Padding a query list with copies of it changes nothing but the list's
    /// length -- and `MemoryHandler::handleNarrowPhase` sizes the subdivision
    /// buffer from that length, as `min(available, 2 * n)` domains.
    int2 disjoint_pair(const MeshArrays& m, const bool is_vf) {
        const auto vlo = [&](const int v) { return std::min(m.V0(v, 0), m.V1(v, 0)); };
        const auto vhi = [&](const int v) { return std::max(m.V0(v, 0), m.V1(v, 0)); };

        const auto elem_lo = [&](const int e) {
            return is_vf ? std::min({vlo(m.F(e, 0)), vlo(m.F(e, 1)), vlo(m.F(e, 2))})
                         : std::min(vlo(m.E(e, 0)), vlo(m.E(e, 1)));
        };
        const Eigen::Index n_elems = is_vf ? m.F.rows() : m.E.rows();

        // The left-most primitive of the first kind, then any element entirely
        // to its right.
        int a = 0;
        double a_hi = std::numeric_limits<double>::infinity();
        const Eigen::Index n_first = is_vf ? m.V0.rows() : m.E.rows();
        for (Eigen::Index i = 0; i < n_first; ++i) {
            const int k = static_cast<int>(i);
            const double hi = is_vf ? vhi(k) : std::max(vhi(m.E(k, 0)), vhi(m.E(k, 1)));
            if (hi < a_hi) {
                a_hi = hi;
                a = k;
            }
        }
        for (Eigen::Index e = 0; e < n_elems; ++e) {
            if (elem_lo(static_cast<int>(e)) > a_hi) return make_int2(a, static_cast<int>(e));
        }
        return make_int2(-1, -1);
    }

    /// Whether the narrow phase's answer depends on the length of the list it is
    /// handed, which it must not.
    ///
    /// Runs one step's curated queries through the unmodified narrow phase twice:
    /// bare, and followed by `pad` copies of a pair rejected on its first test.
    /// A correct search reports the same collisions both ways. In the per-query
    /// build the curated queries found are counted; in the default build the
    /// case's earliest time of impact is compared.
    int run_buffer_probe(const fs::path& data_dir, const std::string& scene, const int step,
                         const std::size_t pad, const std::shared_ptr<smesh::Communicator>& comm) {
        const fs::path dd = data_dir / scene;
        auto cases = scan_cases(dd / "boxes");
        MeshArrays mesh;
        if (!load_mesh_arrays(comm, dd, scene, step, mesh)) return EXIT_FAILURE;
        const DeviceMesh d(mesh);

        int shown = 0;
        for (const CaseFile& c : cases) {
            if (case_step(c.key) != step) continue;
            CaseTruth truth;
            if (!load_case_truth(dd, c, truth)) continue;
            if (!normalize_pairs(truth, c.is_vf, mesh.V0.rows(), mesh.E.rows())) continue;

            const int2 filler = disjoint_pair(mesh, c.is_vf);
            if (filler.x < 0) {
                std::cerr << "error: no disjoint pair to pad " << c.key << " with\n";
                return EXIT_FAILURE;
            }

            std::size_t expected = 0;
            for (std::size_t q = 0; q < truth.size(); ++q) {
                expected += (q < truth.mma.size()) ? (truth.mma[q] != 0)
                                                   : (q < truth.root_toi.size() &&
                                                      sccd::is_finite_bits(truth.root_toi[q]));
            }

            for (const std::size_t n_pad : {std::size_t(0), pad}) {
                std::vector<int2> host(truth.size() + n_pad, filler);
                for (std::size_t q = 0; q < truth.size(); ++q) host[q] = make_int2(truth.c0[q], truth.c1[q]);
                const thrust::device_vector<int2> pairs(host.begin(), host.end());

                auto mh = std::make_shared<scalable_ccd::cuda::MemoryHandler>();
                Collisions collisions;
                scalable_ccd::Scalar toi = 1;
                call_narrow_phase(d, c.is_vf, pairs, mh, collisions, toi);
                cudaDeviceSynchronize();

                std::printf("%-8s pad=%-8zu queries=%-6zu expected=%-6zu", c.key.c_str(), n_pad, truth.size(),
                            expected);
                if (kPerQueryBuild) {
                    const std::vector<double> per = scatter_to_queries(truth, collisions);
                    const Accuracy a = score(truth, per);
                    std::printf(" found=%-6llu fn=%-6llu late=%-4llu", static_cast<unsigned long long>(a.toi_n),
                                static_cast<unsigned long long>(a.fn), static_cast<unsigned long long>(a.toi_late));
                }
                std::printf(" toi=%-.9g gt_earliest=%.9g\n", static_cast<double>(toi), truth.gt_earliest);
                std::fflush(stdout);
            }
            ++shown;
        }
        return shown ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    /// Their own entry point, called as their narrow-phase test calls it, on
    /// the same frame pair. It carries one bound from the vertex-face pass into
    /// the edge-edge pass and returns the step's earliest time of impact, which
    /// must equal the minimum of the harness's two case rows for that step.
    int run_their_ccd(const fs::path& data_dir, const std::string& scene, const int step,
                      const std::shared_ptr<smesh::Communicator>& comm) {
        MeshArrays mesh;
        if (!load_mesh_arrays(comm, data_dir / scene, scene, step, mesh)) return EXIT_FAILURE;

        const auto begin = Clock::now();
#if defined(SCALABLE_CCD_TOI_PER_QUERY)
        Collisions collisions;
        const scalable_ccd::Scalar toi =
            scalable_ccd::cuda::ccd(mesh.V0, mesh.V1, mesh.E, mesh.F, kMinDistance, kMaxIterations, tolerance(),
                                    kAllowZeroToi, collisions, /*memory_limit_GB=*/0);
#else
        const scalable_ccd::Scalar toi =
            scalable_ccd::cuda::ccd(mesh.V0, mesh.V1, mesh.E, mesh.F, kMinDistance, kMaxIterations, tolerance(),
                                    kAllowZeroToi, /*memory_limit_GB=*/0);
#endif
        // Every timer in this file stops behind a device synchronisation, as the
        // three in the per-case path do. Their entry point returns a time of
        // impact and so reads back internally, but a timing claim should not
        // rest on that being true of a version we did not write.
        cudaDeviceSynchronize();
        const double elapsed = ms_since(begin);

        std::printf("scalable_ccd::cuda::ccd  %s step %d  toi=%.9g  %.1f ms\n", scene.c_str(), step,
                    static_cast<double>(toi), elapsed);
        return EXIT_SUCCESS;
    }
#endif

    bool want_device() {
        const char* env = std::getenv("SCCD_BENCH_EXECUTION_SPACE");
        if (!env || env[0] == '\0') env = std::getenv("SCCD_EXECUTION_SPACE");
        if (!env || env[0] == '\0') return false;
        std::string v(env);
        std::transform(v.begin(), v.end(), v.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return v == "cuda" || v == "gpu" || v == "device";
    }

}  // namespace

int main(int argc, char** argv) {
    auto ctx = smesh::initialize(argc, argv);

    if (argc == 2 && std::string(argv[1]) == "--header") {
        std::cout << kCsvHeader << "\n";
        return EXIT_SUCCESS;
    }

    if (argc >= 5 && std::string(argv[1]) == "--buffer-probe") {
#if defined(SCCD_COMPETITOR_WITH_CUDA)
        const std::size_t pad = (argc >= 6) ? static_cast<std::size_t>(std::atoll(argv[5])) : 1000000;
        return run_buffer_probe(fs::path(argv[2]), argv[3], std::atoi(argv[4]), pad, ctx->communicator());
#else
        std::cerr << "error: --buffer-probe needs a CUDA build\n";
        return EXIT_FAILURE;
#endif
    }

    if (argc == 5 && std::string(argv[1]) == "--their-ccd") {
#if defined(SCCD_COMPETITOR_WITH_CUDA)
        return run_their_ccd(fs::path(argv[2]), argv[3], std::atoi(argv[4]), ctx->communicator());
#else
        std::cerr << "error: --their-ccd needs a CUDA build\n";
        return EXIT_FAILURE;
#endif
    }

    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <data-dir> <scene> [<scene> ...]\n"
                  << "       " << argv[0] << " --their-ccd <data-dir> <scene> <step>\n"
                  << "       " << argv[0] << " --buffer-probe <data-dir> <scene> <step> [<pad>]\n";
        return EXIT_FAILURE;
    }

    const fs::path data_dir(argv[1]);

    bool device = want_device();
#if !defined(SCCD_COMPETITOR_WITH_CUDA)
    if (device) {
        std::cerr << "warning: device requested but this harness was built without CUDA; "
                     "Scalable CCD's narrow phase is CUDA only, so the host run reports the "
                     "broad phase alone\n";
        device = false;
    }
#endif

    // The mode names the build, so a diagnostic per-query row can never be read
    // as the library's own.
    const bool per_query = device && kPerQueryBuild;
    const std::string mode = per_query ? "scalable-ccd-perquery"
                             : device  ? "scalable-ccd-device"
                                       : "scalable-ccd-host";
    const std::string broadphase = "sweep";

    std::cout << kCsvHeader << "\n";
    std::cout.setf(std::ios::fmtflags(0), std::ios::floatfield);
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

            MeshArrays mesh;
            if (!load_mesh_arrays(ctx->communicator(), dataset_dir, scene, step, mesh)) {
                ++failures;
                continue;
            }

            Row row;
            row.dataset = scene;
            row.mode = mode;
            row.broadphase = broadphase;
            row.key = c.key;
            row.type = c.is_vf ? "vf" : "ee";
            row.gt_earliest = truth.gt_earliest;
            row.root_n = truth.root_toi.size();

            // The curated pairs are in the dataset's global numbering; the broad
            // phase reports per-list indices. The offsets are the mesh's vertex
            // and edge counts, so this follows the mesh read.
            if (!normalize_pairs(truth, c.is_vf, mesh.V0.rows(), mesh.E.rows())) {
                ++failures;
                continue;
            }

            BroadAccounting acct(truth);
            Collisions collisions;
            scalable_ccd::Scalar case_toi = 1;
            bool ok = false;
            if (device) {
#if defined(SCCD_COMPETITOR_WITH_CUDA)
                const DeviceMesh d(mesh);
                ok = run_case_device(mesh, d, c.is_vf, row, acct, collisions, case_toi);
#endif
            } else {
                std::fprintf(stderr,
                             "error: Scalable CCD is compared on the device only. Its narrow phase "
                             "is CUDA, so its host side is a broad phase with no pipeline behind it "
                             "and not something a caller runs end to end.\n");
                return EXIT_FAILURE;
            }
            if (!ok) {
                ++failures;
                continue;
            }

            row.broad_fn = acct.missing();
            if (row.broad_fn != 0) {
                // A curated pair the broad phase never produced is a missed
                // collision, so name it rather than leaving a count: a numbering
                // that does not line up shows up here as everything being missed.
                std::cerr << "warning: " << scene << "/" << c.key << ": broad phase missed "
                          << row.broad_fn << " of " << truth.size() << " curated pairs;";
                int shown = 0;
                for (std::size_t i = 0; i < truth.size() && shown < 4; ++i) {
                    if (!acct.found[i]) {
                        std::cerr << " (" << truth.c0[i] << "," << truth.c1[i] << ")";
                        ++shown;
                    }
                }
                std::cerr << "\n";
            }

            if (per_query) {
                row.acc = score(truth, scatter_to_queries(truth, collisions));
            }

            // The case's earliest time of impact, the number this comparison is
            // about. A host run has no narrow phase and reports none.
            row.earliest_toi = static_cast<double>(case_toi);
            row.has_narrow = device;
            write_row(std::cout, row);
        }
    }

    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
