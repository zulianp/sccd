// Shared scaffolding for the competitor harnesses.
//
// Everything here exists so that a competitor row and an SCCD row describe the
// same measurement: the same cases, the same geometry, the same numbering for
// the pairs, the same ground truth, and the same columns. Only the call into the
// library under test differs, and that lives in the harness itself.
//
// The dataset layout is the one benchmark/scripts/prepare_data.sh produces:
//
//   <data-dir>/<scene>/boxes/<key>/            c0.int32, c1.int32 -- the pairs
//   <data-dir>/<scene>/queries/<key>.csv       the curated queries themselves
//   <data-dir>/<scene>/frames_raw/<step>/      mesh at that step, raw arrays
//   <data-dir>/<scene>/roots/<key>/toi.float64 ground truth times of impact
//   <data-dir>/<scene>/mma_bool/<key>/mma_bool.uint8   ground truth hit/miss

#ifndef SCCD_COMPETITOR_BENCH_HPP
#define SCCD_COMPETITOR_BENCH_HPP

#include "sccd_base.hpp"
#include "sccd_math.hpp"

#include "smesh_mesh.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace sccd_competitor {

    namespace fs = std::filesystem;

    using idx_t = smesh::idx_t;

    /// The one definition of the result schema, and it has to stay identical to
    /// the `kCsvHeader` in benchmark/bench.exe.cpp. A competitor row that carries
    /// different columns cannot be concatenated with an SCCD row, which is the
    /// whole point of running it. `--header` prints this, and the
    /// `competitor_schema` test diffs it against `sccd_bench --header`.
    inline constexpr const char* kCsvHeader =
        "dataset,mode,broadphase,case,type,queries,prep_ms,broad_ms,narrow_ms,query_narrow_ms,fp,fn,broad_fp,"
        "broad_fn,narrow_ms_s1,toi_n,toi_late,toi_max_late,toi_max_early,toi_med_early,s0_late,s0_margin,"
        "s0_toi,gt_earliest,root_n,s1_min";

    /// The narrow-phase parameters benchmark/bench.exe.cpp measures SCCD with.
    ///
    /// Here for reference only. A competitor's knobs are its own, and passing
    /// these through because the names look alike is how the Scalable CCD
    /// harness came to report no collisions at all: its `max_iterations` is a
    /// budget on box checks, not a subdivision depth, and 69 of them starved the
    /// search. Each harness states the values it uses and why.
    inline constexpr int kSccdNarrowPhaseMaxDepth = 69;
    inline constexpr double kSccdNarrowPhaseTolerance = 3e-8;

    struct CaseFile {
        fs::path dir;
        std::string key;
        bool is_vf = false;
    };

    inline std::uint64_t pair_key(const std::int64_t a, const std::int64_t b) {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) | static_cast<std::uint32_t>(b);
    }

    /// Raw arrays, written by the preparation scripts and read with no header:
    /// the file's size divided by sizeof(T) is the count.
    template <typename T>
    bool read_raw(const fs::path& path, std::vector<T>& values) {
        std::ifstream in(path, std::ios::binary | std::ios::ate);
        if (!in) return false;
        const std::streamsize bytes = in.tellg();
        if (bytes < 0 || static_cast<std::size_t>(bytes) % sizeof(T) != 0) return false;
        values.resize(static_cast<std::size_t>(bytes) / sizeof(T));
        in.seekg(0);
        return static_cast<bool>(in.read(reinterpret_cast<char*>(values.data()), bytes));
    }

    /// The step a case belongs to. A key is a frame index with a query kind
    /// glued on -- `108vf` is the vertex-face queries between frames 108 and 109.
    inline int case_step(const std::string& key) {
        std::size_t end = 0;
        while (end < key.size() && std::isdigit(static_cast<unsigned char>(key[end]))) ++end;
        if (end == 0) return -1;
        return std::atoi(key.substr(0, end).c_str());
    }

    /// Every case of a scene, in the order bench.exe.cpp walks them -- so
    /// `SCCD_BENCH_CASE_BEGIN`/`_END` select the same subset here as there.
    ///
    /// A case is a directory under `boxes/` holding the pair arrays, and it only
    /// counts when the curated queries exist beside it: the `.json` files in the
    /// same directory are the raw download, and the converted form is what the
    /// benchmark reads.
    inline std::vector<CaseFile> scan_cases(const fs::path& boxes_dir) {
        std::vector<CaseFile> cases;
        if (!fs::is_directory(boxes_dir)) return cases;
        const fs::path queries_dir = boxes_dir.parent_path() / "queries";
        for (const auto& entry : fs::directory_iterator(boxes_dir)) {
            if (!entry.is_directory()) continue;
            const fs::path dir = entry.path();
            const std::string key = dir.filename().string();
            if (key.size() < 3) continue;
            if (!fs::exists(dir / "c0.int32") || !fs::exists(dir / "c1.int32")) continue;
            if (!fs::exists(queries_dir / (key + ".csv"))) continue;
            CaseFile c;
            c.dir = dir;
            c.key = key;
            c.is_vf = key.compare(key.size() - 2, 2, "vf") == 0;
            cases.push_back(std::move(c));
        }
        // Lexicographic on the key, which is what bench.exe.cpp does. It has to
        // be the same comparison, not merely a deterministic one: the case range
        // slices this list, so a different order means the same
        // SCCD_BENCH_CASE_BEGIN/END select a different set of cases and the
        // competitor is measured on work SCCD never saw. Sorting by step number
        // instead -- which reads more naturally -- left ten of a hundred and
        // twenty cases in common.
        std::sort(cases.begin(), cases.end(),
                  [](const CaseFile& a, const CaseFile& b) { return a.key < b.key; });
        return cases;
    }

    /// The subsample and the case range from the environment, matching
    /// bench.exe.cpp so a sweep chunk covers the same cases whichever binary runs
    /// it. The two compose in that order there, and must here: MAX_CASES means
    /// "this many spread across the trajectory", and the range then slices
    /// whatever that produced.
    inline void apply_case_range(std::vector<CaseFile>& cases) {
        int max_cases = 0;
        if (const char* v = std::getenv("SCCD_BENCH_MAX_CASES")) max_cases = std::atoi(v);
        if (max_cases > 0 && static_cast<int>(cases.size()) > max_cases) {
            // Spread over each query type separately and hold the two in the
            // proportion the dataset has them, exactly as bench.exe.cpp does.
            //
            // Striding the combined list samples the type as well as the
            // trajectory: keys end in "ee" or "vf", so a scene with equal counts
            // of each sorts into interleaved pairs with edge-edge on every even
            // index, and an even stride returns edge-edge alone. That is bad
            // enough in one binary; across two it is worse, because the SCCD
            // rows and the competitor rows of the same chunk then describe
            // different cases and every head-to-head ratio divides timings of
            // different work.
            std::vector<CaseFile> by_type[2];
            for (const CaseFile& c : cases) {
                by_type[c.is_vf ? 1 : 0].push_back(c);
            }

            const std::size_t total = cases.size();
            int want[2];
            want[1] = static_cast<int>(
                (by_type[1].size() * static_cast<std::size_t>(max_cases) + total / 2) / total);
            want[0] = max_cases - want[1];
            for (int t = 0; t < 2; ++t) {
                const int have = static_cast<int>(by_type[t].size());
                if (want[t] > have) {
                    want[1 - t] += want[t] - have;
                    want[t] = have;
                } else if (want[t] == 0 && have > 0 && want[1 - t] > 1) {
                    want[t] = 1;
                    want[1 - t] -= 1;
                }
            }

            std::vector<CaseFile> subset;
            subset.reserve(static_cast<std::size_t>(max_cases));
            for (int t = 0; t < 2; ++t) {
                for (int k = 0; k < want[t]; ++k) {
                    subset.push_back(by_type[t][static_cast<std::size_t>(
                        static_cast<double>(k) * static_cast<double>(by_type[t].size()) / want[t])]);
                }
            }
            // Back into the one order both binaries slice with.
            std::sort(subset.begin(), subset.end(),
                      [](const CaseFile& a, const CaseFile& b) { return a.key < b.key; });
            cases.swap(subset);
        }

        int begin = 0;
        int end = 0;
        if (const char* v = std::getenv("SCCD_BENCH_CASE_BEGIN")) begin = std::atoi(v);
        if (const char* v = std::getenv("SCCD_BENCH_CASE_END")) end = std::atoi(v);
        if (begin <= 0 && end <= 0) return;
        const int total = static_cast<int>(cases.size());
        const int b = std::max(0, std::min(begin, total));
        const int e = (end > 0) ? std::max(b, std::min(end, total)) : total;
        cases = std::vector<CaseFile>(cases.begin() + b, cases.begin() + e);
    }

    /// The benchmark's canonical edge numbering: every undirected edge of the
    /// face list once, ordered by the second endpoint then the first.
    ///
    /// This is a copy of `benchmark_ordered_edges` in bench.exe.cpp and has to
    /// stay one. It is what makes a competitor's edge-edge pairs comparable to
    /// the curated queries without a mapping table -- the competitor is handed
    /// exactly this array as its edge list, so the indices it returns are already
    /// the benchmark's.
    inline std::array<std::vector<idx_t>, 2> benchmark_ordered_edges(const std::shared_ptr<smesh::Mesh>& mesh) {
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
                if (seen.insert(key).second) keys.push_back(key);
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

    /// Where a scene keeps the mesh for a step.
    ///
    /// The naming is per-scene and this has to agree with `frame_path` in
    /// bench.exe.cpp, or the two harnesses measure different geometry while
    /// reporting the same case key.
    inline fs::path frame_path(const fs::path& dataset_dir, const std::string& dataset, const int step) {
        const fs::path frames = dataset_dir / "frames";
        if (dataset == "cloth-ball") {
            return frames / ("cloth_ball" + std::to_string(step) + ".ply");
        }
        if (dataset == "n-body-simulation") {
            return frames / ("balls16_" + std::to_string(step) + ".ply");
        }
        if (dataset == "cloth-funnel") {
            std::string name = std::to_string(step);
            if (name.size() < 3) name.insert(0, 3 - name.size(), '0');
            return frames / (name + ".ply");
        }
        return frames / (std::to_string(step) + ".ply");
    }

    /// The raw-array form of that mesh, which is what smesh reads.
    ///
    /// bench.exe.cpp converts a missing one by shelling out; here it is required
    /// to exist already. benchmark/scripts/prepare_data.sh converts every frame
    /// of every scene it prepares, so a missing directory means the dataset was
    /// not prepared, and saying so is more useful than converting one frame.
    inline fs::path raw_frame_path(const fs::path& dataset_dir, const std::string& dataset, const int step) {
        return dataset_dir / "frames_raw" / frame_path(dataset_dir, dataset, step).stem();
    }

    /// The curated queries of one case, as the pairs they refer to, plus the
    /// ground truth to score against.
    struct CaseTruth {
        std::vector<std::int32_t> c0;      ///< first element of each query pair
        std::vector<std::int32_t> c1;      ///< second element
        std::vector<double> root_toi;      ///< exact time of impact, NaN for a miss
        std::vector<std::uint8_t> mma;     ///< hit/miss, where it exists
        double gt_earliest = 1.0;          ///< earliest finite root over the case

        std::size_t size() const { return c0.size(); }
    };

    inline bool load_case_truth(const fs::path& dataset_dir, const CaseFile& c, CaseTruth& truth) {
        if (!read_raw(c.dir / "c0.int32", truth.c0) || !read_raw(c.dir / "c1.int32", truth.c1)) {
            std::cerr << "error: missing pair arrays for " << c.key << "\n";
            return false;
        }
        if (truth.c0.size() != truth.c1.size()) {
            std::cerr << "error: query pair arrays disagree for " << c.key << "\n";
            return false;
        }

        // Absent ground truth is not an error here -- the accuracy columns simply
        // have nothing to say -- but it is worth not scoring silently, so the
        // caller can see root_n of zero in the row.
        read_raw(dataset_dir / "roots" / c.key / "toi.float64", truth.root_toi);
        read_raw(dataset_dir / "mma_bool" / c.key / "mma_bool.uint8", truth.mma);

        truth.gt_earliest = 1.0;
        for (const double v : truth.root_toi) {
            if (sccd::is_finite_bits(v) && v < truth.gt_earliest) truth.gt_earliest = v;
        }
        return true;
    }

    /// The curated queries' own coordinates: 24 doubles per query, the four
    /// primitives' points at t=0 then the same four at t=1, from
    /// `queries_raw/<key>.f64`.
    ///
    /// These are the coordinates the exact roots were computed for, and the ones
    /// `sccd_bench` scores its per-pair accuracy on. A per-pair comparison scores
    /// every library on them, so the mesh's float32 rounding is not part of it.
    inline bool read_query_points(const fs::path& dataset_dir, const std::string& key,
                                  std::vector<double>& xyz) {
        if (!read_raw(dataset_dir / "queries_raw" / (key + ".f64"), xyz) || xyz.size() % 24 != 0) {
            std::cerr << "error: missing or malformed queries_raw/" << key
                      << ".f64 -- run benchmark/scripts/prepare_data.sh\n";
            return false;
        }
        return true;
    }

    /// The accuracy columns, computed exactly as bench.exe.cpp computes them.
    ///
    /// `toi[i]` is the competitor's time of impact for query `i`, or 1 for "no
    /// collision in this step". The asymmetry that matters is between `toi_late`
    /// and `toi_max_early`: a time of impact after the true one lets a solver
    /// step through the contact and is a failure, while one before it only costs
    /// accuracy.
    struct Accuracy {
        std::uint64_t fp = 0;
        std::uint64_t fn = 0;
        std::uint64_t toi_n = 0;
        std::uint64_t toi_late = 0;
        double toi_max_late = 0.0;
        double toi_max_early = 0.0;
        double toi_med_early = 0.0;
        double s1_min = 1.0;
    };

    inline Accuracy score(const CaseTruth& truth, const std::vector<double>& toi) {
        Accuracy a;
        std::vector<double> early;
        early.reserve(toi.size());

        for (std::size_t i = 0; i < toi.size(); ++i) {
            bool expected = false;
            if (i < truth.mma.size()) {
                expected = truth.mma[i] != 0;
            } else if (i < truth.root_toi.size()) {
                // Bitwise, not std::isfinite: the roots file encodes "no
                // collision" as NaN, and -ffast-math folds isfinite to true,
                // which would turn every correct miss into a reported fn.
                expected = sccd::is_finite_bits(truth.root_toi[i]);
            }

            const bool found = toi[i] < 1.0;
            a.fp += static_cast<std::uint64_t>(found && !expected);
            a.fn += static_cast<std::uint64_t>(!found && expected);

            if (toi[i] < a.s1_min) a.s1_min = toi[i];

            if (i < truth.root_toi.size() && sccd::is_finite_bits(truth.root_toi[i]) && found) {
                const double ref = truth.root_toi[i];
                const double got = toi[i];
                ++a.toi_n;
                if (got > ref) {
                    ++a.toi_late;
                    if (got - ref > a.toi_max_late) a.toi_max_late = got - ref;
                } else {
                    const double e = ref - got;
                    early.push_back(e);
                    if (e > a.toi_max_early) a.toi_max_early = e;
                }
            }
        }

        if (!early.empty()) {
            std::nth_element(early.begin(), early.begin() + early.size() / 2, early.end());
            a.toi_med_early = early[early.size() / 2];
        }
        return a;
    }

    /// Rewrite the curated pairs from the dataset's global element numbering into
    /// the per-list indices a broad phase actually reports.
    ///
    /// The query files number every element of the mesh in one sequence --
    /// vertices in `[0, n_nodes)`, then edges, then faces -- while a broad phase
    /// returns an index into whichever list it was given. Comparing the two
    /// directly finds nothing, and the failure is silent in the worst way: every
    /// curated pair reads as missed by the broad phase and every query as a false
    /// negative, which looks like the competitor being catastrophically wrong
    /// rather than like a harness that forgot to subtract an offset.
    ///
    /// bench.exe.cpp does this in `normalize_pairs`; this has to agree with it.
    inline bool normalize_pairs(CaseTruth& truth, const bool is_vf, const ptrdiff_t n_nodes,
                                const ptrdiff_t n_edges) {
        const ptrdiff_t edge_offset = n_nodes;
        const ptrdiff_t face_offset = n_nodes + n_edges;

        for (std::size_t i = 0; i < truth.c0.size(); ++i) {
            const ptrdiff_t x = truth.c0[i];
            const ptrdiff_t y = truth.c1[i];
            if (is_vf) {
                // Either orientation appears in the files, so take whichever
                // side is the vertex.
                if (x < n_nodes && y >= face_offset) {
                    truth.c0[i] = static_cast<std::int32_t>(x);
                    truth.c1[i] = static_cast<std::int32_t>(y - face_offset);
                } else if (y < n_nodes && x >= face_offset) {
                    truth.c0[i] = static_cast<std::int32_t>(y);
                    truth.c1[i] = static_cast<std::int32_t>(x - face_offset);
                } else {
                    std::cerr << "error: invalid vertex-face pair (" << x << ", " << y
                              << ") against n_nodes=" << n_nodes << " n_edges=" << n_edges << "\n";
                    return false;
                }
            } else {
                if (x < edge_offset || y < edge_offset || x >= face_offset || y >= face_offset) {
                    std::cerr << "error: invalid edge-edge pair (" << x << ", " << y
                              << ") against n_nodes=" << n_nodes << " n_edges=" << n_edges << "\n";
                    return false;
                }
                truth.c0[i] = static_cast<std::int32_t>(x - edge_offset);
                truth.c1[i] = static_cast<std::int32_t>(y - edge_offset);
            }
        }
        return true;
    }

    /// Which curated pairs a broad phase produced, and how many of the pairs it
    /// produced were not curated ones.
    ///
    /// This is what the `queries`, `broad_fp` and `broad_fn` columns mean in
    /// bench.exe.cpp, and they are easy to misread: `queries` is the number of
    /// candidates the narrow phase was handed, not the number of curated queries,
    /// and `broad_fp` counts candidates that carry no known contact. A case can
    /// therefore have one curated query and ninety thousand candidates -- the
    /// published cloth-ball rows do.
    struct BroadAccounting {
        std::unordered_map<std::uint64_t, std::size_t> slot;
        std::vector<std::uint8_t> found;

        explicit BroadAccounting(const CaseTruth& truth) : found(truth.size(), 0) {
            slot.reserve(truth.size() * 2 + 1);
            for (std::size_t i = 0; i < truth.size(); ++i) {
                slot.emplace(pair_key(truth.c0[i], truth.c1[i]), i);
            }
        }

        /// Account for one candidate pair. The curated list fixes an order for
        /// each pair and a competitor has no reason to use the same one, so both
        /// are tried before calling it a false positive.
        void add(const int a, const int b, std::uint64_t& false_positives) {
            auto it = slot.find(pair_key(a, b));
            if (it == slot.end()) it = slot.find(pair_key(b, a));
            if (it == slot.end()) {
                ++false_positives;
            } else {
                found[it->second] = 1;
            }
        }

        std::uint64_t missing() const {
            std::uint64_t n = 0;
            for (const std::uint8_t f : found) n += (f == 0);
            return n;
        }
    };

    /// Turn the candidate pairs a library reported as colliding into one time of
    /// impact per curated query.
    ///
    /// A query nothing reported is a miss, which is 1. A pair reported more than
    /// once keeps the earliest, because that is the answer a solver would act on.
    /// Both harnesses need this -- one maps a device collision list, the other a
    /// host one -- and two copies of it would drift.
    inline std::vector<double> scatter_to_queries(
        const CaseTruth& truth, const std::vector<std::tuple<int, int, double>>& collisions) {
        std::unordered_map<std::uint64_t, std::size_t> slot;
        slot.reserve(truth.size() * 2 + 1);
        for (std::size_t i = 0; i < truth.size(); ++i) {
            slot.emplace(pair_key(truth.c0[i], truth.c1[i]), i);
        }

        std::vector<double> toi(truth.size(), 1.0);
        for (const auto& [a, b, t] : collisions) {
            // The curated list fixes an order for each pair; a library has no
            // reason to use the same one, so try both.
            auto it = slot.find(pair_key(a, b));
            if (it == slot.end()) it = slot.find(pair_key(b, a));
            if (it == slot.end()) continue;  // a pair outside the curated set
            if (t < toi[it->second]) toi[it->second] = t;
        }
        return toi;
    }

    /// One CSV row, in the schema above.
    struct Row {
        std::string dataset;
        std::string mode;
        std::string broadphase;
        std::string key;
        std::string type;
        std::size_t queries = 0;
        double prep_ms = 0.0;
        double broad_ms = 0.0;
        double narrow_ms = 0.0;
        double query_narrow_ms = 0.0;
        std::uint64_t broad_fp = 0;
        std::uint64_t broad_fn = 0;
        double narrow_ms_s1 = 0.0;
        /// Whether a narrow phase ran at all. A broad-phase-only row has no time
        /// of impact to report, and reporting 1 for it would make `s0_late` fire
        /// on every row -- "no collision found" read as "found one after the true
        /// contact", which is the one failure the column exists to catch. False
        /// leaves the time-of-impact columns neutral instead.
        bool has_narrow = true;
        double earliest_toi = 1.0;
        double gt_earliest = 1.0;
        std::size_t root_n = 0;
        Accuracy acc;
    };

    inline void write_row(std::ostream& out, const Row& r) {
        const int s0_late = (r.has_narrow && r.earliest_toi > r.gt_earliest) ? 1 : 0;
        const double s0_margin = r.has_narrow ? (r.gt_earliest - r.earliest_toi) : 0.0;
        out << r.dataset << ',' << r.mode << ',' << r.broadphase << ',' << r.key << ',' << r.type << ','
            << r.queries << ',' << r.prep_ms << ',' << r.broad_ms << ',' << r.narrow_ms << ','
            << r.query_narrow_ms << ',' << r.acc.fp << ',' << r.acc.fn << ',' << r.broad_fp << ','
            << r.broad_fn << ',' << r.narrow_ms_s1 << ',' << r.acc.toi_n << ',' << r.acc.toi_late << ','
            << r.acc.toi_max_late << ',' << r.acc.toi_max_early << ',' << r.acc.toi_med_early << ','
            << s0_late << ',' << s0_margin << ',' << r.earliest_toi << ',' << r.gt_earliest << ','
            << r.root_n << ',' << r.acc.s1_min << '\n';
    }

}  // namespace sccd_competitor

#endif  // SCCD_COMPETITOR_BENCH_HPP
